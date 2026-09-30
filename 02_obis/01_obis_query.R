#-----------------------------------------------------------------------------
# Author: Paola Gabasa (@paolagab)
# 01_obis_query.R
# Description: Query OBIS occurrence data (parquet dump) via DuckDB and
#              rasterize to a 1° annual grid. Processes years in chunks
#              (N_YEARS_PER_CHUNK) instead of one-by-one (too many queries)
#              or all at once (RAM spike). Defaults calibrated for a laptop
#              with 4 physical cores / ~16GB RAM -- adjust N_THREADS and
#              MEMORY_LIMIT to your own hardware.
#-----------------------------------------------------------------------------

library(DBI)
library(duckdb)
library(dplyr)
library(tidyr)
library(terra)

source("R/data_paths.R")

#----------------------------------------------------
# 0. Hardware-calibrated parameters
#----------------------------------------------------
N_THREADS <- 4          # REAL physical cores (not logical/hyperthreaded)
MEMORY_LIMIT <- "6GB"   # leaves headroom for RStudio + OS; lower to 3-4GB if RAM is tight
N_YEARS_PER_CHUNK <- 5  # chunk size: higher = fewer queries but higher peak memory

#----------------------------------------------------
# 1. Input parquet files
#----------------------------------------------------
input_dir <- rstudioapi::selectDirectory()  # interactive by design: OBIS export location varies
input_dir <- normalizePath(input_dir, winslash = "/")

parquet_files <- list.files(input_dir, pattern = "\\.parquet$", full.names = TRUE)
if (length(parquet_files) == 0) stop("No .parquet files found in: ", input_dir)
print(parquet_files)

if (!dir.exists(query_dir)) dir.create(query_dir, recursive = TRUE)  # query_dir from data_paths.R

years <- 2000:2024

#----------------------------------------------------
# 2. DuckDB connection (hardware-calibrated)
#----------------------------------------------------
con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)


dbExecute(con, sprintf("PRAGMA threads=%d", N_THREADS))
dbExecute(con, sprintf("PRAGMA memory_limit='%s'", MEMORY_LIMIT))

files_sql <- paste(sprintf("'%s'", parquet_files), collapse = ", ")
dbExecute(con, sprintf(
  "CREATE VIEW occ AS SELECT * FROM read_parquet([%s], union_by_name = true);",
  files_sql
))

available_cols <- dbGetQuery(con, "DESCRIBE occ")$column_name
has_uncertainty <- "coordinateUncertaintyInMeters" %in% available_cols
has_basis       <- "basisOfRecord" %in% available_cols

uncertainty_clause <- if (has_uncertainty) {
  "AND (coordinateUncertaintyInMeters IS NULL OR coordinateUncertaintyInMeters <= 100000)"
} else ""
basis_clause <- if (has_basis) {
  "AND (basisOfRecord IS NULL OR LOWER(basisOfRecord) NOT LIKE '%fossil%')"
} else ""

#----------------------------------------------------
# 3. Ocean mask and rasterizing function
#----------------------------------------------------
ocean_mask <- rast(temp_mask)  # temp_mask defined in data_paths.R (from study_area.R output)

obis_to_raster <- function(df, ocean_mask) {
  if (nrow(df) == 0) return(NULL)
  coords <- as.matrix(df[, c("decimalLongitude", "decimalLatitude")])
  v <- vect(coords, type = "points", crs = "EPSG:4326")
  v <- project(v, crs(ocean_mask))
  r <- rasterize(v, ocean_mask, field = 1, fun = "sum", background = 0)
  mask(r, ocean_mask)
}

#----------------------------------------------------
# 4. Loop over BLOCKS of N_YEARS_PER_CHUNK years
#    (2 queries per block, not one per year)
#----------------------------------------------------
year_chunks <- split(years, ceiling(seq_along(years) / N_YEARS_PER_CHUNK))
message(sprintf("Processing %d years in %d block(s) of up to %d years each...",
                length(years), length(year_chunks), N_YEARS_PER_CHUNK))

obis_list <- list()
annual_stats_list <- list()

for (chunk_i in seq_along(year_chunks)) {
  chunk_years <- year_chunks[[chunk_i]]
  y_min <- min(chunk_years); y_max <- max(chunk_years)
  
  message(sprintf("\n[Block %d/%d] Years %d-%d...", chunk_i, length(year_chunks), y_min, y_max))
  
  # --- Raw record count per year, THIS BLOCK ONLY ---
  n_raw_chunk <- dbGetQuery(con, sprintf("
    SELECT date_year, COUNT(*) AS n_raw
    FROM occ
    WHERE date_year BETWEEN %d AND %d
      AND decimalLongitude IS NOT NULL
      AND decimalLatitude IS NOT NULL
    GROUP BY date_year
  ", y_min, y_max))
  
  # --- QC-filtered, deduplicated events, THIS BLOCK ONLY ---
  df_events_chunk <- dbGetQuery(con, sprintf("
    SELECT DISTINCT
      date_year,
      eventDate,
      decimalLongitude,
      decimalLatitude
    FROM occ
    WHERE date_year BETWEEN %d AND %d
      AND decimalLongitude BETWEEN -180 AND 180
      AND decimalLatitude BETWEEN -90 AND 90
      %s
      %s
  ", y_min, y_max, uncertainty_clause, basis_clause))
  
  message(sprintf("  -> %d unique events extracted in this block", nrow(df_events_chunk)))
  
  # --- Block statistics ---
  n_events_chunk <- df_events_chunk %>% count(date_year, name = "n_sampling_events")
  
  stats_chunk <- data.frame(year = chunk_years) %>%
    left_join(n_raw_chunk,    by = c("year" = "date_year")) %>%
    left_join(n_events_chunk, by = c("year" = "date_year")) %>%
    mutate(
      n_raw             = replace_na(n_raw, 0),
      n_sampling_events = replace_na(n_sampling_events, 0),
      n_records_qc      = n_sampling_events
    ) %>%
    rename(n_records_raw = n_raw) %>%
    select(year, n_records_raw, n_records_qc, n_sampling_events)
  
  annual_stats_list[[chunk_i]] <- stats_chunk
  
  # --- Rasterize each year in the block (fast, in-memory) ---
  df_split_chunk <- split(df_events_chunk, df_events_chunk$date_year)
  
  for (yr in chunk_years) {
    df_yr <- df_split_chunk[[as.character(yr)]]
    if (is.null(df_yr) || nrow(df_yr) == 0) {
      message(sprintf("     No records for %d. Skipping...", yr))
      next
    }
    r_year <- obis_to_raster(df_yr, ocean_mask)
    if (!is.null(r_year)) {
      names(r_year) <- paste0("Y", yr)
      obis_list[[as.character(yr)]] <- r_year
    }
  }
  
  # --- Free memory before the next block ---
  rm(df_events_chunk, df_split_chunk, n_raw_chunk, n_events_chunk)
  gc(verbose = FALSE)
}

#----------------------------------------------------
# 5. Consolidation and export
#----------------------------------------------------
obis_list <- obis_list[!sapply(obis_list, is.null)]
if (length(obis_list) == 0) stop("No year produced rasterizable data.")

annual_stats <- bind_rows(annual_stats_list)

obis_stack <- rast(obis_list)ento
terra::plot(obis_stack)  # interactive visual QC only

write.csv(annual_stats, obis_annual_processing_stats, row.names = FALSE)
message("Statistics saved to: ", obis_annual_processing_stats)

writeRaster(obis_stack, filename = obis_count, overwrite = TRUE)
message("Raster saved to: ", obis_count)
























#-----------------------------------------------------------------------------
# Author: Paola Gabasa (@paolagab)
# 01_obis_query.R
# Description: Query OBIS occurrence data (parquet dump) via DuckDB and
#              rasterize to a 1° annual grid. Uses the full historical year
#              range available in the data, capped at 2024 (last complete
#              year -- any partial 2025+ data is dropped). QC applied:
#              excludes fossil and preserved-specimen records, coordinate
#              uncertainty > 100 km, and points falling on land (checked
#              against the ocean mask BEFORE counting sampling events, not
#              just implicitly at rasterization). Processes years in chunks
#              (N_YEARS_PER_CHUNK) instead of one-by-one (too many queries)
#              or all at once (RAM spike). Defaults calibrated for a laptop
#              with 4 physical cores / ~16GB RAM -- adjust N_THREADS and
#              MEMORY_LIMIT to your own hardware.
#-----------------------------------------------------------------------------

library(DBI)
library(duckdb)
library(dplyr)
library(tidyr)
library(terra)

source("R/data_paths.R")

#----------------------------------------------------
# 0. Hardware-calibrated parameters
#----------------------------------------------------
N_THREADS <- 4          # REAL physical cores (not logical/hyperthreaded)
MEMORY_LIMIT <- "6GB"   # leaves headroom for RStudio + OS; lower to 3-4GB if RAM is tight
N_YEARS_PER_CHUNK <- 5  # chunk size for the recent/dense period (after CHUNK_BREAKPOINTS)

# Additional split points (ascending) for coarser chunks over old, sparse
# years -- avoids e.g. splitting min_yr:1950 into dozens of near-empty
# 5-year blocks when the actual record volume there is tiny. Years after
# the last breakpoint are still split into N_YEARS_PER_CHUNK-sized blocks.
CHUNK_BREAKPOINTS <- c(1950, 1990)

LAST_COMPLETE_YEAR <- 2024  # hard cap: excludes any partial/incomplete later years

#----------------------------------------------------
# 1. Input parquet files
#----------------------------------------------------
input_dir <- rstudioapi::selectDirectory()  # interactive by design: OBIS export location varies
input_dir <- normalizePath(input_dir, winslash = "/")

parquet_files <- list.files(input_dir, pattern = "\\.parquet$", full.names = TRUE)
if (length(parquet_files) == 0) stop("No .parquet files found in: ", input_dir)
print(parquet_files)

if (!dir.exists(query_dir)) dir.create(query_dir, recursive = TRUE)  # query_dir from data_paths.R

#----------------------------------------------------
# 2. DuckDB connection (hardware-calibrated)
#----------------------------------------------------
con <- dbConnect(duckdb::duckdb())
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

dbExecute(con, sprintf("PRAGMA threads=%d", N_THREADS))
dbExecute(con, sprintf("PRAGMA memory_limit='%s'", MEMORY_LIMIT))

files_sql <- paste(sprintf("'%s'", parquet_files), collapse = ", ")
dbExecute(con, sprintf(
  "CREATE VIEW occ AS SELECT * FROM read_parquet([%s], union_by_name = true);",
  files_sql
))

# Full historical range available, capped at LAST_COMPLETE_YEAR (drops any
# partial/incomplete more-recent years, e.g. 2025 data collected mid-year).
year_range <- dbGetQuery(con, sprintf(
  "SELECT MIN(date_year) AS min_yr FROM occ WHERE date_year <= %d", LAST_COMPLETE_YEAR
))
years <- year_range$min_yr:LAST_COMPLETE_YEAR
message(sprintf("Processing the full available range: %d-%d.", min(years), max(years)))

available_cols <- dbGetQuery(con, "DESCRIBE occ")$column_name
has_uncertainty <- "coordinateUncertaintyInMeters" %in% available_cols
has_basis       <- "basisOfRecord" %in% available_cols

uncertainty_clause <- if (has_uncertainty) {
  "AND (coordinateUncertaintyInMeters IS NULL OR coordinateUncertaintyInMeters <= 100000)"
} else ""
basis_clause <- if (has_basis) {
  "AND (basisOfRecord IS NULL OR (
     LOWER(basisOfRecord) NOT LIKE '%fossil%'
     AND LOWER(basisOfRecord) NOT LIKE '%specimen%'
   ))"
  # Excludes FossilSpecimen, PreservedSpecimen and LivingSpecimen -- natural
  # history collection records, not field sampling effort. Keeps
  # HumanObservation, MachineObservation, MaterialSample, etc.
} else ""

#----------------------------------------------------
# 3. Ocean mask and rasterizing function
#----------------------------------------------------
ocean_mask <- rast(temp_mask)  # temp_mask defined in data_paths.R (from study_area.R output)

# Drops points falling on land, checked against the ocean mask. Applied to
# the QC'd event table BEFORE counting sampling events, so the reported
# event counts already reflect ocean-only records (not just the final
# raster, which would mask land out anyway but AFTER the stats are computed).
drop_land_points <- function(df, ocean_mask) {
  if (nrow(df) == 0) return(df)
  v <- vect(df, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  v <- project(v, crs(ocean_mask))
  on_ocean <- !is.na(terra::extract(ocean_mask, v)[, 2])
  df[on_ocean, ]
}

obis_to_raster <- function(df, ocean_mask) {
  if (nrow(df) == 0) return(NULL)
  coords <- as.matrix(df[, c("decimalLongitude", "decimalLatitude")])
  v <- vect(coords, type = "points", crs = "EPSG:4326")
  v <- project(v, crs(ocean_mask))
  r <- rasterize(v, ocean_mask, field = 1, fun = "sum", background = 0)
  mask(r, ocean_mask)
}

# Builds year chunks: one coarse chunk per segment defined by
# CHUNK_BREAKPOINTS (for old, sparse years), then N_YEARS_PER_CHUNK-sized
# blocks for everything after the last breakpoint (the recent, dense period).
build_year_chunks <- function(years, breakpoints, chunk_size) {
  bounds <- c(min(years) - 1, sort(breakpoints), max(years))
  chunks <- list()
  for (i in seq_len(length(bounds) - 1)) {
    seg_years <- years[years > bounds[i] & years <= bounds[i + 1]]
    if (length(seg_years) == 0) next
    if (i < length(bounds) - 1) {
      chunks[[length(chunks) + 1]] <- seg_years  # coarse: whole segment as one chunk
    } else {
      chunks <- c(chunks, split(seg_years, ceiling(seq_along(seg_years) / chunk_size)))
    }
  }
  chunks
}

#----------------------------------------------------
# 4. Loop over BLOCKS of N_YEARS_PER_CHUNK years
#    (2 queries per block, not one per year)
#----------------------------------------------------
year_chunks <- build_year_chunks(years, CHUNK_BREAKPOINTS, N_YEARS_PER_CHUNK)
message(sprintf("Processing %d years in %d block(s) (coarse blocks before %s, %d-year blocks after).",
                length(years), length(year_chunks),
                paste(CHUNK_BREAKPOINTS, collapse = "/"), N_YEARS_PER_CHUNK))

obis_list <- list()
annual_stats_list <- list()

for (chunk_i in seq_along(year_chunks)) {
  chunk_years <- year_chunks[[chunk_i]]
  y_min <- min(chunk_years); y_max <- max(chunk_years)
  
  message(sprintf("\n[Block %d/%d] Years %d-%d...", chunk_i, length(year_chunks), y_min, y_max))
  
  # --- Raw record count per year, THIS BLOCK ONLY ---
  n_raw_chunk <- dbGetQuery(con, sprintf("
    SELECT date_year, COUNT(*) AS n_raw
    FROM occ
    WHERE date_year BETWEEN %d AND %d
      AND decimalLongitude IS NOT NULL
      AND decimalLatitude IS NOT NULL
    GROUP BY date_year
  ", y_min, y_max))
  
  # --- QC-filtered, deduplicated events, THIS BLOCK ONLY ---
  df_events_chunk <- dbGetQuery(con, sprintf("
    SELECT DISTINCT
      date_year,
      eventDate,
      decimalLongitude,
      decimalLatitude
    FROM occ
    WHERE date_year BETWEEN %d AND %d
      AND decimalLongitude BETWEEN -180 AND 180
      AND decimalLatitude BETWEEN -90 AND 90
      %s
      %s
  ", y_min, y_max, uncertainty_clause, basis_clause))
  
  n_before_land_filter <- nrow(df_events_chunk)
  df_events_chunk <- drop_land_points(df_events_chunk, ocean_mask)
  message(sprintf("  -> %d unique events after QC, %d after also excluding land points",
                  n_before_land_filter, nrow(df_events_chunk)))
  
  # --- Block statistics ---
  n_events_chunk <- df_events_chunk %>% count(date_year, name = "n_sampling_events")
  
  stats_chunk <- data.frame(year = chunk_years) %>%
    left_join(n_raw_chunk,    by = c("year" = "date_year")) %>%
    left_join(n_events_chunk, by = c("year" = "date_year")) %>%
    mutate(
      n_raw             = replace_na(n_raw, 0),
      n_sampling_events = replace_na(n_sampling_events, 0),
      n_records_qc      = n_sampling_events
    ) %>%
    rename(n_records_raw = n_raw) %>%
    select(year, n_records_raw, n_records_qc, n_sampling_events)
  
  annual_stats_list[[chunk_i]] <- stats_chunk
  
  # --- Rasterize each year in the block (fast, in-memory) ---
  df_split_chunk <- split(df_events_chunk, df_events_chunk$date_year)
  
  for (yr in chunk_years) {
    df_yr <- df_split_chunk[[as.character(yr)]]
    if (is.null(df_yr) || nrow(df_yr) == 0) {
      message(sprintf("     No records for %d. Skipping...", yr))
      next
    }
    r_year <- obis_to_raster(df_yr, ocean_mask)
    if (!is.null(r_year)) {
      names(r_year) <- paste0("Y", yr)
      obis_list[[as.character(yr)]] <- r_year
    }
  }
  
  # --- Free memory before the next block ---
  rm(df_events_chunk, df_split_chunk, n_raw_chunk, n_events_chunk)
  gc(verbose = FALSE)
}


#----------------------------------------------------
# 5. Consolidation and export
#----------------------------------------------------
obis_list <- obis_list[!sapply(obis_list, is.null)]
if (length(obis_list) == 0) stop("No year produced rasterizable data.")

annual_stats <- bind_rows(annual_stats_list) %>% arrange(year)

obis_stack <- rast(obis_list)

# ----------------------------------------------------
# 5b. Restrict to years with gap-free annual reporting
# ----------------------------------------------------
# The temporal window is restricted to years from which continuous
# (gap-free) annual reporting is observed, excluding a small number of
# isolated early records with implausible or unverifiable dates (a few
# early years appear as isolated spikes separated by one or more years
# with zero raw records at all -- consistent with corrupted/mis-parsed
# dates, e.g. a partial "month-day" field read as a 3-4 digit year, rather
# than genuine historical sampling gaps). The cutoff is derived from the
# data itself, not an arbitrary round number: it is the first year after
# which no year has zero raw records.
zero_idx <- which(annual_stats$n_records_raw == 0)
first_continuous_year <- if (length(zero_idx) == 0) {
  min(annual_stats$year)
} else {
  annual_stats$year[max(zero_idx) + 1]
}
message(sprintf("First year of gap-free annual coverage: %d. Restricting the temporal window accordingly.",
                first_continuous_year))

n_years_dropped <- sum(annual_stats$year < first_continuous_year)
n_events_dropped <- sum(annual_stats$n_sampling_events[annual_stats$year < first_continuous_year])
if (n_years_dropped > 0) {
  message(sprintf("Dropping %d early year(s) before %d (%d sampling events).",
                  n_years_dropped, first_continuous_year, n_events_dropped))
}

annual_stats <- annual_stats %>% filter(year >= first_continuous_year)
years <- years[years >= first_continuous_year]

# Match by position against obis_list's own names (set explicitly as
# as.character(yr) when each layer was created), not obis_stack's internal
# layer names -- those have proven unreliable to depend on in this project
# (lost on .nc write/read round-trips elsewhere in the pipeline).
included_years <- as.integer(names(obis_list))
stopifnot(!anyNA(included_years))
kept_layers <- included_years %in% years
obis_stack <- obis_stack[[kept_layers]]
names(obis_stack) <- paste0("Y", included_years[kept_layers])  # restore readable names explicitly

terra::plot(obis_stack)  # interactive visual QC only

write.csv(annual_stats, obis_annual_processing_stats, row.names = FALSE)
message("Statistics saved to: ", obis_annual_processing_stats)

writeRaster(obis_stack, filename = obis_count, overwrite = TRUE)
message("Raster saved to: ", obis_count)

#----------------------------------------------------
# 6. Summary sentence for the manuscript
#----------------------------------------------------
total_sampling_events <- sum(annual_stats$n_sampling_events)

cat(sprintf(
  "\nThis yielded %s distinct sampling events spanning %d-%d.\n",
  format(total_sampling_events, big.mark = ","), min(years), max(years)
))
