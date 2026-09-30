#-----------------------------------------------------------------------------
# 02_eez_high-seas_breakdown.R
# Description: Within the overlap footprint (undersampled OBIS areas x high
#              GFW fishing effort, from 01_overlap.R), downloads drifting
#              longline fishing effort broken down by flag state from GFW,
#              restricted to the overlap region. The query region is built
#              as a coarse 10x10 degree grid intersected with the overlap
#              footprint and dissolved into contiguous polygons, downloaded
#              region-by-region with local caching and retry/backoff, to
#              avoid the HTTP 503 errors a single large, scattered query
#              polygon previously produced.
#
#              Each overlap cell is then classified as falling within an
#              Exclusive Economic Zone (EEZ, + country) or the High Seas
#              (Marine Regions boundaries), with antimeridian-crossing
#              polygons (e.g. Russia, Fiji) split via st_break_antimeridian()
#              before classification, and boundary/sliver-gap cells resolved
#              via nearest-feature assignment. Produces the top N flag
#              states by mean annual fishing effort, ranked separately for
#              EEZ and High Seas (different governance pathways).
#-----------------------------------------------------------------------------



library(sf)
library(terra)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gfwr)      # Global Fishing Watch API



# ----------------------------------------------------
# 1. Load parameters and overlap cells (from 01_overlap.R)
# ----------------------------------------------------

GEARTYPE <- "drifting_longlines"
years <- 2016:2024

TOP_N_FLAGS <- 10

ocean_mask <- rast(temp_mask)
PROJ <- crs(ocean_mask)

overlap_cells_path <- paste(overlap_dir, sprintf("overlap_cells_%s.rds", GEARTYPE), sep = "/")
overlap_sf <- readRDS(overlap_cells_path)
message(sprintf("Overlap cells loaded: %d", nrow(overlap_sf)))

cell_res <- res(ocean_mask)  # c(dx, dy) in meters, Mollweide

# ----------------------------------------------------
# 2. Geometry preparation: Grid generation & grouping
# ----------------------------------------------------
message("Building study area grid for GFW API...")

# 1. Get global bounding box in WGS84
bbox_wgs <- overlap_sf %>%
  st_transform(4326) %>%
  st_bbox()

# 2. Create regular grid over the bounding box
sf_use_s2(FALSE)

# Create 10x10 degree regular grid 
grid_all <- st_make_grid(
  st_as_sfc(bbox_wgs), 
  cellsize = c(10, 10), 
  square = TRUE
)

# 3. Filter cells containing points from the study area
pts_wgs <- overlap_sf %>%
  st_transform(4326) %>%
  st_geometry()

intersect_idx <- st_intersects(grid_all, pts_wgs, sparse = TRUE)
grid_study_area <- grid_all[lengths(intersect_idx) > 0]

# 4. Group contiguous grid cells into continuous polygons (reduces from ~1600 cells to 71 regions)
grid_grouped <- grid_study_area %>%
  st_union() %>%
  st_cast("POLYGON") %>%
  st_sf() %>%
  mutate(grid_id = row_number())

sf_use_s2(TRUE)

# Transform overlap_sf to WGS84 to align with the grid
overlap_wgs <- st_transform(overlap_sf, 4326)

ggplot() +
  # Draw the 71 fused continuous regions (grouped polygons)
  geom_sf(data = grid_grouped, fill = NA, color = "red", linewidth = 0.5) +
  # Draw the original study area cells/points (Overlap)
  geom_sf(data = overlap_wgs, color = "blue", size = 0.6, alpha = 0.5) +
  labs(
    title = sprintf("GFW Query Area (%d grouped polygonal regions)", nrow(grid_grouped)),
    subtitle = "Red: Continuous regions for API | Blue: Original study area (Overlap)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

message(sprintf("Study area successfully divided into %d continuous regions for API query.", nrow(grid_grouped)))

# ----------------------------------------------------
# 2b. Helper function: Download with retries and backoff
# ----------------------------------------------------
download_with_retry <- function(region_poly, y, geartype_filter, max_attempts = 3) {
  attempt <- 1
  res <- NULL
  
  while (attempt <= max_attempts && is.null(res)) {
    tryCatch({
      Sys.sleep(1) # Preventative pause to avoid API rate limits
      
      res <- gfw_ais_fishing_hours(
        spatial_resolution  = "LOW",
        temporal_resolution = "YEARLY",
        group_by            = "FLAGANDGEARTYPE",
        start_date          = paste0(y, "-01-01"),
        end_date            = paste0(y, "-12-31"),
        region_source       = "USER_SHAPEFILE",
        region              = region_poly
      )
      
      if (!is.null(res) && nrow(res) > 0) {
        res <- res %>%
          rename(flag = Flag, geartype = Geartype) %>%
          filter(geartype == geartype_filter) %>%
          mutate(year = y)
      } else {
        res <- data.frame() # Empty dataframe if no records found
      }
    }, error = function(e) {
      message(sprintf("   ⚠️ Attempt %d/%d failed for year %d: %s. Retrying...", 
                      attempt, max_attempts, y, e$message))
      attempt <<- attempt + 1
      Sys.sleep(3 * attempt) # Incremental backoff pause
    })
  }
  return(res)
}

# ----------------------------------------------------
# 3. Download effort by flag (REGION BY REGION WITH CACHE)
# ----------------------------------------------------
dir.create("gfw_cache", showWarnings = FALSE)
message(sprintf("Downloading effort across %d regions...", nrow(grid_grouped)))

effort_by_flag_raw <- map_dfr(seq_len(nrow(grid_grouped)), function(i) {
  cache_file <- sprintf("gfw_cache/region_%d.rds", i)
  
  # Load from cache if region was previously downloaded
  if (file.exists(cache_file)) {
    message(sprintf("Region %d/%d loaded from cache.", i, nrow(grid_grouped)))
    return(readRDS(cache_file))
  }
  
  region_poly <- grid_grouped[i, ]
  message(sprintf("Downloading Region %d/%d...", i, nrow(grid_grouped)))
  
  df_region <- map_dfr(years, function(y) {
    message(sprintf("  -> Year %d...", y))
    download_with_retry(region_poly, y, GEARTYPE)
  })
  
  saveRDS(df_region, cache_file)
  return(df_region)
}) %>% distinct()

gfw_effort_by_flag_raw_path <- paste(overlap_dir, sprintf("gfw_effort_by_flag_raw_%s.rds", GEARTYPE), sep = "/")
saveRDS(effort_by_flag_raw, gfw_effort_by_flag_raw_path)

#effort_by_flag_raw <- readRDS(gfw_effort_by_flag_raw_path)

# ----------------------------------------------------
# 4. Aggregation
# ----------------------------------------------------
message("Aggregating mean annual hours per cell and flag...")

effort_flag_agg <- effort_by_flag_raw %>%
  group_by(Lon, Lat, flag, year) %>%
  summarise(yearly_hours = sum(`Apparent Fishing Hours`, na.rm = TRUE), .groups = "drop") %>%
  complete(nesting(Lon, Lat, flag), year = years, fill = list(yearly_hours = 0)) %>%
  group_by(Lon, Lat, flag) %>%
  summarise(mean_yearly_hours = mean(yearly_hours, na.rm = TRUE), .groups = "drop")
# ----------------------------------------------------
# 5. Rasterize each flag onto the 1-degree grid
# ----------------------------------------------------
message("Rasterizing by flag and filtering to overlap cells...")

n_before <- length(unique(effort_flag_agg$flag))
effort_flag_agg <- effort_flag_agg %>% filter(!is.na(flag) & flag != "")
n_after <- length(unique(effort_flag_agg$flag))

if (n_before != n_after) {
  message(sprintf("Excluded %d NA/blank flag(s) (vessels with no flag reported via AIS).", n_before - n_after))
}

flags <- unique(effort_flag_agg$flag)

# Rebuilt from x/y rather than reusing overlap_sf's own geometry column,
# which is known to contain empty points for a subset of cells.
overlap_sf_clean <- overlap_sf %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("x", "y"), crs = PROJ, remove = FALSE)

flag_cell_values <- map_dfr(flags, function(fl) {
  df_fl <- effort_flag_agg %>% filter(flag == fl)
  if (nrow(df_fl) == 0) return(NULL)
  
  v <- vect(df_fl, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  v <- project(v, PROJ)
  r <- rasterize(v, ocean_mask, field = "mean_yearly_hours", fun = "sum", background = 0)
  r <- mask(r, ocean_mask)
  
  vals <- terra::extract(r, vect(overlap_sf_clean))[, 2]
  data.frame(x = overlap_sf_clean$x, y = overlap_sf_clean$y, flag = fl, mean_yearly_hours = vals)
}) %>%
  filter(mean_yearly_hours > 0)

# ----------------------------------------------------
# 6. Classify each overlap cell: EEZ (+ country) vs High Seas
# ----------------------------------------------------
message("Classifying overlap cells: EEZ (country) vs High Seas...")

sf_use_s2(FALSE)

EEZ_COUNTRY_FIELD <- "SOVEREIGN1"

# 1. Load EEZ, wrap/break at antimeridian, and then reproject
eez_sf <- st_read(eez_shp_path, quiet = TRUE) %>% 
  st_transform(4326) %>%                        # Force WGS84 first
  st_make_valid() %>% 
  st_break_antimeridian(lon_0 = 0) %>%          # Cut at antimeridian (prevents side-to-side artifacts)
  st_collection_extract("POLYGON") %>% 
  st_transform(PROJ) %>%                        # Reproject to Mollweide
  st_make_valid()

# 2. Load High Seas with the same treatment
high_seas_sf <- st_read(high_seas_shp_path, quiet = TRUE) %>% 
  st_transform(4326) %>% 
  st_make_valid() %>% 
  st_break_antimeridian(lon_0 = 0) %>% 
  st_collection_extract("POLYGON") %>% 
  st_transform(PROJ) %>% 
  st_make_valid()

stopifnot(EEZ_COUNTRY_FIELD %in% names(eez_sf))

# 3. Prepare overlap points
overlap_pts <- overlap_sf %>% 
  st_make_valid() %>% 
  mutate(cell_id = row_number())

# 4. Spatial intersections
in_eez <- st_join(overlap_pts, eez_sf[, EEZ_COUNTRY_FIELD], join = st_intersects, left = TRUE) %>%
  st_drop_geometry() %>%
  select(cell_id, eez_country = all_of(EEZ_COUNTRY_FIELD)) %>%
  distinct(cell_id, .keep_all = TRUE)

in_high_seas <- st_join(overlap_pts, high_seas_sf, join = st_intersects, left = TRUE) %>%
  st_drop_geometry() %>%
  distinct(cell_id, .keep_all = TRUE)

in_high_seas$is_high_seas <- !is.na(in_high_seas[[names(high_seas_sf)[1]]])

# 5. Combined classification
jurisdiction <- overlap_pts %>%
  st_drop_geometry() %>%
  select(cell_id, x, y) %>%
  left_join(in_eez, by = "cell_id") %>%
  left_join(in_high_seas %>% select(cell_id, is_high_seas), by = "cell_id") %>%
  mutate(
    jurisdiction = case_when(
      !is.na(eez_country) ~ "EEZ",
      is_high_seas        ~ "High Seas",
      TRUE                ~ NA_character_
    )
  )

# 6. Assign residual boundary/coastal cells by minimum distance
n_unclassified <- sum(is.na(jurisdiction$jurisdiction))

if (n_unclassified > 0) {
  message(sprintf("%d unclassified cells found on boundaries. Resolving by nearest feature...", n_unclassified))
  
  combined_sf <- bind_rows(
    eez_sf %>% transmute(zone_type = "EEZ", eez_country_nn = .data[[EEZ_COUNTRY_FIELD]]),
    high_seas_sf %>% transmute(zone_type = "High Seas", eez_country_nn = NA_character_)
  )
  
  orphan_idx <- which(is.na(jurisdiction$jurisdiction))
  orphan_pts <- st_as_sf(
    jurisdiction[orphan_idx, c("cell_id", "x", "y")],
    coords = c("x", "y"), crs = PROJ, remove = FALSE
  )
  
  nearest_idx <- st_nearest_feature(orphan_pts, combined_sf)
  
  jurisdiction$jurisdiction[orphan_idx] <- combined_sf$zone_type[nearest_idx]
  jurisdiction$eez_country[orphan_idx] <- ifelse(
    combined_sf$zone_type[nearest_idx] == "EEZ",
    combined_sf$eez_country_nn[nearest_idx],
    jurisdiction$eez_country[orphan_idx]
  )
}

sf_use_s2(TRUE)

# Join overlap geometry with the classified jurisdiction table
overlap_classified <- overlap_pts %>%
  left_join(jurisdiction, by = c("cell_id", "x", "y"))

ggplot() +
  # Draw EEZ boundaries in light gray background
  geom_sf(data = eez_sf, fill = "grey90", color = "white", size = 0.2) +
  # Draw classified cells
  geom_sf(data = overlap_classified, aes(color = jurisdiction), size = 0.6, alpha = 0.8) +
  scale_color_manual(
    values = c("EEZ" = "#0A3B5C", "High Seas" = "#7A2021"),
    na.value = "yellow" # Highlight any unclassified cells in bright yellow
  ) +
  labs(
    title = "Jurisdictional Classification Diagnostics",
    subtitle = "Blue: EEZ | Wine: High Seas | Yellow: Unassigned/Orphans",
    color = "Jurisdiction"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

overlap_classified_path <- paste(overlap_dir, sprintf("overlap_classified_%s.rds", GEARTYPE), sep = "/")
saveRDS(overlap_classified, overlap_classified_path)
message("Saved: ", overlap_classified_path)
# ----------------------------------------------------
# 6b. % overlap area (undersampled + high fishing effort) by jurisdiction
# ----------------------------------------------------
message("Computing % of overlap area within EEZ vs High Seas...")

overlap_classified_path <- paste(overlap_dir, sprintf("overlap_classified_%s.rds", GEARTYPE), sep = "/")
if (!file.exists(overlap_classified_path)) {
  stop("overlap_classified file not found. Please run section 6 (jurisdiction classification) first!")
}
overlap_classified <- readRDS(overlap_classified_path)

overlap_area_by_jurisdiction <- overlap_classified %>%
  st_drop_geometry() %>%
  group_by(jurisdiction) %>%
  summarise(
    n_cells  = n(),
    area_km2 = sum(area_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_cells = 100 * n_cells / sum(n_cells),
    pct_area  = 100 * area_km2 / sum(area_km2)
  ) %>%
  arrange(desc(area_km2))

cat("\n========== OVERLAP AREA BY JURISDICTION (", GEARTYPE, ") ==========\n", sep = "")
print(overlap_area_by_jurisdiction)
cat(sprintf(
  "\n%.1f%% of the overlap area falls within EEZs and %.1f%% within High Seas.\n",
  overlap_area_by_jurisdiction$pct_area[overlap_area_by_jurisdiction$jurisdiction == "EEZ"],
  overlap_area_by_jurisdiction$pct_area[overlap_area_by_jurisdiction$jurisdiction == "High Seas"]
))
cat("=========================================================================\n")

overlap_jurisdiction_path <- paste(overlap_dir, sprintf("overlap_area_by_jurisdiction_%s.csv", GEARTYPE), sep = "/")
write.csv(overlap_area_by_jurisdiction, overlap_jurisdiction_path, row.names = FALSE)
message("Saved: ", overlap_jurisdiction_path)

# ----------------------------------------------------
# 7. Combine effort by flag + jurisdiction, and aggregate
# ----------------------------------------------------
flag_jurisdiction <- flag_cell_values %>%
  left_join(jurisdiction %>% select(x, y, jurisdiction), by = c("x", "y"))

flag_summary <- flag_jurisdiction %>%
  group_by(flag, jurisdiction) %>%
  summarise(total_hours = sum(mean_yearly_hours, na.rm = TRUE), .groups = "drop")

flag_summary_path <- paste(overlap_dir, sprintf("flag_eez_summary_raw_%s.csv", GEARTYPE), sep = "/")
write.csv(flag_summary, flag_summary_path, row.names = FALSE)

# ----------------------------------------------------
# 8. Barplot: Top N Flag States with Custom Hex Gradients
# ----------------------------------------------------


rank_top_n <- function(df_jurisdiction, top_n = TOP_N_FLAGS) {
  df_jurisdiction %>%
    arrange(desc(total_hours)) %>%
    slice_head(n = top_n)
}

eez_ranked <- flag_summary %>% filter(jurisdiction == "EEZ")         %>% rank_top_n()
hs_ranked  <- flag_summary %>% filter(jurisdiction == "High Seas") %>% rank_top_n()

order_for_plot <- function(df) {
  df %>% mutate(flag_group = factor(flag, levels = flag[order(total_hours)]))
}

eez_ranked <- order_for_plot(eez_ranked)
hs_ranked  <- order_for_plot(hs_ranked)

n_flags_eez <- flag_summary %>% filter(jurisdiction == "EEZ") %>% nrow()
n_flags_hs  <- flag_summary %>% filter(jurisdiction == "High Seas") %>% nrow()

make_ranked_plot <- function(df, title_label, low_color, high_color, n_total_flags) {
  panel_max <- max(df$total_hours) * 1.1
  
  ggplot(df, aes(x = flag_group, y = total_hours)) +
    geom_col(aes(fill = total_hours), show.legend = FALSE, width = 0.7) +
    scale_fill_gradient(low = low_color, high = high_color) +
    coord_flip(clip = "off") +
    scale_y_continuous(
      limits = c(0, panel_max),
      labels = function(x) round(x / 1e5, 0),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = NULL, 
      y = expression("Mean annual fishing effort" ~ (10^5 ~ h)),
      title = title_label,
      subtitle = sprintf("Top %d of %d active flag states", nrow(df), n_total_flags)
    ) +
    theme_minimal(base_size = 14) +  # Tamaño base general incrementado
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(size = 20, face = "bold"),   # Título más grande
      plot.subtitle = element_text(size = 14, color = "grey30"), # Subtítulo más legible
      axis.text.y = element_text(size = 15, face = "bold"),  # Etiquetas de los países más grandes
      axis.text.x = element_text(size = 15),                 # Números del eje X más grandes
      axis.title.x = element_text(size = 15, face = "bold"), # Título del eje X más grande
      plot.margin = margin(5, 15, 5, 5)
    )
}

# Panel 1: EEZ (Petroleum Blue) and High Seas (Terracotta/Wine)
p_eez <- make_ranked_plot(
  eez_ranked, "EEZ", 
  low_color = "#8DA9C4", high_color = "#0A3B5C", 
  n_flags_eez
)

p_hs <- make_ranked_plot(
  hs_ranked, "High Seas", 
  low_color = "#FFA07A", high_color = "#8B5F65", 
  n_flags_hs
)

# Combine both panels
p_combined <- p_eez + p_hs

# Save image to disk
barplot_path <- paste(overlap_dir, sprintf("6_barplot_flag_by_jurisdiction_%s.png", GEARTYPE), sep = "/")
ggsave(barplot_path, plot = p_combined, device = "png", dpi = 400, width = 30, height = 11, units = "cm")

message("Combined plot successfully saved to: ", barplot_path)
