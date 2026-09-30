#-----------------------------------------------------------------------------
# 03_undersampled_areas.R
# Description: Classifies OBIS sampling into four categories -- never_sampled,
#              low, medium, high -- using terciles of sampling density among
#              sampled cells (following Mestre et al. 2025, Biol. Conserv.,
#              who define "shipless areas" via a first-tercile cutoff on
#              vessel density). Undersampled area = never_sampled + low.
#
#              Unlike Mestre et al., terciles here are computed ONLY on
#              cells with sampling_density > 0 (never_sampled is kept as a
#              separate, more informative category rather than folded into
#              the tercile calculation).
#
#              Reports area globally and by ocean basin.
#-----------------------------------------------------------------------------

source("R/utils.R")
source("R/data_paths.R")

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

# ----------------------------------------------------
# 0. Parameters
# ----------------------------------------------------
DENSITY_METHOD <- "tercile"  # low = p33, high = p66, computed on sampled cells only (see 02a)

# ----------------------------------------------------
# 1. Load already-computed sampling density (from 02a, do not recompute)
# ----------------------------------------------------
sampling_density_tif <- paste(obis_dir, "obis_sampling_density.tif", sep = "/")
density_stack <- rast(sampling_density_tif)

total_events     <- density_stack[["total_events"]]
sampling_density <- density_stack[["sampling_density"]]
years_sampled_recent <- density_stack[["years_sampled_recent"]]

ocean_mask   <- rast(temp_mask)
eff_area_km2 <- rast(eff_area_km2_nc)

# ----------------------------------------------------
# 2. Tercile thresholds (sampled cells only, i.e. density > 0)
# ----------------------------------------------------
density_vals <- values(sampling_density, na.rm = TRUE)[, 1]
density_vals_sampled <- density_vals[density_vals > 0]
log_density_sampled  <- log10(density_vals_sampled)

thr_low_log  <- quantile(log_density_sampled, 1/3, na.rm = TRUE)
thr_high_log <- quantile(log_density_sampled, 2/3, na.rm = TRUE)
thr_low  <- 10^thr_low_log
thr_high <- 10^thr_high_log

message(sprintf(
  "Tercile thresholds: low <= %.6f events/km2 | high > %.6f events/km2 (%.1f%% / %.1f%% of sampled cells)",
  thr_low, thr_high,
  100 * mean(log_density_sampled <= thr_low_log),
  100 * mean(log_density_sampled > thr_high_log)
))

# ----------------------------------------------------
# 3. Classification
# ----------------------------------------------------
message("Classifying cells into never_sampled / low / medium / high...")

df <- as.data.frame(
  c(total_events, sampling_density, years_sampled_recent, eff_area_km2),
  xy = TRUE, na.rm = TRUE
)
names(df) <- c("x", "y", "total_events", "sampling_density", "years_sampled_recent", "area_km2")

df <- df %>%
  mutate(
    sampling_class = case_when(
      total_events == 0             ~ "never_sampled",
      sampling_density <= thr_low   ~ "low",
      sampling_density <= thr_high  ~ "medium",
      TRUE                          ~ "high"
    ),
    is_undersampled = sampling_class %in% c("never_sampled", "low")
  )

cat("\n========== SAMPLING CLASSIFICATION ==========\n")
print(table(df$sampling_class))
cat("================================================\n")

# ----------------------------------------------------
# 4. Class raster
# ----------------------------------------------------
class_lookup <- c(never_sampled = 1L, low = 2L, medium = 3L, high = 4L)
df$class_int <- class_lookup[df$sampling_class]

class_rast <- rast(ocean_mask)
values(class_rast) <- NA
idx <- cellFromXY(class_rast, as.matrix(df[, c("x", "y")]))
class_rast[idx] <- df$class_int
names(class_rast) <- "sampling_class"

undersampled_metrics_tif <- paste(obis_dir, "obis_undersampled_areas.tif", sep = "/")
writeRaster(c(sampling_density, class_rast), undersampled_metrics_tif, overwrite = TRUE)
message("Raster saved to: ", undersampled_metrics_tif)

# ----------------------------------------------------
# 5. SF object (for basin-level stats and later overlap analysis)
# ----------------------------------------------------
PROJ <- crs(ocean_mask)
undersamp_sf <- st_as_sf(df, coords = c("x", "y"), crs = PROJ, remove = FALSE)
undersamp_sf_rds <- paste(obis_dir, "undersamp_sf.rds", sep = "/")
saveRDS(undersamp_sf, undersamp_sf_rds)
message("SF object saved to: ", undersamp_sf_rds)

# ----------------------------------------------------
# 6. Global statistics
# ----------------------------------------------------
total_ocean_km2 <- sum(df$area_km2, na.rm = TRUE)

global_stats <- df %>%
  group_by(sampling_class) %>%
  summarise(n_cells = n(), area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct_ocean = area_km2 / total_ocean_km2 * 100)

undersampled_km2 <- sum(global_stats$area_km2[global_stats$sampling_class %in% c("never_sampled", "low")])

cat("\n========== GLOBAL STATISTICS ==========\n")
cat(sprintf("Total ocean area:      %10.0f km2\n", total_ocean_km2))
cat(sprintf("Undersampled area:     %10.0f km2  (%.1f%% of ocean)\n",
            undersampled_km2, undersampled_km2 / total_ocean_km2 * 100))
cat("\nBreakdown by category:\n")
print(global_stats)
cat("========================================\n")

global_stats_csv <- paste(obis_dir, "sampling_global_stats.csv", sep = "/")
write.csv(global_stats, global_stats_csv, row.names = FALSE)

# ----------------------------------------------------
# 7. Basin-level breakdown (using precomputed basin_lookup from study_area.R)
# ----------------------------------------------------
message("Loading basin lookup table and merging with sampling data...")

# 1. Load the lookup table generated in study_area.R
lookup_path <- paste(temp_dir, "cell_basin_lookup.rds", sep = "/")
if (!file.exists(lookup_path)) {
  stop("basin_lookup file not found. Please run study_area.R first!")
}
basin_lookup <- readRDS(lookup_path)

# 2. Assign cell_id to your current sampling cells dataframe
df$cell_id <- terra::cellFromXY(ocean_mask, as.matrix(df[, c("x", "y")]))

# 3. Merge via instant left_join using cell_id
undersamp_basin <- df %>%
  left_join(basin_lookup %>% select(cell_id, ocean_basin), by = "cell_id")

# Check if any cell missed its basin assignment
n_missing <- sum(is.na(undersamp_basin$ocean_basin))
if (n_missing > 0) {
  warning(sprintf("%d cells missing ocean basin assignment. Assigning to 'Unknown'.", n_missing))
  undersamp_basin$ocean_basin[is.na(undersamp_basin$ocean_basin)] <- "Unknown"
}

# 4. Detailed breakdown statistics by ocean basin and sampling class
basin_stats <- undersamp_basin %>%
  group_by(ocean_basin, sampling_class) %>%
  summarise(n_cells = n(), area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") %>%
  group_by(ocean_basin) %>%
  mutate(basin_total_km2 = sum(area_km2), pct_of_basin = area_km2 / basin_total_km2 * 100) %>%
  ungroup() %>%
  arrange(ocean_basin, sampling_class)

cat("\n========== UNDERSAMPLING BY OCEAN BASIN ==========\n")
print(basin_stats, n = Inf)
cat("=====================================================\n")

basin_stats_csv <- paste(obis_dir, "sampling_stats_by_basin.csv", sep = "/")
write.csv(basin_stats, basin_stats_csv, row.names = FALSE)

# 5. Simplified summary table by ocean basin
basin_summary_simple <- basin_stats %>%
  select(ocean_basin, sampling_class, area_km2) %>%
  pivot_wider(names_from = sampling_class, values_from = area_km2, values_fill = 0) %>%
  mutate(
    across(any_of(c("never_sampled", "low", "medium", "high")), ~ replace_na(., 0)),
    basin_total_km2    = rowSums(pick(any_of(c("never_sampled", "low", "medium", "high"))), na.rm = TRUE),
    total_area_10k_km2 = basin_total_km2 / 1e4,
    pct_never_sampled  = if("never_sampled" %in% names(.)) never_sampled / basin_total_km2 * 100 else 0,
    pct_low            = if("low" %in% names(.)) low / basin_total_km2 * 100 else 0,
    pct_undersampled   = ((if("never_sampled" %in% names(.)) never_sampled else 0) + (if("low" %in% names(.)) low else 0)) / basin_total_km2 * 100,
    pct_medium         = if("medium" %in% names(.)) medium / basin_total_km2 * 100 else 0,
    pct_high           = if("high" %in% names(.)) high / basin_total_km2 * 100 else 0
  ) %>%
  select(ocean_basin, total_area_10k_km2, any_of(c("pct_never_sampled", "pct_low", "pct_undersampled", "pct_medium", "pct_high"))) %>%
  arrange(desc(total_area_10k_km2))

cat("\n========== SIMPLIFIED SUMMARY BY BASIN ==========\n")
print(basin_summary_simple)
cat("====================================================\n")

basin_summary_csv <- paste(obis_dir, "sampling_summary_simple_by_basin.csv", sep = "/")
write.csv(basin_summary_simple, basin_summary_csv, row.names = FALSE)

# ----------------------------------------------------
# 8. Map (never_sampled / low="undersampled" highlighted; medium+high = base ocean colour)
# ----------------------------------------------------
message("Generating undersampled areas map...")

land <- st_read(ne_shp, quiet = TRUE)
land_moll <- st_transform(land, crs = PROJ)

bbox <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = PROJ)

df <- df %>%
  mutate(map_class = ifelse(sampling_class %in% c("medium", "high"), "well_sampled", sampling_class))

cat_colors <- c(never_sampled = "#2C2C54", low = "#8DA9C4", well_sampled = "grey85")
cat_labels <- c(never_sampled = "Never sampled", low = "Undersampled (low tercile)",
                well_sampled = "Medium / high (not flagged)")

map_undersampled <- ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = map_class)) +
  scale_fill_manual(values = cat_colors, labels = cat_labels, name = "Sampling status") +
  geom_sf(data = land_moll, fill = "grey30", color = NA) +
  geom_sf(data = bbox, fill = NA, color = "grey40", linewidth = 0.5) +
  coord_sf(
    xlim = st_bbox(bbox)[c("xmin", "xmax")],
    ylim = st_bbox(bbox)[c("ymin", "ymax")],
    datum = NA, expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "right",
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  ) +
  labs(title = "", x = "", y = "")

map_png <- paste(obis_dir, "map_undersampled_areas.png", sep = "/")
ggsave(map_png, plot = map_undersampled, device = "png", dpi = 400, width = 24, height = 14, units = "cm")

message("Undersampled areas classification complete.")


