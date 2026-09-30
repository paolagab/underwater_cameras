#-----------------------------------------------------------------------------
# 01a_overlap.R
# Description: Computes the spatial overlap between undersampled OBIS areas
#              (03_undersampled_areas.R) and high fishing effort GFW areas
#              (02_high_fe_areas.R). Both rasters share the same
#              1-degree Mollweide grid, so the overlay is a simple cell-by-
#              cell combination, no reprojection/interpolation needed.
#-----------------------------------------------------------------------------


source("R/utils.R")
source("R/data_paths.R")

library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# ----------------------------------------------------
# 0. Parameters
# ----------------------------------------------------
GEARTYPE <- "drifting_longlines"

# ----------------------------------------------------
# 1. Load already-computed classifications (do not recompute)
# ----------------------------------------------------
ocean_mask   <- rast(temp_mask)
eff_area_km2 <- rast(eff_area_km2_nc)
PROJ <- crs(ocean_mask)

undersamp_stack <- rast(paste(obis_dir, "obis_undersampled_areas.tif", sep = "/"))
undersamp_class <- undersamp_stack[["sampling_class"]]
# 03_undersampled_areas.R lookup: 1 = never_sampled, 2 = low, 3 = medium, 4 = high
is_undersampled <- undersamp_class %in% c(1, 2)
names(is_undersampled) <- "is_undersampled"

effort_stack <- rast(paste(fe_dir, sprintf("gfw_high_effort_areas_%s.tif", GEARTYPE), sep = "/"))
effort_class <- effort_stack[["effort_class"]]
# 03b lookup: 1 = no_activity, 2 = low, 3 = medium, 4 = high
is_high_effort <- effort_class == 4
is_any_activity <- effort_class != 1  # "AOO" -- any recorded effort, regardless of intensity
names(is_high_effort)  <- "is_high_effort"
names(is_any_activity) <- "is_any_activity"

# ----------------------------------------------------
# 2. Cell-by-cell overlap
# ----------------------------------------------------
message("Computing undersampled x high-effort overlap...")

overlap <- is_undersampled & is_high_effort
names(overlap) <- "overlap"
overlap <- mask(overlap, ocean_mask)

# Sensitivity variant computed here too, so it can go into the same
# data.frame build below instead of being extracted separately later.
is_medium_or_high <- effort_class %in% c(3, 4)
overlap_relaxed <- is_undersampled & is_medium_or_high
names(overlap_relaxed) <- "overlap_relaxed"
overlap_relaxed <- mask(overlap_relaxed, ocean_mask)

# ----------------------------------------------------
# 3. Single combined data.frame (feeds stats, basin breakdown, and the map)
# ----------------------------------------------------
df <- as.data.frame(
  c(is_undersampled, is_high_effort, is_any_activity, overlap, overlap_relaxed, eff_area_km2),
  xy = TRUE, na.rm = TRUE
)
names(df) <- c("x", "y", "is_undersampled", "is_high_effort", "is_any_activity",
               "overlap", "overlap_relaxed", "area_km2")

total_ocean_km2     <- sum(df$area_km2, na.rm = TRUE)
undersampled_km2    <- sum(df$area_km2[df$is_undersampled], na.rm = TRUE)
high_effort_km2     <- sum(df$area_km2[df$is_high_effort], na.rm = TRUE)
overlap_km2         <- sum(df$area_km2[df$overlap], na.rm = TRUE)

cat("\n========== UNDERSAMPLED (OBIS) x HIGH EFFORT (GFW - ", GEARTYPE, ") ==========\n", sep = "")
cat(sprintf("Total ocean area:        %12.0f km2\n", total_ocean_km2))
cat(sprintf("Undersampled area:       %12.0f km2  (%.1f%% of ocean)\n",
            undersampled_km2, 100 * undersampled_km2 / total_ocean_km2))
cat(sprintf("High effort area (%s): %12.0f km2  (%.1f%% of ocean)\n",
            GEARTYPE, high_effort_km2, 100 * high_effort_km2 / total_ocean_km2))
cat(sprintf("OVERLAP area:            %12.0f km2  (%.1f%% of ocean | %.1f%% of undersampled area | %.1f%% of high-effort area)\n",
            overlap_km2, 100 * overlap_km2 / total_ocean_km2,
            100 * overlap_km2 / undersampled_km2, 100 * overlap_km2 / high_effort_km2))
cat("======================================================================================\n")

overlap_stats <- data.frame(
  metric = c("total_ocean_km2", "undersampled_km2", "high_effort_km2", "overlap_km2",
             "pct_ocean_undersampled", "pct_ocean_high_effort", "pct_ocean_overlap",
             "pct_undersampled_that_is_high_effort", "pct_high_effort_that_is_undersampled"),
  value = c(total_ocean_km2, undersampled_km2, high_effort_km2, overlap_km2,
            100 * undersampled_km2 / total_ocean_km2,
            100 * high_effort_km2 / total_ocean_km2,
            100 * overlap_km2 / total_ocean_km2,
            100 * overlap_km2 / undersampled_km2,
            100 * overlap_km2 / high_effort_km2)
)

overlap_stats_path <- paste(overlap_dir, sprintf("overlap_stats_%s.csv", GEARTYPE), sep = "/")
write.csv(overlap_stats, overlap_stats_path, row.names = FALSE)

# ----------------------------------------------------
# 4. Sensitivity check (DESCRIPTIVE): overlap if "high effort" were relaxed
#    to include the medium tercile too (medium+high vs high-only). Mirrors
#    the earlier p75-vs-p90 sensitivity check, now expressed in terms of
#    the tercile classification instead of an obsolete percentile pair.
# ----------------------------------------------------
overlap_relaxed_km2 <- sum(df$area_km2[df$overlap_relaxed], na.rm = TRUE)

cat(sprintf(
  "\n[DESCRIPTIVE] Relaxing the effort criterion to medium+high tercile (instead of high only)\nwould raise the overlap to %.0f km2 (%.1f%% of ocean), vs %.0f km2 (%.1f%%) with the primary\n(high-only) criterion. Sensitivity check, does not replace the primary result.\n",
  overlap_relaxed_km2, 100 * overlap_relaxed_km2 / total_ocean_km2,
  overlap_km2, 100 * overlap_km2 / total_ocean_km2
))

overlap_stats <- rbind(
  overlap_stats,
  data.frame(metric = c("overlap_medium_or_high_km2_DESCRIPTIVE", "pct_ocean_overlap_medium_or_high_DESCRIPTIVE"),
             value = c(overlap_relaxed_km2, 100 * overlap_relaxed_km2 / total_ocean_km2))
)

# ----------------------------------------------------
# 5. AOO comparison (DESCRIPTIVE): overlap using ANY recorded effort
#    (Area of Occupancy, IUCN Criterion B logic) instead of the high-effort
#    criterion. AOO counts only cells with real activity, no intensity
#    threshold and no assumption of spatial continuity.
# ----------------------------------------------------
AOO_area_km2 <- sum(df$area_km2[df$is_any_activity], na.rm = TRUE)
S_aoo_km2    <- sum(df$area_km2[df$is_undersampled & df$is_any_activity], na.rm = TRUE)
overlap_aoo_pct <- 100 * S_aoo_km2 / undersampled_km2

cat(sprintf(
  "\n[DESCRIPTIVE] AOO (area of occupancy/footprint) for %s: %.0f km2 (%.1f%% of total ocean).\nUsing AOO instead of the high-effort criterion, %.1f%% of the undersampled area\n(%.0f km2 of %.0f km2) falls within the real footprint of fishing activity --\nvs %.1f%% with the primary (high-effort) criterion.\n",
  GEARTYPE, AOO_area_km2, 100 * AOO_area_km2 / total_ocean_km2,
  overlap_aoo_pct, S_aoo_km2, undersampled_km2,
  100 * overlap_km2 / undersampled_km2
))

overlap_stats <- rbind(
  overlap_stats,
  data.frame(metric = c("AOO_area_km2_DESCRIPTIVE", "S_aoo_km2_DESCRIPTIVE", "pct_undersampled_within_AOO_DESCRIPTIVE"),
             value = c(AOO_area_km2, S_aoo_km2, overlap_aoo_pct))
)

write.csv(overlap_stats, overlap_stats_path, row.names = FALSE)

# ----------------------------------------------------
# 6. Basin-level breakdown (using precomputed basin_lookup from study_area.R)
# ----------------------------------------------------
message("Loading basin lookup table and merging with overlap data...")

# 1. Load the lookup table generated in study_area.R
lookup_path <- paste(temp_dir, "cell_basin_lookup.rds", sep = "/")
if (!file.exists(lookup_path)) {
  stop("basin_lookup file not found. Please run study_area.R first!")
}
basin_lookup <- readRDS(lookup_path)

# 2. Assign cell_id to your current dataframe
df$cell_id <- terra::cellFromXY(ocean_mask, as.matrix(df[, c("x", "y")]))

# 3. Merge via instant left_join using cell_id (no heavy spatial intersections needed)
df_basin <- df %>%
  left_join(basin_lookup %>% select(cell_id, ocean_basin), by = "cell_id")

# Check if any cell missed its basin assignment
n_missing <- sum(is.na(df_basin$ocean_basin))
if (n_missing > 0) {
  warning(sprintf("%d cells missing ocean basin assignment. Assigning to 'Unknown'.", n_missing))
  df_basin$ocean_basin[is.na(df_basin$ocean_basin)] <- "Unknown"
}

# 4. Basin-level overlap metrics
basin_overlap <- df_basin %>%
  group_by(ocean_basin) %>%
  summarise(
    basin_area_km2             = sum(area_km2, na.rm = TRUE),
    U_undersampled_km2         = sum(area_km2[is_undersampled], na.rm = TRUE),
    S_overlap_km2              = sum(area_km2[overlap], na.rm = TRUE),
    S_overlap_relaxed_km2      = sum(area_km2[overlap_relaxed], na.rm = TRUE),
    AOO_area_km2_basin         = sum(area_km2[is_any_activity], na.rm = TRUE),
    S_aoo_km2_basin            = sum(area_km2[is_undersampled & is_any_activity], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_ocean_km2              = 100 * basin_area_km2 / sum(basin_area_km2),
    pct_basin_undersampled     = ifelse(basin_area_km2 > 0, U_undersampled_km2 / basin_area_km2 * 100, NA),
    overlap_index_pct          = ifelse(U_undersampled_km2 > 0, S_overlap_km2 / U_undersampled_km2 * 100, NA),
    overlap_index_relaxed_pct  = ifelse(U_undersampled_km2 > 0, S_overlap_relaxed_km2 / U_undersampled_km2 * 100, NA),
    overlap_index_aoo_pct      = ifelse(U_undersampled_km2 > 0, S_aoo_km2_basin / U_undersampled_km2 * 100, NA)
  ) %>%
  arrange(desc(U_undersampled_km2))

cat("\n========== OVERLAP INDEX BY BASIN -", GEARTYPE, "==========\n")
print(basin_overlap)
cat("================================================================\n")

overlap_basin_path <- paste(overlap_dir, sprintf("overlap_index_by_basin_%s.csv", GEARTYPE), sep = "/")
write.csv(basin_overlap, overlap_basin_path, row.names = FALSE)
message("Saved: ", overlap_basin_path)
# ----------------------------------------------------
# 7. Save overlap cells (input for the next step: flag / EEZ / high seas)
# ----------------------------------------------------
overlap_sf <- df %>%
  filter(overlap) %>%
  st_as_sf(coords = c("x", "y"), crs = PROJ, remove = FALSE)

overlap_cells_path <- paste(overlap_dir, sprintf("overlap_cells_%s.rds", GEARTYPE), sep = "/")
saveRDS(overlap_sf, overlap_cells_path)
message("Saved: ", overlap_cells_path)

# ----------------------------------------------------
# 8. Category map (overlap / undersampled only / high effort only / neither)
# ----------------------------------------------------
df <- df %>%
  mutate(
    category = case_when(
      overlap                                  ~ "overlap",
      is_undersampled & !overlap               ~ "undersampled_only",
      is_high_effort & !overlap                ~ "high_effort_only",
      TRUE                                      ~ "neither"
    )
  )

land <- st_read(ne_50m_shp, quiet = TRUE)
land_moll <- st_transform(land, crs = PROJ)

bbox <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = PROJ)

  

cat_colors <- c(
  overlap            = "#D62839",
  undersampled_only  = "#3A86FF",
  high_effort_only   = "#F4A261",
  neither            = "grey88"
)
cat_labels <- c(
  overlap            = "Undersampled + high fishing effort",
  undersampled_only  = "Undersampled only",
  high_effort_only   = sprintf("High fishing effort only (%s)", GEARTYPE),
  neither            = "Neither"
)

map_overlap <- ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = category)) +
  scale_fill_manual(values = cat_colors, labels = cat_labels, name = NULL) +
  geom_sf(data = land_moll, fill = "grey30", color = NA) +
  geom_sf(data = bbox, fill = NA, color = "grey40", linewidth = 0.5) +
  coord_sf(datum = NA) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "right",
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    legend.text      = element_text(size = 8)
  ) +
  labs(title = "", x = "", y = "")

map_overlap_path <- paste(overlap_dir, sprintf("map_overlap_%s.png", GEARTYPE), sep = "/")
ggsave(map_overlap_path, plot = map_overlap, device = "png", dpi = 400, width = 24, height = 14, units = "cm")

message("Map saved: ", map_overlap_path)
