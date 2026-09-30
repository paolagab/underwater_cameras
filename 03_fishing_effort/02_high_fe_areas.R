#-----------------------------------------------------------------------------
# 02_high_fe_areas.R
# Description: Classifies GFW effort density into four categories --
#              no_activity, low, medium, high -- using terciles of effort
#              density among active cells (mirrors 02b_undersampled_areas.R
#              for OBIS). High effort area = high tercile only.
#
#              Reports area globally and by ocean basin, and -- as
#              DESCRIPTIVE context only, not a classification criterion --
#              temporal persistence (years active) within high-effort cells.
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
GEARTYPE <- "drifting_longlines"
THRESHOLD_METHOD <- "tercile"  # high = p66, computed on active cells only
years <- 2016:2024

# ----------------------------------------------------
# 1. Load already-computed density (do not recompute)
# ----------------------------------------------------
ocean_mask     <- rast(temp_mask)
eff_area_km2   <- rast(eff_area_km2_nc)
effort_density <- rast(paste(fe_dir, sprintf("gfw_effort_density_%s.tif", GEARTYPE), sep = "/"))

# ----------------------------------------------------
# 2. Tercile thresholds (active cells only, i.e. density > 0)
# ----------------------------------------------------
effort_vals <- values(effort_density, na.rm = TRUE)[, 1]
effort_vals_active <- effort_vals[effort_vals > 0]
log_effort_active  <- log10(effort_vals_active)

thr_low_log  <- quantile(log_effort_active, 1/3, na.rm = TRUE)
thr_high_log <- quantile(log_effort_active, 2/3, na.rm = TRUE)
thr_low  <- 10^thr_low_log
thr_high <- 10^thr_high_log

message(sprintf(
  "Tercile thresholds (%s): low <= %.6f | high > %.6f hours/km2 (%.1f%% / %.1f%% of active cells)",
  GEARTYPE, thr_low, thr_high,
  100 * mean(log_effort_active <= thr_low_log),
  100 * mean(log_effort_active > thr_high_log)
))

# ----------------------------------------------------
# 3. Classification
# ----------------------------------------------------
message("Classifying cells into no_activity / low / medium / high...")

df <- as.data.frame(c(effort_density, eff_area_km2), xy = TRUE, na.rm = TRUE)
names(df) <- c("x", "y", "effort_density", "area_km2")

df <- df %>%
  mutate(
    effort_class = case_when(
      effort_density == 0           ~ "no_activity",
      effort_density <= thr_low     ~ "low",
      effort_density <= thr_high    ~ "medium",
      TRUE                          ~ "high"
    ),
    is_high_effort = effort_class == "high"
  )

cat("\n========== EFFORT CLASSIFICATION ==========\n")
print(table(df$effort_class))
cat("==============================================\n")

# ----------------------------------------------------
# 4. Class raster
# ----------------------------------------------------
class_lookup <- c(no_activity = 1L, low = 2L, medium = 3L, high = 4L)
df$class_int <- class_lookup[df$effort_class]

class_rast <- rast(ocean_mask)
values(class_rast) <- NA
idx <- cellFromXY(class_rast, as.matrix(df[, c("x", "y")]))
class_rast[idx] <- df$class_int
names(class_rast) <- "effort_class"

high_effort_metrics_tif <- paste(fe_dir, sprintf("gfw_high_effort_areas_%s.tif", GEARTYPE), sep = "/")
writeRaster(c(effort_density, class_rast), high_effort_metrics_tif, overwrite = TRUE)
message("Raster saved to: ", high_effort_metrics_tif)

# ----------------------------------------------------
# 5. SF object (for spatial analysis)
# ----------------------------------------------------
PROJ <- crs(ocean_mask)
effort_sf <- st_as_sf(df, coords = c("x", "y"), crs = PROJ, remove = FALSE)
effort_sf_rds <- paste(fe_dir, sprintf("effort_sf_%s.rds", GEARTYPE), sep = "/")
saveRDS(effort_sf, effort_sf_rds)
message("SF object saved to: ", effort_sf_rds)

# ----------------------------------------------------
# 6. Global statistics
# ----------------------------------------------------
total_ocean_km2 <- sum(df$area_km2, na.rm = TRUE)

global_stats <- df %>%
  group_by(effort_class) %>%
  summarise(n_cells = n(), area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct_ocean = area_km2 / total_ocean_km2 * 100)

high_effort_km2 <- sum(global_stats$area_km2[global_stats$effort_class == "high"])

cat("\n========== GLOBAL STATISTICS ==========\n")
cat(sprintf("Total ocean area:    %10.0f km2\n", total_ocean_km2))
cat(sprintf("High effort area:    %10.0f km2  (%.1f%% of ocean)\n",
            high_effort_km2, high_effort_km2 / total_ocean_km2 * 100))
cat("\nBreakdown by category:\n")
print(global_stats)
cat("========================================\n")

write.csv(global_stats, effort_global_stats_csv, row.names = FALSE)
# ----------------------------------------------------
# 7. Basin-level breakdown (using precomputed basin_lookup from study_area.R)
# ----------------------------------------------------
message("Loading basin lookup table and merging with effort data...")

# 1. Load the lookup table generated in study_area.R
lookup_path <- paste(temp_dir, "cell_basin_lookup.rds", sep = "/")
if (!file.exists(lookup_path)) {
  stop("basin_lookup file not found. Please run study_area.R first!")
}
basin_lookup <- readRDS(lookup_path)

# 2. Assign cell_id to your current active cells data frame
df$cell_id <- terra::cellFromXY(ocean_mask, as.matrix(df[, c("x", "y")]))

# 3. Merge via instant left_join (no heavy spatial reprojections or intersections needed)
df_basin <- df %>%
  left_join(basin_lookup %>% select(cell_id, ocean_basin), by = "cell_id")

# Check if any cell missed its basin assignment
n_missing <- sum(is.na(df_basin$ocean_basin))
if (n_missing > 0) {
  warning(sprintf("%d cells missing ocean basin assignment. Assigning to 'Unknown'.", n_missing))
  df_basin$ocean_basin[is.na(df_basin$ocean_basin)] <- "Unknown"
}

# 4. Basin-level statistics breakdown by effort class and ocean basin
basin_stats <- df_basin %>%
  group_by(ocean_basin, effort_class) %>%
  summarise(
    n_cells = n(),
    area_km2 = sum(area_km2, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n========== BASIN-LEVEL STATISTICS ==========\n")
print(basin_stats)
cat("==============================================\n")

# Save the basin-level statistics table
write.csv(basin_stats, paste(fe_dir, sprintf("gfw_basin_stats_%s.csv", GEARTYPE), sep = "/"), row.names = FALSE)
# ----------------------------------------------------
# 8. Map (no_activity / low+medium ("AOO") / high highlighted)
# ----------------------------------------------------
message("Generating high effort areas map...")

land <- st_read(ne_50m_shp, quiet = TRUE)
land_moll <- st_transform(land, crs = PROJ)

bbox <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = PROJ)

df <- df %>%
  mutate(map_class = case_when(
    effort_class == "no_activity" ~ "no_activity",
    effort_class == "high"        ~ "high",
    TRUE                          ~ "aoo"   # low + medium combined ("some activity")
  ))

cat_colors <- c(no_activity = "transparent", aoo = "#F2C79E", high = "#B5540C")
cat_labels <- c(no_activity = "No recorded activity", aoo = "Active (AOO)",
                high = "High effort (upper tercile)")

map_effort <- ggplot() +
  geom_tile(data = df, aes(x = x, y = y, fill = map_class)) +
  scale_fill_manual(values = cat_colors, labels = cat_labels, name = "Fishing effort") +
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
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  ) +
  labs(title = "", x = "", y = "")

ggsave(map_high_effort_areas_png, plot = map_effort, device = "png", dpi = 400,
       width = 24, height = 14, units = "cm")

# ----------------------------------------------------
# 9. Temporal persistence within high-effort cells (DESCRIPTIVE ONLY)
# ----------------------------------------------------
message("Computing descriptive persistence within high-effort cells...")

effort_raw <- readRDS(gfw_effort_raw_rds) %>% filter(geartype == GEARTYPE)

to_raster_year <- function(df_year, mask_rast) {
  if (nrow(df_year) == 0) return(NULL)
  v <- vect(df_year, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  v <- project(v, crs(mask_rast))
  r <- rasterize(v, mask_rast, field = "yearly_hours", fun = "sum", background = 0)
  mask(r, mask_rast)
}

year_rasters <- list()
for (yr in years) {
  df_yr <- effort_raw %>%
    filter(year == yr) %>%
    group_by(Lon, Lat) %>%
    summarise(yearly_hours = sum(`Apparent Fishing Hours`, na.rm = TRUE), .groups = "drop")
  r_yr <- to_raster_year(df_yr, ocean_mask)
  if (!is.null(r_yr)) {
    names(r_yr) <- paste0("Y", yr)
    year_rasters[[as.character(yr)]] <- r_yr
  }
}

if (length(year_rasters) != length(years)) {
  warning(sprintf("Expected %d years of data, got %d. Check if some year had no activity at all for this gear.",
                  length(years), length(year_rasters)))
}

effort_year_stack <- rast(year_rasters)
years_active <- app(effort_year_stack, fun = function(x) sum(x > 0, na.rm = TRUE))
names(years_active) <- "years_active"
years_active <- mask(years_active, ocean_mask)

df_persist <- as.data.frame(c(class_rast, years_active), xy = TRUE, na.rm = TRUE)
names(df_persist) <- c("x", "y", "class_int", "years_active")
df_high <- df_persist %>% filter(class_int == class_lookup["high"])

cat("\n========== TEMPORAL PERSISTENCE WITHIN HIGH-EFFORT CELLS ==========\n")
cat(sprintf("N high-effort cells: %d\n", nrow(df_high)))
print(quantile(df_high$years_active, probs = c(0.10, 0.25, 0.50, 0.75, 0.90)))
for (k in c(5, 7, 8, 9)) {
  cat(sprintf("%% of high-effort cells active in >= %d of %d years: %.1f%%\n",
              k, length(years), 100 * mean(df_high$years_active >= k)))
}
cat("=======================================================================\n")

message("High effort areas classification complete.")
