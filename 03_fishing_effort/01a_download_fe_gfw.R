#-----------------------------------------------------------------------------
# 01a_download_fe_gfw.R
# Description: Downloads GFW fishing effort and aggregates it to a raster
#              matching the OBIS grid. Focused on drifting longlines (the
#              study's main gear); purse seines optional via GEARTYPES.
#              Map generation lives in a separate script
#              (03d_global_effort_density_map.R).
#-----------------------------------------------------------------------------

source("R/data_paths.R")

library(gfwr)
library(purrr)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

key <- Sys.getenv("GFW_TOKEN")

# ----------------------------------------------------
# 0. Parameters
# ----------------------------------------------------
# Main gear: drifting_longlines. Add "purse_seines" to include seine data
# (e.g. for the FAD trial your supervisor mentioned).
GEARTYPES <- c("drifting_longlines")
# GEARTYPES <- c("drifting_longlines", "seiners")  # <- uncomment to include seine

years <- 2016:2024

# ----------------------------------------------------
# 1. Mask, effective area, projection
# ----------------------------------------------------
ocean_mask   <- rast(temp_mask)
eff_area_km2 <- rast(eff_area_km2_nc)

# ----------------------------------------------------
# 2. Global query polygon for the GFW API
# ----------------------------------------------------
world_polygon <- st_polygon(list(matrix(
  c(-180, -90,
    -180,  90,
    180,  90,
    180, -90,
    -180, -90),
  ncol = 2, byrow = TRUE
))) |>
  st_sfc(crs = 4326) |>
  st_sf(geometry = _)

# ----------------------------------------------------
# 3. Download, by year (all GEARTYPES defined above)
# ----------------------------------------------------
effort_raw <- map_dfr(years, function(y) {
  message("Downloading year: ", y)
  
  gfw_ais_fishing_hours(
    spatial_resolution  = "LOW",
    temporal_resolution = "YEARLY",
    group_by            = "GEARTYPE",
    start_date          = paste0(y, "-01-01"),
    end_date            = paste0(y, "-12-31"),
    region_source       = "USER_SHAPEFILE",
    region              = world_polygon
    # NOTE: if the download fails on authentication, check your installed
    # gfwr version's docs -- get_raster() may need key = gfw_auth()
    # explicitly instead of reading GFW_TOKEN automatically.
  ) |>
    filter(geartype %in% GEARTYPES) |>
    mutate(year = y)
})

saveRDS(effort_raw, gfw_effort_raw_rds)
message("Raw GFW data saved to: ", gfw_effort_raw_rds)

  #----------------------------------------------------
  # 4. Aggregation: mean annual hours per cell
  # ----------------------------------------------------
# IMPORTANT: GFW only returns rows for cell/year combinations with detected
# activity -> a cell active in 3 of 9 years has NO rows for the other 6.
# Averaging directly over the rows present divides by 3 instead of 9,
# inflating "typical annual intensity". The (cell x year) grid is explicitly
# completed with 0 before averaging, so the denominator is always
# length(years).

aggregate_effort <- function(df, all_years) {
  df |>
    group_by(Lon, Lat, year) |>
    summarise(yearly_hours = sum(`Apparent Fishing Hours`, na.rm = TRUE),
              .groups = "drop") |>
    complete(nesting(Lon, Lat), year = all_years, fill = list(yearly_hours = 0)) |>
    group_by(Lon, Lat) |>
    summarise(mean_yearly_hours = mean(yearly_hours, na.rm = TRUE),
              .groups = "drop")
}

to_raster <- function(df, mask_rast, layer_name) {
  v <- vect(df, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  v <- project(v, crs(mask_rast))
  r <- rasterize(v, mask_rast, field = "mean_yearly_hours", fun = "sum", background = 0)
  r <- mask(r, mask_rast)
  names(r) <- layer_name
  r
}

# --- One raster per gear type in GEARTYPES (generic, not hardcoded) ---
effort_by_gear <- map(GEARTYPES, function(g) aggregate_effort(filter(effort_raw, geartype == g), all_years = years))
names(effort_by_gear) <- GEARTYPES

rast_by_gear <- map2(effort_by_gear, GEARTYPES, function(df, g) {
  to_raster(df, ocean_mask, paste0("mean_yearly_hours_", g))
})
names(rast_by_gear) <- GEARTYPES

# --- Combined, ONLY if more than one gear type is active ---
if (length(GEARTYPES) > 1) {
  effort_combined <- aggregate_effort(effort_raw, all_years = years)
  rast_combined <- to_raster(effort_combined, ocean_mask, "mean_yearly_hours_combined")
}

# ----------------------------------------------------
# 5. Save effort rasters and effort DENSITY rasters (effort / effective area)
# ----------------------------------------------------
iwalk(rast_by_gear, function(r, g) {
  writeRaster(r, paste(fe_dir, sprintf("gfw_effort_%s.tif", g), sep = "/"), overwrite = TRUE)
  
  density_g <- ifel(eff_area_km2 > 0, r / eff_area_km2, NA)  # guard against area = 0/NA (avoids Inf)
  writeRaster(density_g, paste(fe_dir, sprintf("gfw_effort_density_%s.tif", g), sep = "/"), overwrite = TRUE)
})

if (length(GEARTYPES) > 1) {
  writeRaster(rast_combined, paste(fe_dir, "gfw_effort_combined.tif", sep = "/"), overwrite = TRUE)
  rast_stack <- c(rast_combined, rast(rast_by_gear))
  writeRaster(rast_stack, paste(fe_dir, "gfw_effort_stack.tif", sep = "/"), overwrite = TRUE)
  
  density_combined <- ifel(eff_area_km2 > 0, rast_combined / eff_area_km2, NA)
  writeRaster(density_combined, paste(fe_dir, "gfw_effort_density_combined.tif", sep = "/"), overwrite = TRUE)
}

message("GFW effort and density rasters saved for: ", paste(GEARTYPES, collapse = ", "))
message("See 03d_global_effort_density_map.R for the corresponding maps.")
