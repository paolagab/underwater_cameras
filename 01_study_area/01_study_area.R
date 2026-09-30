-----------------------------------------------------------------------------
# study_area.R
# Description: Generates a global ocean mask and ancillary spatial layers
#              for subsequent analyses of biodiversity sampling effort and
#              fishing activity. Aggregates GEBCO bathymetry from 15
#              arc-second resolution to 1°, estimates the proportion of
#              ocean surface within each cell, and removes cells with <10%
#              ocean area, continental depressions/inland water bodies, and
#              the Caspian Sea. All layers are reprojected onto a global
#              Mollweide equal-area grid, and effective ocean area (km²) is
#              calculated per analysis cell.
#
#              Outputs: bathymetry, ocean mask, ocean fraction, and
#              effective ocean area rasters (all Mollweide), plus land
#              polygons (Mollweide).
#-----------------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)

source("R/utils.R")
source("R/data_paths.R")

#------------------------------------------------------------------------------
# 0. USER PARAMETERS
#------------------------------------------------------------------------------

# Minimum ocean fraction required for a cell to be retained
exclude_area <- 0.10

# Aggregation factor:
# 15 arc-second GEBCO -> 1° grid (240 × 240 subcells)
fact <- 240

#------------------------------------------------------------------------------
# 1. LOAD INPUT DATA
#------------------------------------------------------------------------------

message("Loading GEBCO bathymetry and land polygons...")

# GEBCO bathymetry
gebco <- rast(gebco_nc)

# Natural Earth land polygons
land <- st_read(ne_shp, quiet = TRUE)
land_v <- vect(land)
crs(land_v) <- "EPSG:4326"

#------------------------------------------------------------------------------
# 2. GENERATE A CLEAN OCEAN MASK
#------------------------------------------------------------------------------

message("Building ocean mask from bathymetry (this can take a while at native GEBCO resolution)...")

# Ocean cells correspond to negative bathymetry values
om <- classify(
  gebco,
  rbind(
    c(-Inf, 0, 1),
    c(0, Inf, NA)
  )
)

# Remove inland depressions/water bodies that may have negative elevation
# (e.g. continental basins) by masking out anything covered by land polygons
om <- mask(om, land_v, inverse = TRUE)

# Apply mask to bathymetry
gebco_mask <- gebco * om

#------------------------------------------------------------------------------
# 3. AGGREGATE TO 1° RESOLUTION
#------------------------------------------------------------------------------

message("Aggregating to 1° resolution...")

# Number of ocean subcells within each 1° cell
om_1d <- aggregate(
  om,
  fact = fact,
  fun = "sum",
  na.rm = TRUE
)

# Fraction of ocean surface per 1° cell
ocean_fraction_1d <- om_1d / (fact^2)

# Mean bathymetry at 1° resolution
gebco_1d <- aggregate(
  gebco_mask,
  fact = fact,
  fun = "mean",
  na.rm = TRUE
)

#------------------------------------------------------------------------------
# 4. REMOVE CELLS WITH <10% OCEAN COVER
#------------------------------------------------------------------------------

message(sprintf("Removing cells with <%.0f%% ocean cover...", exclude_area * 100))

max_om <- global(om_1d, "max", na.rm = TRUE)[[1]]
ncells_ocean <- max_om * exclude_area

om_1d_rec <- classify(
  om_1d,
  rbind(
    c(-Inf, ncells_ocean, NA),
    c(ncells_ocean, Inf, 1)
  )
)

#------------------------------------------------------------------------------
# 5. REMOVE CELLS DOMINATED BY LAND
#------------------------------------------------------------------------------

message("Removing cells dominated by land (>90% land cover)...")

# Temporary 0.1° grid
r_highres <- rast(
  ext(om_1d),
  res = 0.1,
  crs = crs(om_1d)
)

# Rasterize land
rland_highres <- rasterize(
  land_v,
  r_highres,
  field = 1,
  background = 0
)

# Fraction of land cover per 1° cell
rland <- aggregate(
  rland_highres,
  fact = 10,
  fun = "mean"
)

# Remove cells containing >90% land
max_rland <- global(rland, "max", na.rm = TRUE)[[1]]
ncells_land <- max_rland * (1 - exclude_area)

rland_rec <- classify(
  rland,
  rbind(
    c(-Inf, ncells_land, 1),
    c(ncells_land, Inf, NA)
  )
)

#------------------------------------------------------------------------------
# 6. FINAL OCEAN MASK
#------------------------------------------------------------------------------

message("Assembling final ocean mask and excluding the Caspian Sea...")

# NOTE: renamed from "mask" to "final_ocean_mask" -- "mask" is also the name
# of the terra function used later in this script (terra::mask()). R
# resolves name(...) calls to the nearest FUNCTION binding, so the original
# code was not actually broken, but shadowing a core terra function name is
# confusing to read and easy to break in future edits.
final_ocean_mask <- om_1d_rec * rland_rec

# Remove Caspian Sea
ec <- ext(45, 56, 35, 49) + 5
final_ocean_mask[cells(final_ocean_mask, ec)] <- NA

# Apply mask to bathymetry
gebco_1d <- gebco_1d * final_ocean_mask

#------------------------------------------------------------------------------
# 7. CREATE TARGET MOLLWEIDE GRID
#------------------------------------------------------------------------------

r <- rast(
  ext(-18040096, 18040096, -9020048, 9020048),
  resolution = c(100000, 100000),
  crs = PROJ
)

#------------------------------------------------------------------------------
# 8. REPROJECT TO MOLLWEIDE
#------------------------------------------------------------------------------

message("Reprojecting to Mollweide...")

gebco_1d.prj <- project(
  gebco_1d,
  r,
  method = "bilinear"
)

mask.prj <- project(
  final_ocean_mask,
  r,
  method = "near"
)

ocean_fraction.prj <- project(
  ocean_fraction_1d,
  r,
  method = "bilinear"
)

land.prj <- st_transform(
  land,
  crs = PROJ
)

#------------------------------------------------------------------------------
# 9. CALCULATE EFFECTIVE OCEAN AREA
#------------------------------------------------------------------------------

message("Calculating effective ocean area per cell...")

# Area of each Mollweide cell (km²)
area_km2 <- cellSize(
  mask.prj,
  unit = "km"
)

# Ocean area available within each cell
effective_ocean_area <- area_km2 * ocean_fraction.prj

# Apply final ocean mask
effective_ocean_area <- mask(
  effective_ocean_area,
  mask.prj
)


# FIX: bilinear interpolation of ocean_fraction near the antimeridian seam
# (+-180 deg) can collapse to 0/NA even though the ocean mask itself (built
# with nearest-neighbour, more robust at this seam) confirms the cell is
# real ocean. Where that mismatch occurs, fall back to the cell's full
# nominal area rather than letting a real, sampled ocean cell silently
# drop out of downstream analyses via na.rm = TRUE.
n_affected <- global(mask.prj == 1 & (effective_ocean_area == 0 | is.na(effective_ocean_area)),
                     "sum", na.rm = TRUE)[1, 1]
if (n_affected > 0) {
  message(sprintf("Fixing %d ocean cell(s) with effective_ocean_area = 0/NA (antimeridian interpolation artifact)...", n_affected))
}

effective_ocean_area <- ifel(
  mask.prj == 1 & (effective_ocean_area == 0 | is.na(effective_ocean_area)),
  area_km2,
  effective_ocean_area
)

# Diagnostic check (wrapped in message() so it prints even when this script
# is sourced non-interactively, e.g. from a master pipeline script)
total_ocean_area_km2 <- global(effective_ocean_area, "sum", na.rm = TRUE)[1, 1]
message(sprintf("Total effective ocean area: %.0f km²", total_ocean_area_km2))

#------------------------------------------------------------------------------
# 10. OCEAN BASIN ASSIGNMENT 
#------------------------------------------------------------------------------

#(Global Oceans and Seas, GOaS v1 -- Flanders Marine
# Institute). Replaces the earlier IHO Sea Areas v3 workflow: GOaS already
# ships pre-aggregated into 10 named units (no need to grepl-match dozens of
# individual marginal seas, and no centroid-based lat/lon fallback -- which
# also removes the antimeridian-centroid risk that approach carried).

message("Assigning ocean basins using Global Oceans and Seas (GOaS v1)...")

sf_use_s2(FALSE)

goas_path <- paste(raw_dir, "oceans/goas_v01.shp", sep = "/")
goas_sf <- st_read(goas_path, quiet = TRUE) %>%
  st_transform("EPSG:4326") %>%
  st_make_valid()

# Direct lookup, no pattern matching needed -- GOaS's 10 units map
# deterministically to the 5 basins used throughout this project. Mediterranean
# Region and Baltic Sea are treated as Atlantic marginal seas; South China and
# Eastern Archipelagic Seas as a Pacific marginal sea (consistent with the
# grouping used in the earlier IHO-based classification).
basin_map <- c(
  "Southern Ocean"                            = "Southern",
  "South Atlantic Ocean"                      = "Atlantic",
  "North Atlantic Ocean"                      = "Atlantic",
  "Mediterranean Region"                      = "Atlantic",
  "Baltic Sea"                                 = "Atlantic",
  "South Pacific Ocean"                       = "Pacific",
  "North Pacific Ocean"                       = "Pacific",
  "South China and Easter Archipelagic Seas" = "Pacific", #Note that "Eastern" is misspelled in the original file
  "Indian Ocean"                              = "Indian",
  "Arctic Ocean"                               = "Arctic"
)

goas_sf$ocean_basin <- basin_map[goas_sf$name]

# Safety check: catches any GOaS unit whose name doesn't exactly match the
# lookup above (e.g. a version mismatch or an unexpected extra polygon),
# rather than silently producing NA basins downstream.
if (anyNA(goas_sf$ocean_basin)) {
  stop("Unmapped GOaS unit(s): ", paste(unique(goas_sf$name[is.na(goas_sf$ocean_basin)]), collapse = ", "))
}

oceans_sf_wgs84 <- goas_sf %>%
  group_by(ocean_basin) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_make_valid()

saveRDS(oceans_sf_wgs84, temp_ocean_basins)

 

# ----------------------------------------------------
# 11. Asign each grid cell to an ocean basin 
# ----------------------------------------------------
df_cells <- as.data.frame(mask.prj, xy = TRUE, na.rm = TRUE)
colnames(df_cells)[3] <- "mask_val"
df_cells$cell_id <- terra::cellFromXY(mask.prj, as.matrix(df_cells[, c("x", "y")]))

effort_sf_wgs84 <- st_as_sf(df_cells, coords = c("x", "y"), crs = PROJ, remove = FALSE) %>%
  st_transform("EPSG:4326")

effort_basin <- suppressWarnings(
  st_join(effort_sf_wgs84, oceans_sf_wgs84, join = st_intersects, left = TRUE)
)

n_dupes <- sum(duplicated(effort_basin$cell_id))
if (n_dupes > 0) {
  effort_basin <- effort_basin %>% group_by(cell_id) %>% slice(1) %>% ungroup()
}

unassigned_idx <- is.na(effort_basin$ocean_basin)
n_unassigned <- sum(unassigned_idx)
if (n_unassigned > 0) {
  message(sprintf("Resolving %d coastal orphan cells via nearest feature...", n_unassigned))
  orphan_idx <- which(unassigned_idx)
  orphan_pts <- effort_basin[orphan_idx, ]
  nearest_idx <- suppressWarnings(st_nearest_feature(orphan_pts, oceans_sf_wgs84))
  effort_basin$ocean_basin[orphan_idx] <- oceans_sf_wgs84$ocean_basin[nearest_idx]
}

basin_lookup <- effort_basin %>%
  st_drop_geometry() %>%
  select(cell_id, x, y, ocean_basin)

lookup_path <- paste(temp_dir, "cell_basin_lookup.rds", sep = "/")
saveRDS(basin_lookup, lookup_path)
message(sprintf("Saved cell-to-basin lookup table for %d cells to %s", nrow(basin_lookup), lookup_path))

sf_use_s2(TRUE)

ggplot(basin_lookup, aes(x = x, y = y, color = ocean_basin)) +
  geom_point(size = 0.4) +
  coord_fixed() +
  scale_color_manual(values = c(
    Atlantic = "#1f77b4", Pacific = "#ff7f0e", Indian = "#2ca02c",
    Arctic = "#d62728", Southern = "#9467bd"
  )) +
  theme_minimal() +
  labs(color = "Ocean basin")
#------------------------------------------------------------------------------
# 12. EXPORT OUTPUTS
#------------------------------------------------------------------------------

message("Writing outputs...")

writeCDF(
  gebco_1d.prj,
  temp_bathy,
  overwrite = TRUE
)

writeCDF(
  mask.prj,
  temp_mask,
  overwrite = TRUE
)

writeCDF(
  ocean_fraction.prj,
  ocean_fraction_nc,
  overwrite = TRUE
)

writeCDF(
  effective_ocean_area,
  eff_area_km2_nc,
  overwrite = TRUE
)

st_write(
  land.prj,
  file.path(temp_dir, paste0(temp_land, ".shp")),
  delete_layer = TRUE,
  quiet = TRUE
)

message("study_area.R complete.")
