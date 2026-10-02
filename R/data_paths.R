#----------------------------------------------------------------
# data_paths.R
# Defines paths to raw datasets, intermediate products, and outputs
# generated throughout the project. Lives in the R/ project folder,
# alongside utils.R (helper functions used across scripts).
#
#----------------------------------------------------------------

# Projections ----------------------------------------------------------

PROJ <- "+proj=moll +ellps=WGS84"  # Target equal-area projection, Mollweide
GEO  <- "+proj=longlat +ellps=WGS84"  # Geographic (unprojected) reference

# Main directories ----------------------------------------------------

# Machine-specific base path. Replace with a relative/portable path
# before publishing to GitHub (no absolute local paths in the public repo).
user <- "paola"
if (user == "paola") main_path <- "C:/Users/paoga/OneDrive - Universitat de València/underwater_camerasv2/"

raw_dir  <- paste0(main_path, "data/raw")   # raw input data
temp_dir <- paste0(main_path, "data/temp")  # intermediate/temporary products
out_dir  <- paste0(main_path, "data/out")   # final outputs (rasters, tables, figures)

for (d in c(raw_dir, temp_dir, out_dir, fig_dir, tbl_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Raw input data --------------------------------------------------------
# Keep all raw data under the same folder structure, or adjust the
# paths below accordingly.

# Bathymetry
# The GEBCO 2025 Grid, www.gebco.net. GEBCO Compilation Group (2025)
# GEBCO 2025 Grid (doi:10.5285/37c52e96-24ea-67cee063-7086abc05f29)
gebco_nc <- paste(raw_dir, "bathymetry/GEBCO_2025.nc", sep = "/")

# Land mask (Natural Earth, www.naturalearthdata.com)
ne_shp     <- paste(raw_dir, "landmask/ne_10m_land.shp", sep = "/")  # full detail, use for regional/coastal work
ne_50m_shp <- paste(raw_dir, "landmask/ne_50m_land.shp", sep = "/")  # simplified, use for global-scale maps

# Ocean sectors (SeaVoX Salt and Fresh Water Body Gazetteer, v16 2015,
# marineregions.org)
seavox_dir <- paste(raw_dir, "SeaVoX_sea_areas_polygons_v16", sep = "/")
seavox_shp <- "SeaVoX_v16_2015"

# Ocean basin polygons (simplified geometry, used for basin-level
# breakdowns of sampling/effort/overlap statistics)
atlantic <- paste(raw_dir, "oceans/atlantic_simple.shp", sep = "/")
indian   <- paste(raw_dir, "oceans/simple.shp", sep = "/")
pacific  <- paste(raw_dir, "oceans/simple_complete.shp", sep = "/")
arctic   <- paste(raw_dir, "oceans/arctic_simple.shp", sep = "/")
southern <- paste(raw_dir, "oceans/southern_simple.shp", sep = "/")

goas_path <- paste(raw_dir, "oceans/goas_v01.shp", sep = "/")

# EEZ and High Seas polygons (Marine Regions / Flanders Marine Institute)
eez_gpkg       <- paste(raw_dir, "marine_regions/eez_v12.gpkg", sep = "/")
high_seas_gpkg <- paste(raw_dir, "marine_regions/High_Seas_v2.gpkg", sep = "/")
eez_shp_path <- paste(raw_dir, "marine_regions/eez_lowresolution.shp", sep = "/")
high_seas_shp_path <- paste(raw_dir, "marine_regions/High_Seas_v1_geom_corrected.shp", sep = "/")



### Do not edit below this line without checking dependent scripts ###

# Output subfolders -------------------------------------------------
# One subfolder per pipeline component, all under out_dir.

obis_dir    <- paste(out_dir, "obis", sep = "/")             # OBIS sampling outputs
fe_dir      <- paste(out_dir, "fishing_effort", sep = "/")   # GFW fishing effort outputs
overlap_dir <- paste(out_dir, "overlap", sep = "/")          # overlap + flag/EEZ outputs

for (d in c(obis_dir, fe_dir, overlap_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Intermediate (temp) products -----------------------------------------
# Study area / ocean mask construction

temp_bathy   <- paste(temp_dir, "bathy_1d_moll.nc", sep = "/")
temp_mask    <- paste(temp_dir, "mask_1d_moll.nc", sep = "/")
temp_land    <- "land_moll"  # base name, written with st_write() as temp_land.shp
ocean_fraction_nc  <- paste(temp_dir, "ocean_fraction_moll.nc", sep = "/")
eff_area_km2_nc    <- paste(temp_dir, "effective_ocean_area_km2.nc", sep = "/")
temp_ocean_basins <- paste(temp_dir, "gos_oceans_basins.rds", sep = "/")
# OBIS pipeline outputs (obis_dir) --------------------------------------

# Query/rasterize step
query_dir                    <- paste(obis_dir, "query", sep = "/")
obis_count                   <- paste(obis_dir, "obis_count.nc", sep = "/")
obis_annual_processing_stats <- paste(query_dir, "obis_annual_processing_stats.csv", sep = "/")

# Density diagnostics (threshold selection)
density_threshold_diagnostics_csv <- paste(obis_dir, "density_threshold_diagnostics.csv", sep = "/")
density_threshold_diagnostics_png <- paste(obis_dir, "density_distribution_diagnostics.png", sep = "/")

# Sampling classification (low-sampling / never-sampled areas)
# NOTE: filenames below still reflect the earlier median-based approach;
# will be revisited once the tercile-based classification script is updated.
obis_undersampling_metrics <- paste(obis_dir, "obis_undersampling_metrics.tif", sep = "/")
undersamp_sf_rds           <- paste(obis_dir, "undersamp_sf.rds", sep = "/")
undersampling_global_stats <- paste(obis_dir, "undersampling_global_stats.csv", sep = "/")
map_undersampling_png      <- paste(obis_dir, "map_undersampling_categories.png", sep = "/")

# Basin-level breakdown
coldspot_stats_by_basin        <- paste(obis_dir, "coldspot_stats_by_basin.csv", sep = "/")
coldspot_summary_by_basin      <- paste(obis_dir, "coldspot_summary_simple_by_basin.csv", sep = "/")
coldspot_by_basin_barplot_png  <- paste(obis_dir, "coldspot_by_basin_barplot.png", sep = "/")

# Global sampling density map (all years, cumulative)
grid_density_all_time_sf <- paste(obis_dir, "grid_density_all_time_sf.rds", sep = "/")
map_density_all_time_png <- paste(obis_dir, "obis_pelagic_data_log_density_all_time.png", sep = "/")

# GFW / fishing effort pipeline outputs (fe_dir) -------------------------

gfw_effort_raw_rds <- paste(fe_dir, "gfw_effort_raw.rds", sep = "/")
# Per-geartype outputs use sprintf("gfw_effort_%s.tif", geartype) /
# sprintf("gfw_effort_density_%s.tif", geartype) inside 01_download_gfw.R.

effort_threshold_diagnostics_csv <- paste(fe_dir, "effort_threshold_diagnostics.csv", sep = "/")
effort_distribution_diagnostics_png <- paste(fe_dir, "effort_distribution_diagnostics.png", sep = "/")

effort_global_stats_csv    <- paste(fe_dir, "effort_global_stats.csv", sep = "/")
effort_stats_by_basin_csv  <- paste(fe_dir, "effort_stats_by_basin.csv", sep = "/")
effort_summary_by_basin_csv <- paste(fe_dir, "effort_summary_simple_by_basin.csv", sep = "/")

map_high_effort_areas_png  <- paste(fe_dir, "map_high_effort_areas.png", sep = "/")
map_density_effort_all_time_png <- paste(fe_dir, "gfw_effort_density_all_time.png", sep = "/")

# Overlap + flag/EEZ outputs (overlap_dir) --------------------------------

overlap_stats_csv          <- paste(overlap_dir, "overlap_stats.csv", sep = "/")
overlap_index_by_basin_csv <- paste(overlap_dir, "overlap_index_by_basin.csv", sep = "/")
overlap_cells_rds          <- paste(overlap_dir, "overlap_cells.rds", sep = "/")
map_overlap_png            <- paste(overlap_dir, "map_overlap.png", sep = "/")

gfw_effort_by_flag_raw_rds <- paste(overlap_dir, "gfw_effort_by_flag_raw.rds", sep = "/")
flag_eez_summary_raw_csv   <- paste(overlap_dir, "flag_eez_summary_raw.csv", sep = "/")
flag_ranking_eez_csv       <- paste(overlap_dir, "flag_ranking_eez.csv", sep = "/")
flag_ranking_highseas_csv  <- paste(overlap_dir, "flag_ranking_highseas.csv", sep = "/")
barplot_flag_jurisdiction_png <- paste(overlap_dir, "barplot_flag_by_jurisdiction.png", sep = "/")


