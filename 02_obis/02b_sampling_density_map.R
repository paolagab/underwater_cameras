#-----------------------------------------------------------------------------
# 02b_sampling_density_map.R
# Description: Continuous map of cumulative OBIS sampling density (all
#              years combined), complementing the categorical map from
#              02b_undersampled_areas.R. Reuses the density raster already
#              computed and saved by 02a_sampling_density.R (no recompute).
#-----------------------------------------------------------------------------

library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(paletteer)

source("R/utils.R")
source("R/data_paths.R")

# ----------------------------------------------------
# 1. Load already-computed data (do not recompute)
# ----------------------------------------------------
ocean_mask <- rast(temp_mask)
PROJ <- crs(ocean_mask)

sampling_density_tif <- paste(obis_dir, "obis_sampling_density.tif", sep = "/")
sampling_density <- rast(sampling_density_tif)[["sampling_density"]]

# ----------------------------------------------------
# 2. Log10 transform for mapping
# ----------------------------------------------------
message("Preparing log10 density layer...")

# Cells with 0 events -> NA (avoids log10(0) = -Inf and keeps the map clean)
sampling_density[sampling_density == 0] <- NA

log_density_raster <- log10(sampling_density)
names(log_density_raster) <- "log_density"

# ----------------------------------------------------
# 3. Data preparation for ggplot
# ----------------------------------------------------
message("Preparing map layers...")

df_density <- as.data.frame(log_density_raster, xy = TRUE, na.rm = TRUE)
df_blank   <- as.data.frame(ocean_mask, xy = TRUE, na.rm = TRUE)

land <- st_read(ne_50m_shp, quiet = TRUE)  # 50m: enough detail for a global-scale map, renders faster
land_moll <- st_transform(land, crs = PROJ)

bbox <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = PROJ)

# ----------------------------------------------------
# 4. Map
# ----------------------------------------------------
message("Generating map...")

final_raster_map <- ggplot() +
  geom_raster(data = df_blank, aes(x = x, y = y), fill = "grey85") +
  geom_raster(data = df_density, aes(x = x, y = y, fill = log_density)) +
  scale_fill_gradientn(
    colours = rev(paletteer_c("grDevices::Spectral", 20)),
    na.value = "transparent",
    name = expression(log[10] ~ " (events / km" ^ {2} * ")")
  ) +
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
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9)
  ) +
  labs(title = "", x = "", y = "")

# ----------------------------------------------------
# 5. Export
# ----------------------------------------------------
grid_sf_backup <- as.polygons(log_density_raster, values = TRUE) %>% st_as_sf()
saveRDS(grid_sf_backup, grid_density_all_time_sf)

ggsave(
  filename = map_density_all_time_png,
  plot = final_raster_map,
  device = "png",
  dpi = 400,
  width = 24,
  height = 14,
  units = "cm"
)

message("Map saved to: ", map_density_all_time_png)
