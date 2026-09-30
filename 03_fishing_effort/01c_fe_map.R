#-----------------------------------------------------------------------------
# 01c_fe_map.R
# Description: Continuous map of GFW fishing effort density (log10),
#              complementing the categorical map from 03b_high_effort_areas.R.
#              Reuses the density raster already computed and saved by
#              01_download_gfw.R (no recompute). Mirrors the aesthetic of
#              02d_global_sampling_density_map.R for OBIS.
#-----------------------------------------------------------------------------

library(terra)
library(dplyr)
library(sf)
library(ggplot2)
library(paletteer)

source("R/utils.R")
source("R/data_paths.R")

# ----------------------------------------------------
# 0. Parameters
# ----------------------------------------------------
GEARTYPE <- "drifting_longlines"

# ----------------------------------------------------
# 1. Load already-computed data (do not recompute)
# ----------------------------------------------------
ocean_mask <- rast(temp_mask)
PROJ <- crs(ocean_mask)

effort_density <- rast(paste(fe_dir, sprintf("gfw_effort_density_%s.tif", GEARTYPE), sep = "/"))

# ----------------------------------------------------
# 2. Log10 transform for mapping
# ----------------------------------------------------
message("Preparing log10 density layer...")

effort_density[effort_density == 0] <- NA  # avoids log10(0) = -Inf

log_effort_raster <- log10(effort_density)
names(log_effort_raster) <- "log_effort"

# ----------------------------------------------------
# 3. Data preparation for ggplot
# ----------------------------------------------------
message("Preparing map layers...")

df_effort <- as.data.frame(log_effort_raster, xy = TRUE, na.rm = TRUE)
df_blank  <- as.data.frame(ocean_mask, xy = TRUE, na.rm = TRUE)

land <- st_read(ne_50m_shp, quiet = TRUE)
land_moll <- st_transform(land, crs = PROJ)


bbox <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = PROJ)

# ----------------------------------------------------
# 4. Map
# ----------------------------------------------------
message("Generating map...")

final_effort_map <- ggplot() +
  geom_raster(data = df_blank,  aes(x = x, y = y), fill = "transparent") +
  geom_raster(data = df_effort, aes(x = x, y = y, fill = log_effort)) +
  scale_fill_gradientn(
    colours  = paletteer_c("viridis::magma", 20),
    na.value = "transparent",
    name     = expression(log[10] ~ "(mean yearly fishing hours / km" ^ {2} * ")")
  ) +
  geom_sf(data = land_moll, fill = "grey30", color = NA) +
  geom_sf(data = bbox, fill = NA, color = "grey30", linewidth = 0.5) +
  coord_sf(datum = NA) +
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
ggsave(
  filename = map_density_effort_all_time_png,
  plot = final_effort_map,
  device = "png",
  dpi = 400,
  width = 24,
  height = 14,
  units = "cm"
)

message("Map saved to: ", map_density_effort_all_time_png)


