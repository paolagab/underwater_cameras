#-----------------------------------------------------------------------------
# 03_overlap_bivariate_map.R
# Description: "Impact figure" combining undersampled OBIS areas and high
#              GFW fishing effort areas into a single panel.
#
#              Layers (back to front):
#              - AOO footprint (any recorded fishing activity, not just
#                high effort): dot pattern fill. Quiet on purpose --
#                context, not priority signal.
#              - Undersampled / high effort / overlap: explicit discrete
#                categories with solid, forced colours -- NOT alpha
#                blending. An optical blend (semi-transparent layers
#                stacked) was tried first, but the resulting overlap tone
#                could not be controlled precisely and was hard to
#                distinguish at map scale. Assigning the overlap its own
#                fixed colour (#8B668B) guarantees it reads clearly and
#                consistently regardless of the other two colours chosen.
#-----------------------------------------------------------------------------

source("R/data_paths.R")
source("R/utils.R")

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(ggpattern)  # install.packages("ggpattern") if needed; requires ggplot2 >= 4.0.2

# ----------------------------------------------------
# 0. Parameters
# ----------------------------------------------------
GEARTYPE <- "drifting_longlines"
AOO_PATTERN <- "circle"  # "circle" for dots, "stripe" for stripes
CATEGORY_ALPHA <- 0.7

# ----------------------------------------------------
# 1. Load already-computed classifications (do not recompute)
# ----------------------------------------------------
ocean_mask <- rast(temp_mask)
PROJ <- crs(ocean_mask)

undersamp_class <- rast(paste(obis_dir, "obis_undersampled_areas.tif", sep = "/"))[["sampling_class"]]
effort_class    <- rast(paste(fe_dir, sprintf("gfw_high_effort_areas_%s.tif", GEARTYPE), sep = "/"))[["effort_class"]]

# 02b lookup: 1 = never_sampled, 2 = low, 3 = medium, 4 = high
# 03b lookup: 1 = no_activity,   2 = low, 3 = medium, 4 = high

is_undersampled <- undersamp_class %in% c(1, 2)
is_high_effort  <- effort_class == 4
overlap         <- is_undersampled & is_high_effort

# ----------------------------------------------------
# 2. AOO footprint (pattern layer, unchanged)
# ----------------------------------------------------
message("Preparing map layers...")

df_blank <- as.data.frame(ocean_mask, xy = TRUE, na.rm = TRUE)

is_any_activity <- effort_class != 1   # AOO -- any recorded effort
is_any_activity <- mask(is_any_activity, ocean_mask)
aoo_poligonos <- as.polygons(is_any_activity, dissolve = TRUE) |> st_as_sf()
aoo_poligonos <- aoo_poligonos[aoo_poligonos[[1]] == 1, ]  # keep TRUE (any activity) only
aoo_poligonos$layer_label <- "LLD footprint"


# ----------------------------------------------------
# 3. Single categorical layer & valid data separation
# ----------------------------------------------------
category_rast <- ifel(
  overlap, 3,
  ifel(is_undersampled & !is_high_effort, 1,
       ifel(is_high_effort & !is_undersampled, 2, NA))
)
category_rast <- mask(category_rast, ocean_mask)

# Convert ocean_mask to polygons for the white ocean background
ocean_poly <- as.polygons(ocean_mask > -Inf, dissolve = TRUE) %>% 
  st_as_sf()

# Extract only valid categories (drop NA)
df_category <- as.data.frame(category_rast, xy = TRUE, na.rm = TRUE)
names(df_category) <- c("x", "y", "category_int")

category_labels <- c("Undersampled", "High fishing effort", "Overlap")
df_category$category <- factor(category_labels[df_category$category_int], levels = category_labels)

land <- st_read(ne_50m_shp, quiet = TRUE)
land_moll <- st_transform(land, crs = PROJ)

bbox <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs = PROJ)


# ----------------------------------------------------
# 4. Map
# ----------------------------------------------------
message("Generating bivariate impact map...")

category_colors <- c(
  "Undersampled"                = "#8DA9C4",  # same blue as 02b
  "High fishing effort" = "#FFA07A",  # # same dark orange as 03b
  "Overlap"  = "#8B5F65"   # forced overlap color
)


impact_map <- ggplot() +
  
  # Layer 1: Fill the entire bounding box area with grey30 first 
  # (This acts as the land/Caspian Sea base layer underneath)
  geom_sf(data = bbox, fill = "grey30", color = NA) +
  
  # Layer 2: Draw the valid ocean polygon in white OVER the grey box 
  # (This turns the oceans white and leaves land/Caspian Sea as grey30)
  geom_sf(data = ocean_poly, fill = "white", color = NA) +
  
  # Layer 3: AOO footprint -- pattern fill
  geom_sf_pattern(
    data = aoo_poligonos,
    aes(pattern_fill = layer_label),
    pattern = AOO_PATTERN,
    pattern_color = "grey50",
    pattern_density = 0.35,
    pattern_size = 0.03,
    pattern_spacing = 0.012,
    fill = "transparent",
    color = NA
  ) +
  scale_pattern_fill_manual(values = setNames("grey50", unique(aoo_poligonos$layer_label)), name = NULL) +
  
  # Layer 4: Plot valid categories only
  geom_raster(data = df_category, aes(x = x, y = y, fill = category), alpha = CATEGORY_ALPHA) +
  scale_fill_manual(
    values = category_colors, 
    name = NULL
  ) +
  
  # Layer 5: Land shapes and bounding box border on top
  geom_sf(data = land_moll, fill = "grey30", color = NA) +
  geom_sf(data = bbox, fill = NA, color = "grey40", linewidth = 0.5) +
  
  coord_sf(datum = NA) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # Panel background is now clean white
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "right",
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    legend.text      = element_text(size = 10)
  ) +
  labs(title = "", x = "", y = "")


# ----------------------------------------------------
# 5. Export
# ----------------------------------------------------
impact_map_path <- paste(overlap_dir, sprintf("2_map_overlap_bivariate_%s.png", GEARTYPE), sep = "/")
ggsave(impact_map_path, plot = impact_map, device = "png", dpi = 400, width = 24, height = 14, units = "cm")

message("Map saved: ", impact_map_path)
