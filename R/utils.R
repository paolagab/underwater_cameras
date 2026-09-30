#-----------------------------------------------------------------------------
# utils.R
# Description: Shared helper functions used across the project's scripts.
#              bb() builds a densified global bounding box polygon, adapted
#              from March et al. (2020):
#
#              March, D., Boehme, L., Tintoré, J., Vélez-Belchi, P.J., and
#              Godley, B.J. (2020). Towards the integration of animal-borne
#              instruments into global ocean observing systems. Global
#              Change Biology 26, 586-596. https://doi.org/10.1111/gcb.14902
#-----------------------------------------------------------------------------

library(sf)

#' Build a global bounding box polygon, densified and optionally reprojected.
#'
#' A simple 4-corner rectangle transformed into a global projection (e.g.
#' Mollweide) can produce invalid/degenerate geometry with the s2 backend
#' used by sf, especially near the poles and the antimeridian -- long,
#' widely-spaced straight edges do not represent the true curvature of the
#' projected graticule and can crash the R session outright rather than
#' throwing a normal, catchable error. Densifying the edges (adding points
#' every `by` degrees) avoids this.
#'
#' @param xmin,xmax,ymin,ymax Bounding box limits, in degrees (geographic).
#' @param crs Target CRS (PROJ string, EPSG code, or crs object). Defaults
#'   to WGS84 geographic (no reprojection).
#' @param by Point spacing along each edge, in degrees. Default 0.5.
#' @return An sfc polygon in the target CRS.
bb <- function(xmin, xmax, ymin, ymax,
               crs = "+proj=longlat +datum=WGS84 +no_defs", by = 0.5) {
  top    <- data.frame(lon = seq(xmin, xmax, by = by),  lat = ymax)
  bottom <- data.frame(lon = seq(xmax, xmin, by = -by), lat = ymin)
  right  <- data.frame(lon = xmax, lat = seq(ymax, ymin, by = -by))
  left   <- data.frame(lon = xmin, lat = seq(ymin, ymax, by = by))
  coords <- rbind(top, right, bottom, left, top[1, ])
  poly <- st_polygon(list(as.matrix(coords)))
  poly_sf <- st_sfc(poly, crs = 4326)
  st_transform(poly_sf, crs = crs)
}