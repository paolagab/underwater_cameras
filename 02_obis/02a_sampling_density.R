#-----------------------------------------------------------------------------
# 02a_sampling_density.R
# Description: Computes OBIS sampling density (all years combined), explores
#              its distribution to inform the classification threshold used
#              in 02b_undersampled_areas.R, and reports temporal persistence
#              (recency of records) as descriptive context only -- it is not
#              used as a classification criterion.
#
# Threshold candidates evaluated (log10-transformed density, sampled cells
# only -- i.e. sampling_density > 0; never-sampled cells are a separate
# category handled in 02b, not part of this distribution):
#   - Tercile (p33 / p66)        [primary candidate -- see Mestre et al. 2025,
#                                  Biol. Conserv., for precedent: shipless
#                                  areas defined via first-tercile cutoff]
#   - Median                     [robust to residual skew]
#   - Mean                       [geometric mean of density]
#   - Head/Tail Breaks           [natural-break reference]
#   - 2-component Gaussian mixture [crossing point between two populations,
#                                    only meaningful if BIC actually favours G>1]
#-----------------------------------------------------------------------------

source("R/data_paths.R")

library(terra)
library(classInt)
library(mclust)
library(moments)   # skewness()

diagnostics_dir <- paste(obis_dir, "diagnostics", sep = "/")
if (!dir.exists(diagnostics_dir)) dir.create(diagnostics_dir, recursive = TRUE)

sampling_density_tif <- paste(obis_dir, "obis_sampling_density.tif", sep = "/")
persistence_descriptive_csv <- paste(diagnostics_dir, "sampling_persistence_descriptive.csv", sep = "/")

# ----------------------------------------------------
# 0. Parameters
# ----------------------------------------------------
RECENCY_START_YEAR <- 2016   # descriptive cutoff only, aligned with GFW data start

# Derived automatically from the annual stats CSV already saved by
# 01_obis_query.R, rather than hardcoded -- avoids the two files silently
# falling out of sync whenever the processed year range changes (e.g. the
# gap-free-coverage floor applied in 01_obis_query.R).
years_stats <- read.csv(obis_annual_processing_stats)
YEARS_RANGE <- min(years_stats$year):max(years_stats$year)

# ----------------------------------------------------
# 1. Load already-computed inputs (do not recompute)
# ----------------------------------------------------
obis         <- rast(obis_count)      # multiband stack, one layer per year
ocean_mask   <- rast(temp_mask)
eff_area_km2 <- rast(eff_area_km2_nc)

n_years <- nlyr(obis)
if (n_years != length(YEARS_RANGE)) {
  stop(sprintf(
    "OBIS stack has %d layers but obis_annual_processing_stats.csv implies %d years (%d-%d). These should always match -- check that obis_count and the stats CSV came from the same 01_obis_query.R run.",
    n_years, length(YEARS_RANGE), min(YEARS_RANGE), max(YEARS_RANGE)
  ))
}

# The exported .nc does not preserve layer names (Y2000, Y2001...) as
# metadata -- rebuild them from the known year range instead of parsing
# generic "Band1", "Band2"... names.
names(obis) <- paste0("Y", YEARS_RANGE)

message(sprintf("OBIS stack: %d years detected (%d-%d)", n_years, min(YEARS_RANGE), max(YEARS_RANGE)))

# ----------------------------------------------------
# 2. Sampling density (all years combined)
# ----------------------------------------------------
message("Computing cumulative sampling density...")

total_events  <- app(obis, fun = function(x) sum(x, na.rm = TRUE))
years_sampled <- app(obis, fun = function(x) sum(x > 0, na.rm = TRUE))

sampling_density <- ifel(eff_area_km2 > 0, total_events / eff_area_km2, NA)

total_events     <- mask(total_events,     ocean_mask)
years_sampled     <- mask(years_sampled,    ocean_mask)
sampling_density <- mask(sampling_density, ocean_mask)

# Recent-period sampling (>= RECENCY_START_YEAR), descriptive only
recent_idx           <- which(YEARS_RANGE >= RECENCY_START_YEAR)
years_sampled_recent <- app(obis[[recent_idx]], fun = function(x) sum(x > 0, na.rm = TRUE))
years_sampled_recent <- mask(years_sampled_recent, ocean_mask)

names(total_events)         <- "total_events"
names(sampling_density)     <- "sampling_density"
names(years_sampled)        <- "years_sampled"
names(years_sampled_recent) <- "years_sampled_recent"

# Save for 02b to reuse (no need to recompute total_events/density there)
stack_out <- c(total_events, sampling_density, years_sampled, years_sampled_recent)
writeRaster(stack_out, sampling_density_tif, overwrite = TRUE)
message("Sampling density stack saved to: ", sampling_density_tif)

# ----------------------------------------------------
# 3. Distribution shape (raw vs log), sampled cells only
# ----------------------------------------------------
density_vals <- values(sampling_density, na.rm = TRUE)[, 1]
density_vals_sampled <- density_vals[density_vals > 0]
log_density_sampled  <- log10(density_vals_sampled)

n_sampled_cells <- length(density_vals_sampled)
message(sprintf("Cells with some sampling: %d", n_sampled_cells))

cat("\n========== DISTRIBUTION SHAPE ==========\n")
cat(sprintf("Raw scale  -> skewness: %.3f | mean: %.2f | median: %.2f\n",
            skewness(density_vals_sampled), mean(density_vals_sampled), median(density_vals_sampled)))
cat(sprintf("log10 scale -> skewness: %.3f | mean: %.4f | median: %.4f\n",
            skewness(log_density_sampled), mean(log_density_sampled), median(log_density_sampled)))
cat("=========================================\n\n")

# ----------------------------------------------------
# 4. Candidate thresholds (log10 scale, back-transformed to density)
# ----------------------------------------------------
thr_p33    <- quantile(log_density_sampled, 1/3, na.rm = TRUE)
thr_p66    <- quantile(log_density_sampled, 2/3, na.rm = TRUE)
thr_median <- median(log_density_sampled)
thr_mean   <- mean(log_density_sampled)

ht <- classIntervals(log_density_sampled, style = "headtails")
thr_headtails <- ht$brks[2]

bic <- mclustBIC(log_density_sampled, G = 1:4)
cat("\n========== BIC BY NUMBER OF COMPONENTS ==========\n")
print(summary(bic))
cat("====================================================\n")

best_model <- Mclust(log_density_sampled, x = bic)
message(sprintf("Optimal number of components per BIC: %d", best_model$G))

thr_mixture <- NA
if (best_model$G >= 2) {
  mod2 <- Mclust(log_density_sampled, G = 2)
  xgrid <- seq(min(log_density_sampled), max(log_density_sampled), length.out = 3000)
  cls <- predict(mod2, newdata = xgrid)$classification
  low_comp <- which.min(mod2$parameters$mean)
  change_idx <- which(diff(cls == low_comp) != 0)
  if (length(change_idx) > 0) thr_mixture <- xgrid[change_idx[1]]
  cat(sprintf("\n2-component mixture: means = %.3f and %.3f (log10)\n",
              mod2$parameters$mean[1], mod2$parameters$mean[2]))
} else {
  cat("\nBIC does NOT favour a 2+ component mixture: no evidence of a natural\n")
  cat("two-population split. The mixture threshold would not be justified here.\n")
}

pct_below <- function(thr_log) round(100 * mean(log_density_sampled <= thr_log), 1)

comparison <- data.frame(
  method = c("p33_log (tercile, low)", "p66_log (tercile, high)", "median_log", "mean_log",
             "headtails_log", "mixture_2comp"),
  threshold_log10 = c(thr_p33, thr_p66, thr_median, thr_mean, thr_headtails, thr_mixture),
  threshold_density = 10^c(thr_p33, thr_p66, thr_median, thr_mean, thr_headtails, thr_mixture),
  pct_sampled_cells_below = c(
    pct_below(thr_p33), pct_below(thr_p66), pct_below(thr_median), pct_below(thr_mean),
    pct_below(thr_headtails), if (!is.na(thr_mixture)) pct_below(thr_mixture) else NA
  )
)

cat("\n========== THRESHOLD CANDIDATES ==========\n")
print(comparison, row.names = FALSE)
cat("=============================================\n")

write.csv(comparison, density_threshold_diagnostics_csv, row.names = FALSE)

# ----------------------------------------------------
# 5. Histogram with candidates overlaid
# ----------------------------------------------------
png(density_threshold_diagnostics_png, width = 2000, height = 1400, res = 200)

hist(log_density_sampled, breaks = 60, col = "grey80", border = "white",
     main = "log10(sampling density) - cells with events > 0",
     xlab = "log10(events / km2)")

abline(v = thr_p33,       col = "steelblue",  lwd = 2, lty = 1)
abline(v = thr_p66,       col = "darkblue",   lwd = 2, lty = 1)
abline(v = thr_median,    col = "blue",       lwd = 2, lty = 2)
abline(v = thr_mean,      col = "red",        lwd = 2, lty = 2)
abline(v = thr_headtails, col = "darkgreen",  lwd = 2, lty = 3)
if (!is.na(thr_mixture)) abline(v = thr_mixture, col = "purple", lwd = 2, lty = 1)

legend("topright",
       legend = c("p33 (tercile low)", "p66 (tercile high)", "median", "mean", "headtails",
                  if (!is.na(thr_mixture)) "mixture 2 comp." else NULL),
       col = c("steelblue", "darkblue", "blue", "red", "darkgreen",
               if (!is.na(thr_mixture)) "purple" else NULL),
       lwd = 2, lty = c(1, 1, 2, 2, 3, if (!is.na(thr_mixture)) 1 else NULL),
       bty = "n")

dev.off()

# ----------------------------------------------------
# 6. Temporal persistence (DESCRIPTIVE ONLY -- not a classification criterion)
# ----------------------------------------------------
message("Computing descriptive persistence statistics...")

df_persist <- as.data.frame(c(total_events, years_sampled_recent), xy = TRUE, na.rm = TRUE)
names(df_persist) <- c("x", "y", "total_events", "years_sampled_recent")

ever_sampled <- df_persist$total_events > 0
pct_no_recent_records <- 100 * mean(df_persist$years_sampled_recent[ever_sampled] == 0)
n_ever_sampled <- sum(ever_sampled)
n_no_recent <- sum(df_persist$years_sampled_recent[ever_sampled] == 0)

cat(sprintf(
  "\n[DESCRIPTIVE] Of %d cells with some sampling history, %.1f%% (%d cells) have no\nOBIS records since %d. Not used to classify undersampled areas.\n\n",
  n_ever_sampled, pct_no_recent_records, n_no_recent, RECENCY_START_YEAR
))

persistence_summary <- data.frame(
  n_ever_sampled_cells = n_ever_sampled,
  n_no_records_since_recency_year = n_no_recent,
  pct_no_records_since_recency_year = pct_no_recent_records,
  recency_start_year = RECENCY_START_YEAR
)
write.csv(persistence_summary, persistence_descriptive_csv, row.names = FALSE)

message("Diagnostics complete. Review before fixing DENSITY_METHOD in 02b_undersampled_areas.R:")
message(" - ", density_threshold_diagnostics_csv)
message(" - ", density_threshold_diagnostics_png)
message(" - ", persistence_descriptive_csv)
