#-----------------------------------------------------------------------------
# 01b_fe_distribution_diagnostics.R
# Description: Explores the distribution of GFW fishing effort density
#              (log10-transformed) and compares threshold candidates for
#              "high effort" BEFORE committing to one. Same methodological
#              approach as 02a for OBIS density, but looking at the HEAD of
#              the distribution (high values) instead of the tail.
#
#              Primary candidate: tercile (p66 = upper third), following
#              Mestre et al. 2025 (Biol. Conserv.), who use tercile cutoffs
#              on vessel density / species richness for the same kind of
#              conflict-zone mapping. Other methods (p75/p90/p95, Head/Tail
#              Breaks, Gaussian mixture) kept for robustness/documentation.
#-----------------------------------------------------------------------------
source("R/utils.R")
source("R/data_paths.R")

library(terra)
library(classInt)
library(mclust)
library(moments)

# ----------------------------------------------------
# 0. Parameters
# ----------------------------------------------------
GEARTYPE <- "drifting_longlines"  # must match the name used in 01_download_gfw.R

# ----------------------------------------------------
# 1. Load already-computed density (do not recompute)
# ----------------------------------------------------
density_path <- paste(fe_dir, sprintf("gfw_effort_density_%s.tif", GEARTYPE), sep = "/")
effort_density <- rast(density_path)

effort_vals <- values(effort_density, na.rm = TRUE)[, 1]
effort_vals_active <- effort_vals[effort_vals > 0]   # cells with some effort only
log_effort_active <- log10(effort_vals_active)

n_active_cells <- length(effort_vals_active)
message(sprintf("Cells with some fishing effort (%s): %d", GEARTYPE, n_active_cells))

# ----------------------------------------------------
# 2. Distribution shape (raw vs log)
# ----------------------------------------------------
cat("\n========== DISTRIBUTION SHAPE ==========\n")
cat(sprintf("Raw scale   -> skewness: %.3f | mean: %.4f | median: %.4f\n",
            skewness(effort_vals_active), mean(effort_vals_active), median(effort_vals_active)))
cat(sprintf("log10 scale -> skewness: %.3f | mean: %.4f | median: %.4f\n",
            skewness(log_effort_active), mean(log_effort_active), median(log_effort_active)))
cat("===========================================\n\n")

# ----------------------------------------------------
# 3. Threshold candidates (log10, back-transformed)
# ----------------------------------------------------
thr_p33 <- quantile(log_effort_active, 0.33, na.rm = TRUE)  # kept for symmetry/documentation
thr_p66 <- quantile(log_effort_active, 0.66, na.rm = TRUE)  # PRIMARY candidate (tercile, high)
thr_p75 <- quantile(log_effort_active, 0.75, na.rm = TRUE)
thr_p90 <- quantile(log_effort_active, 0.90, na.rm = TRUE)
thr_p95 <- quantile(log_effort_active, 0.95, na.rm = TRUE)

ht <- classIntervals(log_effort_active, style = "headtails")
thr_headtails <- ht$brks[2]

bic <- mclustBIC(log_effort_active, G = 1:4)
cat("\n========== BIC BY NUMBER OF COMPONENTS ==========\n")
print(summary(bic))
cat("====================================================\n")

best_model <- Mclust(log_effort_active, x = bic)
message(sprintf("Optimal number of components per BIC: %d", best_model$G))

thr_mixture <- NA
if (best_model$G >= 2) {
  mod2 <- Mclust(log_effort_active, G = 2)
  xgrid <- seq(min(log_effort_active), max(log_effort_active), length.out = 3000)
  cls <- predict(mod2, newdata = xgrid)$classification
  high_comp <- which.max(mod2$parameters$mean)
  change_idx <- which(diff(cls == high_comp) != 0)
  if (length(change_idx) > 0) thr_mixture <- xgrid[change_idx[1]]
  cat(sprintf("\n2-component mixture: means = %.3f and %.3f (log10)\n",
              mod2$parameters$mean[1], mod2$parameters$mean[2]))
} else {
  cat("\nBIC does NOT favour a 2+ component mixture for this gear -- no evidence\n")
  cat("of a natural two-population split.\n")
}

pct_above <- function(thr_log) round(100 * mean(log_effort_active >= thr_log), 1)

comparison <- data.frame(
  method = c("p33_log (tercile, low)", "p66_log (tercile, high)", "p75_log", "p90_log", "p95_log",
             "headtails_log", "mixture_2comp"),
  threshold_log10 = c(thr_p33, thr_p66, thr_p75, thr_p90, thr_p95, thr_headtails, thr_mixture),
  threshold_hours_km2 = 10^c(thr_p33, thr_p66, thr_p75, thr_p90, thr_p95, thr_headtails, thr_mixture),
  pct_active_cells_above = c(
    pct_above(thr_p33), pct_above(thr_p66), pct_above(thr_p75), pct_above(thr_p90), pct_above(thr_p95),
    pct_above(thr_headtails), if (!is.na(thr_mixture)) pct_above(thr_mixture) else NA
  )
)

cat("\n========== THRESHOLD CANDIDATES ==========\n")
print(comparison, row.names = FALSE)
cat("=============================================\n")

write.csv(comparison, effort_threshold_diagnostics_csv, row.names = FALSE)

# ----------------------------------------------------
# 4. Histogram with candidates overlaid
# ----------------------------------------------------
png(effort_distribution_diagnostics_png, width = 2000, height = 1400, res = 200)

hist(log_effort_active, breaks = 60, col = "grey80", border = "white",
     main = sprintf("log10(effort density) - %s, active cells", GEARTYPE),
     xlab = "log10(mean yearly fishing hours / km2)")

abline(v = thr_p33,       col = "steelblue",  lwd = 2, lty = 1)
abline(v = thr_p66,       col = "darkorange", lwd = 2, lty = 1)
abline(v = thr_p90,       col = "darkred",    lwd = 2, lty = 2)
abline(v = thr_headtails, col = "darkgreen",  lwd = 2, lty = 3)
if (!is.na(thr_mixture)) abline(v = thr_mixture, col = "purple", lwd = 2, lty = 1)

legend("topright",
       legend = c("p33 (tercile low)", "p66 (tercile high)", "p90", "headtails",
                  if (!is.na(thr_mixture)) "mixture 2 comp." else NULL),
       col = c("steelblue", "darkorange", "darkred", "darkgreen",
               if (!is.na(thr_mixture)) "purple" else NULL),
       lwd = 2, lty = c(1, 1, 2, 3, if (!is.na(thr_mixture)) 1 else NULL),
       bty = "n")

dev.off()

message("\nDiagnostics complete. Review before fixing the threshold in 03b_high_effort_areas.R:")
message(" - ", effort_threshold_diagnostics_csv)
message(" - ", effort_distribution_diagnostics_png)
