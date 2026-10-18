#!/usr/bin/env Rscript
# ==============================================================================
# 11_sens_daily_drivers.R -- sensitivity: DAILY vs MONTHLY drivers (averaging check)
# ------------------------------------------------------------------------------
# The stand predictions drive the forests with monthly climatologies. A forest is not
# linear, so the flux at the mean driver is not the mean flux over the days (Jensen).
# Here the same forests are driven, for every day of the climatology years, with that
# day's drivers, and the daily predictions averaged; the ratio to the monthly-driver
# prediction measures the averaging error.
#
# DAILY DRIVERS (CLIM_YEARS):
#   air          the tower 24-h mean of the local day (>= 20 h; 00_soil_temp_climatology.R);
#                not a predictor of the locked forests (dropped 2026-10-17), carried for the
#                with-air forests only
#   soil temp    that script's fitted relation on trailing air (outputs/tables/soil_temp_daily.csv)
#   moisture     the daily water-balance store through the same transfer as the monthly
#                level (outputs/tables/moisture_daily.csv, 02_moisture_climatology.R),
#                x the fixed spatial pattern
# Days without all three are left out; the daily predictions are averaged within each
# calendar month and the months averaged, so each month weighs as in the monthly run.
# MONTHLY DRIVERS: the main drivers, as 01_predict_tree_flux.R / 02_predict_soil_surface.R.
#
# SOIL: every cell of the stand surface. STEMS: a seeded random subsample of N_STEMS
# in-stand inventory stems, 0-2 m band (the band integral of 01_predict_tree_flux.R),
# calibrated, area-weighted. Run sizes from code/lib/run_mode.R: N_STEMS 300 and every
# day in the full run (lite: fewer stems, every 30th day).
#
# Writes outputs/data/sens_daily_drivers.csv, outputs/audit/sens_daily_drivers.txt
# ==============================================================================
source("code/lib/outputs.R"); source("code/lib/run_mode.R")
run_setup <- function(f) {
  L <- readLines(f); m <- grep("^# ---- END OF SETUP", L)[1]; stopifnot(!is.na(m))
  eval(parse(text = L[seq_len(m - 1)]), envir = globalenv())
}
N_STEMS <- RUN$daily_driver_stems; DAY_STEP <- RUN$daily_driver_day_step; SEED <- 20261017
CONV <- 86400 * 365.25 * 16e-6

# ---- daily drivers ---------------------------------------------------------------------
STD <- read.csv("outputs/tables/soil_temp_daily.csv", stringsAsFactors = FALSE)
MD  <- read.csv("outputs/tables/moisture_daily.csv", stringsAsFactors = FALSE)
DD <- merge(STD[, c("date", "air_C", "soil_temp_C", "window_days")], MD[, c("date", "moisture")], by = "date")
DD <- DD[is.finite(DD$air_C) & is.finite(DD$soil_temp_C) & is.finite(DD$moisture), ]
DD <- DD[seq(1, nrow(DD), by = DAY_STEP), ]
DD$month <- as.integer(substr(DD$date, 6, 7))
stopifnot(setequal(DD$month, 1:12))
cat(sprintf("daily drivers: %d days with air, soil temperature and moisture (%s to %s); soil-temperature window %d d\n",
            nrow(DD), min(DD$date), max(DD$date), DD$window_days[1]))
month_mean <- function(v, mo) mean(tapply(v, mo, mean))     # months weigh equally

# ---- soil --------------------------------------------------------------------------------
run_setup("code/06_upscale/02_predict_soil_surface.R")
soil_monthly <- mean(sapply(1:12, function(m)
  mean(soil_predict(shape * CLIM$moisture[CLIM$mo == m], DR$soil_temp_C_mean[DR$month == m][1],
                    DR$air_temp_C_mean[DR$month == m][1], m))))
soil_day <- vapply(seq_len(nrow(DD)), function(i)
  mean(soil_predict(shape * DD$moisture[i], DD$soil_temp_C[i], DD$air_C[i], DD$month[i])), numeric(1))
soil_daily <- month_mean(soil_day, DD$month)
B <- read.csv("outputs/data/canonical_budget.csv", stringsAsFactors = FALSE)
soil_frac <- B$value[B$quantity == "soil_area_fraction"]

# ---- stems (seeded subsample) -------------------------------------------------------------
run_setup("code/06_upscale/01_predict_tree_flux.R")
set.seed(SEED)
smp <- sample(which(INV$in_stand), N_STEMS)
F <- array(NA_real_, c(nrow(INV), K, 12))
for (mo in 1:12) for (k in seq_len(K)) F[smp, k, mo] <- pred_at(edges[k], mo, smp)
w <- INV$A_stem_m2[smp]
band_mon <- band_of(F)[smp, ]                                   # stems x months, uncalibrated
tree_monthly <- sum(rowMeans(band_mon) * INV$cal[smp] * w) / sum(w)
# daily: all days at once per height interval
nd_rows <- rep(smp, times = nrow(DD)); dd_rows <- rep(seq_len(nrow(DD)), each = N_STEMS)
FD <- matrix(NA_real_, length(nd_rows), K)
for (k in seq_len(K)) FD[, k] <- predict_species_marginal(
  function(nd) predict(TreeRF, nd, num.threads = RF_THREADS)$predictions,
  data.frame(dbh_m = INV$dbh_m[nd_rows], soil_moisture_at_tree = INV$vwc_rel[nd_rows] * DD$moisture[dd_rows],
             soil_temp_C = DD$soil_temp_C[dd_rows], air_temp_C = DD$air_C[dd_rows],
             height_cm = edges[k], dead_stem = 0),
  INV$sp[nd_rows], LEVEL_W, cmap)
blen <- band_len_cm[nd_rows]
acc <- pmin(blen, H_LO) * FD[, 1]
for (k in seq_len(K)) acc <- acc + pmax(0, pmin(blen, edges[k + 1]) - max(edges[k], H_LO)) * FD[, k]
band_day <- acc / blen * INV$cal[nd_rows]
stem_day <- tapply(band_day * INV$A_stem_m2[nd_rows], dd_rows, sum) / sum(w)   # area-weighted, per day
tree_daily <- month_mean(as.numeric(stem_day), DD$month)

OUT <- data.frame(term = c("soil (all cells; per m2 ground)", sprintf("stems (%d sampled; 0-2 m band, area-weighted)", N_STEMS)),
                  monthly_drivers = c(soil_monthly * CONV * soil_frac, tree_monthly),
                  daily_drivers   = c(soil_daily * CONV * soil_frac, tree_daily),
                  units = c("mg CH4 m-2 ground yr-1", "nmol m-2 stem s-1"),
                  n_days = nrow(DD), seed = c(NA, SEED))
OUT$ratio_daily_to_monthly <- OUT$daily_drivers / OUT$monthly_drivers
write.csv(OUT, out_path("sens_daily_drivers.csv"), row.names = FALSE)
sink(out_path("sens_daily_drivers.txt", "audit"))
cat(sprintf("Daily vs monthly drivers: %d days (%s to %s), soil temperature from a %d-day trailing-air relation\n",
            nrow(DD), min(DD$date), max(DD$date), DD$window_days[1]))
print(OUT, row.names = FALSE)
sink()
cat(readLines(out_path("sens_daily_drivers.txt", "audit")), sep = "\n")
