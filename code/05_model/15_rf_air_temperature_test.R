#!/usr/bin/env Rscript
# ==============================================================================
# 15_rf_air_temperature_test.R -- sensitivity: both forests WITH air temperature
# ------------------------------------------------------------------------------
# The locked forests (02_rf_models.R) have no air-temperature predictor (dropped
# 2026-10-17). This refits each with the air temperature measured at the deployment
# added (flux_stem / flux_soil air_temp_C) -- everything else as locked: the same rows
# and predictors, num.trees, min.node.size 5, mtry = floor(sqrt(number of predictors)),
# seed 42 -- and reports, for the locked and the with-air forests:
#   OOB R2 and grouped-CV R2 with its SD (the 100 x 5 grouped CV of 04_rf_grouped_cv.R,
#   its own function, evaluated from that file);
#   the stand budget (soil sink per m2 ground, tree band, tree % of soil), driven by the
#   24-h monthly tower climatologies of the main scaling and, for the with-air forests,
#   also by DAYTIME air (the tower mean over the measurement hours: the 5th-95th
#   percentile of the deployments' local start hours, as whole hours, CLIM_YEARS).
#   The locked forests have no air input, so the two drivers give them one value.
# Budgets use the setup sections of 06_upscale/01_predict_tree_flux.R and
# 02_predict_soil_surface.R (up to their END OF SETUP markers), with each forest's own
# out-of-bag calibration and height steps. The locked refits must reproduce the forests
# in outputs/models exactly. Nothing here writes a model to outputs/models.
#
# Writes outputs/data/rf_air_temperature_test.csv, outputs/audit/rf_air_temperature_test.txt
# ==============================================================================
source("code/lib/outputs.R"); source("code/lib/rf_spec.R")
suppressPackageStartupMessages({ library(ranger); library(dplyr) })
load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")
TreeRF_main <- TreeRF; SoilRF_main <- SoilRF
CLIM_YEARS <- 2019:2022
CONV <- 86400 * 365.25 * 16e-6                  # nmol m-2 s-1 -> mg CH4 m-2 yr-1
run_setup <- function(f) {
  L <- readLines(f); m <- grep("^# ---- END OF SETUP", L)[1]; stopifnot(!is.na(m))
  eval(parse(text = L[seq_len(m - 1)]), envir = globalenv())
}

# grouped_cv() exactly as 04_rf_grouped_cv.R defines it
for (e in parse("code/05_model/04_rf_grouped_cv.R")) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
    as.character(e[[2]]) %in% c("NREP", "NFOLD", "r2", "grouped_cv")) eval(e, globalenv())
stopifnot(exists("grouped_cv"), NREP == 100, NFOLD == 5)

fit <- function(X, y) ranger(x = as.data.frame(X), y = y, num.trees = RF_NUM_TREES, min.node.size = 5,
  mtry = floor(sqrt(ncol(X))), importance = "impurity", num.threads = RF_THREADS, oob.error = TRUE, seed = 42)
stopifnot(!"air_temp_C" %in% colnames(X_tree), !"air_temp_C" %in% colnames(X_soil))
Xs <- as.data.frame(X_soil)
SPEC <- list(
  tree_locked  = list(X = X_tree, y = tree_train_complete$y_flux, grp = tree_train_complete$tree_id),
  tree_withair = list(X = cbind(X_tree, air_temp_C = tree_train_complete$air_temp_C), y = tree_train_complete$y_flux,
                      grp = tree_train_complete$tree_id),
  soil_locked  = list(X = Xs, y = soil_train_complete$y_flux, grp = soil_train_complete$site_id),
  soil_withair = list(X = cbind(Xs, air_temp_C = soil_train_complete$air_temp_C), y = soil_train_complete$y_flux,
                      grp = soil_train_complete$site_id))
M <- lapply(SPEC, function(s) fit(s$X, s$y))
stopifnot(identical(M$tree_locked$predictions, TreeRF_main$predictions),
          identical(M$soil_locked$predictions, SoilRF_main$predictions))
CV <- bind_rows(lapply(names(SPEC), function(k) grouped_cv(SPEC[[k]]$X, SPEC[[k]]$y, SPEC[[k]]$grp, M[[k]], k, "unit")))

# ---- daytime air climatology (measurement hours) -----------------------------------------
rt <- function(t) read.csv(file.path("data/tables", paste0(t, ".csv")), colClasses = "character", na.strings = "")
st <- c(rt("flux_stem")$start_time, rt("flux_soil")$start_time)
hr <- as.integer(substr(st, 12, 13)) + as.integer(substr(st, 15, 16)) / 60
hr <- hr[is.finite(hr)]
HR_LO <- floor(quantile(hr, 0.05)); HR_HI <- floor(quantile(hr, 0.95))
w <- rt("weather")
t_beg <- as.POSIXct(substr(sub("T", " ", w$timestamp), 1, 19), tz = "UTC") - 3600   # the hour each record closes
ta <- suppressWarnings(as.numeric(w$air_temp_avg_C))
ok <- is.finite(ta) & as.integer(format(t_beg, "%Y")) %in% CLIM_YEARS & as.integer(format(t_beg, "%H")) >= HR_LO &
      as.integer(format(t_beg, "%H")) <= HR_HI
airday <- tapply(ta[ok], as.integer(format(t_beg[ok], "%m")), mean)
stopifnot(length(airday) == 12)

# ---- stand budget ---------------------------------------------------------------------------
SB <- file.path(tempdir(), "air_test_models"); dir.create(SB, showWarnings = FALSE)
B <- read.csv("outputs/data/canonical_budget.csv", stringsAsFactors = FALSE)
soil_frac <- B$value[B$quantity == "soil_area_fraction"]
budget <- function(treeRF, soilRF) {
  TreeRF <- treeRF; SoilRF <- soilRF
  save(TreeRF, SoilRF, file = file.path(SB, "RF_MODELS.RData"))
  file.copy("outputs/models/TRAINING_DATA.RData", file.path(SB, "TRAINING_DATA.RData"), overwrite = TRUE)
  Sys.setenv(TREE_PRED_SANDBOX = SB)     # 01_predict_tree_flux.R setup reads the forest (and its OOB calibration) from here
  run_setup("code/06_upscale/01_predict_tree_flux.R")
  Sys.unsetenv("TREE_PRED_SANDBOX")
  air24 <- setNames(DR$air_temp_C_mean, DR$month)
  band <- function(air) {
    DR$air_temp_C_mean <<- air[as.character(DR$month)]
    F <- array(NA_real_, c(nrow(INV), K, 12))
    for (mo in 1:12) for (k in seq_len(K)) F[, k, mo] <- pred_at(edges[k], mo)
    b <- rowMeans(band_of(F)) * INV$cal
    sum(b[INV$in_stand] * INV$A_stem_m2[INV$in_stand]) / STAND_AREA_M2 * CONV
  }
  tb <- c(band(air24), band(airday))
  run_setup("code/06_upscale/02_predict_soil_surface.R")
  SoilRF <<- soilRF                       # soil_predict() uses the global SoilRF
  sm <- function(air) mean(sapply(1:12, function(m) mean(soil_predict(shape * CLIM$moisture[CLIM$mo == m],
          DR$soil_temp_C_mean[DR$month == m][1], air[[as.character(m)]], m)))) * CONV * soil_frac
  sb <- c(sm(air24), sm(airday))
  data.frame(drivers = c("24-h monthly", sprintf("daytime %02d-%02d h", HR_LO, HR_HI)),
             soil_mg_m2_yr = sb, tree_band_mg_m2_yr = tb, tree_pct_of_soil = 100 * tb / abs(sb))
}
BUD <- rbind(cbind(spec = "locked (no air temperature)", budget(M$tree_locked, M$soil_locked)),
             cbind(spec = "with air temperature", budget(M$tree_withair, M$soil_withair)))

OUT <- bind_rows(
  CV %>% transmute(section = "skill", spec = model, drivers = NA_character_, oob_r2, grouped_r2, grouped_r2_sd),
  BUD %>% transmute(section = "budget", spec, drivers, soil_mg_m2_yr, tree_band_mg_m2_yr, tree_pct_of_soil))
write.csv(OUT, out_path("rf_air_temperature_test.csv"), row.names = FALSE)
sink(out_path("rf_air_temperature_test.txt", "audit"))
cat("Sensitivity: both forests with the measured air temperature added (locked specification otherwise)\n\n")
print(CV %>% select(model, n_rows, oob_r2, grouped_r2, grouped_r2_sd), row.names = FALSE)
cat(sprintf("\nStand budget (mg CH4 m-2 yr-1; soil per m2 ground). Daytime: tower hours beginning %02d:00-%02d:00, %s\n",
            HR_LO, HR_HI, paste(range(CLIM_YEARS), collapse = "-")))
print(BUD, row.names = FALSE)
sink()
cat(readLines(out_path("rf_air_temperature_test.txt", "audit")), sep = "\n")
