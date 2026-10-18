#!/usr/bin/env Rscript
# ==============================================================================
# 01_predict_tree_flux.R -- per-stem CH4 flux over the 0-2 m band (Fig 9b, budget)
# ------------------------------------------------------------------------------
# Stage C; run from the repository root.
#
# Inventory: outputs/tables/inventory_stems.csv (01_import/08_inventory_build.R).
# Model: the stem forest (species, DBH, soil moisture, soil temperature, height),
#   with a per-species calibration (code/lib/species_calibration.R).
# Drivers: per-stem soil moisture = the stem's relative wetness on the moisture
#   surface x the month's plot-wide level (as the soil grid); soil temperature = the
#   monthly climatology (04_drivers/00_soil_temp_climatology.R).
#
# Height: the forest was trained at 50, 125 and 200 cm, so in height it is a step
# function. The steps are detected (a 1 cm scan of 400 sampled stems), each stem is
# evaluated once per constant interval, and the band mean is the width-weighted
# mean of the intervals: exact. Below 50 cm the 50 cm value is held. Flux above 2 m
# is not predicted here: 07_scaling_full_grid.R applies its vertical scenarios to
# each stem's 2 m value and calibrated band profile, both written here.
#
# TREE_PRED_SANDBOX=<dir>: models are read from <dir> and the tables written there
# (used by 05_model/05_audit_training_population.R to score alternative forests).
#
# Outputs (outputs/tables/, or the sandbox)
#   tree_flux_predictions.csv     per stem: band mean flux, 2 m flux, band area, position
#   tree_band_profile.csv         per stem: calibrated annual-mean flux per height interval
#   tree_band_profile_edges.csv   the interval edges (cm)
#   tree_monthly_stand.csv        stand monthly series (per m2 ground; breast height per m2 bark)
# ==============================================================================
source("code/lib/rf_spec.R")   # RF_NUM_TREES, RF_THREADS
suppressPackageStartupMessages({library(ranger); library(dplyr)})
set.seed(42)
source("code/lib/geometry.R")
source("code/lib/species_levels.R")

SANDBOX   <- Sys.getenv("TREE_PRED_SANDBOX")
MODEL_DIR <- if (nzchar(SANDBOX)) SANDBOX else "outputs/models"
TAB_DIR   <- if (nzchar(SANDBOX)) SANDBOX else "outputs/tables"
load(file.path(MODEL_DIR, "RF_MODELS.RData"))
load(file.path(MODEL_DIR, "TRAINING_DATA.RData"))
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
DR <- rf_workflow_data$DRIVERS
d  <- tree_train_complete
trained <- sort(unique(as.character(d$species_clean)))

INVFILE <- "outputs/tables/inventory_stems.csv"
if (!file.exists(INVFILE))
  stop("missing ", INVFILE, " -- run: Rscript code/01_import/08_inventory_build.R")
INV <- read.csv(INVFILE, stringsAsFactors = FALSE)

H_LO <- 50; H_HI <- 200    # cm, the forest's training range in height

INV <- INV %>%
  filter(is.finite(dbh_m), dbh_m > 0) %>%
  mutate(sp = species_to_model_level(species, trained)) %>%
  as.data.frame()
INV$band_m <- band_length_m(INV$species)   # geometry.R: 2 m; Kalmia 0.75 m
INV$A_stem_m2 <- pi * INV$dbh_m * INV$band_m

cat(sprintf("stems %d | in stand %d | located %d\n",
            nrow(INV), sum(INV$in_stand), sum(INV$located)))
report_species_levels(INV$species, INV$sp, trained)
# Stems with no trained level are predicted under every level and averaged, weighted
# by each level's share of the in-stand 0-2 m band area (code/lib/species_levels.R).
LEVEL_W <- species_level_weights(INV$sp[INV$in_stand], INV$A_stem_m2[INV$in_stand], trained)
NOLEV <- INV$sp == NO_LEVEL
cat(sprintf("  %d stems with no level: band-area weights over %d levels (largest %s %.3f)\n",
            sum(NOLEV), sum(LEVEL_W > 0), names(which.max(LEVEL_W)), max(LEVEL_W)))

# --- moisture: per stem, from the same surface and climatology as the soil term --
# Each stem takes its local relative wetness from the moisture surface (vwc_rel, mean 1
# over the stand) scaled by the month's climatology level, as the soil grid does.
# Unlocated stems take vwc_rel = 1, the stand mean (counts of stems outside the survey
# hull: outputs/data/moisture_surface_coverage.csv).
GRIDFILE <- "outputs/tables/moisture_surface_grid.csv"
CLIMFILE <- "outputs/tables/moisture_climatology_monthly.csv"
for (f in c(GRIDFILE, CLIMFILE))
  if (!file.exists(f)) stop("missing ", f,
    " -- run 04_moisture_surface.R and 02_moisture_climatology.R first")
S    <- read.csv(GRIDFILE, stringsAsFactors = FALSE)
CLIM <- read.csv(CLIMFILE, stringsAsFactors = FALSE)
stopifnot(nrow(CLIM) == 12, all(is.finite(CLIM$moisture)))

# nearest grid cell for each stem; the surface is on a regular PX/PY lattice
gx <- sort(unique(S$PX)); gy <- sort(unique(S$PY))
snap <- function(v, g) g[pmax(1, pmin(length(g), round((v - g[1])/diff(g)[1]) + 1))]
key  <- paste(S$PX, S$PY)
idx  <- match(paste(snap(INV$PX, gx), snap(INV$PY, gy)), key)
INV$vwc_rel <- ifelse(is.na(idx), 1, S$vwc_rel[idx])
INV$vwc_rel[!INV$located | !is.finite(INV$vwc_rel)] <- 1
cat(sprintf("per-stem moisture: %d of %d stems matched a surface cell (rel. %.2f-%.2f)\n",
            sum(!is.na(idx) & INV$located), nrow(INV),
            min(INV$vwc_rel), max(INV$vwc_rel)))

# MOIST[stem, month] = local relative wetness x that month's plot-wide level
MOIST <- outer(INV$vwc_rel, CLIM$moisture[order(CLIM$mo)])
cat(sprintf("monthly stand-mean moisture %.3f-%.3f (climatology, 2019-2022)\n",
            min(colMeans(MOIST)), max(colMeans(MOIST))))

# --- soil temperature: the monthly climatology (00_soil_temp_climatology.R) -----
STFILE <- "outputs/tables/soil_temp_climatology_monthly.csv"
if (!file.exists(STFILE))
  stop("missing ", STFILE, " -- run code/04_drivers/00_soil_temp_climatology.R first")
STCLIM <- read.csv(STFILE, stringsAsFactors = FALSE)
stopifnot(nrow(STCLIM) == 12, all(is.finite(STCLIM$soil_temp_C)))
DR$soil_temp_C_mean <- STCLIM$soil_temp_C[match(DR$month, STCLIM$mo)]

# --- per-species calibration --------------------------------------------------
# A regression forest shrinks predictions toward the conditional mean (low-flux
# species high, high-flux species low), which is right for prediction and wrong for
# a sum. Each species level's prediction is multiplied by observed / out-of-bag
# predicted mean flux of its live training rows; levels with fewer than MIN_CAL_N
# rows keep ratio 1 (code/lib/species_calibration.R). Evaluated under cross-validation
# in 05_model/12_model_family_comparison.R and 08_rf_calibration_sensitivity.R.
source("code/lib/species_calibration.R")   # min sample MIN_CAL_N; see that file
.live <- !d$dead_stem   # calibrate on live stems
cal <- species_calibration(d$species_clean[.live], d$stem_flux_corrected[.live], TreeRF$predictions[.live])
cmap <- setNames(cal$ratio, cal$sp)
INV$cal <- as.numeric(ifelse(is.na(cmap[INV$sp]), 1, cmap[INV$sp]))
cat("\nper-species calibration (observed / out-of-bag predicted):\n")
print(as.data.frame(cal %>% arrange(ratio) %>%
      transmute(sp, n, obs_mean = round(obs_mean,4), oob_mean = round(oob_mean,4),
                ratio = round(ratio,3))), row.names = FALSE)
cat(sprintf("  ratio range %.3f - %.3f across %d levels; %d level(s) with n < %d set to 1\n",
            min(cal$ratio), max(cal$ratio), nrow(cal), sum(cal$n < MIN_CAL_N), MIN_CAL_N))

# Predictions for stems with a level are uncalibrated here (calibrated below through
# INV$cal); a stem with no level gets the weighted mean of the calibrated level
# predictions (its INV$cal is 1), or of the uncalibrated ones with mult = NULL.
pred_at <- function(h, mo, idx = seq_len(nrow(INV)), mult = cmap) predict_species_marginal(
  function(nd) predict(TreeRF, nd, num.threads = RF_THREADS)$predictions,
  data.frame(dbh_m = INV$dbh_m[idx],
             soil_moisture_at_tree = MOIST[idx, mo],
             soil_temp_C = DR$soil_temp_C_mean[mo],       # predictor <- monthly climatology
             air_temp_C = DR$air_temp_C_mean[mo],         # not a predictor of the stem forest; used by the with-air sensitivity
             height_cm = h,
             dead_stem = 0),   # inventory is live stems; ignored unless a variant uses it
  INV$sp[idx], LEVEL_W, mult)

# --- locate the steps, rather than assuming where they are --------------------
# Split thresholds are a property of the forest, so a sample of stems reveals all
# of them. Scanned at 1 cm once, on one month; the result is the exact interval
# structure of the model in height.
set.seed(42)
samp <- sample(nrow(INV), min(400, nrow(INV)))
Hs <- H_LO:H_HI
Ps <- vapply(Hs, function(h) pred_at(h, 7, samp), numeric(length(samp)))
moved <- colSums(abs(Ps[, -1, drop = FALSE] - Ps[, -ncol(Ps), drop = FALSE])) > 0
brk <- Hs[-1][moved]
edges <- c(H_LO, brk, H_HI + 1)          # value is constant on [edges[k], edges[k+1])
K <- length(edges) - 1
cat(sprintf("height steps detected at %s cm -> %d constant intervals\n",
            paste(brk, collapse = ", "), K))

# --- exact band integral ------------------------------------------------------
# Lateral area of a cylinder accumulates uniformly with height (dA/dz = pi*D is
# constant), so the flux-weighted band mean is the height-mean of f times the band
# area. Below H_LO the model is out of training range and the H_LO value is held.
# For Kalmia the band is 0.75 m, entirely below the first step.
band_len_cm <- INV$band_m * 100
band_of <- function(F) {
  band <- matrix(0, nrow(INV), 12)
  for (mo in 1:12) {
    acc <- pmin(band_len_cm, H_LO) * F[, 1, mo]          # 0 .. 50 cm held
    for (k in seq_len(K)) {
      w <- pmax(0, pmin(band_len_cm, edges[k+1]) - max(edges[k], H_LO))
      acc <- acc + w * F[, k, mo]
    }
    band[, mo] <- acc / band_len_cm
  }
  band
}
# ---- END OF SETUP: code/06_upscale/11_sens_daily_drivers.R and 05_model/15_rf_air_temperature_test.R
# ---- evaluate the lines above (inventory, moisture, calibration, pred_at, edges, band_of)
# ---- and drive the same forest with other drivers. Nothing above writes a file.

# --- evaluate every stem once per interval, per month -------------------------
cat(sprintf("predicting %d stems x %d intervals x 12 months (%d calls, not %d)...\n",
            nrow(INV), K, K*12, length(Hs)*12))
F <- array(NA_real_, c(nrow(INV), K, 12))
for (mo in 1:12) for (k in seq_len(K)) F[, k, mo] <- pred_at(edges[k], mo)

band <- band_of(F)
# the same without calibration, for the report below
F_uncal <- F
if (any(NOLEV)) for (mo in 1:12) for (k in seq_len(K))
  F_uncal[NOLEV, k, mo] <- pred_at(edges[k], mo, which(NOLEV), mult = NULL)
# the 2 m anchor for the above-band scenarios: the interval containing H_HI
flux_2m <- rowMeans(F[, K, ])

band_uncal <- band_of(F_uncal)
band    <- band    * INV$cal
flux_2m <- flux_2m * INV$cal

# --- export the calibrated band profile for 07_scaling_full_grid.R -------------
# The grid reads each stem's within-band profile (for the per-stem decay rate of
# exp_band_slope) and its 2 m anchor from here, so both terms of the whole-surface
# total share drivers and calibration.
PROF <- apply(F, c(1, 2), mean) * INV$cal        # n x K, annual mean, CALIBRATED
stopifnot(nrow(PROF) == nrow(INV), ncol(PROF) == K)
# 2 m anchor is the last interval of the profile, so the grid cannot disagree
stopifnot(max(abs(PROF[, K] - flux_2m)) < 1e-12)
PROFOUT <- cbind(
  data.frame(source = INV$source, tag = INV$tag, in_stand = INV$in_stand),
  setNames(as.data.frame(PROF), sprintf("f_i%d", seq_len(K))))
write.csv(PROFOUT, file.path(TAB_DIR, "tree_band_profile.csv"), row.names = FALSE)
write.csv(data.frame(k = seq_len(K), edge_lo_cm = head(edges, -1),
                     edge_hi_cm = tail(edges, -1)),
          file.path(TAB_DIR, "tree_band_profile_edges.csv"), row.names = FALSE)
cat(sprintf("exported calibrated band profile: %d stems x %d intervals (edges %s)\n",
            nrow(PROF), K, paste(edges, collapse = ", ")))

# lon/lat alongside PX/PY, so map panels need no transform of their own
ll <- geo_transforms()$fwd(INV$PX, INV$PY)
OUT <- data.frame(
  source = INV$source, tag = INV$tag, species = INV$species, species_model = INV$sp,
  dbh_m = INV$dbh_m, PX = INV$PX, PY = INV$PY,
  x = ifelse(INV$located, ll$lon, NA_real_), y = ifelse(INV$located, ll$lat, NA_real_),
  located = INV$located,
  in_stand = INV$in_stand, band_m = INV$band_m, A_stem_m2 = INV$A_stem_m2,
  flux_band_nmol_m2_s = rowMeans(band),
  flux_2m_nmol_m2_s   = flux_2m)
OUT$flux_nmol_m2_s <- OUT$flux_band_nmol_m2_s
write.csv(OUT, file.path(TAB_DIR, "tree_flux_predictions.csv"), row.names = FALSE)

# Stand monthly series (the only place it is computed; Figure 9 and the canonical
# budget read it):
#   tree_nmol_m2_s     the 0-2 m band's contribution per m2 of ground: sum(flux x stem
#                      area) over the in-stand stems / stand area (a budget quantity)
#   tree_bh_nmol_m2_s  the area-weighted mean flux at breast height per m2 of bark,
#                      comparable to measured stem fluxes (Figure 9c)
keep <- OUT$in_stand
BH_CM <- 130                                   # breast height, 1.3 m
k_bh <- findInterval(BH_CM, edges, rightmost.closed = TRUE)
bh <- F[, k_bh, ] * INV$cal                    # stems x months, calibrated
w  <- INV$A_stem_m2[keep]
tree_monthly <- data.frame(
  month = 1:12,
  tree_nmol_m2_s    = colSums(band[keep, ] * INV$A_stem_m2[keep]) / STAND_AREA_M2,
  tree_bh_nmol_m2_s = colSums(bh[keep, ] * w) / sum(w))
cat(sprintf("breast height falls in interval %d [%d, %d) cm; mean flux there %.4f nmol m-2 s-1\n",
            k_bh, edges[k_bh], edges[k_bh+1], mean(tree_monthly$tree_bh_nmol_m2_s)))
write.csv(tree_monthly, file.path(TAB_DIR, "tree_monthly_stand.csv"), row.names = FALSE)

CONV <- 86400 * 365.25 * 16e-6
S <- OUT[OUT$in_stand, ]
cat(sprintf("
  stems written  %d  (budget set in_stand = %d)
  band flux      median %.4f   range %.4f - %.4f nmol m-2 s-1
  band budget    %.3f mg CH4 m-2 yr-1  over %.2f ha of censused ground
  uncalibrated   %.3f mg CH4 m-2 yr-1  (%+.1f%% before per-species calibration)
",
  nrow(OUT), nrow(S), median(OUT$flux_band_nmol_m2_s),
  min(OUT$flux_band_nmol_m2_s), max(OUT$flux_band_nmol_m2_s),
  sum(S$flux_band_nmol_m2_s * S$A_stem_m2)/STAND_AREA_M2 * CONV, STAND_AREA_M2/1e4,
  sum(rowMeans(band_uncal)[keep] * INV$A_stem_m2[keep])/STAND_AREA_M2 * CONV,
  100*(sum(rowMeans(band_uncal)[keep]*INV$A_stem_m2[keep]) /
       sum(rowMeans(band)[keep]*INV$A_stem_m2[keep]) - 1)))
