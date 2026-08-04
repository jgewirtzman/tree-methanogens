#!/usr/bin/env Rscript
# ==============================================================================
# predict_tree_flux_current.R
# ------------------------------------------------------------------------------
# Per-stem CH4 flux predictions for the Figure 9 map panel and the tree budget.
#
# INVENTORY: outputs/tables/inventory_stems.csv, built by inventory_build.R
#   from the raw ForestGEO tables. Run that first.
#
#   THIS SCRIPT USED TO READ ITS OWN OUTPUT. It took the tree list from
#   outputs/tables/tree_flux_predictions.csv -- the file it then overwrote -- so
#   the list of mapped stems existed nowhere else in the repo, in a gitignored
#   directory, and budget_canonical.R wrote the SAME path with a different
#   stem set. Running the two in the wrong order deleted 963 stems with no way
#   to recover them. Both now read the inventory and neither writes it.
#
# MODEL: the locked TreeRF (species, dbh_m, soil moisture, soil
#   temperature, air temperature, height).
#
# HEIGHT: the band integral is computed EXACTLY, not on a fine grid.
#   The RF was trained at three heights only (50/125/200 cm, n = 146/143/145), so
#   in height it is a step function: every stem takes 3-4 distinct values across
#   50-200 cm and the transitions are at shared thresholds. Those breakpoints are
#   detected here rather than assumed, then each stem is evaluated once per
#   interval and the intervals are width-weighted. This is exact -- a 1 cm grid
#   computes the same number 38x more slowly, and the 12.5 cm grid used earlier
#   was simply wrong wherever a step fell inside a bin.
#
#   The RF resolves the MEASURED band and nothing else. Flux above 2 m is not its
#   job: it cannot extrapolate (beyond the largest split every leaf returns the
#   same value), so the vertical scenarios in scaling_full_grid.R take over
#   there, applied as shape ratios to each stem's own 2 m value. That anchor and
#   the per-stem lateral area are written out here for the grid to consume.
#
# Output: outputs/tables/tree_flux_predictions.csv  (per-stem; map + budget + grid anchor)
#         outputs/tables/tree_monthly_stand.csv     (stand monthly series)
# ==============================================================================
suppressPackageStartupMessages({library(ranger); library(dplyr)})
set.seed(42)
source("code/lib/geometry.R")
source("code/lib/species_levels.R")

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
DR <- rf_workflow_data$PLACEHOLDER_DRIVERS
d  <- tree_train_complete
trained <- sort(unique(as.character(d$species_clean)))

INVFILE <- "outputs/tables/inventory_stems.csv"
if (!file.exists(INVFILE))
  stop("missing ", INVFILE, " -- run: Rscript code/01_import/inventory_build.R")
INV <- read.csv(INVFILE, stringsAsFactors = FALSE)

STEM_BAND_M   <- 2.00
KALMIA_BAND_M <- 0.75      # shrub; matches compute_geometry() in 02_rf_models.R
H_LO <- 50; H_HI <- 200    # cm, the RF's training range in height

INV <- INV %>%
  filter(is.finite(dbh_m), dbh_m > 0) %>%
  mutate(sp = species_to_model_level(species, trained)) %>%
  as.data.frame()
INV$band_m <- ifelse(INV$species %in% "Kalmia latifolia", KALMIA_BAND_M, STEM_BAND_M)
INV$A_stem_m2 <- pi * INV$dbh_m * INV$band_m

cat(sprintf("stems %d | in stand %d | located %d\n",
            nrow(INV), sum(INV$in_stand), sum(INV$located)))
report_species_levels(INV$species, INV$sp, trained)

# --- moisture: PER STEM, from the same surface and climatology as the soil term -
# This used to be one plot-mean scalar applied to all 8,006 stems, built from raw
# campaign monthly means of soil_moisture_at_tree. That was wrong twice over.
#
# (1) It disagreed with the soil term. predict_soil_surface.R drives the soil
#     flux with the December TPS pattern scaled by the 2019-2022 climatology, cell
#     by cell. The tree term used a different moisture history entirely -- campaign
#     means with 2.5x the seasonal amplitude (0.239 against 0.097), built from as
#     few as ~20 measurements a month, with December having no tree measurements at
#     all. Two terms plotted in one budget figure must share a moisture history.
#
# (2) A single scalar puts every stem on the same side of every split. The forest
#     is a step function in moisture, and with one value per month the whole stand
#     falls off the same cliff at once: at June temperatures, predicted flux jumps
#     from 0.0562 to 0.1204 between moisture 0.27 and 0.28, and the climatology's
#     June level is 0.279. That produced a spurious June maximum 1.6x July's --
#     an artifact of which side of one threshold a single number happened to land.
#
# The stand has a moisture FIELD, not a value. Each stem now takes its local
# relative wetness from the TPS surface (vwc_rel, mean 1 over the stand) and that
# is scaled by the monthly climatology level, exactly as the soil grid is. The
# 51 unlocated stems (0.6%) take vwc_rel = 1, the stand mean. Averaging over 8,006
# stems spread across the field crosses many thresholds instead of one, so the
# month-to-month curve reflects the climatology rather than a split location.
GRIDFILE <- "outputs/tables/moisture_surface_grid.csv"
CLIMFILE <- "outputs/tables/moisture_climatology_monthly.csv"
for (f in c(GRIDFILE, CLIMFILE))
  if (!file.exists(f)) stop("missing ", f,
    " -- run moisture_surface.R and moisture_climatology.R first")
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

# --- soil temperature: multi-year climatology, not campaign means -------------
# DRIVERS$soil_temp_C_mean is the mean of whenever somebody was in the field: 27
# to 624 observations a month from as few as ONE visit, with January and March
# absent entirely and linearly interpolated. soil_temp_climatology.R replaces
# it with a damped/lagged air-temperature model (21-day trailing mean, leave-one-
# date-out CV R2 0.951, RMSE 1.13 C, damping slope 0.683) applied to the full
# 2018-2023 tower record. Air temperature already came from that record.
STFILE <- "outputs/tables/soil_temp_climatology_monthly.csv"
if (!file.exists(STFILE))
  stop("missing ", STFILE, " -- run code/04_drivers/soil_temp_climatology.R first")
STCLIM <- read.csv(STFILE, stringsAsFactors = FALSE)
stopifnot(nrow(STCLIM) == 12, all(is.finite(STCLIM$soil_temp_C)))
DR$soil_temp_C_mean <- STCLIM$soil_temp_C[match(DR$month, STCLIM$mo)]

pred_at <- function(h, mo, idx = seq_len(nrow(INV))) predict(TreeRF, data.frame(
    species = factor(INV$sp[idx], levels = trained), dbh_m = INV$dbh_m[idx],
    soil_moisture_at_tree = MOIST[idx, mo],
    soil_temp_C_mean = DR$soil_temp_C_mean[mo],
    air_temp_C_mean = DR$air_temp_C_mean[mo],
    height_cm = h), num.threads = 1)$predictions

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

# --- evaluate every stem once per interval, per month -------------------------
cat(sprintf("predicting %d stems x %d intervals x 12 months (%d calls, not %d)...\n",
            nrow(INV), K, K*12, length(Hs)*12))
F <- array(NA_real_, c(nrow(INV), K, 12))
for (mo in 1:12) for (k in seq_len(K)) F[, k, mo] <- pred_at(edges[k], mo)

# --- exact band integral ------------------------------------------------------
# Lateral area of a cylinder accumulates uniformly with height (dA/dz = pi*D is
# constant), so the flux-weighted band mean is the height-mean of f times the band
# area. Below H_LO the model is out of training range and the H_LO value is held.
# For Kalmia the band is 0.75 m, which lies entirely below the first step, so this
# general form reproduces the old special case exactly.
band_len_cm <- INV$band_m * 100
band <- matrix(0, nrow(INV), 12)
for (mo in 1:12) {
  acc <- pmin(band_len_cm, H_LO) * F[, 1, mo]          # 0 .. 50 cm held
  for (k in seq_len(K)) {
    w <- pmax(0, pmin(band_len_cm, edges[k+1]) - max(edges[k], H_LO))
    acc <- acc + w * F[, k, mo]
  }
  band[, mo] <- acc / band_len_cm
}
# the 2 m anchor for the above-band scenarios: the interval containing H_HI
flux_2m <- rowMeans(F[, K, ])

# --- per-species calibration --------------------------------------------------
# A regression forest minimises squared error, and that shrinks predictions toward
# the conditional mean: low-flux species come out high and high-flux species low.
# Correct for prediction, wrong for a SUM. Tsuga and Pinus strobus are both
# low-flux and carry 54% of the band area between them, so uncorrected the stand
# tree term is ~20% high (area-weighted ratio 1.205 under cross-validation).
#
# The correction is the ratio of observed to OUT-OF-BAG predicted mean, per
# species level. Using OOB predictions rather than fitted ones is what keeps this
# from being circular. Under repeated 5-fold CV, where the correction is refitted
# on each training fold, it takes the area-weighted sum ratio from 1.205 to 0.992
# while retaining 91% of the forest's skill (R2 0.210 vs 0.230).
#
# It also transports: residual bias across diameter quartiles falls to 0.131 from
# the forest's 0.158, and across moisture quartiles to 0.090 from 0.212. That
# matters because the correction is estimated at the covariate distribution of the
# MEASUREMENTS and applied to inventory stems that are smaller and span all twelve
# months. Calibrating on the predicted value instead -- linear or isotonic -- does
# NOT work (sum ratio 1.11-1.13), because the shrinkage is structured by species
# rather than by predicted magnitude. See model_family_comparison.R.
# NO CLAMP. An earlier version bounded the ratio to [0.2, 5] to guard against
# levels with 1-3 records. On the current model it binds on NO level at all -- every
# ratio falls inside [0.2, 5], so rf_calibration_sensitivity.R measures the clamped
# variant as changing the stand total by exactly 0.00%. It introduced an arbitrary
# parameter that would have had to be defended. The unclamped ratio is what the
# cross-validation in model_family_comparison.R actually evaluated.
cal <- data.frame(sp = as.character(d$species_clean),
                  obs = d$stem_flux_corrected, oob = TreeRF$predictions) %>%
  filter(is.finite(obs), is.finite(oob)) %>%
  group_by(sp) %>%
  summarise(n = dplyr::n(), obs_mean = mean(obs), oob_mean = mean(oob),
            ratio = ifelse(oob_mean > 0, obs_mean/oob_mean, 1), .groups = "drop")
cmap <- setNames(cal$ratio, cal$sp)
INV$cal <- as.numeric(ifelse(is.na(cmap[INV$sp]), 1, cmap[INV$sp]))
cat("\nper-species calibration (observed / out-of-bag predicted):\n")
print(as.data.frame(cal %>% arrange(ratio) %>%
      transmute(sp, n, obs_mean = round(obs_mean,4), oob_mean = round(oob_mean,4),
                ratio = round(ratio,3))), row.names = FALSE)
cat(sprintf("  no clamp applied; ratio range %.3f - %.3f across %d levels\n",
            min(cal$ratio), max(cal$ratio), nrow(cal)))

band_uncal <- band
band    <- band    * INV$cal
flux_2m <- flux_2m * INV$cal

# --- export the CALIBRATED band profile for scaling_full_grid.R -----------
# The grid needs each stem's within-band profile (to fit the per-stem decay rate
# for exp_band_slope, and to report the band shape) and its 2 m anchor. It used to
# rebuild both by re-running the forest itself, and that re-derivation drifted:
#   - it drove the model from DRIVERS$soil_moisture_at_tree and
#     DRIVERS$soil_temp_C_mean, i.e. the campaign means this script replaces above
#     (one plot-mean value per month, December substituted from the soil collars,
#     January and March absent and filled by approx(rule = 2)), and
#   - it applied NO per-species calibration.
# So 74% of the headline scenario -- the extrapolated term, which hangs entirely
# off that anchor -- sat on a different basis from the 26% that came from this
# script's canonical band. A single reported total cannot have its two components
# on different footings, so the profile is exported here and consumed there.
PROF <- apply(F, c(1, 2), mean) * INV$cal        # n x K, annual mean, CALIBRATED
stopifnot(nrow(PROF) == nrow(INV), ncol(PROF) == K)
# 2 m anchor is the last interval of the profile, so the grid cannot disagree
stopifnot(max(abs(PROF[, K] - flux_2m)) < 1e-12)
PROFOUT <- cbind(
  data.frame(source = INV$source, tag = INV$tag, in_stand = INV$in_stand),
  setNames(as.data.frame(PROF), sprintf("f_i%d", seq_len(K))))
write.csv(PROFOUT, "outputs/tables/tree_band_profile.csv", row.names = FALSE)
write.csv(data.frame(k = seq_len(K), edge_lo_cm = head(edges, -1),
                     edge_hi_cm = tail(edges, -1)),
          "outputs/tables/tree_band_profile_edges.csv", row.names = FALSE)
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
write.csv(OUT, "outputs/tables/tree_flux_predictions.csv", row.names = FALSE)

# Stand-level monthly series, per m2 GROUND: sum(flux x stem area) over the
# budget set, divided by censused stand area. This is the tree counterpart of
# the soil monthly series and the only place it is computed -- Figure 9 and the
# canonical budget both read it rather than re-deriving it.
keep <- OUT$in_stand
# Two monthly series, because they answer different questions:
#   tree_nmol_m2_s  the 0-2 m band's contribution per m2 of GROUND -- a budget
#                   quantity, only meaningful once summed and divided by plot area
#   tree_bh_nmol_m2_s  the area-weighted mean flux at BREAST HEIGHT, per m2 of
#                   WOODY SURFACE -- a measured quantity, directly comparable to
#                   the map panel and to any other stem-flux study
# Figure 9's seasonal panel uses the second: a per-ground-area tree series exists
# only because we summed stem area and divided by plot area, and putting it beside
# soil invites reading it as like-for-like when it is diluted by that denominator.
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
write.csv(tree_monthly, "outputs/tables/tree_monthly_stand.csv", row.names = FALSE)

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
