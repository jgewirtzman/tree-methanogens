#!/usr/bin/env Rscript
# ==============================================================================
# rev_predict_tree_flux_current.R
# ------------------------------------------------------------------------------
# Per-stem CH4 flux predictions for the Figure 9 map panel and the tree budget.
#
# INVENTORY: outputs/tables/inventory_stems.csv, built by rev_inventory_build.R
#   from the raw ForestGEO tables. Run that first.
#
#   THIS SCRIPT USED TO READ ITS OWN OUTPUT. It took the tree list from
#   outputs/tables/tree_flux_predictions.csv -- the file it then overwrote -- so
#   the list of mapped stems existed nowhere else in the repo, in a gitignored
#   directory, and rev_budget_canonical.R wrote the SAME path with a different
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
#   same value), so the vertical scenarios in rev_scaling_full_grid.R take over
#   there, applied as shape ratios to each stem's own 2 m value. That anchor and
#   the per-stem lateral area are written out here for the grid to consume.
#
# Output: outputs/tables/tree_flux_predictions.csv  (per-stem; map + budget + grid anchor)
#         outputs/tables/tree_monthly_stand.csv     (stand monthly series)
# ==============================================================================
suppressPackageStartupMessages({library(ranger); library(dplyr)})
set.seed(42)
source("code/revision/rev_geometry.R")

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
DR <- rf_workflow_data$PLACEHOLDER_DRIVERS
d  <- tree_train_complete
trained <- sort(unique(as.character(d$species_clean)))

INVFILE <- "outputs/tables/inventory_stems.csv"
if (!file.exists(INVFILE))
  stop("missing ", INVFILE, " -- run: Rscript code/revision/rev_inventory_build.R")
INV <- read.csv(INVFILE, stringsAsFactors = FALSE)

STEM_BAND_M   <- 2.00
KALMIA_BAND_M <- 0.75      # shrub; matches compute_geometry() in 02_rf_models.R
H_LO <- 50; H_HI <- 200    # cm, the RF's training range in height

INV <- INV %>%
  filter(is.finite(dbh_m), dbh_m > 0) %>%
  mutate(sp = ifelse(!is.na(species) & species %in% trained, species, "SPECIES_OTHER")) %>%
  as.data.frame()
INV$band_m <- ifelse(INV$species %in% "Kalmia latifolia", KALMIA_BAND_M, STEM_BAND_M)
INV$A_stem_m2 <- pi * INV$dbh_m * INV$band_m

cat(sprintf("stems %d | in stand %d | located %d | trained species %d | SPECIES_OTHER %d\n",
            nrow(INV), sum(INV$in_stand), sum(INV$located),
            sum(INV$sp != "SPECIES_OTHER"), sum(INV$sp == "SPECIES_OTHER")))

# monthly drivers, gap-filled the same way as the canonical budget
mm <- d %>% group_by(month) %>%
      summarise(m = mean(soil_moisture_at_tree, na.rm = TRUE), .groups = "drop")
sm <- soil_train_complete %>% group_by(month) %>%
      summarise(ms = mean(soil_moisture_at_site, na.rm = TRUE), .groups = "drop")
DR <- DR %>% left_join(mm, "month") %>% left_join(sm, "month") %>%
      mutate(m = ifelse(is.finite(m), m, ms))
fl <- function(v, mo) approx(mo[is.finite(v)], v[is.finite(v)], xout = mo, rule = 2)$y
DR$m <- fl(DR$m, DR$month); DR$soil_temp_C_mean <- fl(DR$soil_temp_C_mean, DR$month)

pred_at <- function(h, mo, idx = seq_len(nrow(INV))) predict(TreeRF, data.frame(
    species = factor(INV$sp[idx], levels = trained), dbh_m = INV$dbh_m[idx],
    soil_moisture_at_tree = DR$m[mo],
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
# rather than by predicted magnitude. See rev_model_family_comparison.R.
# NO CLAMP. An earlier version bounded the ratio to [0.2, 5] to guard against
# levels with 1-3 records. It bound on exactly one level (GEN_Betula, n = 2, raw
# ratio 0.181), so it changed almost nothing while introducing an arbitrary
# parameter that would have had to be defended. The unclamped ratio is what the
# cross-validation in rev_model_family_comparison.R actually evaluated.
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
