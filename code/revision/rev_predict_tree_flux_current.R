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
# MODEL: the locked TreeRF (species, dbh_m, dbh_within_z, soil moisture, soil
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
  group_by(sp) %>%
  mutate(dbh_within_z = if (n() > 1 && sd(dbh_m, na.rm = TRUE) > 0)
                          as.numeric(scale(dbh_m)) else 0) %>%
  ungroup() %>% as.data.frame()
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
    dbh_within_z = INV$dbh_within_z[idx], soil_moisture_at_tree = DR$m[mo],
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

OUT <- data.frame(
  source = INV$source, tag = INV$tag, species = INV$species, species_model = INV$sp,
  dbh_m = INV$dbh_m, PX = INV$PX, PY = INV$PY, located = INV$located,
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
tree_monthly <- data.frame(
  month = 1:12,
  tree_nmol_m2_s = colSums(band[keep, ] * INV$A_stem_m2[keep]) / STAND_AREA_M2)
write.csv(tree_monthly, "outputs/tables/tree_monthly_stand.csv", row.names = FALSE)

CONV <- 86400 * 365.25 * 16e-6
S <- OUT[OUT$in_stand, ]
cat(sprintf("
  stems written  %d  (budget set in_stand = %d)
  band flux      median %.4f   range %.4f - %.4f nmol m-2 s-1
  band budget    %.3f mg CH4 m-2 yr-1  over %.2f ha of censused ground
",
  nrow(OUT), nrow(S), median(OUT$flux_band_nmol_m2_s),
  min(OUT$flux_band_nmol_m2_s), max(OUT$flux_band_nmol_m2_s),
  sum(S$flux_band_nmol_m2_s * S$A_stem_m2)/STAND_AREA_M2 * CONV, STAND_AREA_M2/1e4))
