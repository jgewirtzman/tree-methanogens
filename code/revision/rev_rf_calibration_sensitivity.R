#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_calibration_sensitivity.R
# ------------------------------------------------------------------------------
# How much of the stand tree term rests on calibration ratios estimated from very
# few measurements, and how wide is the resulting interval?
#
# THE EXPOSURE. Stand totals are calibrated per species level by the ratio of mean
# observed to mean OUT-OF-BAG predicted flux (rev_predict_tree_flux_current.R). The
# ratio is applied unclamped. Most levels are well determined, but after the species
# ladder was fixed at prediction time (2026-07-29) the GEN_Kalmia level -- estimated
# from THREE measurements on three shrubs -- is applied to ~1,900 mountain laurel
# stems. Several genus levels rest on 2-6 records.
#
# This does not re-fit anything. It takes the locked model's OOB predictions and the
# per-stem band fluxes exactly as the budget computes them, and re-weights only the
# calibration ratios, so the reported band is attributable to the ratios alone.
#
# Scenarios:
#   as_shipped        the unclamped ratios the budget uses
#   low_n_to_1        every level with n <= N_SMALL set to 1.0 (no correction)
#   jackknife_lo/hi   every small level at the extreme of its leave-one-out range,
#                     all moved together -- a deliberate worst case, not a CI
#   shrunk            ratio shrunk toward 1 by n/(n + K), the standard remedy
#   clamped           the [0.2, 5] clamp that an earlier version applied
#
# Output: outputs/revision/rf_calibration_sensitivity.csv / .txt
# ==============================================================================
suppressPackageStartupMessages({library(ranger); library(dplyr)})
source("code/revision/rev_geometry.R")
source("code/revision/rev_species_levels.R")

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")

N_SMALL  <- 27      # "few": everything at or below SPECIES_OTHER's own count
SHRINK_K <- 10      # ratio -> 1 + (ratio-1) * n/(n+K); K = the ladder's own threshold
CONV     <- 86400 * 365.25 * 16e-6

d <- tree_train_complete
cal <- data.frame(sp = as.character(d$species_clean),
                  obs = d$stem_flux_corrected, oob = TreeRF$predictions) %>%
  filter(is.finite(obs), is.finite(oob)) %>%
  group_by(sp) %>%
  summarise(n = dplyr::n(), obs_mean = mean(obs), oob_mean = mean(oob),
            ratio = ifelse(mean(oob) > 0, mean(obs)/mean(oob), 1), .groups = "drop") %>%
  as.data.frame()

# leave-one-measurement-out range of each ratio (OOB predictions held fixed)
jk <- t(sapply(cal$sp, function(s) {
  i <- which(as.character(d$species_clean) == s &
             is.finite(d$stem_flux_corrected) & is.finite(TreeRF$predictions))
  if (length(i) < 2) return(c(lo = NA_real_, hi = NA_real_))
  r <- sapply(seq_along(i), function(k) {
    j <- i[-k]; om <- mean(TreeRF$predictions[j])
    if (om > 0) mean(d$stem_flux_corrected[j])/om else 1 })
  c(lo = min(r), hi = max(r))
}))
cal$jk_lo <- jk[, "lo"]; cal$jk_hi <- jk[, "hi"]

# per-stem band flux with the calibration REMOVED, so ratios can be re-applied
TRP <- read.csv("outputs/tables/tree_flux_predictions.csv", stringsAsFactors = FALSE)
TRP <- TRP[TRP$in_stand, ]
rmap <- setNames(cal$ratio, cal$sp)
TRP$ratio_used <- as.numeric(ifelse(is.na(rmap[TRP$species_model]), 1, rmap[TRP$species_model]))
TRP$band_uncal  <- TRP$flux_band_nmol_m2_s / TRP$ratio_used

stand_total <- function(rv) {
  r <- as.numeric(ifelse(is.na(rv[TRP$species_model]), 1, rv[TRP$species_model]))
  sum(TRP$band_uncal * r * TRP$A_stem_m2) / STAND_AREA_M2 * CONV
}
small <- cal$n <= N_SMALL
mk <- function(f) { v <- setNames(cal$ratio, cal$sp); v[small] <- f(cal[small, ]); v }

SC <- data.frame(
  scenario = c("as_shipped","low_n_to_1","jackknife_lo","jackknife_hi","shrunk","clamped"),
  total_mg = c(
    stand_total(setNames(cal$ratio, cal$sp)),
    stand_total(mk(function(x) 1)),
    stand_total(mk(function(x) ifelse(is.na(x$jk_lo), x$ratio, x$jk_lo))),
    stand_total(mk(function(x) ifelse(is.na(x$jk_hi), x$ratio, x$jk_hi))),
    stand_total(setNames(1 + (cal$ratio - 1) * cal$n/(cal$n + SHRINK_K), cal$sp)),
    stand_total(setNames(pmax(pmin(cal$ratio, 5), 0.2), cal$sp))),
  stringsAsFactors = FALSE)
SC$pct_vs_shipped <- 100 * (SC$total_mg / SC$total_mg[1] - 1)

# contribution of each small level, as shipped
contrib <- sapply(cal$sp[small], function(s) {
  i <- TRP$species_model == s
  sum(TRP$flux_band_nmol_m2_s[i] * TRP$A_stem_m2[i]) / STAND_AREA_M2 * CONV })
SM <- data.frame(sp = cal$sp[small], n = cal$n[small],
                 ratio = cal$ratio[small], jk_lo = cal$jk_lo[small], jk_hi = cal$jk_hi[small],
                 stems = sapply(cal$sp[small], function(s) sum(TRP$species_model == s)),
                 mg = as.numeric(contrib),
                 pct_of_tree = 100 * as.numeric(contrib) / SC$total_mg[1],
                 stringsAsFactors = FALSE)
SM <- SM[order(-SM$mg), ]

write.csv(SC, "outputs/revision/rf_calibration_sensitivity.csv", row.names = FALSE)
con <- file("outputs/revision/rf_calibration_sensitivity.txt", "w")
wr <- function(...) cat(sprintf(...), file = con)
wr("CALIBRATION SENSITIVITY OF THE STAND TREE TERM\n")
wr("Ratios re-weighted only; nothing re-fitted. \"few\" = n <= %d measurements.\n\n", N_SMALL)
wr("%-14s %6s %8s %10s %10s %7s %9s %8s\n","level","n","stems","ratio","jk lo-hi","","mg","%% of tree")
for (i in seq_len(nrow(SM))) with(SM[i, ],
  wr("%-14s %6d %8d %10.3f %10.3f %7.3f %9.4f %8.2f\n", sp, n, stems, ratio, jk_lo, jk_hi, mg, pct_of_tree))
wr("\n  all levels with n <= %d together: %.4f mg = %.2f%% of the tree term\n",
   N_SMALL, sum(SM$mg), 100*sum(SM$mg)/SC$total_mg[1])
wr("\n%-16s %12s %12s\n","scenario","total mg","%% vs shipped")
for (i in seq_len(nrow(SC))) with(SC[i, ], wr("%-16s %12.4f %12+.2f\n", scenario, total_mg, pct_vs_shipped))
wr("\njackknife_lo/hi move EVERY small level to its extreme at once, so they bound what\n")
wr("the low-n ratios can do rather than describing a likely outcome.\n")
close(con)
cat(readLines("outputs/revision/rf_calibration_sensitivity.txt"), sep = "\n")
cat("\n  Written: outputs/revision/rf_calibration_sensitivity.{csv,txt}\n")
