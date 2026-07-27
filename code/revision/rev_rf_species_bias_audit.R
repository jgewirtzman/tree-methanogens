#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_species_bias_audit.R
# ------------------------------------------------------------------------------
# For every measured species: how far is the model's prediction from what we
# actually measured, and if it is far, how much of the tree budget does that
# species carry?
#
# THE COMPARISON HAS TO BE MADE TWICE, because the obvious one is misleading.
#
#   MATCHED (the real bias test). Each training measurement has an out-of-bag
#   prediction made at ITS OWN height, month, moisture and temperature. Comparing
#   those to the observed values isolates model error from everything else. This
#   is the number to judge the model on.
#
#   POPULATION (what the budget actually uses). The inventory prediction is a
#   0-2 m band integral averaged over 12 months, applied to every stem of that
#   species including the very small ones. It SHOULD differ from the measured
#   mean, for three reasons that are not model error:
#     - the band mean exceeds the 2 m value, because flux declines with height;
#     - measurements are summer-weighted while the budget is an annual mean;
#     - measured stems are far larger than the median inventory stem.
#   Reporting this without the matched comparison would attribute all three to
#   bias.
#
# The budget sensitivity then rescales each species by its MATCHED ratio -- the
# only defensible correction factor -- and reports what the stand tree term would
# be if every species were unbiased in that sense.
#
# Output: outputs/revision/species_bias_audit.csv / .txt
# ==============================================================================
suppressPackageStartupMessages({library(dplyr); library(ranger)})
load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")
source("code/revision/rev_geometry.R")

CONV <- 86400 * 365.25 * 16e-6

# ---- 1. MATCHED: out-of-bag predictions at the measured conditions -----------
d <- tree_train_complete
stopifnot(length(TreeRF$predictions) == nrow(d))
d$oob <- TreeRF$predictions
d <- d[is.finite(d$stem_flux_corrected) & is.finite(d$oob), ]
M <- d %>% group_by(species) %>%
  summarise(n_meas = dplyr::n(),
            obs_mean = mean(stem_flux_corrected), obs_med = median(stem_flux_corrected),
            oob_mean = mean(oob), oob_med = median(oob),
            matched_ratio = mean(oob)/mean(stem_flux_corrected),
            matched_bias  = mean(oob) - mean(stem_flux_corrected),
            rmse = sqrt(mean((oob - stem_flux_corrected)^2)),
            .groups = "drop") %>%
  filter(!is.na(species))

# ---- 2. POPULATION: what the inventory gets ---------------------------------
TP <- read.csv("outputs/tables/tree_flux_predictions.csv", stringsAsFactors = FALSE)
TP <- TP[TP$in_stand, ]
P <- TP %>% group_by(species) %>%
  summarise(n_stems = dplyr::n(),
            pred_mean = mean(flux_band_nmol_m2_s), pred_med = median(flux_band_nmol_m2_s),
            band_area_m2 = sum(A_stem_m2),
            contrib_mg = sum(flux_band_nmol_m2_s * A_stem_m2)/STAND_AREA_M2 * CONV,
            .groups = "drop")
TOT <- sum(P$contrib_mg)
P$pct_budget <- 100 * P$contrib_mg / TOT

A <- full_join(M, P, "species") %>%
  mutate(pop_ratio = pred_mean/obs_mean) %>%
  arrange(desc(contrib_mg))

cat("=== MATCHED: out-of-bag prediction vs observation, at the same conditions ===\n")
cat("    (this is the model-bias test; ratio 1.00 = unbiased)\n")
print(as.data.frame(A %>% filter(!is.na(n_meas)) %>%
  transmute(species, n_meas, obs_mean = round(obs_mean,4), oob_mean = round(oob_mean,4),
            ratio = round(matched_ratio,3), rmse = round(rmse,3)) %>%
  arrange(desc(abs(log(pmax(A$matched_ratio[!is.na(A$n_meas)], 1e-6)))))), row.names = FALSE)

cat("\n=== POPULATION: inventory prediction vs measured mean ===\n")
cat("    (expected to differ: band integral, annual mean, and much smaller stems)\n")
print(as.data.frame(A %>% filter(!is.na(n_stems)) %>%
  transmute(species, n_stems, band_area_m2 = round(band_area_m2),
            obs_mean = round(obs_mean,4), pred_mean = round(pred_mean,4),
            pop_ratio = round(pop_ratio,2),
            contrib_mg = round(contrib_mg,3), pct_budget = round(pct_budget,1))), row.names = FALSE)

cat(sprintf("\ntotal tree band budget: %.3f mg CH4 m-2 yr-1 over %.2f ha\n", TOT, STAND_AREA_M2/1e4))

# ---- 3. budget sensitivity to the matched bias ------------------------------
cat("\n=== IF EACH SPECIES WERE CORRECTED BY ITS MATCHED RATIO ===\n")
corr <- A %>% filter(!is.na(n_meas), !is.na(n_stems), is.finite(matched_ratio), matched_ratio > 0)
TP2 <- TP %>% left_join(corr %>% select(species, matched_ratio), "species") %>%
  mutate(f = ifelse(is.finite(matched_ratio) & matched_ratio > 0,
                    flux_band_nmol_m2_s / matched_ratio, flux_band_nmol_m2_s))
TOT2 <- sum(TP2$f * TP2$A_stem_m2)/STAND_AREA_M2 * CONV
cat(sprintf("  as modelled          %.3f mg CH4 m-2 yr-1\n", TOT))
cat(sprintf("  bias-corrected       %.3f mg CH4 m-2 yr-1  (%+.1f%%)\n", TOT2, 100*(TOT2-TOT)/TOT))

cat("\n  per-species effect of that correction, largest first:\n")
eff <- TP2 %>% group_by(species) %>%
  summarise(as_modelled = sum(flux_band_nmol_m2_s*A_stem_m2)/STAND_AREA_M2*CONV,
            corrected   = sum(f*A_stem_m2)/STAND_AREA_M2*CONV, .groups = "drop") %>%
  mutate(delta = corrected - as_modelled) %>% arrange(desc(abs(delta)))
print(as.data.frame(eff %>% transmute(species, as_modelled = round(as_modelled,3),
        corrected = round(corrected,3), delta = round(delta,3),
        pct_of_total = round(100*delta/TOT,1)) %>% head(10)), row.names = FALSE)

dir.create("outputs/revision", showWarnings = FALSE, recursive = TRUE)
write.csv(A, "outputs/revision/species_bias_audit.csv", row.names = FALSE)
con <- file("outputs/revision/species_bias_audit.txt", "w")
cat(sprintf("SPECIES BIAS AUDIT\nbuilt %s\ntree band budget %.3f mg CH4 m-2 yr-1; bias-corrected %.3f (%+.1f%%)\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), TOT, TOT2, 100*(TOT2-TOT)/TOT), file = con)
utils::write.table(as.data.frame(A %>% mutate(across(where(is.numeric), ~round(.x,4)))),
                   con, sep = "\t", row.names = FALSE, quote = FALSE)
close(con)
cat("\nwritten: outputs/revision/species_bias_audit.{csv,txt}\n")
