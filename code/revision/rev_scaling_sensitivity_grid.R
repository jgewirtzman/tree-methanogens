#!/usr/bin/env Rscript
# ==============================================================================
# rev_scaling_sensitivity_grid.R
# ------------------------------------------------------------------------------
# Stand-scaled stem-CH4 budget and its sensitivity to the assumptions the data
# cannot constrain.
#
# THE FLUX SIDE IS THE RANDOM FOREST, APPLIED TO THE INVENTORY. For every one of the
# ~7,000 censused stems and every month, TreeRF predicts flux at a grid of heights
# through the measured band, using that tree's species, DBH, within-species DBH
# z-score, and the month's moisture/temperature. Per-tree profiles are integrated over
# height, multiplied by that tree's own stem area, summed over the stand, and divided
# by plot area. Nothing here uses a single stand-average flux profile.
#
# THE RF IS USED ONLY WHERE IT WAS TRAINED (0.5-2 m). A random forest returns a
# constant outside its training range, so an RF "extrapolation" above 2 m would
# silently be the 2 m value. Outside the measured band the parametric forms enter as
# SHAPE RATIOS ONLY -- r(z) = form(z)/form(2 m) -- applied to each tree's own
# RF-predicted 2 m value. Population-level shape, per-tree magnitude.
#
# THREE UNCONSTRAINED CHOICES, crossed:
#
#   A. BELOW 0.5 m  (nothing was measured below 50 cm)
#        floor        hold each tree's 0.5 m value        (conservative)
#        extrapolate  continue the fitted shape downward
#
#   B. ABOVE 2 m  (one tree was ever profiled higher)
#        at_2m        each tree's 2 m value over its whole surface
#        mean_obs     each tree's mean over the measured band, over its whole surface
#        constant     hold the 2 m value above 2 m
#        exp_asym / power / exp_zero / loglinear / linear0   decaying forms
#
#   C. HEIGHT OF THE WOODY SURFACE (m)
#        NO TREE HEIGHTS EXIST IN THIS DATASET -- the inventory carries tag, species,
#        DBH and coordinates only, and the felled oak is the sole tree with a measured
#        height. An earlier version of this script used a DBH-height allometry; that
#        equation was not fitted to these data and not taken from a source, so it has
#        been removed rather than presented as an independent estimate. Canopy height
#        is instead carried as an explicit, labelled sensitivity parameter.
#
#   GEOMETRY: total woody area = WAI x plot area (WAI = 3.07 m2 m-2 ground, Gauci et
#   al. 2024 -- a PER-GROUND-AREA index, not per tree), allocated across stems in
#   proportion to pi*DBH, and assumed uniformly distributed with height to the canopy.
#   Woody surface area was never measured at this site and its vertical distribution is
#   unknown; both are assumptions, not measurements.
#
# The felled oak (0.038 nmol m-2 s-1 at 2 m; 0.034/0.077/0.031 at 6/8/10 m, plus a 10x
# hotspot at 4 m) is the only empirical check above 2 m. It favours the flatter forms,
# but n = 1 tree in a dataset with large tree-to-tree variation grounds the range
# rather than selecting from it.
#
# Output: outputs/revision/scaling_sensitivity_grid.csv
#         outputs/revision/fig_scaling_sensitivity.png
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(ggplot2); library(tidyr)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load("outputs/models/RF_MODELS.RData")          # TreeRF
load("outputs/models/TRAINING_DATA.RData")      # tree_train_complete
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
INV     <- rf_workflow_data$PLACEHOLDER_INVENTORY
DRIVERS <- rf_workflow_data$PLACEHOLDER_DRIVERS

CONV     <- 86400 * 365.25 * 16e-6   # nmol m-2 s-1 -> mg CH4 m-2 yr-1
A_PLOT   <- 200 * 200                # censused core, m2
WAI      <- 3.07                     # m2 woody area per m2 GROUND (not per tree)
SOIL_ANN <- -543.88                  # mg CH4 m-2 yr-1
CANOPY_H <- c(15, 20, 25)            # m; unmeasured -- explicit sensitivity parameter

# --- inventory ---------------------------------------------------------------
fix_dbh <- function(d, s) {
  d <- ifelse(!is.na(d) & d > 3, d/100, d)
  d <- ifelse(grepl("Betula", s) & !is.na(d) & d*100 > 200, d/10, d)
  d <- ifelse(s == "Pinus strobus" & !is.na(d) & d*100 > 230, d/10, d)
  ifelse(s == "Kalmia latifolia" & !is.na(d) & d*100 > 100, d/100,
  ifelse(s == "Kalmia latifolia" & !is.na(d) & d*100 > 10,  d/10, d))
}
INV$dbh <- fix_dbh(INV$dbh_m, INV$species)
INV <- INV[is.finite(INV$dbh) & INV$dbh > 0, ]
INV <- INV %>% group_by(species) %>%
  mutate(dbh_within_z = if (n() > 1 && sd(dbh, na.rm = TRUE) > 0)
           as.numeric(scale(dbh)) else 0) %>% ungroup() %>% as.data.frame()

# per-tree areas. Kalmia is a shrub: 0.75 m of stem, not 2 m.
INV$band_h  <- ifelse(INV$species == "Kalmia latifolia", 0.75, 2)
INV$A_stem  <- pi * INV$dbh * INV$band_h
INV$A_woody <- (WAI * A_PLOT) * (pi * INV$dbh) / sum(pi * INV$dbh)   # allocate ∝ πD

trained <- levels(TreeRF$forest$levels %||% factor(tree_train_complete$species_clean))
if (is.null(trained)) trained <- sort(unique(as.character(tree_train_complete$species_clean)))
INV$sp_model <- ifelse(INV$species %in% trained, INV$species, NA_character_)
cat(sprintf("inventory: %d stems, %d species; %d stems (%.2f%% of stem area) unmeasured species\n",
            nrow(INV), dplyr::n_distinct(INV$species), sum(is.na(INV$sp_model)),
            100*sum(INV$A_stem[is.na(INV$sp_model)])/sum(INV$A_stem)))
cat(sprintf("stem area 0-2 m %.0f m2 (index %.4f) | woody area %.0f m2 (WAI %.2f)\n\n",
            sum(INV$A_stem), sum(INV$A_stem)/A_PLOT, WAI*A_PLOT, WAI))

# --- monthly drivers ---------------------------------------------------------
mm <- tree_train_complete %>% group_by(month) %>%
  summarise(moist = mean(soil_moisture_at_tree, na.rm = TRUE), .groups = "drop")
sm <- soil_train_complete %>% group_by(month) %>%
  summarise(moist_s = mean(soil_moisture_at_site, na.rm = TRUE), .groups = "drop")
DR <- DRIVERS %>% left_join(mm, "month") %>% left_join(sm, "month") %>%
  mutate(moist = ifelse(is.finite(moist), moist, moist_s))
fillm <- function(v, m) approx(m[is.finite(v)], v[is.finite(v)], xout = m, rule = 2)$y
DR$moist            <- fillm(DR$moist, DR$month)
DR$soil_temp_C_mean <- fillm(DR$soil_temp_C_mean, DR$month)

# --- per-tree, per-month, per-height RF predictions ---------------------------
HG <- seq(50, 200, by = 12.5)          # cm, within the RF's training range
predict_at <- function(sp_vec) {
  arr <- array(NA_real_, c(nrow(INV), length(HG), 12))
  for (mo in 1:12) for (j in seq_along(HG)) {
    nd <- data.frame(species = factor(sp_vec, levels = trained),
                     dbh_m = INV$dbh, dbh_within_z = INV$dbh_within_z,
                     soil_moisture_at_tree = DR$moist[mo],
                     soil_temp_C_mean = DR$soil_temp_C_mean[mo],
                     air_temp_C_mean  = DR$air_temp_C_mean[mo],
                     height_cm = HG[j])
    arr[, j, mo] <- predict(TreeRF, nd, num.threads = 1)$predictions
  }
  arr
}
cat("predicting", nrow(INV), "stems x", length(HG), "heights x 12 months ...\n")
# Species not represented in the flux surveys (Kalmia latifolia, Quercus velutina,
# Acer pensylvanicum) go to SPECIES_OTHER, the catch-all level the model was trained
# with for exactly this purpose -- not marginalised over the named species.
P <- predict_at(ifelse(is.na(INV$sp_model), "SPECIES_OTHER", INV$sp_model))
PROF <- apply(P, c(1, 2), mean)        # annual (calendar) mean per tree per height
cat("done. stand-mean profile over the measured band:\n")
print(round(setNames(colMeans(PROF), paste0(HG, "cm")), 4))

# --- shape ratios outside the measured band -----------------------------------
h <- HG/100; fbar <- colMeans(PROF)
lm_log <- lm(log(pmax(fbar,1e-9)) ~ h); lm_lin <- lm(fbar ~ h)
lm_pow <- lm(log(pmax(fbar,1e-9)) ~ log(h))
nls_as <- tryCatch(nls(fbar ~ c + a*exp(-b*h), start = list(c=min(fbar)*.8,a=max(fbar),b=1),
                       control = nls.control(warnOnly = TRUE)), error = function(e) NULL)
SHAPE <- list(
  constant  = function(z) rep(1, length(z)),
  power     = function(z) exp(coef(lm_pow)[2]*(log(pmax(z,.05)) - log(2))),
  exp_zero  = function(z) exp(coef(lm_log)[2]*(z - 2)),
  loglinear = function(z) exp(coef(lm_log)[2]*(z - 2)),
  linear0   = function(z) pmax(coef(lm_lin)[1] + coef(lm_lin)[2]*z, 0) /
                          max(coef(lm_lin)[1] + coef(lm_lin)[2]*2, 1e-9))
if (!is.null(nls_as)) { cf <- coef(nls_as)
  SHAPE$exp_asym <- function(z) (cf["c"] + cf["a"]*exp(-cf["b"]*z)) /
                                (cf["c"] + cf["a"]*exp(-cf["b"]*2)) }

f_2m   <- PROF[, length(HG)]                 # per-tree RF value at 2 m
f_05   <- PROF[, 1]                          # per-tree RF value at 0.5 m
f_band <- rowMeans(PROF)                     # per-tree mean over the measured band

# Below 0.5 m, per tree. This is an extrapolation of the MEASURED BAND trend and is
# deliberately independent of the above-2 m scenario: it is anchored on each tree's own
# 0.5 m value and uses the log-linear slope fitted across 0.5-2 m. (Anchoring it on the
# 2 m value instead would make "extrapolate" reduce near-ground flux below the measured
# 0.5 m value, which inverts the observed height decline.)
b_lo <- coef(lm_log)[2]
below_mean <- function(mode) {
  if (mode == "floor") return(f_05)          # hold the 0.5 m value: conservative
  z <- seq(0.01, 0.5, length.out = 40)
  f_05 * mean(exp(b_lo * (z - 0.5)))         # continue the band trend downward
}

res <- bind_rows(lapply(c("floor","extrapolate"), function(bm)
  bind_rows(lapply(c("at_2m","mean_obs", names(SHAPE)), function(am)
    bind_rows(lapply(CANOPY_H, function(H) {
      shp <- if (am %in% names(SHAPE)) SHAPE[[am]] else SHAPE$constant
      # --- lower stem: measured band (0.5-2 m, all RF heights) + below-0.5 treatment
      band <- (1.5/2)*rowMeans(PROF) + (0.5/2)*below_mean(bm)
      lower_mg <- sum(band * INV$A_stem) / A_PLOT * CONV
      # --- whole woody surface, per tree
      whole <- if (am == "at_2m")    f_2m
        else if (am == "mean_obs")   f_band
        else { zz <- seq(0.01, H, length.out = 300)
               rr <- shp(zz); rr[zz < 0.5 & bm == "floor"] <- shp(0.5)
               f_2m * mean(rr) }
      whole_mg <- sum(whole * INV$A_woody) / A_PLOT * CONV
      data.frame(below = bm, above = am, canopy_H_m = H,
                 lower_2m_mg = lower_mg, whole_surface_mg = whole_mg,
                 pct_of_soil = 100*whole_mg/abs(SOIL_ANN),
                 net_mg = SOIL_ANN + whole_mg)
    }))))))

cat("\n=== STAND-SCALED SENSITIVITY GRID (mg CH4 m-2 yr-1) ===\n\n")
print(res %>% arrange(desc(whole_surface_mg)) %>%
        mutate(across(where(is.numeric), ~round(.x, 2))), row.names = FALSE)
write.csv(res, file.path(outdir, "scaling_sensitivity_grid.csv"), row.names = FALSE)

cat(sprintf("\n  lower stem (0-2 m, RF over inventory) : %.2f to %.2f mg CH4 m-2 yr-1\n",
            min(res$lower_2m_mg), max(res$lower_2m_mg)))
cat(sprintf("  whole woody surface, all scenarios    : %.1f to %.1f mg (%.1f%%-%.1f%% of soil uptake)\n",
            min(res$whole_surface_mg), max(res$whole_surface_mg),
            min(res$pct_of_soil), max(res$pct_of_soil)))
cat(sprintf("  spread %.1f-fold | net budget %.1f to %.1f (sink in %d/%d scenarios)\n",
            max(res$whole_surface_mg)/max(min(res$whole_surface_mg),1e-9),
            min(res$net_mg), max(res$net_mg), sum(res$net_mg < 0), nrow(res)))

p <- res %>% mutate(above = reorder(above, whole_surface_mg)) %>%
  ggplot(aes(above, whole_surface_mg, fill = below)) +
  geom_col(position = position_dodge()) + coord_flip() +
  facet_wrap(~ paste0("canopy height ", canopy_H_m, " m")) +
  scale_fill_manual(values = c(extrapolate = "#B2182B", floor = "#2166AC"),
                    name = "below 0.5 m") +
  labs(title = expression("Stand-scaled tree CH"[4]*" budget under alternative assumptions"),
       subtitle = paste("RF predictions over all censused stems; forms applied as shape ratios only above 2 m.",
                        "\nWoody surface area and its vertical distribution were never measured at this site."),
       x = "treatment above 2 m", y = expression(mg~CH[4]~m^-2~yr^-1)) +
  theme_minimal(base_size = 9)
ggsave(file.path(outdir, "fig_scaling_sensitivity.png"), p, width = 11, height = 5, dpi = 200, bg = "white")
cat("\nWritten: outputs/revision/fig_scaling_sensitivity.png\n")
