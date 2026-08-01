source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# rev_budget_canonical.R  --  THE canonical budget.
# ------------------------------------------------------------------------------
# Every downstream number is computed here and written out; nothing is hardcoded
# in a figure script. This script now ASSEMBLES rather than recomputes: each
# quantity has exactly one producer, and the budget reads it.
#
#   inventory      rev_inventory_build.R        -> inventory_stems.csv
#   tree, per stem rev_predict_tree_flux_current.R -> tree_flux_predictions.csv
#   tree, monthly  rev_predict_tree_flux_current.R -> tree_monthly_stand.csv
#   soil, surface  rev_predict_soil_surface.R   -> soil_surface_{monthly,annual}.csv
#   geometry       rev_geometry.R               -> STAND_AREA_M2
#
# WHAT CHANGED AND WHY
#   * This script used to recompute the tree predictions itself AND write
#     outputs/tables/tree_flux_predictions.csv with a different stem set from the
#     script that also wrote that path, so whichever ran last won. It no longer
#     writes that file at all.
#   * The soil term came from MONTHLY_FLUXES.csv while Figure 9's map drew the
#     interpolated surface -- two estimates of one quantity, 4.7% apart. The
#     surface is now authoritative everywhere (decision 2026-07-27).
#   * Ground area was hardcoded as 200^2. The censused stand is 37,200 m2: the
#     square minus the 7 uncensused 20 m quadrats. See rev_geometry.R, which is
#     the only place that number is defined. (This comment said 37,600 -- a typo
#     that never reached the code, since STAND_AREA_M2 is sourced.)
#   * The measured band was reported as 5.41 here and 5.79 in the scaling grid,
#     and Figure 9 drew one stacked on the other. The grid now consumes
#     `tree_measured_mg_m2_yr` from this file instead of recomputing it.
#
# Output: outputs/revision/canonical_budget.csv   (one row per quantity)
#         outputs/revision/canonical_monthly.csv  (monthly series, both terms)
# ==============================================================================
suppressPackageStartupMessages({library(dplyr)})
source("code/lib/rev_geometry.R")
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

need <- c(inv  = "outputs/tables/inventory_stems.csv",
          tree = "outputs/tables/tree_flux_predictions.csv",
          trm  = "outputs/tables/tree_monthly_stand.csv",
          soim = "outputs/tables/soil_surface_monthly.csv",
          soia = "outputs/tables/soil_surface_annual.csv")
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing inputs:\n  ", paste(miss, collapse = "\n  "),
  "\nrun, in order: rev_inventory_build.R, rev_predict_tree_flux_current.R, rev_predict_soil_surface.R")

INV <- read.csv(need[["inv"]],  stringsAsFactors = FALSE)
TR  <- read.csv(need[["tree"]], stringsAsFactors = FALSE)
TRM <- read.csv(need[["trm"]],  stringsAsFactors = FALSE)
SOM <- read.csv(need[["soim"]], stringsAsFactors = FALSE)
SOA <- read.csv(need[["soia"]], stringsAsFactors = FALSE)
load("outputs/models/RF_MODELS.RData")

SEC_PER_DAY <- 86400; SEC_YR <- 86400*365.25
# 365.25/12, not 30.4: the latter annualises to 364.8 days and put the tree term
# 0.12% below the same quantity computed directly from seconds.
DAYS_PER_MONTH <- 365.25/12
NMOL_TO_MG  <- 16e-6

B_inv  <- INV[INV$in_stand, ]
B_tree <- TR[TR$in_stand, ]
stopifnot(nrow(B_inv) == nrow(B_tree))

basal_area <- sum(B_inv$basal_area_m2)
stem_area  <- sum(B_tree$A_stem_m2)
A_soil     <- STAND_AREA_M2 - basal_area

cat("GEOMETRY\n")
cat(sprintf("  stand %d m2 (%.2f ha) = %d x %d minus notch\n",
            STAND_AREA_M2, STAND_AREA_M2/1e4, PLOT_SIDE_M, PLOT_SIDE_M))
cat(sprintf("  stems %d (%d unlocated, retained) | stem area 0-2 m %.0f m2 (index %.4f)\n",
            nrow(B_inv), sum(!B_inv$located), stem_area, stem_area/STAND_AREA_M2))
cat(sprintf("  basal area %.1f m2 (%.1f m2 ha-1) | soil area %.0f m2 (%.3f of stand)\n",
            basal_area, basal_area/(STAND_AREA_M2/1e4), A_soil, A_soil/STAND_AREA_M2))

# ---- monthly series, both terms per m2 GROUND -------------------------------
soil_monthly <- SOM %>% group_by(month) %>%
                summarise(soil_nmol_m2_s = mean(flux_nmol_m2_s), .groups = "drop")
M <- data.frame(month = 1:12) %>%
  left_join(TRM, "month") %>% left_join(soil_monthly, "month") %>%
  mutate(tree_bh_nmol_m2_s = tree_bh_nmol_m2_s,
         tree_mg_m2_d = tree_nmol_m2_s*SEC_PER_DAY*NMOL_TO_MG,
         # x soil fraction, exactly as the annual term below. Without it this series
         # was per m2 of SOIL while tree_mg_m2_d was per m2 of GROUND, so plot_mg_m2_d
         # summed two different bases -- the very mixing the annual fix removed. It is
         # only 0.45%, but rev_fig09_budget.R plots this series in panel (c) under the
         # label "Soil (per sq.m ground)", which was then untrue.
         soil_mg_m2_d = soil_nmol_m2_s*SEC_PER_DAY*NMOL_TO_MG*(A_soil/STAND_AREA_M2),
         plot_mg_m2_d = tree_mg_m2_d + soil_mg_m2_d)
stopifnot(all(is.finite(M$tree_nmol_m2_s)), all(is.finite(M$soil_nmol_m2_s)))
write.csv(M, out_path("canonical_monthly.csv"), row.names = FALSE)

tree_ann <- sum(M$tree_mg_m2_d)*DAYS_PER_MONTH
# BOTH TERMS ON THE SAME BASIS: per m2 of GROUND.
# The tree term is a sum over stems divided by STAND_AREA_M2, so it is per m2 of ground
# already. The soil term is the spatial MEAN of the predicted surface, i.e. per m2 of
# SOIL, and soil is not the whole stand -- basal area occupies 166.77 m2 of it. Summing
# the two without rescaling adds a per-soil quantity to a per-ground one, which is
# exactly the basis mixing this paper has to be careful about. Multiply by the soil
# fraction (0.99552) to put it on ground. The effect is small, 0.448%, but the repo's
# own legacy code did this (02_rf_models.R: mean_flux * A_soil/A_plot) and
# canonical_budget.csv already exports soil_area_fraction -- so a consumer could
# reasonably have applied it a second time.
soil_ann <- mean(SOA$flux_nmol_m2_s)*SEC_YR*NMOL_TO_MG * (A_soil/STAND_AREA_M2)
tree_ann_2m <- sum(B_tree$flux_2m_nmol_m2_s*B_tree$A_stem_m2)/STAND_AREA_M2*SEC_YR*NMOL_TO_MG

# Grouped-CV skill alongside OOB. OOB is not grouped by tree, and the 1,130 tree rows
# come from only 478 trees (one contributes 13), so OOB scores the model partly on
# trees it has already seen -- it overstates TreeRF by ~2x. Carried here so that any
# script or paragraph quoting model skill has both numbers in front of it and cannot
# reach for the flattering one by default. NA if rev_rf_grouped_cv.R has not run.
GCV <- local({
  f <- "outputs/data/rf_grouped_cv.csv"
  if (!file.exists(f)) { cat("NOTE: rf_grouped_cv.csv absent -- run rev_rf_grouped_cv.R\n")
                         return(c(tree = NA_real_, soil = NA_real_)) }
  g <- read.csv(f, stringsAsFactors = FALSE)
  c(tree = g$grouped_r2[g$model == "TreeRF"][1],
    soil = g$grouped_r2[g$model == "SoilRF"][1])
})

B <- data.frame(
  quantity = c("stand_area_m2","plot_nominal_area_m2","notch_area_m2",
               "stem_area_0_2m_m2","stem_area_index","basal_area_m2","basal_area_m2_ha",
               "soil_area_m2","soil_area_fraction","n_stems","n_stems_unlocated",
               "tree_measured_mg_m2_yr","tree_at_2m_basis_mg_m2_yr","soil_mg_m2_yr",
               "net_measured_mg_m2_yr","tree_pct_of_soil","tree_r2_oob","soil_r2_oob",
               "tree_r2_grouped_cv","soil_r2_grouped_cv"),
  value = c(STAND_AREA_M2, PLOT_AREA_M2, NOTCH_AREA_M2,
            stem_area, stem_area/STAND_AREA_M2, basal_area, basal_area/(STAND_AREA_M2/1e4),
            A_soil, A_soil/STAND_AREA_M2, nrow(B_inv), sum(!B_inv$located),
            tree_ann, tree_ann_2m, soil_ann,
            tree_ann + soil_ann, 100*abs(tree_ann)/abs(soil_ann),
            TreeRF$r.squared, SoilRF$r.squared,
            GCV[["tree"]], GCV[["soil"]]))
write.csv(B, out_path("canonical_budget.csv"), row.names = FALSE)

cat("\nCANONICAL BUDGET\n")
for (i in seq_len(nrow(B))) cat(sprintf("  %-28s %14.4f\n", B$quantity[i], B$value[i]))
cat(sprintf("\n  measured band offsets %.2f%% of soil uptake; net %.1f mg CH4 m-2 yr-1\n",
            100*abs(tree_ann)/abs(soil_ann), tree_ann + soil_ann))
cat("\n  Written: canonical_budget.csv, canonical_monthly.csv\n")
