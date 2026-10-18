#!/usr/bin/env Rscript
# ==============================================================================
# 03_budget_canonical.R -- the canonical stand budget (Fig 9, Results 10)
# ------------------------------------------------------------------------------
# Stage C; run from the repository root. Assembles, does not recompute: each
# quantity has one producer and is read here.
#
#   inventory       01_import/08_inventory_build.R  -> outputs/tables/inventory_stems.csv
#   tree, per stem  01_predict_tree_flux.R          -> tree_flux_predictions.csv
#   tree, monthly   01_predict_tree_flux.R          -> tree_monthly_stand.csv
#   tree, profile   01_predict_tree_flux.R          -> tree_band_profile{,_edges}.csv
#   soil, surface   02_predict_soil_surface.R       -> soil_surface_{monthly,annual}.csv
#   ground area     code/lib/geometry.R             -> STAND_AREA_M2 (censused stand)
#
# Outputs: outputs/data/canonical_budget.csv   (one row per quantity)
#          outputs/data/canonical_monthly.csv  (monthly series, both terms)
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages({library(dplyr)})
source("code/lib/geometry.R")

need <- c(inv  = "outputs/tables/inventory_stems.csv",
          tree = "outputs/tables/tree_flux_predictions.csv",
          trm  = "outputs/tables/tree_monthly_stand.csv",
          prof = "outputs/tables/tree_band_profile.csv",
          edge = "outputs/tables/tree_band_profile_edges.csv",
          soim = "outputs/tables/soil_surface_monthly.csv",
          soia = "outputs/tables/soil_surface_annual.csv")
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing inputs:\n  ", paste(miss, collapse = "\n  "),
  "\nrun, in order: 08_inventory_build.R, 01_predict_tree_flux.R, 02_predict_soil_surface.R")

INV <- read.csv(need[["inv"]],  stringsAsFactors = FALSE)
TR  <- read.csv(need[["tree"]], stringsAsFactors = FALSE)
TRM <- read.csv(need[["trm"]],  stringsAsFactors = FALSE)
PRF <- read.csv(need[["prof"]], stringsAsFactors = FALSE)
EDG <- read.csv(need[["edge"]], stringsAsFactors = FALSE)
SOM <- read.csv(need[["soim"]], stringsAsFactors = FALSE)
SOA <- read.csv(need[["soia"]], stringsAsFactors = FALSE)
load("outputs/models/RF_MODELS.RData")

SEC_PER_DAY <- 86400; SEC_YR <- 86400*365.25
# 365.25/12, so twelve months make one year of seconds
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
  mutate(tree_mg_m2_d = tree_nmol_m2_s*SEC_PER_DAY*NMOL_TO_MG,
         # x soil fraction: per m2 of ground, as the annual term below
         soil_mg_m2_d = soil_nmol_m2_s*SEC_PER_DAY*NMOL_TO_MG*(A_soil/STAND_AREA_M2),
         plot_mg_m2_d = tree_mg_m2_d + soil_mg_m2_d)
stopifnot(all(is.finite(M$tree_nmol_m2_s)), all(is.finite(M$soil_nmol_m2_s)))
write.csv(M, out_path("canonical_monthly.csv"), row.names = FALSE)

tree_ann <- sum(M$tree_mg_m2_d)*DAYS_PER_MONTH
# Both terms per m2 of ground. The tree term is a sum over stems divided by the stand
# area. The soil term is the spatial mean of the surface, i.e. per m2 of soil, so it is
# multiplied by the soil fraction of the stand (basal area excluded).
soil_ann <- mean(SOA$flux_nmol_m2_s)*SEC_YR*NMOL_TO_MG * (A_soil/STAND_AREA_M2)
tree_ann_2m <- sum(B_tree$flux_2m_nmol_m2_s*B_tree$A_stem_m2)/STAND_AREA_M2*SEC_YR*NMOL_TO_MG
# Methods S6: the band without its height profile, each stem's breast-height (1.3 m)
# prediction held from 0 to 2 m (annual mean, calibrated; same rows as the per-stem file)
stopifnot(identical(PRF$source, TR$source), identical(PRF$tag, TR$tag), identical(PRF$in_stand, TR$in_stand))
k_bh <- which(EDG$edge_lo_cm <= 130 & EDG$edge_hi_cm > 130)
tree_ann_bh <- sum(PRF[[sprintf("f_i%d", k_bh)]][TR$in_stand]*B_tree$A_stem_m2)/STAND_AREA_M2*SEC_YR*NMOL_TO_MG

# Grouped-CV skill (whole trees or collars held out; 05_model/04_rf_grouped_cv.R)
# alongside OOB, which scores the stem forest partly on trees it has seen.
GCV <- local({
  f <- "outputs/data/rf_grouped_cv.csv"
  if (!file.exists(f)) stop("missing ", f, " -- run code/05_model/04_rf_grouped_cv.R", call. = FALSE)
  g <- read.csv(f, stringsAsFactors = FALSE)
  c(tree = g$grouped_r2[g$model == "TreeRF"][1],
    soil = g$grouped_r2[g$model == "SoilRF"][1])
})

B <- data.frame(
  quantity = c("stand_area_m2","plot_nominal_area_m2","notch_area_m2",
               "stem_area_0_2m_m2","stem_area_index","basal_area_m2","basal_area_m2_ha",
               "soil_area_m2","soil_area_fraction","n_stems","n_stems_unlocated",
               "tree_measured_mg_m2_yr","tree_at_2m_basis_mg_m2_yr","tree_at_bh_basis_mg_m2_yr","soil_mg_m2_yr",
               "net_measured_mg_m2_yr","tree_pct_of_soil","tree_r2_oob","soil_r2_oob",
               "tree_r2_grouped_cv","soil_r2_grouped_cv"),
  value = c(STAND_AREA_M2, PLOT_AREA_M2, NOTCH_AREA_M2,
            stem_area, stem_area/STAND_AREA_M2, basal_area, basal_area/(STAND_AREA_M2/1e4),
            A_soil, A_soil/STAND_AREA_M2, nrow(B_inv), sum(!B_inv$located),
            tree_ann, tree_ann_2m, tree_ann_bh, soil_ann,
            tree_ann + soil_ann, 100*abs(tree_ann)/abs(soil_ann),
            TreeRF$r.squared, SoilRF$r.squared,
            GCV[["tree"]], GCV[["soil"]]))
write.csv(B, out_path("canonical_budget.csv"), row.names = FALSE)

cat("\nCANONICAL BUDGET\n")
for (i in seq_len(nrow(B))) cat(sprintf("  %-28s %14.4f\n", B$quantity[i], B$value[i]))
cat(sprintf("\n  measured band offsets %.2f%% of soil uptake; net %.1f mg CH4 m-2 yr-1\n",
            100*abs(tree_ann)/abs(soil_ann), tree_ann + soil_ann))
cat("\n  Written: canonical_budget.csv, canonical_monthly.csv\n")
