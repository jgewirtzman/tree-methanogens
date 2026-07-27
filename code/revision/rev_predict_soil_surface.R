#!/usr/bin/env Rscript
# ==============================================================================
# rev_predict_soil_surface.R
# ------------------------------------------------------------------------------
# THE soil flux surface. Monthly and annual, per m2 of ground, clipped to the
# censused stand. Every soil number in the budget and in Figure 9 comes from here.
#
# WHY THIS SCRIPT EXISTS
#   The stand soil flux was being taken from two different places at once:
#   Figure 9's map panel drew the interpolated surface while its budget panel and
#   rev_budget_canonical.R took the soil term from MONTHLY_FLUXES.csv. Those are
#   independent estimates of one quantity and they disagreed by 4.7% (-519.6 vs
#   -543.9 mg CH4 m-2 yr-1). The decision (2026-07-27) is that the interpolated
#   surface is authoritative, because it is the thing that can be clipped to a
#   polygon, and that every panel draws on it.
#
#   The surface also needed regenerating regardless: outputs/models/
#   soil_monthly_predictions.RData dates from 2025-09 and predates the model
#   lock, and code/04_scaling/03_predict_soil_flux.R still builds the old feature
#   set (moisture x temperature interactions, month_sin/month_cos) that the
#   locked SoilRF does not use. Its four predictors are soil_moisture_at_site,
#   soil_temp_C_mean, air_temp_C_mean and month.
#
# CLIPPING IS NOT RENORMALISATION. The soil term is the spatial MEAN of a
#   predicted surface, so restricting the domain does not rescale it the way the
#   tree term rescales -- it changes only because the mean over the new domain
#   differs. The surface as previously averaged ran from PX -41 to +220 and PY -5
#   to +206, so 22% of its cells lay outside the census square entirely, and that
#   outer margin has weaker uptake than the plot interior; including it was
#   diluting the sink.
#
# GRID. The existing surface grid is reused verbatim (its x/y and its VWC), so
#   the map's shape and resolution are unchanged and only the model and the clip
#   differ. Moisture on that grid is VWC in PERCENT; the affine table converts it
#   to the fraction the model was trained on. Note the raster in the workflow
#   RData is on a different scale (already a fraction) -- do not mix them.
#
# Output: outputs/tables/soil_surface_monthly.csv  (cell x month)
#         outputs/tables/soil_surface_annual.csv   (cell, annual mean)
#         outputs/revision/soil_surface_summary.txt
# ==============================================================================
suppressPackageStartupMessages({library(ranger); library(dplyr)})
set.seed(42)
source("code/revision/rev_geometry.R")

load("outputs/models/RF_MODELS.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
DR  <- rf_workflow_data$PLACEHOLDER_DRIVERS
AFF <- read.csv("outputs/tables/MOISTURE_AFFINE_TABLE.csv")

GRIDFILE <- "outputs/tables/soil_flux_extended_annual.csv"
if (!file.exists(GRIDFILE)) stop("missing ", GRIDFILE, " (soil surface grid)")
G <- read.csv(GRIDFILE, stringsAsFactors = FALSE)

# ---- put the grid in plot-local metres and clip to censused ground -----------
tr <- geo_transforms()
p <- tr$inv(G$y, G$x); G$PX <- p$PX; G$PY <- p$PY
G$in_stand <- in_stand(G$PX, G$PY)
cat(sprintf("grid cells %d | inside 200 m square %d | IN STAND %d (%.1f%% kept)\n",
            nrow(G), sum(in_square(G$PX, G$PY)), sum(G$in_stand), 100*mean(G$in_stand)))
S <- G[G$in_stand, ]

# ---- monthly drivers ---------------------------------------------------------
fl <- function(v, mo) approx(mo[is.finite(v)], v[is.finite(v)], xout = mo, rule = 2)$y
DR <- DR %>% arrange(month) %>%
  mutate(soil_temp_C_mean = fl(soil_temp_C_mean, month),
         air_temp_C_mean  = fl(air_temp_C_mean, month))

# ---- monthly moisture: fixed spatial pattern x multi-year monthly level -------
# The plot-wide monthly LEVEL now comes from rev_moisture_climatology.R, which
# averages 2019-2022 rather than tying the seasonal cycle to the single campaign
# year -- a year that had the second-driest Jun-Sep in the record (210 mm against
# 437 mm in 2021). Since the SoilRF learned a moisture x temperature interaction,
# driving it from that one year baked a drought into every month.
#
# The PATTERN still comes from the dense December transect and is treated as
# fixed in time; only its level moves. Scaling is multiplicative, which preserves
# relative wetness and makes the surface match the target mean by construction.
# The 28 soil collars are deliberately NOT used to build the pattern: they span
# PY 32-148 of 200 m, the median plot cell is 31 m from the nearest, and only 29%
# of the plot is within 20 m of one. They set the level (through the climatology's
# de-confounding step), not the shape.
#
# This replaces a per-month affine rescaling whose alpha and beta were fitted
# separately per month on 3-28 observations at R2 = 0.02-0.31, and which was not a
# climatology at all: it jumped from 0.113 in November to 0.355 in December.
CLIMFILE <- "outputs/tables/moisture_climatology_monthly.csv"
if (!file.exists(CLIMFILE))
  stop("missing ", CLIMFILE, " -- run code/revision/rev_moisture_climatology.R first")
CLIM <- read.csv(CLIMFILE, stringsAsFactors = FALSE)
stopifnot(nrow(CLIM) == 12, all(is.finite(CLIM$moisture)))
shape <- S$mean_moisture / mean(S$mean_moisture)     # relative pattern, mean 1
cat(sprintf("moisture: fixed pattern (rel. range %.2f-%.2f) x monthly level %.3f-%.3f\n",
            min(shape), max(shape), min(CLIM$moisture), max(CLIM$moisture)))

# ---- predict per cell per month ---------------------------------------------
NMOL_TO_MG <- 16e-6; SEC_YR <- 86400 * 365.25
out <- vector("list", 12)
for (m in 1:12) {
  moist <- shape * CLIM$moisture[CLIM$mo == m]      # pattern x this month's level
  nd <- data.frame(soil_moisture_at_site = moist,
                   soil_temp_C_mean = DR$soil_temp_C_mean[DR$month == m][1],
                   air_temp_C_mean  = DR$air_temp_C_mean[DR$month == m][1],
                   month = m)
  out[[m]] <- data.frame(cell = seq_len(nrow(S)), PX = S$PX, PY = S$PY,
                         x = S$x, y = S$y, month = m,
                         soil_moisture_at_site = moist,
                         flux_nmol_m2_s = predict(SoilRF, nd, num.threads = 1)$predictions)
}
M <- bind_rows(out)
A <- M %>% group_by(cell, PX, PY, x, y) %>%
     summarise(flux_nmol_m2_s = mean(flux_nmol_m2_s), .groups = "drop") %>%
     mutate(flux_mg_m2_yr = flux_nmol_m2_s * SEC_YR * NMOL_TO_MG)

write.csv(M, "outputs/tables/soil_surface_monthly.csv", row.names = FALSE)
write.csv(A, "outputs/tables/soil_surface_annual.csv",  row.names = FALSE)

# ---- stand-level series: the spatial mean over censused ground ---------------
mon <- M %>% group_by(month) %>%
       summarise(soil_nmol_m2_s = mean(flux_nmol_m2_s), .groups = "drop")
soil_annual_mg <- mean(A$flux_nmol_m2_s) * SEC_YR * NMOL_TO_MG

cat("\n=== STAND SOIL FLUX (spatial mean over censused ground) ===\n")
cat(sprintf("annual  %.2f mg CH4 m-2 yr-1   (%.4f nmol m-2 s-1)\n",
            soil_annual_mg, mean(A$flux_nmol_m2_s)))
print(as.data.frame(mon %>% mutate(soil_nmol_m2_s = round(soil_nmol_m2_s, 4))))

cat("\n-- domain sensitivity, for the record --\n")
dom <- function(sel, lab) {
  gg <- G[sel, ]
  mm <- sapply(1:12, function(m) {
    sh <- gg$mean_moisture / mean(gg$mean_moisture)
    mean(predict(SoilRF, data.frame(
      soil_moisture_at_site = sh * CLIM$moisture[CLIM$mo == m],
      soil_temp_C_mean = DR$soil_temp_C_mean[DR$month == m][1],
      air_temp_C_mean  = DR$air_temp_C_mean[DR$month == m][1],
      month = m), num.threads = 1)$predictions)
  })
  cat(sprintf("  %-34s %6d cells  %9.2f mg m-2 yr-1\n", lab, nrow(gg),
              mean(mm)*SEC_YR*NMOL_TO_MG))
}
dom(rep(TRUE, nrow(G)), "whole surface (old behaviour)")
dom(in_square(G$PX, G$PY), "200 x 200 m square")
dom(G$in_stand, "censused stand (used)")

con <- file("outputs/revision/soil_surface_summary.txt", "w")
cat(sprintf("SOIL SURFACE (rev_predict_soil_surface.R)\nbuilt %s\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = con)
cat(sprintf("grid cells kept (in stand): %d of %d\nstand area: %d m2 (%.2f ha)\n",
            nrow(S), nrow(G), STAND_AREA_M2, STAND_AREA_M2/1e4), file = con)
cat(sprintf("annual soil flux: %.3f mg CH4 m-2 yr-1 (%.5f nmol m-2 s-1)\n",
            soil_annual_mg, mean(A$flux_nmol_m2_s)), file = con)
cat("\nmonthly (nmol m-2 s-1):\n", file = con)
for (i in seq_len(nrow(mon)))
  cat(sprintf("  %2d  %8.4f\n", mon$month[i], mon$soil_nmol_m2_s[i]), file = con)
close(con)
cat("\nwritten: soil_surface_monthly.csv, soil_surface_annual.csv, soil_surface_summary.txt\n")
