#!/usr/bin/env Rscript
# ==============================================================================
# predict_soil_surface.R
# ------------------------------------------------------------------------------
# THE soil flux surface. Monthly and annual, per m2 of ground, clipped to the
# censused stand. Every soil number in the budget and in Figure 9 comes from here.
#
# WHY THIS SCRIPT EXISTS
#   The stand soil flux was being taken from two different places at once:
#   Figure 9's map panel drew the interpolated surface while its budget panel and
#   budget_canonical.R took the soil term from MONTHLY_FLUXES.csv. Those are
#   independent estimates of one quantity and they disagreed by 4.7% (-519.6 vs
#   -543.9 mg CH4 m-2 yr-1). The decision (2026-07-27) is that the interpolated
#   surface is authoritative, because it is the thing that can be clipped to a
#   polygon, and that every panel draws on it.
#
#   The surface also needed regenerating regardless: outputs/models/
#   soil_monthly_predictions.RData dates from 2025-09 and predates the model
#   lock, and code/06_upscale/06_predict_soil_flux.R still builds the old feature
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
# GRID. A regular 2 m plot-local grid over the censused stand, built by
#   moisture_surface.R with a thin-plate spline over the 71 December survey
#   points. This replaces an akima linear interpolation, which fitted flat
#   triangles between points -- the planar facets and creases visible in the old
#   Figure 9 soil panel -- and which was undefined outside the convex hull,
#   leaving 29% of the plot unfillable. See moisture_interpolation.R for the
#   seven-method comparison; note that akima spline, the obvious alternative, was
#   the worst method tested.
#
# Output: outputs/tables/soil_surface_monthly.csv  (cell x month)
#         outputs/tables/soil_surface_annual.csv   (cell, annual mean)
#         outputs/audit/soil_surface_summary.txt
# ==============================================================================
suppressPackageStartupMessages({library(ranger); library(dplyr)})
set.seed(42)
source("code/lib/geometry.R")

load("outputs/models/RF_MODELS.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
DR  <- rf_workflow_data$PLACEHOLDER_DRIVERS
AFF <- read.csv("outputs/tables/MOISTURE_AFFINE_TABLE.csv")

GRIDFILE <- "outputs/tables/moisture_surface_grid.csv"
if (!file.exists(GRIDFILE))
  stop("missing ", GRIDFILE, " -- run code/04_drivers/moisture_surface.R first")
S <- read.csv(GRIDFILE, stringsAsFactors = FALSE)
stopifnot(all(in_stand(S$PX, S$PY)))
cat(sprintf("moisture grid: %d cells on censused ground (%.2f ha)\n",
            nrow(S), nrow(S)*4/1e4))

# ---- monthly drivers ---------------------------------------------------------
fl <- function(v, mo) approx(mo[is.finite(v)], v[is.finite(v)], xout = mo, rule = 2)$y
DR <- DR %>% arrange(month) %>%
  mutate(air_temp_C_mean = fl(air_temp_C_mean, month))
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


# ---- monthly moisture: fixed spatial pattern x multi-year monthly level -------
# The plot-wide monthly LEVEL now comes from moisture_climatology.R, which
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
  stop("missing ", CLIMFILE, " -- run code/04_drivers/moisture_climatology.R first")
CLIM <- read.csv(CLIMFILE, stringsAsFactors = FALSE)
stopifnot(nrow(CLIM) == 12, all(is.finite(CLIM$moisture)))
shape <- S$vwc_rel                                    # already mean 1 over the stand
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

# Domain sensitivity is no longer informative: the pattern is normalised to mean
# 1 over the stand and the level is imposed by the climatology, so the soil mean
# is set by construction rather than emerging from the choice of domain. Clipping
# still governs the map and the tree term.

con <- file("outputs/audit/soil_surface_summary.txt", "w")
cat(sprintf("SOIL SURFACE (predict_soil_surface.R)\nbuilt %s\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = con)
cat(sprintf("grid cells (in stand): %d\nstand area: %d m2 (%.2f ha)\n",
            nrow(S), STAND_AREA_M2, STAND_AREA_M2/1e4), file = con)
cat(sprintf("annual soil flux: %.3f mg CH4 m-2 yr-1 (%.5f nmol m-2 s-1)\n",
            soil_annual_mg, mean(A$flux_nmol_m2_s)), file = con)
cat("\nmonthly (nmol m-2 s-1):\n", file = con)
for (i in seq_len(nrow(mon)))
  cat(sprintf("  %2d  %8.4f\n", mon$month[i], mon$soil_nmol_m2_s[i]), file = con)
close(con)
cat("\nwritten: soil_surface_monthly.csv, soil_surface_annual.csv, soil_surface_summary.txt\n")
