#!/usr/bin/env Rscript
# ==============================================================================
# 02_predict_soil_surface.R -- the soil CH4 flux surface (Fig 9a, budget)
# ------------------------------------------------------------------------------
# Stage C; run from the repository root. Monthly and annual flux per m2 of ground
# on a regular 2 m plot-local grid over the censused stand (the moisture-surface grid
# of 04_drivers/04_moisture_surface.R). Every soil number in the budget and in
# Figure 9 comes from here.
#
# Model: the soil forest (soil moisture, soil temperature, month).
# Drivers per cell and month: soil moisture = the cell's relative wetness on the
#   December 2020 surface (fixed in time) x the month's plot-wide level
#   (04_drivers/02_moisture_climatology.R, 2019-2022); soil temperature = the
#   monthly climatology (04_drivers/00_soil_temp_climatology.R). air_temp_C (the
#   tower 24-h monthly mean) is passed for the with-air sensitivity
#   (05_model/15_rf_air_temperature_test.R) and ignored by the soil forest.
# The soil term is the spatial mean of the surface over the stand.
#
# Outputs: outputs/tables/soil_surface_monthly.csv  (cell x month)
#          outputs/tables/soil_surface_annual.csv   (cell, annual mean)
#          outputs/audit/soil_surface_summary.txt
# ==============================================================================
source("code/lib/rf_spec.R")   # RF_NUM_TREES, RF_THREADS
suppressPackageStartupMessages({library(ranger); library(dplyr)})
set.seed(42)
source("code/lib/geometry.R")

load("outputs/models/RF_MODELS.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
DR  <- rf_workflow_data$DRIVERS

GRIDFILE <- "outputs/tables/moisture_surface_grid.csv"
if (!file.exists(GRIDFILE))
  stop("missing ", GRIDFILE, " -- run code/04_drivers/04_moisture_surface.R first")
S <- read.csv(GRIDFILE, stringsAsFactors = FALSE)
stopifnot(all(in_stand(S$PX, S$PY)))
cat(sprintf("moisture grid: %d cells on censused ground (%.2f ha)\n",
            nrow(S), nrow(S)*4/1e4))

# ---- monthly drivers ---------------------------------------------------------
# air temperature: the tower's 24-h calendar-month means, 2019-2022 (01_load_and_prep_data.R), all 12
DR <- DR %>% arrange(month)
stopifnot(nrow(DR) == 12, all(DR$month == 1:12), all(is.finite(DR$air_temp_C_mean)))
# soil temperature: the monthly climatology (00_soil_temp_climatology.R)
STFILE <- "outputs/tables/soil_temp_climatology_monthly.csv"
if (!file.exists(STFILE))
  stop("missing ", STFILE, " -- run code/04_drivers/00_soil_temp_climatology.R first")
STCLIM <- read.csv(STFILE, stringsAsFactors = FALSE)
stopifnot(nrow(STCLIM) == 12, all(is.finite(STCLIM$soil_temp_C)))
DR$soil_temp_C_mean <- STCLIM$soil_temp_C[match(DR$month, STCLIM$mo)]


# ---- monthly moisture: fixed spatial pattern x multi-year monthly level -------
# The pattern is the December 2020 survey surface, fixed in time; its level moves
# with the 2019-2022 climatology. Scaling is multiplicative, which keeps relative
# wetness and makes the surface mean equal the level. The soil collars set the level
# (through the climatology), not the pattern: they cover too little of the plot.
CLIMFILE <- "outputs/tables/moisture_climatology_monthly.csv"
if (!file.exists(CLIMFILE))
  stop("missing ", CLIMFILE, " -- run code/04_drivers/02_moisture_climatology.R first")
CLIM <- read.csv(CLIMFILE, stringsAsFactors = FALSE)
stopifnot(nrow(CLIM) == 12, all(is.finite(CLIM$moisture)))
shape <- S$vwc_rel                                    # already mean 1 over the stand
cat(sprintf("moisture: fixed pattern (rel. range %.2f-%.2f) x monthly level %.3f-%.3f\n",
            min(shape), max(shape), min(CLIM$moisture), max(CLIM$moisture)))

soil_predict <- function(moist, soil_temp, air_temp, month)    # one month (or day) of drivers
  predict(SoilRF, data.frame(soil_moisture_at_site = moist, soil_temp_C = soil_temp,
                             air_temp_C = air_temp, month = month), num.threads = RF_THREADS)$predictions

# ---- END OF SETUP: code/06_upscale/11_sens_daily_drivers.R and 05_model/15_rf_air_temperature_test.R
# ---- evaluate the lines above (grid, moisture pattern and level, drivers, soil_predict)
# ---- and drive the same forest with other drivers. Nothing above writes a file.

# ---- predict per cell per month ---------------------------------------------
NMOL_TO_MG <- 16e-6; SEC_YR <- 86400 * 365.25
out <- vector("list", 12)
for (m in 1:12) {
  moist <- shape * CLIM$moisture[CLIM$mo == m]      # pattern x this month's level
  out[[m]] <- data.frame(cell = seq_len(nrow(S)), PX = S$PX, PY = S$PY,
                         x = S$x, y = S$y, month = m,
                         soil_moisture_at_site = moist,
                         flux_nmol_m2_s = soil_predict(moist, DR$soil_temp_C_mean[DR$month == m][1],
                                                       DR$air_temp_C_mean[DR$month == m][1], m))
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

con <- file("outputs/audit/soil_surface_summary.txt", "w")
cat(sprintf("SOIL SURFACE (02_predict_soil_surface.R)\nbuilt %s\n\n",
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
