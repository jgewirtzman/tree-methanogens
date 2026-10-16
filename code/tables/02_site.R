#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 02_site.R -- data/tables/site.csv: the stand reference point and plot geometry
# ------------------------------------------------------------------------------
# One row. Plot-local PX/PY (m) map to WGS84 by a great-circle offset from the
# reference point after rotating the plot axes by rotation_rad (the geodetic
# transform stored with the model's .RData; code/lib/geometry.R). Stand area is
# the 37,200 m2 of censused ground (200 x 200 m minus seven uncensused 20 m
# quadrats), the denominator of the budget.
# Writes data/tables/site.csv
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/geometry.R")
OUT <- "data/tables"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

e <- new.env(); load("data/processed/integrated/rf_workflow_input_data_with_2023.RData", envir = e)
w <- e$rf_workflow_data
site <- tibble(
  site_id = "YMF", site_name = "Yale-Myers Forest ForestGEO plot, Connecticut, USA",
  ref_latitude = w$ref_lat, ref_longitude = w$ref_lon, rotation_rad = w$rotation_angle,
  earth_radius_m = w$earth_radius, plot_side_m = PLOT_SIDE_M, quadrat_side_m = GAP_CELL_M,
  plot_area_m2 = PLOT_AREA_M2, uncensused_area_m2 = GAP_AREA_M2, stand_area_m2 = STAND_AREA_M2,
  uncensused_quadrats = paste(sort(UNCENSUSED_QUADRATS), collapse = ";"))

# the budget is on the same basis
B <- "outputs/data/canonical_budget.csv"
if (file.exists(B)) {
  b <- read.csv(B); a <- b$value[b$quantity == "stand_area_m2"]
  if (length(a) && abs(a - site$stand_area_m2) > 1e-6) stop("stand area differs from canonical_budget.csv: ", a)
}
write_csv(site, file.path(OUT, "site.csv"))
cat(sprintf("site.csv: ref %.6f, %.5f; rotation %.6f rad; stand %d m2\n", site$ref_latitude,
            site$ref_longitude, site$rotation_rad, site$stand_area_m2))
