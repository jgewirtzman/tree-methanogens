#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 35_moisture_survey.R -- moisture_survey.csv (data/tables/)
# ------------------------------------------------------------------------------
# The spatial soil-moisture survey of 16 December 2020: one row per flagged point
# (79), from data/raw/inventory/spatial_data/soil_moisture_20201216.csv (handheld
# GPS waypoints, VWC by probe). point_id is the waypoint label.
# Points 788-796 carry one identical GPS position (the receiver did not update);
# shared_position marks them. The moisture surface (code/04_drivers/04_moisture_surface.R)
# keeps the first of them (788), so it is fitted to 71 points.
# GPS elevation (column Elevation, ~820-850, unit not recorded) is not carried;
# collars.csv and the lidar DEM give ground elevation.
#
# Writes data/tables/moisture_survey.csv
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
OUT <- "data/tables"
SRC <- "data/raw/inventory/spatial_data/soil_moisture_20201216.csv"
SURVEY_DATE <- as.Date("2020-12-16")

D <- read.csv(SRC, stringsAsFactors = FALSE)
stopifnot(nrow(D) == 79, !anyDuplicated(D$Label), all(is.finite(D$VWC)), all(is.finite(D$Latitude)))
pos <- paste(round(D$Longitude, 8), round(D$Latitude, 8))
MS <- tibble(point_id = as.character(D$Label), date = SURVEY_DATE, latitude = D$Latitude, longitude = D$Longitude,
             soil_vwc_pct = D$VWC, shared_position = pos %in% pos[duplicated(pos)])
stopifnot(sum(!duplicated(pos)) == 71, sum(MS$shared_position) == 9)
write_csv(MS, file.path(OUT, "moisture_survey.csv"), na = "")
cat(sprintf("moisture_survey: %d points, %d distinct positions, VWC %.1f-%.1f %%\n", nrow(MS), sum(!duplicated(pos)),
            min(MS$soil_vwc_pct), max(MS$soil_vwc_pct)))
