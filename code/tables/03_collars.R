#!/usr/bin/env Rscript
# ==============================================================================
# 03_collars.R -- data/tables/collars.csv: one row per soil collar
# ------------------------------------------------------------------------------
# collar_id is the form the soil flux data use, <Plot letter>_<Plot Tag>
# (WS_1-1); plots.csv labels the same collars WS-1-1. Positions are the handheld
# GPS fixes in plots.csv. Its Elevation column is not metres (566-786 across
# collars whose lidar relief is 20 m), so elevation_m is read from the 1 m 3DEP lidar
# DEM at the GPS position. in_stand follows code/lib/soil_collars.R (the two
# wetland-margin collars outside the censused stand); the GPS position mapped to
# plot coordinates is reported alongside as a check.
# Writes data/tables/collars.csv
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/geometry.R"); source("code/lib/soil_collars.R")
OUT <- "data/tables"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

P <- read.csv("data/raw/inventory/spatial_data/plots.csv", stringsAsFactors = FALSE) %>%
  transmute(collar_id = sub("^([A-Z]+)-", "\\1_", Label), landscape_position = sub("_.*", "", collar_id),
            site = Site, latitude = Latitude, longitude = Longitude)
stopifnot(!anyDuplicated(P$collar_id))

# every collar in the soil flux data has a position, and every position a flux record
S <- read.csv("data/raw/field_data/ipad_data/Cleaned data/soilflux_total.csv", check.names = FALSE)
used <- unique(paste0(trimws(S$`Plot letter`), "_", trimws(S$`Plot Tag`)))
used <- used[!grepl("_$|^_|NA", used)]
if (length(setdiff(used, P$collar_id))) stop("soil flux collars without a position: ", paste(setdiff(used, P$collar_id), collapse = ", "))
if (length(setdiff(P$collar_id, used))) stop("positions without soil flux: ", paste(setdiff(P$collar_id, used), collapse = ", "))

dem <- terra::rast("data/raw/inventory/spatial_data/lidar_3dep_1m_utm18n.tif")
v <- terra::project(terra::vect(P, geom = c("longitude", "latitude"), crs = "EPSG:4326"), terra::crs(dem))
P$elevation_m <- round(terra::extract(dem, v)[, 2], 2)

xy <- geo_transforms()$inv(P$latitude, P$longitude)
P$in_stand <- !P$collar_id %in% OUT_OF_STAND_COLLARS
P$gps_in_stand <- in_stand(xy$PX, xy$PY)

out <- P %>% select(collar_id, landscape_position, latitude, longitude, elevation_m, in_stand)
write_csv(out, file.path(OUT, "collars.csv"))
cat(sprintf("collars.csv: %d collars (%s); %d outside the stand (%s)\n", nrow(out),
            paste(names(table(out$landscape_position)), table(out$landscape_position), collapse = " "),
            sum(!out$in_stand), paste(out$collar_id[!out$in_stand], collapse = ", ")))
cat(sprintf("  elevation %.1f-%.1f m (lidar DEM)\n", min(out$elevation_m), max(out$elevation_m)))
d <- P %>% filter(in_stand != gps_in_stand)
if (nrow(d)) cat(sprintf("  GPS position disagrees with soil_collars.R for: %s (PX, PY = %s)\n",
                         paste(d$collar_id, collapse = ", "),
                         paste(sprintf("%.0f,%.0f", xy$PX[P$in_stand != P$gps_in_stand],
                                       xy$PY[P$in_stand != P$gps_in_stand]), collapse = "; ")))
