# ==============================================================================
# 03_prep_tree_auxfile.R -- goFlux auxiliary file for the monthly semi-rigid stem fluxes
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/02_flux/semirigid/).
#
# One row per tagged-stem deployment with a valid sampling time and geometry.
# UniqueID = Date_PlotTag_PlotLetter_Chamber (a sequence number is appended where a
# stem was measured twice with the same chamber on one day).
#
# Area and chamber volume: the deployment's geometry from 02_join_flux_geometry.R
#   (code/lib/stem_chamber_geometry.R), carried with its source.
# Vtot = chamber volume + TUBING_CM3 + ANALYZER_VOLUME_CM3 (code/lib/chamber_constants.R).
# Tcham: station hourly mean air temperature of the hour containing the start, in true
#   time (logger clock converted by code/lib/weather_time.R, the rule of
#   code/tables/06_weather.R), within 24 h; every deployment must have one (the build
#   stops otherwise). Pcham: 101.325 kPa.
#
# Inputs:  data/processed/flux/flux_with_geometry_fixed.csv (02_join_flux_geometry.R)
#          data/raw/weather/ymf_clean_sorted.csv
# Outputs (data/processed/flux/):
#          auxfile_goFlux_with_weather.csv  (read by 04_goflux_trees.R)
#          auxfile_goFlux_with_weather.txt  (tab-separated copy)
#          auxfile_goFlux_with_weather_formatted_datetime.csv  (+ start.time as text;
#          the base of the assembled monthly stem table)
# ==============================================================================

source("../../lib/project_library.R")   # the project package library (renv.lock), however this script is started
library(dplyr)
library(readr)
library(lubridate)

source("../../lib/chamber_constants.R")  # TUBING_CM3, ANALYZER_VOLUME_CM3
source("../../lib/weather_time.R")       # read_station_air(), station_air_at()

flux_data <- read_csv("../../../data/processed/flux/flux_with_geometry_fixed.csv")

station_air <- read_station_air("../../../data/raw/weather/ymf_clean_sorted.csv")
cat("Station air temperature:", nrow(station_air), "hourly records\n")

prepared <- flux_data %>%
  mutate(
    numeric_time = as.numeric(`Time of sampling`),
    UniqueID_base = paste(Date, `Plot Tag`, `Plot Letter`, Chamber_ID, sep = "_"),
    UniqueID = ave(UniqueID_base, UniqueID_base, FUN = function(x) {
      if (length(x) == 1) x else paste(x, seq_along(x), sep = "_")
    }),
    # sampling time HHMM (1357 or 1357.0) -> HH:MM:00; missing or invalid times -> NA
    time_int = floor(numeric_time),
    hours = floor(time_int / 100),
    minutes = time_int %% 100,
    formatted_time = case_when(
      is.na(numeric_time) ~ NA_character_,
      hours >= 0 & hours <= 23 & minutes >= 0 & minutes <= 59 ~
        sprintf("%02d:%02d:00", hours, minutes),
      TRUE ~ NA_character_
    ),
    datetime_string = case_when(
      !is.na(formatted_time) ~ paste(Date, formatted_time),
      TRUE ~ NA_character_
    ),
    start.time = as.POSIXct(datetime_string, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    Area = surface_area_cm2,                                         # cm2
    Vtot = (volume_cm3 + TUBING_CM3 + ANALYZER_VOLUME_CM3) / 1000    # L
  )

removed_rows <- prepared %>%
  filter(is.na(start.time) | is.na(Area) | is.na(Vtot)) %>%
  select(Date, `Plot Tag`, `Plot Letter`, `Time of sampling`, surface_area_cm2, volume_cm3)

auxfile <- prepared %>%
  filter(!is.na(start.time), !is.na(Area), !is.na(Vtot)) %>%
  mutate(Tcham = station_air_at(start.time, station_air), Pcham = 101.325)
stopifnot(!anyNA(auxfile$Tcham))   # station air temperature within 24 h for every deployment

auxfile <- auxfile %>%
  select(UniqueID, start.time, Area, Vtot, Tcham, Pcham,
         Date, `Plot Tag`, `Plot Letter`, Chamber_ID, `File name`, Notes,
         geometry_source, chamber_dims, D_stem_cm, chamber_volume_cm3 = volume_cm3) %>%
  arrange(start.time)

cat("Deployments:", nrow(flux_data), "->", nrow(auxfile), "in the auxfile | range:",
    as.character(range(auxfile$start.time, na.rm = TRUE)), "\n")
cat("Tcham range:", round(range(auxfile$Tcham, na.rm = TRUE), 1), "degC | duplicate UniqueIDs:",
    sum(duplicated(auxfile$UniqueID)), "\n")
if (nrow(removed_rows) > 0) {
  cat("\nNot in the auxfile (no valid time or no geometry):\n")
  print(removed_rows)
}

write.table(auxfile, "../../../data/processed/flux/auxfile_goFlux_with_weather.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
write_csv(auxfile, "../../../data/processed/flux/auxfile_goFlux_with_weather.csv")
auxfile_readable <- auxfile %>%
  mutate(start.time_formatted = format(start.time, "%Y-%m-%d %H:%M:%S %Z"))
write_csv(auxfile_readable, "../../../data/processed/flux/auxfile_goFlux_with_weather_formatted_datetime.csv")
cat("Wrote data/processed/flux/auxfile_goFlux_with_weather{.csv,.txt,_formatted_datetime.csv}\n")
