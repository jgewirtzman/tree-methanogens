# ==============================================================================
# 03_prep_soil_auxfile.R -- goFlux auxiliary file for the monthly soil-collar fluxes
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/02_flux/semirigid/).
#
# One row per soil deployment on the field sheet that has a file name and a valid
# sampling time. UniqueID = date_plot-letter_file (e.g. 20200615_2Jan-U_13).
#
# Area: SOIL_COLLAR_AREA_CM2 (25.43 cm interior diameter). Vtot as written here uses the
# pre-correction chamber volume SOIL_CHAMBER_VOLUME_PRECORRECTION_L (17.75 L: cap 14.42 L
# + collar headspace 3.33 L) + tubing + analyzer. No fit uses this Vtot: every soil fit
# (04_goflux_soils.R, 05_fix_december.R, 03_merge/01_fix_soil_flux.R) sets Vtot to
# SOIL_VTOT_L (7.53 L collar headspace; code/lib/chamber_constants.R).
# Tcham: station hourly mean air temperature of the hour containing the start, in true
#   time (logger clock converted by code/lib/weather_time.R, the rule of
#   code/tables/06_weather.R), within 24 h; every deployment must have one (the build
#   stops otherwise). Pcham: 101.325 kPa.
#
# Inputs:  data/raw/field_data/ipad_data/Cleaned data/soilflux_total.csv
#          data/raw/weather/ymf_clean_sorted.csv
# Outputs (data/processed/flux/):
#          auxfile_goFlux_soilflux_with_weather.csv  (read by the soil fits)
#          auxfile_goFlux_soilflux_with_weather.txt  (tab-separated copy)
#          auxfile_goFlux_soilflux_with_weather_formatted.csv  (+ start.time as text;
#          the base of the assembled soil tables)
# ==============================================================================

library(dplyr)
library(readr)
library(lubridate)

source("../../lib/chamber_constants.R")  # SOIL_COLLAR_AREA_CM2, SOIL_CHAMBER_VOLUME_PRECORRECTION_L, TUBING_CM3, ANALYZER_VOLUME_CM3
source("../../lib/weather_time.R")       # read_station_air(), station_air_at()

flux_data <- read_csv("../../../data/raw/field_data/ipad_data/Cleaned data/soilflux_total.csv",
                      col_types = cols(
                        Date = col_character(),
                        Site = col_character(),
                        `Plot Tag` = col_character(),
                        `Plot letter` = col_character(),
                        `Person sampling` = col_character(),
                        `Time of sampling` = col_character(),
                        `File name` = col_character(),
                        Notes = col_character()
                      ))

station_air <- read_station_air("../../../data/raw/weather/ymf_clean_sorted.csv")
cat("Station air temperature:", nrow(station_air), "hourly records\n")

total_system_volume_cm3 <- SOIL_CHAMBER_VOLUME_PRECORRECTION_L * 1000 + TUBING_CM3 + ANALYZER_VOLUME_CM3
cat("Area:", SOIL_COLLAR_AREA_CM2, "cm2 | Vtot written:", round(total_system_volume_cm3 / 1000, 3),
    "L (pre-correction chamber volume; the fits use SOIL_VTOT_L =", round(SOIL_VTOT_L, 3), "L)\n")
cat("Field sheet rows:", nrow(flux_data), "| without a file name:", sum(is.na(flux_data$`File name`)), "\n")

flux_data <- flux_data %>%
  mutate(Date_parsed = mdy(Date),                       # M/D/YY
         Date_formatted = format(Date_parsed, "%Y-%m-%d"))

auxfile <- flux_data %>%
  filter(!is.na(`File name`)) %>%                       # no gas record without a file
  mutate(
    # sampling time HHMM (possibly with decimals) -> HH:MM:00; invalid times -> NA
    numeric_time = as.numeric(`Time of sampling`),
    time_int = floor(numeric_time),
    hours = floor(time_int / 100),
    minutes = time_int %% 100,
    formatted_time = case_when(
      hours >= 0 & hours <= 23 & minutes >= 0 & minutes <= 59 ~
        sprintf("%02d:%02d:00", hours, minutes),
      TRUE ~ NA_character_
    ),
    datetime_string = case_when(
      !is.na(formatted_time) ~ paste(Date_formatted, formatted_time),
      TRUE ~ NA_character_
    ),
    start.time = as.POSIXct(datetime_string, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    UniqueID = paste0(gsub("-", "", Date_formatted), "_",
                      `Plot Tag`, "-", `Plot letter`, "_",
                      `File name`),
    Area = SOIL_COLLAR_AREA_CM2,             # cm2
    Vtot = total_system_volume_cm3 / 1000    # L
  ) %>%
  filter(!is.na(start.time), !is.na(Area), !is.na(Vtot)) %>%
  mutate(Tcham = station_air_at(start.time, station_air), Pcham = 101.325)
stopifnot(!anyNA(auxfile$Tcham))   # station air temperature within 24 h for every deployment

auxfile <- auxfile %>%
  select(UniqueID, start.time, Area, Vtot, Tcham, Pcham,
         Date = Date_formatted, Site, `Plot Tag`, `Plot letter`,
         `Person sampling`, `File name`, Notes) %>%
  arrange(start.time)

cat("Deployments:", nrow(auxfile), "| dates:", length(unique(auxfile$Date)),
    "| range:", as.character(range(auxfile$start.time)), "\n")
cat("Duplicate UniqueIDs:", sum(duplicated(auxfile$UniqueID)), "\n")
cat("Rows not in the auxfile (no file name or no valid time):", nrow(flux_data) - nrow(auxfile), "\n")

write.table(auxfile, "../../../data/processed/flux/auxfile_goFlux_soilflux_with_weather.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
write_csv(auxfile, "../../../data/processed/flux/auxfile_goFlux_soilflux_with_weather.csv")
auxfile_readable <- auxfile %>%
  mutate(start.time_formatted = format(start.time, "%Y-%m-%d %H:%M:%S %Z"))
write_csv(auxfile_readable, "../../../data/processed/flux/auxfile_goFlux_soilflux_with_weather_formatted.csv")
cat("Wrote data/processed/flux/auxfile_goFlux_soilflux_with_weather{.csv,.txt,_formatted.csv}\n")
