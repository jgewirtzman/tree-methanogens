# ==============================================================================
# 01_prep_auxfile_2021.R -- goFlux auxiliary file for the 2021 multi-height
#                           static-chamber campaign
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/02_flux/static/). Runs first in
# stage A: 01_import/01_tree_id_consensus.R reads its flux_data_with_dimensions table.
#
# Each deployment on the field sheet gets its chamber's surface area (by chamber
# series: chambers 1-2 = B, 3-4 = C, 5-6 = D), the series' mean measured chamber
# volume, and the tubing and analyzer volumes (Vtot = chamber + tubing + analyzer;
# no drierite). additional_vol.csv lists 70 cm3 for the LGR cell, goFlux's example
# value for the larger benchtop UGGA; the microportable used here is
# ANALYZER_VOLUME_CM3 (code/lib/chamber_constants.R).
#
# Tcham: the field-recorded chamber air temperature; where the sheet has none, the
#   tower air temperature of the hour containing the deployment start, in true time
#   (code/lib/weather_time.R, the rule of every stage-A aux file). Pcham: 101.325 kPa.
# met_air_temp / met_time_diff_hours: the nearest weather-logger record within 6 h of
#   the start (logger clock), kept as a record in the dimensions table; not used for Tcham.
#
# Inputs:  data/raw/field_data/static_chamber_field/Gewirtzman_datasheets_fixing_dates.csv
#          data/raw/field_data/static_chamber_dims/{surface_area,additional_vol,simplified_volume}.csv
#          data/raw/weather/ymf_clean_sorted.csv
# Outputs: data/processed/flux/flux_data_with_dimensions_weather_and_final_airtemp.csv
#          data/processed/flux/goflux_auxfile.csv
# ==============================================================================

library(dplyr)
library(readr)
library(purrr)

source("../../lib/chamber_constants.R")  # ANALYZER_VOLUME_CM3
source("../../lib/weather_time.R")       # read_station_air(), station_air_at()

DIMS    <- "../../../data/raw/field_data/static_chamber_dims"
WEATHER <- "../../../data/raw/weather/ymf_clean_sorted.csv"

flux_data      <- read_csv("../../../data/raw/field_data/static_chamber_field/Gewirtzman_datasheets_fixing_dates.csv",
                           show_col_types = FALSE)
surface_area   <- read_csv(file.path(DIMS, "surface_area.csv"), show_col_types = FALSE)
additional_vol <- read_csv(file.path(DIMS, "additional_vol.csv"), show_col_types = FALSE)
simplified_vol <- read_csv(file.path(DIMS, "simplified_volume.csv"), show_col_types = FALSE)
weather_data   <- read_csv(WEATHER, show_col_types = FALSE)
cat("Field sheet:", nrow(flux_data), "rows; weather records:", nrow(weather_data), "\n")
cat("Chamber IDs on the sheet:",
    paste(sort(unique(flux_data$`Chamber ID`[!is.na(flux_data$`Chamber ID`)])), collapse = ", "), "\n")

# ---- chamber dimensions -------------------------------------------------------
surface_area_clean <- surface_area %>%
  mutate(
    chamber_series = case_when(
      grepl("^A", `Chamber ID`, ignore.case = TRUE) ~ "A",
      grepl("^B", `Chamber ID`, ignore.case = TRUE) ~ "B",
      grepl("^C", `Chamber ID`, ignore.case = TRUE) ~ "C",
      grepl("^D", `Chamber ID`, ignore.case = TRUE) ~ "D",
      TRUE ~ `Chamber ID`
    ),
    surface_area_cm2 = as.numeric(`SA cm2`)
  ) %>%
  select(chamber_series, surface_area_cm2) %>%
  filter(!is.na(surface_area_cm2))

analyzer_volumes <- additional_vol %>%
  mutate(
    analyzer_name = case_when(
      tolower(instrument) == "lgr_mgga" ~ "LGR",
      tolower(instrument) == "lgr" ~ "LGR",
      tolower(instrument) == "picarro" ~ "Picarro",
      TRUE ~ toupper(instrument)
    ),
    analyzer_cell_volume_cm3 = ifelse(analyzer_name == "LGR", ANALYZER_VOLUME_CM3,
                                      as.numeric(analyzer_cell)),
    tubing_volume_cm3 = as.numeric(tubing),
    filter_volume_cm3 = 0
  ) %>%
  select(analyzer_name, analyzer_cell_volume_cm3, tubing_volume_cm3, filter_volume_cm3)

names(simplified_vol) <- gsub("\\n|\\r", "_", names(simplified_vol))
names(simplified_vol) <- gsub("\\s+", "_", names(simplified_vol))
chamber_volumes <- simplified_vol %>%
  filter(!is.na(`Chamber_Alt_ID`) & !is.na(`Total_Volume_(mL)`)) %>%
  mutate(
    chamber_series = `Chamber_Alt_ID`,
    chamber_volume_cm3 = as.numeric(`Total_Volume_(mL)`)  # mL = cm3
  ) %>%
  filter(!is.na(chamber_volume_cm3) & chamber_series %in% c("A", "B", "C", "D")) %>%
  group_by(chamber_series) %>%
  summarise(avg_chamber_volume_cm3 = mean(chamber_volume_cm3, na.rm = TRUE), .groups = "drop")

chamber_dimensions <- tibble(
  chamber_id = c(1, 2, 3, 4, 5, 6),
  chamber_series = c("B", "B", "C", "C", "D", "D"),
  chamber_pair = c("pair_1", "pair_1", "pair_2", "pair_2", "pair_3", "pair_3")
) %>%
  left_join(surface_area_clean, by = "chamber_series") %>%
  left_join(chamber_volumes, by = "chamber_series") %>%
  rename(chamber_volume_cm3 = avg_chamber_volume_cm3)
cat("\nChamber dimensions:\n"); print(chamber_dimensions)

# ---- deployments ----------------------------------------------------------------
result <- flux_data %>%
  mutate(
    # sheet date YYMMDD (e.g. 210806) -> Date
    date_clean = case_when(
      !is.na(Date) & Date >= 210000 & Date <= 220000 ~ {
        year_part <- floor(Date / 10000) + 2000
        month_part <- floor((Date %% 10000) / 100)
        day_part <- Date %% 100
        as.Date(paste(year_part, sprintf("%02d", month_part), sprintf("%02d", day_part), sep = "-"))
      },
      TRUE ~ as.Date(NA)
    ),
    # start time HHMMSS -> POSIXct (UTC label, field clock)
    datetime_flux = case_when(
      !is.na(date_clean) & !is.na(`Flux Start Time`) ~ {
        time_str <- sprintf("%06d", `Flux Start Time`)
        time_formatted <- paste0(substr(time_str, 1, 2), ":",
                                 substr(time_str, 3, 4), ":",
                                 substr(time_str, 5, 6))
        as.POSIXct(paste(as.character(date_clean), time_formatted),
                   format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
      },
      TRUE ~ as.POSIXct(NA)
    ),
    flux_id = paste(format(date_clean, "%Y%m%d"), Plot, Tree_ID, `Measurement height`,
                    `Chamber ID`, `Flux Start Time`, sep = "_"),
    analyzer_clean = "LGR",   # every 2021 deployment used the LGR3 microportable
    start_time = `Flux Start Time`,
    end_time = HMS_End
  ) %>%
  left_join(chamber_dimensions, by = c("Chamber ID" = "chamber_id")) %>%
  left_join(analyzer_volumes, by = c("analyzer_clean" = "analyzer_name"))

# nearest weather-logger record within 6 h (record only; see header)
weather_clean <- weather_data %>%
  mutate(
    datetime_weather = as.POSIXct(TIMESTAMP, tz = "UTC"),
    air_temp_weather = Tair_Avg,
    pressure_weather = 101.325
  ) %>%
  filter(!is.na(datetime_weather) & !is.na(air_temp_weather)) %>%
  select(datetime_weather, air_temp_weather, pressure_weather) %>%
  arrange(datetime_weather)

find_closest_weather <- function(flux_datetime, weather_data) {
  none <- data.frame(air_temp_weather = NA, pressure_weather = NA, time_diff_hours = NA)
  if (is.na(flux_datetime)) return(none)
  time_diffs <- abs(as.numeric(weather_data$datetime_weather - flux_datetime))
  closest_idx <- which.min(time_diffs)
  closest_time_diff_hours <- min(time_diffs) / 3600
  if (closest_time_diff_hours > 6) return(none)
  data.frame(air_temp_weather = weather_data$air_temp_weather[closest_idx],
             pressure_weather = weather_data$pressure_weather[closest_idx],
             time_diff_hours = closest_time_diff_hours)
}
weather_matches <- map_dfr(result$datetime_flux, ~find_closest_weather(.x, weather_clean))

result <- result %>%
  bind_cols(weather_matches) %>%
  rename(met_air_temp = air_temp_weather,
         met_pressure = pressure_weather,
         met_time_diff_hours = time_diff_hours) %>%
  mutate(
    total_system_volume_cm3 = chamber_volume_cm3 + tubing_volume_cm3 + analyzer_cell_volume_cm3 + filter_volume_cm3,
    total_system_volume_L = total_system_volume_cm3 / 1000
  ) %>%
  select(
    date = date_clean,
    datetime = datetime_flux,
    personnel = Personnel,
    plot = Plot,
    tree_id = Tree_ID,
    species = Tree_sp,
    measurement_height = `Measurement height`,
    chamber_id = `Chamber ID`,
    temp_stem = Temp_stem,
    temp_air = Temp_air,
    dbh = DBH,
    start_time,
    end_time,
    analyzer = analyzer_clean,
    gas_analyser_original = Gas_analyser,
    flux_id,
    chamber_series,
    chamber_pair,
    surface_area_cm2,
    chamber_volume_cm3,
    analyzer_cell_volume_cm3,
    tubing_volume_cm3,
    filter_volume_cm3,
    total_system_volume_cm3,
    total_system_volume_L,
    met_air_temp,
    met_pressure,
    met_time_diff_hours
  ) %>%
  arrange(datetime)

# Tcham: field air temperature, else the tower hour (see header)
station_air <- read_station_air(WEATHER)
result <- result %>%
  mutate(
    air_temp_final = ifelse(!is.na(temp_air), temp_air, station_air_at(datetime, station_air)),
    temp_source = case_when(!is.na(temp_air) ~ "field_air",
                            !is.na(air_temp_final) ~ "tower_hour",
                            TRUE ~ NA_character_)
  )

cat("\nChamber air temperature by source:\n"); print(table(result$temp_source, useNA = "ifany"))
cat("Deployments:", nrow(result), "| with area and volume:",
    sum(!is.na(result$surface_area_cm2) & !is.na(result$total_system_volume_L)), "\n")

write_csv(result, "../../../data/processed/flux/flux_data_with_dimensions_weather_and_final_airtemp.csv")
cat("Wrote data/processed/flux/flux_data_with_dimensions_weather_and_final_airtemp.csv\n")

# ---- goFlux auxfile ---------------------------------------------------------------
# No obs.length column: the closure windows are the hand-picked ones
# (data/raw/flux_windows/), and goFlux obs.win() is given its search length as an
# argument (03_goflux_trees_2021.R).
auxfile_goflux <- result %>%
  filter(!is.na(air_temp_final) & !is.na(total_system_volume_L) & !is.na(surface_area_cm2)) %>%
  mutate(
    UniqueID = flux_id,
    start.time = datetime,
    Area = surface_area_cm2,       # cm2
    Vtot = total_system_volume_L,  # L
    Tcham = air_temp_final,        # degC
    Pcham = 101.325                # kPa
  ) %>%
  select(UniqueID, start.time, Area, Vtot, Tcham, Pcham,
         plot, tree_id, species, chamber_id, measurement_height, temp_source) %>%
  arrange(start.time)

cat("\ngoFlux auxfile:", nrow(auxfile_goflux), "deployments\n")
cat("  Tcham range:", round(range(auxfile_goflux$Tcham, na.rm = TRUE), 1), "degC\n")
cat("  Area range:", round(range(auxfile_goflux$Area, na.rm = TRUE), 1), "cm2\n")
cat("  Vtot range:", round(range(auxfile_goflux$Vtot, na.rm = TRUE), 3), "L\n")
write_csv(auxfile_goflux, "../../../data/processed/flux/goflux_auxfile.csv")
cat("Wrote data/processed/flux/goflux_auxfile.csv\n")
