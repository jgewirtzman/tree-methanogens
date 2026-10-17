#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 06_weather.R -- weather.csv (data/tables/)
# ------------------------------------------------------------------------------
# The Yale-Myers weather station, hourly, 2018-09-29 to 2023-05-12 (40,458
# records), from data/raw/weather/ymf_clean_sorted.csv; the columns of the
# archived environmental_timeseries.csv plus the logger's record number and the
# true time. Two columns renamed for their units: solar_rad_avg_kW ->
# solar_rad_avg_W_m2 (the logger's "kW" column peaks near 1000),
# solar_rad_total_MJ -> _MJ_m2.
#
# TIME BASIS (checked here; outputs/audit/weather_build.txt):
#   - timestamp_logger is the logger clock as recorded. Each record closes the
#     hour it summarises (Tair_TMx, the time of the hour's maximum, falls in
#     (timestamp - 1 h, timestamp] for all but 5 records).
#   - No daylight-saving shifts: at every DST change the hour sequence is
#     unbroken (02:00 present each March, no repeated 01:00 each November).
#   - The clock was reset once, on 2021-08-03: the logger stopped between RECORD
#     27251 (stamped 07:00) and 27252 (stamped 12:00).
#   - Rule (code/lib/weather_time.R, also used by the stage-A goFlux aux files for
#     chamber temperature): the logger kept UTC-7 before the reset and UTC-4
#     (daylight time, all year) after it. timestamp = timestamp_logger + 7 h before the reset, + 4 h
#     after, written as local time with the America/New_York offset of that
#     instant (as flux_stem.start_time). Evidence: the centroid of each day's
#     shortwave radiation (the clearer half of the days of each month) against
#     solar noon at the site puts the logger at UTC-7.24 h before the reset and
#     UTC-4.34 h after (medians; drift under 0.1 h a year in either period); the
#     early bias of the centroid (~0.3 h) is common to both periods. On
#     2021-08-03 itself the morning (before the stop) and the evening (after it)
#     both fall at the true hours of sunrise and sunset, so the change is at the
#     stop. The archived file labelled every stamp UTC ("Z"); it is neither.
#
# SCREENING (values set NA here, raw unchanged; every removed value is listed in
# outputs/audit/weather_screened.csv):
#   bounds  air temperatures -40 to 45 C; dew point -50 to 35 C; relative
#           humidity 0-100 %; precipitation 0-100 mm per hour; solar 0-1361 W m-2
#           (an hourly mean cannot exceed the solar constant) and 0-4.9 MJ m-2;
#           station VWC above 0 and at most 0.6 (zero and negative readings are the
#           probe failure of Jun-Sep 2020); ET_ref (logger ETos) -0.1 to 1.5 mm
#           per hour.
#   stuck   wind direction unchanged for 6 h or more. No temperature, dew point,
#           humidity or solar series holds a value long enough to be stuck.
#   WIND SPEED (wind_speed_avg_ms, wind_speed_max_ms) is blanked: implausible
#           throughout (hourly means above 40 m s-1 in over 17,000 hours; a median near
#           7 m s-1 under forest without those; unit not recorded) and unused
#           (04_drivers/01_wb_reference_et.R takes the FAO-56 default of 2 m s-1).
#   Not screened: soil_temp_avg_C (it tracks air temperature, so the probe is not
#   in the soil; 04_drivers/00_soil_temp_climatology.R). The logger computed ET_ref
#   from the wind speed; it is not used.
# USED BY: air_temp_avg_C -- flux lookups (code/lib/flux_times.R tower_air),
#   monthly air drivers (05_model/01_load_and_prep_data.R), daily air
#   (04_drivers/01_wb_reference_et.R, 00_soil_temp_climatology.R); air max/min,
#   relative humidity, solar total, precipitation and VWC -- reference ET and
#   water balance (04_drivers/01_wb_reference_et.R, which sets wind to the
#   FAO-56 default).
#
# (The check against the old archive, data/compiled/environmental_timeseries.csv, was
# removed 2026-10-16: the archive is now a copy of data/tables/.)
# Writes data/tables/weather.csv, outputs/audit/weather_build.txt,
#   outputs/audit/weather_screened.csv
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/weather_time.R")
OUT <- "data/tables"
SRC <- "data/raw/weather/ymf_clean_sorted.csv"
RESET_RECORD <- WEATHER_RESET_RECORD                 # first record stamped after the reset
LOGGER_UTC_H <- WEATHER_LOGGER_UTC_H                 # logger clock = UTC + this
TZ <- "America/New_York"

w <- read_csv(SRC, show_col_types = FALSE, name_repair = "unique_quiet")
ts <- as.POSIXct(w$TIMESTAMP, tz = "UTC")             # clock arithmetic only; not a time zone
stopifnot(nrow(w) == 40458, !anyNA(ts), !anyDuplicated(ts))
gap <- c(NA, diff(as.numeric(ts)) / 3600); rgap <- c(NA, diff(w$RECORD))
i_reset <- which(w$RECORD == RESET_RECORD)
stopifnot(length(i_reset) == 1, gap[i_reset] == 5, rgap[i_reset] == 1)
before <- seq_len(nrow(w)) < i_reset
utc <- weather_logger_to_utc(ts, w$RECORD)                                             # true instant
stopifnot(!is.unsorted(utc, strictly = TRUE))
iso_tz <- function(t) paste0(format(t, "%Y-%m-%dT%H:%M:%S", tz = TZ),
                             sub("(\\d{2})(\\d{2})$", "\\1:\\2", format(t, "%z", tz = TZ)))

W <- w %>% transmute(
  timestamp = iso_tz(utc), timestamp_logger = format(ts, "%Y-%m-%dT%H:%M:%S"), record = RECORD,
  before_clock_reset = before,
  year, month, day, hour, air_temp_C = Tair, air_temp_avg_C = Tair_Avg,
  air_temp_max_C = Tair_Max, air_temp_min_C = Tair_Min, soil_temp_avg_C = Tsoil_Avg,
  relative_humidity_pct = RH, vwc = VWC, vwc_avg = VWC_Avg, precip_mm = Rain_mm_Tot,
  wind_speed_avg_ms = WindSpeed_ms_Avg, wind_speed_max_ms = WindSpeed_ms_Max, wind_dir_deg = WindDir,
  solar_rad_avg_W_m2 = Solar_Rad_kW_Avg, solar_rad_total_MJ_m2 = Solar_Rad_Tot_MJ_Tot,
  dewpoint_C = TdewPointC, ET_ref = ETos)
stopifnot(!anyDuplicated(W$timestamp))

# ---- time-basis evidence ------------------------------------------------------------------
tmx <- as.POSIXct(w$Tair_TMx, format = "%m/%d/%y %H:%M", tz = "UTC")
tmx_lag <- as.numeric(ts - tmx, units = "hours")
dst <- c("2019-03-10", "2019-11-03", "2020-03-08", "2020-11-01", "2021-03-14", "2021-11-07", "2022-03-13", "2022-11-06")
dst_ok <- sapply(dst, function(d) identical(format(ts[format(ts, "%Y-%m-%d") == d], "%H"), sprintf("%02d", 0:23)))
lon <- as.numeric(read_csv(file.path(OUT, "site.csv"), show_col_types = FALSE)$ref_longitude)
solar_noon_utc <- function(d) {                       # hours UTC (NOAA equation of time)
  g <- 2 * pi / 365 * (as.numeric(format(d, "%j")) - 1)
  eot <- 229.18 * (0.000075 + 0.001868 * cos(g) - 0.032077 * sin(g) - 0.014615 * cos(2 * g) - 0.040849 * sin(2 * g))
  (720 - 4 * lon - eot) / 60
}
CEN <- tibble(day = as.Date(format(ts, "%Y-%m-%d")), tc = as.numeric(format(ts, "%H")) - 0.5,   # interval centre
              rad = pmax(w$Solar_Rad_kW_Avg, 0), before) %>%
  group_by(day) %>% filter(n() == 24, n_distinct(before) == 1) %>%
  summarise(cen = sum(rad * tc) / sum(rad), tot = sum(rad), before = before[1], .groups = "drop") %>%
  mutate(off = cen - solar_noon_utc(day)) %>%
  group_by(m = format(day, "%Y-%m")) %>% filter(tot >= median(tot)) %>% ungroup()
off_med <- tapply(CEN$off, CEN$before, median)
drift <- sapply(split(CEN, CEN$before), function(x) coef(lm(off ~ as.numeric(day), x))[[2]] * 365.25)
stopifnot(abs(off_med[["TRUE"]] - LOGGER_UTC_H[["before"]]) < 0.5, abs(off_med[["FALSE"]] - LOGGER_UTC_H[["after"]]) < 0.5)
rd <- W %>% mutate(utc = utc) %>% filter(substr(timestamp, 1, 10) == "2021-08-03", solar_rad_avg_W_m2 > 20)
sun_reset <- c(first = rd$timestamp[1], last = rd$timestamp[nrow(rd)],
               first_logger = rd$timestamp_logger[1], last_logger = rd$timestamp_logger[nrow(rd)])

# ---- screening ---------------------------------------------------------------------------------
BOUNDS <- tribble(~var, ~lo, ~hi,
  "air_temp_C", -40, 45, "air_temp_avg_C", -40, 45, "air_temp_max_C", -40, 45, "air_temp_min_C", -40, 45,
  "dewpoint_C", -50, 35, "relative_humidity_pct", 0, 100, "precip_mm", 0, 100,
  "solar_rad_avg_W_m2", 0, 1361, "solar_rad_total_MJ_m2", 0, 4.9,
  "vwc", 1e-9, 0.6, "vwc_avg", 1e-9, 0.6, "ET_ref", -0.1, 1.5)
run_len <- function(x) { r <- rle(x); rep(r$lengths, r$lengths) }    # length of the run each value sits in
flag <- list()
for (k in seq_len(nrow(BOUNDS))) { v <- BOUNDS$var[k]; x <- W[[v]]
  i <- which(!is.na(x) & (x < BOUNDS$lo[k] | x > BOUNDS$hi[k]))
  flag[[length(flag) + 1]] <- tibble(i, variable = v, rule = sprintf("outside %g to %g", BOUNDS$lo[k], BOUNDS$hi[k])) }
vane <- which(!is.na(W$wind_dir_deg) & run_len(W$wind_dir_deg) >= 6)
flag[[length(flag) + 1]] <- tibble(i = vane, variable = "wind_dir_deg", rule = "wind direction unchanged for 6 h or more")
for (v in c("air_temp_C", "air_temp_avg_C", "air_temp_max_C", "air_temp_min_C", "dewpoint_C"))   # nothing stuck here
  stopifnot(max(run_len(W[[v]])) < 6)
FL <- bind_rows(flag) %>% mutate(value = mapply(function(r, v) W[[v]][r], i, variable)) %>%
  filter(!is.na(value)) %>% distinct(i, variable, .keep_all = TRUE) %>% arrange(variable, i) %>%
  mutate(timestamp = W$timestamp[i], record = W$record[i]) %>% select(timestamp, record, variable, value, rule)
for (v in unique(FL$variable)) W[[v]][match(FL$record[FL$variable == v], W$record)] <- NA
n_wind <- c(avg = sum(!is.na(W$wind_speed_avg_ms)), max = sum(!is.na(W$wind_speed_max_ms)))
W$wind_speed_avg_ms <- NA_real_; W$wind_speed_max_ms <- NA_real_          # blanked (see header)
write_csv(FL, out_path("weather_screened.csv", "audit"))

sink(out_path("weather_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("weather: %d hourly records, %s to %s (true time; logger clock %s to %s)\n", nrow(W), W$timestamp[1],
            W$timestamp[nrow(W)], W$timestamp_logger[1], W$timestamp_logger[nrow(W)]))
cat(sprintf("gaps between consecutive logger stamps: %s\n", paste(names(table(gap)), table(gap), sep = " h x", collapse = "; ")))
cat(sprintf("Tair_TMx within (stamp - 1 h, stamp]: %d of %d records (stamps close the interval; the 5 exceptions sit at the\n  2022-01-31 repeated RECORD and the 2022-06-11 gap)\n",
            sum(tmx_lag >= 0 & tmx_lag < 1, na.rm = TRUE), sum(!is.na(tmx_lag))))
cat(sprintf("DST days with an unbroken 00-23 hour sequence: %d of %d\n", sum(dst_ok), length(dst)))
cat(sprintf("clock reset: logger %s -> %s (logger stopped), RECORD %d -> %d; true time %s -> %s\n",
            W$timestamp_logger[i_reset - 1], W$timestamp_logger[i_reset], W$record[i_reset - 1], W$record[i_reset],
            W$timestamp[i_reset - 1], W$timestamp[i_reset]))
cat(sprintf("logger offset from UTC, daily shortwave centroid vs solar noon at %.3f deg (clearer half of each month's days):\n", lon))
cat(sprintf("  before the reset %.2f h (%d days), after %.2f h (%d days); drift %+.3f / %+.3f h per year\n",
            off_med[["TRUE"]], sum(CEN$before), off_med[["FALSE"]], sum(!CEN$before), drift[["TRUE"]], drift[["FALSE"]]))
cat(sprintf("  by year: %s\n", paste(CEN %>% group_by(y = paste0(substr(m, 1, 4), ifelse(before, "", "*"))) %>%
            summarise(o = sprintf("%s %.2f", y[1], median(off)), .groups = "drop") %>% pull(o), collapse = "; ")))
cat(sprintf("  rule applied: UTC%+d before, UTC%+d after (* = after the reset)\n", LOGGER_UTC_H[["before"]], LOGGER_UTC_H[["after"]]))
cat(sprintf("2021-08-03, radiation > 20 W m-2: first hour ending %s (logger %s, before the stop), last ending %s (logger %s, after it)\n",
            sun_reset[["first"]], sun_reset[["first_logger"]], sun_reset[["last"]], sun_reset[["last_logger"]]))
cat(sprintf("\nscreened (set NA): %d values; by variable and rule:\n", nrow(FL)))
print(as.data.frame(FL %>% count(variable, rule)), row.names = FALSE)
cat(sprintf("\nwind speed blanked: %d hourly means and %d maxima set NA (not screened; unused)\n", n_wind[["avg"]], n_wind[["max"]]))
options(op); sink()
cat(readLines(out_path("weather_build.txt", "audit")), sep = "\n")

write_csv(W, file.path(OUT, "weather.csv"), na = "")
