#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 34_weather.R -- weather.csv (data/tables/)
# ------------------------------------------------------------------------------
# The Yale-Myers weather station, hourly, 2018-09-29 to 2023-05-12 (40,458
# records), from data/raw/weather/ymf_clean_sorted.csv; the columns of the
# archived environmental_timeseries.csv plus the logger's record number. Values
# as logged (not screened: e.g. wind speed up to 289 m/s, vwc down to -8). Two
# columns renamed for their units: solar_rad_avg_kW -> solar_rad_avg_W_m2 (the
# logger's "kW" column peaks near 1000), solar_rad_total_MJ -> _MJ_m2.
#
# TIME BASIS (checked here; outputs/audit/weather_build.txt):
#   - timestamp_logger is the logger clock as recorded. Each record closes the
#     hour it summarises (Tair_TMx, the time of the hour's maximum, falls in
#     (timestamp - 1 h, timestamp] for all but 5 records).
#   - No daylight-saving shifts: at every DST change the hour sequence is
#     unbroken (02:00 present each March, no repeated 01:00 each November).
#   - The clock was reset once, on 2021-08-03: the record stops for 4 h (07:00 ->
#     12:00 with consecutive RECORD numbers 27251 -> 27252, i.e. the logger was
#     stopped -- the 12 one-hour gaps elsewhere are the same kind of stop), and
#     the diurnal cycle moves by ~2.9 h. Solar noon on the logger clock
#     (midpoint of the half-maximum radiation times, interval centres) is 09:32
#     before and 12:26 after; expected 11:48 EST / 12:48 EDT at -72.13 deg. With the
#     same ~0.35 h early bias of that estimate in both periods, the clock kept
#     about UTC-7 before the reset and about UTC-4 (EDT, all year) after it.
#     before_clock_reset marks the records stamped before the reset. The
#     archived file labelled every stamp UTC ("Z"); it is neither UTC nor
#     reliably local time.
#
# Check: values equal data/compiled/environmental_timeseries.csv.
# Writes data/tables/weather.csv, outputs/audit/weather_build.txt
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
OUT <- "data/tables"
SRC <- "data/raw/weather/ymf_clean_sorted.csv"
RESET_RECORD <- 27252                                # first record stamped after the reset

w <- read_csv(SRC, show_col_types = FALSE, name_repair = "unique_quiet")
ts <- as.POSIXct(w$TIMESTAMP, tz = "UTC")             # clock arithmetic only; not a time zone
stopifnot(nrow(w) == 40458, !anyNA(ts), !anyDuplicated(ts))
gap <- c(NA, diff(as.numeric(ts)) / 3600); rgap <- c(NA, diff(w$RECORD))
i_reset <- which(w$RECORD == RESET_RECORD)
stopifnot(length(i_reset) == 1, gap[i_reset] == 5, rgap[i_reset] == 1)

W <- w %>% transmute(
  timestamp_logger = format(ts, "%Y-%m-%dT%H:%M:%S"), record = RECORD,
  before_clock_reset = seq_len(n()) < i_reset,
  year, month, day, hour, air_temp_C = Tair, air_temp_avg_C = Tair_Avg,
  air_temp_max_C = Tair_Max, air_temp_min_C = Tair_Min, soil_temp_avg_C = Tsoil_Avg,
  relative_humidity_pct = RH, vwc = VWC, vwc_avg = VWC_Avg, precip_mm = Rain_mm_Tot,
  wind_speed_avg_ms = WindSpeed_ms_Avg, wind_speed_max_ms = WindSpeed_ms_Max, wind_dir_deg = WindDir,
  solar_rad_avg_W_m2 = Solar_Rad_kW_Avg, solar_rad_total_MJ_m2 = Solar_Rad_Tot_MJ_Tot,
  dewpoint_C = TdewPointC, ET_ref = ETos)

# ---- time-basis evidence ------------------------------------------------------------------
tmx <- as.POSIXct(w$Tair_TMx, format = "%m/%d/%y %H:%M", tz = "UTC")
tmx_lag <- as.numeric(ts - tmx, units = "hours")
dst <- c("2019-03-10", "2019-11-03", "2020-03-08", "2020-11-01", "2021-03-14", "2021-11-07", "2022-03-13", "2022-11-06")
dst_ok <- sapply(dst, function(d) identical(format(ts[format(ts, "%Y-%m-%d") == d], "%H"), sprintf("%02d", 0:23)))
w$hc <- as.numeric(format(ts, "%H")) - 0.5          # interval centre
noon <- function(x) { a <- tapply(x$Solar_Rad_kW_Avg, x$hc, mean); h <- as.numeric(names(a)); hm <- max(a) / 2
  i <- which(a >= hm); r <- approx(a[(min(i) - 1):min(i)], h[(min(i) - 1):min(i)], hm)$y
  f <- approx(a[max(i):(max(i) + 1)], h[max(i):(max(i) + 1)], hm)$y; (r + f) / 2 }
sn <- sapply(split(w, W$before_clock_reset), noon)

# ---- check vs the archived table ---------------------------------------------------------------
CP <- read_csv("data/compiled/environmental_timeseries.csv", show_col_types = FALSE)
stopifnot(nrow(CP) == nrow(W), all(format(CP$timestamp, "%Y-%m-%dT%H:%M:%S", tz = "UTC") == W$timestamp_logger))
names(CP)[names(CP) == "solar_rad_avg_kW"] <- "solar_rad_avg_W_m2"           # renamed: the values are W m-2
names(CP)[names(CP) == "solar_rad_total_MJ"] <- "solar_rad_total_MJ_m2"
vals <- setdiff(names(CP), "timestamp")
bad <- sapply(vals, function(v) sum(!((is.na(CP[[v]]) & is.na(W[[v]])) |
  (!is.na(CP[[v]]) & !is.na(W[[v]]) & abs(CP[[v]] - W[[v]]) <= 1e-12 * pmax(1, abs(W[[v]]))))))   # floating point

sink(out_path("weather_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("weather: %d hourly records, %s to %s (logger clock)\n", nrow(W), W$timestamp_logger[1], W$timestamp_logger[nrow(W)]))
cat(sprintf("gaps between consecutive stamps: %s\n", paste(names(table(gap)), table(gap), sep = " h x", collapse = "; ")))
cat(sprintf("Tair_TMx within (stamp - 1 h, stamp]: %d of %d records (stamps close the interval; the 5 exceptions sit at the\n  2022-01-31 repeated RECORD and the 2022-06-11 gap)\n",
            sum(tmx_lag >= 0 & tmx_lag < 1, na.rm = TRUE), sum(!is.na(tmx_lag))))
cat(sprintf("DST days with an unbroken 00-23 hour sequence: %d of %d\n", sum(dst_ok), length(dst)))
cat(sprintf("clock reset: %s -> %s (logger stopped), RECORD %d -> %d\n", W$timestamp_logger[i_reset - 1], W$timestamp_logger[i_reset],
            W$record[i_reset - 1], W$record[i_reset]))
cat(sprintf("solar noon on the logger clock (half-max midpoint): before reset %.2f h, after %.2f h (expected 11.80 EST, 12.80 EDT)\n",
            sn[["TRUE"]], sn[["FALSE"]]))
cat("\nvs data/compiled/environmental_timeseries.csv -- values differing:\n"); print(bad)
options(op); sink()
cat(readLines(out_path("weather_build.txt", "audit")), sep = "\n")
if (any(bad > 0)) stop("weather does not reproduce environmental_timeseries.csv")

write_csv(W, file.path(OUT, "weather.csv"), na = "")
