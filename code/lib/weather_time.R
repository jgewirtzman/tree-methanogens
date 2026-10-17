# weather_time.R -- met-station logger clock to true time, and station air temperature
# at a deployment. One rule, used by code/tables/06_weather.R (weather.csv timestamp)
# and by the stage-A goFlux aux files (chamber temperature), which run before
# weather.csv exists. Evidence for the rule: header of 06_weather.R.
#   The logger kept UTC-7 until its clock reset on 2021-08-03 (first record after the
#   reset: RECORD 27252) and UTC-4 after it. Each record closes the hour it summarises.
WEATHER_RESET_RECORD <- 27252
WEATHER_LOGGER_UTC_H <- c(before = -7, after = -4)    # logger clock = UTC + this

# logger stamps (POSIXct holding the logger's digits in UTC) and RECORD numbers, in
# record order -> true instant (POSIXct UTC)
weather_logger_to_utc <- function(ts, record) {
  i_reset <- which(record == WEATHER_RESET_RECORD)
  stopifnot(length(i_reset) == 1)
  before <- seq_along(ts) < i_reset
  ts - ifelse(before, WEATHER_LOGGER_UTC_H[["before"]], WEATHER_LOGGER_UTC_H[["after"]]) * 3600
}

# hourly mean air temperature from the raw station file, keyed by the true centre of
# each hour (seconds since the epoch)
read_station_air <- function(path = "data/raw/weather/ymf_clean_sorted.csv") {
  w <- utils::read.csv(path)
  ts <- as.POSIXct(w$TIMESTAMP, tz = "UTC")
  stopifnot(!anyNA(ts), !is.unsorted(ts, strictly = TRUE))
  utc <- weather_logger_to_utc(ts, w$RECORD)
  out <- data.frame(mid = as.numeric(utc) - 1800, tair = as.numeric(w$Tair_Avg))
  out[!is.na(out$tair), ]
}

# station air temperature at local wall-clock times t (POSIXct holding local digits in
# UTC, as the aux files store start.time): the hour whose true centre is nearest, i.e.
# the hour containing t, within max_h hours; NA otherwise
station_air_at <- function(t, wx, max_h = 24) {
  tu <- as.numeric(as.POSIXct(format(t, "%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = "America/New_York"))
  vapply(tu, function(x) {
    if (is.na(x)) return(NA_real_)
    d <- abs(wx$mid - x); j <- which.min(d)
    if (d[j] <= max_h * 3600) wx$tair[j] else NA_real_
  }, numeric(1))
}
