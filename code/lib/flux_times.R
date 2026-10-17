# ==============================================================================
# flux_times.R -- clock and lookup helpers for the flux tables (code/tables/1x_flux_*.R)
# ------------------------------------------------------------------------------
# goFlux start.time and the saved window times end in "Z" but hold the analyzer's
# wall clock. They are read as wall-clock digits (POSIXct in UTC, no shift), corrected
# where the analyzer clock was off, and written as local time with the
# America/New_York offset of that date: 2020-06-15T13:57:00-04:00.
# ==============================================================================
FLUX_TZ <- "America/New_York"

# "2020-06-15T13:57:00.000000Z" -> wall-clock POSIXct (UTC container, digits unchanged)
wall_clock <- function(x) as.POSIXct(substr(sub("T", " ", as.character(x)), 1, 19), tz = "UTC")

# wall-clock POSIXct -> "YYYY-MM-DDTHH:MM:SS-04:00" (offset of that local date and time)
iso_local <- function(t) {
  lt <- as.POSIXct(format(t, "%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = FLUX_TZ)
  z <- sub("(\\d{2})(\\d{2})$", "\\1:\\2", format(lt, "%z"))
  ifelse(is.na(t), NA_character_, paste0(format(lt, "%Y-%m-%dT%H:%M:%S"), z))
}

# hand-picked fitting windows (one per UniqueID unless keep_all)
read_windows <- function(f, keep_all = FALSE) {
  w <- unique(as.data.frame(data.table::fread(f, select = c("UniqueID", "start.time_corr", "end.time_corr"),
                                              colClasses = "character")))
  if (!keep_all) stopifnot(!anyDuplicated(w$UniqueID))
  data.frame(UniqueID = w$UniqueID, win_start = wall_clock(w$start.time_corr), win_end = wall_clock(w$end.time_corr))
}

# "2020-06-15T13:57:00-04:00" -> true instant (POSIXct UTC)
iso_instant <- function(x) as.POSIXct(substr(x, 1, 19), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC") -
  (as.numeric(substr(x, 20, 22)) * 3600 + sign(as.numeric(substr(x, 20, 22))) * as.numeric(substr(x, 24, 25)) * 60)

# met-tower hourly mean air temperature (weather.csv air_temp_avg_C) for a local wall-clock
# time t (as wall_clock() returns): the hour whose centre is nearest t in true time
# (weather.csv timestamp; the hour containing t), within 24 h.
tower_air <- local({
  WX <- NULL
  function(t) {
    if (is.null(WX)) {
      w <- read.csv("data/tables/weather.csv", colClasses = c(timestamp = "character"))
      WX <<- data.frame(mid = as.numeric(iso_instant(w$timestamp)) - 1800, tair = w$air_temp_avg_C)
      WX <<- WX[!is.na(WX$tair), ]; WX <<- WX[order(WX$mid), ]
    }
    tu <- as.numeric(as.POSIXct(format(t, "%Y-%m-%d %H:%M:%S", tz = "UTC"), tz = FLUX_TZ))   # true instant
    vapply(tu, function(x) {
      if (is.na(x)) return(NA_real_)
      d <- abs(WX$mid - x); j <- which.min(d)
      if (d[j] <= 86400) WX$tair[j] else NA_real_
    }, numeric(1))
  }
})

# read.csv, with the columns in `readr_cols` parsed as readr::read_csv parses them:
# the values the model has always used (it read these files with readr). The two
# parsers disagree for some long decimals (readr by up to ~1e-13 relative, and badly
# below 1e-307), so only the columns the model reads are taken from readr.
read_flux_csv <- function(f, readr_cols, ...) {
  d <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, ...)
  txt <- utils::read.csv(f, check.names = FALSE, colClasses = "character", na.strings = character(0), ...)
  stopifnot(all(readr_cols %in% names(d)))
  for (j in readr_cols) {
    v <- suppressWarnings(readr::parse_double(trimws(txt[[j]]), na = c("", "NA")))
    stopifnot(is.double(d[[j]]) || is.integer(d[[j]]), identical(is.na(v), is.na(d[[j]])),
              isTRUE(all.equal(v, as.double(d[[j]]), tolerance = 1e-12)))
    d[[j]] <- v
  }
  d
}
