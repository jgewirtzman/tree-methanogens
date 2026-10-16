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

# met-tower Tair_Avg nearest a time, within 24 h -- the lookup 05_model/01_load_and_prep_data.R
# makes (tower stamps and the given times compared as recorded)
tower_air <- local({
  WX <- NULL
  function(t) {
    if (is.null(WX)) {
      w <- read.csv("data/raw/weather/ymf_clean_sorted.csv")
      WX <<- data.frame(ts = as.POSIXct(w$TIMESTAMP, tz = "UTC"), tair = w$Tair_Avg)
      WX <<- WX[!is.na(WX$tair), ]; WX <<- WX[order(WX$ts), ]
    }
    vapply(seq_along(t), function(i) {
      if (is.na(t[i])) return(NA_real_)
      d <- abs(as.numeric(WX$ts - t[i])); j <- which.min(d)
      if (d[j] <= 86400) WX$tair[j] else NA_real_
    }, numeric(1))
  }
})
