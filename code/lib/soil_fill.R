# soil_fill.R -- soil temperature and moisture for any measurement, taken only from
# data/tables/soil_readings.csv. This is the only place the rule is written; every
# script that needs soil conditions (model prep, drivers, figures, archive) calls
# attach_soil(). The flux tables carry no soil values.
#
# attach_soil(df, id_col, time_col, unit = c("collar", "tree"))
#   df        any data frame
#   id_col    column holding the collar_id (unit "collar") or tree_id (unit "tree")
#   time_col  column holding the measurement time: an ISO 8601 string (the tables'
#             start_time, local time; its first 10 characters are the local date),
#             a Date, or a POSIXct (its date in its own time zone)
#   returns   df plus
#     soil_temp_C        degC
#     soil_vwc_pct       volumetric water content, percent
#     soil_match         how the values were found (first match wins):
#       same_date              the unit's own readings: the mean of all of them dated
#                              within 3 days of the measurement, one of them on its date
#                              (a tree: its holes -- the 2021 survey read some trees'
#                              four holes over two days; a collar: that collar-date's
#                              records; no collar was read twice within 6 days)
#       nearest_date           as same_date, none of them on the measurement date
#       position_same_date     the mean over the collars at the unit's landscape
#                              position (U, I, WD, WS) read that date (a collar's
#                              position is its prefix; a tree's is trees.csv
#                              landscape_position, where it is one of these four)
#       position_nearest_date  that position mean, nearest date within 3 days
#       none                   nothing within 3 days: both values NA
#     soil_reading_date  date of the readings used (own readings: the nearest of them)
# A reading counts if it has a temperature or a VWC; both values come from the matched
# readings (one may be NA). Equal distances in date: the earlier date. The position
# means are over collar readings only (the tree holes of the 2021 survey and the 2023
# field-sheet readings sit at trees, not at a landscape position).
SOIL_FILL_MAX_DAYS <- 3
SOIL_POSITIONS <- c("U", "I", "WD", "WS")

.soil_tables_dir <- function() {
  d <- Filter(dir.exists, c("data/tables", "../data/tables", "../../data/tables", "../../../data/tables"))[1]
  if (is.na(d)) stop("data/tables not found from ", getwd(), call. = FALSE)
  d
}

# one row per unit x date: mean of that unit's readings that date
.soil_unit_dates <- function(sr, id) {
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  sr <- sr[!is.na(sr[[id]]) & !is.na(sr$date) & (!is.na(sr$soil_temp_C) | !is.na(sr$soil_vwc_pct)), ]
  k <- paste(sr[[id]], sr$date)
  data.frame(unit = tapply(sr[[id]], k, `[`, 1), date = as.Date(tapply(as.character(sr$date), k, `[`, 1)),
             soil_temp_C = tapply(sr$soil_temp_C, k, mean_na), soil_vwc_pct = tapply(sr$soil_vwc_pct, k, mean_na),
             row.names = NULL, stringsAsFactors = FALSE)
}

.soil_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(format(x, "%Y-%m-%d")))
  x <- as.character(x); x[!is.na(x) & !grepl("^\\d{4}-\\d{2}-\\d{2}", x)] <- NA
  as.Date(substr(x, 1, 10))
}

attach_soil <- function(df, id_col, time_col, unit = c("collar", "tree")) {
  unit <- match.arg(unit)
  td <- .soil_tables_dir()
  sr <- utils::read.csv(file.path(td, "soil_readings.csv"), stringsAsFactors = FALSE, na.strings = "")
  sr$date <- as.Date(sr$date)
  idc <- if (unit == "collar") "collar_id" else "tree_id"
  own <- sr[!is.na(sr[[idc]]) & !is.na(sr$date) & (!is.na(sr$soil_temp_C) | !is.na(sr$soil_vwc_pct)),
            c(idc, "date", "soil_temp_C", "soil_vwc_pct")]
  names(own)[1] <- "unit"
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  col <- .soil_unit_dates(sr, "collar_id")
  col$position <- sub("_.*", "", col$unit)
  pos <- .soil_unit_dates(transform(col, pid = position), "pid")   # mean over the collars read that date
  ids <- as.character(df[[id_col]])
  d <- .soil_date(df[[time_col]])
  position <- if (unit == "collar") sub("_.*", "", ids) else {
    tr <- utils::read.csv(file.path(td, "trees.csv"), stringsAsFactors = FALSE, na.strings = "")
    tr$landscape_position[match(ids, tr$tree_id)]
  }
  position[!position %in% SOIL_POSITIONS] <- NA
  pick <- function(cand, day) {               # nearest date within the limit; earlier on ties
    if (!nrow(cand)) return(NULL)
    dd <- abs(as.numeric(cand$date - day))
    j <- order(dd, cand$date)[1]
    if (dd[j] <= SOIL_FILL_MAX_DAYS) cand[j, ] else NULL
  }
  n <- nrow(df)
  out <- data.frame(soil_temp_C = rep(NA_real_, n), soil_vwc_pct = NA_real_, soil_match = "none",
                    soil_reading_date = as.Date(NA), stringsAsFactors = FALSE)
  own_by <- split(own, own$unit); pos_by <- split(pos, pos$unit)
  for (i in seq_len(n)) {
    if (is.na(d[i])) next
    hit <- NULL
    if (!is.na(ids[i]) && !is.null(own_by[[ids[i]]])) {   # all own readings within the limit
      o <- own_by[[ids[i]]]; dd <- abs(as.numeric(o$date - d[i])); o <- o[dd <= SOIL_FILL_MAX_DAYS, ]
      if (nrow(o)) {
        dd <- abs(as.numeric(o$date - d[i]))
        hit <- data.frame(soil_temp_C = mean_na(o$soil_temp_C), soil_vwc_pct = mean_na(o$soil_vwc_pct),
                          date = o$date[order(dd, o$date)[1]])
        if (any(o$date == d[i])) hit$date <- d[i]
      }
      pre <- ""
    }
    if (is.null(hit) && !is.na(position[i]) && !is.null(pos_by[[position[i]]])) {
      hit <- pick(pos_by[[position[i]]], d[i]); pre <- "position_"
    }
    if (!is.null(hit)) {
      out$soil_temp_C[i] <- hit$soil_temp_C; out$soil_vwc_pct[i] <- hit$soil_vwc_pct
      out$soil_match[i] <- paste0(pre, if (hit$date == d[i]) "same_date" else "nearest_date")
      out$soil_reading_date[i] <- hit$date
    }
  }
  drop <- intersect(names(out), names(df))
  if (length(drop)) df <- df[, setdiff(names(df), drop), drop = FALSE]
  cbind(df, out)
}
