# soil_collars.R -- soil collars outside the modelled stand; soil readings per collar-date.
# One definition, used by tables/03_collars.R (collars.csv in_stand, which keeps them out
# of the soil model's training rows) and the scripts that report them.
# Both collars sit in the wetland margin outside the censused stand, so the soil model,
# which is applied only inside the stand, is not trained on them.
OUT_OF_STAND_COLLARS <- c("WS_1-1", "WS_1-2")

# Monthly-campaign soil readings at the collars (data/tables/soil_readings.csv), one row
# per collar and date: the mean of that collar-date's records (the iPad sheet and the
# master workbook often both hold it). plot_letter is the
# collar's landscape position (U, I, WD, WS).
collar_soil_readings <- function(path = "data/tables/soil_readings.csv") {
  sr <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "")
  sr <- sr[!is.na(sr$collar_id), ]
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  k <- paste(sr$collar_id, sr$date)
  out <- data.frame(collar_id = tapply(sr$collar_id, k, `[`, 1), date = as.Date(tapply(sr$date, k, `[`, 1)),
                    soil_temp_C = tapply(sr$soil_temp_C, k, mean_na), soil_vwc_pct = tapply(sr$soil_vwc_pct, k, mean_na),
                    n_records = as.integer(tapply(k, k, length)), row.names = NULL)
  out$plot_letter <- sub("_.*", "", out$collar_id)
  out[order(out$date, out$collar_id), ]
}
