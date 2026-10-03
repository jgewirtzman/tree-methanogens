# monthly_cols.R -- one column convention for the monthly-survey flux files
# ------------------------------------------------------------------------------
# The stored monthly stem files carried every goFlux CH4 column twice, as .x and .y
# (identical values; an old merge), and readers asked for the .x copy. A fresh run
# from raw writes plain names, so those readers silently got NA for every tagged
# stem (model prep saw 45 monthly measurements instead of 414). monthly_plain()
# maps either form to plain names; every reader of these files calls it.
monthly_plain <- function(df) {
  nm <- names(df)
  sfx <- grep("^(CH4|CO2)_.*\\.[xy]$", nm, value = TRUE)
  for (n in grep("\\.x$", sfx, value = TRUE)) {
    base <- sub("\\.x$", "", n)
    df[[base]] <- if (base %in% names(df)) dplyr::coalesce(df[[base]], df[[n]]) else df[[n]]
  }
  for (n in grep("\\.y$", sfx, value = TRUE)) {
    base <- sub("\\.y$", "", n)
    if (!base %in% names(df)) df[[base]] <- df[[n]] else df[[base]] <- dplyr::coalesce(df[[base]], df[[n]])
  }
  df[, setdiff(names(df), sfx), drop = FALSE]
}
