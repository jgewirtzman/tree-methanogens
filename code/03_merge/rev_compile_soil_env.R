#!/usr/bin/env Rscript
# ==============================================================================
# rev_compile_soil_env.R
# ------------------------------------------------------------------------------
# Compiles every per-collar soil temperature / moisture sheet in the iPad field
# archive into one tidy table.
#
# WHY
#   The scaling pipeline drew soil temp + VWC from soilmoisture_total.csv, which keys
#   on PLOT rather than collar and stops at Feb 2021 -- so soil temperature had a
#   single value per month, and the March/April/May 2021 campaigns (88 flux
#   measurements) had no measured covariates at all and fell back to monthly medians.
#   The underlying per-date sheets exist for every campaign; they were simply never
#   compiled.
#
# THE FILES ARE MESSY, deliberately handled here rather than by hand:
#   * Numbers exports: real data sits on "Sheet 1 - Table 1", behind an export notice,
#     with a title row above the header row.
#   * Collar is encoded as "U-1-3" (letter-plot-subplot) in the 2021 sheets but split
#     across Plot / Plot letter columns in the 2020 sheets.
#   * The volumetric moisture PERCENT lands in a column called "Soil moisture" in the
#     2020 sheets and "VWC" in the 2021 sheets (where "VWC" in 2020 means the raw
#     sensor period, ~1.4). Chosen here by plausible range, not by column name.
#   * Dates are sometimes blank in-sheet; recovered from the filename.
#
# Output: data/processed/environmental/soil_env_by_collar.csv
#         columns: Date, plot_letter, plot_tag, soil_temp_C, soil_moisture_pct, source
# ==============================================================================

suppressMessages({library(readxl); library(dplyr)})

ROOT <- "."
files <- list.files(file.path(ROOT, "data/raw/field_data/ipad_data"),
                    pattern = "^Soil temp moisture.*\\.xlsx$",
                    recursive = TRUE, full.names = TRUE)
cat("Found", length(files), "per-date soil env sheets\n\n")

# a plausible volumetric water content, in percent
is_pct <- function(v) { v <- suppressWarnings(as.numeric(v)); !is.na(v) & v > 1.6 & v < 80 }

parse_one <- function(f) {
  sheets <- tryCatch(excel_sheets(f), error = function(e) character(0))
  for (s in sheets) {
    d <- suppressMessages(tryCatch(read_excel(f, sheet = s, col_names = FALSE),
                                   error = function(e) NULL))
    if (is.null(d) || nrow(d) < 3 || ncol(d) < 4) next
    d <- as.data.frame(d, stringsAsFactors = FALSE)

    # locate the header row by looking for "Soil temp"
    hdr <- which(apply(d, 1, function(r) any(grepl("soil\\s*temp", r, ignore.case = TRUE))))[1]
    if (is.na(hdr)) next
    nm <- tolower(trimws(as.character(unlist(d[hdr, ]))))
    body <- d[(hdr + 1):nrow(d), , drop = FALSE]
    names(body) <- ifelse(is.na(nm) | nm == "", paste0("v", seq_along(nm)), nm)

    col <- function(pat) { i <- grep(pat, names(body)); if (length(i)) body[[i[1]]] else NA }
    plot_raw  <- col("^plot$")
    letter_raw<- col("plot\\s*letter")
    temp_raw  <- col("soil\\s*temp")
    moist_a   <- col("soil\\s*moisture")
    moist_b   <- col("^vwc")

    if (all(is.na(plot_raw))) next

    # moisture percent: prefer whichever column is in a plausible VWC-% range
    pick_moisture <- function(a, b) {
      a <- suppressWarnings(as.numeric(a)); b <- suppressWarnings(as.numeric(b))
      score <- function(v) if (all(is.na(v))) -1 else mean(is_pct(v), na.rm = TRUE)
      if (score(a) >= score(b)) a else b
    }
    moist <- pick_moisture(moist_a, moist_b)

    # Collar identity. Two encodings, plus one corruption:
    #   2021 sheets : "U-1-3"  (letter-plot-subplot in one cell)
    #   2020 sheets : Plot + "Plot letter" in separate columns -- BUT Numbers silently
    #                 coerced collar labels like "1-3" into DATES, so the Plot column
    #                 arrives as a serial number. These use the Mac 1904 epoch, under
    #                 which 42369..42374 decode cleanly to 1-1..1-6 (the 1899 epoch
    #                 shifts everything by one and yields "12-31").
    pr <- trimws(as.character(plot_raw))
    combined <- grepl("^[A-Za-z]{1,2}-\\d+-\\d+$", pr)
    serial <- suppressWarnings(as.numeric(pr))
    is_ser <- !is.na(serial) & serial > 30000 & serial < 60000
    dec <- as.Date(serial, origin = "1904-01-01")
    pr_dec <- ifelse(is_ser,
                     paste0(as.integer(format(dec, "%m")), "-", as.integer(format(dec, "%d"))),
                     pr)
    letter <- ifelse(combined, toupper(sub("^([A-Za-z]{1,2})-.*", "\\1", pr)),
                     toupper(trimws(as.character(letter_raw))))
    tag <- ifelse(combined, sub("^[A-Za-z]{1,2}-", "", pr), pr_dec)

    # date: in-sheet where present (Excel serial or text), else from the filename
    fdate <- as.Date(sub(".*?(\\d{8}).*", "\\1", basename(f)), format = "%Y%m%d")
    # In-sheet dates are corrupted by the same Numbers coercion and disagree with the
    # campaign by years, so the filename date (YYYYMMDD) is authoritative.
    Date <- fdate

    out <- data.frame(
      Date = Date, plot_letter = letter, plot_tag = tag,
      soil_temp_C = suppressWarnings(as.numeric(temp_raw)),
      soil_moisture_pct = moist,
      source = basename(f), stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(plot_letter), plot_letter %in% c("U","I","WD","WS"),
             grepl("^\\d+-\\d+$", plot_tag),
             !(is.na(soil_temp_C) & is.na(soil_moisture_pct)))
    if (nrow(out) > 0) return(out)
  }
  NULL
}

all_env <- bind_rows(lapply(files, function(f) {
  r <- tryCatch(parse_one(f), error = function(e) { cat("  ! failed:", basename(f), "\n"); NULL })
  if (!is.null(r)) cat(sprintf("  %-46s %3d rows  %s\n", basename(f), nrow(r), format(r$Date[1])))
  r
}))

all_env <- all_env %>%
  filter(!is.na(Date)) %>%
  group_by(Date, plot_letter, plot_tag) %>%
  summarise(soil_temp_C = mean(soil_temp_C, na.rm = TRUE),
            soil_moisture_pct = mean(soil_moisture_pct, na.rm = TRUE),
            source = paste(unique(source), collapse = ";"), .groups = "drop") %>%
  mutate(across(c(soil_temp_C, soil_moisture_pct), ~ifelse(is.nan(.x), NA_real_, .x)))

cat(sprintf("\nCompiled: %d collar-date records across %d dates\n",
            nrow(all_env), length(unique(all_env$Date))))
cat(sprintf("  soil temp non-NA: %d   moisture non-NA: %d\n",
            sum(!is.na(all_env$soil_temp_C)), sum(!is.na(all_env$soil_moisture_pct))))

dir.create(file.path(ROOT, "data/processed/environmental"), showWarnings = FALSE, recursive = TRUE)
outf <- file.path(ROOT, "data/processed/environmental/soil_env_by_collar.csv")
write.csv(all_env, outf, row.names = FALSE)
cat("Written:", outf, "\n\n")

# --- coverage against the soil flux records ---------------------------------
sf <- read.csv(file.path(ROOT, "data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv"),
               stringsAsFactors = FALSE)
sf <- sf[!is.na(sf$CH4_best.flux), ]
sf$k <- paste(as.Date(sf$Date), toupper(trimws(sf$Plot.letter)), trimws(sf$Plot.Tag))
all_env$k <- paste(all_env$Date, all_env$plot_letter, all_env$plot_tag)
cat("=== coverage of soil flux measurements ===\n")
cat(sprintf("  flux records: %d   with compiled env: %d (%.0f%%)\n",
            nrow(sf), sum(sf$k %in% all_env$k), 100 * mean(sf$k %in% all_env$k)))
miss <- sort(unique(as.character(as.Date(sf$Date[!sf$k %in% all_env$k]))))
if (length(miss)) cat("  still missing dates:", paste(miss, collapse = " "), "\n")
