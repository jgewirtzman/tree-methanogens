#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 02b_instrument_precision.R -- analyzer CH4 precision (sigma) per field period
# (Methods S2, "Instrument precision"; read by 03_precision_and_detection.R)
# ------------------------------------------------------------------------------
# sigma = MAD(dx) / sqrt(2), where dx is the difference between consecutive logged
# CH4 (dry, ppb) readings over each field period's ENTIRE raw analyzer record.
#   - R's mad() already applies the 1.4826 normal-consistency factor.
#   - A constant rate of change shifts every difference equally, so a flux cancels
#     in the median-based spread; moves, breaths and chamber changes are isolated
#     large differences that the MAD tolerates.
#   - Differences are taken only between consecutive readings on the same day and
#     no more than 1.5 logging intervals apart (no gaps, no day boundaries).
#   - Logging interval differs by period: 5 s for the 2021 survey, 1 s otherwise.
# Field periods are the raw folders in data/raw/lgr/. The 2023 folder also holds
# the felled-oak day (2022-10-10), which belongs to no survey and is excluded.
#
# goFlux::import2RData(merge = TRUE) merges EVERY file in the working directory's
# RData/ folder, so each period is imported in its own empty temporary directory;
# importing from the repository root would pool all periods.
#
# Supersedes the constants hard-coded in 03_precision_and_detection.R before
# 2026-10-09 (1.200, 1.725, 2.181 ppb), whose generating code was never committed.
# Output: outputs/data/instrument_precision.csv
# ==============================================================================
suppressMessages({ library(goFlux); library(dplyr) })

PERIODS <- data.frame(
  camp   = c("Height+molecular", "Cross-species", "Monthly survey"),
  folder = c("static_2021", "static_2023", "semirigid_2020-2021"),
  drop_days = c("", "2022-10-10", ""),
  stringsAsFactors = FALSE)

root <- getwd()
sigma_for <- function(camp, folder, drop_days) {
  src <- file.path(root, "data/raw/lgr", folder)
  work <- file.path(tempdir(), paste0("prec_", folder)); unlink(work, recursive = TRUE)
  dir.create(file.path(work, "in"), recursive = TRUE)
  old <- setwd(work); on.exit(setwd(old))
  for (z in list.files(src, "\\.zip$", recursive = TRUE, full.names = TRUE)) unzip(z, exdir = "in", junkpaths = TRUE)
  for (f in list.files(src, "\\.txt$", recursive = TRUE, full.names = TRUE)) file.copy(f, "in")
  d <- suppressMessages(suppressWarnings(import2RData(path = "in", instrument = "UGGA", date.format = "mdy",
         timezone = "UTC", keep_all = FALSE, prec = c(0.35, 0.9, 200), merge = TRUE)))
  d <- d %>% arrange(POSIX.time) %>% distinct(POSIX.time, .keep_all = TRUE) %>%
    mutate(day = as.Date(POSIX.time))
  if (nzchar(drop_days)) d <- d %>% filter(!day %in% as.Date(strsplit(drop_days, ";")[[1]]))
  dt <- c(NA, diff(as.numeric(d$POSIX.time)))
  dx <- c(NA, diff(d$CH4dry_ppb))
  same_day <- c(FALSE, diff(as.numeric(d$day)) == 0)
  step <- median(dt, na.rm = TRUE)
  ok <- is.finite(dx) & same_day & dt > 0 & dt <= 1.5 * step
  data.frame(camp = camp, folder = folder, days = length(unique(d$day)), readings = nrow(d),
             logging_interval_s = round(step, 2), differences = sum(ok),
             sigma_ppb = round(mad(dx[ok]) / sqrt(2), 3))
}

res <- do.call(rbind, lapply(seq_len(nrow(PERIODS)), function(i)
  sigma_for(PERIODS$camp[i], PERIODS$folder[i], PERIODS$drop_days[i])))
write.csv(res, out_path("instrument_precision.csv"), row.names = FALSE)
print(res)
