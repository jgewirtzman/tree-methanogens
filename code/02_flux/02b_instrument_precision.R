#!/usr/bin/env Rscript
# ==============================================================================
# 02b_instrument_precision.R -- analyzer precision (sigma) per field day, from the closures
# (Methods S2, "Instrument precision"; read by every goFlux fitting script and by
# 03_precision_and_detection.R)
# ------------------------------------------------------------------------------
# Stage A; workdir: repository root. Runs before the goFlux fits, which use it.
#
# sigma is goFlux::empirical.prec(method = "hadamard"): for each closure, the median
# absolute deviation (x 1.4826) of the second differences of the gas inside the fitted
# window, divided by sqrt(6); a field day's sigma is the median over its closures. A
# second difference cancels any straight line, so the flux does not enter it, and the
# median over closures keeps a few curved or disturbed closures from setting the value.
# Differences never span a change of logging interval (5 s for the 2021 survey, 1 s
# otherwise); a day logged at two intervals gets one row per interval.
#
# One LGR analyzer made every measurement, so a group is a field day x logging interval,
# pooling the stem and soil closures of that day. Computed for CH4 (ppb) and CO2 (ppm).
#
# Inputs:  data/raw/flux_windows/lgr_manual_identification_{2021_multiheight_final,
#            2023_cross_species,semirigid_tree,results_soil}.csv
#          data/raw/flux_windows/untagged/untagged_manID.rds
# Outputs: data/processed/flux/analyzer_sigma_by_day.csv  day, dt_s, closures, sigma
#            (read by the fits and by 03_precision_and_detection.R)
#          outputs/data/instrument_precision.csv  by campaign: median and range of the
#            daily sigma (the values quoted in Methods S2)
# ==============================================================================
source("code/lib/outputs.R")
source("code/lib/flux_windows.R")
suppressMessages({ library(goFlux); library(dplyr) })

WIN <- c("2021_multiheight" = "lgr_manual_identification_2021_multiheight_final.csv",
         "2023_cross_species" = "lgr_manual_identification_2023_cross_species.csv",
         "monthly_stem" = "lgr_manual_identification_semirigid_tree.csv",
         "monthly_soil" = "lgr_manual_identification_results_soil.csv")
COLS <- c("UniqueID", "flag", "POSIX.time", "CH4dry_ppb", "CO2dry_ppm")
W <- bind_rows(
  lapply(names(WIN), function(k) load_windows(window_path(WIN[[k]]))[, COLS] %>% mutate(set = k)),
  readRDS(window_path("untagged/untagged_manID.rds"))[, COLS] %>%
    mutate(POSIX.time = as.POSIXct(POSIX.time, tz = "UTC"), set = "monthly_stem_untagged"))
W$day <- format(W$POSIX.time, "%Y-%m-%d", tz = "UTC")
stopifnot(!anyNA(W$day))
# repeated whole-second timestamps (an analyzer logging slightly faster than 1 Hz) end a run of
# differences in empirical.prec; reported, not removed
cat(sprintf("rows with a repeated timestamp within their closure: %d of %d\n",
            sum(duplicated(W[c("set", "UniqueID", "POSIX.time")])), nrow(W)))
W$UniqueID <- paste(W$set, W$UniqueID)          # an id is unique within its window file

prec_by_day <- function(gas) {
  p <- suppressWarnings(empirical.prec(W, gas, method = "hadamard", by = "day", warn = FALSE))
  # a closure with fewer than 6 second differences has no sigma of its own (it takes its day's);
  # empirical.prec lists such closures as a group row with no interval
  p <- p[is.finite(p$dt_s) & p$n.closures > 0, ]
  p[, c("day", "dt_s", "n.closures", "prec", "ac1", "zero.frac")]
}
CH4 <- prec_by_day("CH4dry_ppb"); CO2 <- prec_by_day("CO2dry_ppm")
S <- full_join(CH4 %>% rename(closures = n.closures, sigma_CH4_ppb = prec, ac1_CH4 = ac1, zero_frac_CH4 = zero.frac),
               CO2 %>% select(day, dt_s, sigma_CO2_ppm = prec), by = c("day", "dt_s")) %>% arrange(day, dt_s)
# the interval a day's closures were logged at (one per day here)
stopifnot(!anyDuplicated(S$day), all(is.finite(S$sigma_CH4_ppb)), all(is.finite(S$sigma_CO2_ppm)),
          all(S$sigma_CH4_ppb > 0), all(S$sigma_CO2_ppm > 0))
write.csv(S, "data/processed/flux/analyzer_sigma_by_day.csv", row.names = FALSE)
cat(sprintf("sigma per field day: %d days, %d closures; CH4 median %.3f ppb (range %.3f-%.3f)\n",
            nrow(S), sum(S$closures), median(S$sigma_CH4_ppb), min(S$sigma_CH4_ppb), max(S$sigma_CH4_ppb)))
if (any(S$ac1_CH4 > -0.5)) cat("days with red noise (lag-1 autocorrelation of second differences > -0.5):",
                               paste(S$day[S$ac1_CH4 > -0.5], collapse = ", "), "\n")

# campaign summary for the text (the 2021 and 2023 surveys; the monthly survey pools stems and soil)
camp_days <- W %>% distinct(set, day) %>%
  mutate(campaign = recode(set, "2021_multiheight" = "Height+molecular", "2023_cross_species" = "Cross-species",
                           .default = "Monthly survey")) %>% distinct(campaign, day)
C <- camp_days %>% inner_join(S, by = "day") %>% group_by(campaign) %>%
  summarise(days = n(), closures = sum(closures), logging_interval_s = paste(sort(unique(dt_s)), collapse = "/"),
            sigma_ppb = round(median(sigma_CH4_ppb), 3), sigma_min_ppb = round(min(sigma_CH4_ppb), 3),
            sigma_max_ppb = round(max(sigma_CH4_ppb), 3), .groups = "drop")
write.csv(C, out_path("instrument_precision.csv"), row.names = FALSE)
print(as.data.frame(C))
