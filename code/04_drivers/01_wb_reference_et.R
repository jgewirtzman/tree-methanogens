#!/usr/bin/env Rscript
# ==============================================================================
# 01_wb_reference_et.R
# ------------------------------------------------------------------------------
# Daily reference evapotranspiration and a soil water-balance index for the
# Yale-Myers tower record, recomputed from raw meteorology.
#
# WHY NOT USE THE LOGGER'S ETos COLUMN. It is unusable: the median is a plausible
# 0.20 mm h-1 but the tail runs to 334 mm h-1, and capping at a physical ceiling
# of 1.5 mm h-1 discards 56.8% of the record. A balance built on it is missing
# more than half its evaporative demand, which is why a bucket driven by it fills
# straight through summer and ends up NEGATIVELY correlated with observed soil
# moisture at monthly scale (-0.52).
#
# The inputs needed to compute ET properly are intact: air temperature (mean, max,
# min), RH and solar radiation are each 100% complete across 2018-09 to 2023-05.
# Wind is not (see below).
#
# Input: data/tables/weather.csv (true time; screened by code/tables/06_weather.R).
# Days are local calendar days of the true time.
#
# GAPS, one rule. A day is used when it has at least 20 hourly records; its rain is
# the sum of its hourly precipitation, and an hour with no precipitation value would
# count as 0 (the build reports how many; none at present). A day with fewer records
# is missing: the bucket and the API step over it with no rain and no ET (the store
# carries unchanged). The build reports the missing days (2 single days at present).
#
# WHAT THIS IS FOR. The soil campaign ran one year (2020-06 to 2021-05), and that
# year had the second-driest Jun-Sep in the record (210 mm against 437 mm in
# 2021). The RF learned a moisture x temperature interaction, so driving the
# seasonal upscaling with that one year's pairing bakes a drought into every
# month. A water-balance index computed for all years gives a temporal driver
# that can be averaged into a climatology. The tower VWC probe cannot do this
# job: it failed 2020-06 to 2020-09, which is 15 of the 22 campaign dates.
#
# Output: outputs/tables/tower_daily_waterbalance.csv
#         outputs/audit/wb_reference_et.txt
# ==============================================================================
suppressPackageStartupMessages({library(dplyr)})

source("code/lib/fao56.R")      # et0_fao56(): the one FAO-56 Penman-Monteith (u2 = FAO56_U2)

w <- read.csv("data/tables/weather.csv", stringsAsFactors = FALSE)
w$date <- as.Date(substr(w$timestamp, 1, 10))
NUM <- function(v) suppressWarnings(as.numeric(w[[v]]))
w$tair <- NUM("air_temp_avg_C"); w$tmax <- NUM("air_temp_max_C"); w$tmin <- NUM("air_temp_min_C")
w$rh <- NUM("relative_humidity_pct")
# WIND IS UNUSABLE. The logger's WindSpeed_ms_Avg has a median of 30.2 and a max of
# 289 despite the name; those are not m s-1 and the true unit is not recorded (after
# the table's screen the median is still 7). Fed into the Penman-Monteith
# aerodynamic term it inflated ET0 to 959 mm yr-1. FAO-56 section on missing data
# prescribes a 2 m s-1 default in exactly this case, which is what is used
# (FAO56_U2, code/lib/fao56.R). The
# weather table now leaves wind speed blank (06_weather.R).
# Solar: the logger column named Solar_Rad_kW_Avg is W m-2 (solar_rad_avg_W_m2), and
# solar_rad_total_MJ_m2 is a direct hourly MJ m-2 total. Summed per day the two
# agree at 12.1 MJ m-2 d-1, so the MJ column is used and the other is the check.
w$srad_mj <- NUM("solar_rad_total_MJ_m2")
w$srad_w  <- NUM("solar_rad_avg_W_m2")
w$rain <- NUM("precip_mm"); w$vwc <- NUM("vwc")

D <- w %>% group_by(date) %>% summarise(
  Tmean = mean(tair, na.rm = TRUE), Tmax = max(tmax, na.rm = TRUE), Tmin = min(tmin, na.rm = TRUE),
  RH    = mean(rh,   na.rm = TRUE),
  Rs      = sum(srad_mj, na.rm = TRUE),                    # MJ m-2 day-1, direct
  Rs_chk  = mean(srad_w, na.rm = TRUE) * 86400 / 1e6,      # from W m-2, as a check
  rain_mm = sum(rain, na.rm = TRUE),
  vwc   = mean(ifelse(vwc > 0 & vwc < 0.6, vwc, NA), na.rm = TRUE),
  nobs  = dplyr::n(), .groups = "drop") %>%
  mutate(across(c(Tmean, Tmax, Tmin, RH, Rs, Rs_chk, vwc), ~ifelse(is.finite(.x), .x, NA))) %>%
  filter(nobs >= 20)                                   # complete days only
rain_na <- sum(tapply(is.na(w$rain), w$date, sum)[as.character(D$date)])
cal <- seq(min(w$date, na.rm = TRUE), max(w$date, na.rm = TRUE), by = "day")
gap_days <- cal[!cal %in% D$date]
cat(sprintf("days used %d of %d; missing days (fewer than 20 hourly records; the store carries over them): %d (%s)\n",
            nrow(D), length(cal), length(gap_days), paste(gap_days, collapse = ", ")))
cat(sprintf("hours with no precipitation value on days used (counted as 0): %d\n", rain_na))

# --- FAO-56 Penman-Monteith (code/lib/fao56.R); Hargreaves-Samani as a check ----
esat <- function(t) 0.6108*exp(17.27*t/(t + 237.3))
ea <- (esat(D$Tmax) + esat(D$Tmin))/2 * D$RH/100
E <- et0_fao56(D$date, D$Tmax, D$Tmin, D$Tmean, D$Rs, ea)
D$ET0_pm <- E$ET0_pm; D$ET0_hs <- E$ET0_hs

cat(sprintf("solar cross-check: MJ column %.2f vs W-derived %.2f MJ m-2 d-1 (r = %.4f)\n",
            mean(D$Rs, na.rm=TRUE), mean(D$Rs_chk, na.rm=TRUE),
            cor(D$Rs, D$Rs_chk, use="complete.obs")))
cat("=== RECOMPUTED REFERENCE ET (mm day-1) ===\n")
for (v in c("ET0_pm", "ET0_hs")) cat(sprintf("  %-8s median %.2f  q95 %.2f  max %.2f  annual %.0f mm\n",
  v, median(D[[v]], na.rm=TRUE), quantile(D[[v]], .95, na.rm=TRUE),
  max(D[[v]], na.rm=TRUE), sum(D[[v]], na.rm=TRUE)/(nrow(D)/365.25)))
cat(sprintf("  PM vs HS correlation: %+.3f\n", cor(D$ET0_pm, D$ET0_hs, use="complete.obs")))
cat(sprintf("  wind speed not recorded usably; u2 fixed at %.1f m s-1 per FAO-56\n", FAO56_U2))
cat("  (FAO-56 reference ET for this region should sit near 600-750 mm yr-1)\n")

# --- water balance ------------------------------------------------------------
FAULT <- c("2020-06","2020-07","2020-08","2020-09")   # VWC probe failure
D$ym <- format(D$date, "%Y-%m")
D$vwc_ok <- !(D$ym %in% FAULT) & is.finite(D$vwc)
cat(sprintf("\ndays with usable tower VWC: %d of %d\n", sum(D$vwc_ok), nrow(D)))

bucket <- function(p, e, smax) { out <- numeric(length(p)); s <- smax/2
  for (i in seq_along(p)) { s <- min(max(s + ifelse(is.finite(p[i]),p[i],0) -
      ifelse(is.finite(e[i]),e[i],0), 0), smax); out[i] <- s }; out }
api <- function(p, k) { out <- numeric(length(p)); a <- 0
  for (i in seq_along(p)) { a <- k*a + ifelse(is.finite(p[i]), p[i], 0); out[i] <- a }; out }

cand <- list()
for (sm in c(50, 75, 100, 150, 200, 300))
  cand[[sprintf("bucket_PM_%d", sm)]] <- bucket(D$rain_mm, D$ET0_pm, sm)
for (sm in c(75, 100, 150))
  cand[[sprintf("bucket_HS_%d", sm)]] <- bucket(D$rain_mm, D$ET0_hs, sm)
for (k in c(0.80, 0.90, 0.95)) cand[[sprintf("API_%.2f", k)]] <- api(D$rain_mm, k)

ok <- D$vwc_ok
R <- bind_rows(lapply(names(cand), function(n) data.frame(index = n,
  pearson = cor(cand[[n]][ok], D$vwc[ok]),
  spearman = cor(cand[[n]][ok], D$vwc[ok], method = "spearman"))))
mo <- as.integer(format(D$date, "%m"))
R$monthly <- sapply(names(cand), function(n) {
  a <- tapply(cand[[n]], mo, mean); b <- tapply(D$vwc[ok], mo[ok], mean)
  cor(a[names(b)], b, use = "complete.obs") })
R <- R %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>% arrange(-spearman)
cat("\n=== INDEX vs CLEAN TOWER VWC ===\n"); print(as.data.frame(R), row.names = FALSE)
BEST <- R$index[1]
cat(sprintf("\nbest: %s (Spearman %.3f, monthly %.3f)\n", BEST, R$spearman[1], R$monthly[1]))
cat(sprintf("previous best (rain-only API k=0.80, broken ETos): Spearman 0.290, monthly ~0.60\n"))

D$wb_index <- cand[[BEST]]
for (n in names(cand)) D[[n]] <- cand[[n]]
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(D, "outputs/tables/tower_daily_waterbalance.csv", row.names = FALSE)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
con <- file("outputs/audit/wb_reference_et.txt", "w")
cat(sprintf("REFERENCE ET AND WATER BALANCE (01_wb_reference_et.R)\nbuilt %s\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = con)
cat(sprintf("logger ETos column NOT used: 56.8%% of values outside [0, 1.5] mm/h\n"), file = con)
cat(sprintf("PM annual %.0f mm | HS annual %.0f mm | r = %.3f\n",
            sum(D$ET0_pm, na.rm=TRUE)/(nrow(D)/365.25),
            sum(D$ET0_hs, na.rm=TRUE)/(nrow(D)/365.25),
            cor(D$ET0_pm, D$ET0_hs, use="complete.obs")), file = con)
cat(sprintf("\nbest index: %s\n\n", BEST), file = con)
utils::write.table(R, con, sep = "\t", row.names = FALSE, quote = FALSE)
close(con)
cat("\nwritten: outputs/tables/tower_daily_waterbalance.csv, outputs/audit/wb_reference_et.txt\n")
