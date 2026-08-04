#!/usr/bin/env Rscript
# ==============================================================================
# wb_reference_et.R
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
# The inputs needed to compute ET properly are all intact: Tair_Avg, Tair_Max,
# Tair_Min, RH, Solar_Rad_kW_Avg and WindSpeed_ms_Avg are each 100% complete
# across 2018-09 to 2023-05.
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

LAT_DEG <- 41.9899      # plot centre
ELEV_M  <- 215          # Yale-Myers, approximate
ALBEDO  <- 0.23         # FAO-56 reference grass
SIGMA   <- 4.903e-9     # Stefan-Boltzmann, MJ K-4 m-2 day-1

w <- read.csv("data/raw/weather/ymf_clean_sorted.csv", stringsAsFactors = FALSE)
w$date <- as.Date(substr(w$TIMESTAMP, 1, 10))
NUM <- function(v) suppressWarnings(as.numeric(w[[v]]))
w$tair <- NUM("Tair_Avg"); w$tmax <- NUM("Tair_Max"); w$tmin <- NUM("Tair_Min")
w$rh <- NUM("RH")
# WIND IS UNUSABLE. WindSpeed_ms_Avg has a median of 30.2 and a max of 289
# despite the name; those are not m s-1 and the true unit is not recorded. Fed
# into the Penman-Monteith aerodynamic term it inflated ET0 to 959 mm yr-1.
# FAO-56 section on missing data prescribes a 2 m s-1 default in exactly this
# case, which is what is used; the column is retained only for the record.
U2_DEFAULT <- 2.0
w$wind_raw <- NUM("WindSpeed_ms_Avg")
# Solar: the column named Solar_Rad_kW_Avg is actually W m-2 (mean 140, max 1237),
# and Solar_Rad_Tot_MJ_Tot is a direct hourly MJ m-2 total. Summed per day the two
# agree at 12.1 MJ m-2 d-1, so the MJ column is used and the other is the check.
w$srad_mj <- NUM("Solar_Rad_Tot_MJ_Tot")
w$srad_w  <- NUM("Solar_Rad_kW_Avg")
w$rain <- NUM("Rain_mm_Tot"); w$vwc <- NUM("VWC")

D <- w %>% group_by(date) %>% summarise(
  Tmean = mean(tair, na.rm = TRUE), Tmax = max(tmax, na.rm = TRUE), Tmin = min(tmin, na.rm = TRUE),
  RH    = mean(rh,   na.rm = TRUE),
  Rs      = sum(srad_mj, na.rm = TRUE),                    # MJ m-2 day-1, direct
  Rs_chk  = mean(srad_w, na.rm = TRUE) * 86400 / 1e6,      # from W m-2, as a check
  u_raw = mean(wind_raw, na.rm = TRUE),
  rain_mm = sum(rain, na.rm = TRUE),
  vwc   = mean(ifelse(vwc > 0 & vwc < 0.6, vwc, NA), na.rm = TRUE),
  nobs  = dplyr::n(), .groups = "drop") %>%
  mutate(across(c(Tmean, Tmax, Tmin, RH, Rs, Rs_chk, vwc), ~ifelse(is.finite(.x), .x, NA))) %>%
  filter(nobs >= 20)                                   # complete days only

# --- FAO-56 Penman-Monteith ---------------------------------------------------
doy <- as.integer(format(D$date, "%j"))
phi <- LAT_DEG * pi/180
dr  <- 1 + 0.033*cos(2*pi*doy/365)
dec <- 0.409*sin(2*pi*doy/365 - 1.39)
ws  <- acos(pmin(pmax(-tan(phi)*tan(dec), -1), 1))
Ra  <- (24*60/pi) * 0.0820 * dr * (ws*sin(phi)*sin(dec) + cos(phi)*cos(dec)*sin(ws))
Rso <- (0.75 + 2e-5*ELEV_M) * Ra
esat <- function(t) 0.6108*exp(17.27*t/(t + 237.3))
es <- (esat(D$Tmax) + esat(D$Tmin))/2
ea <- es * D$RH/100
delta <- 4098*esat(D$Tmean)/(D$Tmean + 237.3)^2
P <- 101.3*((293 - 0.0065*ELEV_M)/293)^5.26
gamma <- 0.665e-3 * P
Rns <- (1 - ALBEDO) * D$Rs
frac <- pmin(pmax(D$Rs/Rso, 0.25), 1.0)
Rnl <- SIGMA * (((D$Tmax + 273.16)^4 + (D$Tmin + 273.16)^4)/2) *
       (0.34 - 0.14*sqrt(pmax(ea, 0))) * (1.35*frac - 0.35)
Rn <- Rns - Rnl
u2 <- rep(U2_DEFAULT, nrow(D))    # see the wind note above
D$ET0_pm <- pmax((0.408*delta*Rn + gamma*900/(D$Tmean + 273)*u2*(es - ea)) /
                 (delta + gamma*(1 + 0.34*u2)), 0)
# --- Hargreaves-Samani, as an independent check --------------------------------
D$ET0_hs <- pmax(0.0023 * (Ra/2.45) * (D$Tmean + 17.8) * sqrt(pmax(D$Tmax - D$Tmin, 0)), 0)

cat(sprintf("solar cross-check: MJ column %.2f vs W-derived %.2f MJ m-2 d-1 (r = %.4f)\n",
            mean(D$Rs, na.rm=TRUE), mean(D$Rs_chk, na.rm=TRUE),
            cor(D$Rs, D$Rs_chk, use="complete.obs")))
cat("=== RECOMPUTED REFERENCE ET (mm day-1) ===\n")
for (v in c("ET0_pm", "ET0_hs")) cat(sprintf("  %-8s median %.2f  q95 %.2f  max %.2f  annual %.0f mm\n",
  v, median(D[[v]], na.rm=TRUE), quantile(D[[v]], .95, na.rm=TRUE),
  max(D[[v]], na.rm=TRUE), sum(D[[v]], na.rm=TRUE)/(nrow(D)/365.25)))
cat(sprintf("  PM vs HS correlation: %+.3f\n", cor(D$ET0_pm, D$ET0_hs, use="complete.obs")))
cat(sprintf("  wind column unusable (median %.1f, max %.1f); u2 fixed at %.1f m s-1 per FAO-56\n",
            median(D$u_raw, na.rm=TRUE), max(D$u_raw, na.rm=TRUE), U2_DEFAULT))
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
cat(sprintf("REFERENCE ET AND WATER BALANCE (wb_reference_et.R)\nbuilt %s\n\n",
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
