#!/usr/bin/env Rscript
# ==============================================================================
# rev_moisture_climatology.R
# ------------------------------------------------------------------------------
# Monthly plot-wide soil moisture, averaged over 2019-2022 instead of tied to the
# single campaign year. Feeds rev_predict_soil_surface.R.
#
# THE PROBLEM. The soil campaign ran one year (2020-06 to 2021-05) and that year
# had the second-driest Jun-Sep in the record (210 mm against 437 mm in 2021).
# The SoilRF learned a moisture x temperature interaction, so driving the
# seasonal upscaling from that one year's pairing bakes a drought into every
# month. The previous route -- a per-month affine rescaling of the December
# transect field, with alpha and beta fitted separately per month on 3-28
# observations at R2 = 0.02-0.31 -- also produced a curve that is not a
# climatology at all: it jumps from 0.113 in November to 0.355 in December.
#
# THE METHOD, in three steps.
#
# 1. DE-CONFOUND THE HANDHELD. Which collars were visited is confounded with the
#    plot-wide level: sample only wetland collars and the day looks wet. Across
#    28 collars x 22 dates, 71.7% of moisture variance is WHICH COLLAR and only
#    16.4% is WHEN. Fitting site + date and reading off the date effects removes
#    the sampling bias; it shrinks apparent day-to-day variation by 38% (SD 0.107
#    -> 0.067). The collar wet/dry ranking is stable enough to justify this
#    (median pairwise Spearman +0.72 across dates).
#
# 2. TIE THAT TO A TOWER WATER BALANCE. A bucket driven by rainfall and the
#    recomputed Penman-Monteith ET (see rev_wb_reference_et.R -- the logger's own
#    ETos, solar and wind columns are all unusable) is regressed on the 22
#    de-confounded campaign levels, then applied to 2019-2022.
#
#    THE INDEX IS SELECTED ON SEASONAL SHAPE, NOT DAILY CORRELATION. Chosen by
#    daily agreement with the tower VWC probe, a short-memory index wins (API
#    k=0.80, 5-day time constant) -- but that probe is shallow and fast, and the
#    resulting climatology puts the annual maximum in December when observations
#    put it in March. Scoring instead on how well the back-predicted monthly
#    curve matches the observed one:
#
#      index                        LOO R2   shape r   peak month (obs = Mar)
#      bucket Smax >= 400 mm         0.789     0.898    Feb
#      bucket Smax = 200 mm          0.676     0.870    Dec
#      API k=0.99 (tau 100 d)        0.551     0.639    Dec
#      API k=0.80 (tau 5 d)          0.334     0.315    Dec
#
#    A large store is what carries seasonal recharge: autumn rain accumulates
#    through a winter when ET is near zero and the profile peaks in early spring.
#    A five-day memory cannot produce that lag.
#
#    WHAT THE STORE IS. A one-line running balance in millimetres of water:
#        S_t = clamp( S_{t-1} + rain_t - ET_t , 0 , Smax )
#    started at Smax/2 and stepped daily. It is not soil moisture; it is a
#    catchment-style water store whose seasonal swing is regressed onto the
#    de-confounded plot moisture, so the regression sets both the units and the
#    absolute level. Rainfall averages 1105 mm yr-1 against 694 mm of reference
#    ET, a surplus of +410 mm yr-1, so the store refills every winter.
#
#    Smax therefore acts as FIELD CAPACITY, not as a free parameter: the store
#    sits at the ceiling through winter and spring (10% of all days, 5% of days
#    in the campaign window) and the seasonal signal is the summer drawdown --
#    to 105 mm in 2020 and 25 mm in 2022, against 225-285 mm in the wetter 2019
#    and 2021. Raising Smax above ~300 mm therefore changes only the winter
#    plateau, which the transfer function is flat across; the Jun-Sep
#    climatology is identical (0.2443) for every Smax from 300 to 1200 mm and
#    the drought correction stays at +0.031. That insensitivity is the reason
#    for choosing 400 rather than tuning it.
#
# 3. SCALE THE SPATIAL FIELD TO THE MONTHLY LEVEL. The relative spatial pattern
#    is taken as fixed in time and comes from the dense December transect, NOT
#    from the collars: the 28 collars span PY 32-148 of 200 m, the median plot
#    cell is 31 m from the nearest one, and only 29% of the plot is within 20 m.
#    They are a clustered sample carrying local structure a plot-scale surface
#    should not reproduce, so they set the LEVEL, not the pattern. Scaling is
#    multiplicative, which preserves relative wetness and makes the surface
#    unbiased in the mean by construction.
#
# Output: outputs/tables/moisture_climatology_monthly.csv
#         outputs/revision/moisture_climatology.txt
# ==============================================================================
suppressPackageStartupMessages({library(dplyr)})

SMAX <- 400            # mm; see note 2 -- insensitive above ~400
YEARS <- 2019:2022

load("outputs/models/TRAINING_DATA.RData")
S <- soil_train_complete %>% mutate(date = as.Date(Date)) %>%
     filter(is.finite(soil_moisture_at_site))

# --- 1. de-confounded plot-wide level per campaign date -----------------------
m0 <- lm(soil_moisture_at_site ~ factor(site_id) + factor(date), data = S)
nd <- expand.grid(site_id = unique(S$site_id), date = unique(S$date))
nd$p <- predict(m0, newdata = data.frame(
  site_id = factor(nd$site_id, levels = levels(factor(S$site_id))),
  date    = factor(nd$date,    levels = levels(factor(S$date)))))
CAMP <- nd %>% group_by(date) %>% summarise(deconf = mean(p), .groups = "drop") %>%
        mutate(date = as.Date(date))
raw <- S %>% group_by(date) %>% summarise(raw = mean(soil_moisture_at_site), .groups = "drop")
cat(sprintf("de-confounding: %d campaign dates | raw SD %.4f -> %.4f (%.0f%% of apparent variation was sampling location)\n",
            nrow(CAMP), sd(raw$raw), sd(CAMP$deconf), 100*(1 - sd(CAMP$deconf)/sd(raw$raw))))

# --- 2. transfer function on the water-balance store --------------------------
WB <- "outputs/tables/tower_daily_waterbalance.csv"
if (!file.exists(WB)) stop("missing ", WB, " -- run code/04_drivers/rev_wb_reference_et.R first")
W <- read.csv(WB, stringsAsFactors = FALSE) %>% arrange(as.Date(date))
W$date <- as.Date(W$date); W$mo <- as.integer(format(W$date, "%m"))
W$yr  <- as.integer(format(W$date, "%Y"))
bucket <- function(p, e, smax) { o <- numeric(length(p)); s <- smax/2
  for (i in seq_along(p)) { s <- min(max(s + ifelse(is.finite(p[i]),p[i],0) -
      ifelse(is.finite(e[i]),e[i],0), 0), smax); o[i] <- s }; o }
W$store <- bucket(W$rain_mm, W$ET0_pm, SMAX)

J <- inner_join(CAMP, W, "date")
fit <- lm(deconf ~ store, data = J)
loo <- sapply(seq_len(nrow(J)), function(i)
  predict(lm(deconf ~ store, data = J[-i, ]), newdata = J[i, ]))
loo_r2 <- 1 - sum((J$deconf - loo)^2)/sum((J$deconf - mean(J$deconf))^2)
cat(sprintf("transfer fit on %d dates: R2 = %.3f | LOO-CV R2 = %.3f | RMSE = %.4f\n",
            nrow(J), summary(fit)$r.squared, loo_r2, sqrt(mean((J$deconf - loo)^2))))
if (loo_r2 < 0.3) warning("transfer function is weak out of sample (LOO R2 < 0.3)")

# --- 3. monthly climatology, and the campaign year for comparison -------------
W$pred <- coef(fit)[1] + coef(fit)[2]*W$store
CLIM <- W %>% filter(yr %in% YEARS) %>% group_by(mo) %>%
  summarise(moisture = mean(pred, na.rm = TRUE),
            between_year_sd = sd(tapply(pred, yr, mean), na.rm = TRUE), .groups = "drop")
CYR <- range(J$date)
CAMPY <- W %>% filter(date >= CYR[1], date <= CYR[2]) %>% group_by(mo) %>%
         summarise(campaign_year = mean(pred), .groups = "drop")
OBS <- CAMP %>% mutate(mo = as.integer(format(date, "%m"))) %>%
       group_by(mo) %>% summarise(observed = mean(deconf), n_dates = dplyr::n(), .groups = "drop")
OUT <- CLIM %>% left_join(CAMPY, "mo") %>% left_join(OBS, "mo") %>%
       mutate(drought_offset = moisture - campaign_year)

sh <- OUT %>% filter(is.finite(observed))
cat(sprintf("back-prediction to the campaign year: shape r = %.3f (peak month %d predicted vs %d observed)\n",
            cor(sh$observed, sh$campaign_year), sh$mo[which.max(sh$campaign_year)], sh$mo[which.max(sh$observed)]))
cat(sprintf("mean drought correction (multi-year minus campaign year): %+.4f m3 m-3; Jun-Sep %+.4f\n",
            mean(OUT$drought_offset), mean(OUT$drought_offset[OUT$mo %in% 6:9])))

cat("\n=== MONTHLY PLOT-WIDE SOIL MOISTURE ===\n")
print(as.data.frame(OUT %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(OUT, "outputs/tables/moisture_climatology_monthly.csv", row.names = FALSE)
dir.create("outputs/revision", showWarnings = FALSE, recursive = TRUE)
con <- file("outputs/revision/moisture_climatology.txt", "w")
cat(sprintf("MOISTURE CLIMATOLOGY (rev_moisture_climatology.R)\nbuilt %s\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = con)
cat(sprintf("index: bucket Smax = %d mm on rainfall and Penman-Monteith ET\n", SMAX), file = con)
cat(sprintf("transfer R2 %.3f | LOO-CV R2 %.3f | shape r %.3f\n",
            summary(fit)$r.squared, loo_r2, cor(sh$observed, sh$campaign_year)), file = con)
cat(sprintf("years averaged: %s\n\n", paste(range(YEARS), collapse = "-")), file = con)
utils::write.table(round(as.data.frame(OUT), 4), con, sep = "\t", row.names = FALSE, quote = FALSE)
close(con)
cat("\nwritten: outputs/tables/moisture_climatology_monthly.csv, outputs/revision/moisture_climatology.txt\n")
