#!/usr/bin/env Rscript
# ==============================================================================
# 00_soil_temp_climatology.R  --  multi-year monthly soil temperature
# ------------------------------------------------------------------------------
# Puts soil temperature on the same footing as moisture and air temperature: a
# multi-year monthly climatology rather than the mean of whenever we happened to
# visit.
#
# WHAT WAS THERE BEFORE. 01_load_and_prep_data.R built soil_temp_C_mean by
# averaging campaign measurements per month (TREE_JULY + TREE_YEAR + SOIL_YEAR).
# That series has 27 to 624 observations a month, JANUARY AND MARCH ARE ABSENT
# ENTIRELY and get linearly interpolated, and each month's value is the
# temperature on the few days somebody was in the field.
#
# WHY NOT THE TOWER. The station logs Tsoil_Avg, but it is not measuring
# soil: r = 0.993 with air temperature, damping ratio 1.03, zero lag at daily
# resolution, seasonal amplitude 24.3 C against air's 23.6 C, and 5% of days below
# -5 C. A probe in the ground damps and lags air by weeks. The column is unusable.
#
# WHAT THIS DOES INSTEAD. Soil temperature at depth is, to good approximation, a
# damped and lagged version of air temperature, so a trailing mean of air
# temperature is the natural predictor. The relation is fitted on the campaign
# measurements -- which ARE genuine soil temperature -- and then applied to the
# full multi-year air record to produce the climatology.
#
# AIR TEMPERATURE: the tower (data/tables/weather.csv). A day is the local day of the
# true time (an hour belongs to the day of its midpoint) with at least 20 hourly
# records, its value the 24-h mean -- the definition of 05_model/01_load_and_prep_data.R.
# A trailing window is a run of consecutive CALENDAR days: the record is laid on the
# full calendar, and a window that spans a day with no value has no value; gaps are
# never closed up. The climatology is the mean over the climatology years
# (CLIM_YEARS, the years of the moisture climatology, 02_moisture_climatology.R).
#
# THE OBSERVATIONS are data/tables/soil_readings.csv: every soil temperature read at a
# collar (monthly campaign 2020-21), a tree hole (2021 survey) or a tree (2023 field
# sheet). (Until 2026-10-16 they were the model's training rows, whose soil
# temperature had itself been gap-filled from a monthly mean: circular.) The script
# reads tables only, so it runs in stage B before the model: 05_model/01 uses its
# monthly values as the soil-temperature driver, for training rows without a reading
# and for the stand predictions alike.
#
# TWO THINGS THAT KEEP THIS HONEST.
#
#   (1) Fitted on DATE MEANS, not individual observations. The readings come from
#       a few dozen distinct dates, so observation-level
#       fitting is pseudoreplicated: it would weight a date by how many collars
#       were read that day and report an R2 inflated by within-date replication.
#
#   (2) Window selected by LEAVE-ONE-DATE-OUT cross-validation, not by in-sample
#       R2. Choosing among 10 candidate windows on the same points that report the
#       fit is exactly the selection bias this project has been correcting
#       elsewhere.
#
# Output: outputs/tables/soil_temp_climatology_monthly.csv
#         outputs/tables/soil_temp_daily.csv  (daily air, trailing air and fitted soil
#           temperature, CLIM_YEARS; read only by the daily-driver sensitivity, 06_upscale/11)
# ==============================================================================
suppressMessages({library(dplyr)})
outdir <- "outputs/tables"

CLIM_YEARS <- 2019:2022
MIN_HOURS <- 20

# ---- daily air temperature from the tower, on the full calendar ---------------
w <- read.csv("data/tables/weather.csv", stringsAsFactors = FALSE)
w$t  <- as.POSIXct(substr(sub("T", " ", w$timestamp), 1, 19), tz = "UTC")   # local wall-clock digits
w$ta <- suppressWarnings(as.numeric(w$air_temp_avg_C))
dd <- w %>% filter(is.finite(ta), !is.na(t)) %>% mutate(date = as.Date(t - 1800)) %>%
  group_by(date) %>% summarise(ta = mean(ta), n = dplyr::n(), .groups = "drop") %>% filter(n >= MIN_HOURS)
AIR <- data.frame(date = seq(min(dd$date), max(dd$date), by = "day"))
AIR$ta <- dd$ta[match(AIR$date, dd$date)]
cat(sprintf("tower air temperature: %d complete days (>= %d h), %s to %s; %d calendar days without a value\n",
            sum(is.finite(AIR$ta)), MIN_HOURS, min(AIR$date), max(AIR$date), sum(!is.finite(AIR$ta))))

# ---- campaign soil temperature, collapsed to one value per date --------------
OBS <- read.csv("data/tables/soil_readings.csv", stringsAsFactors = FALSE) %>%
  transmute(date = as.Date(date), st = suppressWarnings(as.numeric(soil_temp_C))) %>%
  filter(is.finite(st), !is.na(date)) %>%
  group_by(date) %>% summarise(st = mean(st), n_obs = dplyr::n(), .groups = "drop")
cat(sprintf("soil_readings soil temperature: %d dates (%d readings), %s to %s\n",
            nrow(OBS), sum(OBS$n_obs), min(OBS$date), max(OBS$date)))
cat("months represented by DATES (not measurements):\n")
print(table(as.integer(format(OBS$date, "%m"))))

# ---- choose the trailing window by leave-one-date-out CV ---------------------
trail <- function(N) as.numeric(stats::filter(AIR$ta, rep(1/N, N), sides = 1))
WINDOWS <- c(1, 3, 5, 7, 10, 14, 21, 30, 45, 60)
cat("\ntrailing-window selection, leave-one-date-out CV on date means:\n")
score <- sapply(WINDOWS, function(N) {
  A <- AIR; A$sm <- trail(N)
  J <- OBS %>% left_join(A %>% select(date, sm), "date") %>% filter(is.finite(sm))
  if (nrow(J) < 20) return(c(NA, NA, NA))
  pred <- vapply(seq_len(nrow(J)), function(i)
    as.numeric(predict(lm(st ~ sm, J[-i, ]), J[i, ])), numeric(1))
  c(cv_R2 = 1 - sum((J$st - pred)^2)/sum((J$st - mean(J$st))^2),
    cv_RMSE = sqrt(mean((J$st - pred)^2)),
    in_R2 = summary(lm(st ~ sm, J))$r.squared)
})
TAB <- data.frame(window_days = WINDOWS, cv_R2 = score[1, ],
                  cv_RMSE_C = score[2, ], in_sample_R2 = score[3, ])
print(as.data.frame(TAB %>% mutate(across(-window_days, ~round(.x, 3)))), row.names = FALSE)

BEST <- WINDOWS[which.max(TAB$cv_R2)]
AIR$sm <- trail(BEST)
J <- OBS %>% left_join(AIR %>% select(date, sm), "date") %>% filter(is.finite(sm))
FIT <- lm(st ~ sm, J)
cat(sprintf("\nselected window %d days | CV R2 %.3f, RMSE %.2f C on %d dates\n",
            BEST, max(TAB$cv_R2), TAB$cv_RMSE_C[which.max(TAB$cv_R2)], nrow(J)))
cat(sprintf("  soil_temp = %.3f * air_%dd %+.3f\n", coef(FIT)[2], BEST, coef(FIT)[1]))
cat(sprintf("  slope < 1 is the damping a buried sensor must show (tower Tsoil gives 1.03)\n"))

# ---- climatology over the climatology years ------------------------------------
AIR$st_hat <- as.numeric(predict(FIT, AIR))
CLIM <- AIR %>% filter(is.finite(st_hat), as.integer(format(date, "%Y")) %in% CLIM_YEARS) %>%
  mutate(mo = as.integer(format(date, "%m")),
         yr = as.integer(format(date, "%Y"))) %>%
  group_by(mo) %>%
  summarise(soil_temp_C = mean(st_hat), between_year_sd = sd(tapply(st_hat, yr, mean)),
            n_days = dplyr::n(), .groups = "drop")

CAMP <- OBS %>% mutate(mo = as.integer(format(date, "%m"))) %>%
  group_by(mo) %>% summarise(campaign_C = mean(st), n_dates = dplyr::n(), .groups = "drop")
CMP <- CLIM %>% left_join(CAMP, "mo")
cat("\nmonthly soil temperature (C):\n")
print(as.data.frame(CMP %>% mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)

# The disagreement is structured, and that is the argument for the climatology:
# where campaign sampling is dense the two agree, where it is thin they do not.
d <- CMP %>% filter(is.finite(campaign_C)) %>% mutate(diff = abs(soil_temp_C - campaign_C))
cat(sprintf("\nmean |difference| where both exist: %.2f C over %d months\n",
            mean(d$diff), nrow(d)))
cat(sprintf("correlation between |difference| and number of campaign dates: %+.2f\n",
            suppressWarnings(cor(d$diff, d$n_dates))))
cat(sprintf("months the campaign series cannot supply at all: %s\n",
            paste(setdiff(1:12, CAMP$mo), collapse = ", ")))

stopifnot(nrow(CLIM) == 12, all(is.finite(CLIM$soil_temp_C)))
write.csv(CLIM, file.path(outdir, "soil_temp_climatology_monthly.csv"), row.names = FALSE)
write.csv(AIR %>% filter(as.integer(format(date, "%Y")) %in% CLIM_YEARS) %>%
            transmute(date, air_C = ta, air_trailing_C = sm, window_days = BEST, soil_temp_C = st_hat),
          file.path(outdir, "soil_temp_daily.csv"), row.names = FALSE)
cat(sprintf("\nwritten: %s/soil_temp_climatology_monthly.csv\n", outdir))
