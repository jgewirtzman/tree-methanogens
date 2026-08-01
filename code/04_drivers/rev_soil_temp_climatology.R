#!/usr/bin/env Rscript
# ==============================================================================
# rev_soil_temp_climatology.R  --  multi-year monthly soil temperature
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
# WHY NOT THE TOWER. ymf_clean_sorted.csv has Tsoil_Avg, but it is not measuring
# soil: r = 0.993 with air temperature, damping ratio 1.03, zero lag at daily
# resolution, seasonal amplitude 24.3 C against air's 23.6 C, and 5% of days below
# -5 C. A probe in the ground damps and lags air by weeks. The column is unusable.
#
# WHAT THIS DOES INSTEAD. Soil temperature at depth is, to good approximation, a
# damped and lagged version of air temperature, so a trailing mean of air
# temperature is the natural predictor. The relation is fitted on the campaign
# measurements -- which ARE genuine soil temperature -- and then applied to the
# full multi-year tower record to produce the climatology.
#
# TWO THINGS THAT KEEP THIS HONEST.
#
#   (1) Fitted on DATE MEANS, not individual observations. The 1,325 campaign
#       measurements come from only 50 distinct dates, so observation-level
#       fitting is pseudoreplicated: it would weight a date by how many collars
#       were read that day and report an R2 inflated by within-date replication.
#
#   (2) Window selected by LEAVE-ONE-DATE-OUT cross-validation, not by in-sample
#       R2. Choosing among 10 candidate windows on the same points that report the
#       fit is exactly the selection bias this project has been correcting
#       elsewhere.
#
# Output: outputs/tables/soil_temp_climatology_monthly.csv
# ==============================================================================
suppressMessages({library(dplyr)})
outdir <- "outputs/tables"
load("outputs/models/TRAINING_DATA.RData")

# ---- daily air temperature from the tower ------------------------------------
w <- read.csv("data/raw/weather/ymf_clean_sorted.csv", stringsAsFactors = FALSE)
w$date <- as.Date(w$TIMESTAMP)
w$ta   <- suppressWarnings(as.numeric(w$Tair_Avg))
AIR <- w %>% filter(is.finite(ta), !is.na(date)) %>%
  group_by(date) %>% summarise(ta = mean(ta), .groups = "drop") %>% arrange(date)
cat(sprintf("tower air temperature: %d days, %s to %s\n",
            nrow(AIR), min(AIR$date), max(AIR$date)))

# ---- campaign soil temperature, collapsed to one value per date --------------
OBS <- bind_rows(
    tree_train_complete %>% transmute(date = as.Date(Date), st = soil_temp_C_mean),
    soil_train_complete %>% transmute(date = as.Date(Date), st = soil_temp_C_mean)) %>%
  filter(is.finite(st), !is.na(date)) %>%
  group_by(date) %>% summarise(st = mean(st), n_obs = dplyr::n(), .groups = "drop")
cat(sprintf("campaign soil temperature: %d dates (%d measurements), %s to %s\n",
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

# ---- climatology over the whole tower record ---------------------------------
AIR$st_hat <- as.numeric(predict(FIT, AIR))
CLIM <- AIR %>% filter(is.finite(st_hat)) %>%
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
cat(sprintf("\nwritten: %s/soil_temp_climatology_monthly.csv\n", outdir))
