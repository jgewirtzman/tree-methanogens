#!/usr/bin/env Rscript
# ==============================================================================
# 11_flux_soil.R -- data/tables/flux_soil.csv: one row per soil chamber deployment
# ------------------------------------------------------------------------------
# 288 deployments with a CH4 flux, monthly survey 2020-21 (static chambers on PVC
# collars). measurement_id is the goFlux UniqueID; collar_id = <Plot letter>_<Plot Tag>
# (collars.csv).
#
# Sources: data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv (the one
# soil fit, 02_flux/semirigid/04_goflux_soils.R; checked equal to flux_FINAL); windows in
# data/raw/flux_windows/lgr_manual_identification_results_soil.csv; detection from
# data/processed/flux/flux_FINAL.csv; closure screen from
# data/processed/flux/qc_excluded_measurements.csv.
#
# ROW ORDER: the field-sheet order, with the 30 December 2020 deployments last. That is
# the order the soil forest has always been trained on (until 2026-10-18 a separate
# December fit appended them), and a forest's bootstrap depends on row order.
#
# TIMES: the analyzer clock was kept on local time in 2020-21 (as for the stems);
# written as ISO 8601 local time with the America/New_York offset of the date.
# Five 19 May 2021 deployments have two saved windows; the one whose length matches
# the closure length of the fit (flux_FINAL t_sec) is reported.
#
# TEMPERATURES. No soil values: soil conditions are in data/tables/soil_readings.csv
# and are attached by code/lib/soil_fill.R attach_soil(unit = "collar").
#   air_temp_C / air_temp_source  met-tower hourly mean (data/tables/weather.csv, true
#                          time) of the hour containing the logged start, else the
#                          chamber temperature.
#   chamber_temp_C         the Tcham goFlux used (the same station hour, stage A).
#
# FLAGS by key: qc_pass (closure screen, data/processed/flux/qc_excluded_measurements.csv);
# in_rf_training = qc_pass and the collar is in the stand (collars.csv in_stand).
# CH4 flux is parsed as readr parses it (read_flux_csv, code/lib/flux_times.R), the
# values the model has always used.
#
# Writes data/tables/flux_soil.csv
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/flux_times.R")
OUT <- "data/tables"; WD <- "data/raw/flux_windows"
num <- function(x) suppressWarnings(as.numeric(x))

F_SOIL <- "data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv"
s <- read_flux_csv(F_SOIL, "CH4_best.flux") %>% filter(!is.na(CH4_best.flux)) %>%
  arrange(startsWith(UniqueID, "202012"))   # December 2020 last (ROW ORDER, above)
DET <- read.csv("data/processed/flux/flux_FINAL.csv") %>% filter(type == "soil")
QX <- read.csv("data/processed/flux/qc_excluded_measurements.csv", stringsAsFactors = FALSE)
COL <- read.csv(file.path(OUT, "collars.csv"), stringsAsFactors = FALSE)

# ---- windows ----------------------------------------------------------------------
W <- read_windows(file.path(WD, "lgr_manual_identification_results_soil.csv"), keep_all = TRUE) %>%
  left_join(DET %>% select(UniqueID, t_sec), by = "UniqueID") %>%
  mutate(fit = is.na(t_sec) | abs(as.numeric(win_end - win_start, units = "secs") + 1 - t_sec) < 0.5) %>%
  group_by(UniqueID) %>% arrange(!fit, .by_group = TRUE) %>%
  mutate(n_win = n_distinct(paste(win_start, win_end)), ok = any(fit)) %>% slice(1) %>% ungroup()
stopifnot(all(W$ok[W$UniqueID %in% s$UniqueID]))

# ---- environment as the model builds it -----------------------------------------
S <- s %>% transmute(
  UniqueID, date = as.Date(Date), letter = toupper(trimws(`Plot letter`)), tag = as.character(`Plot Tag`),
  collar_id = paste0(letter, "_", tag), st = wall_clock(start.time), Area, Vtot, Tcham, Notes,
  CH4_best.flux, CH4_model, CH4_LM.flux, CH4_LM.r2, CH4_LM.p.val, CH4_HM.flux, CH4_quality.check,
  CO2_best.flux, CO2_model) %>%
  left_join(W %>% select(UniqueID, win_start, win_end), by = "UniqueID")
S <- S %>% mutate(air = tower_air(st),
                  air_temp_source = ifelse(is.na(air), "chamber temperature (no tower record within 24 h)", "met tower, hour of the logged start"),
                  air = coalesce(air, Tcham))
cat(sprintf("flux_soil: %d deployments\n", nrow(S)))
stopifnot(!anyDuplicated(S$UniqueID), all(S$collar_id %in% COL$collar_id),
          all(!is.na(S$st)), all(!is.na(S$win_start)))

d <- match(S$UniqueID, DET$UniqueID)
qx <- match(S$UniqueID, QX$UniqueID)
stopifnot(all(is.na(d) == !is.na(qx)),
          isTRUE(all.equal(S$CH4_best.flux[!is.na(d)], DET$best.flux[d[!is.na(d)]], tolerance = 1e-12)))
in_stand <- COL$in_stand[match(S$collar_id, COL$collar_id)]
out <- S %>% transmute(
  measurement_id = UniqueID, campaign = "monthly_2020_2021", collar_id, chamber_type = "static",
  date, start_time = iso_local(st), window_start = iso_local(win_start), window_end = iso_local(win_end),
  closure_s = DET$t_sec[d], chamber_area_cm2 = Area, system_volume_L = Vtot, chamber_temp_C = Tcham,
  air_temp_C = air, air_temp_source,
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_LM_flux_nmol_m2_s = CH4_LM.flux, CH4_LM_r2 = CH4_LM.r2,
  CH4_LM_pval = CH4_LM.p.val, CH4_HM_flux_nmol_m2_s = CH4_HM.flux, CH4_quality = CH4_quality.check,
  CO2_flux_umol_m2_s = CO2_best.flux, CO2_model,
  CH4_MDF_nmol_m2_s = DET$MDF[d], analyzer_sigma_ppb = DET$sigma[d], detected = DET$detected[d],
  detection_class = DET$class[d],
  qc_pass = is.na(qx), qc_reason = QX$reason[qx],
  in_rf_training = qc_pass & in_stand,
  exclusion_reason = case_when(!qc_pass ~ "chamber not at ambient at closure",
                               !in_stand ~ "collar outside the modelled stand"),
  notes = as.character(Notes))

F_OLD <- file.path(OUT, "flux_soil.csv")
OLD <- if (file.exists(F_OLD)) read.csv(F_OLD, stringsAsFactors = FALSE) else NULL
write_csv(out, file.path(OUT, "flux_soil.csv"), na = "")
cat(sprintf("flux_soil.csv: %d deployments on %d collars; in_rf_training %d; qc fail %d; outside stand %d\n",
            nrow(out), n_distinct(out$collar_id), sum(out$in_rf_training), sum(!out$qc_pass),
            sum(out$exclusion_reason %in% "collar outside the modelled stand")))
cat(sprintf("  windows: %d deployments had two saved windows (the one matching the fit's closure length kept)\n",
            sum(W$n_win[W$UniqueID %in% out$measurement_id] > 1)))
cat(sprintf("  UTC offsets: %s\n", paste(names(table(substr(out$start_time, 20, 25))), table(substr(out$start_time, 20, 25)), collapse = "; ")))
cat(sprintf("  air_temp_source: %s\n", paste(names(table(out$air_temp_source)), table(out$air_temp_source), collapse = "; ")))
if (!is.null(OLD)) {
  oi <- match(out$measurement_id, OLD$measurement_id)
  cat(sprintf("  vs the previous flux_soil.csv: %d not in it; CH4 flux changed in %d; in_rf_training %d -> %d\n", sum(is.na(oi)),
              sum(abs(out$CH4_flux_nmol_m2_s - OLD$CH4_flux_nmol_m2_s[oi]) > 1e-12, na.rm = TRUE),
              sum(OLD$in_rf_training), sum(out$in_rf_training)))
}
