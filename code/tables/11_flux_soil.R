#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 11_flux_soil.R -- data/tables/flux_soil.csv: one row per soil chamber deployment
# ------------------------------------------------------------------------------
# 288 deployments with a CH4 flux, monthly survey 2020-21 (static chambers on PVC
# collars). measurement_id is the goFlux UniqueID; collar_id = <Plot letter>_<Plot Tag>
# (collars.csv).
#
# Sources: data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv (fits
# with the corrected collar geometry, 03_merge/01_fix_soil_flux.R); windows in
# data/raw/flux_windows/lgr_manual_identification_results{,_december}_soil.csv;
# detection from outputs/data/flux_FINAL.csv; closure screen from
# outputs/data/qc_excluded_measurements.csv.
#
# TIMES: the analyzer clock was kept on local time in 2020-21 (as for the stems);
# written as ISO 8601 local time with the America/New_York offset of the date.
# Five 19 May 2021 deployments have two saved windows; the one whose length matches
# the closure length of the fit (flux_FINAL t_sec) is reported.
#
# ENVIRONMENT (the values the model trains on, read by 05_model/01_load_and_prep_data.R; env_source per row):
#   soil temperature, VWC  the collar's own reading that date
#                          (data/processed/environmental/soil_env_by_collar.csv); else
#                          the mean for that plot tag (e.g. 1-3, all landscape
#                          positions) that date in soilmoisture_total.csv.
#   air temperature        met-tower Tair_Avg nearest 00:00 of the date (the model's
#                          lookup), else the chamber temperature (tower at deployment).
#
# FLAGS by key: qc_pass (closure screen); in_rf_training = qc_pass and the collar is
# in the stand (collars.csv in_stand). Checked against the previous stage-B run's
# export (outputs/data/flux_measurements_soil.csv, matched by collar + flux to 8
# decimals; skipped when absent) and data/compiled/flux_soil.csv.
# CH4 flux is parsed as readr parses it (read_flux_csv, code/lib/flux_times.R), the
# values the model has always used.
#
# Writes data/tables/flux_soil.csv
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/flux_times.R")
OUT <- "data/tables"; WD <- "data/raw/flux_windows"
num <- function(x) suppressWarnings(as.numeric(x))

F_SOIL <- "data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv"
s <- read_flux_csv(F_SOIL, "CH4_best.flux") %>% filter(!is.na(CH4_best.flux))
DET <- read.csv("outputs/data/flux_FINAL.csv") %>% filter(type == "soil")
QX <- read.csv("outputs/data/qc_excluded_measurements.csv", stringsAsFactors = FALSE)
COL <- read.csv(file.path(OUT, "collars.csv"), stringsAsFactors = FALSE)

# ---- windows: December 2020 file first, then the main file ----------------------
W <- bind_rows(read_windows(file.path(WD, "lgr_manual_identification_results_december_soil.csv")) %>% mutate(src = 1),
               read_windows(file.path(WD, "lgr_manual_identification_results_soil.csv"), keep_all = TRUE) %>% mutate(src = 2)) %>%
  left_join(DET %>% select(UniqueID, t_sec), by = "UniqueID") %>%
  mutate(fit = is.na(t_sec) | abs(as.numeric(win_end - win_start, units = "secs") + 1 - t_sec) < 0.5) %>%
  group_by(UniqueID) %>% arrange(src, !fit, .by_group = TRUE) %>%
  mutate(n_win = n_distinct(paste(win_start, win_end)), ok = any(fit)) %>% slice(1) %>% ungroup()
stopifnot(all(W$ok[W$UniqueID %in% s$UniqueID]))

# ---- environment as the model builds it -----------------------------------------
SE <- read.csv("data/processed/environmental/soil_env_by_collar.csv", stringsAsFactors = FALSE) %>%
  transmute(date = as.Date(Date), letter = plot_letter, tag = plot_tag, c_temp = soil_temp_C,
            c_vwc = num(soil_moisture_pct) / 100)
stopifnot(!anyDuplicated(SE[c("date", "letter", "tag")]))
SM <- read.csv("data/raw/field_data/ipad_data/Cleaned data/soilmoisture_total.csv", fileEncoding = "UTF-8-BOM",
               stringsAsFactors = FALSE) %>%
  transmute(date = lubridate::mdy(Date, quiet = TRUE), tag = as.character(Plot), st = num(Soil.temp), sv = num(VWC) / 100) %>%
  filter(!is.na(date), !is.na(sv)) %>% group_by(date, tag) %>%
  summarise(p_temp = mean(st, na.rm = TRUE), p_vwc = mean(sv, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(c(p_temp, p_vwc), ~ ifelse(is.nan(.x), NA_real_, .x)))

S <- s %>% transmute(
  UniqueID, date = as.Date(Date), letter = toupper(trimws(`Plot letter`)), tag = as.character(`Plot Tag`),
  collar_id = paste0(letter, "_", tag), st = wall_clock(start.time), Area, Vtot, Tcham, Notes,
  CH4_best.flux, CH4_model, CH4_LM.flux, CH4_LM.r2, CH4_LM.p.val, CH4_HM.flux, CH4_quality.check,
  CO2_best.flux, CO2_model) %>%
  left_join(W %>% select(UniqueID, win_start, win_end), by = "UniqueID") %>%
  left_join(SE, by = c("date", "letter", "tag")) %>% left_join(SM, by = c("date", "tag")) %>%
  mutate(air = coalesce(tower_air(as.POSIXct(as.character(date), tz = "UTC")), Tcham),
         soil_temp_C = coalesce(c_temp, p_temp), soil_vwc = coalesce(c_vwc, p_vwc),
         env_source = paste0(case_when(!is.na(c_temp) | !is.na(c_vwc) ~ "soil: collar reading that date",
                                       !is.na(p_temp) | !is.na(p_vwc) ~ "soil: plot-tag mean that date (soilmoisture_total.csv)",
                                       TRUE ~ "soil: no reading"),
                             "; air: met tower nearest 00:00 of the date"))
stopifnot(nrow(S) == 288, !anyDuplicated(S$UniqueID), all(S$collar_id %in% COL$collar_id),
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
  air_temp_C = air, soil_temp_C, soil_vwc, env_source,
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

# ---- checks against the model's export and the archive --------------------------
# The table is built before the model (stage B reads it), so the export compared is
# that of the previous stage-B run, when there is one.
key8 <- function(x) round(as.numeric(x), 8)
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-9)
tr <- which(out$in_rf_training)
F_FM <- "outputs/data/flux_measurements_soil.csv"
dup_in <- dup_fm <- mm_train <- NA; env_mm <- rep(NA, 3)
if (file.exists(F_FM)) {
  FM <- read.csv(F_FM)
  mk <- paste(out$collar_id, key8(out$CH4_flux_nmol_m2_s)); fk <- paste(FM$site_id, key8(FM$soil_flux_nmol_m2_s))
  fi <- match(mk, fk)
  dup_in <- sum(duplicated(mk)); dup_fm <- sum(duplicated(fk))
  mm_train <- sum(out$in_rf_training != !is.na(fi))
  env_mm <- c(sum(!same(out$air_temp_C[tr], FM$air_temp_C[fi[tr]])), sum(!same(out$soil_temp_C[tr], FM$soil_temp_C[fi[tr]])),
              sum(!same(out$soil_vwc[tr], FM$soil_moisture_abs[fi[tr]])))
  stopifnot(dup_in == 0, dup_fm == 0, mm_train == 0, sum(out$in_rf_training) == nrow(FM))
}
C <- read.csv("data/compiled/flux_soil.csv", stringsAsFactors = FALSE)
ci <- match(out$measurement_id, C$unique_id)
eq <- function(a, b) all(same(as.numeric(a), as.numeric(b)))
stopifnot(!anyNA(ci),
          eq(out$CH4_flux_nmol_m2_s, C$CH4_flux_nmol_m2_s[ci]), eq(out$CH4_MDF_nmol_m2_s, C$CH4_MDF_nmol_m2_s[ci]),
          eq(out$analyzer_sigma_ppb, C$analyzer_sigma_ppb[ci]), eq(out$closure_s, C$closure_s[ci]),
          identical(out$detected, as.logical(C$detected[ci])), identical(out$detection_class, C$detection_class[ci]),
          identical(out$qc_pass, as.logical(C$qc_pass[ci])), identical(out$in_rf_training, as.logical(C$in_rf_training[ci])),
          identical(out$collar_id, C$site_id[ci]), eq(out$chamber_area_cm2, C$chamber_area_cm2[ci]),
          eq(out$system_volume_L, C$chamber_vol_L[ci]), eq(out$CO2_flux_umol_m2_s, C$CO2_flux_umol_m2_s[ci]))
reason_mm <- sum(!is.na(out$exclusion_reason) & out$exclusion_reason != C$exclusion_reason[ci])
stopifnot(reason_mm == 0)

write_csv(out, file.path(OUT, "flux_soil.csv"), na = "")
cat(sprintf("flux_soil.csv: %d deployments on %d collars; in_rf_training %d; qc fail %d; outside stand %d\n",
            nrow(out), n_distinct(out$collar_id), sum(out$in_rf_training), sum(!out$qc_pass),
            sum(out$exclusion_reason %in% "collar outside the modelled stand")))
cat(sprintf("  windows: %d deployments had two saved windows (the one matching the fit's closure length kept)\n",
            sum(W$n_win[W$UniqueID %in% out$measurement_id] > 1)))
cat(sprintf("  UTC offsets: %s\n", paste(names(table(substr(out$start_time, 20, 25))), table(substr(out$start_time, 20, 25)), collapse = "; ")))
cat(sprintf("  flux-value join (collar + flux to 8 dp): %d duplicate keys; %d training mismatches vs the key-based flags\n",
            dup_in, mm_train))
cat(sprintf("  environment vs the model's training rows (%d): mismatches air %d, soil temp %d, VWC %d\n",
            length(tr), env_mm[1], env_mm[2], env_mm[3]))
cat(sprintf("  env source: %s\n", paste(names(table(sub(";.*", "", out$env_source))), table(sub(";.*", "", out$env_source)), collapse = "; ")))
