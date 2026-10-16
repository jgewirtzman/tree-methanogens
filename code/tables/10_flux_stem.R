#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 10_flux_stem.R -- data/tables/flux_stem.csv: one row per stem chamber deployment
# ------------------------------------------------------------------------------
# 1,200 deployments with a CH4 flux: monthly survey 2020-21 (414, incl. 45 untagged
# stems rescued from hand-picked windows), 2021 multi-height (439), 2023
# cross-species (338), 2022 felled oak (9). measurement_id is the goFlux UniqueID.
#
# Sources (stage A): data/processed/flux/{semirigid_tree_final_complete_dataset_with_untagged,
#   tree_flux_2021_multiheight,tree_flux_2023_cross_species}.csv;
#   data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv; closure windows in
#   data/raw/flux_windows/ (untagged: untagged/untagged_manID.rds); detection
#   (MDF, sigma, closure length, class) from outputs/data/flux_FINAL.csv.
# tree_id: tree_tags.csv (source file + column + spelling; species and DBH where
#   one spelling names two trees).
# Untagged monthly stems: the merged file has no start time for them. Start time and
#   window come from their hand-picked window record (untagged_manID.rds, by UniqueID);
#   chamber temperature from untagged_auxfile.csv matched on start time, because 23 of
#   the 45 window IDs carry a species code or _1/_2 suffix from before the species
#   correction and no longer equal the aux IDs. The tree is the merged file's Plot.Tag.
#
# TIMES. goFlux start.time and the window times carry a "Z" but are the analyzer's
# wall clock. Tables hold real local time, ISO 8601 with the America/New_York offset
# of that date (-05:00 EST, -04:00 EDT).
#   monthly 2020-21  analyzer clock kept on local time across both DST changes
#                    (field sampling times vs windows: 1-2 min all year): as recorded.
#   2021 heights     field times were read from the analyzer: as recorded.
#   2023             analyzer clock fast. Offset = Flux Start (System) - Flux Start
#                    (Real) where both were written (25 deployments): per-day median;
#                    days without a pair take the offset linearly interpolated in
#                    time between measured days (held constant before the first).
#                    Offsets outside 0-60 min are entry errors (a 12-h slip) and are
#                    dropped. Per-day offsets: outputs/audit/flux_2023_clock_offsets.csv.
#   felled oak 2022  real_start / real_end from the field sheet (clock 22 min fast);
#                    the fitted window is not archived.
# start_time = goFlux start.time (chamber closure as logged); window_start/_end =
# the hand-picked fitting window (start.time_corr / end.time_corr).
#
# ENVIRONMENT (the values the model trains on, read by 05_model/01_load_and_prep_data.R; env_source per row):
#   monthly   air: met-tower Tair_Avg nearest the logged start (within 24 h; noon of
#             the date where the merged file has no start time), else chamber temp.
#             soil: per-collar readings (data/processed/environmental/soil_env_by_collar.csv)
#             averaged by date and landscape position; nearest date within 3 days.
#   2021      air at the chamber height and tree-mean soil temperature / VWC from
#             data/processed/integrated/merged_tree_dataset_final.csv.
#   2023      field sheet (air temp, soil temp, mean of the VWC probes); tower at
#             00:00 of the date where air temp is missing.
#   felled oak  field sheet air temperature.
#
# FLAGS, by deployment key (no flux-value matching):
#   dead_stem      monthly `dead` column / _dead label; 2021 name containing "dead";
#                  2023 bark-missing or wounding score "dead" (code/lib/dead_stems.R rules).
#   in_rf_training not dead, not the 2021 root-crown deployment, not the felled oak.
# Both are checked against outputs/data/flux_measurements_tree.csv (the previous
# stage-B run's export, matched by campaign + flux to 8 decimals; skipped when absent)
# and data/compiled/flux_stem.csv; flux, MDF and detection against flux_FINAL.csv and
# the archive.
# The columns the model reads (CH4 flux; 2023 air temperature, soil temperature, VWC
# and DBH) are parsed as readr parses them (read_flux_csv, code/lib/flux_times.R).
# dbh_recorded_cm: DBH written on the 2023 field sheet for the deployment (2023 only).
#
# Writes data/tables/flux_stem.csv, outputs/audit/flux_2023_clock_offsets.csv
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/monthly_cols.R"); source("code/lib/flux_times.R")
OUT <- "data/tables"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FD <- "data/processed/flux"; WD <- "data/raw/flux_windows"
num <- function(x) suppressWarnings(as.numeric(x))

# ---- tree_id from tree_tags ----------------------------------------------------
TT <- read.csv(file.path(OUT, "tree_tags.csv"), colClasses = "character", na.strings = "")
tree_of <- function(file, column, tag, species = NA, dbh = NA) {
  t <- TT[TT$source_file == file & TT$source_column == column, ]
  tag <- trimws(as.character(tag)); out <- rep(NA_character_, length(tag))
  for (i in seq_along(tag)) {
    c <- t[t$tag_recorded %in% tag[i], ]
    if (nrow(c) > 1) c <- c[c$species_recorded %in% species[i] & abs(num(c$dbh_recorded_cm) - num(dbh[i])) < 0.05, ]
    if (nrow(c) == 1) out[i] <- c$tree_id
  }
  if (anyNA(out)) stop(file, ": no tree_tags entry for ", paste(unique(tag[is.na(out)]), collapse = ", "), call. = FALSE)
  out
}

# ==============================================================================
# 1. MONTHLY SURVEY 2020-21 (semi-rigid chambers, 1.25 m)
# ==============================================================================
F_MON <- file.path(FD, "semirigid_tree_final_complete_dataset_with_untagged.csv")
mon <- monthly_plain(read_flux_csv(F_MON, "CH4_best.flux")) %>% filter(!is.na(CH4_best.flux))
names(mon) <- gsub(" ", ".", names(mon))
w_mon <- read_windows(file.path(WD, "lgr_manual_identification_semirigid_tree.csv"))
# untagged stems: start time and window from the hand-picked records, by UniqueID
man <- readRDS(file.path(WD, "untagged/untagged_manID.rds")) %>%
  distinct(UniqueID, start.time, start.time_corr, end.time_corr)
stopifnot(!anyDuplicated(man$UniqueID))
man <- man %>% transmute(UniqueID, st_untag = wall_clock(format(start.time, tz = "UTC")),
                         win_start_u = wall_clock(format(start.time_corr, tz = "UTC")),
                         win_end_u = wall_clock(format(end.time_corr, tz = "UTC")))
resc <- read.csv(file.path(FD, "untagged_rescue/untagged_monthly_fluxes.csv"), check.names = FALSE) %>%
  transmute(UniqueID, r_model = CH4_model, r_quality = CH4_quality, r_co2 = CO2_best.flux, r_area = Area, r_vtot = Vtot)
aux_u <- read.csv(file.path(FD, "untagged_rescue/untagged_auxfile.csv")) %>%
  transmute(key = start.time_formatted, r_tcham = Tcham)

M <- mon %>% left_join(w_mon, by = "UniqueID") %>% left_join(man, by = "UniqueID") %>% left_join(resc, by = "UniqueID") %>%
  mutate(untagged = grepl("^UNTAG_", Plot.Tag),
         st = coalesce(wall_clock(start.time), st_untag),
         key = format(st, "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  left_join(aux_u, by = "key") %>%
  mutate(win_start = coalesce(win_start, win_start_u), win_end = coalesce(win_end, win_end_u))
stopifnot(sum(M$untagged) == 45, all(!is.na(M$st)), all(!is.na(M$win_start)))

# environment as the model builds it
M$t_model <- coalesce(wall_clock(M$start.time), as.POSIXct(paste(as.Date(M$Date), "12:00:00"), tz = "UTC"))
SE <- read.csv("data/processed/environmental/soil_env_by_collar.csv", stringsAsFactors = FALSE) %>%
  mutate(Date = as.Date(Date), vwc = num(soil_moisture_pct) / 100) %>% filter(!is.na(Date), !is.na(plot_letter))
SE_pos <- SE %>% group_by(Date, plot_letter) %>%
  summarise(st = mean(soil_temp_C, na.rm = TRUE), sv = mean(vwc, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(c(st, sv), ~ ifelse(is.nan(.x), NA_real_, .x)))
pos_env <- function(d, pos) {                      # nearest date (earliest on ties) within 3 days
  c <- SE_pos[SE_pos$plot_letter == pos, ]
  if (!nrow(c)) return(c(NA, NA))
  dd <- abs(as.numeric(d - c$Date)); j <- order(dd)[1]
  if (dd[j] <= 3) c(c$st[j], c$sv[j]) else c(NA, NA)
}
pe <- t(mapply(pos_env, as.Date(M$t_model), toupper(trimws(M$Plot.Letter))))
tcham_mon <- ifelse(M$untagged, M$r_tcham, num(M$Tcham))
M$air <- coalesce(tower_air(M$t_model), num(M$Tcham))
M$soilT <- pe[, 1]; M$soilV <- pe[, 2]

mon_out <- M %>% transmute(
  measurement_id = UniqueID, campaign = "monthly_2020_2021",
  tree_id = tree_of(F_MON, "Plot.Tag", Plot.Tag),
  height_cm = 125, height_label = NA_character_, chamber_type = "semirigid",
  chamber_id = ifelse(untagged, sub("^UNTAG_[A-Za-z]+_[A-Z]{4}_\\d{8}_(.*?)(_\\d)?$", "\\1", UniqueID),
                      sub("\\.0$", "", as.character(Chamber_ID))),
  st, win_start, win_end,
  chamber_area_cm2 = coalesce(num(Area), r_area), system_volume_L = coalesce(num(Vtot), r_vtot),
  chamber_temp_C = tcham_mon, air_temp_C = air, soil_temp_C = soilT, soil_vwc = soilV,
  env_source = paste0("air: met tower nearest the ", ifelse(is.na(start.time), "date at 12:00", "logged start"),
                      " (chamber temp if none within 24 h); soil: collar readings, mean by landscape position, nearest date within 3 days"),
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model = coalesce(CH4_model, r_model),
  CH4_LM_flux_nmol_m2_s = num(CH4_LM.flux), CH4_LM_r2 = num(CH4_LM.r2), CH4_LM_pval = num(CH4_LM.p.val),
  CH4_HM_flux_nmol_m2_s = num(CH4_HM.flux), CH4_quality = coalesce(CH4_quality.check, r_quality),
  CO2_flux_umol_m2_s = coalesce(num(CO2_best.flux), r_co2), CO2_model,
  bark_missing = NA_character_, wounding = NA_character_, moss_lichen_cover = NA_character_,
  dead_stem = grepl("_dead$", Plot.Tag) | toupper(as.character(dead)) %in% "TRUE",
  untagged_rescue = untagged, root_crown = FALSE, notes = as.character(Notes))
mon_out$chamber_id[mon_out$chamber_id %in% c("na", "NA", "")] <- NA

# ==============================================================================
# 2. 2021 MULTI-HEIGHT (rigid chambers; 0.5, 1.25, 2 m; Kalmia 0.75 m and root crown)
# ==============================================================================
F_H21 <- file.path(FD, "tree_flux_2021_multiheight.csv")
h21 <- read_flux_csv(F_H21, "CH4_best.flux") %>% filter(!is.na(CH4_best.flux)) %>%
  left_join(read_windows(file.path(WD, "lgr_manual_identification_2021_multiheight_final.csv")), by = "UniqueID")
F_MT <- "data/processed/integrated/merged_tree_dataset_final.csv"
MT <- read.csv(F_MT, check.names = FALSE)
MT <- MT[rowSums(!is.na(MT[grep("^CH4_best.flux_", names(MT))])) > 0, ]   # trees with a flux
MT$k <- tree_of(F_MT, "tree_id", MT$tree_id)
h21$tree <- tree_of(F_H21, "tree_id", h21$tree_id)
h21$dead <- grepl("dead", h21$tree_id, ignore.case = TRUE)            # field name, e.g. "Bowser dead red oak"
stopifnot(!anyDuplicated(MT$k))
hl <- trimws(h21$measurement_height)
h21$hcol <- ifelse(hl == "Root Crown", "RootCrown", sub(" cm$", "", hl))
mi <- match(h21$tree, MT$k)
stopifnot(!anyNA(mi))
air21 <- vapply(seq_len(nrow(h21)), function(i) {
  if (h21$hcol[i] == "75" && h21$species[i] == "KALA")       # Kalmia: model reads 75 cm, else root crown
    return(coalesce(MT$Temp_Air_75cm[mi[i]], MT$Temp_Air_RootCrowncm[mi[i]]))
  MT[[paste0("Temp_Air_", h21$hcol[i], "cm")]][mi[i]]
}, numeric(1))
h21_out <- h21 %>% transmute(
  measurement_id = UniqueID, campaign = "2021_multiheight",
  tree_id = tree,
  height_cm = num(sub(" cm$", "", hl)), height_label = ifelse(hl == "Root Crown", "root crown", NA_character_),
  chamber_type = "rigid", chamber_id = as.character(chamber_id),
  st = wall_clock(start.time), win_start, win_end,
  chamber_area_cm2 = Area, system_volume_L = Vtot, chamber_temp_C = Tcham,
  air_temp_C = air21, soil_temp_C = MT$SoilTemp_mean[mi], soil_vwc = MT$VWC_mean[mi] / 100,
  env_source = "air: field air temperature at the chamber height; soil: tree mean of the 2021 survey soil readings (merged_tree_dataset_final.csv)",
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_LM_flux_nmol_m2_s = CH4_LM.flux, CH4_LM_r2 = CH4_LM.r2,
  CH4_LM_pval = CH4_LM.p.val, CH4_HM_flux_nmol_m2_s = CH4_HM.flux, CH4_quality = CH4_quality.check,
  CO2_flux_umol_m2_s = CO2_best.flux, CO2_model,
  bark_missing = NA_character_, wounding = NA_character_, moss_lichen_cover = NA_character_,
  dead_stem = dead, untagged_rescue = FALSE,
  root_crown = hl == "Root Crown", notes = NA_character_)

# ==============================================================================
# 3. 2023 CROSS-SPECIES (rigid chambers, 1.25 m) -- analyzer clock correction
# ==============================================================================
F_23 <- file.path(FD, "tree_flux_2023_cross_species.csv")
y <- read_flux_csv(F_23, c("CH4_best.flux", "air_temp_C", "vwc_mean", "Soil Temp (°C)", "DBH (cm)")) %>% filter(!is.na(CH4_best.flux)) %>%
  left_join(read_windows(file.path(WD, "lgr_manual_identification_2023_cross_species.csv")), by = "UniqueID")
hms_s <- function(x) {                       # "9:42:00", "10;31", "9:51" -> seconds after midnight
  p <- strsplit(gsub(";", ":", trimws(x)), ":")
  vapply(p, function(v) { v <- num(v); if (length(v) < 2 || anyNA(v[1:2])) NA_real_ else
    v[1] * 3600 + v[2] * 60 + ifelse(length(v) > 2 && !is.na(v[3]), v[3], 0) }, numeric(1))
}
y_all <- read.csv(F_23, check.names = FALSE, stringsAsFactors = FALSE)
pair <- y_all %>% transmute(date = as.Date(date_clean), sys = hms_s(`Flux Start (System) Time`),
                            real = hms_s(`Flux Start (Real) Time`)) %>%
  filter(!is.na(sys), !is.na(real)) %>% mutate(offset_min = (sys - real) / 60)
dropped <- pair %>% filter(offset_min <= 0 | offset_min >= 60)
OFF <- pair %>% filter(offset_min > 0, offset_min < 60) %>% group_by(date) %>%
  summarise(n_pairs = n(), offset_min = median(offset_min), .groups = "drop")
DAY <- y_all %>% count(date = as.Date(date_clean), name = "n_rows") %>%
  left_join(y %>% count(date = as.Date(date_clean), name = "n_flux"), by = "date") %>%
  mutate(n_flux = coalesce(n_flux, 0L)) %>% left_join(OFF, by = "date") %>%
  mutate(n_pairs = coalesce(n_pairs, 0L), source = ifelse(n_pairs > 0, "median of that day's pairs", "interpolated"),
         offset_min = ifelse(n_pairs > 0, offset_min,
                             approx(as.numeric(OFF$date), OFF$offset_min, xout = as.numeric(date), rule = 2)$y))
write_csv(DAY, out_path("flux_2023_clock_offsets.csv", "audit"))
off_s <- 60 * DAY$offset_min[match(as.Date(y$date_clean), DAY$date)]
bark <- y[["Bark Missing"]]; bark[bark %in% "0:00:00"] <- "0"      # spreadsheet read 0 as a time
cond <- function(x) { x <- trimws(as.character(x)); x[x %in% c("", "NA")] <- NA; ifelse(tolower(x) == "dead", "dead", x) }
y_air <- coalesce(num(y$air_temp_C), tower_air(as.POSIXct(y$date_clean, tz = "UTC")))
y_out <- y %>% transmute(
  measurement_id = UniqueID, campaign = "2023_cross_species",
  tree_id = tree_of(F_23, "Tree Tag", `Tree Tag`, `Species Code`, `DBH (cm)`),
  height_cm = 125, height_label = NA_character_, chamber_type = "rigid", chamber_id = as.character(`Chamber ID`),
  dbh_recorded_cm = num(`DBH (cm)`),
  st = wall_clock(start.time) - off_s, win_start = win_start - off_s, win_end = win_end - off_s,
  chamber_area_cm2 = Area, system_volume_L = Vtot, chamber_temp_C = Tcham,
  air_temp_C = y_air, soil_temp_C = num(`Soil Temp (°C)`), soil_vwc = num(vwc_mean) / 100,
  env_source = ifelse(is.na(num(air_temp_C)), "air: met tower nearest 00:00 of the date (no field reading); soil: 2023 field sheet",
                      "2023 field sheet (air temperature; soil temperature; mean of the VWC probes)"),
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_LM_flux_nmol_m2_s = CH4_LM.flux, CH4_LM_r2 = CH4_LM.r2,
  CH4_LM_pval = CH4_LM.p.val, CH4_HM_flux_nmol_m2_s = CH4_HM.flux, CH4_quality = CH4_quality.check,
  CO2_flux_umol_m2_s = CO2_best.flux, CO2_model,
  bark_missing = cond(bark), wounding = cond(`Wounding/Holes`), moss_lichen_cover = cond(`Moss/Lichen Cover`),
  untagged_rescue = FALSE, root_crown = FALSE, notes = as.character(Notes)) %>%
  mutate(dead_stem = bark_missing %in% "dead" | wounding %in% "dead")

# ==============================================================================
# 4. 2022 FELLED BLACK OAK (rigid chamber along the felled stem, 0.5-10 m)
# ==============================================================================
F_OAK <- "data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv"
oak <- read.csv(F_OAK, check.names = FALSE, stringsAsFactors = FALSE) %>% filter(!is.na(CH4_best.flux))
h <- trimws(oak$Height_m)
oak_tree <- TT$tree_id[TT$campaign == "2022_felled_oak"][1]
stopifnot(length(oak_tree) == 1, !is.na(oak_tree))
oak_out <- oak %>% transmute(
  measurement_id = UniqueID, campaign = "2022_felled_oak", tree_id = oak_tree,
  height_cm = 100 * num(sub(" .*", "", h)), height_label = ifelse(tolower(h) == "seam", "seam", NA_character_),
  chamber_type = "rigid", chamber_id = Chamber,
  st = as.POSIXct(real_start, tz = "UTC"), win_start = as.POSIXct(NA), win_end = as.POSIXct(NA),
  chamber_area_cm2 = NA_real_, system_volume_L = NA_real_, chamber_temp_C = NA_real_,
  air_temp_C = Air_Temp_C, soil_temp_C = NA_real_, soil_vwc = NA_real_,
  env_source = ifelse(is.na(Air_Temp_C), NA_character_, "felled-oak field sheet (air temperature)"),
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_LM_flux_nmol_m2_s = CH4_LM.flux, CH4_LM_r2 = CH4_LM.r2,
  CH4_LM_pval = CH4_LM.p.val, CH4_HM_flux_nmol_m2_s = CH4_HM.flux, CH4_quality = CH4_quality.check,
  CO2_flux_umol_m2_s = CO2_best.flux, CO2_model,
  bark_missing = NA_character_, wounding = NA_character_, moss_lichen_cover = NA_character_,
  dead_stem = FALSE, untagged_rescue = FALSE, root_crown = FALSE,
  notes = case_when(grepl("restart", h) ~ "repeat of the 4 m deployment after the analyzer battery died",
                    TRUE ~ as.character(Notes)),
  oak_mdf = CH4_MDF, oak_sys = as.POSIXct(sys_start, tz = "UTC"))
oak_offset <- unique(as.numeric(difftime(oak_out$oak_sys, oak_out$st, units = "mins")))

# ==============================================================================
# 5. ASSEMBLE; DETECTION; FLAGS
# ==============================================================================
S <- bind_rows(mon_out, h21_out, y_out, oak_out)
stopifnot(nrow(S) == 1200, !anyDuplicated(S$measurement_id), !anyNA(S$tree_id))
DET <- read.csv("outputs/data/flux_FINAL.csv") %>% filter(type == "stem")
d <- match(S$measurement_id, DET$UniqueID)
stopifnot(all(!is.na(d) | S$campaign == "2022_felled_oak"),
          isTRUE(all.equal(S$CH4_flux_nmol_m2_s[!is.na(d)], DET$best.flux[d[!is.na(d)]], tolerance = 1e-12)))
S <- S %>% mutate(
  closure_s = DET$t_sec[d],
  CH4_MDF_nmol_m2_s = coalesce(DET$MDF[d], oak_mdf), analyzer_sigma_ppb = DET$sigma[d],
  detected = DET$detected[d], detection_class = DET$class[d],
  qc_pass = TRUE, qc_reason = NA_character_,           # no stem deployment failed the closure screen
  in_rf_training = !dead_stem & !root_crown & campaign != "2022_felled_oak",
  exclusion_reason = case_when(campaign == "2022_felled_oak" ~ "felled black oak: analyzed separately",
                               dead_stem ~ "dead stem: the inventory the model is applied to is live stems only",
                               root_crown ~ "root-crown deployment: no stem height"),
  date = as.Date(format(st, "%Y-%m-%d", tz = "UTC")),
  start_time = iso_local(st), window_start = iso_local(win_start), window_end = iso_local(win_end))

# ---- checks against the model's export and the archive --------------------------
# The table is built before the model (stage B reads it), so the export compared is
# that of the previous stage-B run, when there is one.
key8 <- function(x) round(as.numeric(x), 8)
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-9)
surv <- S$campaign != "2022_felled_oak"
tr <- which(surv & S$in_rf_training)
F_FM <- "outputs/data/flux_measurements_tree.csv"
mm_train <- mm_dead <- mm_reason <- integer(0); dup_in <- dup_fm <- NA; env_mm <- rep(NA, 3)
if (file.exists(F_FM)) {
  FM <- read.csv(F_FM)
  mk <- paste(S$campaign, key8(S$CH4_flux_nmol_m2_s)); fk <- paste(FM$campaign, key8(FM$stem_flux_nmol_m2_s))
  dup_in <- sum(duplicated(mk[surv])); dup_fm <- sum(duplicated(fk))
  fi <- match(mk, fk)
  stopifnot(dup_in == 0, dup_fm == 0, all(!is.na(fi[surv])), length(unique(fi[surv])) == nrow(FM))
  mm_train <- which(surv & S$in_rf_training != FM$in_rf_training[fi])
  mm_dead <- which(surv & S$dead_stem != FM$dead_stem[fi])
  mm_reason <- which(surv & !S$in_rf_training & S$exclusion_reason != FM$exclusion_reason[fi])
  # environment: equal to the model's training rows?
  env_mm <- c(sum(!same(S$air_temp_C[tr], FM$air_temp_C[fi[tr]])),
              sum(!same(S$soil_temp_C[tr], FM$soil_temp_C[fi[tr]])),
              sum(!same(S$soil_vwc[tr], FM$soil_moisture_abs[fi[tr]])))
}
# dead_stems.R matches dead flux values across campaigns: would any live deployment collide?
dead_vals <- key8(S$CH4_flux_nmol_m2_s[S$dead_stem])
x_camp <- sum(!S$dead_stem & key8(S$CH4_flux_nmol_m2_s) %in% dead_vals)
C <- read.csv("data/compiled/flux_stem.csv", stringsAsFactors = FALSE)
ci <- match(S$measurement_id, C$unique_id)
stopifnot(!anyNA(ci))
eq <- function(a, b) all(same(as.numeric(a), as.numeric(b)))
stopifnot(eq(S$CH4_flux_nmol_m2_s, C$CH4_flux_nmol_m2_s[ci]), eq(S$CH4_MDF_nmol_m2_s, C$CH4_MDF_nmol_m2_s[ci]),
          eq(S$analyzer_sigma_ppb, C$analyzer_sigma_ppb[ci]), eq(S$closure_s, C$closure_s[ci]),
          identical(S$detected, as.logical(C$detected[ci])), identical(S$detection_class, C$detection_class[ci]),
          eq(S$CH4_LM_r2, C$CH4_LM_r2[ci]), eq(S$CH4_LM_pval, C$CH4_LM_pval[ci]),
          eq(S$CO2_flux_umol_m2_s[!S$untagged_rescue], C$CO2_flux_umol_m2_s[ci][!S$untagged_rescue]),
          eq(S$chamber_area_cm2[surv & !S$untagged_rescue], C$chamber_area_cm2[ci][surv & !S$untagged_rescue]))
stopifnot(length(mm_train) == 0, length(mm_dead) == 0, length(mm_reason) == 0,
          sum(S$in_rf_training) == 1131, sum(S$dead_stem) == 59,
          identical(S$in_rf_training, as.logical(C$in_rf_training[ci])), identical(S$dead_stem, as.logical(C$dead_stem[ci])))

out <- S %>% select(measurement_id, campaign, tree_id, height_cm, height_label, chamber_type, chamber_id,
                    date, start_time, window_start, window_end, closure_s, chamber_area_cm2, system_volume_L,
                    chamber_temp_C, air_temp_C, soil_temp_C, soil_vwc, env_source,
                    CH4_flux_nmol_m2_s, CH4_model, CH4_LM_flux_nmol_m2_s, CH4_LM_r2, CH4_LM_pval,
                    CH4_HM_flux_nmol_m2_s, CH4_quality, CO2_flux_umol_m2_s, CO2_model,
                    CH4_MDF_nmol_m2_s, analyzer_sigma_ppb, detected, detection_class,
                    dbh_recorded_cm, bark_missing, wounding, moss_lichen_cover, qc_pass, qc_reason, dead_stem,
                    untagged_rescue, in_rf_training, exclusion_reason, notes)
write_csv(out, file.path(OUT, "flux_stem.csv"), na = "")

# ---- report -------------------------------------------------------------------
cat(sprintf("flux_stem.csv: %d deployments (%s)\n", nrow(out),
            paste(names(table(out$campaign)), table(out$campaign), collapse = "; ")))
cat(sprintf("  in_rf_training %d; dead_stem %d; untagged_rescue %d\n",
            sum(out$in_rf_training), sum(out$dead_stem), sum(out$untagged_rescue)))
cat(sprintf("  start_time missing: %d; window missing: %d (felled oak: window not archived)\n",
            sum(is.na(out$start_time)), sum(is.na(out$window_start))))
cat("  2023 analyzer clock offsets (min, analyzer minus real):\n")
for (i in seq_len(nrow(DAY))) cat(sprintf("    %s  %5.2f  %s (pairs %d; deployments with flux %d)\n", DAY$date[i],
                                         DAY$offset_min[i], DAY$source[i], DAY$n_pairs[i], DAY$n_flux[i]))
if (nrow(dropped)) cat(sprintf("    dropped pairs (entry errors): %s\n",
                               paste(sprintf("%s %.1f min", dropped$date, dropped$offset_min), collapse = "; ")))
cat(sprintf("  felled oak: analyzer minus real = %s min; real times used\n", paste(oak_offset, collapse = ", ")))
if (file.exists(F_FM)) cat(sprintf("  vs the previous model export, flux-value join (campaign + flux to 8 dp): %d duplicate keys; %d training / %d dead / %d reason mismatches vs the key-based flags\n",
            dup_in, length(mm_train), length(mm_dead), length(mm_reason))) else cat("  no previous model export: flags not compared\n")
cat(sprintf("  dead_stems.R cross-campaign value match: %d live deployments share a dead deployment's flux value\n", x_camp))
cat(sprintf("  environment vs the model's training rows (%d): mismatches air %d, soil temp %d, VWC %d\n",
            length(tr), env_mm[1], env_mm[2], env_mm[3]))
