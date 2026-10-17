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
#   Stage-B tables read: tree_tags, trees, weather (built before this script).
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
# TEMPERATURES (per deployment; no soil values: soil conditions are in
# data/tables/soil_readings.csv and are attached by code/lib/soil_fill.R attach_soil()).
#   air_temp_C / air_temp_source
#     monthly   met-tower hourly mean (weather.csv, true time) of the hour containing the
#               logged start (untagged stems: their hand-picked start), else the
#               chamber temperature.
#     2021      the chamber temperature: the field air temperature at the chamber, else
#               the tower hour containing the start (see chamber_temp_source).
#     2023      field sheet air temperature; the stem temperature measured at the deployment
#               where the sheet has none (the station record ends 2023-05-12).
#     felled oak  field sheet air temperature.
#   chamber_temp_C / chamber_temp_source: the Tcham goFlux used, and where it came from
#     (monthly: station air, the hour containing the start in true time; 2021 (stage A,
#     01_prep_auxfile_2021.R): field air at the chamber, else the station air of the hour
#     containing the start (true time); 2023: the sheet's air
#     temperature, else its stem temperature).
#
# CHAMBER GEOMETRY. chamber_area_cm2 and system_volume_L are those goFlux used.
#   geometry_source / chamber_stem_diameter_cm: monthly survey, the rule of
#   code/lib/stem_chamber_geometry.R (the chamber sheet's row for that stem x month x
#   chamber, else computed for the chamber used that day); 2021 and 2023, the rigid
#   chamber's measured area and volume (data/raw/field_data/static_chamber_dims/).
#
# FLAGS, by deployment key:
#   dead_stem      trees.csv: the tree is dead (the one rule in code/tables/04_trees.R)
#                  and the deployment is on or after its dead_from date.
#   qc_pass        not in outputs/data/qc_excluded_measurements.csv (the closure
#                  screen, code/02_flux/02_qc_c0_screen.R); qc_reason from that file.
#   in_rf_training qc_pass, not dead, not the 2021 root-crown deployment, not the felled oak.
# Flux, MDF and detection are checked against flux_FINAL.csv.
# The columns the model reads (CH4 flux; 2023 air temperature and DBH) are parsed as
# readr parses them (read_flux_csv, code/lib/flux_times.R).
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
# untagged stems: start time and window from the hand-picked records. The records carry
# the click step's IDs; the measurement ID is the auxfile's (corrected species,
# 02_flux/rescue/02_untagged_fluxes.R), matched on the start time as there.
aux_id <- read.csv(file.path(FD, "untagged_rescue/untagged_auxfile.csv")) %>% transmute(key = start.time_formatted, id = UniqueID)
man <- readRDS(file.path(WD, "untagged/untagged_manID.rds")) %>%
  mutate(key = format(as.POSIXct(start.time_formatted, tz = "UTC"), "%Y-%m-%d %H:%M:%S")) %>%
  distinct(UniqueID, key, start.time, start.time_corr, end.time_corr) %>%
  mutate(UniqueID = aux_id$id[match(key, aux_id$key)]) %>% select(-key)
stopifnot(!anyNA(man$UniqueID), !anyDuplicated(man$UniqueID))
man <- man %>% transmute(UniqueID, st_untag = wall_clock(format(start.time, tz = "UTC")),
                         win_start_u = wall_clock(format(start.time_corr, tz = "UTC")),
                         win_end_u = wall_clock(format(end.time_corr, tz = "UTC")))
resc <- read.csv(file.path(FD, "untagged_rescue/untagged_monthly_fluxes.csv"), check.names = FALSE) %>%
  transmute(UniqueID, r_model = CH4_model, r_quality = CH4_quality, r_co2 = CO2_best.flux, r_area = Area, r_vtot = Vtot)
aux_u <- read.csv(file.path(FD, "untagged_rescue/untagged_auxfile.csv")) %>%
  transmute(key = start.time_formatted, r_tcham = Tcham, r_geo = geometry_source, r_dstem = Dstem)
# tagged stems: geometry source and stem diameter from the auxfile, by UniqueID
aux_t <- read.csv(file.path(FD, "auxfile_goFlux_with_weather.csv"), check.names = FALSE) %>%
  transmute(UniqueID, t_area = Area, t_vtot = Vtot, t_geo = geometry_source, t_dstem = D_stem_cm)

M <- mon %>% left_join(w_mon, by = "UniqueID") %>% left_join(man, by = "UniqueID") %>% left_join(resc, by = "UniqueID") %>%
  left_join(aux_t, by = "UniqueID") %>%
  mutate(untagged = grepl("^UNTAG_", Plot.Tag),
         st = coalesce(wall_clock(start.time), st_untag),
         key = format(st, "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  left_join(aux_u, by = "key") %>%
  mutate(win_start = coalesce(win_start, win_start_u), win_end = coalesce(win_end, win_end_u))
stopifnot(sum(M$untagged) == 45, all(!is.na(M$st)), all(!is.na(M$win_start)))
stopifnot(all(as.Date(M$st) == as.Date(M$Date)))
# the geometry goFlux used is the auxfile's
stopifnot(all(M$untagged | (abs(num(M$Area) - M$t_area) < 1e-6 & abs(num(M$Vtot) - M$t_vtot) < 1e-9)))
tcham_mon <- ifelse(M$untagged, M$r_tcham, num(M$Tcham))
M$air <- tower_air(M$st)
M$air_src <- ifelse(is.na(M$air), "chamber temperature (no tower record within 24 h)",
                    paste0("met tower, hour of the ", ifelse(M$untagged, "hand-picked", "logged"), " start"))
M$air <- coalesce(M$air, tcham_mon)

mon_out <- M %>% transmute(
  measurement_id = UniqueID, campaign = "monthly_2020_2021",
  tree_id = tree_of(F_MON, "Plot.Tag", Plot.Tag),
  height_cm = 125, height_label = NA_character_, chamber_type = "semirigid",
  chamber_id = ifelse(untagged, sub("^UNTAG_[A-Za-z]+_[A-Z]{4}_\\d{8}_(.*?)(_\\d)?$", "\\1", UniqueID),
                      sub("\\.0$", "", as.character(Chamber_ID))),
  st, win_start, win_end,
  chamber_area_cm2 = coalesce(num(Area), r_area), system_volume_L = coalesce(num(Vtot), r_vtot),
  geometry_source = coalesce(t_geo, r_geo), chamber_stem_diameter_cm = coalesce(t_dstem, r_dstem),
  chamber_temp_C = tcham_mon, chamber_temp_source = "station air, hour containing the start (true time)",
  air_temp_C = air, air_temp_source = air_src,
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model = coalesce(CH4_model, r_model),
  CH4_LM_flux_nmol_m2_s = num(CH4_LM.flux), CH4_LM_r2 = num(CH4_LM.r2), CH4_LM_pval = num(CH4_LM.p.val),
  CH4_HM_flux_nmol_m2_s = num(CH4_HM.flux), CH4_quality = coalesce(CH4_quality.check, r_quality),
  CO2_flux_umol_m2_s = coalesce(num(CO2_best.flux), r_co2), CO2_model,
  bark_missing = NA_character_, wounding = NA_character_, moss_lichen_cover = NA_character_,
  untagged_rescue = untagged, root_crown = FALSE, notes = as.character(Notes))
mon_out$chamber_id[mon_out$chamber_id %in% c("na", "NA", "")] <- NA

# ==============================================================================
# 2. 2021 MULTI-HEIGHT (rigid chambers; 0.5, 1.25, 2 m; Kalmia 0.75 m and root crown)
# ==============================================================================
F_H21 <- file.path(FD, "tree_flux_2021_multiheight.csv")
h21 <- read_flux_csv(F_H21, "CH4_best.flux") %>% filter(!is.na(CH4_best.flux)) %>%
  left_join(read_windows(file.path(WD, "lgr_manual_identification_2021_multiheight_final.csv")), by = "UniqueID")
h21$tree <- tree_of(F_H21, "tree_id", h21$tree_id)
hl <- trimws(h21$measurement_height)
h21_out <- h21 %>% transmute(
  measurement_id = UniqueID, campaign = "2021_multiheight",
  tree_id = tree,
  height_cm = num(sub(" cm$", "", hl)), height_label = ifelse(hl == "Root Crown", "root crown", NA_character_),
  chamber_type = "rigid", chamber_id = as.character(chamber_id),
  st = wall_clock(start.time), win_start, win_end,
  chamber_area_cm2 = Area, system_volume_L = Vtot,
  geometry_source = "rigid chamber: measured area and volume (static_chamber_dims)", chamber_stem_diameter_cm = NA_real_,
  chamber_temp_C = Tcham,
  chamber_temp_source = recode(temp_source, field_air = "field air temperature at the chamber",
                               tower_hour = "met tower, hour containing the start (no field air temperature)"),
  air_temp_C = Tcham, air_temp_source = ifelse(temp_source == "tower_hour", "met tower, hour of the logged start (no field reading)",
                                               paste("chamber temperature:", chamber_temp_source)),
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_LM_flux_nmol_m2_s = CH4_LM.flux, CH4_LM_r2 = CH4_LM.r2,
  CH4_LM_pval = CH4_LM.p.val, CH4_HM_flux_nmol_m2_s = CH4_HM.flux, CH4_quality = CH4_quality.check,
  CO2_flux_umol_m2_s = CO2_best.flux, CO2_model,
  bark_missing = NA_character_, wounding = NA_character_, moss_lichen_cover = NA_character_,
  untagged_rescue = FALSE,
  root_crown = hl == "Root Crown", notes = NA_character_)
stopifnot(!anyNA(h21_out$chamber_temp_source), all(h21$temp_source %in% c("field_air", "tower_hour")))

# ==============================================================================
# 3. 2023 CROSS-SPECIES (rigid chambers, 1.25 m) -- analyzer clock correction
# ==============================================================================
F_23 <- file.path(FD, "tree_flux_2023_cross_species.csv")
y <- read_flux_csv(F_23, c("CH4_best.flux", "air_temp_C", "DBH (cm)")) %>% filter(!is.na(CH4_best.flux)) %>%
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
# Condition scores are 1-3 (or "dead"). (The two Bark Missing cells that were
# time-formatted in the 2023 sheet, trees 132 and 6352 on 2023-07-12, are blank in the
# raw sheet with a note since 2026-10-17: data/raw/CORRECTIONS.md.)
cond <- function(x) { x <- trimws(as.character(x)); x[x %in% c("", "NA")] <- NA; ifelse(tolower(x) == "dead", "dead", x) }
stopifnot(all(unlist(lapply(y[c("Bark Missing", "Wounding/Holes", "Moss/Lichen Cover", "Visible fungus")], cond)) %in%
                c(NA, "1", "2", "3", "dead")))
y_field <- num(y$air_temp_C)
# air temperature at the deployment: the sheet's air reading, else the stem temperature
# measured there (the 2023 campaign lies after the tower record ends, 2023-05-12); the
# same rule as the chamber temperature (stage A, 02_prep_auxfile_2023.R)
y_air <- coalesce(y_field, num(y$stem_temp_C))
stopifnot(!anyNA(y_air), isTRUE(all.equal(y_air, num(y$Tcham))))
y_out <- y %>% transmute(
  measurement_id = UniqueID, campaign = "2023_cross_species",
  tree_id = tree_of(F_23, "Tree Tag", `Tree Tag`, `Species Code`, `DBH (cm)`),
  height_cm = 125, height_label = NA_character_, chamber_type = "rigid", chamber_id = as.character(`Chamber ID`),
  dbh_recorded_cm = num(`DBH (cm)`),
  st = wall_clock(start.time) - off_s, win_start = win_start - off_s, win_end = win_end - off_s,
  chamber_area_cm2 = Area, system_volume_L = Vtot,
  geometry_source = "rigid chamber: measured area and volume (static_chamber_dims)", chamber_stem_diameter_cm = NA_real_,
  chamber_temp_C = Tcham,
  chamber_temp_source = ifelse(!is.na(y_field), "field sheet air temperature", "field sheet stem temperature (no air temperature)"),
  air_temp_C = y_air,
  air_temp_source = ifelse(!is.na(y_field), "field sheet", "field sheet stem temperature (no air temperature)"),
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_LM_flux_nmol_m2_s = CH4_LM.flux, CH4_LM_r2 = CH4_LM.r2,
  CH4_LM_pval = CH4_LM.p.val, CH4_HM_flux_nmol_m2_s = CH4_HM.flux, CH4_quality = CH4_quality.check,
  CO2_flux_umol_m2_s = CO2_best.flux, CO2_model,
  bark_missing = cond(`Bark Missing`), wounding = cond(`Wounding/Holes`), moss_lichen_cover = cond(`Moss/Lichen Cover`),
  untagged_rescue = FALSE, root_crown = FALSE,
  notes = as.character(Notes))

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
  chamber_area_cm2 = NA_real_, system_volume_L = NA_real_, geometry_source = NA_character_, chamber_stem_diameter_cm = NA_real_,
  chamber_temp_C = NA_real_, chamber_temp_source = NA_character_,
  air_temp_C = Air_Temp_C, air_temp_source = ifelse(is.na(Air_Temp_C), NA_character_, "felled-oak field sheet"),
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_LM_flux_nmol_m2_s = CH4_LM.flux, CH4_LM_r2 = CH4_LM.r2,
  CH4_LM_pval = CH4_LM.p.val, CH4_HM_flux_nmol_m2_s = CH4_HM.flux, CH4_quality = CH4_quality.check,
  CO2_flux_umol_m2_s = CO2_best.flux, CO2_model,
  bark_missing = NA_character_, wounding = NA_character_, moss_lichen_cover = NA_character_,
  untagged_rescue = FALSE, root_crown = FALSE,
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
TR <- read.csv(file.path(OUT, "trees.csv"), colClasses = "character", na.strings = "")
ti <- match(S$tree_id, TR$tree_id)
stopifnot(!anyNA(ti))
QX <- read.csv("outputs/data/qc_excluded_measurements.csv", stringsAsFactors = FALSE)
qx <- match(S$measurement_id, QX$UniqueID)
S <- S %>% mutate(
  closure_s = DET$t_sec[d],
  CH4_MDF_nmol_m2_s = coalesce(DET$MDF[d], oak_mdf), analyzer_sigma_ppb = DET$sigma[d],
  detected = DET$detected[d], detection_class = DET$class[d],
  date = as.Date(format(st, "%Y-%m-%d", tz = "UTC")),
  dead_stem = TR$dead[ti] %in% "TRUE" & date >= as.Date(TR$dead_from[ti]),
  dead_stem = coalesce(dead_stem, FALSE),
  qc_pass = is.na(qx), qc_reason = QX$reason[qx],
  in_rf_training = qc_pass & !dead_stem & !root_crown & campaign != "2022_felled_oak",
  exclusion_reason = case_when(campaign == "2022_felled_oak" ~ "felled black oak: analyzed separately",
                               !qc_pass ~ "chamber not at ambient at closure",
                               dead_stem ~ "dead stem: the inventory the model is applied to is live stems only",
                               root_crown ~ "root-crown deployment: no stem height"),
  start_time = iso_local(st), window_start = iso_local(win_start), window_end = iso_local(win_end))
# monthly stems: the tree's landscape position (trees.csv) is the plot letter of the sheet
pl <- toupper(trimws(M$Plot.Letter))[match(S$measurement_id[S$campaign == "monthly_2020_2021"], M$UniqueID)]
pos_mm <- sum(TR$landscape_position[ti][S$campaign == "monthly_2020_2021"] != pl, na.rm = TRUE)

# ---- previous table (for the report only) ----------------------------------------
F_OLD <- file.path(OUT, "flux_stem.csv")
OLD <- if (file.exists(F_OLD)) read.csv(F_OLD, stringsAsFactors = FALSE) else NULL

out <- S %>% select(measurement_id, campaign, tree_id, height_cm, height_label, chamber_type, chamber_id,
                    date, start_time, window_start, window_end, closure_s, chamber_area_cm2, system_volume_L,
                    geometry_source, chamber_stem_diameter_cm, chamber_temp_C, chamber_temp_source,
                    air_temp_C, air_temp_source,
                    CH4_flux_nmol_m2_s, CH4_model, CH4_LM_flux_nmol_m2_s, CH4_LM_r2, CH4_LM_pval,
                    CH4_HM_flux_nmol_m2_s, CH4_quality, CO2_flux_umol_m2_s, CO2_model,
                    CH4_MDF_nmol_m2_s, analyzer_sigma_ppb, detected, detection_class,
                    dbh_recorded_cm, bark_missing, wounding, moss_lichen_cover, qc_pass, qc_reason, dead_stem,
                    untagged_rescue, in_rf_training, exclusion_reason, notes)
write_csv(out, file.path(OUT, "flux_stem.csv"), na = "")

# ---- report -------------------------------------------------------------------
cat(sprintf("flux_stem.csv: %d deployments (%s)\n", nrow(out),
            paste(names(table(out$campaign)), table(out$campaign), collapse = "; ")))
cat(sprintf("  in_rf_training %d; dead_stem %d; qc_pass FALSE %d; untagged_rescue %d\n",
            sum(out$in_rf_training), sum(out$dead_stem), sum(!out$qc_pass), sum(out$untagged_rescue)))
cat(sprintf("  start_time missing: %d; window missing: %d (felled oak: window not archived)\n",
            sum(is.na(out$start_time)), sum(is.na(out$window_start))))
cat(sprintf("  monthly stems whose trees.csv landscape_position differs from the sheet's plot letter: %d\n", pos_mm))
cat("  geometry_source:\n"); print(as.data.frame(out %>% count(campaign, geometry_source)), row.names = FALSE)
cat("  chamber_temp_source:\n"); print(as.data.frame(out %>% count(campaign, chamber_temp_source)), row.names = FALSE)
cat("  air_temp_source:\n"); print(as.data.frame(out %>% count(campaign, air_temp_source)), row.names = FALSE)
cat("  dead stems (trees.csv rule):\n")
print(as.data.frame(out %>% filter(tree_id %in% TR$tree_id[TR$dead == "TRUE" | !is.na(TR$dead_evidence)]) %>%
  group_by(tree_id, campaign) %>% summarise(n = n(), dead = sum(dead_stem), first = min(date), last = max(date), .groups = "drop")), row.names = FALSE)
cat("  2023 analyzer clock offsets (min, analyzer minus real):\n")
for (i in seq_len(nrow(DAY))) cat(sprintf("    %s  %5.2f  %s (pairs %d; deployments with flux %d)\n", DAY$date[i],
                                         DAY$offset_min[i], DAY$source[i], DAY$n_pairs[i], DAY$n_flux[i]))
if (nrow(dropped)) cat(sprintf("    dropped pairs (entry errors): %s\n",
                               paste(sprintf("%s %.1f min", dropped$date, dropped$offset_min), collapse = "; ")))
cat(sprintf("  felled oak: analyzer minus real = %s min; real times used\n", paste(oak_offset, collapse = ", ")))
if (!is.null(OLD)) {
  oi <- match(out$measurement_id, OLD$measurement_id)
  r <- out$CH4_flux_nmol_m2_s / OLD$CH4_flux_nmol_m2_s[oi]
  cat(sprintf("  vs the previous flux_stem.csv: %d deployments not in it; CH4 flux changed (|ratio-1| > 1e-9) in %d\n",
              sum(is.na(oi)), sum(abs(r - 1) > 1e-9, na.rm = TRUE)))
  for (cp in unique(out$campaign)) { k <- out$campaign == cp & is.finite(r) & abs(r - 1) > 1e-9
    if (any(k)) cat(sprintf("    %-20s n %3d  ratio min %.4f  median %.4f  max %.4f\n", cp, sum(k), min(r[k]), median(r[k]), max(r[k]))) }
  if ("dead_stem" %in% names(OLD)) cat(sprintf("  dead_stem changed: %d (now TRUE: %s; now FALSE: %s)\n",
    sum(out$dead_stem != OLD$dead_stem[oi], na.rm = TRUE),
    paste(out$measurement_id[which(out$dead_stem & !OLD$dead_stem[oi])], collapse = ", "),
    paste(out$measurement_id[which(!out$dead_stem & OLD$dead_stem[oi])], collapse = ", ")))
  cat(sprintf("  in_rf_training: previous %d, now %d\n", sum(OLD$in_rf_training), sum(out$in_rf_training)))
}
