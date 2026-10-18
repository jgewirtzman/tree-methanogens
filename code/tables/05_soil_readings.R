#!/usr/bin/env Rscript
# ==============================================================================
# 05_soil_readings.R -- soil_readings.csv (data/tables/)
# ------------------------------------------------------------------------------
# One row per soil reading, at a tree (tree_id) or at a soil collar (collar_id);
# exactly one of the two is set.
#
# TREES (2021 survey, July-August 2021): one row per hole, four holes per tree
# (N, E, S, W of the stem). Source data/raw/field_data/static_chamber_field/
# Soils_Data_filled.csv -- the field sheet (Soils Data_cleaned.csv) with the tree
# name filled down and the holes shared by two or three trees ("RO8 & RO7",
# "YB11', 12, 13", ...) repeated under each tree; shared_with names the other trees.
#   soil_vwc_pct: mean of the probe readings at the hole (VWC1-3); a hole whose
#     first reading is "saturated" (no other reading) is 60 % -- the value
#     soil_tree_level_means.csv uses; vwc_saturated marks these holes.
#     RO2 hole N's third reading was recorded "3,9" (decimal comma); corrected to 3.9
#     in the raw sheets (data/raw/CORRECTIONS.md, 2026-10-15).
#   mineral_depth_cm: NA where recorded as ">45" or "water" (kept in notes), as in
#     soil_tree_level_means.csv.
#   date: holes whose date was not written (40 holes at 9 trees: PB2, WP7, RO1, RM3,
#     RO103, WP103, RM304, H301, WPBERENICE) take the date of the tree's 2021 flux
#     visit, and date_source says so. The survey read a tree's holes on its flux
#     visit: of the trees with dated holes, the build reports how many were read the
#     same day as the flux and how many within 3 days.
# The tree means in soil_tree_level_means.csv (the field team's file) are the means
# over a tree's holes.
#   Two sets of holes are recorded under a combined name and were never repeated under
#   each tree in the filled sheet: "AB2/3" (note "one sample for AB 2/3") is AB2 and
#   AB3; "WA14/1/XV" is WA14, WA1 and WAX1 -- "XV" read as "X1": WAX1 is the third
#   Lab-plot white ash measured in 2021 and has no soil reading of its own (none of
#   the three has). Their holes are repeated under each tree, with shared_with naming
#   the others, as for the other shared holes (SHARED_NAMES below).
#
# COLLARS (monthly campaign, 2020-21): one row per collar per date per source
# record: the per-date iPad sheets "Soil temp moisture YYYYMMDD*.xlsx" and the
# 2020-21 master workbook data/raw/environmental/Soil_temp_moisture_2020.xlsx.
# Where a collar-date appears in more than one record (both sources, or twice in
# a sheet) each record is a row.
# soilmoisture_total.csv (data/raw/field_data/ipad_data/Cleaned data/), an older
# compilation of the same visits, is checked record by record (section 2b): each of
# its 296 records is either empty or a reading already here. It is never a source of
# readings; a record that matched nothing stops the build.
# Its errors, for the record: 63 readings dated a day early (2020-10-27, 2020-11-05,
# 2021-05-18; the sheets and the fluxes are -28, -06, -19), 22 of them under another
# collar's label (incl. "1-31", "1-0"); its 12 records dated 2020-07-21 are the
# 2020-06-17 readings of the same collars; at six 2020-08-26 collars its VWC column
# holds the probe's raw reading. Every record's match:
# outputs/audit/soil_readings_soilmoisture_total.csv.
# The iPad sheets "Soil temp moisture 20200722.xlsx" and "...20200826.xlsx" were copies
# of the 20200615 sheet (date cell 2020-06-15; the six U_1 readings identical on all
# three dates). Their readings are blanked in the raw files and every compilation of
# them (data/raw/CORRECTIONS.md, 2026-10-16), so U_1 has no reading on those dates.
#
# vwc_source_column: which sheet column soil_vwc_pct was read from. The iPad sheets
# carry two moisture columns ("Soil moisture" and "VWC"); per sheet, the one whose
# values look like a percentage (1.6-80) for more readings is taken (the other
# holds the probe's raw reading); the column chosen is recorded per record.
#
# COLLARS, 2021-03-30: the soil temperatures (the flux sheet's ST column) and the six
# intermediate-collar VWC values (typed in the flux sheet's Notes column) are in that
# day's soil sheet since 2026-10-17 (data/raw/CORRECTIONS.md), so it is read like any
# other date.
#
# TREES, 2023 (data/raw/field_data/static_chamber_field/Compiled YMF Data 2023 -
# Sheet1.csv): one reading per sheet row at the tree measured (tree_tags, with species
# and DBH where a tag names two trees): Soil Temp (degC) and the mean of the VWC
# probes (VWC 1-3; "Saturated" = 60 %, as for the 2021 holes); time_local is the
# real flux start where written. Rows whose tree is not a measured tree are not
# carried (listed in the build report).
#
# Soil temperatures must lie in -5 to 35 C (the build stops otherwise). The one reading
# outside it (44.6 C at tree 748 on 2023-07-07) was blanked in the raw sheet
# (data/raw/CORRECTIONS.md, 2026-10-17).
#
# Checks: tree means reproduce VWC_mean, ORP_mean, SoilTemp_mean,
# OrganicDepth_mean, MineralDepth_mean of soil_tree_level_means.csv (a combined name's
# row against each of its trees).
#
# Writes data/tables/soil_readings.csv, outputs/audit/soil_readings_build.txt,
#   outputs/audit/soil_readings_soilmoisture_total.csv
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages({ library(tidyverse); library(readxl) })
OUT <- "data/tables"
SRC_T <- "data/raw/field_data/static_chamber_field/Soils_Data_filled.csv"
MEANS <- "data/raw/field_data/static_chamber_field/soil_tree_level_means.csv"
SATURATED_VWC <- 60
num <- function(x) suppressWarnings(as.numeric(x))
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-9 * pmax(1, abs(a)))
mean_na <- function(x) { m <- mean(x, na.rm = TRUE); ifelse(is.nan(m), NA_real_, m) }
TT_ALL <- read.csv(file.path(OUT, "tree_tags.csv"), colClasses = "character", na.strings = "")
tags_for <- function(file, column) {
  t <- TT_ALL %>% filter(source_file == file, source_column == column) %>% distinct(tag = trimws(tag_recorded), tree_id, .keep_all = TRUE)
  stopifnot(all(is.na(t$species_recorded)), !anyDuplicated(t$tag)); t
}

# ==============================================================================
# 1. TREE HOLES
# ==============================================================================
f <- read.csv(SRC_T, colClasses = "character", na.strings = character(0), check.names = FALSE)
stopifnot(nrow(f) == 600)
vw <- sapply(c("VWC1", "VWC2", "VWC3"), function(k) num(sub("^(\\d+),(\\d+)$", "\\1.\\2", trimws(f[[k]]))))
sat <- tolower(trimws(f$VWC1)) == "saturated"
min_txt <- trimws(f$`Mineral.Sample.Depth..cm`); dist_txt <- trimws(f$`Distance.from.stem..cm`)
H <- tibble(
  tree_recorded = trimws(f$Tree_ID),
  date = as.Date(trimws(f$Date.Collected), "%m/%d/%y"),
  time_local = ifelse(trimws(f$Time.Collected) == "", NA, sub(";", ":", trimws(f$Time.Collected))),
  hole = gsub("[^NESW]", "", toupper(f$Cardinal.Direction)),
  distance_from_stem_cm = num(dist_txt), hole_depth_cm = num(f$`Total.Hole.depth..cm`),
  organic_depth_cm = num(f$`Organic.Sample.Depth..cm`), mineral_depth_cm = num(min_txt),
  soil_temp_C = num(f$Soil.Temp),
  soil_vwc_pct = ifelse(sat, SATURATED_VWC, apply(vw, 1, mean_na)),
  vwc_n_probes = ifelse(sat, 0L, as.integer(rowSums(!is.na(vw)))),
  vwc_saturated = sat, vwc_source_column = "VWC1-VWC3 (mean of the probe readings)",
  orp_mV = num(f$`ORP..mV`),
  notes = pmap_chr(list(trimws(f$Notes), ifelse(min_txt %in% c(">45", "water"), paste0("mineral depth recorded as '", min_txt, "'"), NA),
                        ifelse(dist_txt == "*wet", "distance recorded as '*wet'", NA),
                        ifelse(grepl(",", f$VWC3), paste0("VWC3 recorded as '", trimws(f$VWC3), "'"), NA)),
                   function(...) { x <- c(...); x <- x[!is.na(x) & x != ""]; if (length(x)) paste(x, collapse = "; ") else NA_character_ }))
stopifnot(all(H$hole %in% c("N", "E", "S", "W")))

tm <- tags_for(MEANS, "Tree_ID")
TREES <- read.csv(file.path(OUT, "trees.csv"), colClasses = "character")
# holes recorded under a combined name: one row per tree (see the header)
SHARED_NAMES <- list("AB2/3" = c("AB2", "AB3"), "WA14/1/XV" = c("WA14", "WA1", "WAX1"))
stopifnot(all(unlist(SHARED_NAMES) %in% TREES$tree_id), all(names(SHARED_NAMES) %in% H$tree_recorded),
          !any(unlist(SHARED_NAMES) %in% tm$tree_id[match(H$tree_recorded, tm$tag)]))
member_of <- function(rec) { i <- match(rec, tm$tag); ifelse(is.na(i), NA_character_, tm$tree_id[i]) }
H <- bind_rows(H %>% filter(!tree_recorded %in% names(SHARED_NAMES)) %>% mutate(tree_id = member_of(tree_recorded)),
               bind_rows(lapply(names(SHARED_NAMES), function(nm)
                 H %>% filter(tree_recorded == nm) %>% slice(rep(seq_len(n()), each = length(SHARED_NAMES[[nm]]))) %>%
                   mutate(tree_id = rep(SHARED_NAMES[[nm]], length.out = n())))))
not_tree <- H %>% filter(is.na(tree_id) | !tree_id %in% TREES$tree_id)
stopifnot(nrow(not_tree) == 0)
# holes repeated under several trees
H <- H %>% group_by(date, time_local, hole, distance_from_stem_cm, hole_depth_cm, organic_depth_cm, mineral_depth_cm,
                    soil_temp_C, soil_vwc_pct, orp_mV) %>%
  mutate(shared_with = map_chr(tree_id, function(t) { o <- setdiff(unique(tree_id), t); if (length(o)) paste(sort(o), collapse = ";") else NA_character_ })) %>%
  ungroup() %>%
  group_by(tree_id, hole) %>% mutate(k = row_number(), n = n()) %>% ungroup() %>%
  mutate(reading_id = paste0("tree_", tree_id, "_", hole, ifelse(n > 1, paste0("_", k), "")), collar_id = NA_character_,
         source = basename(SRC_T))
# holes with no date: the date of the tree's 2021 flux visit (one visit per tree)
F21 <- "data/processed/flux/tree_flux_2021_multiheight.csv"
h21 <- read.csv(F21, check.names = FALSE, stringsAsFactors = FALSE)
t21 <- TT_ALL %>% filter(source_file == F21, source_column == "tree_id") %>% distinct(tag = trimws(tag_recorded), tree_id)
stopifnot(!anyDuplicated(t21$tag))
VIS <- tibble(tree_id = t21$tree_id[match(trimws(h21$tree_id), t21$tag)], date = as.Date(substr(h21$start.time, 1, 10))) %>%
  filter(!is.na(tree_id)) %>% distinct()
vis1 <- VIS %>% group_by(tree_id) %>% filter(n() == 1) %>% ungroup()
dated <- H %>% filter(!is.na(date)) %>% distinct(tree_id, date) %>% inner_join(VIS, by = "tree_id", suffix = c("", "_flux")) %>%
  group_by(tree_id, date) %>% summarise(gap = min(abs(as.numeric(date - date_flux))), .groups = "drop")
H$date_source <- ifelse(is.na(H$date), NA_character_, "recorded")
und <- is.na(H$date)
H$date[und] <- vis1$date[match(H$tree_id[und], vis1$tree_id)]
H$date_source[und & !is.na(H$date)] <- "not recorded: date of the tree's 2021 flux visit"
undated <- H %>% filter(und) %>% count(tree_id, date, name = "holes")

# ==============================================================================
# 2. COLLARS
# ==============================================================================
files <- list.files("data/raw/field_data/ipad_data", pattern = "^Soil temp moisture.*\\.xlsx$", recursive = TRUE, full.names = TRUE)
is_pct <- function(v) { v <- suppressWarnings(as.numeric(v)); !is.na(v) & v > 1.6 & v < 80 }
parse_one <- function(f) {
  sheets <- tryCatch(excel_sheets(f), error = function(e) character(0))
  for (s in sheets) {
    d <- suppressMessages(tryCatch(read_excel(f, sheet = s, col_names = FALSE), error = function(e) NULL))
    if (is.null(d) || nrow(d) < 3 || ncol(d) < 4) next
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    hdr <- which(apply(d, 1, function(r) any(grepl("soil\\s*temp", r, ignore.case = TRUE))))[1]
    if (is.na(hdr)) next
    nm <- tolower(trimws(as.character(unlist(d[hdr, ]))))
    body <- d[(hdr + 1):nrow(d), , drop = FALSE]
    names(body) <- ifelse(is.na(nm) | nm == "", paste0("v", seq_along(nm)), nm)
    col <- function(pat) { i <- grep(pat, names(body)); if (length(i)) body[[i[1]]] else NA }
    plot_raw <- col("^plot$"); letter_raw <- col("plot\\s*letter"); temp_raw <- col("soil\\s*temp")
    moist_a <- col("soil\\s*moisture"); moist_b <- col("^vwc")
    if (all(is.na(plot_raw))) next
    pick_moisture <- function(a, b) {
      a <- suppressWarnings(as.numeric(a)); b <- suppressWarnings(as.numeric(b))
      score <- function(v) if (all(is.na(v))) -1 else mean(is_pct(v), na.rm = TRUE)
      if (score(a) >= score(b)) "a" else "b"
    }
    which_m <- pick_moisture(moist_a, moist_b)
    moist <- suppressWarnings(as.numeric(if (which_m == "a") moist_a else moist_b))
    moist_name <- if (which_m == "a") names(body)[grep("soil\\s*moisture", names(body))[1]] else names(body)[grep("^vwc", names(body))[1]]
    pr <- trimws(as.character(plot_raw))
    combined <- grepl("^[A-Za-z]{1,2}-\\d+-\\d+$", pr)
    serial <- suppressWarnings(as.numeric(pr))
    is_ser <- !is.na(serial) & serial > 30000 & serial < 60000
    dec <- as.Date(serial, origin = "1904-01-01")
    pr_dec <- ifelse(is_ser, paste0(as.integer(format(dec, "%m")), "-", as.integer(format(dec, "%d"))), pr)
    letter <- ifelse(combined, toupper(sub("^([A-Za-z]{1,2})-.*", "\\1", pr)), toupper(trimws(as.character(letter_raw))))
    tag <- ifelse(combined, sub("^[A-Za-z]{1,2}-", "", pr), pr_dec)
    Date <- as.Date(sub(".*?(\\d{8}).*", "\\1", basename(f)), format = "%Y%m%d")
    out <- data.frame(Date = Date, plot_letter = letter, plot_tag = tag,
                      soil_temp_C = suppressWarnings(as.numeric(temp_raw)), soil_moisture_pct = moist,
                      vwc_source_column = paste0("'", moist_name, "' column"),
                      source = basename(f), stringsAsFactors = FALSE) %>%
      filter(!is.na(plot_letter), plot_letter %in% c("U", "I", "WD", "WS"), grepl("^\\d+-\\d+$", plot_tag),
             !(is.na(soil_temp_C) & is.na(soil_moisture_pct)))
    if (nrow(out) > 0) return(out)
  }
  NULL
}
env <- bind_rows(lapply(files, function(f) tryCatch(parse_one(f), error = function(e) NULL)))
MASTER <- "data/raw/environmental/Soil_temp_moisture_2020.xlsx"
m <- as.data.frame(read_excel(MASTER, sheet = 1, .name_repair = "minimal"))[, 1:6]
names(m) <- c("Date", "Site", "Plot", "Subplot", "Soil.temp", "VWC")
site_map <- c(Upland = "U", Intermediate = "I", "Wetland dry" = "WD", "Wetland swamp" = "WS")
env <- bind_rows(env, m %>%
  transmute(Date = as.Date(Date), plot_letter = unname(site_map[trimws(Site)]), plot_tag = paste0(Plot, "-", Subplot),
            soil_temp_C = suppressWarnings(as.numeric(Soil.temp)), soil_moisture_pct = suppressWarnings(as.numeric(VWC)),
            vwc_source_column = "'VWC' column", source = basename(MASTER)) %>%
  filter(!is.na(Date), !is.na(plot_letter), grepl("^\\d+-\\d+$", plot_tag), !(is.na(soil_temp_C) & is.na(soil_moisture_pct)))) %>%
  filter(!is.na(Date))

COLLARS <- read.csv(file.path(OUT, "collars.csv"), colClasses = "character")
C <- env %>% mutate(collar_id = paste0(plot_letter, "_", plot_tag))
stopifnot(all(C$collar_id %in% COLLARS$collar_id))
C <- C %>% transmute(collar_id, date = Date, soil_temp_C, soil_vwc_pct = soil_moisture_pct, vwc_source_column, source)
C$vwc_source_column[is.na(C$soil_vwc_pct)] <- NA

# ==============================================================================
# 2b. soilmoisture_total.csv -- each record matched to a reading above, in order:
#   empty (no temperature, no VWC); the same collar, date and values (NA = NA); the
#   same collar, date and VWC where the record has no temperature; the
#   same collar, date and temperature with the probe's raw reading (< 3) in the VWC
#   column; the same values on the same or the next day under any collar label; the
#   same values at the same collar on another date. An unmatched record stops the build.
# ==============================================================================
SMT_F <- "data/raw/field_data/ipad_data/Cleaned data/soilmoisture_total.csv"
SMT <- read.csv(SMT_F, check.names = FALSE, colClasses = "character", fileEncoding = "UTF-8-BOM")[, 1:6]
names(SMT) <- c("Date", "Plot", "letter", "temp", "moist_raw", "VWC")
SMT <- SMT %>% transmute(smt_row = row_number(), collar_id = paste0(toupper(trimws(letter)), "_", trimws(Plot)),
                         date = as.Date(Date, "%m/%d/%y"), soil_temp_C = num(temp), soil_vwc_pct = num(VWC))
stopifnot(nrow(SMT) == 296, !anyNA(SMT$date))
smt_match <- function(r) {
  if (is.na(r$soil_temp_C) && is.na(r$soil_vwc_pct)) return(c("empty record", NA, NA))
  eqT <- same(C$soil_temp_C, r$soil_temp_C); eqV <- same(C$soil_vwc_pct, r$soil_vwc_pct)
  rules <- list(
    "same collar, date and values" = C$collar_id == r$collar_id & C$date == r$date & eqT & eqV,
    "same collar, date and VWC; the record has no temperature" =
      C$collar_id == r$collar_id & C$date == r$date & is.na(r$soil_temp_C) & eqV,
    "same collar, date and temperature; its VWC column holds the raw probe reading" =
      C$collar_id == r$collar_id & C$date == r$date & eqT & isTRUE(r$soil_vwc_pct < 3),
    "same values a day later" = C$collar_id == r$collar_id & C$date == r$date + 1 & eqT & eqV,
    "same values a day later, under another collar label" = C$date == r$date + 1 & eqT & eqV,
    "same values the same day, under another collar label" = C$date == r$date & eqT & eqV,
    "same collar and values on another date" = C$collar_id == r$collar_id & eqT & eqV)
  for (k in names(rules)) { j <- which(rules[[k]]); if (length(j)) return(c(k, C$collar_id[j[1]], format(C$date[j[1]]))) }
  c("unmatched", NA, NA)
}
SM_M <- bind_cols(SMT, as_tibble(do.call(rbind, lapply(split(SMT, SMT$smt_row), smt_match)), .name_repair = ~ c("match", "match_collar", "match_date")))
stopifnot(all(C$collar_id %in% COLLARS$collar_id))
C <- C %>% group_by(collar_id, date) %>% mutate(k = row_number(), n = n()) %>% ungroup() %>%
  transmute(reading_id = paste0("collar_", collar_id, "_", date, ifelse(n > 1, paste0("_", k), "")),
            tree_id = NA_character_, collar_id, date, soil_temp_C, soil_vwc_pct, vwc_source_column, source,
            date_source = "recorded")

# ==============================================================================
# 2d. TREES, 2023 cross-species survey: one reading per field-sheet row
# ==============================================================================
F_23 <- "data/raw/field_data/static_chamber_field/Compiled YMF Data 2023 - Sheet1.csv"
y <- read.csv(F_23, check.names = FALSE, colClasses = "character", na.strings = character(0))
t23 <- TT_ALL %>% filter(source_file == F_23, source_column == "Tree Tag")
y_tree <- vapply(seq_len(nrow(y)), function(i) {
  c <- t23[trimws(t23$tag_recorded) == trimws(y$`Tree Tag`[i]), ]
  if (nrow(c) > 1) c <- c[c$species_recorded %in% y$`Species Code`[i] & abs(num(c$dbh_recorded_cm) - num(y$`DBH (cm)`[i])) < 0.05, ]
  if (nrow(c) == 1) c$tree_id else NA_character_ }, character(1))
v23 <- sapply(c("VWC 1 (%)", "VWC 2 (%)", "VWC 3 (%)"), function(k) {
  x <- trimws(y[[k]]); ifelse(tolower(x) == "saturated", SATURATED_VWC, num(x)) })
sat23 <- apply(sapply(c("VWC 1 (%)", "VWC 2 (%)", "VWC 3 (%)"), function(k) tolower(trimws(y[[k]])) == "saturated"), 1, any)
Y <- tibble(tree_id = y_tree, tag = y$`Tree Tag`, date = as.Date(y$Date, "%m/%d/%y"),
            time_local = ifelse(grepl("^\\d{1,2}[:;]\\d{2}", trimws(y$`Flux Start (Real) Time`)),
                                substr(gsub(";", ":", trimws(y$`Flux Start (Real) Time`)), 1, 5), NA_character_),
            soil_temp_C = num(y$`Soil Temp (°C)`), soil_vwc_pct = apply(v23, 1, mean_na),
            vwc_n_probes = as.integer(rowSums(!is.na(v23))), vwc_saturated = sat23,
            vwc_source_column = "VWC 1-3 (mean of the probe readings)", notes = NA_character_, source = basename(F_23),
            date_source = "recorded") %>%
  mutate(soil_vwc_pct = ifelse(is.nan(soil_vwc_pct), NA_real_, soil_vwc_pct))
stopifnot(!anyNA(Y$date))
y_not <- Y %>% filter(is.na(tree_id) | !tree_id %in% TREES$tree_id)
Y <- Y %>% filter(!is.na(tree_id), tree_id %in% TREES$tree_id, !(is.na(soil_temp_C) & is.na(soil_vwc_pct))) %>%
  group_by(tree_id, date) %>% mutate(k = row_number(), n = n()) %>% ungroup() %>%
  mutate(reading_id = paste0("tree_", tree_id, "_", date, ifelse(n > 1, paste0("_", k), "")), collar_id = NA_character_) %>%
  select(-k, -n, -tag)
SM_M <- SM_M %>% mutate(match_reading = ifelse(is.na(match_collar), NA, paste0("collar_", match_collar, "_", match_date)))
write_csv(SM_M %>% select(smt_row, collar_id, date, soil_temp_C, soil_vwc_pct, match, match_reading),
          out_path("soil_readings_soilmoisture_total.csv", "audit"), na = "")
if (any(SM_M$match == "unmatched")) stop("soilmoisture_total.csv holds readings not in the field sheets: see ",
                                         out_path("soil_readings_soilmoisture_total.csv", "audit"))
# sheets that repeat another date's readings (copied sheets)
rep_sets <- C %>% filter(!is.na(soil_temp_C)) %>% distinct(collar_id, date, soil_temp_C, soil_vwc_pct) %>%
  group_by(collar_id, soil_temp_C, soil_vwc_pct) %>% filter(n_distinct(date) > 1) %>%
  summarise(dates = paste(sort(unique(format(date))), collapse = " "), .groups = "drop")

# ==============================================================================
# 3. ASSEMBLE
# ==============================================================================
SR <- bind_rows(
  H %>% select(reading_id, tree_id, collar_id, date, date_source, time_local, hole, distance_from_stem_cm, hole_depth_cm,
               organic_depth_cm, mineral_depth_cm, soil_temp_C, soil_vwc_pct, vwc_n_probes, vwc_saturated, vwc_source_column,
               orp_mV, shared_with, source, notes),
  Y, C) %>%
  mutate(hole_depth_cm = as.numeric(hole_depth_cm)) %>%
  arrange(is.na(tree_id), tree_id, collar_id, date, hole)
SOIL_T_RANGE <- c(-5, 35)
bad_t <- which(!is.na(SR$soil_temp_C) & (SR$soil_temp_C < SOIL_T_RANGE[1] | SR$soil_temp_C > SOIL_T_RANGE[2]))
if (length(bad_t)) stop("soil temperatures outside ", SOIL_T_RANGE[1], " to ", SOIL_T_RANGE[2], " C (correct the raw file, ",
                        "data/raw/CORRECTIONS.md): ", paste(SR$reading_id[bad_t], SR$soil_temp_C[bad_t], collapse = "; "))
stopifnot(!anyDuplicated(SR$reading_id), all(xor(is.na(SR$tree_id), is.na(SR$collar_id))))

# ==============================================================================
# 4. CHECKS
# ==============================================================================
# (a) trees vs soil_tree_level_means.csv (the field team's tree means)
tmean <- SR %>% filter(!is.na(tree_id), source == basename(SRC_T)) %>% group_by(tree_id) %>%
  summarise(VWC_mean = mean_na(soil_vwc_pct), ORP_mean = mean_na(orp_mV), SoilTemp_mean = mean_na(soil_temp_C),
            OrganicDepth_mean = mean_na(organic_depth_cm), MineralDepth_mean = mean_na(mineral_depth_cm), .groups = "drop")
VARS <- c("VWC_mean", "ORP_mean", "SoilTemp_mean", "OrganicDepth_mean", "MineralDepth_mean")
SM <- read.csv(MEANS, stringsAsFactors = FALSE) %>% mutate(Tree_ID = trimws(Tree_ID))
SM <- bind_rows(SM %>% filter(!Tree_ID %in% names(SHARED_NAMES)) %>% mutate(tree_id = tm$tree_id[match(Tree_ID, tm$tag)]),
                bind_rows(lapply(names(SHARED_NAMES), function(nm) SM %>% filter(Tree_ID == nm) %>%
                  slice(rep(1, length(SHARED_NAMES[[nm]]))) %>% mutate(tree_id = SHARED_NAMES[[nm]])))) %>%
  filter(!is.na(tree_id), tree_id %in% TREES$tree_id)
diffs <- function(ref, label) {
  j <- full_join(ref %>% select(tree_id, all_of(VARS)), tmean, by = "tree_id", suffix = c("_today", "_table"))
  bind_rows(lapply(VARS, function(v) j %>% filter(!same(.data[[paste0(v, "_today")]], .data[[paste0(v, "_table")]])) %>%
    transmute(against = label, tree_id, variable = v, today = .data[[paste0(v, "_today")]], table = .data[[paste0(v, "_table")]])))
}
D2 <- diffs(SM, "soil_tree_level_means.csv")
# The one expected difference: RMS1's west hole ("aquatic, in stream", no probe reading)
# is recorded saturated (60 %) since 2026-10-17 (CORRECTIONS.md); the field team's
# mean is over the other three holes.
EXPECTED <- tibble(tree_id = "RMS1", variable = "VWC_mean")
unexpected <- D2 %>% anti_join(EXPECTED, by = c("tree_id", "variable"))

# collar-date means (for the report)
cm <- SR %>% filter(!is.na(collar_id)) %>% group_by(collar_id, date) %>%
  summarise(t1 = mean_na(soil_temp_C), m1 = mean_na(soil_vwc_pct), .groups = "drop")

sink(out_path("soil_readings_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("soil_readings: %d rows | trees: %d holes at %d trees | collars: %d records, %d collar-dates, %d collars, %d dates\n",
            nrow(SR), sum(!is.na(SR$tree_id)), n_distinct(SR$tree_id, na.rm = TRUE), sum(!is.na(SR$collar_id)),
            nrow(cm), n_distinct(cm$collar_id), n_distinct(cm$date)))
cat(sprintf("holes under a combined name, repeated under each tree: %s\n",
            paste(sprintf("%s -> %s", names(SHARED_NAMES), sapply(SHARED_NAMES, paste, collapse = ", ")), collapse = "; ")))
cat(sprintf("holes with a saturated VWC probe (set to %g %%): %d; holes shared by several trees: %d rows\n",
            SATURATED_VWC, sum(SR$vwc_saturated, na.rm = TRUE), sum(!is.na(SR$shared_with))))
cat("records per source:\n"); print(as.data.frame(SR %>% count(source)), row.names = FALSE)
cat(sprintf("\n2021 holes with no recorded date, dated to the tree's 2021 flux visit: %d holes at %d trees (%s)\n",
            sum(undated$holes), nrow(undated), paste(sprintf("%s %s", undated$tree_id, undated$date), collapse = "; ")))
cat(sprintf("  trees with dated holes vs their flux visit: %d tree-dates the same day, %d within 3 days, %d further\n",
            sum(dated$gap == 0), sum(dated$gap > 0 & dated$gap <= 3), sum(dated$gap > 3)))
cat(sprintf("2023 field sheet: %d readings carried; %d rows not carried (tree not a measured tree: %s)\n", nrow(Y), nrow(y_not),
            paste(unique(y_not$tag), collapse = ", ")))
cat("VWC column used, by source:\n"); print(as.data.frame(SR %>% count(source, vwc_source_column)), row.names = FALSE)
cat(sprintf("\n(a) tree means vs soil_tree_level_means.csv: %d trees there (in trees.csv), %d values differ\n", nrow(SM), nrow(D2)))
if (nrow(D2)) print(as.data.frame(D2), row.names = FALSE)
cat(sprintf("\n(b) soilmoisture_total.csv (checked only, never a source): %d records; %d unmatched; by match:\n", nrow(SM_M),
            sum(SM_M$match == "unmatched")))
print(as.data.frame(SM_M %>% count(match, date) %>% group_by(match) %>%
  summarise(n = sum(n), dates = paste(format(date), collapse = " "), .groups = "drop")), row.names = FALSE)
cat(sprintf("\nreadings repeated on several dates at the same collar (a copied sheet): %d\n", nrow(rep_sets)))
if (nrow(rep_sets)) print(as.data.frame(rep_sets), row.names = FALSE)
options(op); sink()
cat(readLines(out_path("soil_readings_build.txt", "audit")), sep = "\n")
if (nrow(unexpected)) stop("soil_readings does not reproduce today's means (outputs/audit/soil_readings_build.txt)")

write_csv(SR, file.path(OUT, "soil_readings.csv"), na = "")
