#!/usr/bin/env Rscript
source("code/lib/outputs.R")
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
# The tree means in soil_tree_level_means.csv (which reach the analysis through
# merged_tree_dataset_final.csv) are the means over a tree's holes. Readings at
# trees not in trees.csv (2021 tree list only: AB2/3, WA14/1/XV) are not carried.
#
# COLLARS (monthly campaign, 2020-21): one row per collar per date per source
# record, parsed exactly as code/03_merge/02_compile_soil_env.R (the parser is
# copied below): the per-date iPad sheets "Soil temp moisture YYYYMMDD*.xlsx"
# and the 2020-21 master workbook data/raw/environmental/Soil_temp_moisture_2020.xlsx.
# Where a collar-date appears in more than one record (both sources, or twice in
# a sheet) each record is a row; soil_env_by_collar.csv is their mean.
# soilmoisture_total.csv (data/raw/field_data/ipad_data/Cleaned data/), an older
# compilation of the same visits, is reconciled record by record (section 2b): each
# of its 296 records is either empty or a reading already here, so it adds no row.
# A record that matched nothing would be added (source soilmoisture_total.csv).
# Its errors, for the record: 63 readings dated a day early (2020-10-27, 2020-11-05,
# 2021-05-18; the sheets and the fluxes are -28, -06, -19), 22 of them under another
# collar's label (incl. "1-31", "1-0"); its 12 records dated 2020-07-21 are the
# 2020-06-17 readings of the same collars; at six 2020-08-26 collars its VWC column
# holds the probe's raw reading. Every record's match:
# outputs/audit/soil_readings_soilmoisture_total.csv.
# The iPad sheets "Soil temp moisture 20200722.xlsx" and "...20200826.xlsx" are copies
# of the 20200615 sheet (their date cell reads 2020-06-15; the six U_1 readings are
# identical on all three dates). They are carried as recorded; see the audit file.
#
# Checks: tree means reproduce VWC_mean, ORP_mean, SoilTemp_mean,
# OrganicDepth_mean, MineralDepth_mean of merged_tree_dataset_final.csv (and of
# soil_tree_level_means.csv); collar-date means reproduce soil_env_by_collar.csv.
#
# Writes data/tables/soil_readings.csv, outputs/audit/soil_readings_build.txt,
#   outputs/audit/soil_readings_soilmoisture_total.csv
# ==============================================================================
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
  vwc_saturated = sat,
  orp_mV = num(f$`ORP..mV`),
  notes = pmap_chr(list(trimws(f$Notes), ifelse(min_txt %in% c(">45", "water"), paste0("mineral depth recorded as '", min_txt, "'"), NA),
                        ifelse(dist_txt == "*wet", "distance recorded as '*wet'", NA),
                        ifelse(grepl(",", f$VWC3), paste0("VWC3 recorded as '", trimws(f$VWC3), "'"), NA)),
                   function(...) { x <- c(...); x <- x[!is.na(x) & x != ""]; if (length(x)) paste(x, collapse = "; ") else NA_character_ }))
stopifnot(all(H$hole %in% c("N", "E", "S", "W")))

tm <- tags_for(MEANS, "Tree_ID")
H$tree_id <- tm$tree_id[match(H$tree_recorded, tm$tag)]
TREES <- read.csv(file.path(OUT, "trees.csv"), colClasses = "character")
not_tree <- H %>% filter(is.na(tree_id) | !tree_id %in% TREES$tree_id)
stopifnot(setequal(unique(not_tree$tree_recorded), c("AB2/3", "WA14/1/XV")))
H <- H %>% filter(!is.na(tree_id), tree_id %in% TREES$tree_id)
# holes repeated under several trees
H <- H %>% group_by(date, time_local, hole, distance_from_stem_cm, hole_depth_cm, organic_depth_cm, mineral_depth_cm,
                    soil_temp_C, soil_vwc_pct, orp_mV) %>%
  mutate(shared_with = map_chr(tree_id, function(t) { o <- setdiff(unique(tree_id), t); if (length(o)) paste(sort(o), collapse = ";") else NA_character_ })) %>%
  ungroup() %>%
  group_by(tree_id, hole) %>% mutate(k = row_number(), n = n()) %>% ungroup() %>%
  mutate(reading_id = paste0("tree_", tree_id, "_", hole, ifelse(n > 1, paste0("_", k), "")), collar_id = NA_character_,
         source = basename(SRC_T))

# ==============================================================================
# 2. COLLARS -- parser copied from code/03_merge/02_compile_soil_env.R
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
      if (score(a) >= score(b)) a else b
    }
    moist <- pick_moisture(moist_a, moist_b)
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
            source = basename(MASTER)) %>%
  filter(!is.na(Date), !is.na(plot_letter), grepl("^\\d+-\\d+$", plot_tag), !(is.na(soil_temp_C) & is.na(soil_moisture_pct)))) %>%
  filter(!is.na(Date))

COLLARS <- read.csv(file.path(OUT, "collars.csv"), colClasses = "character")
C <- env %>% mutate(collar_id = paste0(plot_letter, "_", plot_tag))
stopifnot(all(C$collar_id %in% COLLARS$collar_id))
C <- C %>% transmute(collar_id, date = Date, soil_temp_C, soil_vwc_pct = soil_moisture_pct, source)

# ==============================================================================
# 2b. soilmoisture_total.csv -- each record matched to a reading above, in order:
#   empty (no temperature, no VWC); the same collar, date and values (NA = NA); the
#   same collar, date and temperature with the probe's raw reading (< 3) in the VWC
#   column; the same values on the same or the next day under any collar label; the
#   same values at the same collar on another date. Unmatched records are added.
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
    "same collar, date and temperature; its VWC column holds the raw probe reading" =
      C$collar_id == r$collar_id & C$date == r$date & eqT & isTRUE(r$soil_vwc_pct < 3),
    "same values a day later" = C$collar_id == r$collar_id & C$date == r$date + 1 & eqT & eqV,
    "same values a day later, under another collar label" = C$date == r$date + 1 & eqT & eqV,
    "same values the same day, under another collar label" = C$date == r$date & eqT & eqV,
    "same collar and values on another date" = C$collar_id == r$collar_id & eqT & eqV)
  for (k in names(rules)) { j <- which(rules[[k]]); if (length(j)) return(c(k, C$collar_id[j[1]], format(C$date[j[1]]))) }
  c("unmatched: added", NA, NA)
}
SM_M <- bind_cols(SMT, as_tibble(do.call(rbind, lapply(split(SMT, SMT$smt_row), smt_match)), .name_repair = ~ c("match", "match_collar", "match_date")))
C <- bind_rows(C, SM_M %>% filter(match == "unmatched: added") %>%
                 transmute(collar_id, date, soil_temp_C, soil_vwc_pct, source = basename(SMT_F)))
stopifnot(all(C$collar_id %in% COLLARS$collar_id))
C <- C %>% group_by(collar_id, date) %>% mutate(k = row_number(), n = n()) %>% ungroup() %>%
  transmute(reading_id = paste0("collar_", collar_id, "_", date, ifelse(n > 1, paste0("_", k), "")),
            tree_id = NA_character_, collar_id, date, soil_temp_C, soil_vwc_pct, source)
SM_M <- SM_M %>% mutate(match_reading = ifelse(is.na(match_collar), NA, paste0("collar_", match_collar, "_", match_date)))
write_csv(SM_M %>% select(smt_row, collar_id, date, soil_temp_C, soil_vwc_pct, match, match_reading),
          out_path("soil_readings_soilmoisture_total.csv", "audit"), na = "")
# sheets that repeat another date's readings (copied sheets)
rep_sets <- C %>% filter(!is.na(soil_temp_C)) %>% distinct(collar_id, date, soil_temp_C, soil_vwc_pct) %>%
  group_by(collar_id, soil_temp_C, soil_vwc_pct) %>% filter(n_distinct(date) > 1) %>%
  summarise(dates = paste(sort(unique(format(date))), collapse = " "), .groups = "drop")

# ==============================================================================
# 3. ASSEMBLE
# ==============================================================================
SR <- bind_rows(
  H %>% select(reading_id, tree_id, collar_id, date, time_local, hole, distance_from_stem_cm, hole_depth_cm,
               organic_depth_cm, mineral_depth_cm, soil_temp_C, soil_vwc_pct, vwc_n_probes, vwc_saturated, orp_mV,
               shared_with, source, notes),
  C) %>%
  mutate(hole_depth_cm = as.numeric(hole_depth_cm)) %>%
  arrange(is.na(tree_id), tree_id, collar_id, date, hole)
stopifnot(!anyDuplicated(SR$reading_id), all(xor(is.na(SR$tree_id), is.na(SR$collar_id))))

# ==============================================================================
# 4. CHECKS
# ==============================================================================
# (a) trees vs merged_tree_dataset_final.csv and soil_tree_level_means.csv
tmean <- SR %>% filter(!is.na(tree_id)) %>% group_by(tree_id) %>%
  summarise(VWC_mean = mean_na(soil_vwc_pct), ORP_mean = mean_na(orp_mV), SoilTemp_mean = mean_na(soil_temp_C),
            OrganicDepth_mean = mean_na(organic_depth_cm), MineralDepth_mean = mean_na(mineral_depth_cm), .groups = "drop")
VARS <- c("VWC_mean", "ORP_mean", "SoilTemp_mean", "OrganicDepth_mean", "MineralDepth_mean")
MT <- read.csv("data/processed/integrated/merged_tree_dataset_final.csv", check.names = FALSE, stringsAsFactors = FALSE)
tmt <- tags_for("data/processed/integrated/merged_tree_dataset_final.csv", "tree_id")
one <- function(x) { x <- unique(x[!is.na(x)]); if (length(x) > 1) stop("merged rows disagree"); if (length(x)) x else NA_real_ }
MTm <- MT %>% mutate(tree_id = tmt$tree_id[match(trimws(tree_id), tmt$tag)]) %>% filter(!is.na(tree_id)) %>%
  select(tree_id, all_of(VARS)) %>% group_by(tree_id) %>% summarise(across(everything(), one), .groups = "drop") %>%
  filter(if_any(all_of(VARS), ~ !is.na(.x)))
SM <- read.csv(MEANS, stringsAsFactors = FALSE) %>% mutate(tree_id = tm$tree_id[match(trimws(Tree_ID), tm$tag)]) %>%
  filter(!is.na(tree_id), tree_id %in% TREES$tree_id)
diffs <- function(ref, label) {
  j <- full_join(ref %>% select(tree_id, all_of(VARS)), tmean, by = "tree_id", suffix = c("_today", "_table"))
  bind_rows(lapply(VARS, function(v) j %>% filter(!same(.data[[paste0(v, "_today")]], .data[[paste0(v, "_table")]])) %>%
    transmute(against = label, tree_id, variable = v, today = .data[[paste0(v, "_today")]], table = .data[[paste0(v, "_table")]])))
}
D1 <- diffs(MTm, "merged_tree_dataset_final.csv"); D2 <- diffs(SM, "soil_tree_level_means.csv")
EXPECTED <- tibble(tree_id = character(), variable = character())  # none since the RO2 raw fix (CORRECTIONS.md, 2026-10-15)
unexpected <- bind_rows(D1, D2) %>% anti_join(EXPECTED, by = c("tree_id", "variable"))

# (b) collars vs soil_env_by_collar.csv
E <- read.csv("data/processed/environmental/soil_env_by_collar.csv", stringsAsFactors = FALSE) %>%
  transmute(collar_id = paste0(plot_letter, "_", plot_tag), date = as.Date(Date), t0 = soil_temp_C, m0 = soil_moisture_pct)
cm <- SR %>% filter(!is.na(collar_id)) %>% group_by(collar_id, date) %>%
  summarise(t1 = mean_na(soil_temp_C), m1 = mean_na(soil_vwc_pct), .groups = "drop")
cj <- full_join(E, cm, by = c("collar_id", "date"))
cbad <- cj %>% filter(!(same(t0, t1) & same(m0, m1)))

sink(out_path("soil_readings_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("soil_readings: %d rows | trees: %d holes at %d trees | collars: %d records, %d collar-dates, %d collars, %d dates\n",
            nrow(SR), sum(!is.na(SR$tree_id)), n_distinct(SR$tree_id, na.rm = TRUE), sum(!is.na(SR$collar_id)),
            nrow(cm), n_distinct(cm$collar_id), n_distinct(cm$date)))
cat(sprintf("tree holes not carried (tree not in trees.csv): %d (%s)\n", nrow(not_tree), paste(unique(not_tree$tree_recorded), collapse = ", ")))
cat(sprintf("holes with a saturated VWC probe (set to %g %%): %d; holes shared by several trees: %d rows\n",
            SATURATED_VWC, sum(SR$vwc_saturated, na.rm = TRUE), sum(!is.na(SR$shared_with))))
cat(sprintf("collar records per source: %s\n", paste(names(table(C$source == basename(MASTER))), table(C$source == basename(MASTER)), collapse = " ")))
cat(sprintf("\n(a) tree means vs merged_tree_dataset_final.csv: %d trees there, %d trees here, %d values differ\n",
            nrow(MTm), nrow(tmean), nrow(D1)))
cat(sprintf("    tree means vs soil_tree_level_means.csv:     %d trees there (in trees.csv), %d values differ\n", nrow(SM), nrow(D2)))
if (nrow(D1) + nrow(D2)) print(as.data.frame(bind_rows(D1, D2)), row.names = FALSE)
cat(sprintf("\n(b) collar-date means vs soil_env_by_collar.csv: %d there, %d here, %d differ\n", nrow(E), nrow(cm), nrow(cbad)))
if (nrow(cbad)) print(as.data.frame(cbad), row.names = FALSE)
cat(sprintf("\n(c) soilmoisture_total.csv: %d records in; %d rows added (unmatched); by match:\n", nrow(SM_M),
            sum(SM_M$match == "unmatched: added")))
print(as.data.frame(SM_M %>% count(match, date) %>% group_by(match) %>%
  summarise(n = sum(n), dates = paste(format(date), collapse = " "), .groups = "drop")), row.names = FALSE)
cat("\nreadings repeated on several dates at the same collar (copied sheets: 20200722 and 20200826 are the 20200615 sheet):\n")
print(as.data.frame(rep_sets), row.names = FALSE)
options(op); sink()
cat(readLines(out_path("soil_readings_build.txt", "audit")), sep = "\n")
if (nrow(unexpected) || nrow(cbad)) stop("soil_readings does not reproduce today's means (outputs/audit/soil_readings_build.txt)")

write_csv(SR, file.path(OUT, "soil_readings.csv"), na = "")
