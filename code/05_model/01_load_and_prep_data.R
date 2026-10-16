# ==============================================================================
# 01_load_and_prep_data.R -- random-forest workflow inputs from the canonical tables
# ------------------------------------------------------------------------------
# Reads data/tables/ only (built earlier in stage B by code/tables/) plus code/lib:
#   species, site, collars, inventory_stems, trees, tree_tags, flux_stem,
#   flux_soil, internal_gas, weather, moisture_survey, stream_points.
#
# Builds the PLACEHOLDER_* objects 02_rf_models.R reads:
#   TREE_JULY   2021 multi-height campaign + 2023 cross-species campaign
#   TREE_YEAR   monthly survey 2020-21 (semi-rigid chambers)
#   SOIL_YEAR   monthly soil collars 2020-21
#   INVENTORY   2019 census stems in the censused core, with lon/lat
#   MOISTURE_DEC_RASTER / MOISTURE_LOOKUP_XY   December 2020 moisture surface
#   DRIVERS     monthly air and soil temperature
#   TAXONOMY, PLOT_AREA, geodetic functions and constants
#
# Environmental covariates are the per-deployment values of flux_stem / flux_soil
# (air_temp_C, soil_temp_C, soil_vwc; their env_source column says where each came
# from). Run from code/05_model/ (pipeline workdir = script).
#
# "# KNOWN:" marks current behaviour that is a candidate for revision.
#
# Output: data/processed/integrated/rf_workflow_input_data_with_2023.RData
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(akima)
})

cat("=== RF WORKFLOW DATA IMPORT (data/tables) ===\n")
cat("Coordinate System: WGS84 (x=longitude, y=latitude)\n")
cat("Started at:", format(Sys.time()), "\n\n")

TAB <- "../../data/tables"
source("../lib/geometry.R")         # geodetic_transform(), geodetic_inverse()
source("../lib/species_levels.R")   # MODEL_SPECIES

# Every column as text. Numbers: the tables hold write_csv's round-trip decimal
# text, which a correctly rounded parser (jsonlite, i.e. C strtod) reads back to the
# double that was written; R's and readr's parsers misread some long decimals in the
# last digits.
rd <- function(t) read_csv(file.path(TAB, paste0(t, ".csv")), col_types = cols(.default = col_character()),
                           na = "", trim_ws = FALSE, progress = FALSE)
num <- function(x) {
  out <- rep(NA_real_, length(x))
  js <- !is.na(x) & grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", x)
  if (any(js)) out[js] <- as.double(jsonlite::parse_json(paste0("[", paste(x[js], collapse = ","), "]"),
                                                         simplifyVector = TRUE))
  oth <- !is.na(x) & !js
  if (any(oth)) out[oth] <- as.numeric(x[oth])   # Inf, NaN
  out
}
wall <- function(x) as.POSIXct(substr(sub("T", " ", x), 1, 19), tz = "UTC")   # local wall-clock digits

SPECIES  <- rd("species")
SITE     <- rd("site")
COLLARS  <- rd("collars")
STEMS    <- rd("inventory_stems")
TREES    <- rd("trees")
TAGS     <- rd("tree_tags")
FLUX     <- rd("flux_stem")
FLUX_S   <- rd("flux_soil")
GAS      <- rd("internal_gas")
WEATHER  <- rd("weather")
MOIST    <- rd("moisture_survey")
STREAM   <- rd("stream_points")
stopifnot(nrow(SITE) == 1)
# KNOWN: CH4 flux (flux_stem, flux_soil) and the 2023 field-sheet covariates are
#   readr's parse of the processed flux files, as the model has always read them;
#   readr misreads some long decimals (up to ~1e-13 relative).

# tree spellings by source file (tree_tags); the model's tree_id keeps these spellings
F_MERGED <- "data/processed/integrated/merged_tree_dataset_final.csv"
F_H21    <- "data/processed/flux/tree_flux_2021_multiheight.csv"
F_AUX21  <- "data/processed/flux/goflux_auxfile.csv"
F_23     <- "data/processed/flux/tree_flux_2023_cross_species.csv"
F_MON    <- "data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv"
F_YM     <- "data/raw/inventory/spatial_data/YM_trees_measured.csv.csv"
spellings <- function(file, column) TAGS[TAGS$source_file == file & TAGS$source_column == column, ]
# one spelling per tree; where a tree has several, the one equal (case and spaces
# aside) to its spelling in `prefer`
spelling_of <- function(tree, sp, prefer = NULL) {
  norm <- function(z) tolower(gsub(" ", "", z))
  vapply(tree, function(t) {
    c <- unique(sp$tag_recorded[sp$tree_id %in% t])
    if (length(c) > 1 && !is.null(prefer)) c <- c[norm(c) %in% norm(prefer$tag_recorded[prefer$tree_id %in% t])]
    if (length(c) != 1) stop("no unique spelling for tree ", t, call. = FALSE)
    c
  }, character(1), USE.NAMES = FALSE)
}

# =============================================================================
# 1. TAXONOMY
# =============================================================================
TAXONOMY <- data.frame(
  species_code = MODEL_SPECIES,
  species = SPECIES$species[match(MODEL_SPECIES, SPECIES$species_code)],
  genus   = SPECIES$genus[match(MODEL_SPECIES, SPECIES$species_code)],
  family  = SPECIES$family[match(MODEL_SPECIES, SPECIES$species_code)],
  order   = SPECIES$order[match(MODEL_SPECIES, SPECIES$species_code)],
  class   = SPECIES$class[match(MODEL_SPECIES, SPECIES$species_code)])
stopifnot(!anyNA(TAXONOMY))
species_mapping <- setNames(TAXONOMY$species, TAXONOMY$species_code)
cat("Taxonomy:", nrow(TAXONOMY), "species\n")

# =============================================================================
# 2. GEODETIC CONSTANTS (site.csv)
# =============================================================================
ref_lat <- num(SITE$ref_latitude)
ref_lon <- num(SITE$ref_longitude)
rotation_angle <- num(SITE$rotation_rad)
earth_radius <- num(SITE$earth_radius_m)
lat1_rad <- ref_lat * pi/180
lon1_rad <- ref_lon * pi/180
cat("Reference point:", ref_lat, "N,", ref_lon, "W; rotation", rotation_angle * 180/pi, "degrees\n")

# =============================================================================
# 3. INVENTORY: 2019 census stems in the censused core
# =============================================================================
# KNOWN: the model's inventory is the 2019 census only (not the 2018 by-tag survey),
#   first record per tag, raw DBH, species restricted to MODEL_SPECIES, censused core
#   PX/PY 0-200 m; PLOT_AREA is the nominal 200 x 200 m. 02_rf_models.R scales with
#   canonical_inventory() and STAND_AREA_M2 instead; this object feeds the tag joins
#   below. (The two 2019 records the census notes flag as duplicates, tags 3105 and
#   3547, are not in inventory_stems.)
fg19 <- STEMS %>% filter(source == "fg19") %>%
  transmute(Tag = as.integer(tag), PX = num(PX), PY = num(PY), DBH = num(dbh_cm_raw), Species_Code = species_code)
fg19_filtered <- fg19 %>% filter(!is.na(PX), !is.na(PY), !is.na(DBH))
coords_result <- geodetic_transform(fg19_filtered$PX, fg19_filtered$PY,
                                    lat1_rad, lon1_rad, rotation_angle, earth_radius)
fg19_clean <- fg19_filtered %>%
  mutate(
    x = coords_result$longitude,
    y = coords_result$latitude,
    dbh_m = DBH / 100,
    tree_id = as.character(Tag),
    species_code = Species_Code,
    species = species_mapping[Species_Code]
  ) %>%
  filter(!is.na(species)) %>%
  group_by(tree_id) %>%
  slice_head(n = 1) %>%
  ungroup()

INVENTORY <- fg19_clean %>%
  select(tree_id, species, species_code, dbh_m, x, y, PX, PY) %>%
  filter(!is.na(x), !is.na(y))

CENSUS_PX_MAX <- num(SITE$plot_side_m)
CENSUS_PY_MAX <- num(SITE$plot_side_m)
n_before_core <- nrow(INVENTORY)
INVENTORY <- INVENTORY %>% filter(PX >= 0, PX <= CENSUS_PX_MAX, PY >= 0, PY <= CENSUS_PY_MAX)
PLOT_AREA <- CENSUS_PX_MAX * CENSUS_PY_MAX
cat(sprintf("Inventory: %d stems in the censused core (%d outside dropped); plot area %.0f m2\n",
            nrow(INVENTORY), n_before_core - nrow(INVENTORY), PLOT_AREA))

# =============================================================================
# 4. WEATHER (monthly air-temperature drivers)
# =============================================================================
# KNOWN: weather timestamps are the logger clock as recorded (weather.csv
#   timestamp_logger; about UTC-7 before the 2021-08-03 reset, UTC-4 after).
weather_clean <- data.frame(TIMESTAMP = as.POSIXct(sub("T", " ", WEATHER$timestamp_logger), tz = "UTC"),
                            air_temp_C = num(WEATHER$air_temp_avg_C)) %>%
  filter(!is.na(air_temp_C)) %>%
  arrange(TIMESTAMP)
cat("Weather:", nrow(weather_clean), "records,", format(min(weather_clean$TIMESTAMP), "%Y-%m-%d"),
    "to", format(max(weather_clean$TIMESTAMP), "%Y-%m-%d"), "\n\n")

# =============================================================================
# 5. 2021 MULTI-HEIGHT CAMPAIGN (rigid chambers; 50, 125, 200 cm; Kalmia 75 cm)
# =============================================================================
# KNOWN: one row per tree x chamber height, the mean of that tree's deployments at
#   that height (the per-tree layout of merged_tree_dataset_final.csv), not one row
#   per deployment.
# KNOWN: tree_id is the 2021 survey spelling of the tree (merged-file spelling,
#   spaces removed), not trees.tree_id. Where a tree has two survey spellings
#   (SA1, WP6, WP8) the one its fluxes were filed under; soil readings filed under
#   the other spelling are not joined (flux_stem soil_temp_C / soil_vwc NA).
h21 <- FLUX %>% filter(campaign == "2021_multiheight") %>%
  transmute(tree = tree_id,
            hcol = ifelse(height_label %in% "root crown", "RootCrown", as.character(num(height_cm))),
            flux = num(CH4_flux_nmol_m2_s), air = num(air_temp_C),
            soil_temp_C = num(soil_temp_C), soil_moisture_abs = num(soil_vwc),
            start = wall(start_time))
stopifnot(all(h21$hcol %in% c("50", "75", "125", "200", "RootCrown")))
h21$spelling <- spelling_of(h21$tree, spellings(F_MERGED, "tree_id"), spellings(F_H21, "tree_id"))

# KNOWN: the per-height means pass through the merged file's CSV text and R's own
#   number parser (write_csv, then read.csv), reproduced here to the last bit.
csv_round_trip <- function(x)
  as.numeric(strsplit(format_csv(data.frame(v = x), col_names = FALSE), "\n", fixed = TRUE)[[1]])

# KNOWN: species is the merged table's: the code on the 2021 gas sample sheet
#   (internal_gas species_recorded), else one inferred from the name's prefix (its
#   ddPCR-sheet fallback adds nothing for these trees), not trees.species_code
#   (rm105: ACRU vs QURU; 7809, a numeric name without a gas sample, has none and
#   becomes "Unknown NA").
sp21 <- GAS %>% filter(campaign == "2021_survey", !is.na(tree_id)) %>% distinct(tree_id, species_recorded)
stopifnot(!anyDuplicated(sp21$tree_id))
NAME_PREFIX_SPECIES <- c(ab = "FAGR", bb = "BELE", bc = "PRSE", bo = "QUVE", h = "TSCA", ml = "KALA",
                         pb = "BEPA", rm = "ACRU", ro = "QURU", sa = "SAAL", sf = "SAAL", sh = "CAOV",
                         sm = "ACSA", wa = "FRAM", wo = "QUAL", wp = "PIST", yb = "BEAL")
species_from_name <- function(tree_id) vapply(tree_id, function(t) {
  if (t == "bowser") return("QUVE")
  if (t == "cul8r") return("BEAL")
  if (t == "sassafrass") return("SAAL")
  if (grepl("^[0-9]+$", t)) return(NA_character_)
  if (grepl("^h[0-9x-z]", tolower(t))) return("TSCA")
  unname(NAME_PREFIX_SPECIES[tolower(substr(t, 1, 2))])
}, character(1), USE.NAMES = FALSE)
# KNOWN: DBH is the 2021 survey consensus diameter only (trees.dbh_cm where
#   dbh_source is "2021 survey DBH (...)"); trees without one get NA rather than the
#   census or datasheet value.
dbh21 <- TREES %>% filter(grepl("^2021 survey DBH", dbh_source)) %>% transmute(tree_id, dbh = num(dbh_cm))

per_height <- h21 %>% group_by(spelling, tree, hcol) %>%
  summarise(flux = mean(flux, na.rm = TRUE), air = first(air), .groups = "drop") %>%
  mutate(flux = csv_round_trip(flux))
merged_tree_data <- h21 %>% distinct(spelling, tree, soil_temp_C, soil_moisture_abs)
stopifnot(!anyDuplicated(merged_tree_data$spelling))
for (h in c("50", "75", "125", "200", "RootCrown")) {
  p <- per_height[per_height$hcol == h, ]
  merged_tree_data[[paste0("CH4_best.flux_", h, "cm")]] <- p$flux[match(merged_tree_data$spelling, p$spelling)]
  merged_tree_data[[paste0("Temp_Air_", h, "cm")]] <- p$air[match(merged_tree_data$spelling, p$spelling)]
}
merged_tree_data <- merged_tree_data %>%
  mutate(species_id = coalesce(sp21$species_recorded[match(tree, sp21$tree_id)], species_from_name(spelling)),
         dbh = dbh21$dbh[match(tree, dbh21$tree_id)]) %>%
  arrange(spelling) %>%
  mutate(tree_id = gsub(" ", "", spelling)) %>%
  as.data.frame()

# Measurement dates: each tree's first deployment, matched by normalised name.
# KNOWN: dates are matched between the survey spelling and the goFlux aux-file
#   spelling after normalisation; six trees whose spellings differ (BOWSER, H7, PB5,
#   PB501, WP403, WPBERNICE) fall back to 2021-07-30 00:00, and every row of a tree
#   carries the tree's first deployment time.
tree_key <- function(z) toupper(gsub("[^A-Za-z0-9]", "",
                       gsub("prime", "P", gsub("'", "P", as.character(z)))))
first_start <- tapply(as.numeric(h21$start), h21$tree, min)
aux <- spellings(F_AUX21, "tree_id")
aux$first <- first_start[aux$tree_id]
aux <- aux[!is.na(aux$first), ]
JULY_DATES <- tapply(aux$first, tree_key(aux$tag_recorded), min)
base_date <- as.POSIXct("2021-07-30", tz = "UTC")
july_unmatched <- new.env(parent = emptyenv()); july_unmatched$ids <- character(0)
july_date <- function(tree_id) {
  k <- tree_key(tree_id)
  out <- rep(base_date, length(k))
  hit <- k %in% names(JULY_DATES)
  out[hit] <- as.POSIXct(JULY_DATES[k[hit]], origin = "1970-01-01", tz = "UTC")
  if (any(!hit)) july_unmatched$ids <- union(july_unmatched$ids, unique(k[!hit]))
  out
}

# KNOWN: a survey DBH of 3 or less is taken to be in metres already (dbh > 3 -> /100).
make_height_block <- function(df, flux_col, temp_col, height_cm) {
  df %>%
    filter(!is.na(.data[[flux_col]])) %>%
    mutate(
      Date = july_date(tree_id),
      tree_id = as.character(tree_id),
      species_code = species_id,
      species = ifelse(species_id %in% names(species_mapping),
                       species_mapping[species_id], paste("Unknown", species_id)),
      stem_flux_umol_m2_s = .data[[flux_col]],
      air_temp_C = .data[[temp_col]],
      soil_temp_C = soil_temp_C,
      soil_moisture_abs = soil_moisture_abs,
      dbh_m = case_when(is.na(dbh) ~ NA_real_, dbh > 3 ~ dbh / 100, TRUE ~ dbh),
      chamber_type = "rigid",
      measurement_height_cm = height_cm,
      month = month(Date),
      year = year(Date)
    ) %>%
    select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
           air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m, chamber_type,
           measurement_height_cm)
}
TREE_JULY_125 <- make_height_block(merged_tree_data, "CH4_best.flux_125cm", "Temp_Air_125cm", 125)
TREE_JULY_50  <- make_height_block(merged_tree_data, "CH4_best.flux_50cm",  "Temp_Air_50cm",  50)
TREE_JULY_200 <- make_height_block(merged_tree_data, "CH4_best.flux_200cm", "Temp_Air_200cm", 200)

# Kalmia: 75 cm, else root crown
KALMIA_DATA <- merged_tree_data %>%
  filter(species_id == "KALA" | grepl("(?i)kalmia", species_id)) %>%
  filter(!is.na(CH4_best.flux_75cm) | !is.na(CH4_best.flux_RootCrowncm)) %>%
  mutate(
    Date = july_date(tree_id),
    tree_id = as.character(tree_id),
    stem_flux_umol_m2_s = coalesce(CH4_best.flux_75cm, CH4_best.flux_RootCrowncm),
    air_temp_C = coalesce(Temp_Air_75cm, Temp_Air_RootCrowncm),
    soil_temp_C = soil_temp_C,
    soil_moisture_abs = soil_moisture_abs,
    dbh_m = ifelse(!is.na(dbh) & dbh > 0, dbh / 100, NA_real_),
    species_code = "KALA",
    species = "Kalmia latifolia",
    chamber_type = "rigid",
    month = month(Date),
    year = year(Date)
  ) %>%
  select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
         air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m, chamber_type)

# other trees measured only at 75 cm
OTHER_75CM <- merged_tree_data %>%
  filter(is.na(CH4_best.flux_125cm), !is.na(CH4_best.flux_75cm),
         !(species_id %in% c("KALA", "Kalmia"))) %>%
  mutate(
    Date = july_date(tree_id),
    tree_id = as.character(tree_id),
    species_code = species_id,
    species = ifelse(species_id %in% names(species_mapping),
                     species_mapping[species_id], paste("Unknown", species_id)),
    stem_flux_umol_m2_s = CH4_best.flux_75cm,
    air_temp_C = Temp_Air_75cm,
    soil_temp_C = soil_temp_C,
    soil_moisture_abs = soil_moisture_abs,
    dbh_m = case_when(is.na(dbh) ~ NA_real_, dbh > 3 ~ dbh / 100, TRUE ~ dbh),
    chamber_type = "rigid",
    month = month(Date),
    year = year(Date)
  ) %>%
  select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
         air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m, chamber_type)

TREE_JULY_2021 <- bind_rows(TREE_JULY_125, TREE_JULY_50, TREE_JULY_200, KALMIA_DATA, OTHER_75CM)
cat(sprintf("2021 campaign: %d tree x height rows (125 cm %d, 50 cm %d, 200 cm %d, Kalmia %d, other 75 cm %d)\n",
            nrow(TREE_JULY_2021), nrow(TREE_JULY_125), nrow(TREE_JULY_50), nrow(TREE_JULY_200),
            nrow(KALMIA_DATA), nrow(OTHER_75CM)))
cat(sprintf("  dates %s to %s, %d field days\n", format(min(TREE_JULY_2021$Date), "%Y-%m-%d"),
            format(max(TREE_JULY_2021$Date), "%Y-%m-%d"), length(unique(as.Date(TREE_JULY_2021$Date)))))
if (length(july_unmatched$ids))
  cat(sprintf("  %d tree(s) without a matching aux-file name, dated %s: %s\n", length(july_unmatched$ids),
              format(base_date, "%Y-%m-%d"), paste(july_unmatched$ids, collapse = ", ")))

# =============================================================================
# 6. 2023 CROSS-SPECIES CAMPAIGN (rigid chambers, 125 cm)
# =============================================================================
# KNOWN: tree_id is the field-sheet tag as a number, or "2023_" + the cleaned label
#   for ash and untagged stems (untagged numbered in sheet order); a tag written on
#   two trees (1985) gives both the same id.
# KNOWN: species is trees.species_code (equal to the field sheet's code for every
#   2023 deployment); DBH is the census DBH where the tag is a 2019 census tag, else
#   the field sheet's (flux_stem dbh_recorded_cm).
y23 <- FLUX %>% filter(campaign == "2023_cross_species")
tree_2023_clean <- y23 %>%
  transmute(
    tag = spelling_of(tree_id, spellings(F_23, "Tree Tag")),
    tree_tag_num = suppressWarnings(as.numeric(tag)),
    species_code = TREES$species_code[match(tree_id, TREES$tree_id)],
    date_clean = as.POSIXct(date, tz = "UTC"),
    month = month(date_clean),
    year = year(date_clean),
    stem_flux_umol_m2_s = num(CH4_flux_nmol_m2_s),
    dbh_own_m = num(dbh_recorded_cm) / 100,
    soil_moisture_own = num(soil_vwc),
    soil_temp_own = num(soil_temp_C),
    # KNOWN: where the sheet has no air temperature, flux_stem air_temp_C is the
    #   tower reading nearest 00:00 of the date.
    air_temp_C = num(air_temp_C)) %>%
  filter(!is.na(stem_flux_umol_m2_s)) %>%
  mutate(tree_id_2023 = ifelse(!is.na(tree_tag_num), as.character(tree_tag_num),
                               paste0("2023_", gsub("[^A-Za-z0-9]", "", tag))),
         tree_id_2023 = ifelse(duplicated(tree_id_2023) | duplicated(tree_id_2023, fromLast = TRUE),
                               ifelse(is.na(tree_tag_num), paste0(tree_id_2023, "_", ave(seq_along(tree_id_2023), tree_id_2023, FUN = seq_along)),
                                      tree_id_2023), tree_id_2023))

# KNOWN: joined to the census by tag number; a tag shared with another census tree
#   takes that tree's position and DBH.
TREE_JULY_2023 <- tree_2023_clean %>%
  mutate(tree_id = tree_id_2023) %>%
  left_join(INVENTORY %>% select(tree_id, x, y, dbh_m), by = "tree_id") %>%
  mutate(
    species = species_mapping[species_code],
    Date = date_clean,
    chamber_type = "rigid",
    dbh_m = coalesce(dbh_m, dbh_own_m)
  ) %>%
  filter(!is.na(species), !is.na(dbh_m)) %>%
  mutate(soil_temp_C = soil_temp_own, soil_moisture_abs = soil_moisture_own) %>%
  select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
         air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m, chamber_type, x, y)
cat(sprintf("2023 campaign: %d deployments (%d without a census record)\n", nrow(TREE_JULY_2023),
            sum(!tree_2023_clean$tree_id_2023 %in% INVENTORY$tree_id)))

# =============================================================================
# 7. TREE_JULY = 2021 + 2023
# =============================================================================
TREE_JULY_2021 <- TREE_JULY_2021 %>% left_join(INVENTORY %>% select(tree_id, x, y), by = "tree_id")
TREE_JULY <- bind_rows(TREE_JULY_2021, TREE_JULY_2023) %>% arrange(Date)
n_2021 <- sum(TREE_JULY$year == 2021, na.rm = TRUE)
n_2023 <- sum(TREE_JULY$year == 2023, na.rm = TRUE)
cat(sprintf("TREE_JULY: %d rows from %d trees (2021 %d, 2023 %d); with coordinates %d\n\n", nrow(TREE_JULY),
            length(unique(TREE_JULY$tree_id)), n_2021, n_2023, sum(!is.na(TREE_JULY$x) & !is.na(TREE_JULY$y))))

# =============================================================================
# 8. TREE_YEAR: monthly survey 2020-21 (semi-rigid chambers, 125 cm)
# =============================================================================
# KNOWN: tree_id is the survey's Plot Tag (or its digits where only those match a
#   census tag); census species, DBH and position joined by tag number.
# KNOWN: the 45 rescued untagged deployments are dated at 12:00 of their date (their
#   flux_stem air_temp_C is the tower reading nearest that noon), the others at the
#   logged start.
# KNOWN: soil temperature and VWC are flux_stem's: the collar readings averaged by
#   date and landscape position, nearest date within 3 days.
mon <- FLUX %>% filter(campaign == "monthly_2020_2021")
TREE_YEAR <- mon %>%
  transmute(
    measurement_id,
    Date = if_else(untagged_rescue == "TRUE", as.POSIXct(paste(date, "12:00:00"), tz = "UTC"), wall(start_time)),
    tree_id_raw = spelling_of(tree_id, spellings(F_MON, "Plot.Tag")),
    stem_flux_umol_m2_s = num(CH4_flux_nmol_m2_s),
    chamber_type = "semirigid",
    month = month(Date),
    year = year(Date),
    air_temp_C = num(air_temp_C),
    soil_temp_C = num(soil_temp_C),
    soil_moisture_abs = num(soil_vwc))

TREE_YEAR_match1 <- TREE_YEAR %>%
  mutate(tree_id = tree_id_raw) %>%
  left_join(INVENTORY %>% select(tree_id, species, species_code, dbh_m, x, y), by = "tree_id")
TREE_YEAR_match2 <- TREE_YEAR %>%
  mutate(tree_id = gsub("[^0-9]", "", tree_id_raw)) %>%
  filter(!tree_id %in% TREE_YEAR_match1$tree_id[!is.na(TREE_YEAR_match1$species)]) %>%
  left_join(INVENTORY %>% select(tree_id, species, species_code, dbh_m, x, y), by = "tree_id")
TREE_YEAR <- bind_rows(
  TREE_YEAR_match1 %>% filter(!is.na(species)),
  TREE_YEAR_match2 %>% filter(!is.na(species)),
  TREE_YEAR_match1 %>% filter(is.na(species)) %>% mutate(tree_id = tree_id_raw)
) %>%
  distinct() %>%
  # untagged stems: species from the label, UNTAG_<position>_<SPECIES>[_dead]
  mutate(.sp_tag = ifelse(grepl("^UNTAG_", tree_id_raw),
                          sub("^UNTAG_[A-Za-z]+_([A-Z]{4}).*$", "\\1", tree_id_raw), NA_character_),
         species_code = coalesce(species_code, .sp_tag),
         species = coalesce(species, unname(species_mapping[.sp_tag]))) %>%
  select(-.sp_tag)

# The monthly survey's tree list (YM_trees_measured) overrides a census match whose
# species disagrees, and fills a missing DBH.
# KNOWN: matched on the Plot Tag spelling; the list's first entry for a label (3355
#   is listed for two trees).
ym_tags <- spellings(F_YM, "Label")
mon_spelling <- spellings(F_MON, "Plot.Tag")
ym <- TREE_YEAR %>% distinct(tree_id_raw) %>%
  mutate(ym_tree = vapply(tree_id_raw, function(s) {
    c <- unique(ym_tags$tree_id[ym_tags$tag_recorded == s])
    if (length(c) > 1) c <- intersect(c, mon_spelling$tree_id[mon_spelling$tag_recorded == s])
    if (length(c) == 1) c else NA_character_
  }, character(1))) %>%
  filter(!is.na(ym_tree))
i_ym <- match(ym$ym_tree, TREES$tree_id)
ym <- ym %>% transmute(tree_id_raw, ym_code = TREES$species_code[i_ym],
                       ym_dbh_m = ifelse(TREES$dbh_source[i_ym] %in% "monthly survey tree list (YM_trees_measured)",
                                         num(TREES$dbh_cm[i_ym]) / 100, NA_real_))
TREE_YEAR <- TREE_YEAR %>% left_join(ym, by = "tree_id_raw") %>%
  mutate(.fix = !is.na(ym_code) & (is.na(species_code) | species_code != ym_code),
         species_code = ifelse(.fix, ym_code, species_code),
         species = ifelse(.fix, unname(species_mapping[ym_code]), species),
         dbh_m = ifelse(.fix | is.na(dbh_m), coalesce(ym_dbh_m, dbh_m), dbh_m)) %>%
  select(-ym_code, -ym_dbh_m, -.fix)

TREE_YEAR <- TREE_YEAR %>%
  group_by(tree_id, Date) %>%
  slice(1) %>%
  ungroup() %>%
  select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
         air_temp_C, soil_temp_C, soil_moisture_abs, chamber_type, dbh_m, x, y) %>%
  filter(!is.na(stem_flux_umol_m2_s))
cat(sprintf("TREE_YEAR: %d deployments from %d trees; soil moisture for %d; %d without a census match\n\n",
            nrow(TREE_YEAR), length(unique(TREE_YEAR$tree_id)), sum(!is.na(TREE_YEAR$soil_moisture_abs)),
            sum(is.na(TREE_YEAR$species) | is.na(TREE_YEAR$x) | is.na(TREE_YEAR$y))))

# =============================================================================
# 9. SOIL_YEAR: monthly soil collars 2020-21
# =============================================================================
# Deployments that failed the closure screen (flux_soil qc_pass) and the two collars
# outside the stand (collars in_stand: WS_1-1, WS_1-2) are left out.
# KNOWN: air temperature is flux_soil's: the tower reading nearest 00:00 of the date
#   (else chamber temperature), not the deployment time.
# KNOWN: soil temperature and VWC are flux_soil's: the collar's own reading that
#   date, else the plot-tag mean that date (soilmoisture_total.csv).
in_stand_collar <- COLLARS$collar_id[COLLARS$in_stand == "TRUE"]
SOIL_YEAR <- FLUX_S %>%
  filter(qc_pass == "TRUE") %>%
  transmute(
    site_id = collar_id,
    Date = as.POSIXct(as.Date(date), tz = "UTC"),
    month = month(Date),
    year = year(Date),
    soil_flux_umol_m2_s = num(CH4_flux_nmol_m2_s),
    air_temp_C = num(air_temp_C),
    soil_temp_C = num(soil_temp_C),
    soil_moisture_abs = num(soil_vwc)) %>%
  filter(!is.na(soil_flux_umol_m2_s))
n_before_stand <- nrow(SOIL_YEAR)
SOIL_YEAR <- SOIL_YEAR %>% filter(site_id %in% in_stand_collar)
SOIL_YEAR <- SOIL_YEAR %>%
  mutate(x = num(COLLARS$longitude[match(site_id, COLLARS$collar_id)]),
         y = num(COLLARS$latitude[match(site_id, COLLARS$collar_id)])) %>%
  select(site_id, Date, month, year, soil_flux_umol_m2_s, air_temp_C,
         soil_temp_C, soil_moisture_abs, x, y)
stopifnot(!anyNA(SOIL_YEAR$x), !anyNA(SOIL_YEAR$y))
cat(sprintf("SOIL_YEAR: %d deployments on %d collars (%d on out-of-stand collars dropped)\n\n",
            nrow(SOIL_YEAR), length(unique(SOIL_YEAR$site_id)), n_before_stand - nrow(SOIL_YEAR)))

# =============================================================================
# 10. DECEMBER 2020 MOISTURE RASTER
# =============================================================================
# Survey points plus the stream-channel flags at saturation (VWC = 1), interpolated
# (akima, linear) onto a 100 x 100 lon/lat grid.
combined_moisture <- rbind(
  data.frame(Longitude = num(MOIST$longitude), Latitude = num(MOIST$latitude), VWC = num(MOIST$soil_vwc_pct) / 100),
  data.frame(Longitude = num(STREAM$longitude), Latitude = num(STREAM$latitude), VWC = 1.0))
lon_range <- range(combined_moisture$Longitude)
lat_range <- range(combined_moisture$Latitude)
MOISTURE_DEC_RASTER <- interp(
  x = combined_moisture$Longitude,
  y = combined_moisture$Latitude,
  z = combined_moisture$VWC,
  xo = seq(lon_range[1], lon_range[2], length = 100),
  yo = seq(lat_range[1], lat_range[2], length = 100),
  duplicate = "mean"
)

# Lookup function for moisture at any x,y coordinate
moisture_lookup_xy <- function(x, y) {
  moisture_values <- numeric(length(x))
  for(i in 1:length(x)) {
    if(!is.na(x[i]) && !is.na(y[i])) {
      lon_idx <- which.min(abs(MOISTURE_DEC_RASTER$x - x[i]))
      lat_idx <- which.min(abs(MOISTURE_DEC_RASTER$y - y[i]))
      moisture_values[i] <- MOISTURE_DEC_RASTER$z[lon_idx, lat_idx]
    } else {
      moisture_values[i] <- NA
    }
  }
  return(moisture_values)
}
cat("December moisture raster:", length(MOISTURE_DEC_RASTER$x), "x", length(MOISTURE_DEC_RASTER$y), "\n")

# =============================================================================
# 11. MONTHLY DRIVERS
# =============================================================================
# KNOWN: air by calendar month of the logger-clock timestamp, all years pooled; soil
#   temperature by month from the training rows themselves (TREE_JULY, TREE_YEAR,
#   SOIL_YEAR).
air_monthly <- weather_clean %>%
  mutate(month = month(TIMESTAMP)) %>%
  group_by(month) %>%
  summarise(
    air_temp_C_mean = mean(air_temp_C, na.rm = TRUE),
    air_temp_C_sd = sd(air_temp_C, na.rm = TRUE),
    .groups = "drop"
  )
soil_monthly <- bind_rows(
  TREE_JULY %>% select(month, soil_temp_C),
  TREE_YEAR %>% select(month, soil_temp_C),
  SOIL_YEAR %>% select(month, soil_temp_C)
) %>%
  filter(!is.na(soil_temp_C)) %>%
  group_by(month) %>%
  summarise(
    soil_temp_C_mean = mean(soil_temp_C, na.rm = TRUE),
    soil_temp_C_sd = sd(soil_temp_C, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(
    soil_temp_C_mean = ifelse(is.nan(soil_temp_C_mean), NA_real_, soil_temp_C_mean),
    soil_temp_C_sd = ifelse(is.nan(soil_temp_C_sd), NA_real_, soil_temp_C_sd)
  )
DRIVERS <- tibble(month = 1:12) %>%
  left_join(air_monthly, by = "month") %>%
  left_join(soil_monthly, by = "month")
cat("Monthly drivers: air", sum(!is.na(DRIVERS$air_temp_C_mean)), "months, soil",
    sum(!is.na(DRIVERS$soil_temp_C_mean)), "months\n")

# =============================================================================
# 12. CHECKS AND SAVE
# =============================================================================
need <- list(TREE_JULY = c("tree_id", "species", "Date", "stem_flux_umol_m2_s", "x", "y"),
             TREE_YEAR = c("tree_id", "species", "Date", "stem_flux_umol_m2_s", "x", "y"),
             SOIL_YEAR = c("site_id", "Date", "soil_flux_umol_m2_s", "x", "y"),
             INVENTORY = c("tree_id", "species", "dbh_m", "x", "y"))
for (n in names(need)) {
  miss <- setdiff(need[[n]], names(get(n)))
  if (nrow(get(n)) == 0 || length(miss)) stop(n, " is empty or lacks: ", paste(miss, collapse = ", "))
}

rf_workflow_data <- list(
  PLACEHOLDER_TREE_JULY = TREE_JULY,
  PLACEHOLDER_TREE_YEAR = TREE_YEAR,
  PLACEHOLDER_SOIL_YEAR = SOIL_YEAR,
  PLACEHOLDER_INVENTORY = INVENTORY,
  PLACEHOLDER_MOISTURE_DEC_RASTER = MOISTURE_DEC_RASTER,
  PLACEHOLDER_MOISTURE_LOOKUP_XY = moisture_lookup_xy,
  PLACEHOLDER_DRIVERS = DRIVERS,
  PLACEHOLDER_TAXONOMY = TAXONOMY,
  PLACEHOLDER_PLOT_AREA = PLOT_AREA,
  geodetic_transform = geodetic_transform,
  geodetic_inverse = geodetic_inverse,
  ref_lat = ref_lat,
  ref_lon = ref_lon,
  rotation_angle = rotation_angle,
  earth_radius = earth_radius,
  metadata = list(
    coordinate_system = "WGS84 geodetic",
    x_variable = "longitude (decimal degrees)",
    y_variable = "latitude (decimal degrees)",
    creation_date = Sys.time(),
    data_years = unique(c(TREE_JULY$year, TREE_YEAR$year, SOIL_YEAR$year)),
    n_species = length(unique(c(TREE_JULY$species, TREE_YEAR$species))),
    script_version = "2024-12 with 2023 data and geodetic coordinates"
  )
)
save(rf_workflow_data, file = "../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")
cat("\nSaved data/processed/integrated/rf_workflow_input_data_with_2023.RData\n")
cat("Completed at:", format(Sys.time()), "\n")
