# ==============================================================================
# Load and Prep Data for Random Forest Workflow
# ==============================================================================
# Purpose: Loads and compiles all input datasets for the Random Forest workflow
#   (merged tree data, flux, soil, inventory, weather, moisture).
#
# Pipeline stage: 3 — Upscaling
# Run after: 00_harmonization/02_harmonize_all_data.R
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#   - methanogen_tree_flux_complete_dataset.csv (from data/processed/flux/)
#   - semirigid flux datasets (from data/processed/flux/)
#   - ForestGEO inventory CSVs (from data/raw/inventory/)
#   - weather and moisture data (from data/raw/field_data/)
#
# Outputs:
#   - rf_workflow_input_data_with_2023.RData (to data/processed/)
# ==============================================================================

library(dplyr)
library(readr)
library(lubridate)
library(akima)
library(readxl)
library(tidyr)

cat("=== RF WORKFLOW DATA IMPORT - WITH 2023 + GEODETIC ===\n\n")
cat("Coordinate System: WGS84 (x=longitude, y=latitude)\n")
cat("Started at:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. DEFINE FILE PATHS
# =============================================================================

paths <- list(
  merged_tree = "../../data/processed/integrated/merged_tree_dataset_final.csv",
  tree_2023 = "../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv",
  # NOTE (2026 revision): use the dataset that includes the recovered untagged/dead-snag
  # monthly trees (+45 fluxes, +7 trees). The revision analyses all use this file.
  semirigid_tree = "../../data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv",
  soil = "../../data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv",
  inventory = "../../data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_2019.csv",
  met_tower = "../../data/raw/weather/ymf_clean_sorted.csv",
  moisture_extra = "../../data/raw/field_data/ipad_data/Cleaned data/soilmoisture_total.csv",
  # Master per-collar soil temperature + VWC campaign (Date, Site, Plot, Subplot).
  # Resolves to the individual collar, unlike soilmoisture_total.csv which only
  # resolves to the plot, and covers 197 of 288 soil flux records vs 104.
  # Compiled by code/revision/rev_compile_soil_env.R from every per-date iPad sheet
  # PLUS the master workbook: 97% coverage of soil flux records (279/288) with
  # genuinely measured per-collar temperature and VWC, vs 39% before.
  soil_env_collar_csv = "../../data/processed/environmental/soil_env_by_collar.csv",
  soil_env_master = "../../data/raw/environmental/Soil_temp_moisture_2020.xlsx",
  moisture_december = "../../data/raw/inventory/spatial_data/soil_moisture_20201216.csv",
  river = "../../data/raw/inventory/spatial_data/River.xlsx",
  plot_locations = "../../data/raw/inventory/spatial_data/plots.csv"
)

# =============================================================================
# 2. SPECIES MAPPING AND TAXONOMY
# =============================================================================

cat("Setting up species mapping and taxonomy...\n")

species_mapping <- c(
  "ACRU" = "Acer rubrum",
  "ACSA" = "Acer saccharum", 
  "ACPE" = "Acer pensylvanicum",
  "BEAL" = "Betula alleghaniensis",
  "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera",
  "BEPO" = "Betula populifolia",
  "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana",
  "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra",
  "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina",
  "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata",
  "CALA" = "Carya laciniosa",
  "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina",
  "SAAL" = "Sassafras albidum"
)

TAXONOMY <- data.frame(
  species_code = names(species_mapping),
  species = unname(species_mapping),
  genus = sub(" .*", "", unname(species_mapping))
) %>%
  mutate(
    family = case_when(
      genus == "Acer" ~ "Sapindaceae",
      genus == "Betula" ~ "Betulaceae", 
      genus == "Fagus" ~ "Fagaceae",
      genus == "Quercus" ~ "Fagaceae",
      genus == "Fraxinus" ~ "Oleaceae",
      genus == "Pinus" ~ "Pinaceae",
      genus == "Tsuga" ~ "Pinaceae",
      genus == "Carya" ~ "Juglandaceae",
      genus == "Prunus" ~ "Rosaceae",
      genus == "Sassafras" ~ "Lauraceae",
      genus == "Kalmia" ~ "Ericaceae",
      TRUE ~ "Unknown"
    ),
    order = case_when(
      family %in% c("Sapindaceae") ~ "Sapindales",
      family %in% c("Betulaceae", "Fagaceae", "Juglandaceae") ~ "Fagales",
      family == "Oleaceae" ~ "Lamiales",
      family == "Pinaceae" ~ "Pinales",
      family == "Rosaceae" ~ "Rosales",
      family == "Lauraceae" ~ "Laurales",
      family == "Ericaceae" ~ "Ericales",
      TRUE ~ "Unknown"
    ),
    class = ifelse(family == "Pinaceae", "Pinopsida", "Magnoliopsida")
  )

cat("✓ Taxonomy built for", nrow(TAXONOMY), "species\n")

# =============================================================================
# 3. LOAD RAW DATA
# =============================================================================

cat("\nLoading raw data files...\n")

# Load all data files with error handling
load_safe <- function(path, read_func, ...) {
  if (file.exists(path)) {
    return(read_func(path, ...))
  } else {
    warning(paste("File not found:", path))
    return(NULL)
  }
}

merged_tree_data <- load_safe(paths$merged_tree, read.csv)
tree_2023_data <- load_safe(paths$tree_2023, read_csv, show_col_types = FALSE)
semirigid_data <- load_safe(paths$semirigid_tree, read_csv, show_col_types = FALSE)
soil_data <- load_safe(paths$soil, read_csv, show_col_types = FALSE)
fg19 <- load_safe(paths$inventory, read.csv)
ymf_met <- load_safe(paths$met_tower, read.csv)
# NOTE (2026 revision): this CSV carries a UTF-8 BOM, which mangles the first column
# name ("Date") unless declared -- previously it silently failed to load at all.
semirigid_moisture <- load_safe(paths$moisture_extra, read.csv, fileEncoding = "UTF-8-BOM")
plot_locations <- load_safe(paths$plot_locations, read.csv)

# =============================================================================
# CHAMBER INTEGRITY SCREEN (2026 revision)
# -----------------------------------------------------------------------------
# A chamber whose CH4 concentration STARTS far above ambient (~2000 ppb) was not
# properly flushed -- the "flux" then records the chamber venting down, not soil or
# stem exchange. This is the one failure mode present in these data: exactly two
# soil deployments (both 2020-06-17) start at 11,687 and 49,049 ppb and produce
# +31.6 and -176.1 nmol m-2 s-1, the latter ~39x any other uptake observed.
#
# This REPLACES the previous statistical filters, which were doing the wrong thing:
#   * MAD k=8 on soil deleted four genuine wetland-margin emissions (11.6-63.1
#     nmol m-2 s-1, all with normal C0) while keeping the -176 artifact, making the
#     soil sink ~2.7x too strong.
#   * A 1st/99th-percentile trim on tree flux removed the largest emissions -- the
#     stem hotspots this study documents -- costing 33% of the mean.
#   * An R2 >= 0.7 gate is not a quality criterion here at all: rho(R2, |flux|) =
#     0.95, so it simply deletes small fluxes (9x smaller in median magnitude) and
#     makes the soil mean 30% more negative. Both artifacts above pass it easily
#     (R2 = 0.94 and 0.9995).
# Unlike those, a C0 threshold is mechanistic and unbiased with respect to flux size.
# =============================================================================
C0_MAX_PPB <- 5000   # ~2.5x ambient

screen_chamber_integrity <- function(df, c0_col, flux_col, label) {
  if (is.null(df) || !(c0_col %in% names(df))) return(df)
  c0 <- suppressWarnings(as.numeric(df[[c0_col]]))
  bad <- !is.na(c0) & c0 > C0_MAX_PPB & !is.na(suppressWarnings(as.numeric(df[[flux_col]])))
  if (any(bad)) {
    cat(sprintf("  %s: dropping %d deployment(s) with C0 > %d ppb (max %.0f)\n",
                label, sum(bad), C0_MAX_PPB, max(c0[bad])))
    df[[flux_col]][bad] <- NA
  } else {
    cat(sprintf("  %s: no deployments exceed C0 = %d ppb\n", label, C0_MAX_PPB))
  }
  df
}

# The revision dataset was written with syntactic column names (Plot.Tag); the older
# file used literal spaces (`Plot Tag`). Normalise so downstream code is agnostic.
if (!is.null(semirigid_data)) {
  nm <- names(semirigid_data)
  ren <- c("Plot.Tag" = "Plot Tag", "Plot.Letter" = "Plot Letter")
  for (from in names(ren)) if (from %in% nm && !(ren[[from]] %in% nm)) {
    names(semirigid_data)[names(semirigid_data) == from] <- ren[[from]]
  }
  cat("  semirigid columns normalised (", ncol(semirigid_data), "cols,", nrow(semirigid_data), "rows )\n")
}

cat("\nChamber integrity screen:\n")
soil_data      <- screen_chamber_integrity(soil_data,      "CH4_C0",   "CH4_best.flux",   "soil monthly")
semirigid_data <- screen_chamber_integrity(semirigid_data, "CH4_C0.x", "CH4_best.flux.x", "tree monthly")
tree_2023_data <- screen_chamber_integrity(tree_2023_data, "CH4_C0",   "CH4_best.flux",   "tree 2023")

# Check critical data loaded
critical_missing <- c()
if (is.null(merged_tree_data)) critical_missing <- c(critical_missing, "merged_tree")
if (is.null(fg19)) critical_missing <- c(critical_missing, "inventory")
if (is.null(ymf_met)) critical_missing <- c(critical_missing, "met_tower")

if (length(critical_missing) > 0) {
  stop(paste("Critical data files missing:", paste(critical_missing, collapse = ", ")))
}

cat("✓ All critical data loaded\n")
if (!is.null(tree_2023_data)) cat("  ✓ 2023 tree data loaded\n") else cat("  ⚠ 2023 tree data not found\n")
if (!is.null(semirigid_data)) cat("  ✓ Semirigid data loaded\n") else cat("  ⚠ Semirigid data not found\n")
if (!is.null(soil_data)) cat("  ✓ Soil data loaded\n") else cat("  ⚠ Soil data not found\n")

# =============================================================================
# 4. GEODETIC TRANSFORMATION SETUP
# =============================================================================

cat("\nSetting up geodetic transformation...\n")

# Reference point (plot center)
ref_lat <- 41.989211
ref_lon <- -72.13092
lat1_rad <- ref_lat * pi/180
lon1_rad <- ref_lon * pi/180
rotation_angle <- 0.174533  # 10 degrees in radians
earth_radius <- 6371000  # meters

# Forward transformation: plot coords (m) to GPS (decimal degrees)
geodetic_transform <- function(px, py, lat1_rad, lon1_rad, rotation_angle, earth_radius) {
  distance <- sqrt(px^2 + py^2)
  original_bearing <- atan2(px, py)
  rotated_bearing <- original_bearing + rotation_angle
  
  px_rotated <- distance * sin(rotated_bearing)
  py_rotated <- distance * cos(rotated_bearing)
  
  bearing <- atan2(px_rotated, py_rotated)
  distance_m <- sqrt(px_rotated^2 + py_rotated^2)
  
  lat2_rad <- asin(sin(lat1_rad) * cos(distance_m/earth_radius) + 
                     cos(lat1_rad) * sin(distance_m/earth_radius) * cos(bearing))
  
  lon2_rad <- lon1_rad + atan2(sin(bearing) * sin(distance_m/earth_radius) * cos(lat1_rad),
                               cos(distance_m/earth_radius) - sin(lat1_rad) * sin(lat2_rad))
  
  return(list(longitude = lon2_rad * 180/pi, latitude = lat2_rad * 180/pi))
}

# Inverse transformation: GPS to plot coords
geodetic_inverse <- function(lat, lon, lat1_rad, lon1_rad, rotation_angle, earth_radius) {
  lat_rad <- lat * pi/180
  lon_rad <- lon * pi/180
  delta_lon <- lon_rad - lon1_rad
  
  y_comp <- sin(delta_lon) * cos(lat_rad)
  x_comp <- cos(lat1_rad) * sin(lat_rad) - sin(lat1_rad) * cos(lat_rad) * cos(delta_lon)
  bearing <- atan2(y_comp, x_comp)
  
  a <- sin((lat_rad - lat1_rad)/2)^2 + cos(lat1_rad) * cos(lat_rad) * sin(delta_lon/2)^2
  c <- 2 * asin(sqrt(a))
  distance_m <- earth_radius * c
  
  bearing_adjusted <- bearing - rotation_angle
  x <- distance_m * sin(bearing_adjusted)
  y <- distance_m * cos(bearing_adjusted)
  
  return(list(x = x, y = y))
}

cat("✓ Geodetic transformation ready\n")
cat("  Reference point:", ref_lat, "N,", ref_lon, "W\n")
cat("  Rotation angle:", rotation_angle * 180/pi, "degrees\n")

# =============================================================================
# 5. PREPARE INVENTORY WITH GEODETIC GPS COORDINATES
# =============================================================================

cat("\nPreparing inventory with geodetic GPS coordinates...\n")

# Filter valid inventory records
fg19_filtered <- fg19 %>%
  filter(!is.na(PX), !is.na(PY), !is.na(DBH))

cat("  Inventory records with valid coordinates:", nrow(fg19_filtered), "\n")

# Transform plot coordinates to GPS
coords_result <- geodetic_transform(fg19_filtered$PX, fg19_filtered$PY, 
                                    lat1_rad, lon1_rad, rotation_angle, earth_radius)

fg19_clean <- fg19_filtered %>%
  mutate(
    x = coords_result$longitude,  # x = longitude
    y = coords_result$latitude,   # y = latitude
    dbh_m = DBH / 100,
    tree_id = as.character(Tag),  # Keep original Tag format (no "tree_" prefix)
    species_code = Species_Code,
    species = species_mapping[Species_Code]
  ) %>%
  filter(!is.na(species)) %>%
  group_by(tree_id) %>%
  slice_head(n = 1) %>%
  ungroup()

# Report unmapped species codes
unmapped_codes <- unique(fg19_filtered$Species_Code[!fg19_filtered$Species_Code %in% names(species_mapping)])
if (length(unmapped_codes) > 0) {
  cat("  ⚠ Unmapped species codes:", paste(unmapped_codes, collapse = ", "), "\n")
}

INVENTORY <- fg19_clean %>%
  select(tree_id, species, species_code, dbh_m, x, y, PX, PY) %>%
  filter(!is.na(x), !is.na(y))

# =============================================================================
# CENSUSED PLOT DOMAIN (2026 revision)
# -----------------------------------------------------------------------------
# PLOT_AREA is the denominator of the entire tree budget (Phi_tree = sum(F_i*S_i)/A_plot),
# so getting it right matters more than anything else in the scaling.
#
# It was previously the NSEW bounding box of stem lat/lon = 10.23 ha. That is wrong three
# times over: (a) the plot is rotated 10 deg from north, so an axis-aligned box is the wrong
# SHAPE; (b) the census does not actually cover the nominal 350 x 300 m plot -- stems are
# confined to roughly PX 0-200, PY 0-200; and (c) the box was inflated to 10.23 ha by just
# 5 stray stems (0.07%) lying outside the surveyed area, almost certainly coordinate errors.
#
# Evidence for the censused extent: on a 50 m grid every cell with PX<200 & PY<200 holds
# 86-2268 stems, while every cell beyond holds 0-1. Stem density confirms it -- 685 stems/ha
# over the old box is implausibly low for a census whose median DBH is 2.3 cm, whereas
# 1,752 stems/ha over the censused core is right.
#
# Working in the plot's own (PX, PY) metre coordinates also sidesteps the rotation entirely.
# NOTE: this cancels in the whole-woody-surface scenario -- Phi_tree and the stem area index
# both carry 1/A_plot -- but it scales the MEASURED 0-2 m budget directly.
# =============================================================================
CENSUS_PX_MAX <- 200   # m, plot coordinates
CENSUS_PY_MAX <- 200

n_before_core <- nrow(INVENTORY)
INVENTORY <- INVENTORY %>%
  filter(PX >= 0, PX <= CENSUS_PX_MAX, PY >= 0, PY <= CENSUS_PY_MAX)
cat(sprintf("  Censused core: kept %d of %d stems (dropped %d outside PX/PY <= %d m)\n",
            nrow(INVENTORY), n_before_core, n_before_core - nrow(INVENTORY), CENSUS_PX_MAX))

PLOT_AREA <- CENSUS_PX_MAX * CENSUS_PY_MAX

coord_ranges <- list(
  x = range(INVENTORY$x, na.rm = TRUE),
  y = range(INVENTORY$y, na.rm = TRUE)
)

cat("✓ Inventory:", nrow(INVENTORY), "trees with GPS coordinates\n")
cat("  Coordinate ranges: lon [", round(coord_ranges$x[1], 6), "-", round(coord_ranges$x[2], 6), 
    "], lat [", round(coord_ranges$y[1], 6), "-", round(coord_ranges$y[2], 6), "]\n")
cat("  Plot area:", round(PLOT_AREA, 0), "m² (", round(PLOT_AREA/10000, 2), "ha)\n\n")

# =============================================================================
# 6. PREPARE MET TOWER DATA
# =============================================================================

cat("Preparing met tower data...\n")

weather_clean <- ymf_met %>%
  mutate(
    TIMESTAMP = as.POSIXct(TIMESTAMP, tz = "UTC"),
    air_temp_C = Tair_Avg
  ) %>%
  select(TIMESTAMP, air_temp_C) %>%
  filter(!is.na(air_temp_C)) %>%
  arrange(TIMESTAMP)

# Helper function to add met tower temperature
add_met_tower_temp <- function(df, weather_df, time_col = "Date") {
  df$datetime_match <- as.POSIXct(df[[time_col]], tz = "UTC")
  df$air_temp_C_tower <- NA_real_
  
  for(i in 1:nrow(df)) {
    flux_time <- df$datetime_match[i]
    if(!is.na(flux_time)) {
      time_diffs <- abs(as.numeric(weather_df$TIMESTAMP - flux_time))
      if(length(time_diffs) > 0 && min(time_diffs) <= 86400) {  # Within 24 hours
        df$air_temp_C_tower[i] <- weather_df$air_temp_C[which.min(time_diffs)]
      }
    }
  }
  df$datetime_match <- NULL
  return(df)
}

cat("✓ Weather data ready:", nrow(weather_clean), "records\n")
cat("  Date range:", format(min(weather_clean$TIMESTAMP), "%Y-%m-%d"), 
    "to", format(max(weather_clean$TIMESTAMP), "%Y-%m-%d"), "\n\n")

# =============================================================================
# 7. PROCESS 2021 TREE DATA (FROM MERGED FILE)
# =============================================================================

cat("Processing 2021 summer intensive tree data...\n")

# Clean tree IDs
merged_tree_data$tree_id <- gsub(" ", "", merged_tree_data$tree_id)

# Diagnostic: Check tree ID formats
cat("  Sample 2021 tree IDs (first 10):", paste(head(unique(merged_tree_data$tree_id), 10), collapse=", "), "\n")
cat("  Sample inventory tree IDs (first 10):", paste(head(unique(INVENTORY$tree_id), 10), collapse=", "), "\n")

base_date <- as.POSIXct("2021-07-28", tz = "UTC")

# Process trees with 125cm measurements
TREE_JULY_125 <- merged_tree_data %>%
  filter(!is.na(CH4_best.flux_125cm)) %>%
  mutate(
    Date = base_date,
    tree_id = as.character(tree_id),
    species_code = species_id,
    species = ifelse(species_id %in% names(species_mapping), 
                     species_mapping[species_id], 
                     paste("Unknown", species_id)),
    stem_flux_umol_m2_s = CH4_best.flux_125cm,
    air_temp_C = Temp_Air_125cm,
    soil_temp_C = SoilTemp_mean,
    soil_moisture_abs = VWC_mean / 100,
    dbh_m = case_when(
      is.na(dbh) ~ NA_real_,
      dbh > 3 ~ dbh / 100,  # Convert cm to m if needed
      TRUE ~ dbh
    ),
    chamber_type = "rigid",
    measurement_height_cm = 125,
    month = month(Date),
    year = year(Date)
  ) %>%
  select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s, 
         air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m, chamber_type,
         measurement_height_cm)

# =============================================================================
# 2021 HEIGHT CAMPAIGN: 50 cm and 200 cm (2026 revision)
# -----------------------------------------------------------------------------
# Only the 125 cm measurements were previously used, so 291 of the campaign's 439
# fluxes never reached the model. Including all three heights is not just more data:
# the budget integrates flux over the 0-2 m stem band, so a mean ACROSS heights
# estimates that band better than assuming 125 cm is representative of it.
# Height is carried through as a column so it can be inspected or used as a feature.
# =============================================================================
make_height_block <- function(df, flux_col, temp_col, height_cm) {
  if (!(flux_col %in% names(df))) return(NULL)
  df %>%
    filter(!is.na(.data[[flux_col]])) %>%
    mutate(
      Date = base_date,
      tree_id = as.character(tree_id),
      species_code = species_id,
      species = ifelse(species_id %in% names(species_mapping),
                       species_mapping[species_id], paste("Unknown", species_id)),
      stem_flux_umol_m2_s = .data[[flux_col]],
      air_temp_C = if (temp_col %in% names(df)) .data[[temp_col]] else Temp_Air_125cm,
      soil_temp_C = SoilTemp_mean,
      soil_moisture_abs = VWC_mean / 100,
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

TREE_JULY_50  <- make_height_block(merged_tree_data, "CH4_best.flux_50cm",  "Temp_Air_50cm",  50)
TREE_JULY_200 <- make_height_block(merged_tree_data, "CH4_best.flux_200cm", "Temp_Air_200cm", 200)
cat("  50cm chambers:",  if (is.null(TREE_JULY_50)) 0 else nrow(TREE_JULY_50), "\n")
cat("  200cm chambers:", if (is.null(TREE_JULY_200)) 0 else nrow(TREE_JULY_200), "\n")

# Process Kalmia measurements
KALMIA_DATA <- merged_tree_data %>%
  filter(species_id == "KALA" | grepl("(?i)kalmia", species_id)) %>%
  filter(!is.na(CH4_best.flux_75cm) | !is.na(CH4_best.flux_RootCrowncm)) %>%
  mutate(
    Date = base_date,
    tree_id = as.character(tree_id),
    stem_flux_umol_m2_s = coalesce(CH4_best.flux_75cm, CH4_best.flux_RootCrowncm),
    air_temp_C = coalesce(Temp_Air_75cm, Temp_Air_RootCrowncm),
    soil_temp_C = SoilTemp_mean,
    soil_moisture_abs = VWC_mean / 100,
    dbh_m = ifelse(!is.na(dbh) & dbh > 0, dbh / 100, NA_real_),
    species_code = "KALA",
    species = "Kalmia latifolia",
    chamber_type = "rigid",
    month = month(Date),
    year = year(Date)
  ) %>%
  select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
         air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m, chamber_type)

# Process other 75cm measurements
OTHER_75CM <- merged_tree_data %>%
  filter(is.na(CH4_best.flux_125cm), !is.na(CH4_best.flux_75cm),
         !(species_id %in% c("KALA", "Kalmia"))) %>%
  mutate(
    Date = base_date,
    tree_id = as.character(tree_id),
    species_code = species_id,
    species = ifelse(species_id %in% names(species_mapping), 
                     species_mapping[species_id], 
                     paste("Unknown", species_id)),
    stem_flux_umol_m2_s = CH4_best.flux_75cm,
    air_temp_C = Temp_Air_75cm,
    soil_temp_C = SoilTemp_mean,
    soil_moisture_abs = VWC_mean / 100,
    dbh_m = case_when(
      is.na(dbh) ~ NA_real_,
      dbh > 3 ~ dbh / 100,
      TRUE ~ dbh
    ),
    chamber_type = "rigid",
    month = month(Date),
    year = year(Date)
  ) %>%
  select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
         air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m, chamber_type)

# Combine all 2021 July data
TREE_JULY_2021 <- bind_rows(TREE_JULY_125, TREE_JULY_50, TREE_JULY_200, KALMIA_DATA, OTHER_75CM)

cat("✓ 2021 summer data:", nrow(TREE_JULY_2021), "measurements\n")
cat("  125cm chambers:", nrow(TREE_JULY_125), "\n")
cat("  Kalmia:", nrow(KALMIA_DATA), "\n")
cat("  Other 75cm:", nrow(OTHER_75CM), "\n\n")

# =============================================================================
# 8. PROCESS 2023 TREE DATA
# =============================================================================

if (!is.null(tree_2023_data)) {
  cat("Processing 2023 tree data...\n")
  
  # Identify relevant columns
  tree_tag_col <- names(tree_2023_data)[grepl("tree.*tag", names(tree_2023_data), ignore.case = TRUE)][1]
  species_col <- names(tree_2023_data)[grepl("species.*code", names(tree_2023_data), ignore.case = TRUE)][1]
  
  if (!is.na(tree_tag_col) && !is.na(species_col)) {
    tree_2023_clean <- tree_2023_data %>%
      mutate(
        tree_tag_num = suppressWarnings(as.numeric(.data[[tree_tag_col]])),
        species_code = .data[[species_col]],
        date_clean = as.POSIXct(Date, format = "%m/%d/%y", tz = "UTC"),
        month = month(date_clean),
        year = year(date_clean),
        stem_flux_umol_m2_s = CH4_best.flux
      ) %>%
      filter(!is.na(stem_flux_umol_m2_s), !is.na(tree_tag_num))
    
    # Match with inventory
    # =========================================================================
    # NOTE (2026 revision): this block previously required every 2023 tree to match an
    # INVENTORY record, because it took x, y AND dbh from the inventory and then filtered
    # on !is.na(x). That discarded 93 measurements on trees that were never censused --
    # even though the 2023 file carries its OWN DBH (338/338), VWC (~330), soil temp
    # (338) and air temp (337). It also threw away those measured covariates, setting
    # soil_moisture_abs to NA outright.
    #
    # Training data does not need to be inventory stems: the model is FIT on measured
    # trees and APPLIED to the inventory. (The 2021 height campaign makes this obvious --
    # only 1 of its 235 trees is a censused stem, yet it trains the model fine.)
    # Inventory coordinates are still joined where available, since they let a row
    # contribute to the moisture affine calibration, but they are no longer required.
    # =========================================================================
    # Resolve the campaign's own covariate columns by pattern -- their literal names
    # contain spaces, parentheses and a degree symbol, so they cannot be referenced directly.
    .pick <- function(df, pat) { h <- grep(pat, names(df), value = TRUE, ignore.case = TRUE)[1]
                                 if (is.na(h)) NA_real_ else suppressWarnings(as.numeric(df[[h]])) }
    tree_2023_clean$dbh_own_m         <- .pick(tree_2023_clean, "^DBH") / 100
    tree_2023_clean$soil_moisture_own <- .pick(tree_2023_clean, "^vwc_mean$") / 100
    tree_2023_clean$soil_temp_own     <- .pick(tree_2023_clean, "^Soil.?Temp")
    tree_2023_clean$air_temp_own      <- .pick(tree_2023_clean, "^air_temp_C$")

    tree_2023_matched <- tree_2023_clean %>%
      mutate(tree_id = as.character(tree_tag_num)) %>%
      left_join(INVENTORY %>% select(tree_id, x, y, dbh_m), by = "tree_id") %>%
      mutate(
        species = species_mapping[species_code],
        Date = date_clean,
        chamber_type = "rigid",
        # prefer the inventory DBH (censused, consistent) but fall back to the
        # campaign's own field measurement rather than dropping the tree
        dbh_m = dplyr::coalesce(dbh_m, dbh_own_m)
      ) %>%
      filter(!is.na(species), !is.na(dbh_m))

    # Add met tower temperature (fallback where the field reading is missing)
    tree_2023_matched <- add_met_tower_temp(tree_2023_matched, weather_clean, "Date")

    TREE_JULY_2023 <- tree_2023_matched %>%
      mutate(
        air_temp_C        = dplyr::coalesce(air_temp_own, air_temp_C_tower),
        soil_temp_C       = soil_temp_own,
        soil_moisture_abs = soil_moisture_own
      ) %>%
      select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
             air_temp_C, soil_temp_C, soil_moisture_abs,
             dbh_m, chamber_type, x, y)
    
    # Check for unmatched trees
    unmatched_2023 <- tree_2023_clean %>%
      mutate(tree_id = as.character(tree_tag_num)) %>%  # Same format as above
      anti_join(INVENTORY, by = "tree_id")
    
    cat("✓ 2023 data:", nrow(TREE_JULY_2023), "measurements matched to inventory\n")
    if (nrow(unmatched_2023) > 0) {
      cat("  ⚠", nrow(unmatched_2023), "measurements unmatched (no inventory record)\n")
    }
    cat("\n")
  } else {
    cat("  ⚠ Could not identify required columns in 2023 data\n")
    TREE_JULY_2023 <- NULL
  }
} else {
  TREE_JULY_2023 <- NULL
  cat("⚠ 2023 tree data not available\n\n")
}

# =============================================================================
# 9. ADD GPS COORDINATES AND MERGE 2021 + 2023
# =============================================================================

cat("Adding GPS coordinates and merging datasets...\n")

# Add coordinates to 2021 data
TREE_JULY_2021 <- TREE_JULY_2021 %>%
  left_join(INVENTORY %>% select(tree_id, x, y), by = "tree_id")

# Combine all July data
if (!is.null(TREE_JULY_2023)) {
  TREE_JULY <- bind_rows(TREE_JULY_2021, TREE_JULY_2023) %>%
    arrange(Date)
} else {
  TREE_JULY <- TREE_JULY_2021
}

# Report summary
n_2021 <- sum(TREE_JULY$year == 2021, na.rm=TRUE)
n_2023 <- sum(TREE_JULY$year == 2023, na.rm=TRUE)
n_with_coords <- sum(!is.na(TREE_JULY$x) & !is.na(TREE_JULY$y))

cat("✓ TREE_JULY:", nrow(TREE_JULY), "measurements from",
    length(unique(TREE_JULY$tree_id)), "trees\n")
cat("  2021:", n_2021, "measurements\n")
if (n_2023 > 0) cat("  2023:", n_2023, "measurements\n")
cat("  With coordinates:", n_with_coords, "of", nrow(TREE_JULY), "\n\n")

# =============================================================================
# 10. PREPARE MOISTURE CALIBRATION DATA
# =============================================================================

cat("Preparing moisture calibration data...\n")

# =============================================================================
# PER-COLLAR SOIL TEMPERATURE AND MOISTURE (2026 revision)
# -----------------------------------------------------------------------------
# The pipeline previously drew these only from soilmoisture_total.csv, which keys on
# PLOT (1 or 2) rather than collar, so every collar in a plot received the same value
# and soil temperature ended up with exactly ONE distinct value per month across the
# whole dataset -- no spatial information at all.
#
# Soil_temp_moisture_2020.xlsx records Site + Plot + Subplot, i.e. the individual
# collar, on 19 dates. Matched on date + landscape position + collar it covers 197 of
# 288 soil flux records (68%) with genuinely measured temperature and VWC, against 104
# (39%) before. The remaining 91 fall on four dates when the environmental campaign did
# not run; those keep the monthly-mean fallback.
# =============================================================================
soil_env_collar <- NULL
if (!is.null(paths$soil_env_collar_csv) && file.exists(paths$soil_env_collar_csv)) {
  soil_env_collar <- read.csv(paths$soil_env_collar_csv, stringsAsFactors = FALSE) %>%
    mutate(Date = as.Date(Date),
           soil_moisture_abs = suppressWarnings(as.numeric(soil_moisture_pct)) / 100) %>%
    filter(!is.na(Date), !is.na(plot_letter)) %>%
    select(Date, plot_letter, plot_tag, soil_temp_C, soil_moisture_abs)
  cat("\u2713 Per-collar soil env (compiled):", nrow(soil_env_collar), "records on",
      length(unique(soil_env_collar$Date)), "dates\n")
} else if (!is.null(paths$soil_env_master) && file.exists(paths$soil_env_master)) {
  suppressMessages(library(readxl))
  .env <- read_excel(paths$soil_env_master, sheet = 1)
  names(.env) <- make.names(names(.env))
  .site_map <- c("Upland" = "U", "Intermediate" = "I",
                 "Wetland dry" = "WD", "Wetland swamp" = "WS")
  soil_env_collar <- .env %>%
    mutate(
      Date        = as.Date(Date),
      plot_letter = unname(.site_map[trimws(Site)]),
      plot_tag    = paste0(Plot, "-", Subplot),
      soil_temp_C = suppressWarnings(as.numeric(Soil.temp)),
      soil_moisture_abs = suppressWarnings(as.numeric(VWC)) / 100
    ) %>%
    filter(!is.na(Date), !is.na(plot_letter), !is.na(soil_temp_C)) %>%
    select(Date, plot_letter, plot_tag, soil_temp_C, soil_moisture_abs)
  cat("✓ Per-collar soil env data:", nrow(soil_env_collar), "records on",
      length(unique(soil_env_collar$Date)), "dates\n")
}

if (!is.null(semirigid_moisture)) {
  moisture_clean <- semirigid_moisture %>%
    mutate(
      Date = suppressWarnings(mdy(Date)),
      site_id = as.character(Plot),
      plot_tag = as.character(Plot),
      plot_letter = toupper(trimws(Plot.letter)),
      soil_temp_C = as.numeric(Soil.temp),
      soil_moisture_abs = as.numeric(VWC) / 100
    ) %>%
    filter(!is.na(Date), !is.na(soil_moisture_abs)) %>%
    select(Date, site_id, plot_tag, plot_letter, soil_temp_C, soil_moisture_abs)
  
  cat("✓ Moisture calibration data ready:", nrow(moisture_clean), "records\n\n")
} else {
  moisture_clean <- data.frame()
  cat("⚠ Moisture calibration data not available\n\n")
}

# =============================================================================
# 11. PREPARE MONTHLY TREE DATA (SEMIRIGID)
# =============================================================================

ENV_MATCH_DAYS <- 3   # sensitivity: 1 day gives ~76% coverage, 3 days ~89%

if (!is.null(semirigid_data)) {
  cat("Preparing semirigid tree data...\n")
  
  # Diagnostic: Check semirigid tree ID formats
  cat("  Sample semirigid Plot Tags (first 20):", paste(head(unique(semirigid_data$`Plot Tag`), 20), collapse=", "), "\n")
  
  TREE_YEAR <- semirigid_data %>%
    mutate(
      Date = as.POSIXct(start.time, tz = "UTC"),
      tree_id_raw = as.character(`Plot Tag`),
      plot_letter = toupper(trimws(`Plot Letter`)),
      stem_flux_umol_m2_s = CH4_best.flux.x,
      chamber_type = "semirigid",
      month = month(Date),
      year = year(Date)
    )
  
  # Try multiple ID matching strategies
  # Strategy 1: Direct match
  TREE_YEAR_match1 <- TREE_YEAR %>%
    mutate(tree_id = tree_id_raw) %>%
    left_join(INVENTORY %>% select(tree_id, species, species_code, dbh_m, x, y), by = "tree_id")
  
  # Strategy 2: Remove any non-numeric characters and try again
  TREE_YEAR_match2 <- TREE_YEAR %>%
    mutate(tree_id = gsub("[^0-9]", "", tree_id_raw)) %>%
    filter(!tree_id %in% TREE_YEAR_match1$tree_id[!is.na(TREE_YEAR_match1$species)]) %>%
    left_join(INVENTORY %>% select(tree_id, species, species_code, dbh_m, x, y), by = "tree_id")
  
  # Combine successful matches
  TREE_YEAR <- bind_rows(
    TREE_YEAR_match1 %>% filter(!is.na(species)),
    TREE_YEAR_match2 %>% filter(!is.na(species)),
    TREE_YEAR_match1 %>% filter(is.na(species)) %>% mutate(tree_id = tree_id_raw)
  ) %>%
    distinct()
  
  # Add met tower temperature
  TREE_YEAR <- add_met_tower_temp(TREE_YEAR, weather_clean, "Date")
  
  # IMPROVED MOISTURE MATCHING
  # Add moisture data with temporal proximity matching
  if (nrow(moisture_clean) > 0) {
    cat("  Matching moisture data with temporal proximity (±7 days)...\n")
    
    TREE_YEAR <- TREE_YEAR %>%
      mutate(
        # Create date columns for matching
        Date_only = as.Date(Date),
        # Preserve any existing values
        air_temp_C_original = coalesce(air_temp_C_tower, Tcham),
        soil_temp_C_original = if("soil_temp_C" %in% names(.)) soil_temp_C else NA_real_,
        soil_moisture_abs_original = if("soil_moisture_abs" %in% names(.)) soil_moisture_abs else NA_real_
      ) %>%
      # Join with moisture data allowing many-to-many
      # NOTE (2026 revision): joins the COMPILED per-collar env archive (291 records,
      # 22 dates) rather than soilmoisture_total.csv (103 records, plot-level, ends
      # Feb 2021). Trees and soil collars were visited on adjacent days of the same
      # campaign, so the closest-match logic below recovers most of them; monthly-tree
      # soil temp/moisture coverage rises from 32% to ~76% at +/-1 day.
      left_join(
        (if (!is.null(soil_env_collar) && nrow(soil_env_collar) > 0) {
            soil_env_collar %>%
              group_by(Date, plot_letter) %>%
              summarise(soil_temp_C = mean(soil_temp_C, na.rm = TRUE),
                        soil_moisture_abs = mean(soil_moisture_abs, na.rm = TRUE),
                        .groups = "drop") %>%
              mutate(across(c(soil_temp_C, soil_moisture_abs), ~ifelse(is.nan(.x), NA_real_, .x)))
          } else moisture_clean) %>%
          mutate(moisture_date = as.Date(Date)) %>%
          select(moisture_date, plot_letter, soil_temp_C_moisture = soil_temp_C, 
                 soil_moisture_abs_moisture = soil_moisture_abs),
        by = "plot_letter",
        relationship = "many-to-many"
      ) %>%
      mutate(
        # Calculate date difference
        date_diff = abs(as.numeric(Date_only - moisture_date)),
        # Same campaign visit. Trees and collars were measured 0-3 days apart; a
        # wider window would borrow conditions from a different visit.
        soil_temp_C_matched = ifelse(date_diff <= ENV_MATCH_DAYS, soil_temp_C_moisture, NA),
        soil_moisture_abs_matched = ifelse(date_diff <= ENV_MATCH_DAYS, soil_moisture_abs_moisture, NA)
      ) %>%
      # Keep only the closest match per tree-date combination
      group_by(tree_id, Date) %>%
      arrange(date_diff) %>%
      slice(1) %>%
      ungroup() %>%
      # Final column selection and naming
      mutate(
        air_temp_C = air_temp_C_original,
        soil_temp_C = coalesce(soil_temp_C_matched, soil_temp_C_original),
        soil_moisture_abs = coalesce(soil_moisture_abs_matched, soil_moisture_abs_original)
      ) %>%
      select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s, 
             air_temp_C, soil_temp_C, soil_moisture_abs, chamber_type, dbh_m, x, y)
    
    # Report matching success
    n_with_moisture <- sum(!is.na(TREE_YEAR$soil_moisture_abs))
    cat("    Matched moisture for", n_with_moisture, "of", nrow(TREE_YEAR), "measurements\n")
  } else {
    # No moisture data available
    TREE_YEAR <- TREE_YEAR %>%
      mutate(
        air_temp_C = coalesce(air_temp_C_tower, Tcham),
        soil_temp_C = NA_real_,
        soil_moisture_abs = NA_real_
      ) %>%
      select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s, 
             air_temp_C, soil_temp_C, soil_moisture_abs, chamber_type, dbh_m, x, y)
  }
  
  # Final filter for valid flux measurements
  TREE_YEAR <- TREE_YEAR %>%
    filter(!is.na(stem_flux_umol_m2_s))
  
  # Check for unmatched semirigid trees
  unmatched_semirigid <- TREE_YEAR %>%
    filter(is.na(species) | is.na(x) | is.na(y))
  
  cat("✓ Semirigid data:", nrow(TREE_YEAR), "measurements\n")
  if (nrow(unmatched_semirigid) > 0) {
    cat("  ⚠", nrow(unmatched_semirigid), "measurements without inventory match\n")
    cat("    Note: Semirigid chambers may use different tree labeling system\n")
    cat("    For RF workflow, observed moisture values will be used for training\n")
  }
  cat("\n")
  
} else {
  TREE_YEAR <- data.frame()
  cat("⚠ Semirigid data not available\n\n")
}

# =============================================================================
# 12. PREPARE SOIL DATA
# =============================================================================

# NOTE (2026 revision): these paths were stale after the data reorganisation, so
# plot_locations silently failed to load and SOIL_YEAR came out EMPTY -- the saved
# RData predated the breakage and hid it. Fail loudly instead.
if (is.null(plot_locations)) stop("plot_locations failed to load - check paths$plot_locations")
if (is.null(soil_data))      stop("soil_data failed to load - check paths$soil")

if (!is.null(soil_data) && !is.null(plot_locations)) {
  cat("Preparing soil data...\n")
  
  SOIL_YEAR <- soil_data %>%
    mutate(
      Date = as.POSIXct(Date, tz = "UTC"),
      plot_tag = `Plot Tag`,
      plot_letter = toupper(trimws(`Plot letter`)),
      site_id = paste0(plot_letter, "_", plot_tag),
      soil_flux_umol_m2_s = CH4_best.flux,
      month = month(Date),
      year = year(Date)
    ) %>%
    add_met_tower_temp(weather_clean, "Date") %>%
    mutate(air_temp_C = coalesce(air_temp_C_tower, Tcham)) %>%
    select(site_id, Date, month, year, plot_letter, plot_tag, soil_flux_umol_m2_s, air_temp_C) %>%
    filter(!is.na(soil_flux_umol_m2_s))

  # =========================================================================
  # OUT-OF-STAND COLLARS (2026 revision)
  # -------------------------------------------------------------------------
  # The soil campaign spans a designed moisture gradient (U -> I -> WD -> WS). Two
  # WS collars sit in OPEN WETLAND, not forest: WS-1-1 and WS-1-2 are the only
  # collars with zero censused stems within 10 m (nearest tree 17.0 and 16.1 m;
  # every other collar has 11-502). They are a different ecosystem, and they
  # produce ALL five soil fluxes above 5 nmol m-2 s-1 (mean +3.9 and +12.6, max
  # +63.1, vs <= +0.4 for every other collar).
  #
  # They are excluded from the STAND budget because the budget upscales to the
  # censused forest, which is also what the tree side upscales to -- like footprints
  # on both sides. This is a spatial criterion, not a statistical one: the other wet
  # collars (WS-1-3, WS-2-1/2/3, all WD) sit inside the forest and are KEPT, since
  # the plot genuinely contains wet patches we want the model to capture.
  #
  # This replaces the old MAD k=8 filter, which deleted four genuine wetland-margin
  # emissions while retaining the one real artifact, making the sink ~2.7x too strong.
  # =========================================================================
  OUT_OF_STAND_COLLARS <- c("WS_1-1", "WS_1-2")
  n_before_stand <- nrow(SOIL_YEAR)
  SOIL_YEAR <- SOIL_YEAR %>% filter(!(site_id %in% OUT_OF_STAND_COLLARS))
  cat(sprintf("  Out-of-stand collars removed: %s (%d of %d measurements)\n",
              paste(OUT_OF_STAND_COLLARS, collapse = ", "),
              n_before_stand - nrow(SOIL_YEAR), n_before_stand))

  # Per-collar measurements first (date + landscape position + collar), so each collar
  # gets ITS OWN temperature and moisture rather than a plot- or month-level constant.
  if (!is.null(soil_env_collar) && nrow(soil_env_collar) > 0) {
    SOIL_YEAR <- SOIL_YEAR %>%
      mutate(Date_only = as.Date(Date)) %>%
      left_join(soil_env_collar %>%
                  rename(soil_temp_collar = soil_temp_C,
                         soil_moisture_collar = soil_moisture_abs),
                by = c("Date_only" = "Date", "plot_letter", "plot_tag"))
    cat(sprintf("  Per-collar soil env matched: %d of %d measurements (%.0f%%)\n",
                sum(!is.na(SOIL_YEAR$soil_temp_collar)), nrow(SOIL_YEAR),
                100 * mean(!is.na(SOIL_YEAR$soil_temp_collar))))
  } else {
    SOIL_YEAR$soil_temp_collar <- NA_real_
    SOIL_YEAR$soil_moisture_collar <- NA_real_
  }

  # Add moisture data if available
  if (nrow(moisture_clean) > 0) {
    SOIL_YEAR <- SOIL_YEAR %>%
      mutate(Date_only = as.Date(Date)) %>%
      left_join(
        moisture_clean %>%
          mutate(Date_only = as.Date(Date), plot_tag_match = plot_tag) %>%
          group_by(Date_only, plot_tag_match) %>%
          summarise(
            soil_temp_C = mean(soil_temp_C, na.rm = TRUE),
            soil_moisture_abs = mean(soil_moisture_abs, na.rm = TRUE),
            .groups = "drop"
          ),
        by = c("Date_only", "plot_tag" = "plot_tag_match")
      ) %>%
      mutate(
        # collar-level measurement takes precedence over the plot-level average
        soil_temp_C = dplyr::coalesce(soil_temp_collar, soil_temp_C),
        soil_moisture_abs = dplyr::coalesce(soil_moisture_collar, soil_moisture_abs),
        soil_temp_C = ifelse(is.nan(soil_temp_C), NA_real_, soil_temp_C),
        soil_moisture_abs = ifelse(is.nan(soil_moisture_abs), NA_real_, soil_moisture_abs)
      ) %>%
      select(-Date_only)
  } else {
    SOIL_YEAR$soil_temp_C <- NA_real_
    SOIL_YEAR$soil_moisture_abs <- NA_real_
  }
  
  # Add plot coordinates
  SOIL_YEAR <- SOIL_YEAR %>%
    mutate(plot_num = as.integer(sub("-.*", "", plot_tag)))
  
  # =========================================================================
  # PER-COLLAR COORDINATES (2026 revision)
  # -------------------------------------------------------------------------
  # plots.csv holds a row per collar (Label = e.g. "WS-1-2") with its own lat/lon.
  # These were being averaged to a single point per PLOT, so all 30 collars collapsed
  # onto just two coordinate pairs. The moisture surface is sampled at (x, y), so
  # every collar in a plot received an identical moisture value and the model could
  # not distinguish an upland collar from a wetland one spatially at all.
  # Match on the collar label first; the plot/letter averages remain as fallbacks.
  # =========================================================================
  collar_coords <- plot_locations %>%
    mutate(collar_key = toupper(gsub("[^A-Za-z0-9]", "", Label))) %>%
    group_by(collar_key) %>%
    summarise(collar_lat = mean(Latitude, na.rm = TRUE),
              collar_lon = mean(Longitude, na.rm = TRUE), .groups = "drop")

  SOIL_YEAR <- SOIL_YEAR %>%
    mutate(collar_key = toupper(gsub("[^A-Za-z0-9]", "", site_id))) %>%
    left_join(collar_coords, by = "collar_key")
  cat(sprintf("  Per-collar coordinates matched: %d of %d measurements (%d distinct positions)\n",
              sum(!is.na(SOIL_YEAR$collar_lat)), nrow(SOIL_YEAR),
              length(unique(paste(SOIL_YEAR$collar_lat, SOIL_YEAR$collar_lon)))))

  # Get plot-level coordinates
  plot_coords_soil <- plot_locations %>%
    group_by(Plot) %>%
    summarise(
      plot_lat = mean(Latitude, na.rm = TRUE),
      plot_lon = mean(Longitude, na.rm = TRUE),
      .groups = "drop"
    )
  
  SOIL_YEAR <- SOIL_YEAR %>%
    left_join(plot_coords_soil, by = c("plot_num" = "Plot"))
  
  # Get letter-level coordinates as fallback
  plot_coords_letter <- plot_locations %>%
    mutate(plot_letter = toupper(trimws(Site))) %>%
    group_by(plot_letter) %>%
    summarise(
      letter_lat = mean(Latitude, na.rm = TRUE),
      letter_lon = mean(Longitude, na.rm = TRUE),
      .groups = "drop"
    )
  
  SOIL_YEAR <- SOIL_YEAR %>%
    left_join(plot_coords_letter, by = "plot_letter") %>%
    mutate(
      final_lat = coalesce(collar_lat, plot_lat, letter_lat),
      final_lon = coalesce(collar_lon, plot_lon, letter_lon),
      x = final_lon,  # x = longitude
      y = final_lat   # y = latitude
    )
  
  # Final selection
  SOIL_YEAR <- SOIL_YEAR %>%
    select(site_id, Date, month, year, soil_flux_umol_m2_s, air_temp_C, 
           soil_temp_C, soil_moisture_abs, x, y)
  
  n_with_coords_soil <- sum(!is.na(SOIL_YEAR$x) & !is.na(SOIL_YEAR$y))
  cat("✓ Soil data prepared:", nrow(SOIL_YEAR), "measurements from",
      length(unique(SOIL_YEAR$site_id)), "sites\n")
  cat("  With coordinates:", n_with_coords_soil, "of", nrow(SOIL_YEAR), "\n\n")
  
} else {
  SOIL_YEAR <- data.frame()
  cat("⚠ Soil data or plot locations not available\n\n")
}

# =============================================================================
# 13. CREATE DECEMBER MOISTURE RASTER
# =============================================================================

cat("Creating December moisture raster...\n")

if(file.exists(paths$moisture_december) && file.exists(paths$river)) {
  moisture_dec <- read.csv(paths$moisture_december)
  river_data <- read_excel(paths$river)
  
  moisture_dec_clean <- moisture_dec %>%
    filter(!is.na(VWC)) %>%
    mutate(VWC_frac = VWC / 100)
  
  # Combine moisture and river data
  combined_moisture <- rbind(
    data.frame(
      Longitude = moisture_dec_clean$Longitude,
      Latitude = moisture_dec_clean$Latitude,
      VWC = moisture_dec_clean$VWC_frac
    ),
    data.frame(
      Longitude = river_data$Longitude,
      Latitude = river_data$Latitude,
      VWC = 1.0  # River = saturated
    )
  )
  
  lon_range <- range(combined_moisture$Longitude)
  lat_range <- range(combined_moisture$Latitude)
  
  # Create interpolated raster
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
  
  cat("✓ December moisture raster created\n")
  cat("  Grid size:", length(MOISTURE_DEC_RASTER$x), "×", length(MOISTURE_DEC_RASTER$y), "\n")
} else {
  MOISTURE_DEC_RASTER <- NULL
  moisture_lookup_xy <- function(x, y) { rep(NA_real_, length(x)) }
  cat("⚠ December moisture files not found - raster not created\n")
}

# =============================================================================
# 14. CREATE MONTHLY DRIVERS
# =============================================================================

cat("\nCreating monthly climate drivers...\n")

# Air temperature by month
air_monthly <- weather_clean %>%
  mutate(month = month(TIMESTAMP)) %>%
  group_by(month) %>%
  summarise(
    air_temp_C_mean = mean(air_temp_C, na.rm = TRUE),
    air_temp_C_sd = sd(air_temp_C, na.rm = TRUE),
    .groups = "drop"
  )

# Soil temperature by month (from all available sources)
soil_monthly <- bind_rows(
  TREE_JULY %>% select(month, soil_temp_C),
  if(exists("TREE_YEAR") && nrow(TREE_YEAR) > 0) TREE_YEAR %>% select(month, soil_temp_C) else NULL,
  if(exists("SOIL_YEAR") && nrow(SOIL_YEAR) > 0) SOIL_YEAR %>% select(month, soil_temp_C) else NULL
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

# Combine drivers
DRIVERS <- tibble(month = 1:12) %>%
  left_join(air_monthly, by = "month") %>%
  left_join(soil_monthly, by = "month")

cat("✓ Monthly drivers created\n")
cat("  Months with air temp data:", sum(!is.na(DRIVERS$air_temp_C_mean)), "\n")
cat("  Months with soil temp data:", sum(!is.na(DRIVERS$soil_temp_C_mean)), "\n")

# =============================================================================
# 15. FINAL VALIDATION AND SUMMARY
# =============================================================================

cat("\n=== VALIDATION ===\n")

# Validation function
validate_columns <- function(df, required_cols, df_name) {
  if (is.null(df) || nrow(df) == 0) {
    cat("⚠", df_name, "is empty or NULL\n")
    return(FALSE)
  }
  
  missing <- setdiff(required_cols, names(df))
  if(length(missing) > 0) {
    cat("⚠", df_name, "missing columns:", paste(missing, collapse=", "), "\n")
    return(FALSE)
  } else {
    cat("✓", df_name, "has all required columns\n")
    return(TRUE)
  }
}

# Validate main datasets
validate_columns(TREE_JULY, 
                 c("tree_id", "species", "Date", "stem_flux_umol_m2_s", "x", "y"),
                 "TREE_JULY")

if (exists("TREE_YEAR") && nrow(TREE_YEAR) > 0) {
  validate_columns(TREE_YEAR,
                   c("tree_id", "species", "Date", "stem_flux_umol_m2_s", "x", "y"),
                   "TREE_YEAR")
}

if (exists("SOIL_YEAR") && nrow(SOIL_YEAR) > 0) {
  validate_columns(SOIL_YEAR,
                   c("site_id", "Date", "soil_flux_umol_m2_s", "x", "y"),
                   "SOIL_YEAR")
}

validate_columns(INVENTORY,
                 c("tree_id", "species", "dbh_m", "x", "y"),
                 "INVENTORY")

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("COORDINATE SYSTEM: WGS84 geodetic (x=longitude, y=latitude)\n")
cat("Reference point:", ref_lat, "N,", ref_lon, "W\n")
cat("Rotation angle:", rotation_angle * 180/pi, "degrees\n\n")

cat("TREE_JULY:", nrow(TREE_JULY), "measurements from",
    length(unique(TREE_JULY$tree_id)), "trees\n")
if (exists("n_2021")) cat("  2021:", n_2021, "measurements\n")
if (exists("n_2023") && n_2023 > 0) cat("  2023:", n_2023, "measurements\n")

if (exists("TREE_YEAR") && nrow(TREE_YEAR) > 0) {
  cat("TREE_YEAR:", nrow(TREE_YEAR), "measurements from",
      length(unique(TREE_YEAR$tree_id)), "trees\n")
}

if (exists("SOIL_YEAR") && nrow(SOIL_YEAR) > 0) {
  cat("SOIL_YEAR:", nrow(SOIL_YEAR), "measurements from",
      length(unique(SOIL_YEAR$site_id)), "sites\n")
}

cat("INVENTORY:", nrow(INVENTORY), "trees with GPS coordinates\n")
cat("Plot area:", round(PLOT_AREA/10000, 2), "ha\n")

# Data quality summary
cat("\n=== DATA QUALITY ===\n")

# Check coordinate coverage
tree_july_coords <- sum(!is.na(TREE_JULY$x) & !is.na(TREE_JULY$y))
cat("TREE_JULY with coordinates:", tree_july_coords, "of", nrow(TREE_JULY),
    "(", round(100*tree_july_coords/nrow(TREE_JULY), 1), "%)\n")

if (exists("TREE_YEAR") && nrow(TREE_YEAR) > 0) {
  tree_year_coords <- sum(!is.na(TREE_YEAR$x) & !is.na(TREE_YEAR$y))
  cat("TREE_YEAR with coordinates:", tree_year_coords, "of", nrow(TREE_YEAR),
      "(", round(100*tree_year_coords/nrow(TREE_YEAR), 1), "%)\n")
}

if (exists("SOIL_YEAR") && nrow(SOIL_YEAR) > 0) {
  soil_year_coords <- sum(!is.na(SOIL_YEAR$x) & !is.na(SOIL_YEAR$y))
  cat("SOIL_YEAR with coordinates:", soil_year_coords, "of", nrow(SOIL_YEAR),
      "(", round(100*soil_year_coords/nrow(SOIL_YEAR), 1), "%)\n")
}

# =============================================================================
# 16. SAVE OUTPUT
# =============================================================================

cat("\nSaving datasets...\n")

# Ensure empty dataframes if not created
if (!exists("TREE_YEAR")) TREE_YEAR <- data.frame()
if (!exists("SOIL_YEAR")) SOIL_YEAR <- data.frame()

# Create output list
rf_workflow_data <- list(
  # Main datasets
  PLACEHOLDER_TREE_JULY = TREE_JULY,
  PLACEHOLDER_TREE_YEAR = TREE_YEAR,
  PLACEHOLDER_SOIL_YEAR = SOIL_YEAR,
  PLACEHOLDER_INVENTORY = INVENTORY,
  
  # Moisture raster
  PLACEHOLDER_MOISTURE_DEC_RASTER = MOISTURE_DEC_RASTER,
  PLACEHOLDER_MOISTURE_LOOKUP_XY = moisture_lookup_xy,
  
  # Climate drivers
  PLACEHOLDER_DRIVERS = DRIVERS,
  
  # Taxonomy
  PLACEHOLDER_TAXONOMY = TAXONOMY,
  
  # Plot information
  PLACEHOLDER_PLOT_AREA = PLOT_AREA,
  
  # Geodetic transformation functions and parameters
  geodetic_transform = geodetic_transform,
  geodetic_inverse = geodetic_inverse,
  ref_lat = ref_lat,
  ref_lon = ref_lon,
  rotation_angle = rotation_angle,
  earth_radius = earth_radius,
  
  # Metadata
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

# Save the data
save(rf_workflow_data, file = "../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")

cat("✓ All datasets saved to: rf_workflow_input_data_with_2023.RData\n")
cat("\n=== DATA IMPORT COMPLETE ===\n")
cat("Completed at:", format(Sys.time()), "\n")
cat("\nReady to use with RF workflow script\n")
cat("Remember: All spatial data uses geodetic GPS coordinates\n")
cat("         x = longitude (decimal degrees)\n")
cat("         y = latitude (decimal degrees)\n")

# Final notes
cat("\nNOTES FOR DOWNSTREAM ANALYSIS:\n")
cat("1. Coordinate system is WGS84 geodetic (not projected meters)\n")
cat("2. Use geodetic_transform() to convert plot coords to GPS\n")
cat("3. Use geodetic_inverse() to convert GPS to plot coords\n")
cat("4. moisture_lookup_xy() function works with GPS coordinates\n")
cat("5. Check for missing coordinates before spatial analysis\n")