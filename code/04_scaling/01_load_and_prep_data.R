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
  semirigid_tree = "../../data/processed/flux/semirigid_tree_final_complete_dataset.csv",
  soil = "../../data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv",
  inventory = "../../data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_2019.csv",
  met_tower = "../../data/raw/weather/ymf_clean_sorted.csv",
  moisture_extra = "../../data/raw/ipad_data/Cleaned data/soilmoisture_total.csv",
  moisture_december = "../../data/raw/spatial_data/soil_moisture_20201216.csv",
  river = "../../data/raw/spatial_data/River.xlsx",
  plot_locations = "../../data/raw/spatial_data/plots.csv"
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
semirigid_moisture <- load_safe(paths$moisture_extra, read.csv)
plot_locations <- load_safe(paths$plot_locations, read.csv)

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

# Validate coordinate transformation
coord_ranges <- list(
  x = range(INVENTORY$x, na.rm = TRUE),
  y = range(INVENTORY$y, na.rm = TRUE)
)

# Calculate plot area
lon_dist_m <- (coord_ranges$x[2] - coord_ranges$x[1]) * 111320 * cos(mean(INVENTORY$y, na.rm=TRUE) * pi/180)
lat_dist_m <- (coord_ranges$y[2] - coord_ranges$y[1]) * 111320
PLOT_AREA <- lon_dist_m * lat_dist_m

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
    month = month(Date),
    year = year(Date)
  ) %>%
  select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s, 
         air_temp_C, soil_temp_C, soil_moisture_abs, dbh_m, chamber_type)

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
TREE_JULY_2021 <- bind_rows(TREE_JULY_125, KALMIA_DATA, OTHER_75CM)

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
    tree_2023_matched <- tree_2023_clean %>%
      mutate(tree_id = as.character(tree_tag_num)) %>%  # Use tag number directly, no prefix
      left_join(INVENTORY %>% select(tree_id, x, y, dbh_m), by = "tree_id") %>%
      mutate(
        species = species_mapping[species_code],
        Date = date_clean,
        chamber_type = "rigid"
      ) %>%
      filter(!is.na(x), !is.na(y), !is.na(species))
    
    # Add met tower temperature
    tree_2023_matched <- add_met_tower_temp(tree_2023_matched, weather_clean, "Date")
    
    # Create final 2023 dataset
    TREE_JULY_2023 <- tree_2023_matched %>%
      mutate(
        soil_temp_C = NA_real_, 
        soil_moisture_abs = NA_real_
      ) %>%
      select(tree_id, species, species_code, Date, month, year, stem_flux_umol_m2_s,
             air_temp_C = air_temp_C_tower, soil_temp_C, soil_moisture_abs, 
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
      left_join(
        moisture_clean %>%
          mutate(moisture_date = as.Date(Date)) %>%
          select(moisture_date, plot_letter, soil_temp_C_moisture = soil_temp_C, 
                 soil_moisture_abs_moisture = soil_moisture_abs),
        by = "plot_letter",
        relationship = "many-to-many"
      ) %>%
      mutate(
        # Calculate date difference
        date_diff = abs(as.numeric(Date_only - moisture_date)),
        # Only keep moisture within 7 days
        soil_temp_C_matched = ifelse(date_diff <= 7, soil_temp_C_moisture, NA),
        soil_moisture_abs_matched = ifelse(date_diff <= 7, soil_moisture_abs_moisture, NA)
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
      final_lat = coalesce(plot_lat, letter_lat),
      final_lon = coalesce(plot_lon, letter_lon),
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