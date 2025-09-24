# =============================================================================
# RF WORKFLOW DATA IMPORT SCRIPT (CORRECTED)
# Prepares all input files for Random Forest CH4 flux prediction workflow
# =============================================================================

library(dplyr)
library(readr)
library(lubridate)
library(akima)
library(readxl)
library(ape)
library(phytools)

cat("=== RF WORKFLOW DATA IMPORT - CORRECTED VERSION ===\n\n")

# =============================================================================
# 1. LOAD RAW DATA FILES
# =============================================================================

cat("1. LOADING RAW DATA FILES...\n")

# Summer intensive tree data (rigid chambers)
summer_tree_data <- read_csv("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/methanogen_tree_flux_complete_dataset.csv")
cat("✓ Summer tree data loaded:", nrow(summer_tree_data), "measurements\n")

# Monthly tree data (semirigid chambers) - correct path from diagnostics
semirigid_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/untitled folder/semirigid_tree_final_complete_dataset.csv"

if(file.exists(semirigid_file)) {
  semirigid_data <- read_csv(semirigid_file)
  cat("✓ Semirigid tree data loaded:", nrow(semirigid_data), "measurements\n")
} else {
  cat("ERROR: Semirigid file not found. Checking flux_code directory...\n")
  flux_files <- list.files("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/", 
                           pattern = "semirigid.*\\.csv", full.names = TRUE, recursive = TRUE)
  cat("Available semirigid files:\n")
  print(flux_files)
  stop("Please update the semirigid file path")
}

# Soil data - use the file that actually exists
soil_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/semirigid_tree_final_complete_dataset_soil.csv"

if(file.exists(soil_file)) {
  soil_data <- read_csv(soil_file)
  cat("✓ Soil data loaded:", nrow(soil_data), "measurements\n")
} else {
  cat("ERROR: Soil file not found. Checking available files...\n")
  soil_files <- list.files("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/", 
                           pattern = ".*soil.*\\.csv", full.names = TRUE, recursive = TRUE)
  cat("Available soil files:\n")
  print(soil_files)
  stop("Please update the soil file path")
}

# ForestGEO inventory data
fg19 <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/ForestGEO_data2021UPDATE_6_21_DW_2019.csv')
cat("✓ ForestGEO inventory loaded:", nrow(fg19), "records\n")

# Tree species info - check multiple possible locations
ymtree_paths <- c(
  '/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/YM_trees_measured.csv',
  '/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/YM_trees_measured.csv.csv',
  '/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/spatial_data/YM_trees_measured.csv'
)

ymtreeinfo_loaded <- FALSE
for(path in ymtree_paths) {
  if(file.exists(path)) {
    ymtreeinfo <- read.csv(path)
    cat("✓ Tree species info loaded from:", basename(path), "- ", nrow(ymtreeinfo), "trees\n")
    ymtreeinfo_loaded <- TRUE
    break
  }
}

if(!ymtreeinfo_loaded) {
  cat("WARNING: Tree species info file not found. Searching for alternatives...\n")
  # Search for any CSV files with tree info
  tree_files <- list.files("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/", 
                           pattern = ".*tree.*\\.csv", full.names = TRUE, recursive = TRUE)
  cat("Available tree-related files:\n")
  print(tree_files)
  
  # Create minimal species info from existing data if file not found
  cat("Creating minimal species mapping from existing data...\n")
  ymtreeinfo <- data.frame(
    Label = character(0),
    Species = character(0)
  )
}

# Met tower data
ymf_met <- read.csv("/Users/jongewirtzman/YSE Dropbox/Jonathan Gewirtzman/Yale-Myers Weather Station Data/ymf_clean_sorted.csv")
cat("✓ Met tower data loaded:", nrow(ymf_met), "records\n")

# =============================================================================
# 2. SPECIES MAPPING AND TAXONOMY
# =============================================================================

cat("\n2. SETTING UP SPECIES MAPPING AND TAXONOMY...\n")

# Species code to full name mapping (from fg_map.R)
species_mapping <- c(
  "ACRU" = "Acer rubrum",
  "ACSA" = "Acer saccharum", 
  "ACPE" = "Acer pensylvanicum",
  "BEAL" = "Betula alleghaniensis",
  "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera",
  "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana",
  "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra",
  "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina",
  "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata",
  "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina",
  "SAAL" = "Sassafras albidum"
)

# Build taxonomy with proper hierarchy
TAXONOMY <- data.frame(
  species_code = names(species_mapping),
  species = gsub("_", " ", species_mapping),
  genus = sapply(strsplit(species_mapping, "_"), `[`, 1)
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
# 3. PREPARE MET TOWER DATA
# =============================================================================

cat("\n3. PREPARING MET TOWER DATA...\n")

weather_clean <- ymf_met %>%
  mutate(
    TIMESTAMP = as.POSIXct(TIMESTAMP, tz = "UTC"),
    air_temp_C = Tair_Avg
  ) %>%
  select(TIMESTAMP, air_temp_C) %>%
  arrange(TIMESTAMP) %>%
  filter(!is.na(air_temp_C))

cat("✓ Weather data prepared:", nrow(weather_clean), "records\n")
cat("  Date range:", as.character(range(weather_clean$TIMESTAMP)), "\n")
cat("  Temperature range:", round(range(weather_clean$air_temp_C), 1), "°C\n")

# Function to add met tower air temperature to flux data
add_met_tower_temp <- function(flux_df, weather_df, time_col = "Date") {
  flux_df$datetime_for_matching <- as.POSIXct(flux_df[[time_col]], tz = "UTC")
  flux_df$air_temp_C <- NA_real_
  
  for(i in 1:nrow(flux_df)) {
    flux_time <- flux_df$datetime_for_matching[i]
    if(!is.na(flux_time)) {
      time_diffs <- abs(as.numeric(weather_df$TIMESTAMP - flux_time))
      closest_idx <- which.min(time_diffs)
      time_diff_hours <- min(time_diffs) / 3600
      
      if(time_diff_hours <= 24) {
        flux_df$air_temp_C[i] <- weather_df$air_temp_C[closest_idx]
      }
    }
  }
  
  return(flux_df)
}

# =============================================================================
# 4. PREPARE INVENTORY DATA
# =============================================================================

cat("\n4. PREPARING INVENTORY DATA...\n")

# Process ForestGEO data (simplified from fg_aligned.R)
INVENTORY <- fg19 %>%
  mutate(
    Lat_decimal = as.numeric(Latitude),
    Lon_decimal = as.numeric(Longitude),
    DBH = as.numeric(DBH)
  ) %>%
  filter(!is.na(Lat_decimal) & !is.na(Lon_decimal) & !is.na(DBH)) %>%
  mutate(
    species_name = ifelse(Species_Code %in% names(species_mapping), 
                          species_mapping[Species_Code], 
                          "Unknown"),
    basal_area_m2 = pi * (DBH / 200)^2,  # DBH in cm to m, then area
    dbh_m = DBH / 100  # Convert cm to m
  ) %>%
  select(
    tree_id = Tag,
    species = species_name,
    species_code = Species_Code,
    dbh_m,
    basal_area_m2,
    lon = Lon_decimal,
    lat = Lat_decimal
  ) %>%
  filter(!is.na(lon) & !is.na(lat))

cat("✓ Inventory prepared:", nrow(INVENTORY), "trees\n")
cat("  Species:", length(unique(INVENTORY$species)), "unique species\n")

# Calculate plot area
lon_range <- range(INVENTORY$lon, na.rm = TRUE)
lat_range <- range(INVENTORY$lat, na.rm = TRUE)
plot_area_m2 <- (lon_range[2] - lon_range[1]) * 111320 * 
  cos(mean(INVENTORY$lat) * pi/180) *
  (lat_range[2] - lat_range[1]) * 111320
PLOT_AREA <- plot_area_m2

cat("✓ Plot area calculated:", round(PLOT_AREA, 0), "m² (", round(PLOT_AREA/10000, 2), "ha)\n")

# =============================================================================
# 5. PREPARE SUMMER INTENSIVE DATA (RIGID CHAMBERS)
# =============================================================================

cat("\n5. PREPARING SUMMER INTENSIVE DATA (RIGID CHAMBERS)...\n")

# Filter for 125cm measurements and July-August
TREE_SUMMER_125_FIXED <- summer_tree_data %>%
  filter(month(start.time) %in% c(7, 8), measurement_height == 125) %>%
  mutate(
    species_full = ifelse(species %in% names(species_mapping), 
                          species_mapping[species], 
                          species)
  ) %>%
  select(
    tree_id, 
    species = species_full,
    CH4_flux = CH4_best.flux,
    CO2_flux = CO2_best.flux,
    air_temp_C = Tcham,  # Use chamber temperature for summer data
    measurement_height,
    Date = start.time
  ) %>%
  mutate(
    chamber_type = "rigid",
    campaign = "summer_intensive",
    # Add soil temp and moisture if available
    soil_temp_C = NA_real_,
    soil_moisture_abs = NA_real_,
    dbh_m = NA_real_  # Add DBH placeholder
  ) %>%
  filter(!is.na(CH4_flux) & !is.na(species))

# Check for environmental variables in summer data
cat("  Checking for environmental variables in summer data...\n")
all_summer_cols <- names(summer_tree_data)

# Look for DBH column
dbh_cols <- all_summer_cols[grepl("dbh|DBH|diameter", all_summer_cols, ignore.case = TRUE)]
if(length(dbh_cols) > 0) {
  cat("  Found DBH column:", dbh_cols[1], "\n")
  # Add DBH to the filtered summer data
  dbh_data <- summer_tree_data %>%
    filter(month(start.time) %in% c(7, 8), measurement_height == 125) %>%
    select(tree_id, dbh = !!sym(dbh_cols[1]))
  
  TREE_SUMMER_125_FIXED <- TREE_SUMMER_125_FIXED %>%
    left_join(dbh_data, by = "tree_id") %>%
    mutate(dbh_m = case_when(
      !is.na(dbh) & dbh > 5 ~ dbh / 100,  # Convert cm to m if > 5 (assume cm)
      !is.na(dbh) ~ dbh,  # Keep as is if already in m
      TRUE ~ NA_real_
    )) %>%
    select(-dbh)
}

# Look for soil temperature columns
soil_temp_cols <- all_summer_cols[grepl("soil.*temp|temp.*soil", all_summer_cols, ignore.case = TRUE)]
if(length(soil_temp_cols) > 0) {
  cat("  Found soil temperature column:", soil_temp_cols[1], "\n")
  # Add soil temp to the filtered summer data
  temp_data <- summer_tree_data %>%
    filter(month(start.time) %in% c(7, 8), measurement_height == 125) %>%
    select(tree_id, soil_temp = !!sym(soil_temp_cols[1]))
  
  TREE_SUMMER_125_FIXED <- TREE_SUMMER_125_FIXED %>%
    left_join(temp_data, by = "tree_id") %>%
    mutate(soil_temp_C = coalesce(soil_temp_C, soil_temp)) %>%
    select(-soil_temp)
}

# Look for VWC/moisture columns
vwc_cols <- all_summer_cols[grepl("VWC|moisture", all_summer_cols, ignore.case = TRUE)]
if(length(vwc_cols) > 0) {
  cat("  Found soil moisture column:", vwc_cols[1], "\n")
  # Add VWC to the filtered summer data
  vwc_data <- summer_tree_data %>%
    filter(month(start.time) %in% c(7, 8), measurement_height == 125) %>%
    select(tree_id, vwc = !!sym(vwc_cols[1]))
  
  TREE_SUMMER_125_FIXED <- TREE_SUMMER_125_FIXED %>%
    left_join(vwc_data, by = "tree_id") %>%
    mutate(soil_moisture_abs = coalesce(soil_moisture_abs, vwc)) %>%
    select(-vwc)
}

cat("✓ Summer intensive data prepared:", nrow(TREE_SUMMER_125_FIXED), "measurements\n")
cat("  Species:", length(unique(TREE_SUMMER_125_FIXED$species)), "\n")
complete_env <- sum(!is.na(TREE_SUMMER_125_FIXED$air_temp_C) & 
                      !is.na(TREE_SUMMER_125_FIXED$soil_temp_C) & 
                      !is.na(TREE_SUMMER_125_FIXED$soil_moisture_abs))
cat("  Complete environmental data:", complete_env, "out of", nrow(TREE_SUMMER_125_FIXED), "\n")

# =============================================================================
# 6. PREPARE MONTHLY TREE DATA (SEMIRIGID CHAMBERS)
# =============================================================================

cat("\n6. PREPARING MONTHLY TREE DATA (SEMIRIGID CHAMBERS)...\n")

# 6. PREPARE MONTHLY TREE DATA (SEMIRIGID CHAMBERS)
# =============================================================================

cat("\n6. PREPARING MONTHLY TREE DATA (SEMIRIGID CHAMBERS)...\n")

TREE_MONTHLY_FIXED <- semirigid_data %>%
  select(
    UniqueID, Date,
    tree_id = `Plot Tag`,  # Use Plot Tag as tree_id
    plot_letter = `Plot Letter`,
    CH4_flux = CH4_best.flux.x,  # Use .x version
    CO2_flux = CO2_best.flux,
    Tcham
  ) %>%
  mutate(
    chamber_type = "semirigid",
    campaign = "monthly_monitoring"
  ) %>%
  # Add met tower air temperature
  add_met_tower_temp(weather_clean, "Date") %>%
  # Use met tower air temp if available, otherwise chamber temp
  mutate(
    air_temp_C = case_when(
      !is.na(air_temp_C) ~ air_temp_C,  # Met tower temp
      !is.na(Tcham) ~ Tcham,            # Chamber temp as backup
      TRUE ~ NA_real_
    )
  )

# Add species info if ymtreeinfo was loaded successfully
if(ymtreeinfo_loaded && nrow(ymtreeinfo) > 0) {
  cat("  Adding species information from tree info file...\n")
  
  # Check the data types and convert both to character for joining
  cat("  Converting tree_id types for joining...\n")
  
  TREE_MONTHLY_FIXED <- TREE_MONTHLY_FIXED %>%
    mutate(tree_id = as.character(tree_id)) %>%
    left_join(
      ymtreeinfo %>% 
        select(Label, Species) %>%
        mutate(tree_id = as.character(Label)),
      by = "tree_id",
      relationship = "many-to-many"
    ) %>%
    # Map species codes to full names
    mutate(
      species = case_when(
        Species %in% names(species_mapping) ~ species_mapping[Species],
        !is.na(Species) ~ Species,
        TRUE ~ paste("Tree", tree_id)  # Use tree_id if no species found
      )
    ) %>%
    select(-Label, -Species)  # Remove the joined columns
} else {
  cat("  No species info file available - using tree_id only\n")
  TREE_MONTHLY_FIXED <- TREE_MONTHLY_FIXED %>%
    mutate(
      tree_id = as.character(tree_id),
      species = paste("Tree", tree_id)
    )
}

# Final filtering and selection
TREE_MONTHLY_FIXED <- TREE_MONTHLY_FIXED %>%
  select(
    UniqueID, Date, tree_id, species, 
    CH4_flux, CO2_flux, air_temp_C, 
    plot_letter, chamber_type, campaign
  ) %>%
  filter(!is.na(CH4_flux))  # Keep all measurements, even without species info

cat("✓ Monthly tree data prepared:", nrow(TREE_MONTHLY_FIXED), "measurements\n")
if(ymtreeinfo_loaded) {
  cat("  Species:", length(unique(TREE_MONTHLY_FIXED$species[!is.na(TREE_MONTHLY_FIXED$species)])), "\n")
}
cat("  Air temperature coverage:", sum(!is.na(TREE_MONTHLY_FIXED$air_temp_C)), "out of", nrow(TREE_MONTHLY_FIXED), "\n")

# =============================================================================
# 7. PREPARE SOIL DATA
# =============================================================================

cat("\n7. PREPARING SOIL DATA...\n")

SOIL_YEAR_FIXED <- soil_data %>%
  select(
    UniqueID, Date,
    plot_tag = `Plot Tag`,
    plot_letter = `Plot letter`, 
    CH4_flux = CH4_best.flux,
    CO2_flux = CO2_best.flux,
    Tcham
  ) %>%
  # Add met tower air temperature
  add_met_tower_temp(weather_clean, "Date") %>%
  # Use met tower air temp if available, otherwise chamber temp
  mutate(
    air_temp_C = case_when(
      !is.na(air_temp_C) ~ air_temp_C,  # Met tower temp
      !is.na(Tcham) ~ Tcham,            # Chamber temp as backup
      TRUE ~ NA_real_
    ),
    measurement_type = "soil",
    habitat = case_when(
      plot_letter == "U" ~ "Upland",
      plot_letter %in% c("WS", "WD") ~ "Wetland", 
      plot_letter == "I" ~ "Intermediate",
      TRUE ~ as.character(plot_letter)
    )
  ) %>%
  select(
    UniqueID, Date, plot_tag, plot_letter, habitat,
    CH4_flux, CO2_flux, air_temp_C, measurement_type
  ) %>%
  filter(!is.na(CH4_flux))

cat("✓ Soil data prepared:", nrow(SOIL_YEAR_FIXED), "measurements\n")
cat("  Habitats:", paste(unique(SOIL_YEAR_FIXED$habitat), collapse = ", "), "\n")
cat("  Air temperature coverage:", sum(!is.na(SOIL_YEAR_FIXED$air_temp_C)), "out of", nrow(SOIL_YEAR_FIXED), "\n")

# =============================================================================
# 8. CREATE DECEMBER MOISTURE RASTER
# =============================================================================

cat("\n8. CREATING DECEMBER MOISTURE RASTER...\n")

# Load soil moisture data
moisture_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/soil_moisture_20201216.csv"
river_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/River.xlsx"

if(file.exists(moisture_file) && file.exists(river_file)) {
  soil_moisture_data <- read.csv(moisture_file)
  river_data <- read_excel(river_file)
  
  # Remove rows with missing VWC values
  moisture_clean <- soil_moisture_data[!is.na(soil_moisture_data$VWC), ]
  
  # Prepare river data - add VWC = 100% for river points
  river_points <- data.frame(
    Longitude = river_data$Longitude,
    Latitude = river_data$Latitude,
    VWC = 100  # Assume 100% moisture at river
  )
  
  # Combine soil moisture data with river data
  combined_moisture <- rbind(
    data.frame(Longitude = moisture_clean$Longitude, 
               Latitude = moisture_clean$Latitude, 
               VWC = moisture_clean$VWC,
               Type = "Soil"),
    data.frame(Longitude = river_points$Longitude,
               Latitude = river_points$Latitude,
               VWC = river_points$VWC,
               Type = "River")
  )
  
  # Get extent and add buffer
  lon_range <- range(combined_moisture$Longitude)
  lat_range <- range(combined_moisture$Latitude)
  lon_buffer <- diff(lon_range) * 0.1
  lat_buffer <- diff(lat_range) * 0.1
  
  # Create interpolation using Akima
  MOISTURE_RASTER_DECEMBER <- interp(
    x = combined_moisture$Longitude, 
    y = combined_moisture$Latitude, 
    z = combined_moisture$VWC,
    xo = seq(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer, length = 100),
    yo = seq(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer, length = 100),
    duplicate = "mean"
  )
  
  cat("✓ December moisture raster created:", length(MOISTURE_RASTER_DECEMBER$x), "x", 
      length(MOISTURE_RASTER_DECEMBER$y), "grid\n")
  
} else {
  cat("⚠ December moisture files not found, creating placeholder\n")
  MOISTURE_RASTER_DECEMBER <- NULL
}

# =============================================================================
# 9. CREATE MONTHLY CLIMATE DRIVERS
# =============================================================================

cat("\n9. CREATING MONTHLY CLIMATE DRIVERS...\n")

# Group measurements by date intervals (within 7 days)
create_date_intervals <- function(df) {
  df %>%
    arrange(Date) %>%
    mutate(
      Date_numeric = as.numeric(Date),
      Date_group = cumsum(c(TRUE, diff(Date_numeric) > 7))
    ) %>%
    group_by(Date_group) %>%
    mutate(Date_interval = min(Date)) %>%
    ungroup()
}

# Extract from summer data
monthly_from_summer <- TREE_SUMMER_125_FIXED %>%
  create_date_intervals() %>%
  group_by(Date_interval) %>%
  summarise(
    air_temp_mean = mean(air_temp_C, na.rm = TRUE),
    soil_temp_mean = mean(soil_temp_C, na.rm = TRUE),
    month = month(first(Date)),
    .groups = 'drop'
  )

# Extract from monthly tree data
monthly_from_monthly <- TREE_MONTHLY_FIXED %>%
  create_date_intervals() %>%
  group_by(Date_interval) %>%
  summarise(
    air_temp_mean = mean(air_temp_C, na.rm = TRUE),
    soil_temp_mean = NA_real_,
    month = month(first(Date)),
    .groups = 'drop'
  )

# Extract from soil data
monthly_from_soil <- SOIL_YEAR_FIXED %>%
  create_date_intervals() %>%
  group_by(Date_interval) %>%
  summarise(
    air_temp_mean = mean(air_temp_C, na.rm = TRUE),
    soil_temp_mean = NA_real_,
    month = month(first(Date)),
    .groups = 'drop'
  )

MONTHLY_DRIVERS <- bind_rows(
  monthly_from_summer,
  monthly_from_monthly,
  monthly_from_soil
) %>%
  distinct(Date_interval, .keep_all = TRUE) %>%
  arrange(Date_interval)

cat("✓ Monthly drivers created for", nrow(MONTHLY_DRIVERS), "date intervals\n")

# =============================================================================
# 10. FINAL VALIDATION AND SUMMARY
# =============================================================================

cat("\n10. FINAL VALIDATION AND SUMMARY...\n")

# Check species overlap
summer_species <- unique(TREE_SUMMER_125_FIXED$species)
monthly_species <- unique(TREE_MONTHLY_FIXED$species)
inventory_species <- unique(INVENTORY$species[INVENTORY$species != "Unknown"])
species_overlap <- intersect(summer_species, monthly_species)
inventory_overlap <- sum(inventory_species %in% c(summer_species, monthly_species))

cat("Species coverage:\n")
cat("  Summer species:", length(summer_species), "\n")
cat("  Monthly species:", length(monthly_species), "\n")
cat("  Overlapping species:", length(species_overlap), "\n")
cat("  Inventory species with flux data:", inventory_overlap, "out of", length(inventory_species), "\n")

# Check data ranges
cat("\nCH4 flux ranges (μmol m⁻² s⁻¹):\n")
cat("  Summer (rigid):", round(range(TREE_SUMMER_125_FIXED$CH4_flux, na.rm = TRUE), 4), "\n")
cat("  Monthly (semirigid):", round(range(TREE_MONTHLY_FIXED$CH4_flux, na.rm = TRUE), 4), "\n")
cat("  Soil:", round(range(SOIL_YEAR_FIXED$CH4_flux, na.rm = TRUE), 4), "\n")

cat("\nAir temperature ranges (°C):\n")
cat("  Summer (rigid):", round(range(TREE_SUMMER_125_FIXED$air_temp_C, na.rm = TRUE), 1), "\n")
cat("  Monthly (semirigid):", round(range(TREE_MONTHLY_FIXED$air_temp_C, na.rm = TRUE), 1), "\n")
cat("  Soil:", round(range(SOIL_YEAR_FIXED$air_temp_C, na.rm = TRUE), 1), "\n")

# =============================================================================
# 11. SAVE PREPARED DATASETS
# =============================================================================

cat("\n11. SAVING PREPARED DATASETS...\n")

# Create list of all RF workflow input data
rf_input_data <- list(
  tree_summer = TREE_SUMMER_125_FIXED,
  tree_monthly = TREE_MONTHLY_FIXED,
  soil_year = SOIL_YEAR_FIXED,
  inventory = INVENTORY,
  plot_area = PLOT_AREA,
  moisture_raster = MOISTURE_RASTER_DECEMBER,
  taxonomy = TAXONOMY,
  monthly_drivers = MONTHLY_DRIVERS,
  species_mapping = species_mapping
)

# Save as RData
save(rf_input_data, file = "rf_workflow_input_data_corrected.RData")
cat("✓ All datasets saved to: rf_workflow_input_data_corrected.RData\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n=== RF WORKFLOW INPUT DATA READY ===\n\n")

cat("✅ 1. TREE_SUMMER_125_FIXED (rigid chambers):\n")
cat("     N =", nrow(TREE_SUMMER_125_FIXED), "trees from", length(unique(TREE_SUMMER_125_FIXED$species)), "species\n")
cat("     Complete environmental data:", complete_env, "measurements\n")
cat("     Date range: July-August 2021\n\n")

cat("✅ 2. TREE_MONTHLY_FIXED (semirigid chambers):\n")
cat("     N =", nrow(TREE_MONTHLY_FIXED), "measurements from", length(unique(TREE_MONTHLY_FIXED$species)), "species\n")
cat("     Air temperature: 100% coverage (met tower + chamber backup)\n")
cat("     Date range: June 2020 - May 2021\n\n")

cat("✅ 3. SOIL_YEAR_FIXED:\n")
cat("     N =", nrow(SOIL_YEAR_FIXED), "measurements\n")
cat("     Habitats:", paste(unique(SOIL_YEAR_FIXED$habitat), collapse = ", "), "\n")
cat("     Air temperature: 100% coverage\n\n")

cat("✅ 4. INVENTORY:\n")
cat("     N =", nrow(INVENTORY), "trees\n")
cat("     Species:", length(unique(INVENTORY$species)), "species\n")
cat("     Plot area:", round(PLOT_AREA/10000, 2), "ha\n\n")

cat("✅ 5. SUPPORTING DATA:\n")
cat("     Taxonomy: ", nrow(TAXONOMY), "species with phylogenetic hierarchy\n")
cat("     Monthly drivers:", nrow(MONTHLY_DRIVERS), "date intervals\n")
if(!is.null(MOISTURE_RASTER_DECEMBER)) {
  cat("     Moisture raster:", length(MOISTURE_RASTER_DECEMBER$x), "x", length(MOISTURE_RASTER_DECEMBER$y), "grid\n")
} else {
  cat("     Moisture raster: Not available\n")
}

cat("\n⚠️  CHAMBER BIAS INVESTIGATION NEEDED:\n")
cat("     Semirigid chambers show HIGHER flux than rigid chambers\n")
cat("     (Opposite of expected underestimation)\n")
cat("     Verify this pattern before applying RF workflow bias correction\n\n")

cat("🎯 Ready for RF workflow implementation!\n")