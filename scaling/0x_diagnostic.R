# =============================================================================
# DIAGNOSTIC CODE FOR ENVIRONMENTAL DATA AND SPECIES MAPPING
# =============================================================================

cat("=== DETAILED ENVIRONMENTAL DATA AUDIT ===\n")

# 1. Check what environmental data we actually have
cat("1. SUMMER INTENSIVE (RIGID) ENVIRONMENTAL DATA:\n")
if(exists("TREE_SUMMER_125")) {
  cat("Available columns:\n")
  env_cols <- names(TREE_SUMMER_125)[grepl("temp|moisture|VWC|soil", names(TREE_SUMMER_125), ignore.case = TRUE)]
  print(env_cols)
  
  # Check if soil temp and moisture are in the summer data
  has_soil_temp <- any(grepl("soil.*temp", names(TREE_SUMMER_125), ignore.case = TRUE))
  has_soil_moisture <- any(grepl("VWC|moisture", names(TREE_SUMMER_125), ignore.case = TRUE))
  
  cat("Has soil temperature:", has_soil_temp, "\n")
  cat("Has soil moisture:", has_soil_moisture, "\n")
  
  if(has_soil_temp) {
    soil_temp_cols <- names(TREE_SUMMER_125)[grepl("soil.*temp", names(TREE_SUMMER_125), ignore.case = TRUE)]
    cat("Soil temp columns:", paste(soil_temp_cols, collapse = ", "), "\n")
  }
  
  if(has_soil_moisture) {
    moisture_cols <- names(TREE_SUMMER_125)[grepl("VWC|moisture", names(TREE_SUMMER_125), ignore.case = TRUE)]
    cat("Soil moisture columns:", paste(moisture_cols, collapse = ", "), "\n")
  }
}

cat("\n2. YEARLONG (SEMIRIGID) ENVIRONMENTAL DATA:\n")
if(exists("semirigid_data")) {
  # Check if soil data is in the tree file or separate soil file
  cat("Tree measurements file environmental columns:\n")
  tree_env_cols <- names(semirigid_data)[grepl("temp|moisture|VWC", names(semirigid_data), ignore.case = TRUE)]
  print(tree_env_cols)
}

# Check soil dataset for environmental data
if(exists("soil_data")) {
  cat("Soil measurements file environmental columns:\n")
  soil_env_cols <- names(soil_data)[grepl("temp|moisture|VWC", names(soil_data), ignore.case = TRUE)]
  print(soil_env_cols)
}

# =============================================================================
# SPECIES NAME MAPPING FIX
# =============================================================================

cat("\n=== SPECIES NAME MAPPING ISSUE ===\n")

# Show the species name formats in each dataset
cat("Summer flux species (sample):\n")
print(head(unique(TREE_SUMMER_125$species), 10))

cat("\nMonthly flux species (sample):\n") 
print(head(unique(TREE_MONTHLY$species), 10))

cat("\nInventory species (sample):\n")
print(head(unique(INVENTORY$species), 10))

# Create species mapping from your fg_map.R code
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

# Apply mapping to make species names consistent
if(exists("TREE_SUMMER_125")) {
  TREE_SUMMER_125 <- TREE_SUMMER_125 %>%
    mutate(species_full = ifelse(species %in% names(species_mapping), 
                                 species_mapping[species], 
                                 species))
}

if(exists("TREE_MONTHLY")) {
  TREE_MONTHLY <- TREE_MONTHLY %>%
    mutate(species_full = ifelse(species %in% names(species_mapping), 
                                 species_mapping[species], 
                                 species))
}

# Re-check species overlap after mapping
cat("\n=== SPECIES OVERLAP AFTER MAPPING ===\n")
if(exists("TREE_SUMMER_125") && exists("TREE_MONTHLY") && exists("INVENTORY")) {
  summer_species_full <- unique(TREE_SUMMER_125$species_full[!is.na(TREE_SUMMER_125$species_full)])
  monthly_species_full <- unique(TREE_MONTHLY$species_full[!is.na(TREE_MONTHLY$species_full)])
  inventory_species_clean <- unique(INVENTORY$species[!is.na(INVENTORY$species)])
  
  cat("Summer species (mapped):", length(summer_species_full), "\n")
  cat("Monthly species (mapped):", length(monthly_species_full), "\n")
  cat("Inventory species:", length(inventory_species_clean), "\n")
  
  all_flux_species <- unique(c(summer_species_full, monthly_species_full))
  overlap <- sum(inventory_species_clean %in% all_flux_species)
  
  cat("Inventory species with flux data:", overlap, "out of", length(inventory_species_clean), "\n")
  
  if(overlap > 0) {
    overlapping_species <- inventory_species_clean[inventory_species_clean %in% all_flux_species]
    cat("Overlapping species:\n")
    print(overlapping_species)
  }
}

# =============================================================================
# MET TOWER DATA LOADING
# =============================================================================

cat("\n=== LOADING MET TOWER DATA ===\n")

# Load met tower data from your soil_auxfile.R path
met_file <- "/Users/jongewirtzman/YSE Dropbox/Jonathan Gewirtzman/Yale-Myers Weather Station Data/ymf_clean_sorted.csv"

if(file.exists(met_file)) {
  ymf_met <- read.csv(met_file)
  cat("Met tower data loaded:", nrow(ymf_met), "records\n")
  
  # Check column names
  cat("Met tower columns:\n")
  print(names(ymf_met))
  
  # Prepare weather data (following your soil_auxfile.R approach)
  weather_clean <- ymf_met %>%
    mutate(
      TIMESTAMP = as.POSIXct(TIMESTAMP, tz = "UTC"),
      air_temp_C = Tair_Avg  # Use average air temperature
    ) %>%
    select(TIMESTAMP, air_temp_C) %>%
    arrange(TIMESTAMP) %>%
    filter(!is.na(air_temp_C))
  
  cat("Weather data range:", as.character(range(weather_clean$TIMESTAMP)), "\n")
  cat("Temperature range:", round(range(weather_clean$air_temp_C, na.rm = TRUE), 1), "°C\n")
  
} else {
  cat("ERROR: Met tower file not found at:", met_file, "\n")
  cat("Please check the file path\n")
}

# =============================================================================
# SEMIRIGID DATA COLUMN SELECTION
# =============================================================================

cat("\n=== SEMIRIGID CH4 FLUX COLUMN SELECTION ===\n")

if(exists("semirigid_data")) {
  # Check both CH4 flux columns
  ch4_x_available <- "CH4_best.flux.x" %in% names(semirigid_data)
  ch4_y_available <- "CH4_best.flux.y" %in% names(semirigid_data)
  
  cat("CH4_best.flux.x available:", ch4_x_available, "\n")
  cat("CH4_best.flux.y available:", ch4_y_available, "\n")
  
  if(ch4_x_available) {
    cat("CH4_best.flux.x summary:\n")
    print(summary(semirigid_data$CH4_best.flux.x))
    missing_x <- sum(is.na(semirigid_data$CH4_best.flux.x))
    cat("Missing values in .x:", missing_x, "\n")
  }
  
  if(ch4_y_available) {
    cat("CH4_best.flux.y summary:\n") 
    print(summary(semirigid_data$CH4_best.flux.y))
    missing_y <- sum(is.na(semirigid_data$CH4_best.flux.y))
    cat("Missing values in .y:", missing_y, "\n")
  }
  
  # Use CH4_best.flux.x as you specified
  cat("Using CH4_best.flux.x for monthly tree measurements\n")
}












# Check if environmental data exists in original source files
cat("=== CHECKING ORIGINAL SEMIRIGID FILES ===\n")

# Check what columns are actually in the semirigid files
if(exists("semirigid_data")) {
  cat("All columns in semirigid tree data:\n")
  all_cols <- names(semirigid_data)
  temp_cols <- all_cols[grepl("temp|Temp|TEMP", all_cols)]
  moisture_cols <- all_cols[grepl("VWC|moisture|Moisture", all_cols)]
  
  cat("Temperature-related columns:\n")
  print(temp_cols)
  cat("Moisture-related columns:\n") 
  print(moisture_cols)
  
  # Check if Tcham is the air temperature column
  if("Tcham" %in% names(semirigid_data)) {
    cat("Tcham (chamber temp) summary:\n")
    print(summary(semirigid_data$Tcham))
  }
}

# Check soil data file
if(exists("soil_data")) {
  cat("\nAll columns in soil data:\n")
  soil_all_cols <- names(soil_data)
  soil_temp_cols <- soil_all_cols[grepl("temp|Temp|TEMP", soil_all_cols)]
  soil_moisture_cols <- soil_all_cols[grepl("VWC|moisture|Moisture", soil_all_cols)]
  
  cat("Soil temperature-related columns:\n")
  print(soil_temp_cols)
  cat("Soil moisture-related columns:\n")
  print(soil_moisture_cols)
}

# Function to match air temperature from met tower to flux measurements
add_met_tower_temp <- function(flux_df, weather_df, time_col = "Date") {
  # Convert flux dates to datetime for matching
  flux_df$datetime_for_matching <- as.POSIXct(flux_df[[time_col]], tz = "UTC")
  
  # For each flux measurement, find closest weather observation
  flux_df$air_temp_C <- NA_real_
  
  for(i in 1:nrow(flux_df)) {
    flux_time <- flux_df$datetime_for_matching[i]
    if(!is.na(flux_time)) {
      # Find closest weather data (within 24 hours)
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

# Apply to semirigid data if weather data is available
if(exists("weather_clean") && exists("TREE_MONTHLY")) {
  cat("Adding met tower air temperature to monthly tree data...\n")
  TREE_MONTHLY <- add_met_tower_temp(TREE_MONTHLY, weather_clean, "Date")
  
  # Check results
  matched_temp <- sum(!is.na(TREE_MONTHLY$air_temp_C))
  cat("Air temperature matched for", matched_temp, "out of", nrow(TREE_MONTHLY), "measurements\n")
}

if(exists("weather_clean") && exists("SOIL_YEAR")) {
  cat("Adding met tower air temperature to soil data...\n")
  SOIL_YEAR <- add_met_tower_temp(SOIL_YEAR, weather_clean, "Date")
  
  matched_temp_soil <- sum(!is.na(SOIL_YEAR$air_temp_C))
  cat("Air temperature matched for", matched_temp_soil, "out of", nrow(SOIL_YEAR), "measurements\n")
}

# Update datasets with fixed species names and environmental data
cat("=== UPDATING DATASETS WITH CORRECTED DATA ===\n")

# Use species_full (mapped names) instead of species
TREE_SUMMER_125_FIXED <- TREE_SUMMER_125 %>%
  select(
    tree_id, 
    species = species_full,  # Use mapped full names
    dbh,
    CH4_flux = CH4_flux,
    CO2_flux = CO2_flux,
    air_temp_C = air_temp,   # Rename for consistency
    soil_temp_C = soil_temp, 
    soil_moisture_abs = VWC,
    plot,
    measurement_height
  ) %>%
  mutate(
    chamber_type = "rigid",
    campaign = "summer_intensive",
    Date = as.Date("2021-07-30")  # Representative date
  ) %>%
  filter(!is.na(CH4_flux) & !is.na(species))

cat("TREE_SUMMER_125_FIXED prepared:", nrow(TREE_SUMMER_125_FIXED), "measurements\n")

# Monthly tree data with met tower air temp
TREE_MONTHLY_FIXED <- TREE_MONTHLY %>%
  select(
    UniqueID, Date,
    tree_id,  
    species = species_full,  # Use mapped full names
    CH4_flux = CH4_best.flux.x,  # Use .x column
    CO2_flux = CO2_best.flux,
    air_temp_C,  # From met tower
    plot_letter
  ) %>%
  mutate(
    chamber_type = "semirigid",
    campaign = "monthly_monitoring"
  ) %>%
  filter(!is.na(CH4_flux) & !is.na(species))

cat("TREE_MONTHLY_FIXED prepared:", nrow(TREE_MONTHLY_FIXED), "measurements\n")

# Print summary of fixed datasets
cat("\n=== FIXED DATASETS SUMMARY ===\n")
cat("Summer intensive:\n")
cat("- Species:", length(unique(TREE_SUMMER_125_FIXED$species)), "\n")
cat("- Measurements:", nrow(TREE_SUMMER_125_FIXED), "\n")
cat("- Environmental data complete:", 
    sum(!is.na(TREE_SUMMER_125_FIXED$air_temp_C) & 
          !is.na(TREE_SUMMER_125_FIXED$soil_temp_C) & 
          !is.na(TREE_SUMMER_125_FIXED$soil_moisture_abs)), "\n")

cat("\nMonthly monitoring:\n")
cat("- Species:", length(unique(TREE_MONTHLY_FIXED$species)), "\n") 
cat("- Measurements:", nrow(TREE_MONTHLY_FIXED), "\n")
cat("- Air temp from met tower:", sum(!is.na(TREE_MONTHLY_FIXED$air_temp_C)), "\n")

# Check what columns are actually in TREE_MONTHLY
cat("=== CHECKING TREE_MONTHLY COLUMN NAMES ===\n")
if(exists("TREE_MONTHLY")) {
  cat("All columns in TREE_MONTHLY:\n")
  print(names(TREE_MONTHLY))
  
  # Look for CH4 flux columns specifically
  ch4_cols <- names(TREE_MONTHLY)[grepl("CH4", names(TREE_MONTHLY))]
  cat("\nCH4-related columns:\n")
  print(ch4_cols)
  
  # Look for CO2 flux columns
  co2_cols <- names(TREE_MONTHLY)[grepl("CO2", names(TREE_MONTHLY))]
  cat("\nCO2-related columns:\n")
  print(co2_cols)
}

# Also check the original semirigid_data to see what columns it has
cat("\n=== CHECKING ORIGINAL SEMIRIGID_DATA COLUMNS ===\n")
if(exists("semirigid_data")) {
  cat("CH4 columns in original semirigid_data:\n")
  semirigid_ch4_cols <- names(semirigid_data)[grepl("CH4", names(semirigid_data))]
  print(semirigid_ch4_cols)
}
# Fix the monthly tree data preparation with correct column names
cat("=== FIXING MONTHLY TREE DATA PREPARATION ===\n")

# Go back to original semirigid_data and prepare correctly
if(exists("semirigid_data")) {
  TREE_MONTHLY_FIXED <- semirigid_data %>%
    select(
      UniqueID, Date,
      tree_id,  
      plot_letter,
      CH4_flux = CH4_best.flux.x,  # This should exist in semirigid_data
      CO2_flux = CO2_best.flux,
      Tcham  # Air temperature from chamber
    ) %>%
    # Add species info
    left_join(
      ymtreeinfo %>% 
        select(Label, Species) %>%
        mutate(tree_id = as.character(Label)),
      by = "tree_id"
    ) %>%
    # Map species codes to full names
    mutate(
      species = case_when(
        Species %in% names(species_mapping) ~ species_mapping[Species],
        !is.na(Species) ~ Species,
        TRUE ~ NA_character_
      ),
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
    ) %>%
    select(
      UniqueID, Date, tree_id, species, 
      CH4_flux, CO2_flux, air_temp_C, 
      plot_letter, chamber_type, campaign
    ) %>%
    filter(!is.na(CH4_flux) & !is.na(species))
  
  cat("TREE_MONTHLY_FIXED prepared:", nrow(TREE_MONTHLY_FIXED), "measurements\n")
}

# Also fix the soil data preparation
cat("=== FIXING SOIL DATA PREPARATION ===\n")

if(exists("soil_data")) {
  SOIL_YEAR_FIXED <- soil_data %>%
    select(
      UniqueID, Date,
      plot_tag = `Plot Tag`,
      plot_letter = `Plot letter`, 
      CH4_flux = CH4_best.flux,
      CO2_flux = CO2_best.flux,
      Tcham  # Air temperature from chamber
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
  
  cat("SOIL_YEAR_FIXED prepared:", nrow(SOIL_YEAR_FIXED), "measurements\n")
}

# Print summary of all fixed datasets
cat("\n=== FINAL FIXED DATASETS SUMMARY ===\n")

if(exists("TREE_SUMMER_125_FIXED")) {
  cat("1. TREE_SUMMER_125_FIXED (rigid chambers):\n")
  cat("   Measurements:", nrow(TREE_SUMMER_125_FIXED), "\n")
  cat("   Species:", length(unique(TREE_SUMMER_125_FIXED$species)), "\n")
  cat("   Complete environmental data:", 
      sum(!is.na(TREE_SUMMER_125_FIXED$air_temp_C) & 
            !is.na(TREE_SUMMER_125_FIXED$soil_temp_C) & 
            !is.na(TREE_SUMMER_125_FIXED$soil_moisture_abs)), "\n")
}

if(exists("TREE_MONTHLY_FIXED")) {
  cat("\n2. TREE_MONTHLY_FIXED (semirigid chambers):\n")
  cat("   Measurements:", nrow(TREE_MONTHLY_FIXED), "\n")
  cat("   Species:", length(unique(TREE_MONTHLY_FIXED$species)), "\n")
  cat("   Air temperature data:", sum(!is.na(TREE_MONTHLY_FIXED$air_temp_C)), "\n")
}

if(exists("SOIL_YEAR_FIXED")) {
  cat("\n3. SOIL_YEAR_FIXED:\n")
  cat("   Measurements:", nrow(SOIL_YEAR_FIXED), "\n")
  cat("   Habitats:", paste(unique(SOIL_YEAR_FIXED$habitat), collapse = ", "), "\n")
  cat("   Air temperature data:", sum(!is.na(SOIL_YEAR_FIXED$air_temp_C)), "\n")
}







# Check what columns are actually in TREE_MONTHLY
cat("=== CHECKING TREE_MONTHLY COLUMN NAMES ===\n")
if(exists("TREE_MONTHLY")) {
  cat("All columns in TREE_MONTHLY:\n")
  print(names(TREE_MONTHLY))
  
  # Look for CH4 flux columns specifically
  ch4_cols <- names(TREE_MONTHLY)[grepl("CH4", names(TREE_MONTHLY))]
  cat("\nCH4-related columns:\n")
  print(ch4_cols)
}

# Fix the monthly tree data preparation with correct column names
cat("=== FIXING MONTHLY TREE DATA WITH CORRECT COLUMNS ===\n")

# Create TREE_MONTHLY_FIXED from the correct source with correct column names
TREE_MONTHLY_FIXED <- TREE_MONTHLY %>%
  select(
    UniqueID, Date,
    tree_id,  
    species = species_full,  # Use mapped full names
    CH4_flux,  # This should be the column name that exists
    # CO2_flux = CO2_flux,  # Check if this exists
    air_temp_C,  # From met tower (already added)
    plot_letter
  ) %>%
  mutate(
    chamber_type = "semirigid",
    campaign = "monthly_monitoring"
  ) %>%
  filter(!is.na(CH4_flux) & !is.na(species))

cat("TREE_MONTHLY_FIXED prepared:", nrow(TREE_MONTHLY_FIXED), "measurements\n")

# Also fix the soil data
cat("=== FIXING SOIL DATA ===\n")

SOIL_YEAR_FIXED <- SOIL_YEAR %>%
  select(
    UniqueID, Date,
    plot_tag, plot_letter, habitat,
    CH4_flux,  # This should be the correct column name
    air_temp_C,  # From met tower (already added)
    measurement_type
  ) %>%
  filter(!is.na(CH4_flux))

cat("SOIL_YEAR_FIXED prepared:", nrow(SOIL_YEAR_FIXED), "measurements\n")

# Print summary of all fixed datasets
cat("\n=== FINAL FIXED DATASETS SUMMARY ===\n")

cat("1. TREE_SUMMER_125_FIXED (rigid chambers):\n")
cat("   Measurements:", nrow(TREE_SUMMER_125_FIXED), "\n")
cat("   Species:", length(unique(TREE_SUMMER_125_FIXED$species)), "\n")
cat("   Complete environmental data:", 
    sum(!is.na(TREE_SUMMER_125_FIXED$air_temp_C) & 
          !is.na(TREE_SUMMER_125_FIXED$soil_temp_C) & 
          !is.na(TREE_SUMMER_125_FIXED$soil_moisture_abs)), "\n")

cat("\n2. TREE_MONTHLY_FIXED (semirigid chambers):\n")
cat("   Measurements:", nrow(TREE_MONTHLY_FIXED), "\n")
cat("   Species:", length(unique(TREE_MONTHLY_FIXED$species)), "\n")
cat("   Air temperature data:", sum(!is.na(TREE_MONTHLY_FIXED$air_temp_C)), "\n")

cat("\n3. SOIL_YEAR_FIXED:\n")
cat("   Measurements:", nrow(SOIL_YEAR_FIXED), "\n")
cat("   Habitats:", paste(unique(SOIL_YEAR_FIXED$habitat), collapse = ", "), "\n")
cat("   Air temperature data:", sum(!is.na(SOIL_YEAR_FIXED$air_temp_C)), "\n")

# Check data ranges for validation
cat("\n=== DATA VALIDATION ===\n")

cat("CH4 flux ranges:\n")
cat("Summer (rigid):", round(range(TREE_SUMMER_125_FIXED$CH4_flux, na.rm = TRUE), 4), "\n")
cat("Monthly (semirigid):", round(range(TREE_MONTHLY_FIXED$CH4_flux, na.rm = TRUE), 4), "\n")
cat("Soil:", round(range(SOIL_YEAR_FIXED$CH4_flux, na.rm = TRUE), 4), "\n")

cat("\nAir temperature ranges:\n")
cat("Summer (rigid):", round(range(TREE_SUMMER_125_FIXED$air_temp_C, na.rm = TRUE), 1), "°C\n")
cat("Monthly (semirigid):", round(range(TREE_MONTHLY_FIXED$air_temp_C, na.rm = TRUE), 1), "°C\n")
cat("Soil:", round(range(SOIL_YEAR_FIXED$air_temp_C, na.rm = TRUE), 1), "°C\n")

# Check species overlap between datasets
all_species_summer <- unique(TREE_SUMMER_125_FIXED$species)
all_species_monthly <- unique(TREE_MONTHLY_FIXED$species)
species_overlap <- intersect(all_species_summer, all_species_monthly)

cat("\nSpecies overlap between summer and monthly:\n")
cat("Summer species:", length(all_species_summer), "\n")
cat("Monthly species:", length(all_species_monthly), "\n")
cat("Overlapping species:", length(species_overlap), "\n")
print(species_overlap)

