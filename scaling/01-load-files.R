# =============================================================================
# FILE 1: TREE JULY CAMPAIGN DATA (Rigid Chambers)
# =============================================================================

library(dplyr)
library(readr)
library(lubridate)

# Load the July tree flux data
cat("=== LOADING JULY TREE FLUX DATA ===\n")
july_tree_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/methanogen_tree_flux_complete_dataset.csv"

if(file.exists(july_tree_file)) {
  july_tree_data <- read_csv(july_tree_file)
  
  cat("File loaded successfully\n")
  cat("Dimensions:", nrow(july_tree_data), "rows x", ncol(july_tree_data), "columns\n\n")
  
  # Check column names
  cat("=== COLUMN NAMES ===\n")
  print(names(july_tree_data))
  
  # Check for environmental variables
  cat("\n=== CHECKING FOR ENVIRONMENTAL VARIABLES ===\n")
  env_cols <- c("air_temp", "soil_temp", "moisture", "VWC", "temperature", "temp")
  matching_cols <- names(july_tree_data)[grepl(paste(env_cols, collapse = "|"), 
                                               names(july_tree_data), 
                                               ignore.case = TRUE)]
  cat("Potential environmental columns found:\n")
  print(matching_cols)
  
  # Check for flux columns
  cat("\n=== CHECKING FOR FLUX COLUMNS ===\n")
  flux_cols <- names(july_tree_data)[grepl("CH4|CO2|flux", names(july_tree_data), ignore.case = TRUE)]
  cat("Flux columns found:\n")
  print(flux_cols)
  
  # Check date range
  cat("\n=== DATE INFORMATION ===\n")
  if("Date" %in% names(july_tree_data)) {
    july_tree_data <- july_tree_data %>% mutate(Date = as.Date(Date))
    cat("Date range:", as.character(range(july_tree_data$Date, na.rm = TRUE)), "\n")
    
    # Check which dates are in July
    july_dates <- july_tree_data %>% 
      filter(month(Date) == 7) %>% 
      summarise(n = n(), dates = toString(unique(Date)))
    cat("July measurements:", july_dates$n, "\n")
    cat("July dates:", july_dates$dates, "\n")
  }
  
  # Check for tree identifiers
  cat("\n=== TREE IDENTIFIERS ===\n")
  id_cols <- names(july_tree_data)[grepl("tree|tag|id", names(july_tree_data), ignore.case = TRUE)]
  cat("Tree ID columns found:\n")
  print(id_cols)
  
  # Check for species information
  cat("\n=== SPECIES INFORMATION ===\n")
  species_cols <- names(july_tree_data)[grepl("species", names(july_tree_data), ignore.case = TRUE)]
  cat("Species columns found:\n")
  print(species_cols)
  if(length(species_cols) > 0) {
    cat("Unique species codes:\n")
    print(unique(july_tree_data[[species_cols[1]]]))
  }
  
  # Check for measurement height
  cat("\n=== MEASUREMENT HEIGHT ===\n")
  height_cols <- names(july_tree_data)[grepl("height|measurement", names(july_tree_data), ignore.case = TRUE)]
  cat("Height-related columns found:\n")
  print(height_cols)
  
  # Show first few rows
  cat("\n=== FIRST 3 ROWS OF DATA ===\n")
  print(head(july_tree_data, 3))
  
} else {
  cat("ERROR: File not found at:", july_tree_file, "\n")
  cat("Please check the file path\n")
}
# =============================================================================
# FILE 1: TREE SUMMER CAMPAIGN DATA (Rigid Chambers - July/August)
# =============================================================================

library(dplyr)
library(readr)
library(lubridate)

# Load the summer tree flux data
cat("=== LOADING SUMMER (JULY-AUGUST) TREE FLUX DATA ===\n")
summer_tree_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/methanogen_tree_flux_complete_dataset.csv"

if(file.exists(summer_tree_file)) {
  summer_tree_data <- read_csv(summer_tree_file)
  
  cat("File loaded successfully\n")
  cat("Dimensions:", nrow(summer_tree_data), "rows x", ncol(summer_tree_data), "columns\n\n")
  
  # Check column names
  cat("=== COLUMN NAMES ===\n")
  print(names(summer_tree_data))
  
  # Check for environmental variables
  cat("\n=== CHECKING FOR ENVIRONMENTAL VARIABLES ===\n")
  env_cols <- c("air_temp", "soil_temp", "moisture", "VWC", "temperature", "temp", "Tcham", "Pcham")
  matching_cols <- names(summer_tree_data)[grepl(paste(env_cols, collapse = "|"), 
                                                 names(summer_tree_data), 
                                                 ignore.case = TRUE)]
  cat("Potential environmental columns found:\n")
  print(matching_cols)
  
  # Check for flux columns
  cat("\n=== CHECKING FOR FLUX COLUMNS ===\n")
  flux_cols <- names(summer_tree_data)[grepl("CH4|CO2|flux", names(summer_tree_data), ignore.case = TRUE)]
  cat("Flux columns found:\n")
  print(flux_cols)
  
  # Check date range and filter for July-August
  cat("\n=== DATE INFORMATION ===\n")
  if("Date" %in% names(summer_tree_data)) {
    summer_tree_data <- summer_tree_data %>% mutate(Date = as.Date(Date))
    cat("Full date range:", as.character(range(summer_tree_data$Date, na.rm = TRUE)), "\n")
    
    # Check July-August measurements
    summer_months_data <- summer_tree_data %>% 
      filter(month(Date) %in% c(7, 8))
    
    cat("\nJuly-August measurements:", nrow(summer_months_data), "out of", nrow(summer_tree_data), "total\n")
    
    # Breakdown by month
    monthly_summary <- summer_tree_data %>%
      mutate(Month = month(Date, label = TRUE)) %>%
      count(Month)
    cat("\nMeasurements by month:\n")
    print(monthly_summary)
  }
  
  # Check for tree identifiers
  cat("\n=== TREE IDENTIFIERS ===\n")
  id_cols <- names(summer_tree_data)[grepl("tree|tag|id", names(summer_tree_data), ignore.case = TRUE)]
  cat("Tree ID columns found:\n")
  print(id_cols)
  
  # Check for species information
  cat("\n=== SPECIES INFORMATION ===\n")
  species_cols <- names(summer_tree_data)[grepl("species", names(summer_tree_data), ignore.case = TRUE)]
  cat("Species columns found:\n")
  print(species_cols)
  if(length(species_cols) > 0) {
    cat("\nUnique species codes:\n")
    print(sort(unique(summer_tree_data[[species_cols[1]]])))
  }
  
  # Check for measurement height
  cat("\n=== MEASUREMENT HEIGHT ===\n")
  height_cols <- names(summer_tree_data)[grepl("height|measurement", names(summer_tree_data), ignore.case = TRUE)]
  cat("Height-related columns found:\n")
  print(height_cols)
  if("measurement_height" %in% names(summer_tree_data)) {
    cat("Height values present:\n")
    print(table(summer_tree_data$measurement_height))
  }
  
  # Check for coordinate columns
  cat("\n=== COORDINATE COLUMNS ===\n")
  coord_cols <- names(summer_tree_data)[grepl("lat|lon|coord|x|y", names(summer_tree_data), ignore.case = TRUE)]
  cat("Coordinate columns found:\n")
  print(coord_cols)
  
  # Check for plot/site information
  cat("\n=== PLOT/SITE INFORMATION ===\n")
  plot_cols <- names(summer_tree_data)[grepl("plot|site", names(summer_tree_data), ignore.case = TRUE)]
  cat("Plot/site columns found:\n")
  print(plot_cols)
  
  # Show structure of first few rows (transposed for readability)
  cat("\n=== DATA STRUCTURE (First 2 rows, transposed) ===\n")
  print(t(head(summer_tree_data, 2)))
  
} else {
  cat("ERROR: File not found at:", summer_tree_file, "\n")
  cat("Please check the file path\n")
}
# =============================================================================
# FILE 1 PROCESSING: Extract dates and check for missing environmental data
# =============================================================================

# Extract date from start.time and check date range
cat("\n=== PROCESSING SUMMER TREE DATA ===\n")

summer_tree_data <- summer_tree_data %>%
  mutate(
    Date = as.Date(start.time),
    Month = month(Date),
    Year = year(Date)
  )

# Check date range
cat("Date range:", as.character(range(summer_tree_data$Date, na.rm = TRUE)), "\n")

# Filter for July-August only
summer_tree_filtered <- summer_tree_data %>%
  filter(Month %in% c(7, 8))

cat("July-August measurements:", nrow(summer_tree_filtered), "out of", nrow(summer_tree_data), "total\n")

# Check breakdown by month and year
date_summary <- summer_tree_data %>%
  count(Year, Month) %>%
  arrange(Year, Month)
cat("\nMeasurements by year-month:\n")
print(date_summary)

# Check environmental data availability
cat("\n=== ENVIRONMENTAL DATA CHECK ===\n")
cat("Tcham (air temp) available:", sum(!is.na(summer_tree_data$Tcham)), "out of", nrow(summer_tree_data), "\n")
cat("Tcham range:", range(summer_tree_data$Tcham, na.rm = TRUE), "°C\n")

# Check if we have soil moisture/temp elsewhere
cat("\nMissing from this dataset:\n")
cat("- Soil temperature\n")
cat("- Soil moisture (VWC)\n")
cat("- Coordinates (lat/lon)\n")

# Check plot information
cat("\n=== PLOT/LOCATION INFORMATION ===\n")
plot_summary <- summer_tree_data %>%
  count(plot) %>%
  arrange(desc(n))
cat("Unique plots:\n")
print(plot_summary)

# Check species distribution
cat("\n=== SPECIES DISTRIBUTION ===\n")
species_summary <- summer_tree_data %>%
  filter(Month %in% c(7, 8)) %>%
  count(species) %>%
  arrange(desc(n))
print(species_summary)

# Check for coordinate extraction from plot names
cat("\n=== CHECKING FOR COORDINATE PATTERNS ===\n")
# Some plot names might contain coordinates or plot identifiers
unique_plots <- unique(summer_tree_data$plot)
cat("Sample plot names:\n")
print(head(unique_plots, 10))

# Save the filtered summer data
TREE_SUMMER <- summer_tree_filtered %>%
  select(
    UniqueID, Date, plot, tree_id, species, 
    measurement_height, 
    CH4_best.flux, CO2_best.flux,
    CH4_LM.r2, CH4_HM.r2, CO2_LM.r2, CO2_HM.r2,
    Tcham, Pcham, chamber_id
  ) %>%
  mutate(
    chamber_type = "rigid",
    dataset = "summer_intensive"
  )

cat("\n=== SUMMER TREE DATA PREPARED ===\n")
cat("Saved as TREE_SUMMER with", nrow(TREE_SUMMER), "rows\n")
cat("Columns retained:", paste(names(TREE_SUMMER), collapse = ", "), "\n")

# Check if we need to load additional environmental data
cat("\n=== NEXT STEPS ===\n")
cat("Need to obtain from other sources:\n")
cat("1. Soil temperature - from field measurements or weather station\n")
cat("2. Soil moisture (VWC) - from field measurements\n")
cat("3. Geographic coordinates - from plot mapping files\n")
# =============================================================================
# CHECK FOR COMPLETE DATASET WITH ENVIRONMENTAL DATA
# =============================================================================

cat("=== LOOKING FOR COMPLETE DATASET WITH ENVIRONMENTAL DATA ===\n")

# Try the merged dataset from height effect script
merged_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/merged_tree_dataset_final.csv"

if(file.exists(merged_file)) {
  merged_data <- read_csv(merged_file)
  
  cat("Merged dataset found and loaded!\n")
  cat("Dimensions:", nrow(merged_data), "rows x", ncol(merged_data), "columns\n\n")
  
  # Check for environmental columns
  cat("=== CHECKING FOR ENVIRONMENTAL VARIABLES ===\n")
  env_patterns <- c("VWC", "moisture", "soil.*temp", "air.*temp", "temp", "humidity")
  env_cols <- names(merged_data)[grepl(paste(env_patterns, collapse = "|"), 
                                       names(merged_data), 
                                       ignore.case = TRUE)]
  cat("Environmental columns found:\n")
  print(env_cols)
  
  # Check for coordinate columns
  cat("\n=== CHECKING FOR COORDINATES ===\n")
  coord_patterns <- c("lat", "lon", "longitude", "latitude", "coord")
  coord_cols <- names(merged_data)[grepl(paste(coord_patterns, collapse = "|"), 
                                         names(merged_data), 
                                         ignore.case = TRUE)]
  cat("Coordinate columns found:\n")
  print(coord_cols)
  
  # Check column names
  cat("\n=== ALL COLUMN NAMES ===\n")
  print(names(merged_data))
  
  # Check if this has tree IDs that match
  cat("\n=== CHECKING TREE IDENTIFIERS ===\n")
  id_cols <- names(merged_data)[grepl("tree|tag|id", names(merged_data), ignore.case = TRUE)]
  cat("ID columns found:\n")
  print(id_cols)
  
} else {
  cat("Merged dataset not found at expected location.\n")
  cat("Need to check for other complete dataset files.\n")
}

# Also check if there's a complete flux dataset with environmental data
cat("\n=== CHECKING FOR OTHER COMPLETE DATASETS ===\n")

# List all CSV files in the flux_code directory
flux_dir <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/"
csv_files <- list.files(flux_dir, pattern = "*.csv", full.names = FALSE)

cat("CSV files in flux_code directory:\n")
print(csv_files)

# Look for files with "complete" or "final" in the name
complete_files <- csv_files[grepl("complete|final", csv_files, ignore.case = TRUE)]
cat("\nPotential complete dataset files:\n")
print(complete_files)
# =============================================================================
# FILE 2: MONTHLY TREE DATA (Semirigid Chambers)
# =============================================================================

cat("\n=== LOADING MONTHLY TREE DATA (SEMIRIGID) ===\n")

# Check if semirigid tree dataset exists in flux_code
semirigid_files <- list.files(
  "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/",
  pattern = "semirigid.*tree.*\\.csv",
  full.names = TRUE
)

cat("Semirigid tree files found:\n")
print(basename(semirigid_files))

# Try the untitled folder mentioned in your scripts
semirigid_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/untitled folder/semirigid_tree_final_complete_dataset.csv"

if(file.exists(semirigid_file)) {
  semirigid_data <- read_csv(semirigid_file)
  
  cat("\nSemirigid dataset loaded!\n")
  cat("Dimensions:", nrow(semirigid_data), "rows x", ncol(semirigid_data), "columns\n\n")
  
  # Check for date information
  cat("=== DATE INFORMATION ===\n")
  date_cols <- names(semirigid_data)[grepl("date|Date|time|Time", names(semirigid_data), ignore.case = TRUE)]
  cat("Date/time columns found:\n")
  print(date_cols)
  
  if("Date" %in% names(semirigid_data)) {
    semirigid_data <- semirigid_data %>% mutate(Date = as.Date(Date))
    cat("\nDate range:", as.character(range(semirigid_data$Date, na.rm = TRUE)), "\n")
    
    # Check monthly distribution
    monthly_dist <- semirigid_data %>%
      mutate(YearMonth = format(Date, "%Y-%m")) %>%
      count(YearMonth) %>%
      arrange(YearMonth)
    cat("\nMeasurements by year-month:\n")
    print(monthly_dist)
  }
  
  # Check for environmental variables
  cat("\n=== ENVIRONMENTAL VARIABLES ===\n")
  env_patterns <- c("VWC", "moisture", "soil.*temp", "air.*temp", "Tcham", "temp")
  env_cols <- names(semirigid_data)[grepl(paste(env_patterns, collapse = "|"), 
                                          names(semirigid_data), 
                                          ignore.case = TRUE)]
  cat("Environmental columns found:\n")
  print(env_cols)
  
  # Check for tree and species info
  cat("\n=== TREE AND SPECIES INFO ===\n")
  tree_cols <- names(semirigid_data)[grepl("tree|tag|species", names(semirigid_data), ignore.case = TRUE)]
  cat("Tree/species columns:\n")
  print(tree_cols)
  
  # Check for flux columns
  cat("\n=== FLUX COLUMNS ===\n")
  flux_cols <- names(semirigid_data)[grepl("CH4.*flux|CO2.*flux", names(semirigid_data))]
  cat("Main flux columns:\n")
  print(head(flux_cols, 10))
  
  # Sample first few rows of key columns
  if(all(c("Date", "tree_id", "species") %in% names(semirigid_data))) {
    cat("\n=== SAMPLE DATA ===\n")
    sample_data <- semirigid_data %>%
      select(Date, tree_id, species, matches("CH4_best.flux"), matches("VWC|moisture|temp")) %>%
      head(3)
    print(sample_data)
  }
  
} else {
  cat("Semirigid file not found at expected location.\n")
  cat("Checking alternative locations...\n")
}

# Now let's match the merged dataset with flux measurements to get complete data
cat("\n=== CREATING COMPLETE DATASET WITH ENVIRONMENTAL DATA ===\n")

# The merged_data has environmental variables but aggregated by tree
# The summer_tree_data has individual measurements
# Let's see if we can match them

if(exists("merged_data") && exists("summer_tree_data")) {
  # Check tree IDs in both datasets
  cat("Tree IDs in summer flux data:", length(unique(summer_tree_data$tree_id)), "unique trees\n")
  cat("Tree IDs in merged environmental data:", length(unique(merged_data$tree_id)), "unique trees\n")
  
  # Check overlap
  common_trees <- intersect(summer_tree_data$tree_id, merged_data$tree_id)
  cat("Trees in common:", length(common_trees), "\n")
  
  if(length(common_trees) > 0) {
    cat("Sample of common tree IDs:", head(common_trees, 10), "\n")
  }
}
# =============================================================================
# EXAMINE SEMIRIGID DATA MORE CAREFULLY
# =============================================================================

cat("=== EXAMINING SEMIRIGID (MONTHLY) DATA STRUCTURE ===\n")

# Check all column names to find tree/species identifiers
cat("\nAll column names in semirigid data:\n")
print(names(semirigid_data))

# Check what Plot Tag and Plot Letter contain
cat("\n=== PLOT IDENTIFIERS ===\n")
cat("Unique Plot Tags:\n")
print(unique(semirigid_data$`Plot Tag`))

cat("\nUnique Plot Letters:\n")
print(unique(semirigid_data$`Plot Letter`))

# Check if there's tree information in other columns
cat("\n=== CHECKING FOR TREE INFO IN OTHER COLUMNS ===\n")
if("Notes" %in% names(semirigid_data)) {
  sample_notes <- head(unique(semirigid_data$Notes), 20)
  cat("Sample Notes (might contain tree info):\n")
  print(sample_notes)
}

# Check UniqueID structure - might contain tree info
cat("\n=== UNIQUE ID STRUCTURE ===\n")
sample_ids <- head(semirigid_data$UniqueID, 10)
cat("Sample UniqueIDs:\n")
print(sample_ids)

# Parse UniqueID to extract components
semirigid_data <- semirigid_data %>%
  mutate(
    # Try to extract tree info from UniqueID or Plot Tag
    tree_id = as.character(`Plot Tag`),
    plot_letter = `Plot Letter`
  )

# Now let's prepare both datasets for the RF workflow
cat("\n=== PREPARING DATASETS FOR RF WORKFLOW ===\n")

# 1. SUMMER (RIGID) CHAMBERS - July-August intensive
TREE_SUMMER <- summer_tree_data %>%
  filter(month(Date) %in% c(7, 8)) %>%
  select(
    UniqueID, Date, plot, tree_id, species,
    measurement_height,
    CH4_best.flux,  # This is at various heights
    CO2_best.flux,
    CH4_LM.r2, CH4_HM.r2, 
    CO2_LM.r2, CO2_HM.r2,
    Tcham  # Air temperature from chamber
  ) %>%
  mutate(
    chamber_type = "rigid",
    campaign = "summer_intensive"
  )

cat("TREE_SUMMER prepared:", nrow(TREE_SUMMER), "measurements\n")

# For the merged data with flux by height, extract the 125cm measurements
TREE_SUMMER_125 <- merged_data %>%
  select(
    tree_id, 
    species = species_id,
    dbh,
    CH4_flux = CH4_best.flux_125cm,  # 125cm height flux
    CO2_flux = CO2_best.flux_125cm,
    air_temp = Temp_Air_125cm,
    VWC = VWC_mean,
    soil_temp = SoilTemp_mean,
    plot
  ) %>%
  mutate(
    measurement_height = 125,
    chamber_type = "rigid",
    campaign = "summer_intensive"
  ) %>%
  filter(!is.na(CH4_flux))  # Keep only rows with flux data

cat("TREE_SUMMER_125 prepared:", nrow(TREE_SUMMER_125), "trees with 125cm flux\n")

# 2. MONTHLY (SEMIRIGID) CHAMBERS
TREE_MONTHLY <- semirigid_data %>%
  select(
    UniqueID, Date,
    tree_id,  # Created from Plot Tag
    plot_letter,
    CH4_flux = CH4_best.flux.x,  
    CO2_flux = CO2_best.flux,
    CH4_LM_r2 = CH4_LM.r2.x,
    CH4_HM_r2 = CH4_HM.r2.x,
    CO2_LM_r2 = CO2_LM.r2,
    CO2_HM_r2 = CO2_HM.r2,
    Tcham
  ) %>%
  mutate(
    chamber_type = "semirigid",
    campaign = "monthly_monitoring"
  )

cat("TREE_MONTHLY prepared:", nrow(TREE_MONTHLY), "measurements\n")

# Check date ranges
cat("\n=== DATE RANGES ===\n")
cat("Summer intensive (rigid):", 
    as.character(range(TREE_SUMMER$Date, na.rm = TRUE)), "\n")
cat("Monthly monitoring (semirigid):", 
    as.character(range(TREE_MONTHLY$Date, na.rm = TRUE)), "\n")

# Save these as our first two input files
cat("\n=== FILES PREPARED ===\n")
cat("1. TREE_SUMMER_125: Summer intensive with environmental data\n")
cat("2. TREE_MONTHLY: Monthly monitoring (needs environmental data)\n")

# =============================================================================
# JOIN SPECIES INFO TO MONTHLY TREE DATA
# =============================================================================

cat("=== ADDING SPECIES INFO TO MONTHLY TREES ===\n")

# Load the tree info
ymtreeinfo <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/YM_trees_measured.csv.csv')

cat("Tree info loaded:", nrow(ymtreeinfo), "trees\n")
cat("Species in monthly trees:", unique(ymtreeinfo$Species), "\n")

# Join species info to monthly data
TREE_MONTHLY <- TREE_MONTHLY %>%
  left_join(
    ymtreeinfo %>% 
      select(Label, Species, Site, D_stem, Latitude, Longitude) %>%
      mutate(tree_id = as.character(Label)),
    by = "tree_id"
  ) %>%
  rename(
    species = Species,
    dbh = D_stem,
    lat = Latitude,
    lon = Longitude
  )

# Check the join
cat("\n=== JOIN RESULTS ===\n")
cat("Monthly measurements with species info:", sum(!is.na(TREE_MONTHLY$species)), "out of", nrow(TREE_MONTHLY), "\n")
cat("Unique species in monthly data:", paste(unique(TREE_MONTHLY$species[!is.na(TREE_MONTHLY$species)]), collapse = ", "), "\n")

# =============================================================================
# FIX SOIL DATA LOADING
# =============================================================================

cat("=== FIXING SOIL DATA LOADING ===\n")

# Check actual column names in soil data
cat("Actual column names in soil data:\n")
soil_cols <- names(soil_data)
print(soil_cols[grepl("Plot|plot|Tag|tag", soil_cols, ignore.case = TRUE)])

# Prepare soil dataset with correct column names
SOIL_YEAR <- soil_data %>%
  select(
    UniqueID, Date,
    plot_tag = `Plot Tag`,      # Without the period
    plot_letter = `Plot letter`, # Without the period
    CH4_flux = CH4_best.flux,
    CO2_flux = CO2_best.flux,
    CH4_LM_r2 = CH4_LM.r2,
    CH4_HM_r2 = CH4_HM.r2,
    CO2_LM_r2 = CO2_LM.r2,
    CO2_HM_r2 = CO2_HM.r2,
    air_temp = Tcham
  ) %>%
  mutate(
    measurement_type = "soil",
    habitat = case_when(
      plot_letter == "U" ~ "Upland",
      plot_letter %in% c("WS", "WD") ~ "Wetland",
      plot_letter == "I" ~ "Intermediate",
      TRUE ~ as.character(plot_letter)
    )
  )

cat("SOIL_YEAR prepared:", nrow(SOIL_YEAR), "measurements\n")

# =============================================================================
# LOAD FORESTGEO DATA DIRECTLY FROM fg_aligned.R CODE
# =============================================================================

cat("\n=== LOADING FORESTGEO INVENTORY ===\n")

# Load the data files directly
fg19 <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/ForestGEO_data2021UPDATE_6_21_DW_2019.csv')
fgplot <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/ForestGEO_data2021UPDATE_6_21_DW_byplot.csv')
fgtag <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/ForestGEO_data2021UPDATE_6_21_DW_bytag.csv')

cat("ForestGEO files loaded\n")

# Quick processing to get the essential inventory data
# Using simplified version of fg_aligned.R logic
library(tidyr)

# Extract coordinates from fg19
fg19_coords <- fg19 %>%
  select(Tag, Species_Code, PX, PY, Latitude, Longitude, DBH, Status, Quadrat, Sub_Quadrat) %>%
  rename(Tag_ID = Tag, Species = Species_Code) %>%
  mutate(
    Lat_decimal = as.numeric(Latitude),
    Lon_decimal = as.numeric(Longitude),
    DBH = as.numeric(DBH)
  ) %>%
  filter(!is.na(Lat_decimal) & !is.na(Lon_decimal) & !is.na(DBH))

# Simple species mapping
species_mapping <- c(
  "ACRU" = "Acer rubrum", "ACSA" = "Acer saccharum", 
  "BEAL" = "Betula alleghaniensis", "BELE" = "Betula lenta",
  "FAGR" = "Fagus grandifolia", "FRAM" = "Fraxinus americana",
  "PIST" = "Pinus strobus", "QURU" = "Quercus rubra",
  "QUAL" = "Quercus alba", "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata", "PRSE" = "Prunus serotina"
)

# Create inventory
INVENTORY <- fg19_coords %>%
  mutate(
    species_name = ifelse(Species %in% names(species_mapping), 
                          species_mapping[Species], 
                          "Unknown"),
    basal_area_m2 = pi * (DBH / 200)^2  # DBH in cm to m, then area
  ) %>%
  select(
    tree_id = Tag_ID,
    species = species_name,
    species_code = Species,
    dbh = DBH,
    basal_area_m2,
    lon = Lon_decimal,
    lat = Lat_decimal
  ) %>%
  filter(!is.na(lon) & !is.na(lat))

cat("INVENTORY prepared:", nrow(INVENTORY), "trees\n")
cat("Species in inventory:", length(unique(INVENTORY$species)), "unique species\n")

# Calculate plot area
lon_range <- range(INVENTORY$lon, na.rm = TRUE)
lat_range <- range(INVENTORY$lat, na.rm = TRUE)
plot_area_m2 <- (lon_range[2] - lon_range[1]) * 111320 * 
  cos(mean(INVENTORY$lat) * pi/180) *
  (lat_range[2] - lat_range[1]) * 111320
PLOT_AREA <- plot_area_m2

cat("Plot area:", round(PLOT_AREA, 0), "m² (", round(PLOT_AREA/10000, 2), "ha)\n")

# =============================================================================
# FILE 5: DECEMBER MOISTURE RASTER
# =============================================================================

cat("\n=== LOADING DECEMBER MOISTURE RASTER ===\n")

# Source the map.R script to get the extended moisture interpolation
# First check if we can load pre-saved moisture data
moisture_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/spatial_data/soil_moisture_20201216.csv"

if(file.exists(moisture_file)) {
  december_moisture <- read.csv(moisture_file)
  cat("December moisture data loaded:", nrow(december_moisture), "points\n")
  
  # Check what we have
  cat("Columns:", paste(names(december_moisture), collapse = ", "), "\n")
  cat("VWC range:", range(december_moisture$VWC, na.rm = TRUE), "%\n")
  
  # We'll need to interpolate this to create a raster
  # For now, save the raw points
  MOISTURE_DECEMBER_POINTS <- december_moisture %>%
    select(Longitude, Latitude, VWC) %>%
    filter(!is.na(VWC))
  
  cat("December moisture points prepared:", nrow(MOISTURE_DECEMBER_POINTS), "\n")
} else {
  cat("December moisture file not found\n")
}



# =============================================================================
# FILE 6: DECEMBER MOISTURE RASTER (from extended interpolation)
# =============================================================================

cat("=== CREATING DECEMBER MOISTURE RASTER ===\n")

# We need to create the extended moisture interpolation from your map.R code
# This uses the extended Akima interpolation that covers the full ForestGEO area
# =============================================================================
# FILE 6: DECEMBER MOISTURE RASTER (CORRECTED - from extended interpolation)
# =============================================================================

cat("=== CREATING DECEMBER MOISTURE RASTER (CORRECTED) ===\n")

library(akima)
library(readxl)
library(lubridate)

# Load the correct moisture data file from your mapping scripts
moisture_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/soil_moisture_20201216.csv"

# Check if file exists first
if(!file.exists(moisture_file)) {
  cat("ERROR: Moisture file not found at:", moisture_file, "\n")
  cat("Please check the file path\n")
} else {
  
  # Read the soil moisture data (this matches your mapping script)
  soil_moisture_data <- read.csv(moisture_file)
  cat("Soil moisture data loaded:", nrow(soil_moisture_data), "points\n")
  
  # Read the river data (also from your mapping script)
  river_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/Yale Myers Methane Project/spatial_data/River.xlsx"
  
  if(file.exists(river_file)) {
    river_data <- read_excel(river_file)
    cat("River data loaded:", nrow(river_data), "points\n")
    
    # Remove rows with missing VWC values (from your mapping script)
    moisture_clean <- soil_moisture_data[!is.na(soil_moisture_data$VWC), ]
    
    # Prepare river data - add VWC = 100% for river points (from your mapping script)
    river_points <- data.frame(
      Longitude = river_data$Longitude,
      Latitude = river_data$Latitude,
      VWC = 100  # Assume 100% moisture at river
    )
    
    # Combine soil moisture data with river data for interpolation (exactly as in your mapping script)
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
    
    cat("Combined moisture points:", nrow(combined_moisture), "\n")
    cat("VWC range:", round(min(combined_moisture$VWC), 2), "to", round(max(combined_moisture$VWC), 2), "\n")
    
    # Get extent of combined data and add buffer (from your mapping script)
    lon_range <- range(combined_moisture$Longitude)
    lat_range <- range(combined_moisture$Latitude)
    lon_buffer <- diff(lon_range) * 0.1  # 10% buffer
    lat_buffer <- diff(lat_range) * 0.1
    
    # Create interpolation using akima package with combined data (soil + river)
    # This is the CORRECT approach from your mapping scripts
    MOISTURE_RASTER_DECEMBER <- interp(
      x = combined_moisture$Longitude, 
      y = combined_moisture$Latitude, 
      z = combined_moisture$VWC,
      xo = seq(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer, length = 100),
      yo = seq(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer, length = 100),
      duplicate = "mean"
    )
    
    cat("Moisture raster created successfully: ", length(MOISTURE_RASTER_DECEMBER$x), "x", 
        length(MOISTURE_RASTER_DECEMBER$y), "grid\n")
    
    # Convert to dataframe for checking (optional)
    moisture_df_check <- expand.grid(Longitude = MOISTURE_RASTER_DECEMBER$x, 
                                     Latitude = MOISTURE_RASTER_DECEMBER$y)
    moisture_df_check$VWC_predicted <- as.vector(MOISTURE_RASTER_DECEMBER$z)
    moisture_df_check <- moisture_df_check[!is.na(moisture_df_check$VWC_predicted), ]
    
    cat("Interpolated moisture grid points:", nrow(moisture_df_check), "\n")
    cat("Interpolated VWC range:", round(range(moisture_df_check$VWC_predicted), 2), "\n")
    
  } else {
    cat("ERROR: River file not found at:", river_file, "\n")
  }
}

# Note: This creates the same extended moisture interpolation that you use successfully 
# in your mapping scripts, which includes both soil moisture measurements AND river 
# influence (rivers set to 100% VWC). This is much more accurate than using only 
# the December 2020 point measurements.

cat("Moisture raster created: ", length(MOISTURE_RASTER_DECEMBER$x), "x", 
    length(MOISTURE_RASTER_DECEMBER$y), "grid\n")

# =============================================================================
# FILE 7: SPECIES TAXONOMY (using phylogenetic tree)
# =============================================================================

cat("\n=== BUILDING SPECIES TAXONOMY FROM PHYLOGENY ===\n")

library(ape)
library(phytools)

# Load phylogenetic tree
tree_file <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/metadata/PhytoPhylo"

if(file.exists(tree_file)) {
  phylo_tree <- read.tree(tree_file)
  cat("Phylogenetic tree loaded with", length(phylo_tree$tip.label), "tips\n")
  
  # Get species in our data
  all_species_codes <- unique(c(
    TREE_SUMMER_125$species[!is.na(TREE_SUMMER_125$species)],
    TREE_MONTHLY$species[!is.na(TREE_MONTHLY$species)]
  ))
  
  # Species code to binomial mapping
  species_mapping <- c(
    "ACRU" = "Acer_rubrum", "ACSA" = "Acer_saccharum", 
    "BEAL" = "Betula_alleghaniensis", "BELE" = "Betula_lenta",
    "BEPA" = "Betula_papyrifera",
    "FAGR" = "Fagus_grandifolia", "FRAM" = "Fraxinus_americana",
    "PIST" = "Pinus_strobus", "QURU" = "Quercus_rubra",
    "QUAL" = "Quercus_alba", "QUVE" = "Quercus_velutina",
    "TSCA" = "Tsuga_canadensis", "TSLA" = "Tsuga_canadensis",
    "CAOV" = "Carya_ovata", "PRSE" = "Prunus_serotina",
    "SAAL" = "Sassafras_albidum", "KALA" = "Kalmia_latifolia"
  )
  
  # Build taxonomy with proper hierarchy
  TAXONOMY <- data.frame(
    species_code = names(species_mapping),
    species = gsub("_", " ", species_mapping),
    genus = sapply(strsplit(species_mapping, "_"), `[`, 1)
  )
  
  # Add family based on genus
  TAXONOMY$family <- case_when(
    TAXONOMY$genus == "Acer" ~ "Sapindaceae",
    TAXONOMY$genus == "Betula" ~ "Betulaceae", 
    TAXONOMY$genus == "Fagus" ~ "Fagaceae",
    TAXONOMY$genus == "Quercus" ~ "Fagaceae",
    TAXONOMY$genus == "Fraxinus" ~ "Oleaceae",
    TAXONOMY$genus == "Pinus" ~ "Pinaceae",
    TAXONOMY$genus == "Tsuga" ~ "Pinaceae",
    TAXONOMY$genus == "Carya" ~ "Juglandaceae",
    TAXONOMY$genus == "Prunus" ~ "Rosaceae",
    TAXONOMY$genus == "Sassafras" ~ "Lauraceae",
    TAXONOMY$genus == "Kalmia" ~ "Ericaceae",
    TRUE ~ "Unknown"
  )
  
  # Add order
  TAXONOMY$order <- case_when(
    TAXONOMY$family %in% c("Sapindaceae") ~ "Sapindales",
    TAXONOMY$family %in% c("Betulaceae", "Fagaceae", "Juglandaceae") ~ "Fagales",
    TAXONOMY$family == "Oleaceae" ~ "Lamiales",
    TAXONOMY$family == "Pinaceae" ~ "Pinales",
    TAXONOMY$family == "Rosaceae" ~ "Rosales",
    TAXONOMY$family == "Lauraceae" ~ "Laurales",
    TAXONOMY$family == "Ericaceae" ~ "Ericales",
    TRUE ~ "Unknown"
  )
  
  # Add class
  TAXONOMY$class <- ifelse(TAXONOMY$family == "Pinaceae", "Pinopsida", "Magnoliopsida")
  
  cat("Taxonomy built for", nrow(TAXONOMY), "species\n")
  print(TAXONOMY)
  
} else {
  cat("Phylogenetic tree file not found\n")
  # Create basic taxonomy without phylogeny
  TAXONOMY <- data.frame(species_code = character(), species = character())
}

# =============================================================================
# FILE 8: MONTHLY CLIMATE DRIVERS (from measurements)
# =============================================================================

cat("\n=== EXTRACTING MONTHLY CLIMATE DRIVERS ===\n")

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

# Extract from all datasets
monthly_from_summer <- TREE_SUMMER_125 %>%
  mutate(Date = as.Date("2021-07-30")) %>%
  summarise(
    Date_interval = unique(Date),
    air_temp_mean = mean(air_temp, na.rm = TRUE),
    soil_temp_mean = mean(soil_temp, na.rm = TRUE),
    month = 7
  )

monthly_from_monthly <- TREE_MONTHLY %>%
  create_date_intervals() %>%
  group_by(Date_interval) %>%
  summarise(
    air_temp_mean = mean(Tcham, na.rm = TRUE),
    soil_temp_mean = NA_real_,
    month = month(first(Date)),
    .groups = 'drop'
  )

monthly_from_soil <- SOIL_YEAR %>%
  create_date_intervals() %>%
  group_by(Date_interval) %>%
  summarise(
    air_temp_mean = mean(air_temp, na.rm = TRUE),
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

cat("Monthly drivers created for", nrow(MONTHLY_DRIVERS), "date intervals\n")
print(MONTHLY_DRIVERS)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n=== ALL INPUT FILES READY FOR RF WORKFLOW ===\n")
cat("✓ 1. TREE_SUMMER_125:", nrow(TREE_SUMMER_125), "trees\n")
cat("✓ 2. TREE_MONTHLY:", nrow(TREE_MONTHLY), "measurements\n")
cat("✓ 3. SOIL_YEAR:", nrow(SOIL_YEAR), "measurements\n")
cat("✓ 4. INVENTORY:", nrow(INVENTORY), "trees\n")
cat("✓ 5. PLOT_AREA:", round(PLOT_AREA, 0), "m²\n")
cat("✓ 6. MOISTURE_RASTER_DECEMBER: ", length(MOISTURE_RASTER_DECEMBER$x), "x", 
    length(MOISTURE_RASTER_DECEMBER$y), "grid\n")
cat("✓ 7. TAXONOMY:", nrow(TAXONOMY), "species with phylogenetic hierarchy\n")
cat("✓ 8. MONTHLY_DRIVERS:", nrow(MONTHLY_DRIVERS), "date intervals\n")




# Save all prepared datasets for the RF workflow
cat("\n=== SAVING ALL PREPARED DATASETS ===\n")

# Create a list of all datasets
rf_input_data <- list(
  tree_summer = TREE_SUMMER_125,
  tree_monthly = TREE_MONTHLY,
  soil_year = SOIL_YEAR,
  inventory = INVENTORY,
  plot_area = PLOT_AREA,
  moisture_raster = MOISTURE_RASTER_DECEMBER,
  taxonomy = TAXONOMY,
  monthly_drivers = MONTHLY_DRIVERS
)

# Save as RData for easy loading
save(rf_input_data, file = "rf_workflow_input_data.RData")

# Print final summary
cat("\n=== RF WORKFLOW INPUT DATA READY ===\n\n")

cat("1. TREE SUMMER (rigid chambers):\n")
cat("   - N =", nrow(TREE_SUMMER_125), "trees at 125cm height\n")
cat("   - Has VWC and soil temp\n")
cat("   - Date: July-August 2021\n\n")

cat("2. TREE MONTHLY (semirigid chambers):\n")
cat("   - N =", nrow(TREE_MONTHLY), "measurements\n")
cat("   - Species available for", sum(!is.na(TREE_MONTHLY$species)), "measurements\n")
cat("   - Date range: June 2020 - May 2021\n\n")

cat("3. SOIL YEAR:\n")
cat("   - N =", nrow(SOIL_YEAR), "measurements\n")
cat("   - Habitats:", paste(unique(SOIL_YEAR$habitat), collapse = ", "), "\n\n")

cat("4. INVENTORY:\n")
cat("   - N =", nrow(INVENTORY), "trees\n")
cat("   - Species:", length(unique(INVENTORY$species)), "\n")
cat("   - For scaling up from measured subset\n\n")

cat("5. PLOT AREA:\n")
cat("   -", round(PLOT_AREA, 0), "m² (", round(PLOT_AREA/10000, 2), "ha)\n\n")

cat("6. MOISTURE RASTER:\n")
cat("   - December 2020 interpolation\n")
cat("   - Grid:", length(MOISTURE_RASTER_DECEMBER$x), "x", length(MOISTURE_RASTER_DECEMBER$y), "\n\n")

cat("7. TAXONOMY:\n")
cat("   -", nrow(TAXONOMY), "species with phylogenetic hierarchy\n")
cat("   - Levels: species, genus, family, order, class\n\n")

cat("8. MONTHLY DRIVERS:\n")
cat("   -", nrow(MONTHLY_DRIVERS), "measurement periods\n")
cat("   - Air temp available for all\n")
cat("   - Soil temp available for summer only\n\n")





# Check environmental variables in summer intensive data
cat("=== SUMMER INTENSIVE ENVIRONMENTAL VARIABLES ===\n")
if(exists("summer_tree_data")) {
  env_cols_summer <- names(summer_tree_data)[grepl("temp|moisture|VWC|Tcham", names(summer_tree_data), ignore.case = TRUE)]
  cat("Environmental columns in summer data:\n")
  print(env_cols_summer)
  
  # Check for missing values
  for(col in env_cols_summer) {
    missing_pct <- round(100 * sum(is.na(summer_tree_data[[col]])) / nrow(summer_tree_data), 1)
    cat(col, "- Missing:", missing_pct, "%\n")
  }
}

# Check environmental variables in yearlong semirigid data
cat("\n=== YEARLONG ENVIRONMENTAL VARIABLES ===\n")
if(exists("semirigid_data")) {
  env_cols_year <- names(semirigid_data)[grepl("temp|moisture|VWC|Tcham", names(semirigid_data), ignore.case = TRUE)]
  cat("Environmental columns in yearlong data:\n")
  print(env_cols_year)
}

# Check species coverage
cat("=== SPECIES COVERAGE CHECK ===\n")
summer_species <- unique(TREE_SUMMER_125$species[!is.na(TREE_SUMMER_125$species)])
monthly_species <- unique(TREE_MONTHLY$species[!is.na(TREE_MONTHLY$species)])
inventory_species <- unique(INVENTORY$species[!is.na(INVENTORY$species)])

cat("Summer intensive species:", length(summer_species), "\n")
cat("Monthly monitoring species:", length(monthly_species), "\n") 
cat("Inventory species:", length(inventory_species), "\n")

all_measured_species <- unique(c(summer_species, monthly_species))
cat("Total species with flux measurements:", length(all_measured_species), "\n")
cat("Inventory species with flux data:", sum(inventory_species %in% all_measured_species), "out of", length(inventory_species), "\n")
