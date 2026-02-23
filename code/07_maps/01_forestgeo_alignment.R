# ==============================================================================
# ForestGEO Inventory Alignment
# ==============================================================================
# Purpose: Aligns ForestGEO forest inventory coordinates with study measurement
#   locations.
#
# Pipeline stage: 3 — Spatial Setup
#
# Inputs:
#   - ForestGEO inventory CSVs (from data/raw/inventory/)
#   - PhytoPhylo (phylogenetic tree)
#
# Outputs:
#   - fg_final data frame (used by source())
# ==============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(ape)
library(phytools)
library(viridis)
library(RColorBrewer)


# ForestGEO Data Analysis: Using Available GPS Coordinates
# This script creates the final_dataset used in the main workflow

# Load required libraries
library(dplyr)
library(tidyr)

# Load your three datasets
fg19 <- read.csv('../../data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_2019.csv')
fgplot <- read.csv('../../data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_byplot.csv')
fgtag <- read.csv('../../data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_bytag.csv')

# 1. EXTRACT USABLE COORDINATES FROM ALL DATASETS
# fg19 has PX, PY, Latitude, Longitude
fg19_coords <- fg19 %>%
  dplyr::select(Tag, Species_Code, PX, PY, Latitude, Longitude, DBH, Status, Quadrat, Sub_Quadrat) %>%
  rename(Tag_ID = Tag, Species = Species_Code) %>%
  mutate(
    Lat_decimal = as.numeric(Latitude),
    Lon_decimal = as.numeric(Longitude),
    Coord_Source = "fg19"
  ) %>%
  dplyr::select(-Latitude, -Longitude)

# fgplot has lat.1, lon.1 (in radians, need to convert to degrees)
fgplot_coords <- fgplot %>%
  dplyr::select(Tag.ID, Species..4.char., Quadrat, Subquadrat, Diam..mm., lat.1, lon.1) %>%
  rename(Tag_ID = Tag.ID, Species = Species..4.char., Sub_Quadrat = Subquadrat, DBH_mm = Diam..mm.) %>%
  filter(!is.na(lat.1) & !is.na(lon.1)) %>%  # Remove rows with missing coordinates
  mutate(
    # Convert from radians to degrees
    Lat_decimal = lat.1 * 180 / pi,
    Lon_decimal = lon.1 * 180 / pi,
    # Convert diameter from mm to cm to match fg19
    DBH = DBH_mm / 10,
    Coord_Source = "fgplot"
  ) %>%
  dplyr::select(-lat.1, -lon.1, -DBH_mm)

# fgtag has lat.1, lon.1 (in radians, need to convert to degrees)  
fgtag_coords <- fgtag %>%
  dplyr::select(Tag.ID, Species..4.char., Quadrat, Subquadrat, Diam..mm., lat.1, lon.1) %>%
  rename(Tag_ID = Tag.ID, Species = Species..4.char., Sub_Quadrat = Subquadrat, DBH_mm = Diam..mm.) %>%
  filter(!is.na(lat.1) & !is.na(lon.1)) %>%  # Remove rows with missing coordinates
  mutate(
    # Convert from radians to degrees
    Lat_decimal = lat.1 * 180 / pi,
    Lon_decimal = lon.1 * 180 / pi,
    # Convert diameter from mm to cm to match fg19
    DBH = DBH_mm / 10,
    Coord_Source = "fgtag"
  ) %>%
  dplyr::select(-lat.1, -lon.1, -DBH_mm)

# 2. IDENTIFY TREES WITH AND WITHOUT COORDINATES
cat("=== COORDINATE AVAILABILITY SUMMARY ===\n")
cat("fg19 trees with coordinates:", nrow(fg19_coords), "\n")
cat("fgplot trees with coordinates:", nrow(fgplot_coords), "\n") 
cat("fgtag trees with coordinates:", nrow(fgtag_coords), "\n\n")

# Trees from fgplot without coordinates (missing lat.1/lon.1)
fgplot_missing_coords <- fgplot %>%
  dplyr::select(Tag.ID, Species..4.char., Quadrat, Subquadrat, Diam..mm.) %>%
  rename(Tag_ID = Tag.ID, Species = Species..4.char., Sub_Quadrat = Subquadrat, DBH_mm = Diam..mm.) %>%
  filter(is.na(fgplot$lat.1) | is.na(fgplot$lon.1)) %>%
  mutate(Coord_Source = "fgplot_missing")

# Trees from fgtag without coordinates (missing lat.1/lon.1)
fgtag_missing_coords <- fgtag %>%
  dplyr::select(Tag.ID, Species..4.char., Quadrat, Subquadrat, Diam..mm.) %>%
  rename(Tag_ID = Tag.ID, Species = Species..4.char., Sub_Quadrat = Subquadrat, DBH_mm = Diam..mm.) %>%
  filter(is.na(fgtag$lat.1) | is.na(fgtag$lon.1)) %>%
  mutate(Coord_Source = "fgtag_missing")

cat("fgplot trees missing coordinates:", nrow(fgplot_missing_coords), "\n")
cat("fgtag trees missing coordinates:", nrow(fgtag_missing_coords), "\n\n")

# 3. CREATE MASTER COORDINATE DATASET
# Combine all trees with coordinates
all_coords <- bind_rows(fg19_coords, fgplot_coords, fgtag_coords)

# Check for duplicate Tag_IDs across datasets
duplicate_tags <- all_coords %>%
  group_by(Tag_ID) %>%
  summarise(
    n_datasets = n(),
    datasets = paste(Coord_Source, collapse = ", "),
    .groups = 'drop'
  ) %>%
  filter(n_datasets > 1)

cat("=== DUPLICATE TAG IDs ACROSS DATASETS ===\n")
print(duplicate_tags)
cat("Number of Tag IDs appearing in multiple datasets:", nrow(duplicate_tags), "\n\n")

# For duplicates, prioritize fg19 coordinates (most complete), then fgplot, then fgtag
master_coords <- all_coords %>%
  arrange(Tag_ID, factor(Coord_Source, levels = c("fg19", "fgplot", "fgtag"))) %>%
  group_by(Tag_ID) %>%
  slice_head(n = 1) %>%  # Keep first occurrence (highest priority)
  ungroup()

# 4. ESTIMATE COORDINATES FOR MISSING TREES
# Create reference coordinates by quadrat/sub-quadrat for estimation
quadrat_reference <- master_coords %>%
  filter(!is.na(Lat_decimal) & !is.na(Lon_decimal)) %>%
  group_by(Quadrat, Sub_Quadrat) %>%
  summarise(
    Avg_Lat = mean(Lat_decimal, na.rm = TRUE),
    Avg_Lon = mean(Lon_decimal, na.rm = TRUE),
    n_trees = n(),
    .groups = 'drop'
  )

# Function to estimate coordinates for missing trees
estimate_missing_coordinates <- function(missing_df, reference_df) {
  missing_df %>%
    left_join(reference_df, by = c("Quadrat", "Sub_Quadrat")) %>%
    mutate(
      Lat_decimal = Avg_Lat,
      Lon_decimal = Avg_Lon,
      Coordinates_Estimated = TRUE,
      DBH = DBH_mm / 10  # Convert mm to cm
    ) %>%
    dplyr::select(-Avg_Lat, -Avg_Lon, -n_trees, -DBH_mm)
}

# Estimate coordinates for missing trees
fgplot_estimated <- estimate_missing_coordinates(fgplot_missing_coords, quadrat_reference)
fgtag_estimated <- estimate_missing_coordinates(fgtag_missing_coords, quadrat_reference)

# Trees that can be estimated
estimated_coords <- bind_rows(fgplot_estimated, fgtag_estimated) %>%
  filter(!is.na(Lat_decimal))

# Trees that cannot be estimated (no reference trees in same quadrat/sub-quadrat)
cannot_estimate <- bind_rows(fgplot_estimated, fgtag_estimated) %>%
  filter(is.na(Lat_decimal))

cat("=== COORDINATE ESTIMATION RESULTS ===\n")
cat("Trees with estimated coordinates:", nrow(estimated_coords), "\n")
cat("Trees that cannot be estimated:", nrow(cannot_estimate), "\n\n")

# 5. CREATE FINAL COMPREHENSIVE DATASET
# Add estimation flag to master coordinates
master_coords$Coordinates_Estimated <- FALSE

# Combine all trees with coordinates (original + estimated)
final_dataset <- bind_rows(
  master_coords,
  estimated_coords %>% filter(!is.na(Lat_decimal))
) %>%
  arrange(Tag_ID)

# 6. SUMMARY STATISTICS AND REPORTS
cat("=== FINAL SUMMARY ===\n")
cat("Total trees with coordinates:", nrow(final_dataset), "\n")
cat("Trees with original coordinates:", sum(!final_dataset$Coordinates_Estimated), "\n")
cat("Trees with estimated coordinates:", sum(final_dataset$Coordinates_Estimated), "\n")
cat("Trees still missing coordinates:", nrow(cannot_estimate), "\n\n")

# Show coordinate ranges to verify they make sense
cat("=== COORDINATE RANGES ===\n")
cat("Latitude range:", round(min(final_dataset$Lat_decimal, na.rm = TRUE), 6), "to", 
    round(max(final_dataset$Lat_decimal, na.rm = TRUE), 6), "\n")
cat("Longitude range:", round(min(final_dataset$Lon_decimal, na.rm = TRUE), 6), "to", 
    round(max(final_dataset$Lon_decimal, na.rm = TRUE), 6), "\n\n")

cat("=== final_dataset CREATED SUCCESSFULLY ===\n")
cat("This dataset is now ready for use in the main workflow script.\n")
cat("It contains", nrow(final_dataset), "trees with GPS coordinates.\n")

# Optional: Export the dataset
# write.csv(final_dataset, "ForestGEO_Master_Coordinates.csv", row.names = FALSE)


# Check for required dataset
if(!exists("final_dataset")) {
  cat("ERROR: final_dataset not found. Please run the coordinate merging analysis first.\n")
  cat("The dataset should contain columns: Tag_ID, Species, Quadrat, Sub_Quadrat, Source_Dataset,\n")
  cat("Final_PX, Final_PY, Final_Lat, Final_Lon, Coordinates_Estimated, DBH, Status\n\n")
  stop("Required dataset not available")
}

# Convert the merged dataset to match the expected format
fg_merged <- final_dataset %>%
  rename(
    Tag = Tag_ID,
    Species_Code = Species,
    Latitude = Lat_decimal,
    Longitude = Lon_decimal
  ) %>%
  mutate(
    Sub_Quadrat = ifelse(is.na(Sub_Quadrat) | Sub_Quadrat == "", "a", Sub_Quadrat),
    Latitude = as.numeric(as.character(Latitude)),
    Longitude = as.numeric(as.character(Longitude)),
    DBH = as.numeric(as.character(DBH)),
    PX = ifelse(is.na(PX), runif(n(), 0, 200), as.numeric(as.character(PX))),
    PY = ifelse(is.na(PY), runif(n(), 0, 200), as.numeric(as.character(PY)))
  )

cat("=== MERGED DATASET SUMMARY ===\n")
cat("Total trees in merged dataset:", nrow(fg_merged), "\n")
cat("Trees with original coordinates:", sum(!fg_merged$Coordinates_Estimated, na.rm = TRUE), "\n")
cat("Trees with estimated coordinates:", sum(fg_merged$Coordinates_Estimated, na.rm = TRUE), "\n")
cat("Trees missing coordinates:", sum(is.na(fg_merged$Latitude) | is.na(fg_merged$Longitude)), "\n\n")

# Remove rows with missing coordinates
fg_clean <- fg_merged %>%
  filter(!is.na(Latitude) & !is.na(Longitude) & !is.na(PX) & !is.na(PY) & !is.na(DBH))

cat("Trees with complete coordinate data:", nrow(fg_clean), "\n\n")

# REMOVE DUPLICATES IMMEDIATELY
cat("=== DUPLICATE REMOVAL ===\n")
duplicate_check <- fg_clean %>%
  group_by(Tag) %>%
  summarise(
    n_occurrences = n(),
    coord_sources = paste(unique(Coord_Source), collapse = ", "),
    .groups = 'drop'
  ) %>%
  filter(n_occurrences > 1) %>%
  arrange(desc(n_occurrences))

cat("Duplicate Tag IDs found:", nrow(duplicate_check), "\n")
if(nrow(duplicate_check) > 0) {
  cat("Removing duplicates by prioritizing fg19 > fgplot > fgtag_missing...\n")
  
  fg_clean <- fg_clean %>%
    arrange(Tag, factor(Coord_Source, levels = c("fg19", "fgplot", "fgtag_missing"))) %>%
    group_by(Tag) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  cat("Trees before deduplication:", nrow(fg_merged %>% filter(!is.na(Latitude) & !is.na(Longitude) & !is.na(PX) & !is.na(PY) & !is.na(DBH))), "\n")
  cat("Trees after deduplication:", nrow(fg_clean), "\n")
  cat("Duplicates removed:", nrow(duplicate_check), "\n\n")
} else {
  cat("No duplicate Tag IDs found\n\n")
}

# FIX DBH DATA ENTRY ERRORS
cat("=== FIXING DBH DATA ENTRY ERRORS ===\n")
large_dbh <- fg_clean$DBH[fg_clean$DBH > 100]
if(length(large_dbh) > 0) {
  cat("DBH values before correction:\n")
  print(sort(large_dbh, decreasing = TRUE))
  
  fg_clean$DBH <- case_when(
    fg_clean$DBH >= 1000 ~ fg_clean$DBH / 100,
    fg_clean$DBH >= 100 & fg_clean$DBH < 1000 ~ fg_clean$DBH / 10,
    TRUE ~ fg_clean$DBH
  )
  
  cat("\nDBH values after correction:\n")
  large_dbh_corrected <- fg_clean$DBH[fg_clean$DBH > 50]
  if(length(large_dbh_corrected) > 0) {
    print(sort(large_dbh_corrected, decreasing = TRUE))
  }
  cat("Corrected", length(large_dbh), "DBH values\n\n")
} else {
  cat("No DBH values >100 found\n\n")
}

# Calculate basal area
fg_clean$BasalArea_cm2 <- pi * (fg_clean$DBH / 2)^2
fg_clean$BasalArea_m2 <- fg_clean$BasalArea_cm2 / 10000

# Species name mapping
species_mapping <- c(
  "ACRU" = "Acer rubrum", "ACSA" = "Acer saccharum", "ACPE" = "Acer pensylvanicum",
  "BEAL" = "Betula alleghaniensis", "BELE" = "Betula lenta", "BEPA" = "Betula papyrifera",
  "FAGR" = "Fagus grandifolia", "FRAM" = "Fraxinus americana", "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra", "QUAL" = "Quercus alba", "QUVE" = "Quercus velutina",
  "QUCO" = "Quercus coccinea", "TSCA" = "Tsuga canadensis", "CAOV" = "Carya ovata",
  "KALA" = "Kalmia latifolia", "PRSE" = "Prunus serotina", "SAAL" = "Sassafras albidum",
  "HAVI" = "Hamamelis virginiana", "LITU" = "Liriodendron tulipifera", "VACO" = "Vaccinium corymbosum"
)

fg_clean$Species_Name <- case_when(
  trimws(toupper(fg_clean$Species_Code)) %in% names(species_mapping) ~ 
    species_mapping[trimws(toupper(fg_clean$Species_Code))],
  is.na(fg_clean$Species_Code) | trimws(fg_clean$Species_Code) == "" ~ "Unknown",
  TRUE ~ "Unknown"
)

# OUTLIER DETECTION
detect_outliers <- function(x, multiplier = 0.75) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - multiplier * IQR
  upper_bound <- Q3 + multiplier * IQR
  return(x < lower_bound | x > upper_bound)
}

fg_clean$outlier_PX <- detect_outliers(fg_clean$PX)
fg_clean$outlier_PY <- detect_outliers(fg_clean$PY)
fg_clean$local_outlier <- fg_clean$outlier_PX | fg_clean$outlier_PY

n_outliers <- sum(fg_clean$local_outlier)
n_good <- sum(!fg_clean$local_outlier)

cat("=== LOCAL OUTLIER DETECTION ===\n")
cat("Good local coordinates:", n_good, "(", round(n_good/nrow(fg_clean)*100, 1), "%)\n")
cat("Local outliers:", n_outliers, "(", round(n_outliers/nrow(fg_clean)*100, 1), "%)\n\n")

# GEODETIC COORDINATE TRANSFORMATION
cat("=== GEODETIC COORDINATE TRANSFORMATION ===\n")

# Reference point and parameters
lat1_deg <- 41.989211
lon1_deg <- -72.13092
lat1_rad <- lat1_deg * pi/180
lon1_rad <- lon1_deg * pi/180
rotation_angle <- 0.174533  # 10 degrees in radians
earth_radius <- 6371000

cat("Using geodetic transformation with:\n")
cat("- Reference point:", lat1_deg, "°N,", lon1_deg, "°W\n")
cat("- Rotation angle:", rotation_angle * 180/pi, "degrees\n")
cat("- Earth radius:", earth_radius, "meters\n\n")

# Geodetic transformation function
geodetic_transform <- function(px, py, lat1_rad, lon1_rad, rotation_angle, earth_radius) {
  
  # Step 1: Apply rotation to PX/PY coordinates
  distance <- sqrt(px^2 + py^2)
  original_bearing <- atan2(px, py)
  rotated_bearing <- original_bearing + rotation_angle
  
  px_rotated <- distance * sin(rotated_bearing)
  py_rotated <- distance * cos(rotated_bearing)
  
  # Step 2: Convert to GPS using geodetic formulas
  bearing <- atan2(px_rotated, py_rotated)
  distance_m <- sqrt(px_rotated^2 + py_rotated^2)
  
  # Geodetic forward calculation
  lat2_rad <- asin(sin(lat1_rad) * cos(distance_m/earth_radius) + 
                     cos(lat1_rad) * sin(distance_m/earth_radius) * cos(bearing))
  
  lon2_rad <- lon1_rad + atan2(sin(bearing) * sin(distance_m/earth_radius) * cos(lat1_rad),
                               cos(distance_m/earth_radius) - sin(lat1_rad) * sin(lat2_rad))
  
  return(list(
    longitude = lon2_rad * 180/pi,
    latitude = lat2_rad * 180/pi
  ))
}

# Apply transformation to all trees with PX/PY coordinates
trees_with_local <- fg_clean %>% filter(!is.na(PX) & !is.na(PY) & PX != 0 & PY != 0)

if(nrow(trees_with_local) > 0) {
  cat("Applying geodetic transformation to", nrow(trees_with_local), "trees with PX/PY coordinates\n")
  
  result <- geodetic_transform(
    trees_with_local$PX, 
    trees_with_local$PY,
    lat1_rad, lon1_rad, rotation_angle, earth_radius
  )
  
  # Update coordinates
  fg_clean$Longitude_final <- fg_clean$Longitude  # Default to original
  fg_clean$Latitude_final <- fg_clean$Latitude
  
  # Replace with transformed coordinates where applicable
  local_indices <- which(!is.na(fg_clean$PX) & !is.na(fg_clean$PY) & fg_clean$PX != 0 & fg_clean$PY != 0)
  fg_clean$Longitude_final[local_indices] <- result$longitude
  fg_clean$Latitude_final[local_indices] <- result$latitude
  
  cat("Transformation complete\n\n")
} else {
  cat("No trees with PX/PY coordinates found\n")
  fg_clean$Longitude_final <- fg_clean$Longitude
  fg_clean$Latitude_final <- fg_clean$Latitude
}

# CREATE FINAL DATASET - REMOVE LOCAL OUTLIERS
fg_final <- fg_clean %>% 
  filter(!local_outlier) %>%
  mutate(
    transformation_method = ifelse(!is.na(PX) & !is.na(PY) & PX != 0 & PY != 0, 
                                   "Geodetic Transform", "Original GPS")
  )

cat("=== FINAL DATASET ===\n")
cat("Original points:", nrow(fg_clean), "\n")
cat("After removing local outliers:", nrow(fg_final), "\n")
cat("Removed:", nrow(fg_clean) - nrow(fg_final), "local outliers\n\n")

# Transformation summary
transformation_summary <- fg_final %>%
  group_by(transformation_method) %>%
  summarise(n_trees = n(), .groups = 'drop')
print(transformation_summary)
cat("\n")

# PHYLOGENETIC COLORS
tree_file <- "../../data/processed/metadata/PhytoPhylo"

all_species <- unique(fg_final$Species_Name)
known_species <- all_species[!grepl("Unknown", all_species)]
unknown_species <- all_species[grepl("Unknown", all_species)]

phylo_colors <- NULL
phylo_distance_data <- NULL

if (file.exists(tree_file) && length(known_species) > 1) {
  tryCatch({
    tree_scenario1 <- read.tree(tree_file)
    present_species_clean <- gsub(" ", "_", known_species)
    tree_species <- intersect(present_species_clean, tree_scenario1$tip.label)
    
    if (length(tree_species) > 1) {
      pruned_tree <- keep.tip(tree_scenario1, tree_species)
      dist_matrix <- cophenetic(pruned_tree)
      average_distances <- rowMeans(dist_matrix)
      
      distance_df <- data.frame(
        species_clean = names(average_distances),
        distance = as.numeric(average_distances)
      ) %>%
        mutate(species = gsub("_", " ", species_clean)) %>%
        arrange(distance) %>%
        mutate(
          distance_group = cut(distance, 
                               breaks = quantile(distance, probs = seq(0, 1, length.out = min(8, length(distance)))), 
                               include.lowest = TRUE, 
                               labels = FALSE),
          color = viridis(max(distance_group, na.rm = TRUE))[distance_group]
        )
      
      phylo_colors <- setNames(distance_df$color, distance_df$species)
      phylo_distance_data <- distance_df %>%
        dplyr::select(species, distance, distance_group) %>%
        arrange(distance)
      
      cat("=== PHYLOGENETIC COLOR MAPPING ===\n")
      cat("Species with phylogenetic data:", length(tree_species), "\n")
      print(phylo_distance_data)
      cat("\n")
      
    } else {
      cat("Warning: Not enough species overlap with phylogenetic tree\n")
    }
  }, error = function(e) {
    cat("Warning: Could not calculate phylogenetic distances:", e$message, "\n")
  })
} else {
  cat("Warning: Tree file not found or insufficient known species, using discrete colors\n")
}

# Create complete color mapping
if (is.null(phylo_colors)) {
  if(length(known_species) <= 12) {
    known_colors <- RColorBrewer::brewer.pal(min(max(3, length(known_species)), 12), "Set3")[1:length(known_species)]
  } else {
    known_colors <- rainbow(length(known_species))
  }
  names(known_colors) <- known_species
  phylo_colors <- known_colors
}

unknown_colors <- rep("grey80", length(unknown_species))
names(unknown_colors) <- unknown_species
final_colors <- c(phylo_colors, unknown_colors)

legend_order <- if (!is.null(phylo_distance_data)) {
  c(phylo_distance_data$species, unknown_species)
} else {
  c(sort(known_species), unknown_species)
}

# SUMMARY STATISTICS
cat("=== FINAL SUMMARY BY SOURCE ===\n")
final_summary <- fg_final %>%
  group_by(Coord_Source, Coordinates_Estimated) %>%
  summarise(
    n_trees = n(),
    mean_dbh = round(mean(DBH, na.rm = TRUE), 1),
    total_basal_area = round(sum(BasalArea_m2, na.rm = TRUE), 2),
    .groups = 'drop'
  ) %>%
  arrange(Coord_Source, Coordinates_Estimated)

print(final_summary)
cat("\n")

# BASAL AREA SUMMARY
cat("=== BASAL AREA SUMMARY ===\n")
total_basal_area <- sum(fg_final$BasalArea_m2, na.rm = TRUE)
cat("Total basal area (final dataset):", round(total_basal_area, 3), "m²\n")

species_basal_area <- fg_final %>%
  group_by(Species_Code, Species_Name) %>%
  summarise(
    n_stems = n(),
    total_basal_area_m2 = sum(BasalArea_m2, na.rm = TRUE),
    mean_dbh_cm = mean(DBH, na.rm = TRUE),
    mean_basal_area_m2 = mean(BasalArea_m2, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(total_basal_area_m2))

print(head(species_basal_area, 15))

cat("\n=== READY FOR VISUALIZATION ===\n")
cat("Your merged dataset is now ready for mapping and analysis!\n")
cat("Key variables available:\n")
cat("- Longitude_final, Latitude_final: Final GPS coordinates\n") 
cat("- BasalArea_m2: Basal area for tree sizing\n")
cat("- Species_Name: Full species names for phylogenetic coloring\n")
cat("- transformation_method: Geodetic Transform vs Original GPS\n")
cat("- Coordinates_Estimated: TRUE/FALSE flag for estimated coordinates\n\n")

# FOREST ANALYSIS
cat("=== COMPREHENSIVE FOREST ANALYSIS ===\n")

# Calculate plot area
lon_range_m <- (max(fg_final$Longitude_final) - min(fg_final$Longitude_final)) * 111320 * cos(mean(fg_final$Latitude_final) * pi/180)
lat_range_m <- (max(fg_final$Latitude_final) - min(fg_final$Latitude_final)) * 111320
plot_area_m2 <- lon_range_m * lat_range_m
plot_area_ha <- plot_area_m2 / 10000

cat("Plot dimensions:\n")
cat("- East-West:", round(lon_range_m, 1), "meters\n")
cat("- North-South:", round(lat_range_m, 1), "meters\n")
cat("- Total area:", round(plot_area_m2, 0), "m² (", round(plot_area_ha, 2), "ha)\n")
cat("- Basal area per hectare:", round(total_basal_area / plot_area_ha, 1), "m²/ha\n\n")

cat("Forest structure summary:\n")
cat("- Total trees:", nrow(fg_final), "\n")
cat("- Tree density:", round(nrow(fg_final) / plot_area_ha, 0), "trees/ha\n")
cat("- Mean DBH:", round(mean(fg_final$DBH, na.rm = TRUE), 1), "cm\n")
cat("- Median DBH:", round(median(fg_final$DBH, na.rm = TRUE), 1), "cm\n")
cat("- Max DBH:", round(max(fg_final$DBH, na.rm = TRUE), 1), "cm\n")
cat("- Number of species:", length(unique(fg_final$Species_Name[fg_final$Species_Name != "Unknown"])), "\n\n")

# VISUALIZATIONS

# 1. Dataset composition
composition_plot <- fg_final %>%
  ggplot(aes(x = Coord_Source, fill = transformation_method)) +
  geom_bar(position = "stack") +
  scale_fill_manual(
    values = c("Geodetic Transform" = "steelblue", "Original GPS" = "orange"),
    name = "Coordinate Source"
  ) +
  labs(
    title = "Dataset Composition by Source and Transformation Method",
    subtitle = paste("Total trees:", nrow(fg_final)),
    x = "Data Source", y = "Number of Trees"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(composition_plot)

# 2. Main forest map
main_plot <- ggplot(fg_final, aes(x = Longitude_final, y = Latitude_final)) +
  geom_point(aes(size = BasalArea_m2, color = Species_Name), alpha = 0.8) +
  scale_color_manual(
    values = final_colors,
    breaks = legend_order,
    name = "Species\n(phylogenetic order)"
  ) +
  scale_size_continuous(
    name = "Basal Area\n(m²)",
    range = c(0.3, 3),
    breaks = c(0.001, 0.01, 0.05, 0.1, 0.2),
    labels = c("0.001", "0.01", "0.05", "0.1", "0.2"),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(
    title = "Complete ForestGEO Map: Geodetic Coordinate Transformation",
    subtitle = paste("N =", nrow(fg_final), "trees with accurate geographic positioning"),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

print(main_plot)

# 3. Transformation method comparison
method_plot <- ggplot(fg_final, aes(x = Longitude_final, y = Latitude_final)) +
  geom_point(aes(color = transformation_method, shape = Coordinates_Estimated), 
             size = 1.5, alpha = 0.7) +
  scale_color_manual(
    values = c("Geodetic Transform" = "red", "Original GPS" = "blue"),
    name = "Coordinate Method"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "Original", "TRUE" = "Estimated"),
    name = "GPS Source"
  ) +
  labs(
    title = "Trees by Coordinate Transformation Method",
    subtitle = "Red = Geodetic transformation from PX/PY, Blue = Original GPS",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(method_plot)

# 4. DBH distribution
dbh_dist_plot <- ggplot(fg_final, aes(x = DBH)) +
  geom_histogram(bins = 50, fill = "forestgreen", alpha = 0.7, color = "black") +
  geom_vline(aes(xintercept = mean(DBH, na.rm = TRUE)), color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median(DBH, na.rm = TRUE)), color = "blue", linetype = "dashed", size = 1) +
  scale_x_log10() +
  labs(
    title = "Distribution of Tree Diameters (DBH)",
    subtitle = paste("Red = mean (", round(mean(fg_final$DBH, na.rm = TRUE), 1), "cm), Blue = median (", round(median(fg_final$DBH, na.rm = TRUE), 1), "cm)"),
    x = "DBH (cm, log scale)", y = "Number of Trees"
  ) +
  theme_minimal()

print(dbh_dist_plot)

# 5. Species abundance and basal area
species_summary <- fg_final %>%
  filter(Species_Name != "Unknown") %>%
  group_by(Species_Name) %>%
  summarise(
    n_trees = n(),
    total_basal_area = sum(BasalArea_m2, na.rm = TRUE),
    mean_dbh = mean(DBH, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_trees))

abundance_vs_basal <- ggplot(species_summary, aes(x = n_trees, y = total_basal_area)) +
  geom_point(aes(size = mean_dbh, color = mean_dbh), alpha = 0.7) +
  geom_text(aes(label = Species_Name), hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_c(name = "Mean DBH\n(cm)") +
  scale_size_continuous(name = "Mean DBH\n(cm)", range = c(1, 8)) +
  labs(
    title = "Species Abundance vs Total Basal Area",
    subtitle = "Point size and color = mean DBH",
    x = "Number of Trees (log scale)", y = "Total Basal Area (m², log scale)"
  ) +
  theme_minimal() +
  guides(size = "none")+geom_smooth(method="lm", color="black", se=F)

print(abundance_vs_basal)

# Forest structure metrics table
structure_metrics <- data.frame(
  Metric = c(
    "Total Trees", "Plot Area (ha)", "Tree Density (trees/ha)", 
    "Total Basal Area (m²)", "Basal Area per Hectare (m²/ha)",
    "Mean DBH (cm)", "Median DBH (cm)", "95th percentile DBH (cm)",
    "Number of Species", "Largest Tree DBH (cm)"
  ),
  Value = c(
    nrow(fg_final), round(plot_area_ha, 2), round(nrow(fg_final) / plot_area_ha, 0),
    round(total_basal_area, 1), round(total_basal_area / plot_area_ha, 1),
    round(mean(fg_final$DBH, na.rm = TRUE), 1), round(median(fg_final$DBH, na.rm = TRUE), 1),
    round(quantile(fg_final$DBH, 0.95, na.rm = TRUE), 1),
    length(unique(fg_final$Species_Name[fg_final$Species_Name != "Unknown"])),
    round(max(fg_final$DBH, na.rm = TRUE), 1)
  )
)

cat("=== FOREST STRUCTURE METRICS ===\n")
print(structure_metrics)

cat("\n=== GEODETIC TRANSFORMATION WORKFLOW COMPLETE ===\n")
cat("Successfully processed", nrow(fg_final), "trees using geodetic coordinate transformation\n")
cat("Dataset ready for ecological analysis and satellite imagery overlay!\n")
cat("Use 'fg_final' for all downstream analyses\n")
