fg19<-read.csv("/Users/jongewirtzman/Downloads/ForestGEO_data2021UPDATE_6_21_DW_2019.csv")
fgplot<-read.csv("/Users/jongewirtzman/Downloads/ForestGEO_data2021UPDATE_6_21_DW_bytag.csv")
fgtag<-read.csv("/Users/jongewirtzman/Downloads/ForestGEO_data2021UPDATE_6_21_DW_byplot.csv")

str(fg19)
head(fg19)

str(fgplot)
head(fgplot)

str(fgtag)
head(fgtag)


# ForestGEO Data Analysis: Using Available GPS Coordinates

# Load required libraries
library(dplyr)
library(tidyr)

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

# Display sample of trees with estimated coordinates
cat("=== SAMPLE OF TREES WITH ESTIMATED COORDINATES ===\n")
print(head(estimated_coords %>% dplyr::select(Tag_ID, Species, Quadrat, Sub_Quadrat, Lat_decimal, Lon_decimal, Coord_Source), 10))

cat("\n=== TREES THAT CANNOT BE ESTIMATED ===\n")
print(head(cannot_estimate %>% dplyr::select(Tag_ID, Species, Quadrat, Sub_Quadrat, Coord_Source), 10))

# Optional: Export results
# write.csv(final_dataset, "ForestGEO_Master_Coordinates.csv", row.names = FALSE)
# write.csv(estimated_coords, "Trees_with_Estimated_Coordinates.csv", row.names = FALSE)
# write.csv(cannot_estimate, "Trees_Cannot_Estimate_Coordinates.csv", row.names = FALSE)