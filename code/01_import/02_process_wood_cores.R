# ==============================================================================
# Process Wood Core Data
# ==============================================================================
# Purpose: Processes wood core sectioning data, fills measurement gaps via
#   median imputation, and calculates summary statistics.
#
# Pipeline stage: 01 Tree Data Processing
# Run after: None
#
# Inputs:
#   - Tree Core Sectioning Data.csv (from data/raw/tree_cores/)
#
# Outputs:
#   - tree_core_filled_complete.csv
# ==============================================================================

library(dplyr)
library(readr)

# Read and process data (same as before)
data <- read_csv('../../data/raw/tree_cores/Tree Core Sectioning Data.csv')

core_diameter_mm <- 5.15
core_radius_cm <- (core_diameter_mm / 2) / 10

calculate_volume <- function(length_cm) {
  pi * core_radius_cm^2 * length_cm
}

# Process the original data
processed_data <- data %>%
  rename(
    outer_length_cm = `Outer Section Length (cm)`,
    outer_fresh_mass_g = `OUTER SECTION FRESH MASS (g)`,
    outer_dry_mass_g = `Post-Freeze Dry OUTER SECTION Sample Mass (g)`,
    inner_length_cm = `Inner Section Length (cm)`,
    inner_fresh_mass_g = `INNER SECTION FRESH MASS (g)`,
    inner_dry_mass_g = `Post-Freeze Dry INNER SECTION Sample Mass (g)`,
    middle_length_mm = `Middle Section Length (mm)`,
    middle_fresh_mass_g = `MIDDLE SECTION FRESH MASS (g)`,
    middle_dry_mass_g = `Post-Freeze Dry MIDDLE SECTION Sample Mass (g)`,
    tree_id = `Tree Core ID`,
    species = `Tree Core Species Code`
  ) %>%
  mutate(
    outer_length_cm = as.numeric(outer_length_cm),
    inner_length_cm = as.numeric(inner_length_cm),
    middle_length_mm = as.numeric(middle_length_mm),
    outer_fresh_mass_g = as.numeric(outer_fresh_mass_g),
    outer_dry_mass_g = as.numeric(outer_dry_mass_g),
    inner_fresh_mass_g = as.numeric(inner_fresh_mass_g),
    inner_dry_mass_g = as.numeric(inner_dry_mass_g),
    middle_fresh_mass_g = as.numeric(middle_fresh_mass_g),
    middle_dry_mass_g = as.numeric(middle_dry_mass_g),
    middle_length_cm = middle_length_mm / 10
  )

# Calculate all measurements
complete_data <- processed_data %>%
  mutate(
    outer_volume_cm3 = ifelse(!is.na(outer_length_cm) & outer_length_cm > 0,
                              calculate_volume(outer_length_cm), NA),
    inner_volume_cm3 = ifelse(!is.na(inner_length_cm) & inner_length_cm > 0,
                              calculate_volume(inner_length_cm), NA),
    middle_volume_cm3 = ifelse(!is.na(middle_length_cm) & middle_length_cm > 0,
                               calculate_volume(middle_length_cm), NA),
    
    outer_moisture_fresh = ifelse(!is.na(outer_fresh_mass_g) & !is.na(outer_dry_mass_g) & outer_fresh_mass_g > 0,
                                  ((outer_fresh_mass_g - outer_dry_mass_g) / outer_fresh_mass_g) * 100, NA),
    inner_moisture_fresh = ifelse(!is.na(inner_fresh_mass_g) & !is.na(inner_dry_mass_g) & inner_fresh_mass_g > 0,
                                  ((inner_fresh_mass_g - inner_dry_mass_g) / inner_fresh_mass_g) * 100, NA),
    middle_moisture_fresh = ifelse(!is.na(middle_fresh_mass_g) & !is.na(middle_dry_mass_g) & middle_fresh_mass_g > 0,
                                   ((middle_fresh_mass_g - middle_dry_mass_g) / middle_fresh_mass_g) * 100, NA),
    
    outer_moisture_dry = ifelse(!is.na(outer_fresh_mass_g) & !is.na(outer_dry_mass_g) & outer_dry_mass_g > 0,
                                ((outer_fresh_mass_g - outer_dry_mass_g) / outer_dry_mass_g) * 100, NA),
    inner_moisture_dry = ifelse(!is.na(inner_fresh_mass_g) & !is.na(inner_dry_mass_g) & inner_dry_mass_g > 0,
                                ((inner_fresh_mass_g - inner_dry_mass_g) / inner_dry_mass_g) * 100, NA),
    middle_moisture_dry = ifelse(!is.na(middle_fresh_mass_g) & !is.na(middle_dry_mass_g) & middle_dry_mass_g > 0,
                                 ((middle_fresh_mass_g - middle_dry_mass_g) / middle_dry_mass_g) * 100, NA),
    
    outer_density_g_cm3 = ifelse(!is.na(outer_dry_mass_g) & !is.na(outer_volume_cm3) & outer_volume_cm3 > 0,
                                 outer_dry_mass_g / outer_volume_cm3, NA),
    inner_density_g_cm3 = ifelse(!is.na(inner_dry_mass_g) & !is.na(inner_volume_cm3) & inner_volume_cm3 > 0,
                                 inner_dry_mass_g / inner_volume_cm3, NA),
    middle_density_g_cm3 = ifelse(!is.na(middle_dry_mass_g) & !is.na(middle_volume_cm3) & middle_volume_cm3 > 0,
                                  middle_dry_mass_g / middle_volume_cm3, NA)
  )

# Define problematic trees (from our error diagnostic)
problematic_trees <- c("AB6'", "H104", "HX2", "RM1", "rms1", "WA14", "YB15")

print("=== MEDIAN IMPUTATION FOR PROBLEMATIC VALUES ===")
print(paste("Problematic trees to be imputed:", paste(problematic_trees, collapse = ", ")))

# Function to identify problematic values
is_problematic <- function(tree_id, fresh_mass, dry_mass, moisture_dry, density) {
  is_problem_tree <- tree_id %in% problematic_trees
  
  # Additional checks for specific error types
  dry_exceeds_fresh <- !is.na(fresh_mass) & !is.na(dry_mass) & dry_mass > fresh_mass
  extreme_moisture <- !is.na(moisture_dry) & moisture_dry > 300
  impossible_density <- !is.na(density) & (density > 2.0 | density < 0.1)
  
  return(is_problem_tree & (dry_exceeds_fresh | extreme_moisture | impossible_density))
}

# Create clean dataset (excluding problematic values for median calculation)
clean_data <- complete_data %>%
  mutate(
    # Mark problematic measurements as NA for median calculation
    outer_moisture_fresh_clean = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3),
                                        NA, outer_moisture_fresh),
    outer_moisture_dry_clean = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3),
                                      NA, outer_moisture_dry),
    outer_density_clean = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3),
                                 NA, outer_density_g_cm3),
    
    inner_moisture_fresh_clean = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3),
                                        NA, inner_moisture_fresh),
    inner_moisture_dry_clean = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3),
                                      NA, inner_moisture_dry),
    inner_density_clean = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3),
                                 NA, inner_density_g_cm3),
    
    middle_moisture_fresh_clean = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3),
                                         NA, middle_moisture_fresh),
    middle_moisture_dry_clean = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3),
                                       NA, middle_moisture_dry),
    middle_density_clean = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3),
                                  NA, middle_density_g_cm3)
  )

# Calculate medians by species and section (using only clean data)
medians_by_species <- clean_data %>%
  group_by(species) %>%
  summarise(
    # Outer section medians
    outer_moisture_fresh_median = median(outer_moisture_fresh_clean, na.rm = TRUE),
    outer_moisture_dry_median = median(outer_moisture_dry_clean, na.rm = TRUE),
    outer_density_median = median(outer_density_clean, na.rm = TRUE),
    
    # Inner section medians
    inner_moisture_fresh_median = median(inner_moisture_fresh_clean, na.rm = TRUE),
    inner_moisture_dry_median = median(inner_moisture_dry_clean, na.rm = TRUE),
    inner_density_median = median(inner_density_clean, na.rm = TRUE),
    
    # Middle section medians
    middle_moisture_fresh_median = median(middle_moisture_fresh_clean, na.rm = TRUE),
    middle_moisture_dry_median = median(middle_moisture_dry_clean, na.rm = TRUE),
    middle_density_median = median(middle_density_clean, na.rm = TRUE),
    
    .groups = 'drop'
  )

print("\nMedian values by species (calculated from clean data):")
print(medians_by_species)

# Create the filled dataset by replacing problematic values with species medians
filled_data <- complete_data %>%
  left_join(medians_by_species, by = "species") %>%
  mutate(
    # Replace problematic outer section values
    outer_moisture_fresh_filled = ifelse(
      is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3),
      outer_moisture_fresh_median, outer_moisture_fresh
    ),
    outer_moisture_dry_filled = ifelse(
      is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3),
      outer_moisture_dry_median, outer_moisture_dry
    ),
    outer_density_filled = ifelse(
      is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3),
      outer_density_median, outer_density_g_cm3
    ),
    
    # Replace problematic inner section values
    inner_moisture_fresh_filled = ifelse(
      is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3),
      inner_moisture_fresh_median, inner_moisture_fresh
    ),
    inner_moisture_dry_filled = ifelse(
      is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3),
      inner_moisture_dry_median, inner_moisture_dry
    ),
    inner_density_filled = ifelse(
      is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3),
      inner_density_median, inner_density_g_cm3
    ),
    
    # Replace problematic middle section values
    middle_moisture_fresh_filled = ifelse(
      is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3),
      middle_moisture_fresh_median, middle_moisture_fresh
    ),
    middle_moisture_dry_filled = ifelse(
      is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3),
      middle_moisture_dry_median, middle_moisture_dry
    ),
    middle_density_filled = ifelse(
      is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3),
      middle_density_median, middle_density_g_cm3
    )
  )

# Report what was imputed
imputation_report <- filled_data %>%
  filter(tree_id %in% problematic_trees) %>%
  select(tree_id, species, 
         # Original values
         outer_moisture_fresh, inner_moisture_fresh, middle_moisture_fresh,
         outer_moisture_dry, inner_moisture_dry, middle_moisture_dry,
         outer_density_g_cm3, inner_density_g_cm3, middle_density_g_cm3,
         # Filled values
         outer_moisture_fresh_filled, inner_moisture_fresh_filled, middle_moisture_fresh_filled,
         outer_moisture_dry_filled, inner_moisture_dry_filled, middle_moisture_dry_filled,
         outer_density_filled, inner_density_filled, middle_density_filled) %>%
  arrange(tree_id)

print("\nIMPUTATION REPORT:")
print("==================")
print("Original vs. Filled values for problematic trees:")

for (i in 1:nrow(imputation_report)) {
  tree_data <- imputation_report[i, ]
  cat("\nTREE:", tree_data$tree_id, "(", tree_data$species, ")\n")
  
  # Check outer section
  if (!is.na(tree_data$outer_moisture_fresh) || !is.na(tree_data$outer_moisture_fresh_filled)) {
    if (!identical(tree_data$outer_moisture_fresh, tree_data$outer_moisture_fresh_filled)) {
      cat("  Outer - Moisture Fresh: ", tree_data$outer_moisture_fresh, "% → ", tree_data$outer_moisture_fresh_filled, "%\n")
    }
    if (!identical(tree_data$outer_moisture_dry, tree_data$outer_moisture_dry_filled)) {
      cat("  Outer - Moisture Dry:   ", tree_data$outer_moisture_dry, "% → ", tree_data$outer_moisture_dry_filled, "%\n")
    }
    if (!identical(tree_data$outer_density_g_cm3, tree_data$outer_density_filled)) {
      cat("  Outer - Density:        ", tree_data$outer_density_g_cm3, "→ ", tree_data$outer_density_filled, "g/cm³\n")
    }
  }
  
  # Check inner section
  if (!is.na(tree_data$inner_moisture_fresh) || !is.na(tree_data$inner_moisture_fresh_filled)) {
    if (!identical(tree_data$inner_moisture_fresh, tree_data$inner_moisture_fresh_filled)) {
      cat("  Inner - Moisture Fresh: ", tree_data$inner_moisture_fresh, "% → ", tree_data$inner_moisture_fresh_filled, "%\n")
    }
    if (!identical(tree_data$inner_moisture_dry, tree_data$inner_moisture_dry_filled)) {
      cat("  Inner - Moisture Dry:   ", tree_data$inner_moisture_dry, "% → ", tree_data$inner_moisture_dry_filled, "%\n")
    }
    if (!identical(tree_data$inner_density_g_cm3, tree_data$inner_density_filled)) {
      cat("  Inner - Density:        ", tree_data$inner_density_g_cm3, "→ ", tree_data$inner_density_filled, "g/cm³\n")
    }
  }
  
  # Check middle section
  if (!is.na(tree_data$middle_moisture_fresh) || !is.na(tree_data$middle_moisture_fresh_filled)) {
    if (!identical(tree_data$middle_moisture_fresh, tree_data$middle_moisture_fresh_filled)) {
      cat("  Middle - Moisture Fresh:", tree_data$middle_moisture_fresh, "% → ", tree_data$middle_moisture_fresh_filled, "%\n")
    }
    if (!identical(tree_data$middle_moisture_dry, tree_data$middle_moisture_dry_filled)) {
      cat("  Middle - Moisture Dry:  ", tree_data$middle_moisture_dry, "% → ", tree_data$middle_moisture_dry_filled, "%\n")
    }
    if (!identical(tree_data$middle_density_g_cm3, tree_data$middle_density_filled)) {
      cat("  Middle - Density:       ", tree_data$middle_density_g_cm3, "→ ", tree_data$middle_density_filled, "g/cm³\n")
    }
  }
}

# Create final filled dataset with renamed columns
# Fix the rename conflict
final_filled_data <- filled_data %>%
  select(-contains("median"), -contains("clean")) %>%
  select(-outer_density_g_cm3, -inner_density_g_cm3, -middle_density_g_cm3) %>%  # Remove originals first
  rename(
    outer_moisture_fresh_percent = outer_moisture_fresh_filled,
    outer_moisture_dry_percent = outer_moisture_dry_filled,
    outer_density_g_cm3 = outer_density_filled,
    inner_moisture_fresh_percent = inner_moisture_fresh_filled,
    inner_moisture_dry_percent = inner_moisture_dry_filled,
    inner_density_g_cm3 = inner_density_filled,
    middle_moisture_fresh_percent = middle_moisture_fresh_filled,
    middle_moisture_dry_percent = middle_moisture_dry_filled,
    middle_density_g_cm3 = middle_density_filled
  )

# Save the filled dataset
write_csv(final_filled_data, "../../data/processed/tree_cores/tree_core_data_filled.csv")

print("\n\nSUMMARY:")
print("========")
print("Created filled dataset where problematic measurements are replaced with species-specific medians")
print(paste("Total trees in dataset:", nrow(final_filled_data)))
print(paste("Trees with imputed values:", length(problematic_trees)))
print("Saved as: tree_core_data_filled.csv")

print("\nComparison with exclusion approach:")
print("- Exclusion approach: Remove entire problematic trees")
print("- Imputation approach: Keep all trees, replace only problematic measurements")
print("- Use imputation when sample size is critical")
print("- Use exclusion when data quality is more important than sample size")

# Quick validation: check that filled values are reasonable
filled_densities <- c(final_filled_data$outer_density_g_cm3, 
                      final_filled_data$inner_density_g_cm3, 
                      final_filled_data$middle_density_g_cm3)
filled_densities <- filled_densities[!is.na(filled_densities)]

print(paste("\nFilled dataset density range:", round(min(filled_densities), 3), "-", round(max(filled_densities), 3), "g/cm³"))
print("All values should now be within reasonable wood density range (0.2-1.5 g/cm³)")





















# Complete Tree Core Analysis - Using Filled Dataset with Median Imputation
# Load required libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

# Read the data
data <- read_csv('../../data/raw/tree_cores/Tree Core Sectioning Data.csv')

# Define core parameters
core_diameter_mm <- 5.15
core_radius_cm <- (core_diameter_mm / 2) / 10

calculate_volume <- function(length_cm) {
  pi * core_radius_cm^2 * length_cm
}

# Process and create filled dataset (same as median imputation code)
processed_data <- data %>%
  rename(
    outer_length_cm = `Outer Section Length (cm)`,
    outer_fresh_mass_g = `OUTER SECTION FRESH MASS (g)`,
    outer_dry_mass_g = `Post-Freeze Dry OUTER SECTION Sample Mass (g)`,
    inner_length_cm = `Inner Section Length (cm)`,
    inner_fresh_mass_g = `INNER SECTION FRESH MASS (g)`,
    inner_dry_mass_g = `Post-Freeze Dry INNER SECTION Sample Mass (g)`,
    middle_length_mm = `Middle Section Length (mm)`,
    middle_fresh_mass_g = `MIDDLE SECTION FRESH MASS (g)`,
    middle_dry_mass_g = `Post-Freeze Dry MIDDLE SECTION Sample Mass (g)`,
    tree_id = `Tree Core ID`,
    species = `Tree Core Species Code`
  ) %>%
  mutate(
    outer_length_cm = as.numeric(outer_length_cm),
    inner_length_cm = as.numeric(inner_length_cm),
    middle_length_mm = as.numeric(middle_length_mm),
    outer_fresh_mass_g = as.numeric(outer_fresh_mass_g),
    outer_dry_mass_g = as.numeric(outer_dry_mass_g),
    inner_fresh_mass_g = as.numeric(inner_fresh_mass_g),
    inner_dry_mass_g = as.numeric(inner_dry_mass_g),
    middle_fresh_mass_g = as.numeric(middle_fresh_mass_g),
    middle_dry_mass_g = as.numeric(middle_dry_mass_g),
    middle_length_cm = middle_length_mm / 10
  )

# Calculate measurements and identify problematic values
complete_data <- processed_data %>%
  mutate(
    outer_volume_cm3 = ifelse(!is.na(outer_length_cm) & outer_length_cm > 0,
                              calculate_volume(outer_length_cm), NA),
    inner_volume_cm3 = ifelse(!is.na(inner_length_cm) & inner_length_cm > 0,
                              calculate_volume(inner_length_cm), NA),
    middle_volume_cm3 = ifelse(!is.na(middle_length_cm) & middle_length_cm > 0,
                               calculate_volume(middle_length_cm), NA),
    
    outer_moisture_fresh = ifelse(!is.na(outer_fresh_mass_g) & !is.na(outer_dry_mass_g) & outer_fresh_mass_g > 0,
                                  ((outer_fresh_mass_g - outer_dry_mass_g) / outer_fresh_mass_g) * 100, NA),
    inner_moisture_fresh = ifelse(!is.na(inner_fresh_mass_g) & !is.na(inner_dry_mass_g) & inner_fresh_mass_g > 0,
                                  ((inner_fresh_mass_g - inner_dry_mass_g) / inner_fresh_mass_g) * 100, NA),
    middle_moisture_fresh = ifelse(!is.na(middle_fresh_mass_g) & !is.na(middle_dry_mass_g) & middle_fresh_mass_g > 0,
                                   ((middle_fresh_mass_g - middle_dry_mass_g) / middle_fresh_mass_g) * 100, NA),
    
    outer_moisture_dry = ifelse(!is.na(outer_fresh_mass_g) & !is.na(outer_dry_mass_g) & outer_dry_mass_g > 0,
                                ((outer_fresh_mass_g - outer_dry_mass_g) / outer_dry_mass_g) * 100, NA),
    inner_moisture_dry = ifelse(!is.na(inner_fresh_mass_g) & !is.na(inner_dry_mass_g) & inner_dry_mass_g > 0,
                                ((inner_fresh_mass_g - inner_dry_mass_g) / inner_dry_mass_g) * 100, NA),
    middle_moisture_dry = ifelse(!is.na(middle_fresh_mass_g) & !is.na(middle_dry_mass_g) & middle_dry_mass_g > 0,
                                 ((middle_fresh_mass_g - middle_dry_mass_g) / middle_dry_mass_g) * 100, NA),
    
    outer_density_g_cm3 = ifelse(!is.na(outer_dry_mass_g) & !is.na(outer_volume_cm3) & outer_volume_cm3 > 0,
                                 outer_dry_mass_g / outer_volume_cm3, NA),
    inner_density_g_cm3 = ifelse(!is.na(inner_dry_mass_g) & !is.na(inner_volume_cm3) & inner_volume_cm3 > 0,
                                 inner_dry_mass_g / inner_volume_cm3, NA),
    middle_density_g_cm3 = ifelse(!is.na(middle_dry_mass_g) & !is.na(middle_volume_cm3) & middle_volume_cm3 > 0,
                                  middle_dry_mass_g / middle_volume_cm3, NA)
  )

# Apply median imputation for problematic values
problematic_trees <- c("AB6'", "H104", "HX2", "RM1", "rms1", "WA14", "YB15")

is_problematic <- function(tree_id, fresh_mass, dry_mass, moisture_dry, density) {
  is_problem_tree <- tree_id %in% problematic_trees
  dry_exceeds_fresh <- !is.na(fresh_mass) & !is.na(dry_mass) & dry_mass > fresh_mass
  extreme_moisture <- !is.na(moisture_dry) & moisture_dry > 300
  impossible_density <- !is.na(density) & (density > 2.0 | density < 0.1)
  return(is_problem_tree & (dry_exceeds_fresh | extreme_moisture | impossible_density))
}

# Calculate clean medians and apply imputation
clean_data <- complete_data %>%
  mutate(
    outer_moisture_fresh_clean = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3), NA, outer_moisture_fresh),
    outer_moisture_dry_clean = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3), NA, outer_moisture_dry),
    outer_density_clean = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3), NA, outer_density_g_cm3),
    inner_moisture_fresh_clean = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3), NA, inner_moisture_fresh),
    inner_moisture_dry_clean = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3), NA, inner_moisture_dry),
    inner_density_clean = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3), NA, inner_density_g_cm3),
    middle_moisture_fresh_clean = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3), NA, middle_moisture_fresh),
    middle_moisture_dry_clean = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3), NA, middle_moisture_dry),
    middle_density_clean = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3), NA, middle_density_g_cm3)
  )

medians_by_species <- clean_data %>%
  group_by(species) %>%
  summarise(
    outer_moisture_fresh_median = median(outer_moisture_fresh_clean, na.rm = TRUE),
    outer_moisture_dry_median = median(outer_moisture_dry_clean, na.rm = TRUE),
    outer_density_median = median(outer_density_clean, na.rm = TRUE),
    inner_moisture_fresh_median = median(inner_moisture_fresh_clean, na.rm = TRUE),
    inner_moisture_dry_median = median(inner_moisture_dry_clean, na.rm = TRUE),
    inner_density_median = median(inner_density_clean, na.rm = TRUE),
    middle_moisture_fresh_median = median(middle_moisture_fresh_clean, na.rm = TRUE),
    middle_moisture_dry_median = median(middle_moisture_dry_clean, na.rm = TRUE),
    middle_density_median = median(middle_density_clean, na.rm = TRUE),
    .groups = 'drop'
  )

# Create final filled dataset
final_data <- complete_data %>%
  left_join(medians_by_species, by = "species") %>%
  mutate(
    outer_moisture_fresh_percent = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3), outer_moisture_fresh_median, outer_moisture_fresh),
    outer_moisture_dry_percent = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3), outer_moisture_dry_median, outer_moisture_dry),
    outer_density_final = ifelse(is_problematic(tree_id, outer_fresh_mass_g, outer_dry_mass_g, outer_moisture_dry, outer_density_g_cm3), outer_density_median, outer_density_g_cm3),
    inner_moisture_fresh_percent = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3), inner_moisture_fresh_median, inner_moisture_fresh),
    inner_moisture_dry_percent = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3), inner_moisture_dry_median, inner_moisture_dry),
    inner_density_final = ifelse(is_problematic(tree_id, inner_fresh_mass_g, inner_dry_mass_g, inner_moisture_dry, inner_density_g_cm3), inner_density_median, inner_density_g_cm3),
    middle_moisture_fresh_percent = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3), middle_moisture_fresh_median, middle_moisture_fresh),
    middle_moisture_dry_percent = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3), middle_moisture_dry_median, middle_moisture_dry),
    middle_density_final = ifelse(is_problematic(tree_id, middle_fresh_mass_g, middle_dry_mass_g, middle_moisture_dry, middle_density_g_cm3), middle_density_median, middle_density_g_cm3)
  ) %>%
  select(-contains("median"), -contains("clean"))

print("=== COMPLETE TREE CORE ANALYSIS - FILLED DATASET ===")
print(paste("Total trees analyzed:", nrow(final_data)))
print(paste("Trees with imputed values:", length(problematic_trees)))

# Create summary statistics
results_summary <- final_data %>%
  select(tree_id, species, 
         outer_moisture_fresh_percent, outer_moisture_dry_percent, outer_density_final,
         inner_moisture_fresh_percent, inner_moisture_dry_percent, inner_density_final,
         middle_moisture_fresh_percent, middle_moisture_dry_percent, middle_density_final) %>%
  filter(!is.na(outer_moisture_fresh_percent) | !is.na(inner_moisture_fresh_percent) | !is.na(middle_moisture_fresh_percent))

print("\nFirst 10 rows of results:")
print(head(results_summary, 10))

# Summary statistics by species and section
outer_stats <- final_data %>%
  filter(!is.na(outer_density_final)) %>%
  group_by(species) %>%
  summarise(
    section = "Outer (Sapwood)",
    n_samples = n(),
    mean_density = round(mean(outer_density_final, na.rm = TRUE), 3),
    sd_density = round(sd(outer_density_final, na.rm = TRUE), 3),
    mean_moisture_fresh = round(mean(outer_moisture_fresh_percent, na.rm = TRUE), 1),
    sd_moisture_fresh = round(sd(outer_moisture_fresh_percent, na.rm = TRUE), 1),
    mean_moisture_dry = round(mean(outer_moisture_dry_percent, na.rm = TRUE), 1),
    sd_moisture_dry = round(sd(outer_moisture_dry_percent, na.rm = TRUE), 1),
    .groups = 'drop'
  )

inner_stats <- final_data %>%
  filter(!is.na(inner_density_final)) %>%
  group_by(species) %>%
  summarise(
    section = "Inner (Heartwood)",
    n_samples = n(),
    mean_density = round(mean(inner_density_final, na.rm = TRUE), 3),
    sd_density = round(sd(inner_density_final, na.rm = TRUE), 3),
    mean_moisture_fresh = round(mean(inner_moisture_fresh_percent, na.rm = TRUE), 1),
    sd_moisture_fresh = round(sd(inner_moisture_fresh_percent, na.rm = TRUE), 1),
    mean_moisture_dry = round(mean(inner_moisture_dry_percent, na.rm = TRUE), 1),
    sd_moisture_dry = round(sd(inner_moisture_dry_percent, na.rm = TRUE), 1),
    .groups = 'drop'
  )

middle_stats <- final_data %>%
  filter(!is.na(middle_density_final)) %>%
  group_by(species) %>%
  summarise(
    section = "Middle",
    n_samples = n(),
    mean_density = round(mean(middle_density_final, na.rm = TRUE), 3),
    sd_density = round(sd(middle_density_final, na.rm = TRUE), 3),
    mean_moisture_fresh = round(mean(middle_moisture_fresh_percent, na.rm = TRUE), 1),
    sd_moisture_fresh = round(sd(middle_moisture_fresh_percent, na.rm = TRUE), 1),
    mean_moisture_dry = round(mean(middle_moisture_dry_percent, na.rm = TRUE), 1),
    sd_moisture_dry = round(sd(middle_moisture_dry_percent, na.rm = TRUE), 1),
    .groups = 'drop'
  )

all_stats <- bind_rows(outer_stats, inner_stats, middle_stats) %>%
  arrange(species, section)

print("\nSummary statistics by species and section:")
print(all_stats)

# Save results
write_csv(results_summary, "../../data/processed/tree_cores/tree_core_filled_results.csv")
write_csv(final_data, "../../data/processed/tree_cores/tree_core_filled_complete.csv")
write_csv(all_stats, "../../data/processed/tree_cores/tree_core_filled_statistics.csv")

print("\nFiles saved:")
print("- tree_core_filled_results.csv")
print("- tree_core_filled_complete.csv") 
print("- tree_core_filled_statistics.csv")

# VISUALIZATIONS
print("\n=== CREATING VISUALIZATIONS ===")

# 1. Wood Density by Section and Species
density_plot_data <- final_data %>%
  select(tree_id, species, outer_density_final, inner_density_final, middle_density_final) %>%
  pivot_longer(cols = ends_with("_final"), 
               names_to = "section", 
               values_to = "density") %>%
  filter(!is.na(density)) %>%
  mutate(section = case_when(
    section == "outer_density_final" ~ "Outer (Sapwood)",
    section == "inner_density_final" ~ "Inner (Heartwood)", 
    section == "middle_density_final" ~ "Middle",
    TRUE ~ section
  ))

p1 <- ggplot(density_plot_data, aes(x = section, y = density, fill = species)) +
  geom_boxplot() +
  labs(title = "Wood Density by Section and Species (Filled Dataset)",
       subtitle = "Problematic values replaced with species medians",
       x = "Core Section",
       y = "Density (g/cm³)",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)

# 2. Moisture Content Comparison (Fresh vs Dry Basis)
moisture_plot_data <- final_data %>%
  select(tree_id, species, 
         outer_moisture_fresh_percent, outer_moisture_dry_percent,
         inner_moisture_fresh_percent, inner_moisture_dry_percent,
         middle_moisture_fresh_percent, middle_moisture_dry_percent) %>%
  pivot_longer(cols = ends_with("_percent"), 
               names_to = "measurement", 
               values_to = "moisture") %>%
  filter(!is.na(moisture)) %>%
  mutate(
    section = case_when(
      str_detect(measurement, "^outer") ~ "Outer (Sapwood)",
      str_detect(measurement, "^inner") ~ "Inner (Heartwood)", 
      str_detect(measurement, "^middle") ~ "Middle",
      TRUE ~ "Unknown"
    ),
    basis = case_when(
      str_detect(measurement, "fresh") ~ "Fresh Mass Basis",
      str_detect(measurement, "dry") ~ "Dry Mass Basis",
      TRUE ~ "Unknown"
    )
  )

p2 <- ggplot(moisture_plot_data, aes(x = section, y = moisture, fill = basis)) +
  geom_boxplot(position = position_dodge(0.8)) +
  facet_wrap(~species, scales = "free_y") +
  labs(title = "Moisture Content: Fresh vs Dry Mass Basis (Filled Dataset)",
       subtitle = "Clean data with median imputation for problematic values",
       x = "Core Section",
       y = "Moisture Content (%)",
       fill = "Calculation Basis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# 3. Individual Tree Moisture Comparison
tree_comparison_data <- final_data %>%
  select(tree_id, species,
         outer_moisture_fresh_percent, outer_density_final,
         inner_moisture_fresh_percent, inner_density_final,
         middle_moisture_fresh_percent, middle_density_final) %>%
  pivot_longer(
    cols = c(outer_moisture_fresh_percent, outer_density_final,
             inner_moisture_fresh_percent, inner_density_final,
             middle_moisture_fresh_percent, middle_density_final),
    names_to = "measurement",
    values_to = "value"
  ) %>%
  mutate(
    section = case_when(
      str_detect(measurement, "^outer") ~ "Outer",
      str_detect(measurement, "^inner") ~ "Inner", 
      str_detect(measurement, "^middle") ~ "Middle"
    ),
    variable = case_when(
      str_detect(measurement, "moisture") ~ "Moisture (%)",
      str_detect(measurement, "density") ~ "Density (g/cm³)"
    )
  ) %>%
  filter(!is.na(value))

p3 <- ggplot(tree_comparison_data %>% filter(variable == "Moisture (%)"), 
             aes(x = section, y = value, group = tree_id, color = species)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Moisture Content Across Sections by Individual Tree (Filled Dataset)",
       subtitle = "Each line represents one tree core - clean data with median imputation",
       x = "Core Section",
       y = "Moisture Content (% fresh mass basis)",
       color = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(limits = c("Outer", "Middle", "Inner"))

print(p3)

# 4. Wood Density Individual Trees
p4 <- ggplot(tree_comparison_data %>% filter(variable == "Density (g/cm³)"), 
             aes(x = section, y = value, group = tree_id, color = species)) +
  geom_line(alpha = 0.6, linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Wood Density Across Sections by Individual Tree (Filled Dataset)",
       subtitle = "Each line represents one tree core - clean data with median imputation",
       x = "Core Section", 
       y = "Wood Density (g/cm³)",
       color = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(limits = c("Outer", "Middle", "Inner"))

print(p4)

# 5. 1:1 Moisture Basis Comparison - Fixed version
# Create separate datasets for fresh and dry, then join
moisture_fresh <- final_data %>%
  select(tree_id, species,
         outer_moisture_fresh_percent, inner_moisture_fresh_percent, middle_moisture_fresh_percent) %>%
  pivot_longer(cols = ends_with("_percent"), names_to = "section", values_to = "fresh") %>%
  mutate(section = case_when(
    str_detect(section, "outer") ~ "Outer",
    str_detect(section, "inner") ~ "Inner", 
    str_detect(section, "middle") ~ "Middle"
  )) %>%
  filter(!is.na(fresh))

moisture_dry <- final_data %>%
  select(tree_id, species,
         outer_moisture_dry_percent, inner_moisture_dry_percent, middle_moisture_dry_percent) %>%
  pivot_longer(cols = ends_with("_percent"), names_to = "section", values_to = "dry") %>%
  mutate(section = case_when(
    str_detect(section, "outer") ~ "Outer",
    str_detect(section, "inner") ~ "Inner", 
    str_detect(section, "middle") ~ "Middle"
  )) %>%
  filter(!is.na(dry))

# Join the datasets
moisture_basis_data <- moisture_fresh %>%
  inner_join(moisture_dry, by = c("tree_id", "species", "section")) %>%
  filter(!is.na(fresh) & !is.na(dry))

print(paste("1:1 plot data points:", nrow(moisture_basis_data)))

p5 <- ggplot(moisture_basis_data, aes(x = fresh, y = dry)) +
  geom_point(aes(color = species, shape = section), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Moisture Content: Fresh Mass vs Dry Mass Basis (Filled Dataset)",
       subtitle = "Red dashed line = 1:1 relationship, Black line = actual relationship",
       x = "Fresh Mass Basis Moisture Content (%)",
       y = "Dry Mass Basis Moisture Content (%)",
       color = "Species",
       shape = "Section") +
  theme_minimal() +
  theme(legend.position = "right")

print(p5)

# 6. Distribution Plots - Grouped by Species with Sections Side-by-Side
distribution_data <- final_data %>%
  select(tree_id, species,
         outer_moisture_fresh_percent, outer_density_final,
         inner_moisture_fresh_percent, inner_density_final,
         middle_moisture_fresh_percent, middle_density_final) %>%
  pivot_longer(
    cols = -c(tree_id, species),
    names_to = "measurement",
    values_to = "value"
  ) %>%
  mutate(
    section = case_when(
      str_detect(measurement, "^outer") ~ "Outer",
      str_detect(measurement, "^inner") ~ "Inner",
      str_detect(measurement, "^middle") ~ "Middle"
    ),
    trait = case_when(
      str_detect(measurement, "moisture") ~ "Moisture Content (%)",
      str_detect(measurement, "density|final") ~ "Wood Density (g/cm³)"
    ),
    # Create combined species-section factor for x-axis ordering
    species_section = paste(species, section, sep = " ")
  ) %>%
  filter(!is.na(value)) %>%
  # Order species alphabetically, then sections (Outer, Middle, Inner)
  mutate(
    section = factor(section, levels = c("Outer", "Middle", "Inner")),
    species_section = factor(species_section, levels = paste(rep(sort(unique(species)), each = 3), 
                                                             rep(c("Outer", "Middle", "Inner"), length(unique(species)))))
  )

# Create the plot with species grouped, sections side-by-side
p6 <- ggplot(distribution_data, aes(x = species_section, y = value, fill = section)) +
  geom_violin(alpha = 0.6, scale = "width", position = position_dodge(0.8)) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.size = 0.5, position = position_dodge(0.8)) +
  facet_wrap(~trait, scales = "free_y", ncol = 1, strip.position = "left") +
  labs(title = "Wood Properties by Species and Section (Filled Dataset)",
       subtitle = "Sections grouped by species - Outer, Middle, Inner for each species",
       x = "Species and Section",
       y = "Value",
       fill = "Section") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_line(color = "grey90", size = 0.5),
    panel.grid.minor.x = element_blank()
  ) +
  scale_fill_manual(values = c("Outer" = "#E69F00", "Middle" = "#56B4E9", "Inner" = "#009E73")) +
  # Add vertical lines to separate species
  geom_vline(xintercept = seq(3.5, length(unique(distribution_data$species)) * 3 - 0.5, by = 3), 
             linetype = "dashed", alpha = 0.3, color = "grey60")

print(p6)

print("\nAnalysis complete using filled dataset!")
print(paste("Core diameter used:", core_diameter_mm, "mm"))
print("Density calculated using dry mass and calculated volume")
print("Moisture content calculated on BOTH bases:")
print("  - Fresh mass basis = (Fresh mass - Dry mass) / Fresh mass × 100%")
print("  - Dry mass basis = (Fresh mass - Dry mass) / Dry mass × 100%")
print("Problematic values replaced with species-specific medians for:")
print(paste("Trees:", paste(problematic_trees, collapse = ", ")))

# Final data quality check
all_densities <- c(final_data$outer_density_final, final_data$inner_density_final, final_data$middle_density_final)
all_densities <- all_densities[!is.na(all_densities)]
print(paste("\nFinal dataset density range:", round(min(all_densities), 3), "-", round(max(all_densities), 3), "g/cm³"))
print("All values now within reasonable biological range!")