# ==============================================================================
# Investigate Duplicate Samples
# ==============================================================================
# Purpose: Identifies and analyzes duplicate samples across data sources.
# ==============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Load the mapping and create lookup (same as before)
mapping <- read_csv("../../data/processed/tree_data/tree_id_comprehensive_mapping.csv")
mapping <- mapping %>% mutate(across(everything(), ~ifelse(.x == "NA", NA, .x)))

create_comprehensive_lookup <- function(mapping) {
  name_columns <- names(mapping)[grepl("^name_in_|^variant_", names(mapping))]
  name_columns <- c(name_columns, "primary_id")
  
  lookup_list <- list()
  for (col in name_columns) {
    if (col %in% names(mapping)) {
      temp_lookup <- mapping %>%
        filter(!is.na(.data[[col]])) %>%
        select(Tree_ID_normalized, original_name = all_of(col)) %>%
        mutate(original_name = as.character(original_name))
      lookup_list[[col]] <- temp_lookup
    }
  }
  
  comprehensive_lookup <- bind_rows(lookup_list) %>%
    distinct() %>%
    mutate(
      original_name_lower = tolower(trimws(original_name)),
      Tree_ID_normalized = tolower(trimws(Tree_ID_normalized))
    )
  return(comprehensive_lookup)
}

tree_lookup <- create_comprehensive_lookup(mapping)

standardize_tree_id <- function(tree_ids) {
  clean_ids <- tolower(trimws(as.character(tree_ids)))
  standardized <- tree_lookup$Tree_ID_normalized[match(clean_ids, tree_lookup$original_name_lower)]
  result <- ifelse(is.na(standardized), clean_ids, standardized)
  return(result)
}

# =============================================================================
# 1. DBH DATA DUPLICATES
# =============================================================================

cat("=== CHECKING DBH DUPLICATES ===\n")

dbh <- read_csv("../../data/processed/tree_data/tree_dbh_consensus_comprehensive.csv")

# Check for duplicates BEFORE standardization
dbh_raw_dups <- dbh %>%
  count(normalized_id) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("Raw DBH duplicates (before standardization):\n")
print(dbh_raw_dups)

if (nrow(dbh_raw_dups) > 0) {
  cat("\nExample DBH duplicate records:\n")
  example_id <- dbh_raw_dups$normalized_id[1]
  print(dbh %>% filter(normalized_id == example_id))
}

# Check for duplicates AFTER standardization
dbh_with_std <- dbh %>%
  mutate(tree_id = standardize_tree_id(normalized_id)) %>%
  filter(!is.na(tree_id), !is.na(DBH))

dbh_std_dups <- dbh_with_std %>%
  count(tree_id) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("\nDBH duplicates after standardization:\n")
print(dbh_std_dups)

if (nrow(dbh_std_dups) > 0) {
  cat("\nExample standardized DBH duplicate records:\n")
  example_std_id <- dbh_std_dups$tree_id[1]
  print(dbh_with_std %>% filter(tree_id == example_std_id) %>% select(tree_id, normalized_id, Plot, DBH, source))
}

# =============================================================================
# 2. GAS DATA DUPLICATES  
# =============================================================================

cat("\n=== CHECKING GAS DUPLICATES ===\n")

gas <- read_csv("../../data/processed/internal_gas/sample_data_only.csv")

# Check raw duplicates
gas_raw_dups <- gas %>%
  count(Tree.ID) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("Raw gas duplicates:\n")
print(gas_raw_dups)

# Check after standardization
gas_with_std <- gas %>%
  mutate(tree_id = standardize_tree_id(Tree.ID)) %>%
  filter(!is.na(tree_id))

gas_std_dups <- gas_with_std %>%
  count(tree_id) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("\nGas duplicates after standardization:\n")
print(gas_std_dups)

if (nrow(gas_std_dups) > 0) {
  cat("\nExample gas duplicate records:\n")
  example_gas_id <- gas_std_dups$tree_id[1]
  print(gas_with_std %>% filter(tree_id == example_gas_id) %>% select(tree_id, Tree.ID, Species.ID, CO2_concentration, CH4_concentration))
}

# =============================================================================
# 3. FLUX DATA DUPLICATES
# =============================================================================

cat("\n=== CHECKING FLUX DUPLICATES ===\n")

flux <- read_csv("../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv")

# Check raw duplicates
flux_raw_dups <- flux %>%
  count(tree_id) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("Raw flux duplicates (by tree_id):\n")
print(head(flux_raw_dups, 10))

# Check after standardization
flux_with_std <- flux %>%
  mutate(tree_id_std = standardize_tree_id(tree_id)) %>%
  filter(!is.na(tree_id_std), !is.na(measurement_height))

flux_std_dups <- flux_with_std %>%
  count(tree_id_std) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("\nFlux duplicates after standardization:\n")
print(head(flux_std_dups, 10))

# Check for duplicates within tree + height combinations
flux_height_dups <- flux_with_std %>%
  count(tree_id_std, measurement_height) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("\nFlux duplicates by tree + height:\n")
print(head(flux_height_dups, 10))

if (nrow(flux_height_dups) > 0) {
  cat("\nExample flux duplicate records (same tree + height):\n")
  example_flux_tree <- flux_height_dups$tree_id_std[1]
  example_flux_height <- flux_height_dups$measurement_height[1]
  print(flux_with_std %>% 
          filter(tree_id_std == example_flux_tree, measurement_height == example_flux_height) %>%
          select(tree_id_std, tree_id, plot, species, measurement_height, CO2_best.flux, CH4_best.flux, start.time))
}

# =============================================================================
# 4. SOIL DATA DUPLICATES
# =============================================================================

cat("\n=== CHECKING SOIL DUPLICATES ===\n")

soil <- read_csv("../../data/raw/field_data/static_chamber_field/soil_tree_level_means.csv")

# Check raw duplicates
soil_raw_dups <- soil %>%
  count(Tree_ID) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("Raw soil duplicates:\n")
print(soil_raw_dups)

# Check after standardization
soil_with_std <- soil %>%
  mutate(tree_id = standardize_tree_id(Tree_ID)) %>%
  filter(!is.na(tree_id))

soil_std_dups <- soil_with_std %>%
  count(tree_id) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("\nSoil duplicates after standardization:\n")
print(soil_std_dups)

# =============================================================================
# 5. WOOD DATA DUPLICATES
# =============================================================================

cat("\n=== CHECKING WOOD DUPLICATES ===\n")

wood <- read_csv("../../data/processed/tree_cores/tree_core_filled_complete.csv")

# Check raw duplicates
wood_raw_dups <- wood %>%
  count(tree_id) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("Raw wood duplicates:\n")
print(head(wood_raw_dups, 10))

# Check after standardization
wood_with_std <- wood %>%
  mutate(tree_id_std = standardize_tree_id(tree_id)) %>%
  filter(!is.na(tree_id_std))

wood_std_dups <- wood_with_std %>%
  count(tree_id_std) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("\nWood duplicates after standardization:\n")
print(head(wood_std_dups, 10))

if (nrow(wood_std_dups) > 0) {
  cat("\nExample wood duplicate records:\n")
  example_wood_id <- wood_std_dups$tree_id_std[1]
  print(wood_with_std %>% 
          filter(tree_id_std == example_wood_id) %>%
          select(tree_id_std, tree_id, species, Date, outer_density_final, inner_density_final))
}

# =============================================================================
# 6. ddPCR DATA DUPLICATES
# =============================================================================

cat("\n=== CHECKING ddPCR DUPLICATES ===\n")

ddpcr <- read_csv("../../data/processed/molecular/processed_ddpcr_data.csv")

# Extract tree ID and check duplicates
ddpcr_with_extracted <- ddpcr %>%
  mutate(tree_id_extracted = sub(".*_", "", sub("\\n.*", "", `Inner.Core.Sample.ID`)))

ddpcr_raw_dups <- ddpcr_with_extracted %>%
  count(tree_id_extracted) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("Raw ddPCR duplicates (by extracted tree ID):\n")
print(head(ddpcr_raw_dups, 10))

# Check after standardization  
ddpcr_with_std <- ddpcr_with_extracted %>%
  mutate(tree_id_std = standardize_tree_id(tree_id_extracted)) %>%
  filter(!is.na(tree_id_std), !is.na(Target), !is.na(core_type), !is.na(analysis_type))

ddpcr_std_dups <- ddpcr_with_std %>%
  count(tree_id_std) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("\nddPCR duplicates after standardization:\n")
print(head(ddpcr_std_dups, 10))

# Check for duplicates within tree + target + core_type + analysis_type combinations
ddpcr_combo_dups <- ddpcr_with_std %>%
  count(tree_id_std, Target, core_type, analysis_type) %>%
  filter(n > 1) %>%
  arrange(desc(n))

cat("\nddPCR duplicates by tree + target + core_type + analysis_type:\n")
print(head(ddpcr_combo_dups, 10))

if (nrow(ddpcr_combo_dups) > 0) {
  cat("\nExample ddPCR duplicate records (same tree + target + core + analysis):\n")
  example_ddpcr <- ddpcr_combo_dups[1,]
  print(ddpcr_with_std %>% 
          filter(tree_id_std == example_ddpcr$tree_id_std, 
                 Target == example_ddpcr$Target,
                 core_type == example_ddpcr$core_type,
                 analysis_type == example_ddpcr$analysis_type) %>%
          select(tree_id_std, tree_id_extracted, Target, core_type, analysis_type, `Conc.copies.µL.`, Date.Weighed))
}

cat("\n=== DUPLICATE INVESTIGATION COMPLETE ===\n")
cat("Review the output above to understand why duplicates exist before deciding how to handle them.\n")
