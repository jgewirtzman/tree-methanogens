# ==============================================================================
# Harmonize All Data
# ==============================================================================
# Purpose: Master data integration merging tree IDs, DBH, flux, soil, wood
#   cores, internal gas, and ddPCR gene data into a single dataset.
#
# Pipeline stage: 02 Integration
# Run after: All Stage 1 processing scripts
#
# Inputs:
#   - tree_id_comprehensive_mapping.csv (from 03_tree_data/)
#   - tree_dbh_consensus_comprehensive.csv (from 03_tree_data/)
#   - sample_data_only.csv (from 03_tree_data/)
#   - methanogen_tree_flux_complete_dataset.csv (from 01_flux_processing/static/)
#   - soil_tree_level_means.csv (from 00_harmonization/)
#   - tree_core_filled_complete.csv (from 03_tree_data/)
#   - processed_ddpcr_data.csv (from 02_ddpcr/)
#
# Outputs:
#   - merged_tree_dataset_final.csv (KEY OUTPUT for downstream analysis)
# ==============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Load the mapping file and create lookup
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
        dplyr::select(Tree_ID_normalized, original_name = all_of(col)) %>%
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

processed_datasets <- list()

# =============================================================================
# PROCESS FILE 1: DBH DATA
# =============================================================================

cat("=== PROCESSING DBH DATA ===\n")

dbh <- read_csv("../../data/processed/tree_data/tree_dbh_consensus_comprehensive.csv")

dbh_processed <- dbh %>%
  dplyr::select(
    tree_id_raw = normalized_id,
    plot_name = Plot,
    dbh = DBH,
    source
  ) %>%
  # Standardize tree IDs
  mutate(
    tree_id = standardize_tree_id(tree_id_raw),
    source = "dbh_data"
  ) %>%
  # Remove untagged trees and rows with missing data
  filter(!is.na(tree_id), !is.na(dbh), tree_id != "untagged") %>%
  # Average duplicates
  group_by(tree_id) %>%
  summarise(
    plot_name = first(plot_name),
    dbh = mean(dbh, na.rm = TRUE),
    source = first(source),
    .groups = 'drop'
  )

processed_datasets[["dbh"]] <- dbh_processed
cat("DBH data: ", nrow(dbh_processed), " records for ", length(unique(dbh_processed$tree_id)), " unique trees\n")

# =============================================================================
# PROCESS FILE 2: GAS CONCENTRATION DATA  
# =============================================================================

cat("=== PROCESSING GAS CONCENTRATION DATA ===\n")

gas <- read_csv("../../data/processed/internal_gas/sample_data_only.csv")

gas_processed <- gas %>%
  dplyr::select(
    tree_id_raw = Tree.ID,
    species_id = Species.ID,
    CO2_concentration,
    CH4_concentration, 
    N2O_concentration,
    O2_concentration
  ) %>%
  # Standardize tree IDs
  mutate(
    tree_id = standardize_tree_id(tree_id_raw),
    source = "gas_concentrations"
  ) %>%
  # Remove untagged trees and rows with missing data
  filter(!is.na(tree_id), tree_id != "untagged") %>%
  dplyr::select(-tree_id_raw)

processed_datasets[["gas_concentrations"]] <- gas_processed
cat("Gas data: ", nrow(gas_processed), " records for ", length(unique(gas_processed$tree_id)), " unique trees\n")

# =============================================================================
# PROCESS FILE 3: FLUX DATA (with PB4? fix)
# =============================================================================

cat("=== PROCESSING FLUX DATA ===\n")

flux <- read_csv("../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv")

flux_processed <- flux %>%
  dplyr::select(
    tree_id_raw = tree_id,
    plot_name = plot,
    species_name = species,
    measurement_height,
    CO2_best.flux,
    CH4_best.flux,
    Temp_Air = Tcham
  ) %>%
  # Fix the PB4? to pb5 issue
  mutate(
    tree_id_raw = case_when(
      tree_id_raw == "PB4?" ~ "pb5",
      TRUE ~ tree_id_raw
    )
  ) %>%
  # Standardize tree IDs
  mutate(
    tree_id = standardize_tree_id(tree_id_raw)
  ) %>%
  # Remove untagged trees and rows with missing data
  filter(!is.na(tree_id), !is.na(measurement_height), tree_id != "untagged") %>%
  dplyr::select(-tree_id_raw)

# Convert to wide format with duplicate handling
flux_wide <- flux_processed %>%
  # Clean up measurement heights for column names
  mutate(
    height_clean = case_when(
      measurement_height == "Root Crown" ~ "RootCrown",
      measurement_height == "75 cm" ~ "75",
      TRUE ~ as.character(measurement_height)
    ),
    height_suffix = paste0("_", height_clean, "cm")
  ) %>%
  # Group by tree_id and height, average duplicates
  group_by(tree_id, height_suffix) %>%
  summarise(
    plot_name = first(plot_name),
    species_name = first(species_name),
    CO2_best.flux = mean(CO2_best.flux, na.rm = TRUE),
    CH4_best.flux = mean(CH4_best.flux, na.rm = TRUE), 
    Temp_Air = mean(Temp_Air, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Pivot wider
  pivot_wider(
    id_cols = c(tree_id, plot_name, species_name),
    names_from = height_suffix,
    values_from = c(CO2_best.flux, CH4_best.flux, Temp_Air),
    names_sep = ""
  ) %>%
  mutate(source = "flux_data")

processed_datasets[["flux"]] <- flux_wide
cat("Flux data: ", nrow(flux_wide), " records for ", length(unique(flux_wide$tree_id)), " unique trees\n")

# =============================================================================
# PROCESS FILE 4: SOIL DATA
# =============================================================================

cat("=== PROCESSING SOIL DATA ===\n")

soil <- read_csv("../../data/raw/field_data/static_chamber_field/soil_tree_level_means.csv")

soil_processed <- soil %>%
  dplyr::select(
    tree_id_raw = Tree_ID,
    VWC_mean,
    ORP_mean, 
    SoilTemp_mean,
    OrganicDepth_mean,
    MineralDepth_mean
  ) %>%
  # Standardize tree IDs
  mutate(
    tree_id = standardize_tree_id(tree_id_raw),
    source = "soil_data"
  ) %>%
  # Remove untagged trees and rows with missing data
  filter(!is.na(tree_id), tree_id != "untagged") %>%
  dplyr::select(-tree_id_raw)

processed_datasets[["soil"]] <- soil_processed
cat("Soil data: ", nrow(soil_processed), " records for ", length(unique(soil_processed$tree_id)), " unique trees\n")

# =============================================================================
# PROCESS FILE 5: WOOD/TREE CORE DATA
# =============================================================================

cat("=== PROCESSING WOOD/TREE CORE DATA ===\n")

wood <- read_csv("../../data/processed/tree_cores/tree_core_filled_complete.csv")

wood_processed <- wood %>%
  dplyr::select(
    tree_id_raw = tree_id,
    species_name = species,
    outer_moisture_fresh_percent,
    inner_moisture_fresh_percent,
    middle_moisture_fresh_percent,
    outer_moisture_dry_percent,
    inner_moisture_dry_percent,
    middle_moisture_dry_percent,
    outer_density_final,
    inner_density_final,
    middle_density_final
  ) %>%
  # Standardize tree IDs
  mutate(
    tree_id = standardize_tree_id(tree_id_raw),
    source = "wood_core_data"
  ) %>%
  # Remove untagged trees and rows with missing data
  filter(!is.na(tree_id), tree_id != "untagged") %>%
  # Average duplicates
  group_by(tree_id) %>%
  summarise(
    species_name = first(species_name),
    across(ends_with("_percent"), ~ mean(.x, na.rm = TRUE)),
    across(ends_with("_final"), ~ mean(.x, na.rm = TRUE)),
    source = first(source),
    .groups = 'drop'
  )

processed_datasets[["wood"]] <- wood_processed
cat("Wood data: ", nrow(wood_processed), " records for ", length(unique(wood_processed$tree_id)), " unique trees\n")


# =============================================================================
# PROCESS FILE 6: ddPCR DATA (with unit conversion and verification)
# =============================================================================

cat("=== PROCESSING ddPCR DATA ===\n")

ddpcr <- read_csv("../../data/processed/molecular/processed_ddpcr_data.csv")

# First, let's examine the data structure to verify sample mass correspondence
cat("Examining sample mass structure...\n")
sample_mass_check <- ddpcr %>%
  dplyr::select(`Inner.Core.Sample.ID`, Target, core_type, analysis_type, `Sample.Mass.Added.to.Tube..mg.`) %>%
  group_by(`Inner.Core.Sample.ID`) %>%
  summarise(
    unique_masses = n_distinct(`Sample.Mass.Added.to.Tube..mg.`, na.rm = TRUE),
    mass_values = paste(unique(`Sample.Mass.Added.to.Tube..mg.`), collapse = ", "),
    .groups = 'drop'
  )

cat("Sample mass verification:\n")
cat("Samples with multiple different masses: ", sum(sample_mass_check$unique_masses > 1, na.rm = TRUE), "\n")
if(sum(sample_mass_check$unique_masses > 1, na.rm = TRUE) > 0) {
  cat("WARNING: Some samples have multiple different sample masses!\n")
  print(sample_mass_check %>% filter(unique_masses > 1))
}

ddpcr_processed <- ddpcr %>%
  # Extract tree_id from Inner.Core.Sample.ID
  mutate(
    tree_id_raw = sub(".*_", "", sub("\\n.*", "", `Inner.Core.Sample.ID`))
  ) %>%
  dplyr::select(
    tree_id_raw,
    species_name = species,
    Target,
    core_type,
    analysis_type,
    concentration_original = `Conc.copies.µL.`,
    sample_mass = `Sample.Mass.Added.to.Tube..mg.`
  ) %>%
  # Standardize tree IDs
  mutate(
    tree_id = standardize_tree_id(tree_id_raw)
  ) %>%
  # Remove untagged trees and rows with missing data
  filter(!is.na(tree_id), !is.na(Target), !is.na(core_type), !is.na(analysis_type), 
         !is.na(concentration_original), !is.na(sample_mass), tree_id != "untagged") %>%
  # Apply unit conversion: copies/µL to copies/g dry wood
  mutate(
    concentration_per_g = concentration_original * 75 / sample_mass * 1000
  ) %>%
  dplyr::select(-tree_id_raw, -concentration_original, -sample_mass) %>%
  rename(concentration = concentration_per_g)

# Convert to wide format with duplicate averaging
ddpcr_wide <- ddpcr_processed %>%
  # Create combined column name - keeping original structure
  mutate(
    column_name = paste(Target, core_type, analysis_type, sep = "_")
  ) %>%
  # Average duplicates within same tree+target+core+analysis combination
  group_by(tree_id, species_name, column_name) %>%
  summarise(concentration = mean(concentration, na.rm = TRUE), .groups = 'drop') %>%
  # Pivot wider - keeping original column naming pattern with ddpcr prefix
  pivot_wider(
    id_cols = c(tree_id, species_name),
    names_from = column_name,
    values_from = concentration,
    names_prefix = "ddpcr_"
  ) %>%
  mutate(source = "ddpcr_data")

processed_datasets[["ddpcr"]] <- ddpcr_wide
cat("ddPCR data: ", nrow(ddpcr_wide), " records for ", length(unique(ddpcr_wide$tree_id)), " unique trees\n")
cat("Units converted from copies/µL to copies/g dry wood\n")
cat("Column names preserved in original format\n")

# =============================================================================
# FINAL MERGE: Combine all processed datasets
# =============================================================================

cat("\n=== MERGING ALL DATASETS ===\n")

# Start with DBH data as the base
merged_data <- processed_datasets[["dbh"]]
cat("Starting with DBH data: ", nrow(merged_data), " trees\n")

# Merge each dataset sequentially
for(dataset_name in names(processed_datasets)[-1]) {
  dataset <- processed_datasets[[dataset_name]]
  before_count <- nrow(merged_data)
  
  merged_data <- merged_data %>%
    full_join(dataset, by = "tree_id")
  
  after_count <- nrow(merged_data)
  cat("Added ", dataset_name, " - Trees before: ", before_count, ", Trees after: ", after_count, "\n")
}

# Clean up the merged dataset
merged_data_final <- merged_data %>%
  # Remove duplicate source columns
  dplyr::select(-matches("source\\.[xy]$")) %>%
  # Handle duplicate plot and species columns
  {
    if ("plot_name" %in% names(.)) {
      mutate(., plot_final = plot_name)
    } else if ("plot_name.x" %in% names(.)) {
      mutate(., plot_final = coalesce(plot_name.x, plot_name.y))
    } else {
      mutate(., plot_final = NA_character_)
    }
  } %>%
  {
    if ("species_name" %in% names(.)) {
      mutate(., species_final = species_name)
    } else if ("species_name.x" %in% names(.)) {
      mutate(., species_final = coalesce(species_name.x, species_name.y))
    } else {
      mutate(., species_final = NA_character_)
    }
  } %>%
  # Remove duplicate columns
  dplyr::select(-matches("\\.[xy]$"), -matches("plot_name"), -matches("species_name")) %>%
  # Rename final columns
  rename(
    plot = plot_final,
    species = species_final
  ) %>%
  # Arrange by tree_id
  arrange(tree_id)

cat("\n=== MERGE COMPLETE ===\n")
cat("Final dataset dimensions: ", nrow(merged_data_final), " rows × ", ncol(merged_data_final), " columns\n")
cat("Total unique trees: ", length(unique(merged_data_final$tree_id)), "\n")

# Summary of data availability
cat("\nData availability summary:\n")
cat("Trees with DBH data: ", sum(!is.na(merged_data_final$dbh)), "\n")
cat("Trees with gas concentration data: ", sum(!is.na(merged_data_final$CO2_concentration)), "\n")
cat("Trees with flux data: ", sum(!is.na(merged_data_final$CO2_best.flux_50cm)), "\n") 
cat("Trees with soil data: ", sum(!is.na(merged_data_final$VWC_mean)), "\n")
cat("Trees with wood core data: ", sum(!is.na(merged_data_final$outer_density_final)), "\n")
cat("Trees with ddPCR data: ", sum(!is.na(merged_data_final$ddpcr_mcra_probe_Mineral_loose)), "\n")

# Save the final merged dataset
write_csv(merged_data_final, "../../data/processed/integrated/merged_tree_dataset_final.csv")
cat("\nFinal merged dataset saved as: 'data/processed/integrated/merged_tree_dataset_final.csv'\n")

# Show sample of final data
cat("\nSample of final dataset:\n")
print(merged_data_final %>% 
        dplyr::select(tree_id, dbh, CO2_concentration, CO2_best.flux_50cm, VWC_mean, outer_density_final) %>% 
        head(10))

cat("\nPIPELINE COMPLETE!\n")
cat("Changes made:\n")
cat("- Removed all 'untagged' trees from all datasets\n")
cat("- Fixed 'PB4?' to 'pb5' in flux data\n") 
cat("- Averaged remaining duplicates as planned\n")



















# =============================================================================
# SPECIES COMPLETION: Fill missing species based on tree_id patterns
# =============================================================================

# Species abbreviation mapping based on first two letters of tree_id
species_mapping <- c(
  "ab" = "FAGR",     # American beech
  "bb" = "BELE",     # Black birch  
  "bc" = "PRSE",     # Black cherry
  "bo" = "QUVE",     # Black oak
  "h" = "TSCA",      # Hemlock
  "ml" = "KALA",     # Mountain laurel
  "pb" = "BEPA",     # Paper birch
  "rm" = "ACRU",     # Red maple
  "ro" = "QURU",     # Red oak
  "sa" = "SAAL",     # Sassafras
  "sf" = "SAAL",     # Sassafras
  "sh" = "CAOV",     # Shagbark hickory
  "sm" = "ACSA",     # Sugar maple
  "wa" = "FRAM",     # White ash
  "wo" = "QUAL",     # White oak
  "wp" = "PIST",     # White pine
  "yb" = "BEAL"      # Yellow birch
)

# Function to determine species abbreviation from tree_id
determine_species_from_id <- function(tree_id) {
  # Handle special cases first
  if (tree_id == "bowser") return("QUVE")
  if (tree_id == "cul8r") return("BEAL")
  if (tree_id == "sassafrass") return("SAAL")
  
  # Handle numeric-only IDs
  if (grepl("^[0-9]+$", tree_id)) return(NA_character_)
  
  # Handle hemlock trees that start with 'h' followed by numbers or letters
  if (grepl("^h[0-9x-z]", tolower(tree_id))) return("TSCA")
  
  # Extract first two letters, convert to lowercase
  prefix <- tolower(substr(tree_id, 1, 2))
  
  # Look up species by prefix
  species <- species_mapping[prefix]
  return(ifelse(is.na(species), NA_character_, species))
}

# Fill in missing species values by overwriting species_id
merged_data_final <- merged_data_final %>%
  mutate(
    species_from_id = sapply(tree_id, determine_species_from_id),
    species_id = coalesce(species_id, species, species_from_id)
  ) %>%
  dplyr::select(-species_from_id, -species)

cat("\n=== SPECIES COMPLETION ===\n")
cat("Trees with species information: ", sum(!is.na(merged_data_final$species_id)), "\n")

# Re-save the final merged dataset with species completion
write_csv(merged_data_final, "../../data/processed/integrated/merged_tree_dataset_final.csv")
cat("Updated dataset saved with species completion\n")