# ==============================================================================
# Confirm ddPCR Import Integrity
# ==============================================================================
# Purpose: Validates processed ddPCR import by comparing individual
#   concentration values against original files.
#
# Pipeline stage: 01 Molecular Processing
# Run after: 02_check_completeness.R
#
# Inputs:
#   - processed_ddpcr_data.csv
#   - Original ddPCR files (for comparison)
#
# Outputs:
#   - Validation report (console output)
# ==============================================================================

library(dplyr)
library(tidyr)

# Assume files are already loaded: ddpcr_full, ddpcr_transposed, merged_final
# And mapping functions are available

cat("=== DETAILED VALUE-LEVEL COMPARISON ===\n")

# =============================================================================
# PROCESS ORIGINAL FILES TO MATCH PIPELINE FORMAT
# =============================================================================

# Process full ddPCR data exactly like the pipeline does
ddpcr_full_pipeline_format <- ddpcr_full %>%
  # Extract tree_id from Inner.Core.Sample.ID (same as pipeline)
  mutate(
    tree_id_raw = sub(".*_", "", sub("\\n.*", "", `Inner.Core.Sample.ID`)),
    tree_id = standardize_tree_id(tree_id_raw)
  ) %>%
  # Filter exactly like the pipeline
  filter(!is.na(tree_id), !is.na(Target), !is.na(core_type), tree_id != "untagged") %>%
  select(tree_id, Target, core_type, concentration = matches("Conc\\.copies")) %>%
  # Remove rows with missing concentration data
  filter(!is.na(concentration))

cat("Original ddPCR data (pipeline format):\n")
cat("Total measurements:", nrow(ddpcr_full_pipeline_format), "\n")
cat("Unique trees:", length(unique(ddpcr_full_pipeline_format$tree_id)), "\n")
cat("Unique targets:", paste(unique(ddpcr_full_pipeline_format$Target), collapse = ", "), "\n")
cat("Unique core types:", paste(unique(ddpcr_full_pipeline_format$core_type), collapse = ", "), "\n")

# =============================================================================
# CONVERT FINAL MERGED DATA BACK TO LONG FORMAT
# =============================================================================

# Convert final merged data back to long format for comparison
ddpcr_final_long <- merged_final %>%
  select(tree_id, starts_with("ddpcr_")) %>%
  # Convert to long format
  pivot_longer(
    cols = starts_with("ddpcr_"),
    names_to = "column_name",
    values_to = "concentration",
    values_drop_na = TRUE
  ) %>%
  # Parse the column names to extract target, core_type, analysis_type
  mutate(
    # Remove ddpcr_ prefix
    clean_name = str_remove(column_name, "^ddpcr_"),
    # Split by underscores - need to handle complex target names
    parts = str_split(clean_name, "_"),
    # Extract components - this is tricky because target names can have underscores
    analysis_type = map_chr(parts, ~ tail(.x, 1)),  # Last part is analysis type
    core_type = map_chr(parts, ~ .x[length(.x) - 1]),  # Second to last is core type
    # Target is everything except the last two parts
    Target = map_chr(parts, ~ paste(head(.x, length(.x) - 2), collapse = "_"))
  ) %>%
  select(tree_id, Target, core_type, analysis_type, concentration) %>%
  filter(!is.na(concentration))

cat("\nFinal merged data (converted to long):\n")
cat("Total measurements:", nrow(ddpcr_final_long), "\n")
cat("Unique trees:", length(unique(ddpcr_final_long$tree_id)), "\n")
cat("Unique targets:", paste(unique(ddpcr_final_long$Target), collapse = ", "), "\n")
cat("Unique core types:", paste(unique(ddpcr_final_long$core_type), collapse = ", "), "\n")
cat("Unique analysis types:", paste(unique(ddpcr_final_long$analysis_type), collapse = ", "), "\n")

# =============================================================================
# CREATE COMPARISON KEY AND CHECK FOR MISSING VALUES
# =============================================================================

# Create comparison keys for original data
ddpcr_original_keys <- ddpcr_full_pipeline_format %>%
  # Note: Original data doesn't have analysis_type distinction, so we need to handle this
  mutate(
    comparison_key = paste(tree_id, Target, core_type, sep = "_")
  ) %>%
  select(comparison_key, tree_id, Target, core_type, concentration)

# Create comparison keys for final data - we need to account for loose/strict splits
ddpcr_final_keys <- ddpcr_final_long %>%
  mutate(
    comparison_key = paste(tree_id, Target, core_type, sep = "_")
  ) %>%
  # Group by the key and check if we have both loose and strict values
  group_by(comparison_key, tree_id, Target, core_type) %>%
  summarise(
    n_analysis_types = n(),
    has_loose = any(analysis_type == "loose"),
    has_strict = any(analysis_type == "strict"),
    loose_conc = ifelse(has_loose, concentration[analysis_type == "loose"][1], NA),
    strict_conc = ifelse(has_strict, concentration[analysis_type == "strict"][1], NA),
    .groups = 'drop'
  )

cat("\n=== COMPARISON RESULTS ===\n")

# Check which original measurements are missing from final
missing_measurements <- ddpcr_original_keys %>%
  filter(!comparison_key %in% ddpcr_final_keys$comparison_key)

cat("Original measurements missing from final:", nrow(missing_measurements), "\n")
if (nrow(missing_measurements) > 0) {
  cat("Missing measurements:\n")
  print(missing_measurements %>% select(tree_id, Target, core_type, concentration))
}

# Check which final measurements weren't in original
extra_measurements <- ddpcr_final_keys %>%
  filter(!comparison_key %in% ddpcr_original_keys$comparison_key)

cat("\nFinal measurements not in original:", nrow(extra_measurements), "\n")
if (nrow(extra_measurements) > 0) {
  cat("Extra measurements (first 10):\n")
  print(head(extra_measurements %>% select(tree_id, Target, core_type, n_analysis_types), 10))
}

# =============================================================================
# DETAILED COMPARISON FOR SPECIFIC TREE-TARGET COMBINATIONS
# =============================================================================

cat("\n=== DETAILED VALUE COMPARISON FOR SAMPLE MEASUREMENTS ===\n")

# Compare actual concentration values for matching measurements
value_comparison <- ddpcr_original_keys %>%
  inner_join(
    ddpcr_final_keys %>% select(comparison_key, loose_conc, strict_conc),
    by = "comparison_key"
  ) %>%
  mutate(
    # Check if original value matches either loose or strict
    matches_loose = abs(concentration - loose_conc) < 0.001,
    matches_strict = abs(concentration - strict_conc) < 0.001,
    matches_either = matches_loose | matches_strict
  )

# Count matching vs non-matching values
matching_count <- sum(value_comparison$matches_either, na.rm = TRUE)
total_comparison <- nrow(value_comparison)

cat("Measurements with matching values:", matching_count, "/", total_comparison, "\n")
cat("Percentage match:", round((matching_count / total_comparison) * 100, 1), "%\n")

# Show examples of non-matching values
non_matching <- value_comparison %>%
  filter(!matches_either) %>%
  select(tree_id, Target, core_type, original = concentration, loose = loose_conc, strict = strict_conc)

if (nrow(non_matching) > 0) {
  cat("\nNon-matching values (first 10):\n")
  print(head(non_matching, 10))
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== SUMMARY ===\n")
cat("Original data points:", nrow(ddpcr_original_keys), "\n")
cat("Final data points:", nrow(ddpcr_final_long), "\n")
cat("Missing from final:", nrow(missing_measurements), "\n")
cat("Added in final:", nrow(extra_measurements), "\n")
cat("Value matches:", matching_count, "/", total_comparison, "\n")

if (nrow(missing_measurements) == 0) {
  cat("\nGOOD: No original ddPCR measurements were lost.\n")
} else {
  cat("\nWARNING: Some original ddPCR measurements are missing from final data.\n")
}

if (nrow(non_matching) > 0) {
  cat("WARNING: Some concentration values don't match between original and final.\n")
} else {
  cat("GOOD: All concentration values match between original and final data.\n")
}