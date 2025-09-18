# Simple join: Flux data + Surface Area + Volume
# Links flux observations to chamber surface area and volume only

library(readxl)
library(writexl)
library(dplyr)

# Read flux data
flux_data <- read_excel("/Users/jongewirtzman/My Drive/Matthes_Lab/stem-CH4-flux/yale-forest-ch4/ipad_data/Cleaned data/treeflux_total.xlsx", sheet = "Sheet1")

# Read geometry data
geometry_data <- read_excel("combined_corrected_results.xlsx", sheet = "all_data")
geometry_data <- results

# Clean and examine flux identifiers
flux_clean <- flux_data %>%
  mutate(
    Sample_ID = as.character(`Plot Tag`),
    Chamber_ID = as.character(Chamber)
  ) %>%
  filter(!is.na(Sample_ID) & !is.na(Chamber_ID))

# Clean and examine geometry identifiers  
geometry_clean <- geometry_data %>%
  mutate(
    Sample_ID = as.character(Sample),
    Chamber_ID = as.character(`Chamber used`)
  ) %>%
  filter(!is.na(Sample_ID) & !is.na(Chamber_ID))

# Diagnostic analysis
cat("IDENTIFIER COMPARISON:\n")
cat("======================\n\n")

cat("FLUX DATA:\n")
cat("Unique Sample IDs:", n_distinct(flux_clean$Sample_ID), "\n")
cat("Unique Chamber IDs:", n_distinct(flux_clean$Chamber_ID), "\n")
cat("Sample ID examples:", paste(head(sort(unique(flux_clean$Sample_ID)), 10), collapse = ", "), "\n")
cat("Chamber ID examples:", paste(head(sort(unique(flux_clean$Chamber_ID)), 10), collapse = ", "), "\n\n")

cat("GEOMETRY DATA:\n")
cat("Unique Sample IDs:", n_distinct(geometry_clean$Sample_ID), "\n") 
cat("Unique Chamber IDs:", n_distinct(geometry_clean$Chamber_ID), "\n")
cat("Sample ID examples:", paste(head(sort(unique(geometry_clean$Sample_ID)), 10), collapse = ", "), "\n")
cat("Chamber ID examples:", paste(head(sort(unique(geometry_clean$Chamber_ID)), 10), collapse = ", "), "\n\n")

# Find overlaps
flux_samples <- unique(flux_clean$Sample_ID)
geom_samples <- unique(geometry_clean$Sample_ID)
flux_chambers <- unique(flux_clean$Chamber_ID)
geom_chambers <- unique(geometry_clean$Chamber_ID)

common_samples <- intersect(flux_samples, geom_samples)
common_chambers <- intersect(flux_chambers, geom_chambers)

cat("OVERLAP ANALYSIS:\n")
cat("Common Sample IDs:", length(common_samples), "out of", length(flux_samples), "flux samples\n")
cat("Common Chamber IDs:", length(common_chambers), "out of", length(flux_chambers), "flux chambers\n")
cat("Common samples:", paste(head(sort(common_samples), 20), collapse = ", "), "\n")
cat("Common chambers:", paste(head(sort(common_chambers), 10), collapse = ", "), "\n\n")

# Try different join strategies

# Strategy 1: Sample ID only (ignore chamber mismatch)
geometry_by_sample <- geometry_clean %>%
  mutate(
    surface_area_cm2 = as.numeric(Sc_cm2_corrected),
    volume_cm3 = as.numeric(Vc_cm3_corrected)
  ) %>%
  group_by(Sample_ID) %>%
  summarise(
    surface_area_cm2 = mean(surface_area_cm2, na.rm = TRUE),
    volume_cm3 = mean(volume_cm3, na.rm = TRUE),
    n_chamber_configs = n_distinct(Chamber_ID),
    chambers_used = paste(unique(Chamber_ID), collapse = ", "),
    .groups = 'drop'
  )

join_sample_only <- flux_clean %>%
  left_join(geometry_by_sample, by = "Sample_ID") %>%
  mutate(has_geometry = !is.na(surface_area_cm2))

success_rate_sample <- round(sum(join_sample_only$has_geometry) / nrow(join_sample_only) * 100, 1)

cat("STRATEGY 1 - Sample ID only:\n")
cat("Success rate:", success_rate_sample, "%\n")
cat("Records with geometry:", sum(join_sample_only$has_geometry), "out of", nrow(join_sample_only), "\n\n")

# Strategy 2: Fuzzy matching on Sample ID
# Convert numeric sample IDs and try matching
flux_clean_numeric <- flux_clean %>%
  mutate(Sample_ID_numeric = as.numeric(Sample_ID))

geometry_clean_numeric <- geometry_clean %>%
  mutate(Sample_ID_numeric = as.numeric(Sample_ID))

geometry_numeric <- geometry_clean_numeric %>%
  filter(!is.na(Sample_ID_numeric)) %>%
  mutate(
    surface_area_cm2 = as.numeric(Sc_cm2_corrected),
    volume_cm3 = as.numeric(Vc_cm3_corrected)
  ) %>%
  group_by(Sample_ID_numeric) %>%
  summarise(
    surface_area_cm2 = mean(surface_area_cm2, na.rm = TRUE),
    volume_cm3 = mean(volume_cm3, na.rm = TRUE),
    .groups = 'drop'
  )

join_numeric <- flux_clean_numeric %>%
  filter(!is.na(Sample_ID_numeric)) %>%
  left_join(geometry_numeric, by = "Sample_ID_numeric") %>%
  mutate(has_geometry = !is.na(surface_area_cm2))

success_rate_numeric <- round(sum(join_numeric$has_geometry) / nrow(join_numeric) * 100, 1)

cat("STRATEGY 2 - Numeric Sample ID matching:\n")
cat("Success rate:", success_rate_numeric, "%\n")
cat("Records with geometry:", sum(join_numeric$has_geometry), "out of", nrow(join_numeric), "\n\n")

# Choose the best strategy
if (success_rate_sample >= success_rate_numeric) {
  final_joined <- join_sample_only
  best_method <- "Sample ID only"
  best_rate <- success_rate_sample
} else {
  final_joined <- join_numeric
  best_method <- "Numeric Sample ID"
  best_rate <- success_rate_numeric
}

cat("BEST STRATEGY:", best_method, "with", best_rate, "% success rate\n\n")

# Show successful joins
successful_joins <- final_joined %>%
  filter(has_geometry) %>%
  select(Date, Sample_ID, Chamber_ID, surface_area_cm2, volume_cm3)

cat("SUCCESSFUL JOINS (first 20):\n")
print(head(successful_joins, 20))

# Show failed joins for diagnosis
failed_joins <- final_joined %>%
  filter(!has_geometry) %>%
  select(Date, Sample_ID, Chamber_ID) %>%
  head(10)

cat("\nFAILED JOINS (sample):\n")
print(failed_joins)

# Export the best result
write_csv(final_joined, "flux_with_geometry_fixed.csv")
#cat("\nExported fixed join as 'flux_with_geometry_fixed.xlsx'\n")

# Summary of issues and recommendations
cat("\nDIAGNOSIS AND RECOMMENDATIONS:\n")
cat("================================\n")

if (best_rate < 50) {
  cat("⚠️  Low success rate suggests major identifier mismatches.\n")
  cat("Possible issues:\n")
  cat("- Sample IDs use different formats between datasets\n")
  cat("- Chamber IDs recorded inconsistently\n") 
  cat("- Temporal mismatch (flux data from different time period)\n")
  cat("- Data entry errors\n\n")
  cat("Recommendations:\n")
  cat("- Review sample ID formats in both datasets\n")
  cat("- Check if Plot Tag in flux data matches Sample in geometry data\n")
  cat("- Consider manual inspection of identifier patterns\n")
} else {
  cat("✅ Good success rate achieved with", best_method, "strategy\n")
  cat("Dataset ready for flux calculations\n")
}