# Data Completeness Analysis
# Calculate percentage of non-missing data for each tree and overall

library(dplyr)

# Data Completeness Analysis
# Calculate percentage of non-missing data for each tree and overall

library(dplyr)

# Calculate completeness for each row (tree)
completeness_by_tree <- merged_data_final %>%
  mutate(
    total_columns = ncol(.) - 1,  # Exclude tree_id column
    non_missing = rowSums(!is.na(select(., -tree_id))),
    percent_complete = round((non_missing / total_columns) * 100, 1)
  ) %>%
  select(tree_id, total_columns, non_missing, percent_complete) %>%
  arrange(desc(percent_complete))

# Display summary statistics
cat("=== DATA COMPLETENESS SUMMARY ===\n")
cat("Total trees:", nrow(completeness_by_tree), "\n")
cat("Total possible data points:", nrow(merged_data_final) * (ncol(merged_data_final) - 1), "\n")
cat("Actual data points:", sum(completeness_by_tree$non_missing), "\n")
cat("Overall completeness:", round(mean(completeness_by_tree$percent_complete), 1), "%\n\n")

cat("Completeness distribution:\n")
cat("Min:", min(completeness_by_tree$percent_complete), "%\n")
cat("25th percentile:", quantile(completeness_by_tree$percent_complete, 0.25), "%\n")
cat("Median:", median(completeness_by_tree$percent_complete), "%\n")
cat("75th percentile:", quantile(completeness_by_tree$percent_complete, 0.75), "%\n")
cat("Max:", max(completeness_by_tree$percent_complete), "%\n\n")

# Show trees with highest completeness
cat("=== TOP 10 MOST COMPLETE TREES ===\n")
print(head(completeness_by_tree, 10))

cat("\n=== BOTTOM 10 LEAST COMPLETE TREES ===\n")
print(tail(completeness_by_tree, 10))

# Completeness by ranges
completeness_ranges <- completeness_by_tree %>%
  mutate(
    completeness_range = case_when(
      percent_complete >= 80 ~ "80-100%",
      percent_complete >= 60 ~ "60-79%", 
      percent_complete >= 40 ~ "40-59%",
      percent_complete >= 20 ~ "20-39%",
      TRUE ~ "0-19%"
    )
  ) %>%
  count(completeness_range, name = "n_trees") %>%
  mutate(percentage = round((n_trees / nrow(completeness_by_tree)) * 100, 1))

cat("\n=== TREES BY COMPLETENESS RANGES ===\n")
print(completeness_ranges)

# Completeness by data type
data_type_completeness <- merged_data_final %>%
  summarise(
    dbh_complete = sum(!is.na(dbh)),
    gas_complete = sum(!is.na(CO2_concentration)),
    flux_complete = sum(!is.na(CO2_best.flux_50cm)),
    soil_complete = sum(!is.na(VWC_mean)),
    wood_complete = sum(!is.na(outer_density_final)),
    ddpcr_complete = sum(!is.na(ddpcr_mcra_probe_Mineral_loose)),
    total_trees = n()
  ) %>%
  pivot_longer(
    cols = ends_with("_complete"),
    names_to = "data_type",
    values_to = "n_trees"
  ) %>%
  mutate(
    data_type = gsub("_complete", "", data_type),
    percentage = round((n_trees / total_trees) * 100, 1)
  ) %>%
  arrange(desc(percentage))

cat("\n=== COMPLETENESS BY DATA TYPE ===\n")
print(data_type_completeness)

# Trees with all major data types
trees_with_all_data <- merged_data_final %>%
  filter(
    !is.na(dbh),
    !is.na(CO2_concentration),
    !is.na(CO2_best.flux_50cm),
    !is.na(VWC_mean),
    !is.na(outer_density_final),
    !is.na(ddpcr_mcra_probe_Mineral_loose)
  )

cat("\n=== TREES WITH ALL MAJOR DATA TYPES ===\n")
cat("Trees with complete data across all 6 major datasets:", nrow(trees_with_all_data), "\n")
cat("Percentage of total trees:", round((nrow(trees_with_all_data) / nrow(merged_data_final)) * 100, 1), "%\n")

if (nrow(trees_with_all_data) > 0) {
  cat("\nTrees with complete data:\n")
  print(trees_with_all_data$tree_id)
}

# Check for trees with very low completeness (potential data issues)
low_completeness_trees <- completeness_by_tree %>%
  filter(percent_complete < 10)

if (nrow(low_completeness_trees) > 0) {
  cat("\n=== TREES WITH <10% COMPLETENESS (POTENTIAL ISSUES) ===\n")
  print(low_completeness_trees)
}