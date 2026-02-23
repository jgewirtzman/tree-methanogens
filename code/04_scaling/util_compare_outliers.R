# ==============================================================================
# Compare Outlier Removal Strategies
# ==============================================================================
# Purpose: Diagnostic comparison of model predictions with and without outlier
#   removal.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(ranger)
library(patchwork)

cat("\n=== OUTLIER REMOVAL COMPARISON DIAGNOSTIC ===\n")
cat("Comparing three strategies: None, 1% Tails, MAD (k=8)\n\n")

# -----------------------------------------------------------------------------
# SETUP: Load data and define functions
# -----------------------------------------------------------------------------

# Load the original data
load("../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")

# Extract and prepare base datasets with unit conversion
TREE_JULY_BASE <- rf_workflow_data$PLACEHOLDER_TREE_JULY %>%
  mutate(stem_flux_umol_m2_s = stem_flux_umol_m2_s / 1000)  # Convert nmol to μmol

TREE_YEAR_BASE <- rf_workflow_data$PLACEHOLDER_TREE_YEAR %>%
  mutate(stem_flux_umol_m2_s = stem_flux_umol_m2_s / 1000)

SOIL_YEAR_BASE <- rf_workflow_data$PLACEHOLDER_SOIL_YEAR %>%
  mutate(soil_flux_umol_m2_s = soil_flux_umol_m2_s / 1000)

# MAD-based outlier removal function
remove_outliers_mad <- function(x, k = 8) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  
  if(mad_val == 0 || is.na(mad_val)) {
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    
    if(iqr == 0) return(!is.na(x))
    
    lower <- q1 - 3 * iqr
    upper <- q3 + 3 * iqr
    return(x >= lower & x <= upper)
  } else {
    lower <- med - k * mad_val
    upper <- med + k * mad_val
    return(x >= lower & x <= upper)
  }
}

# -----------------------------------------------------------------------------
# FUNCTION: Process data with specified outlier method
# -----------------------------------------------------------------------------

process_with_outlier_method <- function(method = c("none", "percentile", "mad")) {
  
  cat("Processing with method:", method, "\n")
  
  # Combine tree data
  tree_combined <- bind_rows(
    TREE_JULY_BASE %>% mutate(chamber_type = "rigid"),
    TREE_YEAR_BASE %>% mutate(chamber_type = coalesce(chamber_type, "semirigid"))
  ) %>%
    filter(!is.na(stem_flux_umol_m2_s))
  
  # Transform to asinh scale
  tree_combined$z <- asinh(tree_combined$stem_flux_umol_m2_s)
  
  # Apply outlier removal based on method
  if(method == "percentile") {
    # Remove top and bottom 1%
    tree_percentiles <- quantile(tree_combined$z, probs = c(0.01, 0.99), na.rm = TRUE)
    outlier_mask_tree <- tree_combined$z >= tree_percentiles[1] & 
      tree_combined$z <= tree_percentiles[2]
    tree_combined <- tree_combined[outlier_mask_tree, ]
    
    # Same for soil
    SOIL_YEAR <- SOIL_YEAR_BASE
    SOIL_YEAR$z_soil <- asinh(SOIL_YEAR$soil_flux_umol_m2_s)
    soil_percentiles <- quantile(SOIL_YEAR$z_soil, probs = c(0.01, 0.99), na.rm = TRUE)
    outlier_mask_soil <- SOIL_YEAR$z_soil >= soil_percentiles[1] & 
      SOIL_YEAR$z_soil <= soil_percentiles[2]
    SOIL_YEAR <- SOIL_YEAR[outlier_mask_soil, ]
    
  } else if(method == "mad") {
    # MAD-based removal
    outlier_mask_tree <- remove_outliers_mad(tree_combined$z, k = 8)
    tree_combined <- tree_combined[outlier_mask_tree, ]
    
    SOIL_YEAR <- SOIL_YEAR_BASE
    SOIL_YEAR$z_soil <- asinh(SOIL_YEAR$soil_flux_umol_m2_s)
    outlier_mask_soil <- remove_outliers_mad(SOIL_YEAR$z_soil, k = 8)
    SOIL_YEAR <- SOIL_YEAR[outlier_mask_soil, ]
    
  } else {
    # No outlier removal
    SOIL_YEAR <- SOIL_YEAR_BASE
  }
  
  # Clean up temporary columns
  if("z_soil" %in% names(SOIL_YEAR)) {
    SOIL_YEAR$z_soil <- NULL
  }
  
  n_tree <- nrow(tree_combined)
  n_soil <- nrow(SOIL_YEAR)
  
  # Return processed datasets
  return(list(
    tree_data = tree_combined,
    soil_data = SOIL_YEAR,
    n_tree = n_tree,
    n_soil = n_soil,
    method = method
  ))
}

# -----------------------------------------------------------------------------
# FUNCTION: Simplified RF training and monthly prediction
# -----------------------------------------------------------------------------

train_and_predict_monthly <- function(processed_data) {
  
  tree_data <- processed_data$tree_data
  soil_data <- processed_data$soil_data
  
  # Simple feature engineering for trees
  tree_train <- tree_data %>%
    mutate(
      stem_flux_corrected = stem_flux_umol_m2_s,
      y_asinh = asinh(stem_flux_corrected),
      month_sin = sin(2 * pi * month / 12),
      month_cos = cos(2 * pi * month / 12),
      chamber_rigid = as.numeric(chamber_type == "rigid"),
      dbh_m = coalesce(dbh_m, median(dbh_m, na.rm = TRUE)),
      air_temp_C = coalesce(air_temp_C, median(air_temp_C, na.rm = TRUE)),
      soil_temp_C = coalesce(soil_temp_C, median(soil_temp_C, na.rm = TRUE)),
      soil_moisture_abs = coalesce(soil_moisture_abs, median(soil_moisture_abs, na.rm = TRUE))
    )
  
  # Tree features
  X_tree <- tree_train %>%
    select(dbh_m, air_temp_C, soil_temp_C, soil_moisture_abs, 
           month_sin, month_cos, chamber_rigid) %>%
    as.matrix()
  
  y_tree <- tree_train$y_asinh
  
  # Train tree RF
  tree_rf <- ranger(
    x = X_tree,
    y = y_tree,
    num.trees = 200,
    min.node.size = 5,
    mtry = floor(sqrt(ncol(X_tree))),
    num.threads = 1
  )
  
  # Simple feature engineering for soil
  soil_train <- soil_data %>%
    mutate(
      y_asinh = asinh(soil_flux_umol_m2_s),
      month_sin = sin(2 * pi * month / 12),
      month_cos = cos(2 * pi * month / 12),
      air_temp_C = coalesce(air_temp_C, median(air_temp_C, na.rm = TRUE)),
      soil_temp_C = coalesce(soil_temp_C, median(soil_temp_C, na.rm = TRUE)),
      soil_moisture_abs = coalesce(soil_moisture_abs, median(soil_moisture_abs, na.rm = TRUE))
    )
  
  # Soil features
  X_soil <- soil_train %>%
    select(air_temp_C, soil_temp_C, soil_moisture_abs, month_sin, month_cos) %>%
    as.matrix()
  
  y_soil <- soil_train$y_asinh
  
  # Train soil RF
  soil_rf <- ranger(
    x = X_soil,
    y = y_soil,
    num.trees = 200,
    min.node.size = 5,
    mtry = floor(sqrt(ncol(X_soil))),
    num.threads = 1
  )
  
  # Monthly predictions (simplified - using training data means per month)
  monthly_results <- map_df(1:12, function(m) {
    
    # Get monthly statistics from training data
    tree_month <- tree_train %>% filter(month == m)
    soil_month <- soil_train %>% filter(month == m)
    
    # Use mean of predictions for the month (simplified approach)
    tree_flux_mean <- ifelse(nrow(tree_month) > 0,
                             mean(sinh(tree_rf$predictions[tree_train$month == m]), na.rm = TRUE),
                             0)
    
    soil_flux_mean <- ifelse(nrow(soil_month) > 0,
                             mean(sinh(soil_rf$predictions[soil_train$month == m]), na.rm = TRUE),
                             0)
    
    tibble(
      month = m,
      tree_flux_umol_m2_s = tree_flux_mean,
      soil_flux_umol_m2_s = soil_flux_mean,
      tree_n = nrow(tree_month),
      soil_n = nrow(soil_month)
    )
  })
  
  # Add model performance metrics
  monthly_results$tree_r2 <- tree_rf$r.squared
  monthly_results$soil_r2 <- soil_rf$r.squared
  monthly_results$method <- processed_data$method
  
  return(monthly_results)
}

# -----------------------------------------------------------------------------
# MAIN COMPARISON
# -----------------------------------------------------------------------------

# Process data with each method
results_none <- process_with_outlier_method("none")
results_pct <- process_with_outlier_method("percentile")
results_mad <- process_with_outlier_method("mad")

# Get monthly predictions for each
monthly_none <- train_and_predict_monthly(results_none)
monthly_pct <- train_and_predict_monthly(results_pct)
monthly_mad <- train_and_predict_monthly(results_mad)

# Combine results
all_monthly <- bind_rows(
  monthly_none %>% mutate(method = "No Removal"),
  monthly_pct %>% mutate(method = "1% Tails"),
  monthly_mad %>% mutate(method = "MAD (k=8)")
)

# -----------------------------------------------------------------------------
# VISUALIZATION
# -----------------------------------------------------------------------------

# Plot 1: Tree flux comparison by month
p1 <- ggplot(all_monthly, aes(x = month, y = tree_flux_umol_m2_s * 1000, 
                              color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Tree CH₄ Flux: Impact of Outlier Removal",
       x = "Month",
       y = expression(paste("Mean Flux (nmol m"^-2, " s"^-1, ")")),
       color = "Outlier Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plot 2: Soil flux comparison by month
p2 <- ggplot(all_monthly, aes(x = month, y = soil_flux_umol_m2_s * 1000, 
                              color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12) +
  labs(title = "Soil CH₄ Flux: Impact of Outlier Removal",
       x = "Month",
       y = expression(paste("Mean Flux (nmol m"^-2, " s"^-1, ")")),
       color = "Outlier Method") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine plots
combined_plot <- p1 / p2
ggsave("../../outputs/figures/DIAGNOSTIC_outlier_comparison_monthly.png", combined_plot, 
       width = 10, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# SUMMARY STATISTICS
# -----------------------------------------------------------------------------

cat("\n=== OUTLIER REMOVAL IMPACT SUMMARY ===\n\n")

# Data retention
cat("DATA RETENTION:\n")
cat("Method         | Tree Obs | Soil Obs\n")
cat("---------------|----------|----------\n")
cat(sprintf("No Removal     | %8d | %8d\n", results_none$n_tree, results_none$n_soil))
cat(sprintf("1%% Tails      | %8d | %8d\n", results_pct$n_tree, results_pct$n_soil))
cat(sprintf("MAD (k=8)      | %8d | %8d\n", results_mad$n_tree, results_mad$n_soil))

cat("\nPERCENT RETAINED:\n")
cat(sprintf("1%% Tails:  Trees = %.1f%%, Soil = %.1f%%\n", 
            100 * results_pct$n_tree / results_none$n_tree,
            100 * results_pct$n_soil / results_none$n_soil))
cat(sprintf("MAD (k=8): Trees = %.1f%%, Soil = %.1f%%\n",
            100 * results_mad$n_tree / results_none$n_tree,
            100 * results_mad$n_soil / results_none$n_soil))

# Annual totals comparison
annual_comparison <- all_monthly %>%
  group_by(method) %>%
  summarise(
    annual_tree_nmol = sum(tree_flux_umol_m2_s * 1000) * 30.4,
    annual_soil_nmol = sum(soil_flux_umol_m2_s * 1000) * 30.4,
    tree_r2 = first(tree_r2),
    soil_r2 = first(soil_r2),
    .groups = "drop"
  )

cat("\nANNUAL TOTALS (nmol m⁻² yr⁻¹):\n")
print(annual_comparison)

# Relative differences
baseline <- annual_comparison %>% filter(method == "No Removal")
cat("\nRELATIVE DIFFERENCE FROM NO REMOVAL:\n")
for(m in c("1% Tails", "MAD (k=8)")) {
  comp <- annual_comparison %>% filter(method == m)
  tree_diff <- 100 * (comp$annual_tree_nmol - baseline$annual_tree_nmol) / baseline$annual_tree_nmol
  soil_diff <- 100 * (comp$annual_soil_nmol - baseline$annual_soil_nmol) / baseline$annual_soil_nmol
  cat(sprintf("%s: Tree = %+.1f%%, Soil = %+.1f%%\n", m, tree_diff, soil_diff))
}

# Model performance comparison
cat("\nMODEL PERFORMANCE (R²):\n")
cat("Method         | Tree R² | Soil R²\n")
cat("---------------|---------|--------\n")
for(m in unique(annual_comparison$method)) {
  row <- annual_comparison %>% filter(method == m)
  cat(sprintf("%-14s | %7.3f | %7.3f\n", m, row$tree_r2, row$soil_r2))
}

# Create difference plot
all_monthly_wide <- all_monthly %>%
  select(month, method, tree_flux_umol_m2_s, soil_flux_umol_m2_s) %>%
  pivot_wider(names_from = method, values_from = c(tree_flux_umol_m2_s, soil_flux_umol_m2_s))

# Calculate differences from no removal baseline
diff_data <- all_monthly %>%
  filter(method != "No Removal") %>%
  left_join(
    all_monthly %>% 
      filter(method == "No Removal") %>%
      select(month, tree_base = tree_flux_umol_m2_s, soil_base = soil_flux_umol_m2_s),
    by = "month"
  ) %>%
  mutate(
    tree_diff_pct = 100 * (tree_flux_umol_m2_s - tree_base) / tree_base,
    soil_diff_pct = 100 * (soil_flux_umol_m2_s - soil_base) / soil_base
  )

# Plot percent differences
p3 <- ggplot(diff_data, aes(x = month)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(aes(y = tree_diff_pct, color = method, linetype = "Tree"), size = 1) +
  geom_line(aes(y = soil_diff_pct, color = method, linetype = "Soil"), size = 1) +
  geom_point(aes(y = tree_diff_pct, color = method), size = 2) +
  geom_point(aes(y = soil_diff_pct, color = method), size = 2, shape = 17) +
  scale_x_continuous(breaks = 1:12) +
  scale_linetype_manual(values = c("Tree" = "solid", "Soil" = "dashed")) +
  labs(title = "Percent Change from No Outlier Removal",
       x = "Month",
       y = "Percent Change (%)",
       color = "Method",
       linetype = "Flux Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("../../outputs/figures/DIAGNOSTIC_outlier_percent_change.png", p3, width = 10, height = 6, dpi = 300)

# Save detailed results
write.csv(all_monthly, "../../outputs/tables/DIAGNOSTIC_monthly_outlier_comparison.csv", row.names = FALSE)
write.csv(annual_comparison, "../../outputs/tables/DIAGNOSTIC_annual_outlier_comparison.csv", row.names = FALSE)

cat("\n✓ Diagnostic comparison complete\n")
cat("  Saved: DIAGNOSTIC_outlier_comparison_monthly.png\n")
cat("  Saved: DIAGNOSTIC_outlier_percent_change.png\n")
cat("  Saved: DIAGNOSTIC_monthly_outlier_comparison.csv\n")
cat("  Saved: DIAGNOSTIC_annual_outlier_comparison.csv\n")