# ==============================================================================
# RF Diagnostic Plots
# ==============================================================================
# Purpose: Diagnostic visualizations of RF model performance, feature
#   importance, and seasonal patterns.
#
# Pipeline stage: 3 — Upscaling
# Run after: 02_rf_models.R
#
# Inputs:
#   - model objects in memory (from 02_rf_models)
#
# Outputs:
#   - diagnostic PNGs (to outputs/figures/)
# ==============================================================================

# Fix 1: RF IMPORTANCE PLOTS WITH FEATURE NAMES
# Replace the existing QC5 section with this:

cat("\nGenerating feature importance plots with names...\n")

# Tree importance with proper names
tree_imp_raw <- importance(TreeRF)
tree_importance <- data.frame(
  feature = names(tree_imp_raw),  # This gets the actual feature names
  Importance = as.numeric(tree_imp_raw)
) %>%
  arrange(desc(Importance)) %>%
  head(20)  # Top 20 features

p5a <- ggplot(tree_importance, aes(x = reorder(feature, Importance), y = Importance)) +
  geom_col(fill = "darkblue") +
  coord_flip() +
  labs(title = "Tree Model - Top 20 Feature Importances",
       x = "Feature", y = "Importance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))  # Ensure text is readable

ggsave("../../outputs/figures/QC_TREE_IMPORTANCE.png", p5a, width = 10, height = 8)

# Soil importance with proper names
soil_imp_raw <- importance(SoilRF)
soil_importance <- data.frame(
  feature = names(soil_imp_raw),  # This gets the actual feature names
  Importance = as.numeric(soil_imp_raw)
) %>%
  arrange(desc(Importance))

p5b <- ggplot(soil_importance, aes(x = reorder(feature, Importance), y = Importance)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Soil Model - Feature Importances",
       x = "Feature", y = "Importance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

ggsave("../../outputs/figures/QC_SOIL_IMPORTANCE.png", p5b, width = 10, height = 6)

cat("✓ Feature importance plots with names saved\n")

# =============================================================================
# Fix 2: MONTHLY FLUX PLOT (nmol m-2 s-1)
# Add this as a new plot after your existing QC plots
# =============================================================================

cat("\nGenerating monthly flux predictions plot...\n")

# Extract the per-area fluxes from your results
# These are the mean fluxes before scaling by area fractions
monthly_flux_plot_data <- bind_rows(
  tree_results %>% 
    dplyr::select(month, mean_flux) %>%
    mutate(
      flux_nmol = mean_flux * 1000,  # Convert μmol to nmol
      type = "Tree stems"
    ),
  soil_results %>%
    dplyr::select(month, mean_flux) %>%
    mutate(
      flux_nmol = mean_flux * 1000,  # Convert μmol to nmol
      type = "Soil"
    )
)

# Create the monthly flux plot
p_monthly <- ggplot(monthly_flux_plot_data, aes(x = month, y = flux_nmol, color = type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = c("Tree stems" = "forestgreen", "Soil" = "brown")) +
  labs(title = "Monthly CH₄ Flux Predictions (per unit area)",
       x = "Month",
       y = expression(paste("CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
       color = "Source",
       caption = "Note: These are per-area fluxes, not scaled by plot fractions") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("../../outputs/figures/MONTHLY_FLUX_PREDICTIONS.png", p_monthly, width = 10, height = 6)

# Also create a version with both on same panel but different y-axes due to scale differences
p_monthly_dual <- ggplot(monthly_flux_plot_data, aes(x = month, y = flux_nmol)) +
  geom_line(aes(color = type), size = 1.2) +
  geom_point(aes(color = type), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~ type, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = c("Tree stems" = "forestgreen", "Soil" = "brown")) +
  labs(title = "Monthly CH₄ Flux Predictions by Source",
       x = "Month",
       y = expression(paste("CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
       caption = "Note: Different y-axis scales; soil shows uptake (negative), trees show emission (positive)") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave("../../outputs/figures/MONTHLY_FLUX_PREDICTIONS_DUAL.png", p_monthly_dual, width = 8, height = 10)

# Print summary statistics for the plot
cat("\nMonthly flux summary (nmol m-2 s-1):\n")
monthly_flux_summary <- monthly_flux_plot_data %>%
  group_by(type) %>%
  summarise(
    mean_flux = mean(flux_nmol),
    min_flux = min(flux_nmol),
    max_flux = max(flux_nmol),
    peak_month = month[which.max(abs(flux_nmol))],
    .groups = "drop"
  )
print(monthly_flux_summary)

cat("✓ Monthly flux plots saved\n")

# =============================================================================
# COMBINED OBSERVATION VS PREDICTION SCATTER PLOT
# Add this after your other QC plots
# =============================================================================

cat("\nGenerating combined observation vs prediction plot...\n")

# Combine tree and soil observations and predictions
combined_obs_pred <- bind_rows(
  # Tree data
  data.frame(
    observed = sinh(y_tree),  # Back-transform from asinh
    predicted = sinh(TreeRF$predictions),  # Back-transform predictions
    source = "Trees",
    chamber_type = tree_combined$I_monthly[complete_rows]  # Get chamber type for coloring
  ),
  # Soil data
  data.frame(
    observed = sinh(y_soil),  # Back-transform
    predicted = sinh(SoilRF$predictions),  # Back-transform
    source = "Soil",
    chamber_type = 2  # Dummy value for soil (for consistent coloring)
  )
) %>%
  mutate(
    # Convert to nmol for display
    observed_nmol = observed * 1000,
    predicted_nmol = predicted * 1000
  )

# Calculate overall R² and statistics
lm_combined <- lm(predicted_nmol ~ observed_nmol, data = combined_obs_pred)
r2_combined <- summary(lm_combined)$r.squared
rmse_combined <- sqrt(mean((combined_obs_pred$predicted_nmol - combined_obs_pred$observed_nmol)^2, na.rm = TRUE))

# Calculate R² by source for annotation
r2_by_source <- combined_obs_pred %>%
  group_by(source) %>%
  summarise(
    r2 = cor(observed_nmol, predicted_nmol, use = "complete.obs")^2,
    n = n(),
    .groups = "drop"
  )

# Create the combined scatter plot
p_combined <- ggplot(combined_obs_pred, aes(x = observed_nmol, y = predicted_nmol)) +
  # Points colored by source
  geom_point(aes(color = source, shape = source), alpha = 0.6, size = 2) +
  # Overall linear fit
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", size = 1) +
  # 1:1 line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  # Axis transformations for better visualization
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  # Colors and shapes
  scale_color_manual(values = c("Trees" = "forestgreen", "Soil" = "brown")) +
  scale_shape_manual(values = c("Trees" = 16, "Soil" = 17)) +
  # Labels
  labs(
    title = "Combined Observation vs Prediction: Trees and Soil",
    subtitle = bquote(Overall~R^2 == .(round(r2_combined, 3))~";"~RMSE == .(round(rmse_combined, 2))~"nmol m"^-2~"s"^-1),
    x = expression(paste("Observed CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
    color = "Source",
    shape = "Source"
  ) +
  # Add R² annotations for each source
  annotate("text", 
           x = min(combined_obs_pred$observed_nmol[combined_obs_pred$source == "Trees"], na.rm = TRUE),
           y = max(combined_obs_pred$predicted_nmol[combined_obs_pred$source == "Trees"], na.rm = TRUE),
           label = paste0("Trees R² = ", round(r2_by_source$r2[r2_by_source$source == "Trees"], 3)),
           hjust = 0, vjust = 1, color = "forestgreen", size = 4) +
  annotate("text",
           x = max(combined_obs_pred$observed_nmol[combined_obs_pred$source == "Soil"], na.rm = TRUE),
           y = min(combined_obs_pred$predicted_nmol[combined_obs_pred$source == "Soil"], na.rm = TRUE),
           label = paste0("Soil R² = ", round(r2_by_source$r2[r2_by_source$source == "Soil"], 3)),
           hjust = 1, vjust = 0, color = "brown", size = 4) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("../../outputs/figures/COMBINED_OBS_VS_PRED.png", p_combined, width = 9, height = 8)

# Alternative version with facets for clarity
p_combined_facet <- ggplot(combined_obs_pred, aes(x = observed_nmol, y = predicted_nmol)) +
  geom_point(aes(color = source), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ source, scales = "free", ncol = 2) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_color_manual(values = c("Trees" = "forestgreen", "Soil" = "brown")) +
  labs(
    title = "Observation vs Prediction by Source",
    x = expression(paste("Observed CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("../../outputs/figures/COMBINED_OBS_VS_PRED_FACET.png", p_combined_facet, width = 12, height = 6)

# Print summary statistics
cat("\nCombined Model Performance:\n")
cat("  Overall R²:", round(r2_combined, 3), "\n")
cat("  Overall RMSE:", round(rmse_combined, 2), "nmol/m²/s\n")
cat("  Tree R²:", round(r2_by_source$r2[r2_by_source$source == "Trees"], 3), 
    "(n =", r2_by_source$n[r2_by_source$source == "Trees"], ")\n")
cat("  Soil R²:", round(r2_by_source$r2[r2_by_source$source == "Soil"], 3),
    "(n =", r2_by_source$n[r2_by_source$source == "Soil"], ")\n")

# Calculate percentage of variance explained
cat("\nVariance explained:\n")
cat("  Trees:", round(100 * r2_by_source$r2[r2_by_source$source == "Trees"], 1), "%\n")
cat("  Soil:", round(100 * r2_by_source$r2[r2_by_source$source == "Soil"], 1), "%\n")
cat("  Combined:", round(100 * r2_combined, 1), "%\n")

cat("✓ Combined observation vs prediction plots saved\n")



# =============================================================================
# DIAGNOSTIC PLOTS: FIXED VERSION
# Handles tree ID issue and moisture binning error
# =============================================================================

library(ggridges)  # For ridge plots - install.packages("ggridges") if needed

cat("\n=== DIAGNOSTIC DISTRIBUTIONS (FIXED) ===\n")

# =============================================================================
# 1. TREE FLUX DISTRIBUTIONS BY SPECIES (EXCLUDING TREE IDS)
# =============================================================================

cat("\nAnalyzing tree flux distributions by species...\n")

# Combine observed and predicted values by species
tree_species_data <- data.frame(
  species = tree_train$species[complete_rows],
  observed = sinh(y_tree) * 1000,  # Convert to nmol
  predicted = sinh(TreeRF$predictions) * 1000,  # Convert to nmol
  stringsAsFactors = FALSE
) %>%
  # IMPORTANT: Filter out entries that are tree IDs, not species names
  filter(!is.na(species) & 
           species != "" & 
           !grepl("^Tree [0-9]", species) &  # Remove "Tree 1234" style entries
           !grepl("^[0-9]+-[0-9]+", species))  # Remove "3-1" style entries

# Get actual species with enough observations
species_counts <- table(tree_species_data$species)
common_species <- names(species_counts[species_counts >= 5])  # At least 5 observations

cat("Found", length(common_species), "species with ≥5 observations\n")
cat("Species:", paste(common_species, collapse = ", "), "\n")

# Filter to common species
tree_species_filtered <- tree_species_data %>%
  filter(species %in% common_species)

if(nrow(tree_species_filtered) > 0) {
  
  # Reshape for plotting
  tree_species_long <- tree_species_filtered %>%
    pivot_longer(cols = c(observed, predicted),
                 names_to = "type",
                 values_to = "flux_nmol") %>%
    mutate(type = factor(type, levels = c("observed", "predicted")))
  
  # Calculate summary stats by species
  species_summary <- tree_species_long %>%
    group_by(species, type) %>%
    summarise(
      mean = mean(flux_nmol, na.rm = TRUE),
      median = median(flux_nmol, na.rm = TRUE),
      sd = sd(flux_nmol, na.rm = TRUE),
      n = n()/2,  # Divide by 2 because we doubled rows in pivot_longer
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = type, 
                values_from = c(mean, median, sd),
                names_sep = "_")
  
  cat("\nSpecies-level statistics (nmol/m²/s):\n")
  print(species_summary)
  
  # PLOT 1: Box plots by species
  p_species_box <- ggplot(tree_species_long, aes(x = species, y = flux_nmol, fill = type)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
    scale_fill_manual(values = c("observed" = "skyblue", "predicted" = "coral"),
                      labels = c("Observed", "Predicted")) +
    labs(title = "Tree CH₄ Flux Distribution by Species",
         subtitle = "Showing only identified species (tree IDs excluded)",
         x = "Species",
         y = expression(paste("CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
         fill = "Data Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "bottom") +
    coord_cartesian(ylim = quantile(tree_species_long$flux_nmol, c(0.01, 0.99), na.rm = TRUE))
  
  ggsave("../../outputs/figures/DIAGNOSTIC_SPECIES_BOXPLOT.png", p_species_box, width = 12, height = 7)
  
  # PLOT 2: Scatter plot of observed vs predicted by species
  p_species_scatter <- ggplot(tree_species_filtered, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = FALSE, size = 0.8, color = "blue") +
    facet_wrap(~ species, scales = "free", ncol = 3) +
    labs(title = "Observed vs Predicted by Species",
         x = expression(paste("Observed CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
         y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))) +
    theme_minimal()
  
  ggsave("../../outputs/figures/DIAGNOSTIC_SPECIES_SCATTER.png", p_species_scatter, width = 12, height = 10)
  
} else {
  cat("\nWARNING: No valid species found after filtering out tree IDs\n")
  cat("This suggests species information may not be properly linked\n")
}

# =============================================================================
# 2. FLUX DISTRIBUTIONS BY SOIL MOISTURE (FIXED BINNING)
# =============================================================================

cat("\nAnalyzing flux distributions by soil moisture...\n")

# Combine tree and soil data for moisture analysis
moisture_data <- bind_rows(
  # Tree data
  data.frame(
    source = "Trees",
    moisture = tree_train$soil_moisture_at_tree[complete_rows],
    observed = sinh(y_tree) * 1000,  # nmol
    predicted = sinh(TreeRF$predictions) * 1000
  ),
  # Soil data
  data.frame(
    source = "Soil",
    moisture = soil_train$soil_moisture_at_site[complete_rows_soil],
    observed = sinh(y_soil) * 1000,
    predicted = sinh(SoilRF$predictions) * 1000
  )
) %>%
  filter(!is.na(moisture))

# Check moisture distribution
cat("\nMoisture range:", range(moisture_data$moisture), "\n")
cat("Unique moisture values:", length(unique(moisture_data$moisture)), "\n")

# Create moisture bins - handle case where there aren't enough unique values
if(length(unique(moisture_data$moisture)) >= 10) {
  # Use deciles if we have enough unique values
  moisture_breaks <- unique(quantile(moisture_data$moisture, probs = seq(0, 1, by = 0.1), na.rm = TRUE))
  n_bins <- length(moisture_breaks) - 1
  bin_labels <- paste0("D", 1:n_bins)
} else {
  # Use fewer bins based on unique values
  n_bins <- min(5, length(unique(moisture_data$moisture)))
  moisture_breaks <- unique(quantile(moisture_data$moisture, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  bin_labels <- paste0("B", 1:n_bins)
  cat("Using", n_bins, "bins due to limited moisture variation\n")
}

moisture_data <- moisture_data %>%
  mutate(
    moisture_bin = cut(moisture,
                       breaks = moisture_breaks,
                       labels = bin_labels,
                       include.lowest = TRUE)
  )

# Calculate summary stats by moisture bin
moisture_summary <- moisture_data %>%
  group_by(source, moisture_bin) %>%
  summarise(
    moisture_mean = mean(moisture, na.rm = TRUE),
    moisture_min = min(moisture, na.rm = TRUE),
    moisture_max = max(moisture, na.rm = TRUE),
    obs_mean = mean(observed, na.rm = TRUE),
    obs_sd = sd(observed, na.rm = TRUE),
    pred_mean = mean(predicted, na.rm = TRUE),
    pred_sd = sd(predicted, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(moisture_bin))  # Remove any NA bins

cat("\nMoisture bin statistics:\n")
print(moisture_summary)

# Reshape for plotting
moisture_long <- moisture_data %>%
  filter(!is.na(moisture_bin)) %>%
  pivot_longer(cols = c(observed, predicted),
               names_to = "type",
               values_to = "flux_nmol") %>%
  mutate(type = factor(type, levels = c("observed", "predicted")))

# PLOT 3: Box plots by moisture bin
p_moisture_box <- ggplot(moisture_long, aes(x = moisture_bin, y = flux_nmol, fill = type)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
  facet_wrap(~ source, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("observed" = "lightblue", "predicted" = "orange"),
                    labels = c("Observed", "Predicted")) +
  labs(title = "CH₄ Flux Distribution by Soil Moisture",
       subtitle = paste("Binned into", n_bins, "groups from dry to wet"),
       x = "Moisture Bin",
       y = expression(paste("CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
       fill = "Data Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("../../outputs/figures/DIAGNOSTIC_MOISTURE_BOXPLOT.png", p_moisture_box, width = 10, height = 10)

# PLOT 4: Mean flux vs moisture bin with error bars
p_moisture_means <- ggplot(moisture_summary, aes(x = moisture_bin)) +
  geom_point(aes(y = obs_mean, color = "Observed"), size = 3) +
  geom_errorbar(aes(ymin = obs_mean - obs_sd, ymax = obs_mean + obs_sd, color = "Observed"),
                width = 0.2, alpha = 0.5) +
  geom_point(aes(y = pred_mean, color = "Predicted"), size = 3) +
  geom_errorbar(aes(ymin = pred_mean - pred_sd, ymax = pred_mean + pred_sd, color = "Predicted"),
                width = 0.2, alpha = 0.5) +
  geom_line(aes(y = obs_mean, color = "Observed", group = 1), size = 1) +
  geom_line(aes(y = pred_mean, color = "Predicted", group = 1), size = 1, linetype = "dashed") +
  facet_wrap(~ source, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  labs(title = "Mean CH₄ Flux by Moisture Level",
       subtitle = "Error bars show ±1 SD",
       x = paste("Moisture Bin (", bin_labels[1], "= driest,", bin_labels[n_bins], "= wettest)"),
       y = expression(paste("Mean CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
       color = "Data Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("../../outputs/figures/DIAGNOSTIC_MOISTURE_MEANS.png", p_moisture_means, width = 10, height = 8)

# PLOT 5: Continuous moisture relationship
p_moisture_continuous <- ggplot(moisture_data, aes(x = moisture)) +
  geom_point(aes(y = observed, color = "Observed"), alpha = 0.5, size = 1.5) +
  geom_point(aes(y = predicted, color = "Predicted"), alpha = 0.5, size = 1.5) +
  geom_smooth(aes(y = observed, color = "Observed"), method = "loess", se = TRUE, size = 1) +
  geom_smooth(aes(y = predicted, color = "Predicted"), method = "loess", se = TRUE, size = 1) +
  facet_wrap(~ source, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  labs(title = "CH₄ Flux vs Soil Moisture (continuous)",
       x = "Soil Moisture (m³/m³)",
       y = expression(paste("CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
       color = "Data Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("../../outputs/figures/DIAGNOSTIC_MOISTURE_CONTINUOUS.png", p_moisture_continuous, width = 10, height = 8)

# Calculate correlation by groups
cat("\nCorrelation analysis:\n")

# By species (if we have valid species)
if(nrow(tree_species_filtered) > 0) {
  species_cor <- tree_species_filtered %>%
    group_by(species) %>%
    summarise(
      correlation = cor(observed, predicted, use = "complete.obs"),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(correlation))
  
  cat("\nCorrelation by species:\n")
  print(species_cor)
}

# By moisture bin
moisture_cor <- moisture_data %>%
  filter(!is.na(moisture_bin)) %>%
  group_by(source, moisture_bin) %>%
  summarise(
    correlation = cor(observed, predicted, use = "complete.obs"),
    n = n(),
    moisture_range = paste0("[", round(min(moisture), 3), "-", round(max(moisture), 3), "]"),
    .groups = "drop"
  )

cat("\nCorrelation by moisture bin:\n")
print(moisture_cor)

cat("\n✓ Diagnostic distribution plots complete\n")

# Check if training data actually shows temp relationship
temp_flux_relationship <- tree_train %>%
  mutate(temp_bin = cut(air_temp_C_mean, breaks = 5)) %>%
  group_by(temp_bin) %>%
  summarise(
    n = n(),
    mean_flux = mean(stem_flux_corrected),
    median_flux = median(stem_flux_corrected),
    .groups = "drop"
  )

print(temp_flux_relationship)

# Plot it
library(ggplot2)
ggplot(tree_train, aes(x = air_temp_C_mean, y = stem_flux_corrected)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "red") +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  labs(title = "Training Data: Temperature vs Flux",
       subtitle = "Does flux increase with temperature?")

