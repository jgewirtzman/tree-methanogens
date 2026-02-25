# ==============================================================================
# RF Publication Plots (Figures 1, 9, S6)
# ==============================================================================
# Purpose: Publication-quality RF model performance and seasonal prediction
#   figures.
#
# Pipeline stage: 4 — Visualization
# Run after: 02_rf_models.R
#
# Outputs:
#   - figure PDFs/PNGs (to outputs/figures/)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)  # For combining plots - install if needed
library(scales)
library(ranger)

## Load required data from upstream pipeline
if (!file.exists("outputs/models/RF_MODELS.RData") || !file.exists("outputs/models/TRAINING_DATA.RData")) {
  stop("Required data files not found. Run 02_rf_models.R first to generate:\n",
       "  - outputs/models/RF_MODELS.RData\n",
       "  - outputs/models/TRAINING_DATA.RData")
}
load("outputs/models/RF_MODELS.RData")           # TreeRF, SoilRF
load("outputs/models/TRAINING_DATA.RData")        # tree_train_complete, X_tree, X_soil, soil_train_complete

## Load monthly predictions (for panels g/h and Figure 2)
load("outputs/models/tree_monthly_predictions.RData")
tree_monthly_raw <- monthly_predictions
load("outputs/models/soil_monthly_predictions.RData")
soil_monthly_raw <- monthly_predictions

## Summarise by month for tree_results / soil_results
tree_results <- tree_monthly_raw %>%
  group_by(month) %>%
  summarise(mean_flux = mean(flux_umol_m2_s, na.rm = TRUE), .groups = "drop")

soil_results <- soil_monthly_raw %>%
  group_by(month) %>%
  summarise(mean_flux = mean(flux_umol_m2_s, na.rm = TRUE), .groups = "drop")

## Build monthly_results (combined plot-level) for Figure 2
monthly_results <- tree_results %>%
  rename(Phi_tree_umol_m2_s = mean_flux) %>%
  left_join(
    soil_results %>% rename(Phi_soil_umol_m2_s = mean_flux),
    by = "month"
  ) %>%
  mutate(Phi_plot_umol_m2_s = Phi_tree_umol_m2_s + Phi_soil_umol_m2_s)

cat("\n=== GENERATING PUBLICATION FIGURES ===\n")

# Set publication theme
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.position = "bottom",
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold")
  )

# =============================================================================
# FIGURE 1: MODEL PERFORMANCE (MAIN RESULT)
# Combined obs vs pred for trees and soil
# =============================================================================

cat("\nCreating Figure 1: Model Performance...\n")

# Prepare data
perf_data <- bind_rows(
  tree_train_complete %>%
    mutate(
      observed_nmol = stem_flux_corrected * 1000,
      predicted_nmol = pred_flux * 1000,
      source = "Tree stems"
    ) %>%
    dplyr::select(observed_nmol, predicted_nmol, source),
  
  soil_train_complete %>%
    mutate(
      observed_nmol = soil_flux_umol_m2_s * 1000,
      predicted_nmol = pred_flux * 1000,
      source = "Soil"
    ) %>%
    dplyr::select(observed_nmol, predicted_nmol, source)
)

# Calculate statistics for each panel
stats_by_source <- perf_data %>%
  group_by(source) %>%
  summarise(
    r2 = cor(observed_nmol, predicted_nmol, use = "complete.obs")^2,
    rmse = sqrt(mean((predicted_nmol - observed_nmol)^2, na.rm = TRUE)),
    n = n(),
    .groups = "drop"
  )

# Create labels with statistics
stats_labels <- stats_by_source %>%
  mutate(
    label = sprintf("R² = %.3f\nRMSE = %.1f\nn = %d", r2, rmse, n),
    x = ifelse(source == "Tree stems", 
               min(perf_data$observed_nmol[perf_data$source == "Tree stems"], na.rm=TRUE),
               max(perf_data$observed_nmol[perf_data$source == "Soil"], na.rm=TRUE) * 0.5),
    y = ifelse(source == "Tree stems",
               max(perf_data$predicted_nmol[perf_data$source == "Tree stems"], na.rm=TRUE) * 0.8,
               max(perf_data$predicted_nmol[perf_data$source == "Soil"], na.rm=TRUE) * 0.8)
  )

fig1 <- ggplot(perf_data, aes(x = observed_nmol, y = predicted_nmol)) +
  geom_point(aes(color = source), alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.8, alpha = 0.2) +
  geom_text(data = stats_labels, aes(x = x, y = y, label = label),
            hjust = 0, vjust = 1, size = 3.5, fontface = "bold") +
  facet_wrap(~ source, scales = "free", ncol = 2) +
  scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_color_manual(values = c("Tree stems" = "#2E7D32", "Soil" = "#8B4513")) +
  labs(
    x = expression(paste("Observed CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))
  ) +
  theme_pub +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"))

# ggsave("outputs/figures/Figure1_Model_Performance.pdf", fig1, width = 7, height = 3.5, dpi = 300)
# ggsave("outputs/figures/Figure1_Model_Performance.png", fig1, width = 7, height = 3.5, dpi = 300)

# =============================================================================
# FIGURE 2: SEASONAL PATTERNS
# Monthly predictions showing tree vs soil contributions
# =============================================================================

cat("\nCreating Figure 2: Seasonal Patterns...\n")

# Prepare monthly data
monthly_plot_data <- monthly_results %>%
  mutate(
    tree_nmol = Phi_tree_umol_m2_s * 1000,
    soil_nmol = Phi_soil_umol_m2_s * 1000,
    total_nmol = Phi_plot_umol_m2_s * 1000
  ) %>%
  dplyr::select(month, tree_nmol, soil_nmol, total_nmol) %>%
  pivot_longer(cols = c(tree_nmol, soil_nmol, total_nmol),
               names_to = "component",
               values_to = "flux_nmol") %>%
  mutate(
    component = factor(component,
                       levels = c("tree_nmol", "soil_nmol", "total_nmol"),
                       labels = c("Tree stems", "Soil", "Total plot"))
  )

fig2 <- ggplot(monthly_plot_data, aes(x = month, y = flux_nmol, color = component)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(
    values = c("Tree stems" = "#2E7D32", "Soil" = "#8B4513", "Total plot" = "black"),
    name = ""
  ) +
  labs(
    x = "Month",
    y = expression(paste("Plot-level CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))
  ) +
  theme_pub +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = "black", size = 0.3))

# ggsave("outputs/figures/Figure2_Seasonal_Patterns.pdf", fig2, width = 7, height = 4, dpi = 300)
# ggsave("outputs/figures/Figure2_Seasonal_Patterns.png", fig2, width = 7, height = 4, dpi = 300)

# =============================================================================
# FIGURE 3: FEATURE IMPORTANCE (TOP 10 EACH)
# Side-by-side for tree and soil models
# =============================================================================

cat("\nCreating Figure 3: Feature Importance...\n")

# Tree importance (top 10)
tree_imp <- data.frame(
  feature = names(importance(TreeRF)),
  importance = as.numeric(importance(TreeRF))
) %>%
  arrange(desc(importance)) %>%
  head(10) %>%
  mutate(
    model = "Tree stems",
    feature_clean = gsub("species_factor", "", feature),
    feature_clean = gsub("\\.soil_moisture_at_tree", " × Moisture", feature_clean)
  )

# Soil importance (all features since fewer)
soil_imp <- data.frame(
  feature = names(importance(SoilRF)),
  importance = as.numeric(importance(SoilRF))
) %>%
  arrange(desc(importance)) %>%
  mutate(
    model = "Soil",
    feature_clean = feature
  )

# Combine
imp_combined <- bind_rows(tree_imp, soil_imp) %>%
  group_by(model) %>%
  mutate(importance_scaled = importance / max(importance) * 100) %>%
  ungroup()

fig3 <- ggplot(imp_combined, aes(x = reorder(feature_clean, importance), 
                                 y = importance_scaled)) +
  geom_col(aes(fill = model)) +
  coord_flip() +
  facet_wrap(~ model, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Tree stems" = "#2E7D32", "Soil" = "#8B4513")) +
  labs(
    x = "",
    y = "Relative importance (%)"
  ) +
  theme_pub +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 9))

# ggsave("outputs/figures/Figure3_Feature_Importance.pdf", fig3, width = 8, height = 5, dpi = 300)
# ggsave("outputs/figures/Figure3_Feature_Importance.png", fig3, width = 8, height = 5, dpi = 300)

# =============================================================================
# OPTIONAL FIGURE 4: RESIDUAL DIAGNOSTICS
# QQ plots and residuals vs fitted
# =============================================================================

cat("\nCreating Figure 4 (Supplementary): Residual Diagnostics...\n")

# Prepare residuals
tree_resid <- data.frame(
  fitted = TreeRF$predictions,
  residual = tree_train_complete$y_asinh - TreeRF$predictions,
  model = "Tree stems"
)

soil_resid <- data.frame(
  fitted = SoilRF$predictions,
  residual = soil_train_complete$y_asinh - SoilRF$predictions,
  model = "Soil"
)

resid_data <- bind_rows(tree_resid, soil_resid)

# Residuals vs fitted
p1 <- ggplot(resid_data, aes(x = fitted, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = model), alpha = 0.4, size = 1) +
  geom_smooth(aes(color = model), method = "loess", se = TRUE, size = 0.8) +
  facet_wrap(~ model, scales = "free", ncol = 2) +
  scale_color_manual(values = c("Tree stems" = "#2E7D32", "Soil" = "#8B4513")) +
  labs(x = "Fitted values (asinh-transformed)",
       y = "Residuals",
       title = "(a) Residuals vs Fitted") +
  theme_pub +
  theme(legend.position = "none", plot.title = element_text(size = 10))

# QQ plots
p2 <- ggplot(resid_data, aes(sample = residual)) +
  stat_qq(aes(color = model), alpha = 0.4, size = 1) +
  stat_qq_line(aes(color = model), size = 0.8) +
  facet_wrap(~ model, scales = "free", ncol = 2) +
  scale_color_manual(values = c("Tree stems" = "#2E7D32", "Soil" = "#8B4513")) +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles",
       title = "(b) Normal Q-Q Plot") +
  theme_pub +
  theme(legend.position = "none", plot.title = element_text(size = 10))

fig4 <- p1 / p2
# ggsave("outputs/figures/FigureS1_Residual_Diagnostics.pdf", fig4, width = 7, height = 6, dpi = 300)
# ggsave("outputs/figures/FigureS1_Residual_Diagnostics.png", fig4, width = 7, height = 6, dpi = 300)

# =============================================================================
# SUMMARY TABLE FOR MANUSCRIPT
# =============================================================================

cat("\nCreating summary statistics table...\n")

summary_table <- data.frame(
  Component = c("Tree stems", "Soil", "Total plot"),
  `R²` = c(
    round(TreeRF$r.squared, 3),
    round(SoilRF$r.squared, 3),
    NA
  ),
  `RMSE (nmol/m²/s)` = c(
    round(sqrt(TreeRF$prediction.error) * 1000, 1),
    round(sqrt(SoilRF$prediction.error) * 1000, 1),
    NA
  ),
  `Training obs` = c(
    nrow(tree_train_complete),
    nrow(soil_train_complete),
    NA
  ),
  `Annual flux (mg CH4/m²/yr)` = if (exists("annual_summary")) c(
    round(annual_summary$annual_tree_mg_m2, 1),
    round(annual_summary$annual_soil_mg_m2, 1),
    round(annual_summary$annual_plot_mg_m2, 1)
  ) else c(NA, NA, NA),
  check.names = FALSE
)

write.csv(summary_table, "outputs/tables/Table1_Model_Summary.csv", row.names = FALSE)
print(summary_table)

cat("\n=== PUBLICATION FIGURES COMPLETE ===\n")
cat("\nGenerated files:\n")
cat("  • Figure1_Model_Performance.pdf/png - Main model validation\n")
cat("  • Figure2_Seasonal_Patterns.pdf/png - Monthly flux predictions\n")
cat("  • Figure3_Feature_Importance.pdf/png - Key predictors\n")
cat("  • FigureS1_Residual_Diagnostics.pdf/png - Supplementary diagnostics\n")
cat("  • Table1_Model_Summary.csv - Summary statistics\n")
cat("\nRECOMMENDATION FOR PAPER:\n")
cat("  Main text: Figures 1-2 (performance + seasonal patterns)\n")
cat("  Supplement: Figure 3 (feature importance) + Figure S1 (diagnostics)\n")

# =============================================================================
# COMPREHENSIVE FIGURE: PERFORMANCE + IMPORTANCE + PARTIAL DEPENDENCE
# All in one layout for both tree and soil models
# =============================================================================

library(ggplot2)
library(patchwork)
library(DescTools)
library(viridis)
library(hexbin)

cat("\n=== CREATING COMPREHENSIVE MODEL FIGURE ===\n")

# =============================================================================
# PANEL A & B: MODEL PERFORMANCE WITH DENSITY AND CCC
# =============================================================================

cat("\nPreparing performance panels...\n")

# Prepare data WITHOUT filtering negatives
perf_data <- bind_rows(
  tree_train_complete %>%
    mutate(
      observed_nmol = stem_flux_corrected * 1000,
      predicted_nmol = pred_flux * 1000,
      source = "Tree stems"
    ) %>%
    dplyr::select(observed_nmol, predicted_nmol, source),
  
  soil_train_complete %>%
    mutate(
      observed_nmol = soil_flux_umol_m2_s * 1000,
      predicted_nmol = pred_flux * 1000,
      source = "Soil"
    ) %>%
    dplyr::select(observed_nmol, predicted_nmol, source)
)

# Calculate statistics - FIX for CCC being a list
stats_by_source <- perf_data %>%
  group_by(source) %>%
  summarise(
    r2 = cor(observed_nmol, predicted_nmol, use = "complete.obs")^2,
    pearson_r = cor(observed_nmol, predicted_nmol, use = "complete.obs"),
    ccc_result = list(CCC(observed_nmol, predicted_nmol, na.rm = TRUE)),
    rmse = sqrt(mean((predicted_nmol - observed_nmol)^2, na.rm = TRUE)),
    bias = mean(predicted_nmol - observed_nmol, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(ccc = sapply(ccc_result, function(x) x$rho.c[1])) %>%
  dplyr::select(-ccc_result)

print(stats_by_source)

# Separate data for each source
perf_tree <- perf_data %>% filter(source == "Tree stems")
perf_soil <- perf_data %>% filter(source == "Soil")
stats_tree <- stats_by_source %>% filter(source == "Tree stems")
stats_soil <- stats_by_source %>% filter(source == "Soil")

# Panel A: Tree performance
pA <- ggplot(perf_tree, aes(x = observed_nmol, y = predicted_nmol)) +
  geom_hex(bins = 40) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.8, alpha = 0.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
           label = sprintf("CCC = %.3f\nR² = %.3f\nn = %d", 
                           stats_tree$ccc, stats_tree$r2, stats_tree$n),
           size = 3, fontface = "bold") +
  scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_fill_viridis_c(name = "Count", trans = "log10", option = "magma") +
  labs(x = expression(paste("Observed (nmol m"^-2, " s"^-1, ")")),
       y = expression(paste("Predicted (nmol m"^-2, " s"^-1, ")")),
       title = "(a) Tree stems") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        legend.key.height = unit(0.8, "cm"),
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

# Panel B: Soil performance
pB <- ggplot(perf_soil, aes(x = observed_nmol, y = predicted_nmol)) +
  geom_hex(bins = 40) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.8, alpha = 0.8) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
           label = sprintf("CCC = %.3f\nR² = %.3f\nn = %d", 
                           stats_soil$ccc, stats_soil$r2, stats_soil$n),
           size = 3, fontface = "bold") +
  scale_x_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_fill_viridis_c(name = "Count", trans = "log10", option = "magma") +
  labs(x = expression(paste("Observed (nmol m"^-2, " s"^-1, ")")),
       y = expression(paste("Predicted (nmol m"^-2, " s"^-1, ")")),
       title = "(b) Soil") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        legend.key.height = unit(0.8, "cm"),
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1)

# =============================================================================
# PANEL C & D: FEATURE IMPORTANCE
# =============================================================================

cat("\nPreparing feature importance panels...\n")

# Tree importance (top 8)
tree_imp <- data.frame(
  feature = names(importance(TreeRF)),
  importance = as.numeric(importance(TreeRF))
) %>%
  arrange(desc(importance)) %>%
  head(8) %>%
  mutate(
    feature_clean = gsub("species_factor", "", feature),
    feature_clean = gsub("\\.soil_moisture_at_tree", " × Moist", feature_clean),
    feature_clean = gsub("_", " ", feature_clean)
  )

pC <- ggplot(tree_imp, aes(x = reorder(feature_clean, importance), y = importance)) +
  geom_col(fill = "#2E7D32") +
  coord_flip() +
  labs(x = "", y = "Importance", title = "(c) Tree feature importance") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

# Soil importance (all features)
soil_imp <- data.frame(
  feature = names(importance(SoilRF)),
  importance = as.numeric(importance(SoilRF))
) %>%
  arrange(desc(importance)) %>%
  mutate(feature_clean = gsub("_", " ", feature))

pD <- ggplot(soil_imp, aes(x = reorder(feature_clean, importance), y = importance)) +
  geom_col(fill = "#8B4513") +
  coord_flip() +
  labs(x = "", y = "Importance", title = "(d) Soil feature importance") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 8))

# =============================================================================
# PANEL E & F: PARTIAL DEPENDENCE PLOTS
# =============================================================================

cat("\nComputing partial dependence plots...\n")

# Function to compute partial dependence
compute_partial_dependence <- function(model, data, var_name, n_points = 50) {
  var_range <- quantile(data[[var_name]], probs = c(0.05, 0.95), na.rm = TRUE)
  var_seq <- seq(var_range[1], var_range[2], length.out = n_points)
  
  predictions <- sapply(var_seq, function(val) {
    data_copy <- data
    data_copy[[var_name]] <- val
    mean(predict(model, data_copy)$predictions, na.rm = TRUE)
  })
  
  data.frame(
    x = var_seq,
    y = sinh(predictions) * 1000  # Convert back to nmol
  )
}

# Tree: moisture dependence
pd_tree_moisture <- compute_partial_dependence(
  TreeRF, 
  as.data.frame(X_tree), 
  "soil_moisture_at_tree",
  n_points = 50
)

pE <- ggplot(pd_tree_moisture, aes(x = x, y = y)) +
  geom_line(color = "#2E7D32", size = 1.2) +
  geom_ribbon(aes(ymin = y * 0.9, ymax = y * 1.1), alpha = 0.2, fill = "#2E7D32") +
  labs(x = "Soil moisture (m³/m³)",
       y = expression(paste("Predicted flux (nmol m"^-2, " s"^-1, ")")),
       title = "(e) Tree partial dependence: moisture") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

# Soil: temperature dependence
pd_soil_temp <- compute_partial_dependence(
  SoilRF,
  as.data.frame(X_soil),
  "soil_temp_C_mean",
  n_points = 50
)

pF <- ggplot(pd_soil_temp, aes(x = x, y = y)) +
  geom_line(color = "#8B4513", size = 1.2) +
  geom_ribbon(aes(ymin = y * 0.9, ymax = y * 1.1), alpha = 0.2, fill = "#8B4513") +
  labs(x = "Soil temperature (°C)",
       y = expression(paste("Predicted flux (nmol m"^-2, " s"^-1, ")")),
       title = "(f) Soil partial dependence: temperature") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

# =============================================================================
# COMBINE ALL PANELS
# =============================================================================

cat("\nCombining panels into final figure...\n")

# Layout: 
# Row 1: Performance (A, B)
# Row 2: Importance (C, D)  
# Row 3: Partial dependence (E, F)

layout <- "
AABBCC
AABBCC
DDDEEE
DDDEEE
FFFFGG
FFFFGG
"

final_fig <- pA + pB + pC + pD + pE + pF +
  plot_layout(design = layout, heights = c(1, 1, 1, 1, 1, 1))

# ggsave("outputs/figures/Figure_Comprehensive_Model_Analysis.pdf", final_fig, 
#        width = 12, height = 10, dpi = 300)
# ggsave("outputs/figures/Figure_Comprehensive_Model_Analysis.png", final_fig, 
#        width = 12, height = 10, dpi = 300)

cat("\n✓ Comprehensive figure saved\n")
cat("\nFigure includes:\n")
cat("  Panels A-B: Model performance (obs vs pred) with CCC\n")
cat("  Panels C-D: Top feature importances\n")
cat("  Panels E-F: Partial dependence for key predictors\n")
cat("\nStatistics summary:\n")
print(stats_by_source)













# =============================================================================
# PARTIAL DEPENDENCE FOR TOP 6 NON-INTERACTION FEATURES
# =============================================================================

library(ggplot2)
library(patchwork)

cat("\n=== IDENTIFYING TOP NON-INTERACTION FEATURES ===\n")

# Get importance and filter out interactions
tree_imp_all <- data.frame(
  feature = names(importance(TreeRF)),
  importance = as.numeric(importance(TreeRF))
) %>%
  arrange(desc(importance))

soil_imp_all <- data.frame(
  feature = names(importance(SoilRF)),
  importance = as.numeric(importance(SoilRF))
) %>%
  arrange(desc(importance))

# Filter out interaction features (contain "x" or "Moist" suffix for species×moisture)
tree_top6 <- tree_imp_all %>%
  filter(!grepl("\\.soil_moisture_at_tree$", feature)) %>%  # Remove species×moisture interactions
  head(6)

soil_top6 <- soil_imp_all %>%
  filter(!grepl("moisture_x|x_", feature)) %>%  # Remove interaction terms
  head(6)

cat("\nTop 6 non-interaction features for TREE model:\n")
print(tree_top6)

cat("\nTop 6 non-interaction features for SOIL model:\n")
print(soil_top6)

# =============================================================================
# COMPUTE PARTIAL DEPENDENCE FOR ALL TOP 6
# =============================================================================

compute_pd <- function(model, data, var_name, n_points = 50) {
  # Check if variable exists in data
  if(!var_name %in% names(data)) {
    cat("Warning: Variable", var_name, "not found in data\n")
    return(NULL)
  }
  
  var_range <- quantile(data[[var_name]], probs = c(0.05, 0.95), na.rm = TRUE)
  var_seq <- seq(var_range[1], var_range[2], length.out = n_points)
  
  pred_matrix <- sapply(var_seq, function(val) {
    data_copy <- data
    data_copy[[var_name]] <- val
    predict(model, data_copy)$predictions
  })
  
  mean_pred <- apply(pred_matrix, 2, mean, na.rm = TRUE)
  sd_pred <- apply(pred_matrix, 2, sd, na.rm = TRUE)
  
  data.frame(
    x = var_seq,
    y_mean = sinh(mean_pred) * 1000,
    y_lower = sinh(mean_pred - sd_pred) * 1000,
    y_upper = sinh(mean_pred + sd_pred) * 1000,
    feature = var_name
  )
}

# Tree partial dependence
tree_pd_list <- lapply(tree_top6$feature, function(feat) {
  compute_pd(TreeRF, as.data.frame(X_tree), feat)
})
tree_pd_combined <- bind_rows(tree_pd_list)

# Soil partial dependence
soil_pd_list <- lapply(soil_top6$feature, function(feat) {
  compute_pd(SoilRF, as.data.frame(X_soil), feat)
})
soil_pd_combined <- bind_rows(soil_pd_list)

# Create clean feature names for plotting
clean_feature_name <- function(name) {
  name %>%
    gsub("species_factor", "", .) %>%
    gsub("_", " ", .) %>%
    gsub("dbh within z", "DBH (within-species)", .) %>%
    gsub("soil moisture at tree", "Soil moisture", .) %>%
    gsub("soil moisture at site", "Soil moisture", .) %>%
    gsub("soil temp C mean", "Soil temperature", .) %>%
    gsub("air temp C mean", "Air temperature", .) %>%
    gsub("taxon prior asinh", "Taxon prior", .) %>%
    gsub("chamber rigid", "Chamber (rigid)", .) %>%
    gsub("month sin", "Month (sine)", .) %>%
    gsub("month cos", "Month (cosine)", .)
}

tree_pd_combined$feature_clean <- clean_feature_name(tree_pd_combined$feature)
soil_pd_combined$feature_clean <- clean_feature_name(soil_pd_combined$feature)

# =============================================================================
# CREATE PARTIAL DEPENDENCE PLOTS
# =============================================================================

# Tree PD plots (2x3 grid)
p_tree_pd <- ggplot(tree_pd_combined, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#2E7D32") +
  geom_line(color = "#2E7D32", size = 1) +
  facet_wrap(~ feature_clean, scales = "free", ncol = 3) +
  labs(x = "",
       y = expression(paste("Partial effect (nmol m"^-2, " s"^-1, ")")),
       title = "(e) Tree stem partial dependence") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    strip.text = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

# Soil PD plots (2x3 grid)
p_soil_pd <- ggplot(soil_pd_combined, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#8B4513") +
  geom_line(color = "#8B4513", size = 1) +
  facet_wrap(~ feature_clean, scales = "free", ncol = 3) +
  labs(x = "",
       y = expression(paste("Partial effect (nmol m"^-2, " s"^-1, ")")),
       title = "(f) Soil partial dependence") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    strip.text = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# FINAL COMPREHENSIVE FIGURE
# =============================================================================

cat("\nCreating final comprehensive figure...\n")

# Layout:
# Row 1: Tree performance | Tree importance | [Tree PD spans both rows]
# Row 2: Soil performance | Soil importance | 

final_fig <- (pA | pC | p_tree_pd) / 
  (pB | pD | p_soil_pd) +
  plot_layout(widths = c(1, 0.8, 1.8))

# ggsave("outputs/figures/Figure_Comprehensive_Top6PD.pdf", final_fig, 
#        width = 16, height = 8, dpi = 300)
# ggsave("outputs/figures/Figure_Comprehensive_Top6PD.png", final_fig, 
#        width = 16, height = 8, dpi = 300)

# Alternative: Vertical layout for narrow journals
final_fig_vertical <- (pA | pB) / 
  (pC | pD) / 
  p_tree_pd / 
  p_soil_pd +
  plot_layout(heights = c(1, 0.8, 1, 1))

# ggsave("outputs/figures/Figure_Comprehensive_Top6PD_vertical.pdf", final_fig_vertical, 
#        width = 10, height = 14, dpi = 300)
# ggsave("outputs/figures/Figure_Comprehensive_Top6PD_vertical.png", final_fig_vertical, 
#        width = 10, height = 14, dpi = 300)

cat("\n✓ Comprehensive figure with top 6 non-interaction features created\n")
cat("\nGenerated two layouts:\n")
cat("  1. Horizontal (16×8): Performance | Importance | 6 PD plots\n")
cat("  2. Vertical (10×14): All panels stacked\n")

# Print summary
cat("\n=== PARTIAL DEPENDENCE SUMMARY ===\n")
cat("\nTREE MODEL - Top 6 features:\n")
for(i in 1:nrow(tree_top6)) {
  cat(sprintf("  %d. %s (importance: %.1f)\n", 
              i, tree_top6$feature[i], tree_top6$importance[i]))
}

cat("\nSOIL MODEL - Top 6 features:\n")
for(i in 1:nrow(soil_top6)) {
  cat(sprintf("  %d. %s (importance: %.1f)\n", 
              i, soil_top6$feature[i], soil_top6$importance[i]))
}















# =============================================================================
# PARTIAL DEPENDENCE: 4 FEATURES IN 2x2 GRID FOR EACH MODEL
# =============================================================================

library(ggplot2)
library(patchwork)

cat("\n=== COMPUTING PARTIAL DEPENDENCE FOR 4 KEY FEATURES ===\n")

# Function to compute PD with uncertainty
compute_pd <- function(model, data, var_name, n_points = 50) {
  if(!var_name %in% names(data)) {
    cat("Warning:", var_name, "not found\n")
    return(NULL)
  }
  
  var_range <- quantile(data[[var_name]], probs = c(0.05, 0.95), na.rm = TRUE)
  var_seq <- seq(var_range[1], var_range[2], length.out = n_points)
  
  pred_matrix <- sapply(var_seq, function(val) {
    data_copy <- data
    data_copy[[var_name]] <- val
    predict(model, data_copy)$predictions
  })
  
  mean_pred <- apply(pred_matrix, 2, mean, na.rm = TRUE)
  sd_pred <- apply(pred_matrix, 2, sd, na.rm = TRUE)
  
  data.frame(
    x = var_seq,
    y_mean = sinh(mean_pred) * 1000,
    y_lower = sinh(mean_pred - sd_pred) * 1000,
    y_upper = sinh(mean_pred + sd_pred) * 1000
  )
}

# =============================================================================
# TREE PARTIAL DEPENDENCE (4 features)
# =============================================================================

pd_tree_air <- compute_pd(TreeRF, as.data.frame(X_tree), "air_temp_C_mean")
pd_tree_soil_moist <- compute_pd(TreeRF, as.data.frame(X_tree), "soil_moisture_at_tree")
pd_tree_soil_temp <- compute_pd(TreeRF, as.data.frame(X_tree), "soil_temp_C_mean")
pd_tree_dbh <- compute_pd(TreeRF, as.data.frame(X_tree), "dbh_within_z")

# Create plots
p_tree_air <- ggplot(pd_tree_air, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#2E7D32") +
  geom_line(color = "#2E7D32", size = 1) +
  labs(x = "Air temperature (°C)", y = "",
       title = "(e) Tree partial dependence") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

p_tree_moist <- ggplot(pd_tree_soil_moist, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#2E7D32") +
  geom_line(color = "#2E7D32", size = 1) +
  labs(x = "Soil moisture (m³/m³)", y = "") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

p_tree_stemp <- ggplot(pd_tree_soil_temp, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#2E7D32") +
  geom_line(color = "#2E7D32", size = 1) +
  labs(x = "Soil temperature (°C)", y = "") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

p_tree_dbh <- ggplot(pd_tree_dbh, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#2E7D32") +
  geom_line(color = "#2E7D32", size = 1) +
  labs(x = "DBH (within-species z-score)", y = "") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

# Combine tree PD into 2x2 — title on top-left plot survives nesting
tree_pd_panel <- (p_tree_air | p_tree_moist) / (p_tree_stemp | p_tree_dbh)

# =============================================================================
# SOIL PARTIAL DEPENDENCE (4 features, no DBH)
# =============================================================================

pd_soil_air <- compute_pd(SoilRF, as.data.frame(X_soil), "air_temp_C_mean")
pd_soil_moist <- compute_pd(SoilRF, as.data.frame(X_soil), "soil_moisture_at_site")
pd_soil_temp <- compute_pd(SoilRF, as.data.frame(X_soil), "soil_temp_C_mean")
pd_soil_SI <- compute_pd(SoilRF, as.data.frame(X_soil), "SI")

# Create plots
p_soil_air <- ggplot(pd_soil_air, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#8B4513") +
  geom_line(color = "#8B4513", size = 1) +
  labs(x = "Air temperature (°C)", y = "",
       title = "(f) Soil partial dependence") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

p_soil_moist <- ggplot(pd_soil_moist, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#8B4513") +
  geom_line(color = "#8B4513", size = 1) +
  labs(x = "Soil moisture (m³/m³)", y = "") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

p_soil_stemp <- ggplot(pd_soil_temp, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#8B4513") +
  geom_line(color = "#8B4513", size = 1) +
  labs(x = "Soil temperature (°C)", y = "") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

p_soil_SI <- ggplot(pd_soil_SI, aes(x = x, y = y_mean)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2, fill = "#8B4513") +
  geom_line(color = "#8B4513", size = 1) +
  labs(x = "Seasonal index", y = "") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

# Combine soil PD into 2x2 — title on top-left plot survives nesting
soil_pd_panel <- (p_soil_air | p_soil_moist) / (p_soil_stemp | p_soil_SI)

# =============================================================================
# FINAL FIGURE: PERFORMANCE | IMPORTANCE | PD (2x2)
# =============================================================================

final_fig <- (pA | pC | tree_pd_panel) / 
  (pB | pD | soil_pd_panel) +
  plot_layout(widths = c(1, 0.8, 1.2))

# ggsave("outputs/figures/Figure_Final_2x2PD.pdf", final_fig, width = 15, height = 8, dpi = 300)
# ggsave("outputs/figures/Figure_Final_2x2PD.png", final_fig, width = 15, height = 8, dpi = 300)

cat("\n✓ Final figure with 2x2 partial dependence panels created\n")
cat("\nStructure:\n")
cat("  Row 1: Tree performance | Tree importance | Tree PD (2×2 grid)\n")
cat("  Row 2: Soil performance | Soil importance | Soil PD (2×2 grid)\n")
cat("\nTree PD features: Air temp, Soil moisture, Soil temp, DBH\n")
cat("Soil PD features: Air temp, Soil moisture, Soil temp, Seasonal index\n")











# =============================================================================
# ADD MONTHLY FLUX PANELS (ROW 3)
# =============================================================================

cat("\n=== CREATING MONTHLY FLUX PANELS ===\n")

# Prepare monthly data for plotting
monthly_tree <- tree_results %>%
  mutate(flux_nmol = mean_flux * 1000) %>%
  dplyr::select(month, flux_nmol)

monthly_soil <- soil_results %>%
  mutate(flux_nmol = mean_flux * 1000) %>%
  dplyr::select(month, flux_nmol)

# Tree monthly flux
p_tree_monthly <- ggplot(monthly_tree, aes(x = month, y = flux_nmol)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_line(color = "#2E7D32", size = 1.2) +
  geom_point(color = "#2E7D32", size = 2.5) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  labs(x = "Month",
       y = expression(paste("CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
       title = "(g) Tree monthly predictions") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Soil monthly flux
p_soil_monthly <- ggplot(monthly_soil, aes(x = month, y = flux_nmol)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_line(color = "#8B4513", size = 1.2) +
  geom_point(color = "#8B4513", size = 2.5) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  labs(x = "Month",
       y = expression(paste("CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
       title = "(h) Soil monthly predictions") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# =============================================================================
# FINAL FIGURE: 3 ROWS
# Row 1: Performance | Importance | PD (2x2)
# Row 2: [same for soil]
# Row 3: Monthly flux (tree) | Monthly flux (soil) | [empty or combined]
# =============================================================================

# Option 1: Monthly plots span their respective columns
final_fig_3row <- (pA | pC | tree_pd_panel) / 
  (pB | pD | soil_pd_panel) /
  (p_tree_monthly | p_soil_monthly | plot_spacer()) +
  plot_layout(widths = c(1, 0.8, 1.2),
              heights = c(1, 1, 0.6))

# ggsave("outputs/figures/Figure_Final_3Rows.pdf", final_fig_3row, 
#        width = 15, height = 11, dpi = 300)
# ggsave("outputs/figures/Figure_Final_3Rows.png", final_fig_3row, 
#        width = 15, height = 11, dpi = 300)

# Option 2: Monthly plots together in one wider panel
p_monthly_combined <- (p_tree_monthly | p_soil_monthly)

final_fig_3row_alt <- (pA | pC | tree_pd_panel) / 
  (pB | pD | soil_pd_panel) /
  (p_monthly_combined) +
  plot_layout(widths = c(1, 0.8, 1.2),
              heights = c(1, 1, 0.6))

# ggsave("outputs/figures/Figure_Final_3Rows_alt.pdf", final_fig_3row_alt, 
#        width = 15, height = 11, dpi = 300)
ggsave("outputs/figures/supplementary/figS15_rf_predictions.png", final_fig_3row_alt,
       width = 12, height = 8.8, dpi = 300)

# Option 3: Single combined monthly plot (both on same panel)
monthly_combined_data <- bind_rows(
  monthly_tree %>% mutate(source = "Tree stems"),
  monthly_soil %>% mutate(source = "Soil")
)

p_monthly_single <- ggplot(monthly_combined_data, aes(x = month, y = flux_nmol, 
                                                      color = source)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  scale_color_manual(values = c("Tree stems" = "#2E7D32", "Soil" = "#8B4513"),
                     name = "") +
  labs(x = "Month",
       y = expression(paste("CH"[4], " flux (nmol m"^-2, " s"^-1, ")")),
       title = "(g) Monthly flux predictions") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3)
  )

final_fig_3row_single <- (pA | pC | tree_pd_panel) / 
  (pB | pD | soil_pd_panel) /
  (plot_spacer() | p_monthly_single | plot_spacer()) +
  plot_layout(widths = c(1, 0.8, 1.2),
              heights = c(1, 1, 0.6))

# ggsave("outputs/figures/Figure_Final_3Rows_single.pdf", final_fig_3row_single, 
#        width = 15, height = 11, dpi = 300)
# ggsave("outputs/figures/Figure_Final_3Rows_single.png", final_fig_3row_single, 
#        width = 15, height = 11, dpi = 300)

cat("\n✓ Three versions created:\n")
cat("  1. Figure_Final_3Rows.pdf - Monthly plots in left two columns\n")
cat("  2. Figure_Final_3Rows_alt.pdf - Monthly plots span full width\n")
cat("  3. Figure_Final_3Rows_single.pdf - Combined monthly plot (centered)\n")
cat("\nFinal structure:\n")
cat("  Row 1: Tree performance | Tree importance | Tree PD (2×2)\n")
cat("  Row 2: Soil performance | Soil importance | Soil PD (2×2)\n")
cat("  Row 3: Monthly flux predictions\n")