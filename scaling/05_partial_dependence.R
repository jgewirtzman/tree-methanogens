# =============================================================================
# PARTIAL DEPENDENCE PLOTS
# Showing how each predictor affects predictions while holding others constant
# =============================================================================

library(pdp)  # For partial dependence plots - install.packages("pdp") if needed
library(gridExtra)  # For arranging multiple plots

cat("\n=== GENERATING PARTIAL DEPENDENCE PLOTS ===\n")

# =============================================================================
# TREE MODEL PARTIAL DEPENDENCE
# =============================================================================

cat("\nCalculating partial dependence for Tree model...\n")

# Select key features to analyze (avoid too many for computational reasons)
tree_key_features <- c(
  "soil_moisture_at_tree",
  "air_temp_C_mean", 
  "soil_temp_C_mean",
  "dbh_m",
  "SI",
  "month_sin",
  "month_cos",
  "chamber_rigid"
)

# Store PD plots for trees
tree_pd_plots <- list()

for(feature in tree_key_features) {
  if(feature %in% colnames(X_tree)) {
    cat("  Computing PD for", feature, "...\n")
    
    # Calculate partial dependence
    pd_result <- partial(
      TreeRF,
      pred.var = feature,
      train = X_tree,
      n.trees = TreeRF$num.trees,
      grid.resolution = 50,
      plot = FALSE
    )
    
    # Transform back from asinh scale
    pd_result$yhat_original <- sinh(pd_result$yhat)
    pd_result$yhat_nmol <- pd_result$yhat_original * 1000  # Convert to nmol
    
    # Create plot
    p <- ggplot(pd_result, aes_string(x = feature, y = "yhat_nmol")) +
      geom_line(color = "darkgreen", size = 1.2) +
      geom_ribbon(aes(ymin = yhat_nmol * 0.95, ymax = yhat_nmol * 1.05), 
                  alpha = 0.2, fill = "green") +
      labs(title = paste("Partial Dependence:", feature),
           x = feature,
           y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
    
    tree_pd_plots[[feature]] <- p
  }
}

# Combine tree PD plots
n_plots <- length(tree_pd_plots)
n_cols <- 3
n_rows <- ceiling(n_plots / n_cols)

png("PD_TREE_MODEL.png", width = 12, height = 4 * n_rows)
do.call(grid.arrange, c(tree_pd_plots, ncol = n_cols))
dev.off()

cat("✓ Tree model PD plots saved to PD_TREE_MODEL.png\n")

# =============================================================================
# SOIL MODEL PARTIAL DEPENDENCE
# =============================================================================

cat("\nCalculating partial dependence for Soil model...\n")

soil_key_features <- c(
  "soil_moisture_at_site",
  "soil_temp_C_mean",
  "air_temp_C_mean",
  "SI",
  "month_sin",
  "month_cos"
)

soil_pd_plots <- list()

for(feature in soil_key_features) {
  if(feature %in% colnames(X_soil)) {
    cat("  Computing PD for", feature, "...\n")
    
    pd_result <- partial(
      SoilRF,
      pred.var = feature,
      train = X_soil,
      n.trees = SoilRF$num.trees,
      grid.resolution = 50,
      plot = FALSE
    )
    
    # Transform back from asinh scale
    pd_result$yhat_original <- sinh(pd_result$yhat)
    pd_result$yhat_nmol <- pd_result$yhat_original * 1000
    
    p <- ggplot(pd_result, aes_string(x = feature, y = "yhat_nmol")) +
      geom_line(color = "brown", size = 1.2) +
      geom_ribbon(aes(ymin = yhat_nmol * 0.95, ymax = yhat_nmol * 1.05), 
                  alpha = 0.2, fill = "tan") +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      labs(title = paste("Partial Dependence:", feature),
           x = feature,
           y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
    
    soil_pd_plots[[feature]] <- p
  }
}

# Combine soil PD plots
n_plots_soil <- length(soil_pd_plots)
n_rows_soil <- ceiling(n_plots_soil / n_cols)

png("PD_SOIL_MODEL.png", width = 12, height = 4 * n_rows_soil)
do.call(grid.arrange, c(soil_pd_plots, ncol = n_cols))
dev.off()

cat("✓ Soil model PD plots saved to PD_SOIL_MODEL.png\n")

# =============================================================================
# 2D PARTIAL DEPENDENCE - INTERACTION EFFECTS
# =============================================================================

cat("\nCalculating 2D partial dependence for key interactions...\n")

# Tree model: Moisture × Temperature interaction
if(all(c("soil_moisture_at_tree", "air_temp_C_mean") %in% colnames(X_tree))) {
  cat("  Computing Tree: moisture × temperature interaction...\n")
  
  pd_2d_tree <- partial(
    TreeRF,
    pred.var = c("soil_moisture_at_tree", "air_temp_C_mean"),
    train = X_tree,
    n.trees = TreeRF$num.trees,
    grid.resolution = 25,
    plot = FALSE
  )
  
  # Transform predictions
  pd_2d_tree$yhat_nmol <- sinh(pd_2d_tree$yhat) * 1000
  
  p_tree_2d <- ggplot(pd_2d_tree, aes(x = soil_moisture_at_tree, y = air_temp_C_mean, z = yhat_nmol)) +
    geom_tile(aes(fill = yhat_nmol)) +
    geom_contour(color = "white", alpha = 0.5) +
    scale_fill_viridis_c(option = "plasma") +
    labs(title = "Tree Model: Moisture × Temperature Interaction",
         x = "Soil Moisture (m³/m³)",
         y = "Air Temperature (°C)",
         fill = "CH₄ flux\n(nmol/m²/s)") +
    theme_minimal()
  
  ggsave("PD_2D_TREE_MOISTURE_TEMP.png", p_tree_2d, width = 10, height = 8)
}

# Soil model: Moisture × Soil Temperature interaction
if(all(c("soil_moisture_at_site", "soil_temp_C_mean") %in% colnames(X_soil))) {
  cat("  Computing Soil: moisture × temperature interaction...\n")
  
  pd_2d_soil <- partial(
    SoilRF,
    pred.var = c("soil_moisture_at_site", "soil_temp_C_mean"),
    train = X_soil,
    n.trees = SoilRF$num.trees,
    grid.resolution = 25,
    plot = FALSE
  )
  
  pd_2d_soil$yhat_nmol <- sinh(pd_2d_soil$yhat) * 1000
  
  p_soil_2d <- ggplot(pd_2d_soil, aes(x = soil_moisture_at_site, y = soil_temp_C_mean, z = yhat_nmol)) +
    geom_tile(aes(fill = yhat_nmol)) +
    geom_contour(color = "white", alpha = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = "Soil Model: Moisture × Soil Temperature Interaction",
         x = "Soil Moisture (m³/m³)",
         y = "Soil Temperature (°C)",
         fill = "CH₄ flux\n(nmol/m²/s)") +
    theme_minimal()
  
  ggsave("PD_2D_SOIL_MOISTURE_TEMP.png", p_soil_2d, width = 10, height = 8)
}

cat("✓ 2D partial dependence plots saved\n")

# =============================================================================
# SPECIALIZED PD: MONTH EFFECTS (RECONSTRUCTED FROM SIN/COS)
# =============================================================================

cat("\nReconstructing month effects from cyclic features...\n")

# Create synthetic data varying only month
month_grid <- data.frame(
  month = 1:12,
  month_sin = sin(2 * pi * (1:12) / 12),
  month_cos = cos(2 * pi * (1:12) / 12)
)

# For trees
if(all(c("month_sin", "month_cos") %in% colnames(X_tree))) {
  # Create base dataset with median values for all other features
  base_tree <- X_tree[1, , drop = FALSE]
  for(col in colnames(base_tree)) {
    if(is.numeric(base_tree[[col]])) {
      base_tree[[col]] <- median(X_tree[[col]], na.rm = TRUE)
    }
  }
  
  # Predict for each month
  tree_month_pd <- data.frame(month = 1:12)
  for(i in 1:12) {
    test_data <- base_tree
    test_data$month_sin <- month_grid$month_sin[i]
    test_data$month_cos <- month_grid$month_cos[i]
    pred <- predict(TreeRF, test_data)$predictions
    tree_month_pd$yhat_nmol[i] <- sinh(pred) * 1000
  }
  
  p_tree_month <- ggplot(tree_month_pd, aes(x = month, y = yhat_nmol)) +
    geom_line(color = "darkgreen", size = 1.5) +
    geom_point(color = "darkgreen", size = 3) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "Tree Model: Pure Month Effect (from cyclic features)",
         subtitle = "All other variables held at median values",
         x = "Month",
         y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))) +
    theme_minimal()
  
  ggsave("PD_TREE_MONTH_EFFECT.png", p_tree_month, width = 10, height = 6)
}

# For soil
if(all(c("month_sin", "month_cos") %in% colnames(X_soil))) {
  base_soil <- X_soil[1, , drop = FALSE]
  for(col in colnames(base_soil)) {
    if(is.numeric(base_soil[[col]])) {
      base_soil[[col]] <- median(X_soil[[col]], na.rm = TRUE)
    }
  }
  
  soil_month_pd <- data.frame(month = 1:12)
  for(i in 1:12) {
    test_data <- base_soil
    test_data$month_sin <- month_grid$month_sin[i]
    test_data$month_cos <- month_grid$month_cos[i]
    pred <- predict(SoilRF, test_data)$predictions
    soil_month_pd$yhat_nmol[i] <- sinh(pred) * 1000
  }
  
  p_soil_month <- ggplot(soil_month_pd, aes(x = month, y = yhat_nmol)) +
    geom_line(color = "brown", size = 1.5) +
    geom_point(color = "brown", size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "Soil Model: Pure Month Effect (from cyclic features)",
         subtitle = "All other variables held at median values",
         x = "Month",
         y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))) +
    theme_minimal()
  
  ggsave("PD_SOIL_MONTH_EFFECT.png", p_soil_month, width = 10, height = 6)
}

# =============================================================================
# ICE PLOTS (Individual Conditional Expectation) for selected features
# =============================================================================

cat("\nGenerating ICE plots for key features...\n")

# Tree model - moisture ICE plot
if("soil_moisture_at_tree" %in% colnames(X_tree)) {
  cat("  Computing ICE for tree moisture...\n")
  
  # Sample subset for ICE (too many lines make plot unreadable)
  n_ice <- min(100, nrow(X_tree))
  ice_sample <- sample(1:nrow(X_tree), n_ice)
  
  ice_tree_moisture <- partial(
    TreeRF,
    pred.var = "soil_moisture_at_tree",
    train = X_tree[ice_sample, ],
    n.trees = TreeRF$num.trees,
    ice = TRUE,
    grid.resolution = 50,
    plot = FALSE
  )
  
  ice_tree_moisture$yhat_nmol <- sinh(ice_tree_moisture$yhat) * 1000
  
  # Calculate PD line (average of ICE curves)
  pd_line <- ice_tree_moisture %>%
    group_by(soil_moisture_at_tree) %>%
    summarise(pd_yhat = mean(yhat_nmol), .groups = "drop")
  
  p_ice <- ggplot() +
    geom_line(data = ice_tree_moisture, 
              aes(x = soil_moisture_at_tree, y = yhat_nmol, group = yhat.id),
              alpha = 0.2, color = "gray50") +
    geom_line(data = pd_line,
              aes(x = soil_moisture_at_tree, y = pd_yhat),
              color = "darkgreen", size = 2) +
    labs(title = "ICE Plot: Tree Moisture Effect",
         subtitle = paste("Gray lines: individual predictions (n=", n_ice, "), Green: average (PD)"),
         x = "Soil Moisture (m³/m³)",
         y = expression(paste("Predicted CH"[4], " flux (nmol m"^-2, " s"^-1, ")"))) +
    theme_minimal()
  
  ggsave("ICE_TREE_MOISTURE.png", p_ice, width = 10, height = 7)
}

cat("\n✓ Partial dependence analysis complete\n")
cat("Generated plots:\n")
cat("  - PD_TREE_MODEL.png: 1D partial dependence for tree predictors\n")
cat("  - PD_SOIL_MODEL.png: 1D partial dependence for soil predictors\n")
cat("  - PD_2D_*.png: 2D interaction plots\n")
cat("  - PD_*_MONTH_EFFECT.png: Reconstructed month effects\n")
cat("  - ICE_TREE_MOISTURE.png: Individual conditional expectation plot\n")