# ==============================================================================
# Seasonal Flux Maps (Figure 9)
# ==============================================================================
# Purpose: Creates combined tree+soil seasonal flux maps.
#
# Pipeline stage: 4 — Visualization
# Run after: 03_predict_tree_flux.R, 03_predict_soil_flux.R
#
# Inputs:
#   - tree_monthly_predictions.RData (from 03_predict_tree_flux)
#   - soil_monthly_predictions.RData (from 03_predict_soil_flux)
#
# Outputs:
#   - combined_tree_soil_flux_seasonal.png/pdf (to outputs/figures/)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)

cat("=== CREATING COMBINED FLUX FIGURE ===\n")
cat("Starting at:", format(Sys.time()), "\n\n")

# =============================================================================
# LOAD SAVED PREDICTION DATA
# =============================================================================

cat("Loading saved prediction data...\n")

# Load soil predictions
load("../../outputs/models/soil_monthly_predictions.RData")
soil_predictions <- monthly_predictions
rm(monthly_predictions)  # Clear to avoid confusion

# Load tree predictions  
load("../../outputs/models/tree_monthly_predictions.RData")
tree_predictions <- monthly_predictions
rm(monthly_predictions)

cat("  Soil predictions loaded: ", nrow(soil_predictions), " points\n")
cat("  Tree predictions loaded: ", nrow(tree_predictions), " points\n\n")

# =============================================================================
# SELECT TARGET MONTHS
# =============================================================================

# Target months: Dec (12), April (4), July (7), Sept (9)
target_months <- c(12, 4, 7, 9)
month_labels <- c("December", "April", "July", "September")

# Check available months
available_soil <- unique(soil_predictions$month)
available_tree <- unique(tree_predictions$month)
available_both <- intersect(available_soil, available_tree)

cat("Available months in both datasets:", paste(sort(available_both), collapse = ", "), "\n")

# Use intersection of target and available months
display_months <- intersect(target_months, available_both)

if (length(display_months) == 0) {
  stop("No target months available in both datasets!")
}

cat("Will display months:", paste(display_months, collapse = ", "), "\n\n")

# =============================================================================
# CREATE PLOTTING FUNCTIONS
# =============================================================================

# Function to create soil month map
create_soil_month_map <- function(month_val, pred_df, shared_limits = NULL, show_legend = TRUE) {
  month_data <- pred_df %>% filter(month == month_val)
  
  # Use shared limits if provided, otherwise calculate
  if (is.null(shared_limits)) {
    flux_limits <- quantile(pred_df$flux_nmol_m2_s, c(0.02, 0.98), na.rm = TRUE)
  } else {
    flux_limits <- shared_limits
  }
  
  ggplot(month_data, aes(x = x, y = y, fill = flux_nmol_m2_s)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white", 
      high = "red",
      midpoint = 0,
      limits = flux_limits,
      oob = scales::squish,
      name = "Soil Flux\n(nmol/m²/s)"
    ) +
    coord_equal() +
    labs(
      title = NULL,
      subtitle = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(5, 2, 5, 2, "pt"),  # Reduce left/right margins
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = if(show_legend) "right" else "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}

# Function to create tree month map  
create_tree_month_map <- function(month_val, pred_df, shared_limits = NULL, show_legend = TRUE) {
  month_data <- pred_df %>% filter(month == month_val)
  
  # Use shared limits if provided
  if (is.null(shared_limits)) {
    flux_limits <- quantile(pred_df$flux_nmol_m2_s, c(0.02, 0.98), na.rm = TRUE)
  } else {
    flux_limits <- shared_limits
  }
  
  # Month name for title
  month_name <- month.name[month_val]
  
  ggplot(month_data, aes(x = x, y = y)) +
    geom_point(aes(size = BasalArea_m2, color = flux_nmol_m2_s), 
               alpha = 0.7, stroke = 0.1) +
    scale_color_viridis_c(
      option = "plasma",
      limits = flux_limits,
      oob = scales::squish,
      name = "Tree Flux\n(nmol/m²/s)"
    ) +
    scale_size_continuous(
      name = "Basal Area",
      range = c(0.3, 2.5),
      guide = "none"
    ) +
    coord_equal() +
    labs(
      title = month_name,
      subtitle = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_blank(),
      plot.margin = margin(5, 2, 5, 2, "pt"),  # Reduce left/right margins
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = if(show_legend) "right" else "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}

# =============================================================================
# CREATE INDIVIDUAL PLOTS
# =============================================================================

cat("Creating individual month plots...\n")

# Calculate shared color limits for better comparison
soil_limits <- quantile(soil_predictions$flux_nmol_m2_s, c(0.02, 0.98), na.rm = TRUE)
tree_limits <- quantile(tree_predictions$flux_nmol_m2_s, c(0.02, 0.98), na.rm = TRUE)

cat("  Soil flux range for color scale:", round(soil_limits, 2), "\n")
cat("  Tree flux range for color scale:", round(tree_limits, 2), "\n")

# Create plots for each month
tree_plots <- list()
soil_plots <- list()

for (i in 1:length(display_months)) {
  m <- display_months[i]
  cat("  Creating plots for month", m, "(", month.name[m], ")...\n")
  
  # Only show legend on the last plot of each row
  show_legend <- (i == length(display_months))
  
  tree_plots[[i]] <- create_tree_month_map(m, tree_predictions, tree_limits, show_legend = show_legend)
  soil_plots[[i]] <- create_soil_month_map(m, soil_predictions, soil_limits, show_legend = show_legend)
}

# =============================================================================
# COMBINE PLOTS - HORIZONTAL LAYOUT
# =============================================================================

cat("\nCombining plots into final figure...\n")

# Create the combined layout based on number of available months
# Back to: rows = data types (tree top, soil bottom), columns = months
if (length(display_months) == 4) {
  # Ideal case: all 4 months available
  tree_row <- tree_plots[[1]] | tree_plots[[2]] | tree_plots[[3]] | tree_plots[[4]]
  soil_row <- soil_plots[[1]] | soil_plots[[2]] | soil_plots[[3]] | soil_plots[[4]]
  
} else if (length(display_months) == 3) {
  # 3 months: add empty plot
  empty_plot <- ggplot() + theme_void()
  tree_row <- tree_plots[[1]] | tree_plots[[2]] | tree_plots[[3]] | empty_plot
  soil_row <- soil_plots[[1]] | soil_plots[[2]] | soil_plots[[3]] | empty_plot
  
} else if (length(display_months) == 2) {
  # 2 months: duplicate to fill space
  tree_row <- tree_plots[[1]] | tree_plots[[2]] | tree_plots[[1]] | tree_plots[[2]]
  soil_row <- soil_plots[[1]] | soil_plots[[2]] | soil_plots[[1]] | soil_plots[[2]]
  
} else {
  # 1 month: show just once in each row
  tree_row <- tree_plots[[1]]
  soil_row <- soil_plots[[1]]
}

# Combine top (tree) and bottom (soil) rows
combined_figure <- tree_row / soil_row

# No annotations - clean figure

# =============================================================================
# SAVE FIGURE
# =============================================================================

cat("\nSaving combined figure...\n")

# Save as high-quality PNG
# Wide format for horizontal layout
ggsave("../../outputs/figures/combined_tree_soil_flux_seasonal.png", 
       combined_figure, 
       width = 16, 
       height = 8, 
       dpi = 300,
       bg = "white")

cat("  Saved: combined_tree_soil_flux_seasonal.png\n")

# Also save as PDF for publication
ggsave("../../outputs/figures/combined_tree_soil_flux_seasonal.pdf", 
       combined_figure, 
       width = 16, 
       height = 8,
       bg = "white")

cat("  Saved: combined_tree_soil_flux_seasonal.pdf\n")

# =============================================================================
# SUMMARY STATISTICS FOR DISPLAYED MONTHS
# =============================================================================

cat("\n=== SUMMARY FOR DISPLAYED MONTHS ===\n")

for (m in display_months) {
  cat("\n", month.name[m], ":\n", sep = "")
  
  # Soil stats
  soil_month <- soil_predictions %>% filter(month == m)
  cat("  Soil - Mean flux:", round(mean(soil_month$flux_nmol_m2_s), 3), 
      "nmol/m²/s | Uptake:", round(100 * mean(soil_month$uptake), 1), "%\n")
  
  # Tree stats  
  tree_month <- tree_predictions %>% filter(month == m)
  cat("  Tree - Mean flux:", round(mean(tree_month$flux_nmol_m2_s), 3),
      "nmol/m²/s | N trees:", length(unique(tree_month$tree_id)), "\n")
}

cat("\n=== COMBINED FIGURE COMPLETE ===\n")
cat("Figure shows", length(display_months), "months in horizontal layout\n")
cat("Tree flux (top row) and soil flux (bottom row)\n")
cat("Finished at:", format(Sys.time()), "\n")