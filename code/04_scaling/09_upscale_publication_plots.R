# ==============================================================================
# Upscaled Flux Publication Plots
# ==============================================================================
# Purpose: Publication-quality plots of upscaled tree and soil flux results
#   with species breakdown.
#
# Pipeline stage: 4 — Visualization
# Run after: 03_predict_tree_flux.R, 03_predict_soil_flux.R
#
# Inputs:
#   - tree_monthly_predictions.RData (from 03_predict_tree_flux)
#   - soil_monthly_predictions.RData (from 03_predict_soil_flux)
#
# Outputs:
#   - combined flux figure PNGs/PDFs (to outputs/figures/)
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
rm(monthly_predictions)

# Load tree predictions  
load("../../outputs/models/tree_monthly_predictions.RData")
tree_predictions <- monthly_predictions
rm(monthly_predictions)

cat("  Soil predictions loaded: ", nrow(soil_predictions), " points\n")
cat("  Tree predictions loaded: ", nrow(tree_predictions), " points\n\n")

# =============================================================================
# SELECT TARGET MONTHS
# =============================================================================

target_months <- c(12, 4, 7, 9)  # Dec, Apr, Jul, Sep

available_soil <- unique(soil_predictions$month)
available_tree <- unique(tree_predictions$month)
available_both <- intersect(available_soil, available_tree)

cat("Available months in both datasets:", paste(sort(available_both), collapse = ", "), "\n")

display_months <- intersect(target_months, available_both)
if (length(display_months) == 0) stop("No target months available in both datasets!")

cat("Will display months:", paste(display_months, collapse = ", "), "\n\n")

# =============================================================================
# ROTATED MINIMUM BOUNDING BOX (TREE) + CLIP SOIL TO THAT EXTENT
# =============================================================================

calculate_minimum_bounding_box <- function(points, buffer_pct = 0.01) {
  x <- points$x; y <- points$y
  angles <- seq(0, 179, by = 1) * pi / 180
  min_area <- Inf; best_angle <- 0; best_box <- NULL
  
  for (angle in angles) {
    cos_a <- cos(angle); sin_a <- sin(angle)
    x_rot <- x * cos_a - y * sin_a
    y_rot <- x * sin_a + y * cos_a
    xr <- range(x_rot); yr <- range(y_rot)
    area <- diff(xr) * diff(yr)
    if (area < min_area) {
      min_area <- area; best_angle <- angle
      xb <- diff(xr) * buffer_pct; yb <- diff(yr) * buffer_pct
      xr_b <- xr + c(-xb, xb); yr_b <- yr + c(-yb, yb)
      corners_rot <- expand.grid(x = xr_b, y = yr_b)
      corners_orig <- data.frame(
        x = corners_rot$x * cos(-angle) - corners_rot$y * sin(-angle),
        y = corners_rot$x * sin(-angle) + corners_rot$y * cos(-angle)
      )
      best_box <- list(
        corners = corners_orig,
        angle = best_angle * 180 / pi,
        area = min_area,
        rotated_ranges = list(x = xr_b, y = yr_b)
      )
    }
  }
  best_box
}

point_in_rotated_box <- function(test_points, box_info) {
  angle <- box_info$angle * pi / 180
  cos_a <- cos(angle); sin_a <- sin(angle)
  x_rot <- test_points$x * cos_a - test_points$y * sin_a
  y_rot <- test_points$x * sin_a + test_points$y * cos_a
  x_in <- x_rot >= box_info$rotated_ranges$x[1] & x_rot <= box_info$rotated_ranges$x[2]
  y_in <- y_rot >= box_info$rotated_ranges$y[1] & y_rot <= box_info$rotated_ranges$y[2]
  x_in & y_in
}

cat("Calculating minimum rotated bounding box for tree locations...\n")
tree_coords <- data.frame(x = tree_predictions$x, y = tree_predictions$y)
tree_bbox <- calculate_minimum_bounding_box(tree_coords, buffer_pct = 0.02)
cat("  Optimal rotation angle:", round(tree_bbox$angle, 2), "degrees\n")
cat("  Bounding box area:", round(tree_bbox$area, 2), "square units\n")

cat("Clipping soil predictions to tree bounding box...\n")
soil_coords <- data.frame(x = soil_predictions$x, y = soil_predictions$y)
inside_bbox <- point_in_rotated_box(soil_coords, tree_bbox)
soil_predictions_clipped <- soil_predictions[inside_bbox, ]
cat("  Original soil points:", nrow(soil_predictions), "\n")
cat("  Clipped soil points:", nrow(soil_predictions_clipped), "\n\n")

# Shared spatial extent after clipping
x_range <- range(c(soil_predictions_clipped$x, tree_predictions$x), na.rm = TRUE)
y_range <- range(c(soil_predictions_clipped$y, tree_predictions$y), na.rm = TRUE)
cat("Shared spatial extent:\n")
cat("  x range:", paste(round(x_range, 3), collapse = " to "), "\n")
cat("  y range:", paste(round(y_range, 3), collapse = " to "), "\n\n")

# =============================================================================
# PLOTTING FUNCTIONS (titles + subtitles inside panels; default gridlines)
#   - subtitles carry summary stats (Tree: mean; Soil: mean + % uptake)
#   - legends shared per row via plot_layout(guides='collect')
#   - coord_equal with shared x/y limits for consistent extents
# =============================================================================

create_soil_month_map <- function(month_val, pred_df, shared_limits = NULL,
                                  xlim = NULL, ylim = NULL, show_legend = TRUE) {
  month_data <- pred_df %>% dplyr::filter(month == month_val)
  flux_limits <- if (is.null(shared_limits)) {
    quantile(pred_df$flux_nmol_m2_s, c(0.02, 0.98), na.rm = TRUE)
  } else shared_limits
  
  mean_flux <- mean(month_data$flux_nmol_m2_s, na.rm = TRUE)
  pct_uptake <- 100 * mean(month_data$uptake, na.rm = TRUE)
  
  ggplot(month_data, aes(x = x, y = y, fill = flux_nmol_m2_s)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = flux_limits, oob = scales::squish,
      name = "Soil Flux\n(nmol/m²/s)"
    ) +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    coord_equal(xlim = xlim, ylim = ylim, clip = "off") +
    labs(
      # no title -> removes month name from soil plots
      subtitle = sprintf("Mean: %.2f\nUptake: %.1f%%", mean_flux, pct_uptake)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),                # ensure no leftover space
      plot.subtitle = element_text(size = 9, hjust = 0.5, margin = margin(b = 2)),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}


create_tree_month_map <- function(month_val, pred_df, shared_limits = NULL,
                                  xlim = NULL, ylim = NULL, show_legend = TRUE) {
  month_data <- pred_df %>% dplyr::filter(month == month_val)
  flux_limits <- if (is.null(shared_limits)) {
    quantile(pred_df$flux_nmol_m2_s, c(0.02, 0.98), na.rm = TRUE)
  } else shared_limits
  
  mean_flux <- mean(month_data$flux_nmol_m2_s, na.rm = TRUE)
  
  ggplot(month_data, aes(x = x, y = y)) +
    geom_point(aes(size = BasalArea_m2, color = flux_nmol_m2_s),
               alpha = 0.7, stroke = 0.1) +
    scale_color_viridis_c(
      option = "plasma",
      limits = flux_limits, oob = scales::squish,
      name = "Tree Flux\n(nmol/m²/s)"
    ) +
    scale_size_continuous(name = "Basal Area", range = c(0.3, 2.5), guide = "none") +
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    coord_equal(xlim = xlim, ylim = ylim, clip = "off") +
    labs(
      title = month.name[month_val],
      subtitle = sprintf("Mean: %.2f", mean_flux)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, margin = margin(b = 2)),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}

# =============================================================================
# INDIVIDUAL PANELS (shared color limits; only first in each row shows legend)
# =============================================================================

cat("Creating individual month plots...\n")

soil_limits <- quantile(soil_predictions_clipped$flux_nmol_m2_s, c(0.02, 0.98), na.rm = TRUE)
tree_limits <- quantile(tree_predictions$flux_nmol_m2_s, c(0.02, 0.98), na.rm = TRUE)

cat("  Soil flux range for color scale:", paste(round(soil_limits, 2), collapse = " to "), "\n")
cat("  Tree flux range for color scale:", paste(round(tree_limits, 2), collapse = " to "), "\n")

tree_plots <- vector("list", length(display_months))
soil_plots <- vector("list", length(display_months))

for (i in seq_along(display_months)) {
  m <- display_months[i]
  cat("  Creating plots for month", m, "(", month.name[m], ")...\n")
  tree_plots[[i]] <- create_tree_month_map(m, tree_predictions, tree_limits,
                                           xlim = x_range, ylim = y_range,
                                           show_legend = (i == 1))
  soil_plots[[i]] <- create_soil_month_map(m, soil_predictions_clipped, soil_limits,
                                           xlim = x_range, ylim = y_range,
                                           show_legend = (i == 1))
}

# =============================================================================
# COMBINE PLOTS — CLEAN ROWS + SHARED LEGENDS PER ROW
# =============================================================================

cat("\nCombining plots into final figure (shared legends per row)...\n")

assemble_row <- function(plot_list) {
  n <- length(plot_list)
  if (n == 4) {
    row_core <- plot_list[[1]] | plot_list[[2]] | plot_list[[3]] | plot_list[[4]]
  } else if (n == 3) {
    empty_plot <- ggplot() + theme_void()
    row_core <- plot_list[[1]] | plot_list[[2]] | plot_list[[3]] | empty_plot
  } else if (n == 2) {
    row_core <- plot_list[[1]] | plot_list[[2]] | plot_list[[1]] | plot_list[[2]]
  } else {
    row_core <- plot_list[[1]]
  }
  # Collect only the (single) legend from the first panel
  row_core + plot_layout(guides = "collect")
}

tree_row <- assemble_row(tree_plots)
soil_row <- assemble_row(soil_plots)

combined_figure <- (tree_row / soil_row) +
  plot_annotation(
    title = "Seasonal CH₄ Flux Patterns: Trees and Soil",
    subtitle = paste("Months displayed:", paste(month.name[display_months], collapse = ", ")),
    caption = "Top: Tree flux (color = flux, size = basal area) • Bottom: Soil flux (interpolated raster)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 0.5)
    )
  )

# =============================================================================
# SAVE MAIN FIGURE
# =============================================================================

cat("\nSaving combined figure...\n")
# ggsave("../../outputs/figures/combined_tree_soil_flux_seasonal.png",
#        combined_figure, width = 16, height = 8, dpi = 300, bg = "white")
cat("  Saved: combined_tree_soil_flux_seasonal.png\n")

# ggsave("../../outputs/figures/combined_tree_soil_flux_seasonal.pdf",
#        combined_figure, width = 16, height = 8, bg = "white")
cat("  Saved: combined_tree_soil_flux_seasonal.pdf\n")

# =============================================================================
# SUMMARY STATISTICS FOR DISPLAYED MONTHS
# =============================================================================

cat("\n=== SUMMARY FOR DISPLAYED MONTHS ===\n")
for (m in display_months) {
  cat("\n", month.name[m], ":\n", sep = "")
  soil_month <- soil_predictions_clipped %>% dplyr::filter(month == m)
  cat("  Soil - Mean flux:", round(mean(soil_month$flux_nmol_m2_s, na.rm = TRUE), 3),
      "nmol/m²/s | Uptake:", round(100 * mean(soil_month$uptake, na.rm = TRUE), 1), "%\n")
  tree_month <- tree_predictions %>% dplyr::filter(month == m)
  cat("  Tree - Mean flux:", round(mean(tree_month$flux_nmol_m2_s, na.rm = TRUE), 3),
      "nmol/m²/s | N trees:", length(unique(tree_month$tree_id)), "\n")
}
cat("\n=== COMBINED FIGURE COMPLETE ===\n")
cat("Finished at:", format(Sys.time()), "\n")

# =============================================================================
# FACETED ANNUAL FLUX BY SPECIES (WITH FULL NAMES) + SPATIAL, BOXPLOTS, ETC.
#   Assumes 'species_counts' and 'annual_avg' already exist in env
#   Uses shared x/y extents for spatial faceting
# =============================================================================

cat("\nCreating faceted annual flux distribution by species...\n")

species_mapping <- c(
  "ACRU" = "Acer rubrum", "ACSA" = "Acer saccharum", "ACPE" = "Acer pensylvanicum",
  "BEAL" = "Betula alleghaniensis", "BELE" = "Betula lenta", "BEPA" = "Betula papyrifera",
  "BEPO" = "Betula populifolia", "FAGR" = "Fagus grandifolia", "FRAM" = "Fraxinus americana",
  "PIST" = "Pinus strobus", "QURU" = "Quercus rubra", "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina", "TSCA" = "Tsuga canadensis", "CAOV" = "Carya ovata",
  "CALA" = "Carya laciniosa", "KALA" = "Kalmia latifolia", "PRSE" = "Prunus serotina",
  "SAAL" = "Sassafras albidum"
)

species_for_dist <- species_counts %>%
  dplyr::filter(n > 50, species != "KALA") %>%
  dplyr::pull(species)
cat("  Species with n>50 (excluding Kalmia latifolia):", length(species_for_dist), "species\n")

species_dist_data <- annual_avg %>%
  dplyr::filter(species %in% species_for_dist) %>%
  dplyr::mutate(
    species_full = ifelse(species %in% names(species_mapping), species_mapping[species], species)
  ) %>%
  dplyr::group_by(species, species_full) %>%
  dplyr::mutate(n_trees = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(species_label = paste0(species_full, "\n(n=", n_trees, ")"))

# Histogram + density
p_species_dist <- ggplot(species_dist_data, aes(x = mean_flux_nmol)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_density(color = "darkred", linewidth = 0.8, alpha = 0.5) +
  facet_wrap(~ species_label, scales = "free_y", ncol = 4) +
  labs(
    title = "Distribution of Mean Annual CH₄ Flux by Species",
    subtitle = paste(length(species_for_dist), "species with n>50 (Kalmia latifolia excluded)"),
    x = "Mean Annual Flux (nmol/m²/s)", y = "Density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8, face = "italic"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6)
  )

# ggsave("../../outputs/figures/tree_flux_distribution_by_species.png", p_species_dist, width = 14, height = 10, dpi = 150)
cat("  Saved: tree_flux_distribution_by_species.png\n")

# Boxplot + jitter
species_order <- species_dist_data %>%
  dplyr::group_by(species_full) %>%
  dplyr::summarise(mean_flux_sp = mean(mean_flux_nmol), .groups = "drop") %>%
  dplyr::arrange(mean_flux_sp)

species_dist_data$species_full_ordered <- factor(species_dist_data$species_full, levels = species_order$species_full)

p_species_box <- ggplot(species_dist_data, aes(x = mean_flux_nmol, y = species_full_ordered)) +
  geom_boxplot(aes(fill = species_full_ordered), alpha = 0.7, show.legend = FALSE) +
  geom_point(aes(color = species_full_ordered), position = position_jitter(height = 0.2), size = 1, alpha = 0.5, show.legend = FALSE) +
  scale_fill_viridis_d(option = "plasma") +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Mean Annual CH₄ Flux Distribution by Species",
    subtitle = paste(length(species_for_dist), "species with n>50 (Kalmia latifolia excluded) • Points are trees"),
    x = "Mean Annual Flux (nmol/m²/s)", y = "Species"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, face = "italic"))

# ggsave("../../outputs/figures/tree_flux_boxplot_by_species.png", p_species_box, width = 10, height = 8, dpi = 150)
cat("  Saved: tree_flux_boxplot_by_species.png\n")

# Spatial faceting with shared XY extent
p_species_spatial <- ggplot(species_dist_data, aes(x = x, y = y)) +
  geom_point(aes(color = mean_flux_nmol, size = BasalArea_m2), alpha = 0.7) +
  facet_wrap(~ species_full, ncol = 4) +
  scale_color_viridis_c(option = "plasma", name = "Annual Mean Tree\nFlux (nmol/m²/s)") +
  scale_size_continuous(range = c(0.5, 2), guide = "none") +
  scale_x_continuous(limits = x_range, expand = c(0, 0)) +
  scale_y_continuous(limits = y_range, expand = c(0, 0)) +
  coord_equal(xlim = x_range, ylim = y_range) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8, face = "italic"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8)
  )

# ggsave("../../outputs/figures/tree_flux_spatial_faceted_species.png", p_species_spatial, width = 14, height = 10, dpi = 150)
cat("  Saved: tree_flux_spatial_faceted_species.png\n")

# Species stats table
faceted_species_stats <- species_dist_data %>%
  dplyr::group_by(species, species_full) %>%
  dplyr::summarise(
    n_trees = dplyr::n(),
    mean_flux = mean(mean_flux_nmol),
    median_flux = median(mean_flux_nmol),
    sd_flux = sd(mean_flux_nmol),
    cv_flux = sd_flux / mean_flux * 100,
    min_flux = min(mean_flux_nmol),
    max_flux = max(mean_flux_nmol),
    q25_flux = quantile(mean_flux_nmol, 0.25),
    q75_flux = quantile(mean_flux_nmol, 0.75),
    iqr_flux = q75_flux - q25_flux,
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(mean_flux))

cat("\nStatistics for species with n>50 (excluding Kalmia latifolia):\n")
print(faceted_species_stats %>% dplyr::select(species_full, n_trees, mean_flux, median_flux, cv_flux))
write.csv(faceted_species_stats, "../../outputs/tables/tree_flux_faceted_species_stats.csv", row.names = FALSE)

# Bar + error
faceted_species_stats$species_full_ordered <- factor(faceted_species_stats$species_full, levels = species_order$species_full)

p_species_bar <- ggplot(faceted_species_stats, aes(x = species_full_ordered, y = mean_flux)) +
  geom_bar(stat = "identity", aes(fill = mean_flux), show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_flux - sd_flux, ymax = mean_flux + sd_flux), width = 0.3, alpha = 0.5) +
  scale_fill_viridis_c(option = "plasma") +
  coord_flip() +
  labs(
    title = "Mean Annual CH₄ Flux by Species (±1 SD)",
    subtitle = paste(length(species_for_dist), "species with n>50 (Kalmia latifolia excluded)"),
    x = "Species", y = "Mean Annual Flux (nmol/m²/s)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, face = "italic"))

# ggsave("../../outputs/figures/tree_flux_barplot_by_species.png", p_species_bar, width = 10, height = 8, dpi = 150)
cat("  Saved: tree_flux_barplot_by_species.png\n")

# Violin as option
p_species_violin <- ggplot(species_dist_data, aes(x = species_full_ordered, y = mean_flux_nmol)) +
  geom_violin(aes(fill = species_full_ordered), alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(width = 0.1, alpha = 0.8, outlier.size = 0.5) +
  scale_fill_viridis_d(option = "plasma") +
  coord_flip() +
  labs(
    title = "Distribution of Mean Annual CH₄ Flux by Species",
    subtitle = paste(length(species_for_dist), "species with n>50 (Kalmia latifolia excluded)"),
    x = "Species", y = "Mean Annual Flux (nmol/m²/s)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9, face = "italic"))

# ggsave("../../outputs/figures/tree_flux_violin_by_species.png", p_species_violin, width = 10, height = 8, dpi = 150)
cat("  Saved: tree_flux_violin_by_species.png\n")

# =============================================================================
# MEGA COMBINED PLOT (top: trees, mid: soil, bottom: spatial species)
# =============================================================================

mega_combined_plot <- (tree_row / soil_row / p_species_spatial) +
  plot_layout(heights = c(1, 1, 1.5)) +
  plot_annotation(
    #title = "Seasonal & Species-Level CH₄ Flux Overview",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
mega_combined_plot

ggsave("../../outputs/figures/main/fig9_upscaled_flux_seasonal.png",
       mega_combined_plot, width = 12, height = 12, dpi = 300, bg = "white")

cat("\n✓ All figures saved successfully!\n")
