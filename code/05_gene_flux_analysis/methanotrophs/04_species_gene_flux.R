# ==============================================================================
# Species-Level Gene-Flux Correlations (Figures S7, S8, S11)
# ==============================================================================
# Purpose: Species-level gene-flux correlation analysis and mixed models.
#
# Pipeline stage: 3 — Analysis
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#
# UPSTREAM DEPENDENCIES (must be run first, or objects must exist in env):
#   - 05_radial_gene_plots.R  -> p_mcra_result, p_sum_result, p_overlay
#   - 01_summary_stats.R      -> analysis_mcra, analysis_methanotroph,
#                                 analysis_ratio, cor_area_mcra,
#                                 cor_area_methanotroph, pearson_ratio,
#                                 species_comparison_data, y_limit
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(ggrepel)
library(patchwork)

# ============================================================
# FUNCTION TO SPLIT SPECIES NAMES INTO TWO LINES
# ============================================================

break_species_name <- function(name) {
  parts <- strsplit(name, " ")[[1]]
  if (length(parts) == 2) {
    return(paste0(parts[1], "\n", parts[2]))
  }
  return(name)
}

# ============================================================
# RE-CREATE RADIAL PLOTS WITH TWO-LINE SPECIES NAMES
# Requires: p_mcra_result, p_sum_result, p_overlay from
#   05_radial_gene_plots.R
# ============================================================

if (!exists("p_mcra_result") || !exists("p_sum_result") || !exists("p_overlay")) {
  cat("WARNING: Radial plot objects (p_mcra_result, p_sum_result, p_overlay) not found.\n")
  cat("  Run 05_radial_gene_plots.R first, or source it to populate these objects.\n")
  cat("  Skipping radial plot section.\n\n")
} else {

# mcrA plot with two-line species names
p_mcra_data <- p_mcra_result$plot$data
p_mcra_data$species_label <- sapply(as.character(p_mcra_data$species_label), break_species_name)

p_mcra <- ggplot(p_mcra_data, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = bquote("log"[10] ~ "(mcrA)"),
    option = "inferno",
    direction = 1,
    breaks = function(x) c(min(x), mean(c(min(x), max(x))), max(x)),
    labels = function(x) sprintf("%.1f", x)
  ) +
  geom_circle(
    data = p_mcra_data %>% dplyr::distinct(species_label, R),
    aes(x0 = 0, y0 = 0, r = R),
    inherit.aes = FALSE, color = "black", linewidth = 0.3
  ) +
  facet_wrap(~ species_label, ncol = 5) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic", size = 9),
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.6, "cm"),  # Smaller legend
    legend.key.width = unit(0.4, "cm"),   # Smaller legend
    legend.title.align = 0,  # Left-align title
    legend.justification = "left",  # Left-align legend
    panel.spacing = unit(0.5, "lines"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Sum plot (pmoA+mmoX) with two-line species names
p_sum_data <- p_sum_result$plot$data
p_sum_data$species_label <- sapply(as.character(p_sum_data$species_label), break_species_name)

p_sum <- ggplot(p_sum_data, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = bquote("log"[10] ~ "(pmoA+mmoX)"),
    option = "mako",
    direction = 1,
    breaks = function(x) c(min(x), mean(c(min(x), max(x))), max(x)),
    labels = function(x) sprintf("%.1f", x)
  ) +
  geom_circle(
    data = p_sum_data %>% dplyr::distinct(species_label, R),
    aes(x0 = 0, y0 = 0, r = R),
    inherit.aes = FALSE, color = "black", linewidth = 0.3
  ) +
  facet_wrap(~ species_label, ncol = 5) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic", size = 9),
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.6, "cm"),  # Smaller legend
    legend.key.width = unit(0.4, "cm"),   # Smaller legend
    legend.title.align = 0,  # Left-align title
    legend.justification = "left",  # Left-align legend
    panel.spacing = unit(0.5, "lines"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Overlay plot with two-line species names
overlay_data_modified <- p_overlay$data
overlay_data_modified$species_label <- sapply(as.character(overlay_data_modified$species_label), break_species_name)

p_overlay <- ggplot(overlay_data_modified, aes(x = x, y = y)) +
  geom_raster(aes(fill = Clog_methan), alpha = 0.6, interpolate = TRUE) +
  scale_fill_gradient(
    low = "white", 
    high = "#1E88E5",
    name = bquote("log"[10] ~ "(pmoA+mmoX)"),
    breaks = function(x) {
      min_val <- min(overlay_data_modified$Clog_methan)
      max_val <- max(overlay_data_modified$Clog_methan)
      c(min_val, mean(c(min_val, max_val)), max_val)
    },
    labels = function(x) sprintf("%.1f", x)
  ) +
  new_scale_fill() +
  geom_raster(aes(fill = Clog_mcra), alpha = 0.6, interpolate = TRUE) +
  scale_fill_gradient(
    low = "white", 
    high = "#E53935",
    name = bquote("log"[10] ~ "(mcrA)"),
    breaks = function(x) {
      min_val <- min(overlay_data_modified$Clog_mcra)
      max_val <- max(overlay_data_modified$Clog_mcra)
      c(min_val, mean(c(min_val, max_val)), max_val)
    },
    labels = function(x) sprintf("%.1f", x)
  ) +
  geom_circle(
    data = overlay_data_modified %>% dplyr::distinct(species_label, R),
    aes(x0 = 0, y0 = 0, r = R),
    inherit.aes = FALSE, color = "black", linewidth = 0.5
  ) +
  facet_wrap(~ species_label, ncol = 5) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic", size = 9),
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.6, "cm"),  # Smaller legend
    legend.key.width = unit(0.4, "cm"),   # Smaller legend
    legend.title.align = 0,  # Left-align title
    legend.justification = "left",  # Left-align legend
    legend.box = "vertical",  # Stack legends vertically
    legend.box.just = "left",  # Left-align legend box
    panel.spacing = unit(0.5, "lines"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Combine radial plots vertically
combined_plot <- p_mcra / p_sum / p_overlay

}  # end if (radial plot objects exist)

# ============================================================
# SPECIES COMPARISON 2x2 PLOTS
# Requires: analysis_mcra, analysis_methanotroph, analysis_ratio,
#   cor_area_mcra, cor_area_methanotroph, pearson_ratio,
#   species_comparison_data, y_limit from 01_summary_stats.R
# ============================================================

required_analysis_objs <- c("analysis_mcra", "analysis_methanotroph",
                            "analysis_ratio", "cor_area_mcra",
                            "cor_area_methanotroph", "pearson_ratio",
                            "species_comparison_data", "y_limit")
missing_objs <- required_analysis_objs[!sapply(required_analysis_objs, exists)]

if (length(missing_objs) > 0) {
  cat("WARNING: Missing upstream objects for species comparison plots:\n")
  cat("  ", paste(missing_objs, collapse = ", "), "\n")
  cat("  Run 01_summary_stats.R first to populate these objects.\n")
  cat("  Skipping species comparison and final figure sections.\n\n")
} else {

# Panel 1: Species - mcrA
p_species_mcra_2x2 <- ggplot(analysis_mcra,
                             aes(x = log10(median_mcra + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#C03221",
              fill = "#F5C9C3", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#C03221") +
  geom_text_repel(aes(label = species), size = 3.5, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R\u00B2 = %.3f\np = %.3f",
                           cor_area_mcra$estimate^2,
                           cor_area_mcra$p.value),
           hjust = 1.05, vjust = 1.2, size = 5,
           fill = "white", alpha = 0.9) +
  coord_cartesian(clip = "off") +
  labs(x = expression("log"[10]*" median mcrA"),
       y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_classic(base_size = 11.7) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Panel 2: Species - Methanotrophs
p_species_methanotroph_2x2 <- ggplot(analysis_methanotroph,
                                     aes(x = log10(median_methanotroph + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#4A6FA5",
              fill = "#C5D5E8", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#4A6FA5") +
  geom_text_repel(aes(label = species), size = 3.5, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R\u00B2 = %.3f\np = %.3f",
                           cor_area_methanotroph$estimate^2,
                           cor_area_methanotroph$p.value),
           hjust = 1.05, vjust = 1.2, size = 5,
           fill = "white", alpha = 0.9) +
  coord_cartesian(clip = "off") +
  labs(x = expression("log"[10]*" median (pmoA+mmoX)"),
       y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_classic(base_size = 11.7) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Panel 3: Species - Ratio (with adjusted x-axis title position)
p_species_ratio_2x2 <- ggplot(analysis_ratio,
                              aes(x = median_log_ratio, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#6B5B95",
              fill = "#D7D2E0", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#6B5B95") +
  geom_text_repel(aes(label = species), size = 3.5, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.5) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R\u00B2 = %.3f\np = %.3f",
                           pearson_ratio$estimate^2,
                           pearson_ratio$p.value),
           hjust = 1.05, vjust = 1.2, size = 5,
           fill = "white", alpha = 0.9) +
  coord_cartesian(clip = "off") +
  labs(x = expression("log"[10]*" ratio"),
       y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_classic(base_size = 11.7) +
  theme(
    axis.title = element_text(size = 16),
    axis.title.x = element_text(vjust = 10),  # Move x-axis title up
    axis.text = element_text(size = 14),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Panel 4: Model Comparison (with line breaks in x-axis labels)
p_species_comparison_2x2 <- ggplot(species_comparison_data,
                                   aes(x = Model, y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.2f", R2)),
            vjust = -0.3, size = 5) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#285238"),
                    labels = c("NS", "p < 0.05"),
                    name = "") +
  scale_x_discrete(labels = c("mcrA" = "mcrA", 
                              "pmoA" = "pmoA",
                              "mmoX" = "mmoX",
                              "pmoA+mmoX" = "pmoA+\nmmoX",
                              "Ratio" = "Ratio")) +
  ylim(0, y_limit) +
  labs(x = "", y = expression("Model R"^2)) +
  theme_classic(base_size = 11.7) +
  theme(
    axis.title = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.position = "top",
    legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
    legend.margin = ggplot2::margin(t = 2, r = 0, b = 0, l = 0),
    plot.margin = ggplot2::margin(t = 5, r = 5, b = 2, l = 5)
  )

# Create 2x2 layout
species_2x2_layout <- (p_species_mcra_2x2 | p_species_methanotroph_2x2) /
  (p_species_ratio_2x2 | p_species_comparison_2x2)

}  # end if (analysis objects exist)

# ============================================================
# FINAL SIDE-BY-SIDE LAYOUT
# Requires both combined_plot (radial) and species_2x2_layout
# ============================================================

if (exists("combined_plot") && exists("species_2x2_layout")) {

side_by_side <- (combined_plot | species_2x2_layout) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")",
                  theme = theme(plot.tag = element_text(size = 11, face = "bold")))

side_by_side


# ggsave("outputs/figures/Figure_Radial_And_Species_Comparison.pdf",
#        side_by_side,
#        width = 16,
#        height = 10,
#        limitsize = FALSE)

ggsave("outputs/figures/main/fig8_radial_species_comparison.png",
       side_by_side,
       width = 16.5,
       height = 10,
       dpi = 300,
       limitsize = FALSE)

} else {
  cat("WARNING: Cannot create final figure - missing combined_plot and/or species_2x2_layout.\n")
  cat("  Ensure both 05_radial_gene_plots.R and 01_summary_stats.R have been run.\n")
}

