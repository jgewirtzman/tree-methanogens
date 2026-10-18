# ==============================================================================
# fig08_radial-species.R -- Figure 8 (radial gene distributions and species-level
#                           gene-flux relationships) and Figure S24 (mcrA vs
#                           methanotrophs)
# ------------------------------------------------------------------------------
# Figure 8
#   (a-c) species-averaged stem cross-sections: mcrA, pmoA + mmoX, and their log10
#         ratio (code/07_molecular/helper_radial_gene_plots.R); every cross-section at
#         the same size
#   (d-f) species median flux against species median mcrA, pmoA + mmoX and the
#         mcrA : methanotroph ratio; (g) R2 of the five species-level gene models
#         (code/07_molecular/helper_scale_dependent_gene_patterns.R)
# Figure S24: log10 (pmoA + mmoX) against log10 mcrA, per tree and per species.
# Reads data/tables/ through the two helpers (run from the repository root).
# Outputs: outputs/figures/original/main/fig8_radial_species_comparison.png,
#          outputs/figures/original/supplementary/figS14_mcra_vs_methanotroph.png
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(ggrepel)
library(patchwork)

source("code/07_molecular/helper_scale_dependent_gene_patterns.R")
source("code/07_molecular/helper_radial_gene_plots.R")
RADIAL_NCOL <- 10
# every cross-section drawn at the same size (not scaled to DBH)
unit_circle <- function(d) dplyr::mutate(d, x = x / R, y = y / R, R = 1)

# species names on two lines
break_species_name <- function(name) {
  parts <- strsplit(name, " ")[[1]]
  if (length(parts) == 2) {
    return(paste0(parts[1], "\n", parts[2]))
  }
  return(name)
}

# ============================================================
# (a-c) RADIAL CROSS-SECTIONS
# ============================================================
p_mcra_data <- mcra_cloud
p_mcra_data$species_label <- sapply(as.character(p_mcra_data$species_label), break_species_name)
p_mcra_data <- unit_circle(p_mcra_data)

p_mcra <- ggplot(p_mcra_data, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = bquote(atop("log"[10], "(" * italic(mcrA) * ")")),
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
  facet_wrap(~ species_label, ncol = RADIAL_NCOL) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic", size = 8.6), strip.clip = "off",
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm"),  # Smaller legend
    legend.key.width = unit(0.4, "cm"),   # Smaller legend
    legend.title.align = 0,  # Left-align title
    legend.justification = "left",  # Left-align legend
    legend.box.spacing = unit(2, "pt"),
    panel.spacing = unit(0.5, "lines"),
    panel.background = element_rect(fill = "white", color = NA)
  )

p_sum_data <- sum_cloud
p_sum_data$species_label <- sapply(as.character(p_sum_data$species_label), break_species_name)
p_sum_data <- unit_circle(p_sum_data)

p_sum <- ggplot(p_sum_data, aes(x = x, y = y, fill = Clog)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = bquote(atop("log"[10], "(" * italic(pmoA) * " + " * italic(mmoX) * ")")),
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
  facet_wrap(~ species_label, ncol = RADIAL_NCOL) +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic", size = 8.6), strip.clip = "off",
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm"),  # Smaller legend
    legend.key.width = unit(0.4, "cm"),   # Smaller legend
    legend.title.align = 0,  # Left-align title
    legend.justification = "left",  # Left-align legend
    legend.box.spacing = unit(2, "pt"),
    panel.spacing = unit(0.5, "lines"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# (c) one diverging scale: log10(mcrA / (pmoA + mmoX)); red = methanogen-dominated,
# white = balanced, blue = methanotroph-dominated
overlay_data_modified <- overlay_cloud
overlay_data_modified$species_label <- sapply(as.character(overlay_data_modified$species_label), break_species_name)
overlay_data_modified <- unit_circle(overlay_data_modified)
rat <- overlay_data_modified %>% mutate(ratio = Clog_mcra - Clog_methan)
lim <- max(abs(range(rat$ratio, na.rm = TRUE)))
p_overlay <- ggplot(rat, aes(x = x, y = y, fill = ratio)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                       name = bquote(atop("log"[10] ~ "(" * italic(mcrA) * " :", "methanotroph)")),
                       labels = function(x) sub("-", "\u2212", sprintf("%.1f", x))) +
  geom_circle(data = rat %>% dplyr::distinct(species_label, R), aes(x0 = 0, y0 = 0, r = R),
              inherit.aes = FALSE, color = "black", linewidth = 0.3) +
  facet_wrap(~ species_label, ncol = RADIAL_NCOL) + coord_equal() +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), strip.text = element_text(face = "italic", size = 8.6), strip.clip = "off",
        strip.background = element_blank(), legend.position = "right",
        legend.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.key.height = unit(0.6, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.justification = "left", panel.spacing = unit(0.5, "lines"),
        panel.background = element_rect(fill = "white", color = NA))

# ============================================================
# (d-g) SPECIES-LEVEL GENE-FLUX PANELS
# ============================================================
# Panel 1: Species - mcrA
p_species_mcra_2x2 <- ggplot(analysis_mcra,
                             aes(x = log10(median_mcra + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#C03221",
              fill = "#F5C9C3", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#C03221") +
  geom_text_repel(aes(label = species), size = 3.5, fontface = "italic", seed = 42, min.segment.length = 0.1, force = 8, force_pull = 0.5,
                  max.time = 5, max.iter = 1e5, box.padding = 0.35, point.padding = 0.25, max.overlaps = Inf) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = log10(c(10, 15, 20, 30, 50) * 1e3), labels = c("10", "15", "20", "30", "50")) +
  labs(x = expression(atop("Median "*italic(mcrA), "("*10^3*" copies g"^-1*")")),
       y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_classic(base_size = 11.7) +
  theme(
    axis.title = element_text(size = 12.5),
    axis.text = element_text(size = 12), plot.subtitle = element_text(size = 12),
    plot.margin = margin(5, 30, 5, 5),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Panel 2: Species - Methanotrophs
p_species_methanotroph_2x2 <- ggplot(analysis_methanotroph,
                                     aes(x = log10(median_methanotroph + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#4A6FA5",
              fill = "#C5D5E8", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#4A6FA5") +
  geom_text_repel(aes(label = species), size = 3.5, fontface = "italic", seed = 42, min.segment.length = 0.1, force = 8, force_pull = 0.5,
                  max.time = 5, max.iter = 1e5, box.padding = 0.35, point.padding = 0.25, max.overlaps = Inf) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = log10(c(30, 50, 100, 200, 400) * 1e3), labels = c("30", "50", "100", "200", "400")) +
  labs(x = expression(atop("Median "*italic(pmoA)*" + "*italic(mmoX), "("*10^3*" copies g"^-1*")")),
       y = NULL) +
  theme_classic(base_size = 11.7) +
  theme(
    axis.title = element_text(size = 12.5),
    axis.text = element_text(size = 12), plot.subtitle = element_text(size = 12),
    plot.margin = margin(5, 30, 5, 5),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Panel 3: Species - Ratio (with adjusted x-axis title position)
p_species_ratio_2x2 <- ggplot(analysis_ratio,
                              aes(x = median_log_ratio, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#6B5B95",
              fill = "#D7D2E0", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#6B5B95") +
  geom_text_repel(aes(label = species), size = 3.5, fontface = "italic", seed = 42, min.segment.length = 0.1, force = 8, force_pull = 0.5,
                  max.time = 5, max.iter = 1e5, box.padding = 0.35, point.padding = 0.25, max.overlaps = Inf) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.5) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = log10(c(0.05, 0.1, 0.2, 0.5, 1)), labels = c("0.05", "0.1", "0.2", "0.5", "1")) +
  labs(x = expression(atop(italic(mcrA)*" : methanotroph", "ratio")),
       y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_classic(base_size = 11.7) +
  theme(
    axis.title = element_text(size = 12.5),
    axis.text = element_text(size = 12), plot.subtitle = element_text(size = 12),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# ============================================================
# FIGURE 8: cross-sections on top, species panels in a row below
# ============================================================
bm <- theme(plot.margin = margin(16, 8, 5, 5), axis.title = element_text(size = 11.5))
# (g) as horizontal bars: model names read horizontally without adding height under the panel
g_h <- ggplot(species_comparison_data, aes(x = R2, y = Model, fill = Significant)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", R2)), hjust = -0.15, size = 4) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#285238"), labels = c("NS", "p < 0.05"), name = "") +
  scale_y_discrete(labels = c("mcrA" = expression(italic(mcrA)), "pmoA" = expression(italic(pmoA)),
                              "mmoX" = expression(italic(mmoX)),
                              "pmoA+mmoX" = expression(italic(pmoA)*"+"*italic(mmoX)), "Ratio" = "Ratio")) +
  scale_x_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression("Model R"^2), y = NULL) +
  theme_classic(base_size = 11.7) +
  theme(axis.text = element_text(size = 12), legend.position = "bottom", legend.text = element_text(size = 12),
        legend.margin = margin(0, 0, 0, 0))
bottom_row <- wrap_plots(p_species_mcra_2x2 + bm, p_species_methanotroph_2x2 + bm,
                        p_species_ratio_2x2 + labs(y = NULL) + bm, g_h + bm,
                        nrow = 1, widths = c(1, 1, 1, 1))
no_strip <- theme(strip.text = element_blank())
combined_plot <- p_mcra / (p_sum + no_strip) / (p_overlay + no_strip)   # species named once, above row (a)
combined_plot <- combined_plot & theme(legend.key.height = unit(0.26, "cm"), legend.title = element_text(size = 10.5),
                                       legend.text = element_text(size = 9.5))
side_by_side <- wrap_plots(combined_plot, bottom_row, ncol = 1, heights = c(2.6, 2.2))
side_by_side <- side_by_side +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")",
                  theme = theme(plot.tag = element_text(size = 11, face = "bold"))) &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave("outputs/figures/original/main/fig8_radial_species_comparison.png",
       side_by_side, width = 14, height = 9.5, dpi = 300, limitsize = FALSE)

# ============================================================
# FIGURE S24: mcrA vs (pmoA + mmoX), tree and species level
# ============================================================
# ------------------------------------------------------------
# INDIVIDUAL TREE LEVEL: mcrA vs methanotroph_total
# ------------------------------------------------------------

# Calculate correlation
cor_tree_mcra_methanotroph <- cor.test(
  tree_level_complete$log_tree_mcra,
  tree_level_complete$log_methanotroph
)

# Create plot (NO trendline - not significant)
p_tree_mcra_vs_methanotroph <- ggplot(tree_level_complete,
                                      aes(y = log_methanotroph,
                                          x = log_tree_mcra,
                                          color = species)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = species_palette, name = "Species") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("r = %.3f\np = %.3f",
                           cor_tree_mcra_methanotroph$estimate,
                           cor_tree_mcra_methanotroph$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  labs(
    title = "Individual Tree Level",
    y = expression(log[10]~(pmoA+mmoX)),
    x = expression(log[10]~mcrA)
  ) +
  theme_pub_gene +
  theme(
    legend.position = "right",
    legend.text = element_text(face = "italic", size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ------------------------------------------------------------
# SPECIES LEVEL: median mcrA vs median methanotroph
# ------------------------------------------------------------

# Combine mcrA and methanotroph data for species-level analysis
species_mcra_methanotroph <- analysis_mcra %>%
  dplyr::select(species, species_id, median_mcra, q25_mcra, q75_mcra) %>%
  left_join(
    analysis_methanotroph %>% 
      dplyr::select(species, median_methanotroph, q25_methanotroph, q75_methanotroph),
    by = "species"
  )

# Calculate correlation
cor_species_mcra_methanotroph <- cor.test(
  log10(species_mcra_methanotroph$median_mcra + 1),
  log10(species_mcra_methanotroph$median_methanotroph + 1)
)

# Create plot (NO trendline - not significant)
p_species_mcra_vs_methanotroph <- ggplot(species_mcra_methanotroph,
                                         aes(y = log10(median_methanotroph + 1),
                                             x = log10(median_mcra + 1),
                                             color = species)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = species_palette) +
  geom_text_repel(aes(label = species), size = 3, fontface = "italic", seed = 42,
                  box.padding = 0.3, max.overlaps = 20) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("r = %.3f\np = %.3f",
                           cor_species_mcra_methanotroph$estimate,
                           cor_species_mcra_methanotroph$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  labs(
    title = "Species Level",
    y = expression(log[10]~median~(pmoA+mmoX)),
    x = expression(log[10]~median~mcrA)
  ) +
  theme_pub_gene +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# ------------------------------------------------------------
# COMBINE AND SAVE
# ------------------------------------------------------------

# Side-by-side layout with panel labels
mcra_vs_methanotroph_layout <- (p_tree_mcra_vs_methanotroph |
  p_species_mcra_vs_methanotroph) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")",
                  theme = theme(plot.tag = element_text(size = 11, face = "bold")))

ggsave("outputs/figures/original/supplementary/figS14_mcra_vs_methanotroph.png", mcra_vs_methanotroph_layout,
       width = 13, height = 5.5, dpi = 300)


# Print summary
cat("\n### mcrA vs METHANOTROPH CORRELATION SUMMARY ###\n")
cat("\nIndividual Tree Level:\n")
cat(sprintf("  Pearson r = %.3f\n", cor_tree_mcra_methanotroph$estimate))
cat(sprintf("  p-value = %.4f\n", cor_tree_mcra_methanotroph$p.value))
cat(sprintf("  n = %d trees\n", nrow(tree_level_complete)))

cat("\nSpecies Level:\n")
cat(sprintf("  Pearson r = %.3f\n", cor_species_mcra_methanotroph$estimate))
cat(sprintf("  p-value = %.4f\n", cor_species_mcra_methanotroph$p.value))
cat(sprintf("  n = %d species\n", nrow(species_mcra_methanotroph)))
