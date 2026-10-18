# ==============================================================================
# figS15-16_internal-gas.R -- Figures S16 and S17: internal stem gases
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
# Figure S16  internal CH4 by species (quasirandom points coloured by O2; mean +/- SE),
#             species ordered by mean CH4, pseudo-log axis
# Figure S17  (a) CO2 against CH4, (b) O2 against CH4, (c) heartwood mcrA against CH4,
#             (d) stem flux at 50, 125 and 200 cm against CH4; linear fits with R2 and
#             p; species coloured in phylogenetic order (code/lib/phylo_palette.R)
# Inputs: data/tables/ via code/lib/tables_read.R: survey_2021_trees() (internal_gas
#   2021 whole-stem samples, flux_stem 2021 heights, ddpcr + samples, trees, species);
#   the PhytoPhylo megatree (data/raw/external/phylogeny/PhytoPhylo)
# Outputs: outputs/figures/original/supplementary/figS7_internal_gas_beeswarm.png  (S16)
#          outputs/figures/original/supplementary/figS8_internal_gas_profiles.png  (S17)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(ggbeeswarm); library(patchwork); library(broom); library(cowplot)
})
source("code/lib/tables_read.R")
source("code/lib/phylo_palette.R")
tree_all <- survey_2021_trees()   # 2021 survey per tree; species from trees.csv

# Axis labels: ppm below 1%, % above
format_mixed_labels <- function(x) {
  pct <- x / 10000
  ifelse(pct < 1,
         sprintf("%g ppm", x),      # Show as ppm for values < 1%
         sprintf("%.1f%%", pct))    # Show as % for values >= 1%
}

# ==============================================================================
# Figure S16: internal CH4 by species
# ==============================================================================
ridge_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(O2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_name(species_id)) %>%
  filter(!is.na(Species_Latin)) %>%
  mutate(Species_Latin = reorder(Species_Latin, CH4_concentration, mean, na.rm = TRUE))

species_stats <- ridge_data %>%
  group_by(Species_Latin) %>%
  summarise(
    mean_ch4 = mean(CH4_concentration, na.rm = TRUE),
    se_ch4 = sd(CH4_concentration, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

p_beeswarm <- ggplot(ridge_data, aes(x = Species_Latin, y = CH4_concentration)) +
  geom_quasirandom(aes(color = O2_concentration), 
                   size = 2.5, alpha = 0.75, width = 0.3) +
  geom_pointrange(data = species_stats,
                  aes(x = Species_Latin, y = mean_ch4,
                      ymin = mean_ch4 - se_ch4, ymax = mean_ch4 + se_ch4),
                  color = "black", size = 0.5, linewidth = 1) +
  scale_color_gradientn(
    colors = rev(c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")),
    values = scales::rescale(c(0, 100000, 160000, 190000, 200000, 220000)),
    name = expression(O[2]~concentration),
    labels = function(x) sprintf("%.1f%%", x / 10000),
    guide = guide_colorbar(
      barwidth = 10,
      barheight = 0.8,
      title.position = "top",
      title.hjust = 0.5,
      ticks.colour = "black"
    )
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = format_mixed_labels
  ) +
  coord_flip() +
  labs(y = expression(CH[4]~concentration), x = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom"
  )
ggsave("outputs/figures/original/supplementary/figS7_internal_gas_beeswarm.png",
       p_beeswarm, width = 10, height = 8, dpi = 300)

# ==============================================================================
# Figure S17: internal gases, mcrA and flux
# ==============================================================================
# Prepare data for CH4 conc vs CO2 conc
conc_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(CO2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_name(species_id)) %>%
  filter(!is.na(Species_Latin))

# Prepare data for CH4 conc vs O2 conc
o2_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(O2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_name(species_id)) %>%
  filter(!is.na(Species_Latin))

# Prepare data for heartwood mcrA vs CH4 conc
mcra_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(ddpcr_mcra_probe_Inner_loose) & !is.na(species_id)) %>%
  mutate(
    Species_Latin = species_name(species_id),
    mcra_copies = ddpcr_mcra_probe_Inner_loose + 1
  ) %>%
  filter(!is.na(Species_Latin))

# Prepare data for CH4 flux vs CH4 conc
flux_conc_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_name(species_id)) %>%
  filter(!is.na(Species_Latin)) %>%
  pivot_longer(
    cols = c(CH4_best.flux_50cm, CH4_best.flux_125cm, CH4_best.flux_200cm),
    names_to = "height",
    values_to = "CH4_flux",
    values_drop_na = TRUE
  ) %>%
  mutate(
    height_clean = case_when(
      height == "CH4_best.flux_50cm" ~ "50 cm",
      height == "CH4_best.flux_125cm" ~ "125 cm",
      height == "CH4_best.flux_200cm" ~ "200 cm",
      TRUE ~ height
    ),
    height_clean = factor(height_clean, levels = c("50 cm", "125 cm", "200 cm"))
  )


# species colours in phylogenetic order; names with "Unknown" grey and last
phylo <- phylo_species_colours(c(conc_data$Species_Latin, o2_data$Species_Latin,
                                 flux_conc_data$Species_Latin, mcra_data$Species_Latin))

# Calculate correlations - CO2
pearson1 <- cor.test(log10(conc_data$CH4_concentration + 1), 
                     log10(conc_data$CO2_concentration + 1), method = "pearson")
model1 <- lm(log10(CH4_concentration + 1) ~ log10(CO2_concentration + 1), data = conc_data)
stats1 <- glance(model1)
label1 <- sprintf("atop(R^2 == %.3f, p~'%s')",
                  stats1$r.squared,
                  ifelse(pearson1$p.value < 0.001, "< 0.001", sprintf("= %.3f", pearson1$p.value)))

# Calculate correlations - O2
pearson2 <- cor.test(log10(o2_data$CH4_concentration + 1), 
                     log10(o2_data$O2_concentration + 1), method = "pearson")
model2 <- lm(log10(CH4_concentration + 1) ~ log10(O2_concentration + 1), data = o2_data)
stats2 <- glance(model2)
label2 <- sprintf("atop(R^2 == %.3f, p~'%s')",
                  stats2$r.squared,
                  ifelse(pearson2$p.value < 0.001, "< 0.001", sprintf("= %.3f", pearson2$p.value)))

# Calculate correlations - mcrA
pearson_mcra <- cor.test(log10(mcra_data$mcra_copies), 
                         log10(mcra_data$CH4_concentration + 1), method = "pearson")
model_mcra <- lm(log10(CH4_concentration + 1) ~ log10(mcra_copies), data = mcra_data)
stats_mcra <- glance(model_mcra)
label_mcra <- sprintf("atop(R^2 == %.3f, p~'%s')",
                      stats_mcra$r.squared,
                      ifelse(pearson_mcra$p.value < 0.001, "< 0.001", sprintf("= %.3f", pearson_mcra$p.value)))

# Calculate stats for each height
flux_stats <- flux_conc_data %>%
  group_by(height_clean) %>%
  do({
    data <- .
    model <- lm(log10(abs(CH4_flux) + 0.01) ~ log10(CH4_concentration + 1), data = data)
    pearson_test <- cor.test(log10(data$CH4_concentration + 1), 
                             log10(abs(data$CH4_flux) + 0.01), method = "pearson")
    stats <- glance(model)
    data.frame(
      height_clean = unique(data$height_clean),
      label = sprintf("atop(R^2 == %.3f, p~'%s')",
                      stats$r.squared,
                      ifelse(pearson_test$p.value < 0.001, "< 0.001", sprintf("= %.3f", pearson_test$p.value)))
    )
  }) %>%
  ungroup()


# ---- one species order and palette for every panel ----
species_levels <- phylo$order
species_palette <- phylo$colours[species_levels]
names(species_palette) <- species_levels

# Coerce Species_Latin to the same factor in EVERY dataset
conc_data$Species_Latin     <- factor(conc_data$Species_Latin, levels = species_levels)
o2_data$Species_Latin       <- factor(o2_data$Species_Latin,   levels = species_levels)
mcra_data$Species_Latin     <- factor(mcra_data$Species_Latin, levels = species_levels)
flux_conc_data$Species_Latin<- factor(flux_conc_data$Species_Latin, levels = species_levels)

# Panels without legends; the species legend is taken from (a) and placed at the right
p1_core <- ggplot(conc_data, aes(x = CO2_concentration, y = CH4_concentration)) +
  geom_point(aes(color = Species_Latin), size = 2.5, alpha = 0.7, show.legend = TRUE) +  # legend source
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, label = label1, hjust = -0.1, vjust = 1.5, parse = TRUE) +
  scale_color_manual(values = species_palette, limits = species_levels, breaks = species_levels,
                     drop = FALSE, name = "Species\n(phylogenetic order)") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                     breaks = c(0, 1000, 10000, 50000, 100000, 200000),
                     labels = format_mixed_labels) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                     breaks = c(0, 10, 100, 1000, 10000, 100000),
                     labels = format_mixed_labels) +
  labs(x = expression(CO[2]~concentration), y = expression(CH[4]~concentration)) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 8, face = "italic"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# Extract legend from p1_core
legend_grob <- cowplot::get_legend(p1_core)

# Now make no-legend versions of all plots
p1 <- p1_core + theme(legend.position = "none")

p2 <- ggplot(o2_data, aes(x = O2_concentration, y = CH4_concentration)) +
  geom_point(aes(color = Species_Latin), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, label = label2, hjust = -0.1, vjust = 1.5, parse = TRUE) +
  scale_color_manual(values = species_palette, limits = species_levels, breaks = species_levels, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 220000, by = 50000),
                     labels = function(x) sprintf("%.1f%%", x / 10000)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                     breaks = c(0, 10, 100, 1000, 10000, 100000),
                     labels = format_mixed_labels) +
  labs(x = expression(O[2]~concentration), y = expression(CH[4]~concentration)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

p3 <- ggplot(mcra_data, aes(x = mcra_copies, y = CH4_concentration)) +
  geom_point(aes(color = Species_Latin), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1, show.legend = FALSE) +
  annotate("text", x = 0, y = Inf, label = label_mcra, hjust = -0.1, vjust = 1.5, parse = TRUE) +
  scale_color_manual(values = species_palette, limits = species_levels, breaks = species_levels, drop = FALSE) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                     breaks = c(0, 10, 100, 1000, 10000, 100000),
                     labels = format_mixed_labels) +
  labs(x = expression("Heartwood "*italic(mcrA)*" (copies g"^-1*")"),
       y = expression(CH[4]~concentration)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

p4 <- ggplot(flux_conc_data, aes(x = CH4_concentration, y = CH4_flux)) +
  geom_point(aes(color = Species_Latin), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", show.legend = FALSE) +
  geom_text(data = flux_stats, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, show.legend = FALSE, parse = TRUE) +
  facet_wrap(~ height_clean, ncol = 3) +
  scale_color_manual(values = species_palette, limits = species_levels, breaks = species_levels, drop = FALSE) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                     breaks = c(0, 10, 100, 1000, 10000, 100000),
                     labels = format_mixed_labels) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.01),
                     breaks = c(-0.1, 0, 0.1, 0.5, 1, 2, 5)) +
  labs(x = expression(CH[4]~concentration),
       y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

# Tag theme on each plot (ggplot2 4.x)
tag_theme <- theme(plot.tag = element_text(size = 11, face = "bold"))
p1 <- p1 + tag_theme
p2 <- p2 + tag_theme
p3 <- p3 + tag_theme
p4 <- p4 + tag_theme

# Arrange the plots with patchwork (no legends inside)
core_grid <- ((p1 | p2 | p3) / p4) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")")

# Place the extracted legend to the right of the whole patchwork
final_plot <- cowplot::ggdraw() + cowplot::draw_plot(cowplot::plot_grid(core_grid, legend_grob, ncol = 2, rel_widths = c(1, 0.22))) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("outputs/figures/original/supplementary/figS8_internal_gas_profiles.png",
       final_plot, width = 12, height = 7.5, dpi = 300, bg = "white")
