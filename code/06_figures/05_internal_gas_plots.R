# ==============================================================================
# Internal Gas Concentration Plots (Figures S4-S5)
# ==============================================================================
# Purpose: Internal gas concentration visualizations showing CH4/CO2 inside
#   tree stems.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
# ==============================================================================

library(tidyverse)

tree_all<-read.csv('data/processed/integrated/merged_tree_dataset_final.csv')

# Species name mapping
species_mapping <- c(
  "ACRU" = "Acer rubrum",
  "ACSA" = "Acer saccharum", 
  "BEAL" = "Betula alleghaniensis",
  "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera",
  "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana",
  "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra",
  "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata",
  "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina",
  "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina",
  "SAAL" = "Sassafras albidum"
)

library(ggplot2)
library(ggridges)
library(dplyr)
library(viridis)


library(ggplot2)
library(ggridges)
library(dplyr)
library(viridis)

# Prepare data - filter for valid CH4 and O2 concentrations
ridge_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(O2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_mapping[species_id]) %>%
  filter(!is.na(Species_Latin)) %>%
  group_by(Species_Latin) %>%
  #filter(n() > 2) %>%  # Only keep species with >2 observations
  ungroup() %>%
  mutate(Species_Latin = reorder(Species_Latin, CH4_concentration, mean, na.rm = TRUE))

# Calculate species statistics for mean points
species_stats <- ridge_data %>%
  group_by(Species_Latin) %>%
  summarise(
    mean_ch4 = mean(CH4_concentration, na.rm = TRUE),
    se_ch4 = sd(CH4_concentration, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

# Create ridgeline plot
p_ch4_ridges <- ggplot(ridge_data, aes(x = CH4_concentration, y = Species_Latin)) +
  geom_density_ridges(
    fill = "lightgray",
    alpha = 0.6,
    scale = 1.5,
    rel_min_height = 0.01,
    color = "gray30",
    linewidth = 0.3
  ) +
  geom_point(
    aes(color = O2_concentration),
    position = position_jitter(width = 0, height = 0.15),
    size = 1.5,
    alpha = 0.6
  ) +
  geom_pointrange(
    data = species_stats,
    aes(x = mean_ch4, 
        xmin = mean_ch4 - se_ch4, 
        xmax = mean_ch4 + se_ch4,
        y = Species_Latin),
    color = "black",
    size = 0.75,
    linewidth = 0.75,
    fatten = 2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(
    name = expression(O[2]~concentration),
    option = "plasma",
    labels = function(x) sprintf("%.1f%%", x / 10000),  # Convert ppm to %
    guide = guide_colorbar(
      position = "bottom",
      barwidth = 10,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = function(x) sprintf("%.1f%%", x / 10000)  # Convert ppm to %
  ) +
  labs(
    x = expression(CH[4]~concentration),
    y = ""
  ) +
  theme_ridges(grid = TRUE) +
  theme(
    axis.title.x = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "gray95", linewidth = 0.2)
  )

print(p_ch4_ridges)



library(ggplot2)
library(dplyr)
library(viridis)

# Custom labeling function for x-axis
format_mixed_labels <- function(x) {
  pct <- x / 10000
  ifelse(pct < 1,
         sprintf("%g ppm", x),      # Show as ppm for values < 1%
         sprintf("%.1f%%", pct))    # Show as % for values >= 1%
}

p_facet <- ggplot(ridge_data, aes(x = CH4_concentration)) +
  geom_density(fill = "lightgray", alpha = 0.6) +
  geom_point(aes(y = 0, color = O2_concentration), 
             position = position_jitter(height = 0.01),
             size = 2, alpha = 0.8) +
  geom_vline(data = species_stats, aes(xintercept = mean_ch4),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_gradientn(
    colors = c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"),
    values = scales::rescale(c(0, 120000, 150000, 200000, 205000, 220000)),
    name = expression(O[2]~concentration),
    labels = function(x) sprintf("%.1f%%", x / 10000)
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = format_mixed_labels
  ) +
  facet_wrap(~Species_Latin, ncol = 3, scales = "free_y") +
  labs(x = expression(CH[4]~concentration), y = "Density") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "italic", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
p_facet



p_violin <- ggplot(ridge_data, aes(x = Species_Latin, y = CH4_concentration)) +
  geom_violin(fill = "lightgray", alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = O2_concentration), 
              width = 0.15, size = 2.5, alpha = 0.7) +
  scale_color_gradientn(
    colors = c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"),
    values = scales::rescale(c(0, 120000, 150000, 200000, 205000, 220000)),
    name = expression(O[2]~concentration),
    labels = function(x) sprintf("%.1f%%", x / 10000)
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
    legend.position = "bottom"
  )
p_violin


library(ggbeeswarm)

p_beeswarm <- ggplot(ridge_data, aes(x = Species_Latin, y = CH4_concentration)) +
  geom_quasirandom(aes(color = O2_concentration), 
                   size = 2.5, alpha = 0.75, width = 0.3) +
  geom_pointrange(data = species_stats,
                  aes(x = Species_Latin, y = mean_ch4,
                      ymin = mean_ch4 - se_ch4, ymax = mean_ch4 + se_ch4),
                  color = "black", size = 0.5, linewidth = 1) +
  scale_color_gradientn(
    colors = rev(c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")),
    values = scales::rescale(c(0, 170000, 180000, 190000, 200000, 220000)),
    name = expression(O[2]~concentration),
    labels = function(x) sprintf("%.1f%%", x / 10000)
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
    legend.position = "bottom"
  )
p_beeswarm




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
p_beeswarm

# Save beeswarm figure (Fig S4)
ggsave("outputs/figures/supplementary/figS4_internal_gas_beeswarm.png",
       p_beeswarm, width = 10, height = 8, dpi = 300)











library(ggplot2)
library(dplyr)
library(patchwork)
library(broom)
library(ape)
library(phytools)
library(viridis)
library(RColorBrewer)

# [Keep all the data preparation code the same up to the plots]

# Prepare data for CH4 conc vs CO2 conc
conc_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(CO2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_mapping[species_id]) %>%
  filter(!is.na(Species_Latin))

# Prepare data for CH4 conc vs O2 conc
o2_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(O2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_mapping[species_id]) %>%
  filter(!is.na(Species_Latin))

# Prepare data for heartwood mcrA vs CH4 conc
mcra_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(ddpcr_mcra_probe_Inner_loose) & !is.na(species_id)) %>%
  mutate(
    Species_Latin = species_mapping[species_id],
    mcra_copies = ddpcr_mcra_probe_Inner_loose + 1
  ) %>%
  filter(!is.na(Species_Latin))

# Prepare data for CH4 flux vs CH4 conc
flux_conc_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_mapping[species_id]) %>%
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

format_mixed_labels <- function(x) {
  pct <- x / 10000
  ifelse(pct < 1, sprintf("%g ppm", x), sprintf("%.1f%%", pct))
}

# --- PHYLOGENETIC COLOR CALCULATION ---
tree_file <- "data/processed/metadata/PhytoPhylo"

all_species <- unique(c(conc_data$Species_Latin, o2_data$Species_Latin, 
                        flux_conc_data$Species_Latin, mcra_data$Species_Latin))
known_species <- all_species[!grepl("Unknown", all_species)]
unknown_species <- all_species[grepl("Unknown", all_species)]

phylo_colors <- NULL
phylo_distance_data <- NULL

if (file.exists(tree_file) && length(known_species) > 1) {
  tryCatch({
    tree_scenario1 <- read.tree(tree_file)
    present_species_clean <- gsub(" ", "_", known_species)
    tree_species <- intersect(present_species_clean, tree_scenario1$tip.label)
    
    if (length(tree_species) > 1) {
      pruned_tree <- keep.tip(tree_scenario1, tree_species)
      dist_matrix <- cophenetic(pruned_tree)
      average_distances <- rowMeans(dist_matrix)
      
      distance_df <- data.frame(
        species_clean = names(average_distances),
        distance = as.numeric(average_distances)
      ) %>%
        mutate(species = gsub("_", " ", species_clean)) %>%
        arrange(distance) %>%
        mutate(
          distance_group = cut(distance, 
                               breaks = quantile(distance, probs = seq(0, 1, length.out = min(8, length(distance)))), 
                               include.lowest = TRUE, 
                               labels = FALSE),
          color = viridis(max(distance_group, na.rm = TRUE))[distance_group]
        )
      
      phylo_colors <- setNames(distance_df$color, distance_df$species)
      phylo_distance_data <- distance_df %>%
        dplyr::select(species, distance, distance_group) %>%
        arrange(distance)
    }
  }, error = function(e) {
    cat("Warning: Could not calculate phylogenetic distances:", e$message, "\n")
  })
}

if (is.null(phylo_colors)) {
  if(length(known_species) <= 12) {
    known_colors <- RColorBrewer::brewer.pal(min(max(3, length(known_species)), 12), "Set3")[1:length(known_species)]
  } else {
    known_colors <- rainbow(length(known_species))
  }
  names(known_colors) <- known_species
  phylo_colors <- known_colors
}

unknown_colors <- rep("grey80", length(unknown_species))
names(unknown_colors) <- unknown_species
final_colors <- c(phylo_colors, unknown_colors)

legend_order <- if (!is.null(phylo_distance_data)) {
  c(phylo_distance_data$species, unknown_species)
} else {
  c(sort(known_species), unknown_species)
}

# Calculate correlations - CO2 (MODIFIED: only R² and p-value)
pearson1 <- cor.test(log10(conc_data$CH4_concentration + 1), 
                     log10(conc_data$CO2_concentration + 1), method = "pearson")
model1 <- lm(log10(CH4_concentration + 1) ~ log10(CO2_concentration + 1), data = conc_data)
stats1 <- glance(model1)
label1 <- sprintf("atop(R^2 == %.3f, p~'%s')",
                  stats1$r.squared,
                  ifelse(pearson1$p.value < 0.001, "< 0.001", sprintf("= %.3f", pearson1$p.value)))

# Calculate correlations - O2 (MODIFIED: only R² and p-value)
pearson2 <- cor.test(log10(o2_data$CH4_concentration + 1), 
                     log10(o2_data$O2_concentration + 1), method = "pearson")
model2 <- lm(log10(CH4_concentration + 1) ~ log10(O2_concentration + 1), data = o2_data)
stats2 <- glance(model2)
label2 <- sprintf("atop(R^2 == %.3f, p~'%s')",
                  stats2$r.squared,
                  ifelse(pearson2$p.value < 0.001, "< 0.001", sprintf("= %.3f", pearson2$p.value)))

# Calculate correlations - mcrA (MODIFIED: only R² and p-value)
pearson_mcra <- cor.test(log10(mcra_data$mcra_copies), 
                         log10(mcra_data$CH4_concentration + 1), method = "pearson")
model_mcra <- lm(log10(CH4_concentration + 1) ~ log10(mcra_copies), data = mcra_data)
stats_mcra <- glance(model_mcra)
label_mcra <- sprintf("atop(R^2 == %.3f, p~'%s')",
                      stats_mcra$r.squared,
                      ifelse(pearson_mcra$p.value < 0.001, "< 0.001", sprintf("= %.3f", pearson_mcra$p.value)))

# Calculate stats for each height (MODIFIED: only R² and p-value)
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


# ---- FORCE IDENTICAL LEVELS + PALETTE ----
species_levels <- legend_order                  # from your code above
species_palette <- final_colors[species_levels] # ensure names align
names(species_palette) <- species_levels

# Coerce Species_Latin to the same factor in EVERY dataset
conc_data$Species_Latin     <- factor(conc_data$Species_Latin, levels = species_levels)
o2_data$Species_Latin       <- factor(o2_data$Species_Latin,   levels = species_levels)
mcra_data$Species_Latin     <- factor(mcra_data$Species_Latin, levels = species_levels)
flux_conc_data$Species_Latin<- factor(flux_conc_data$Species_Latin, levels = species_levels)


library(cowplot)  # add this

# --- keep the factor/palette forcing from the previous message ---

# (Re)build plots exactly as before BUT with legends OFF in all panels
p1_core <- ggplot(conc_data, aes(x = CO2_concentration, y = CH4_concentration)) +
  geom_point(aes(color = Species_Latin), size = 2.5, alpha = 0.7, show.legend = TRUE) +  # we'll use this to pull legend
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
  labs(x = expression("Heartwood mcrA (copies g"^-1*")"),
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

# Add tag theme to individual plots for ggplot2 4.x compatibility
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
final_plot <- cowplot::plot_grid(core_grid, legend_grob, ncol = 2, rel_widths = c(1, 0.22))

print(final_plot)

# Save combined internal gas figure (Fig S5)
ggsave("outputs/figures/supplementary/figS5_internal_gas_profiles.png",
       final_plot, width = 12, height = 7.5, dpi = 300)