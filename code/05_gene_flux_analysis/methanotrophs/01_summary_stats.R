# ==============================================================================
# Methanotroph Summary Statistics
# ==============================================================================
# Purpose: Summary statistics for methanotroph gene abundances across species
#   and compartments.
#
# Pipeline stage: 3 — Analysis
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
# ==============================================================================

cat("\n============================================================\n")
cat("CREATING 3x2 SIDE-BY-SIDE LAYOUT FIGURE\n")
cat("Using COMPLETE GENE DATA (filtered for all genes present)\n")
cat("============================================================\n")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ============================================================
# PSEUDOLOG TRANSFORMATION FUNCTIONS
# ============================================================

# For individual tree level - stretched near zero
pseudolog10_individual <- function(x) {
  asinh(x / 0.2) / log(10)
}

inv_pseudolog10_individual <- function(x) {
  0.2 * sinh(x * log(10))
}

# For species level - original scaling
pseudolog10 <- function(x) {
  asinh(x / 2) / log(10)
}

inv_pseudolog10 <- function(x) {
  2 * sinh(x * log(10))
}

# ============================================================
# DATA PREPARATION (from ratio analysis code)
# ============================================================

# Load data (ensure paths are correct for your system)
ymf2023 <- read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv")
ymf2021 <- read.csv('data/processed/integrated/merged_tree_dataset_final.csv')

# Species mapping
species_mapping <- c(
  "ACRU" = "Acer rubrum", "ACSA" = "Acer saccharum", 
  "BEAL" = "Betula alleghaniensis", "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera", "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana", "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra", "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata", "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina", "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina", "SAAL" = "Sassafras albidum"
)

# Area-weighted function
area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1*r2 + r2^2
  
  gene_outer + (gene_inner - gene_outer) * (S / (3 * R^2))
}

# Prepare long format data for all genes
prepare_long_all_genes <- function(df) {
  df %>%
    dplyr::select(tree_id, starts_with("ddpcr_")) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    pivot_longer(
      cols = starts_with("ddpcr_"),
      names_to = "measurement_type",
      values_to = "gene_copies"
    ) %>%
    filter(!is.na(gene_copies)) %>%
    separate(
      measurement_type,
      into = c("method", "gene", "part1", "part2", "part3"),
      sep = "_", extra = "merge", fill = "right"
    ) %>%
    mutate(
      is_probe = (part1 == "probe"),
      location = if_else(is_probe, part2, part1),
      stringency = if_else(is_probe, part3, part2),
      location = str_remove(location, "probe"),
      location = case_when(
        location %in% c("Inner","inner") ~ "Inner",
        location %in% c("Outer","outer") ~ "Outer",
        TRUE ~ location
      ),
      gene = case_when(
        gene == "mcra" ~ "mcrA",
        gene == "mmox" ~ "mmoX",
        gene == "pmoa" ~ "pmoA",
        TRUE ~ gene
      ),
      sample_type = case_when(
        location == "Inner" ~ "Heartwood",
        location == "Outer" ~ "Sapwood",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(stringency == "loose", 
           !is.na(sample_type),
           sample_type %in% c("Heartwood", "Sapwood"))
}

# Create area-weighted gene dataset
long_all <- prepare_long_all_genes(ymf2021)

tree_genes_weighted <- long_all %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
  filter(is.finite(dbh)) %>%
  group_by(tree_id, species_id, dbh, gene) %>%
  summarise(
    gene_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
    gene_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(gene_inner), is.finite(gene_outer)) %>%
  mutate(
    gene_area_weighted = mapply(area_weighted_gene, gene_inner, gene_outer, dbh),
    species = species_mapping[species_id]
  ) %>%
  dplyr::select(tree_id, species_id, species, gene, gene_area_weighted) %>%
  pivot_wider(names_from = gene, values_from = gene_area_weighted)

# Add methanotroph total
tree_genes_weighted <- tree_genes_weighted %>%
  mutate(
    methanotroph_total = case_when(
      is.na(pmoA) & !is.na(mmoX) ~ mmoX,
      !is.na(pmoA) & is.na(mmoX) ~ pmoA,
      !is.na(pmoA) & !is.na(mmoX) ~ pmoA + mmoX,
      TRUE ~ NA_real_
    )
  )

# ============================================================
# CRITICAL FILTER: KEEP ONLY TREES WITH ALL GENES PRESENT
# ============================================================

cat("\n=== FILTERING FOR COMPLETE GENE DATA ===\n")
cat("Before filtering:\n")
cat("  Total trees:", nrow(tree_genes_weighted), "\n")

tree_genes_complete <- tree_genes_weighted %>%
  filter(!is.na(mcrA), !is.na(pmoA), !is.na(mmoX), !is.na(methanotroph_total))

cat("After filtering for complete genes:\n")
cat("  Trees with all genes:", nrow(tree_genes_complete), "\n")
cat("  Species represented:", n_distinct(tree_genes_complete$species), "\n")

# Check which species have enough trees
species_counts <- tree_genes_complete %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

cat("\nTrees per species (complete gene data):\n")
print(species_counts)

# ============================================================
# FLUX DATA AND MERGING WITH COMPLETE GENE DATA
# ============================================================

# Flux data
flux_all <- bind_rows(
  ymf2023 %>%
    dplyr::select(species_id = Species.Code, CH4_flux = CH4_best.flux) %>%
    mutate(species = species_mapping[species_id]),
  ymf2021 %>%
    dplyr::select(species_id, CH4_flux = CH4_best.flux_125cm) %>%
    mutate(species = species_mapping[species_id])
) %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux), !is.na(species))

flux_by_species <- flux_all %>%
  group_by(species, species_id) %>%
  summarise(
    n_flux = n(),
    median_flux = median(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# Tree-level data with flux - USING COMPLETE GENE DATA ONLY
tree_level_complete <- tree_genes_complete %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
    by = "tree_id"
  ) %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux)) %>%
  mutate(
    log_tree_mcra = log10(mcrA + 1),
    log_tree_mmox = log10(mmoX + 1),
    log_pmoa = log10(pmoA + 1),
    log_methanotroph = log10(methanotroph_total + 1),
    ratio_mcra_methanotroph = (mcrA + 1) / (methanotroph_total + 1),
    log_ratio = log10(ratio_mcra_methanotroph)
  )

cat("\n=== FINAL TREE-LEVEL DATA ===\n")
cat("Trees with complete genes AND flux:", nrow(tree_level_complete), "\n")

# ============================================================
# SPECIES-LEVEL AGGREGATIONS - ALL FROM SAME FILTERED DATASET
# ============================================================

# Calculate species-level summaries for each gene from the SAME filtered dataset
analysis_mcra <- tree_level_complete %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_mcra = median(mcrA, na.rm = TRUE),
    q25_mcra = quantile(mcrA, 0.25, na.rm = TRUE),
    q75_mcra = quantile(mcrA, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

analysis_pmoa <- tree_level_complete %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_pmoa = median(pmoA, na.rm = TRUE),
    q25_pmoa = quantile(pmoA, 0.25, na.rm = TRUE),
    q75_pmoa = quantile(pmoA, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

analysis_mmox <- tree_level_complete %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_mmox = median(mmoX, na.rm = TRUE),
    q25_mmox = quantile(mmoX, 0.25, na.rm = TRUE),
    q75_mmox = quantile(mmoX, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

analysis_methanotroph <- tree_level_complete %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_methanotroph = median(methanotroph_total, na.rm = TRUE),
    q25_methanotroph = quantile(methanotroph_total, 0.25, na.rm = TRUE),
    q75_methanotroph = quantile(methanotroph_total, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

analysis_ratio <- tree_level_complete %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_ratio = median(ratio_mcra_methanotroph, na.rm = TRUE),
    median_log_ratio = median(log_ratio, na.rm = TRUE),
    q25_ratio = quantile(ratio_mcra_methanotroph, 0.25, na.rm = TRUE),
    q75_ratio = quantile(ratio_mcra_methanotroph, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

cat("\n=== SPECIES-LEVEL SAMPLE SIZES ===\n")
cat("Species in mcrA analysis:", nrow(analysis_mcra), "\n")
cat("Species in pmoA analysis:", nrow(analysis_pmoa), "\n")
cat("Species in mmoX analysis:", nrow(analysis_mmox), "\n")
cat("Species in methanotroph analysis:", nrow(analysis_methanotroph), "\n")
cat("Species in ratio analysis:", nrow(analysis_ratio), "\n")

# Correlations for species-level
cor_area_mcra <- cor.test(log10(analysis_mcra$median_mcra + 1), 
                          analysis_mcra$median_flux)
cor_area_pmoa <- cor.test(log10(analysis_pmoa$median_pmoa + 1), 
                          analysis_pmoa$median_flux)
cor_area_mmox <- cor.test(log10(analysis_mmox$median_mmox + 1), 
                          analysis_mmox$median_flux)
cor_area_methanotroph <- cor.test(log10(analysis_methanotroph$median_methanotroph + 1),
                                  analysis_methanotroph$median_flux)
pearson_ratio <- cor.test(analysis_ratio$median_log_ratio,
                          analysis_ratio$median_flux)

# ============================================================
# SHARED PUBLICATION INFRASTRUCTURE
# ============================================================

library(grid)

# Consistent species color palette (shared across S8 left-column and S11)
species_in_data <- sort(unique(tree_level_complete$species))
species_palette <- setNames(
  scales::viridis_pal()(length(species_in_data)),
  species_in_data
)

# Shared theme for all scatter/bar panels
theme_pub_gene <- theme_classic(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.tag = element_text(size = 11, face = "bold")
  )

# Y-axis labels with units
y_lab_individual <- expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))
y_lab_species    <- expression(Median~CH[4]~flux~(nmol~m^{-2}~s^{-1}))

# ============================================================
# LEFT COLUMN: INDIVIDUAL TREE LEVEL - PSEUDOLOG FLUX
# ============================================================

# Panel A1: Individual - mcrA
tree_lm_mcra <- lm(CH4_flux ~ species + log_tree_mcra, data = tree_level_complete)
tree_mcra_r2 <- summary(tree_lm_mcra)$r.squared
tree_mcra_p <- summary(tree_lm_mcra)$coefficients["log_tree_mcra", "Pr(>|t|)"]

p_tree_mcra <- ggplot(tree_level_complete,
                      aes(x = log_tree_mcra, y = pseudolog10_individual(CH4_flux), color = species)) +
  geom_point(size = 1.95, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C",
              fill = "#FADBD8", alpha = 0.2, linewidth = 1) +
  geom_hline(yintercept = pseudolog10_individual(0), linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", tree_mcra_r2, tree_mcra_p),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  scale_y_continuous(
    breaks = pseudolog10_individual(c(0, 0.1, 1)),
    labels = c("0", "0.1", "1")
  ) +
  scale_color_manual(values = species_palette, name = "Species") +
  labs(tag = "(a)", x = expression(log[10]~mcrA), y = y_lab_individual) +
  theme_pub_gene +
  theme(legend.position = "none")

# Panel A2: Individual - pmoA
tree_lm_pmoa <- lm(CH4_flux ~ species + log_pmoa, data = tree_level_complete)
tree_pmoa_r2 <- summary(tree_lm_pmoa)$r.squared
tree_pmoa_p <- summary(tree_lm_pmoa)$coefficients["log_pmoa", "Pr(>|t|)"]

p_tree_pmoa <- ggplot(tree_level_complete,
                      aes(x = log_pmoa, y = pseudolog10_individual(CH4_flux), color = species)) +
  geom_point(size = 1.95, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#3498DB",
              fill = "#D6EAF8", alpha = 0.2, linewidth = 1) +
  geom_hline(yintercept = pseudolog10_individual(0), linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", tree_pmoa_r2, tree_pmoa_p),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  scale_y_continuous(
    breaks = pseudolog10_individual(c(0, 0.1, 1)),
    labels = c("0", "0.1", "1")
  ) +
  scale_color_manual(values = species_palette, name = "Species") +
  labs(tag = "(c)", x = expression(log[10]~pmoA), y = "") +
  theme_pub_gene +
  theme(legend.position = "none")

# Panel A3: Individual - mmoX
tree_lm_mmox <- lm(CH4_flux ~ species + log_tree_mmox, data = tree_level_complete)
tree_mmox_r2 <- summary(tree_lm_mmox)$r.squared
tree_mmox_p <- summary(tree_lm_mmox)$coefficients["log_tree_mmox", "Pr(>|t|)"]

p_tree_mmox <- ggplot(tree_level_complete,
                      aes(x = log_tree_mmox, y = pseudolog10_individual(CH4_flux), color = species)) +
  geom_point(size = 1.95, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#5DADE2",
              fill = "#D6EAF8", alpha = 0.2, linewidth = 1) +
  geom_hline(yintercept = pseudolog10_individual(0), linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", tree_mmox_r2, tree_mmox_p),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  scale_y_continuous(
    breaks = pseudolog10_individual(c(0, 0.1, 1)),
    labels = c("0", "0.1", "1")
  ) +
  scale_color_manual(values = species_palette, name = "Species") +
  labs(tag = "(e)", x = expression(log[10]~mmoX), y = "") +
  theme_pub_gene +
  theme(legend.position = "none")

# Panel A4: Individual - Methanotrophs
tree_lm_methanotroph <- lm(CH4_flux ~ species + log_methanotroph, data = tree_level_complete)
tree_methanotroph_r2 <- summary(tree_lm_methanotroph)$r.squared
tree_methanotroph_p <- summary(tree_lm_methanotroph)$coefficients["log_methanotroph", "Pr(>|t|)"]

p_tree_methanotroph <- ggplot(tree_level_complete,
                              aes(x = log_methanotroph, y = pseudolog10_individual(CH4_flux), color = species)) +
  geom_point(size = 1.95, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#1F618D",
              fill = "#D6EAF8", alpha = 0.2, linewidth = 1) +
  geom_hline(yintercept = pseudolog10_individual(0), linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", tree_methanotroph_r2, tree_methanotroph_p),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  scale_y_continuous(
    breaks = pseudolog10_individual(c(0, 0.1, 1)),
    labels = c("0", "0.1", "1")
  ) +
  scale_color_manual(values = species_palette, name = "Species") +
  labs(tag = "(g)", x = expression(log[10]~(pmoA+mmoX)), y = "") +
  theme_pub_gene +
  theme(legend.position = "none")

# Panel A5: Individual - Ratio
tree_lm_ratio <- lm(CH4_flux ~ species + log_ratio, data = tree_level_complete)
tree_ratio_r2 <- summary(tree_lm_ratio)$r.squared
tree_ratio_p <- summary(tree_lm_ratio)$coefficients["log_ratio", "Pr(>|t|)"]

p_tree_ratio <- ggplot(tree_level_complete,
                       aes(x = log_ratio, y = pseudolog10_individual(CH4_flux), color = species)) +
  geom_point(size = 1.95, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#9B59B6",
              fill = "#EBDEF0", alpha = 0.2, linewidth = 1) +
  geom_hline(yintercept = pseudolog10_individual(0), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.5) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", tree_ratio_r2, tree_ratio_p),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  scale_y_continuous(
    breaks = pseudolog10_individual(c(0, 0.1, 1)),
    labels = c("0", "0.1", "1")
  ) +
  scale_color_manual(values = species_palette, name = "Species") +
  labs(tag = "(i)", x = expression(log[10]~ratio), y = "") +
  theme_pub_gene +
  theme(legend.position = "none")

# ============================================================
# RIGHT COLUMN: SPECIES LEVEL - LINEAR FLUX (WITH ERROR BARS)
# ============================================================

# Panel B1: Species - mcrA (with x error bars)
p_species_mcra <- ggplot(analysis_mcra,
                         aes(x = log10(median_mcra + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#E74C3C",
              fill = "#FADBD8", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#E74C3C") +
  geom_text_repel(aes(label = species), size = 2.6, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", cor_area_mcra$estimate^2, cor_area_mcra$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  labs(tag = "(b)", x = expression(log[10]~median~mcrA), y = y_lab_species) +
  theme_pub_gene

# Panel B2: Species - pmoA (with x error bars)
p_species_pmoa <- ggplot(analysis_pmoa,
                         aes(x = log10(median_pmoa + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#3498DB",
              fill = "#D6EAF8", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#3498DB") +
  geom_text_repel(aes(label = species), size = 2.6, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", cor_area_pmoa$estimate^2, cor_area_pmoa$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  labs(tag = "(d)", x = expression(log[10]~median~pmoA), y = "") +
  theme_pub_gene

# Panel B3: Species - mmoX (with x error bars)
p_species_mmox <- ggplot(analysis_mmox,
                         aes(x = log10(median_mmox + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#5DADE2",
              fill = "#D6EAF8", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#5DADE2") +
  geom_text_repel(aes(label = species), size = 2.6, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", cor_area_mmox$estimate^2, cor_area_mmox$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  labs(tag = "(f)", x = expression(log[10]~median~mmoX), y = "") +
  theme_pub_gene

# Panel B4: Species - Methanotrophs (with x error bars)
p_species_methanotroph <- ggplot(analysis_methanotroph,
                                 aes(x = log10(median_methanotroph + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#1F618D",
              fill = "#D6EAF8", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#1F618D") +
  geom_text_repel(aes(label = species), size = 2.6, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", cor_area_methanotroph$estimate^2, cor_area_methanotroph$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  labs(tag = "(h)", x = expression(log[10]~median~(pmoA+mmoX)), y = "") +
  theme_pub_gene

# Panel B5: Species - Ratio
p_species_ratio <- ggplot(analysis_ratio,
                          aes(x = median_log_ratio, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#9B59B6",
              fill = "#EBDEF0", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#9B59B6") +
  geom_text_repel(aes(label = species), size = 2.6, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.5) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f", pearson_ratio$estimate^2, pearson_ratio$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.25, fill = "white", alpha = 0.9) +
  labs(tag = "(j)", x = expression(log[10]~ratio), y = "") +
  theme_pub_gene

# ============================================================
# BOTTOM ROW: MODEL COMPARISONS (A6 and B6)
# ============================================================

# Panel A6: Individual level comparison
tree_comparison_data <- data.frame(
  Model = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  R2 = c(tree_mcra_r2, tree_pmoa_r2, tree_mmox_r2, 
         tree_methanotroph_r2, tree_ratio_r2),
  P_value = c(tree_mcra_p, tree_pmoa_p, tree_mmox_p, 
              tree_methanotroph_p, tree_ratio_p)
) %>%
  mutate(
    Significant = P_value < 0.05,
    Model = factor(Model, levels = Model[order(R2)])
  )

# Panel B6: Species level comparison
species_comparison_data <- data.frame(
  Model = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  R2 = c(cor_area_mcra$estimate^2,
         cor_area_pmoa$estimate^2,
         cor_area_mmox$estimate^2,
         cor_area_methanotroph$estimate^2,
         pearson_ratio$estimate^2),
  P_value = c(cor_area_mcra$p.value,
              cor_area_pmoa$p.value,
              cor_area_mmox$p.value,
              cor_area_methanotroph$p.value,
              pearson_ratio$p.value)
) %>%
  mutate(
    Significant = P_value < 0.05,
    Model = factor(Model, levels = Model[order(R2)])
  )

# Calculate shared y-axis limit
max_r2 <- max(c(tree_comparison_data$R2, species_comparison_data$R2))
y_limit <- max_r2 * 1.2

p_tree_comparison <- ggplot(tree_comparison_data,
                            aes(x = Model, y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            vjust = -0.3, size = 3.25) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#27AE60"),
                    labels = c("NS", "p < 0.05"),
                    name = "") +
  ylim(0, y_limit) +
  labs(tag = "(k)", x = "", y = expression(R^2)) +
  theme_pub_gene +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10))

p_species_comparison <- ggplot(species_comparison_data,
                               aes(x = Model, y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            vjust = -0.3, size = 3.25) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#27AE60"),
                    labels = c("NS", "p < 0.05"),
                    name = "") +
  ylim(0, y_limit) +
  labs(tag = "(l)", x = "", y = "") +
  theme_pub_gene +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10))

# ============================================================
# COMBINE LAYOUT (6x2 with column headers)
# ============================================================

# Column headers
header_individual <- wrap_elements(full =
  textGrob("Individual Tree Level", gp = gpar(fontface = "bold", fontsize = 13)))
header_species <- wrap_elements(full =
  textGrob("Species Level", gp = gpar(fontface = "bold", fontsize = 13)))

# Create a legend-only panel as a ggplot (cowplot::get_legend doesn't scale in patchwork)
legend_panel <- ggplot(tree_level_complete,
    aes(x = log_tree_mcra, y = CH4_flux, color = species)) +
  geom_point(alpha = 0) +
  scale_color_manual(values = species_palette, name = "Species") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(size = 3, alpha = 1))) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(face = "italic", size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0))

# Build left and right columns with headers
left_column <- header_individual /
  p_tree_mcra / p_tree_pmoa / p_tree_mmox /
  p_tree_methanotroph / p_tree_ratio / p_tree_comparison +
  plot_layout(ncol = 1, heights = c(0.06, 1, 1, 1, 1, 1, 0.85))

right_column <- header_species /
  p_species_mcra / p_species_pmoa / p_species_mmox /
  p_species_methanotroph / p_species_ratio / p_species_comparison +
  plot_layout(ncol = 1, heights = c(0.06, 1, 1, 1, 1, 1, 0.85))

# Combine: two columns side-by-side, legend row at bottom
main_grid <- left_column | right_column
combined_layout <- main_grid / legend_panel +
  plot_layout(heights = c(1, 0.08))

ggsave("outputs/figures/supplementary/figS10_scale_dependent_genes.png", combined_layout,
       width = 12, height = 15, dpi = 300)

# ============================================================
# INDIVIDUAL-ONLY LAYOUT (2 columns x 3 rows)
# ============================================================

individual_layout <- (p_tree_mcra | p_tree_pmoa) /
  (p_tree_mmox | p_tree_methanotroph) /
  (p_tree_ratio | p_tree_comparison)

# ggsave("outputs/figures/Figure_Individual_Only_Complete_Genes.pdf", individual_layout,
#        width = 8, height = 10)
# ggsave("outputs/figures/Figure_Individual_Only_Complete_Genes.png", individual_layout,
#        width = 8, height = 10, dpi = 300)

# ============================================================
# SPECIES-ONLY LAYOUT (2 columns x 3 rows)
# ============================================================

species_layout <- (p_species_mcra | p_species_pmoa) /
  (p_species_mmox | p_species_methanotroph) /
  (p_species_ratio | p_species_comparison)

# ggsave("outputs/figures/Figure_Species_Only_Complete_Genes.pdf", species_layout,
#        width = 8, height = 10)
# ggsave("outputs/figures/Figure_Species_Only_Complete_Genes.png", species_layout,
#        width = 8, height = 10, dpi = 300)

# ============================================================
# SPECIES-ONLY 2x2 LAYOUT (mcrA, methanotrophs, ratio, comparison)
# ============================================================

# Panel 1: Species - mcrA (custom red color)
p_species_mcra_2x2 <- ggplot(analysis_mcra,
                             aes(x = log10(median_mcra + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#C03221",
              fill = "#F5C9C3", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#C03221") +
  geom_text_repel(aes(label = species), size = 2.6, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           cor_area_mcra$estimate^2,
                           cor_area_mcra$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.25,
           fill = "white", alpha = 0.9) +
  labs(x = expression("log"[10]*" median mcrA"),
       y = expression("Median CH"[4]*" flux")) +
  theme_classic(base_size = 11.7) +
  theme(
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Panel 2: Species - Methanotrophs (custom blue color)
p_species_methanotroph_2x2 <- ggplot(analysis_methanotroph,
                                     aes(x = log10(median_methanotroph + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#4A6FA5",
              fill = "#C5D5E8", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#4A6FA5") +
  geom_text_repel(aes(label = species), size = 2.6, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           cor_area_methanotroph$estimate^2,
                           cor_area_methanotroph$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.25,
           fill = "white", alpha = 0.9) +
  labs(x = expression("log"[10]*" median (pmoA+mmoX)"),
       y = "") +
  theme_classic(base_size = 11.7) +
  theme(
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Panel 3: Species - Ratio (custom purple color)
p_species_ratio_2x2 <- ggplot(analysis_ratio,
                              aes(x = median_log_ratio, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "#6B5B95",
              fill = "#D7D2E0", alpha = 0.2, linewidth = 1) +
  geom_point(size = 3.25, alpha = 0.85, color = "#6B5B95") +
  geom_text_repel(aes(label = species), size = 2.6, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.5) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           pearson_ratio$estimate^2,
                           pearson_ratio$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.25,
           fill = "white", alpha = 0.9) +
  labs(x = expression("log"[10]*" ratio"),
       y = expression("Median CH"[4]*" flux")) +
  theme_classic(base_size = 11.7) +
  theme(
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Panel 4: Model Comparison (custom green for significance)
p_species_comparison_2x2 <- ggplot(species_comparison_data,
                                   aes(x = Model, y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            vjust = -0.3, size = 3.25) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#285238"),
                    labels = c("NS", "p < 0.05"),
                    name = "") +
  ylim(0, y_limit) +
  labs(x = "", y = "") +
  theme_classic(base_size = 11.7) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10.4))

# Create 2x2 layout
species_2x2_layout <- (p_species_mcra_2x2 | p_species_methanotroph_2x2) /
  (p_species_ratio_2x2 | p_species_comparison_2x2)

# ggsave("outputs/figures/Figure_Species_2x2_Complete_Genes.pdf", species_2x2_layout,
#        width = 8, height = 8)
# ggsave("outputs/figures/Figure_Species_2x2_Complete_Genes.png", species_2x2_layout,
#        width = 8, height = 8, dpi = 300)

cat("\n============================================================\n")
cat("ALL FIGURES CREATED (COMPLETE GENE DATA)\n")
cat("============================================================\n")
cat("\nFiles saved:\n")
cat("  1. Combined layout (side-by-side):\n")
cat("     - Figure_Scale_Dependent_Complete_Genes.pdf\n")
cat("     - Figure_Scale_Dependent_Complete_Genes.png\n")
cat("  2. Individual tree level only (2x3):\n")
cat("     - Figure_Individual_Only_Complete_Genes.pdf\n")
cat("     - Figure_Individual_Only_Complete_Genes.png\n")
cat("  3. Species level only (2x3):\n")
cat("     - Figure_Species_Only_Complete_Genes.pdf\n")
cat("     - Figure_Species_Only_Complete_Genes.png\n")
cat("  4. Species level 2x2 (key metrics):\n")
cat("     - Figure_Species_2x2_Complete_Genes.pdf\n")
cat("     - Figure_Species_2x2_Complete_Genes.png\n")

cat("\n### SUMMARY TABLE ###\n")
summary_table <- rbind(
  tree_comparison_data %>% mutate(Scale = "Individual"),
  species_comparison_data %>% mutate(Scale = "Species")
) %>%
  dplyr::select(Scale, Model, R2, P_value, Significant) %>%
  arrange(Scale, desc(R2))

print(summary_table)

write.csv(summary_table, "outputs/tables/Scale_Dependent_Complete_Genes_Summary.csv", row.names = FALSE)

cat("\n### DATA FILTERING SUMMARY ###\n")
cat("Trees with complete gene data:", nrow(tree_level_complete), "\n")
cat("Species analyzed:", nrow(analysis_ratio), "\n")
cat("All analyses use the SAME filtered dataset\n")

cat("\n============================================================\n")

















# ============================================================
# ADDITIONAL ANALYSIS: mcrA vs. (pmoA+mmoX) 
# At both individual tree and species levels
# ============================================================

cat("\n============================================================\n")
cat("CREATING mcrA vs METHANOTROPH PLOTS\n")
cat("============================================================\n")

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
  geom_text_repel(aes(label = species), size = 3, fontface = "italic",
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

# Save combined figure
# ggsave("outputs/figures/Figure_mcrA_vs_Methanotroph.pdf", mcra_vs_methanotroph_layout,
#        width = 12, height = 5)
ggsave("outputs/figures/supplementary/figS13_mcra_vs_methanotroph.png", mcra_vs_methanotroph_layout,
       width = 13, height = 5.5, dpi = 300)

# Save individual plots
# ggsave("outputs/figures/Figure_mcrA_vs_Methanotroph_Individual.pdf",
#        p_tree_mcra_vs_methanotroph,
#        width = 6, height = 5)
# ggsave("outputs/figures/Figure_mcrA_vs_Methanotroph_Species.pdf",
#        p_species_mcra_vs_methanotroph,
#        width = 6, height = 5)

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

cat("\nFigures saved:\n")
cat("  - Figure_mcrA_vs_Methanotroph.pdf (combined)\n")
cat("  - Figure_mcrA_vs_Methanotroph.png (combined)\n")
cat("  - Figure_mcrA_vs_Methanotroph_Individual.pdf\n")
cat("  - Figure_mcrA_vs_Methanotroph_Species.pdf\n")

cat("\n============================================================\n")
