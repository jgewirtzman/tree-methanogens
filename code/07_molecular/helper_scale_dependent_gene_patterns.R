# ==============================================================================
# helper_scale_dependent_gene_patterns.R -- species-level gene-flux summaries
# ------------------------------------------------------------------------------
# SOURCED, never run alone (from the repository root). Builds, from the canonical
# tables (code/lib/species_gene_data.R):
#   tree_level_complete, flux_all, ...      (species_gene_data.R)
#   analysis_mcra, _pmoa, _mmox, _methanotroph, _ratio
#                                           per species (>= 5 trees with every gene
#                                           and >= 5 fluxes): median and quartiles of
#                                           the gene, median flux
#   cor_area_mcra, _pmoa, _mmox, _methanotroph, pearson_ratio
#                                           species-level Pearson tests (log10 gene)
#   tree_comparison_data, species_comparison_data, y_limit
#                                           R2 and p of the five gene models at tree
#                                           (flux ~ species + gene) and species level
#   species_palette, theme_pub_gene         shared plot styling
# Sourced by fig04, fig08, figS21 and 09_tables_stats/14-20.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

source("code/lib/species_gene_data.R")

cat("\n=== FILTERING FOR COMPLETE GENE DATA ===\n")
cat("Before filtering:\n")
cat("  Total trees:", nrow(tree_genes_weighted), "\n")
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
# SHARED PLOT STYLING
# ============================================================
library(grid)

# species colours (Fig S24)
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

# ============================================================
# GENE MODELS: tree level (flux ~ species + gene) and species level
# ============================================================
tree_lm_mcra <- lm(CH4_flux ~ species + log_tree_mcra, data = tree_level_complete)
tree_mcra_r2 <- summary(tree_lm_mcra)$r.squared
tree_mcra_p <- summary(tree_lm_mcra)$coefficients["log_tree_mcra", "Pr(>|t|)"]
tree_lm_pmoa <- lm(CH4_flux ~ species + log_pmoa, data = tree_level_complete)
tree_pmoa_r2 <- summary(tree_lm_pmoa)$r.squared
tree_pmoa_p <- summary(tree_lm_pmoa)$coefficients["log_pmoa", "Pr(>|t|)"]
tree_lm_mmox <- lm(CH4_flux ~ species + log_tree_mmox, data = tree_level_complete)
tree_mmox_r2 <- summary(tree_lm_mmox)$r.squared
tree_mmox_p <- summary(tree_lm_mmox)$coefficients["log_tree_mmox", "Pr(>|t|)"]
tree_lm_methanotroph <- lm(CH4_flux ~ species + log_methanotroph, data = tree_level_complete)
tree_methanotroph_r2 <- summary(tree_lm_methanotroph)$r.squared
tree_methanotroph_p <- summary(tree_lm_methanotroph)$coefficients["log_methanotroph", "Pr(>|t|)"]
tree_lm_ratio <- lm(CH4_flux ~ species + log_ratio, data = tree_level_complete)
tree_ratio_r2 <- summary(tree_lm_ratio)$r.squared
tree_ratio_p <- summary(tree_lm_ratio)$coefficients["log_ratio", "Pr(>|t|)"]

# R2 and p of each gene model, tree level
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

# species level
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
