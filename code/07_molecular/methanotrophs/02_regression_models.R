# ==============================================================================
# Methanotroph Regression Models
# ==============================================================================
# Purpose: Multiple regression models linking pmoA/mmoxY to CH4 flux.
#
# Pipeline stage: 3 — Analysis
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#   - flux data (from data/processed/flux/)
# ==============================================================================

cat("\n============================================================\n")
cat("CREATING 3x2 SIDE-BY-SIDE LAYOUT FIGURE\n")
cat("Using ratio analysis data preparation\n")
cat("============================================================\n")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ============================================================
# PSEUDOLOG TRANSFORMATION FUNCTIONS
# ============================================================

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
ymf2023 <- read.csv("../../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv")
ymf2021 <- read.csv('../../../data/processed/integrated/merged_tree_dataset_final.csv')

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

# Tree-level data with flux
tree_level_ratio <- tree_genes_weighted %>%
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

# Species-level aggregations
analysis_mcra <- tree_genes_weighted %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm), 
            by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.na(mcrA)) %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_mcra = median(mcrA, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

analysis_methanotroph <- tree_genes_weighted %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm), 
            by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.na(methanotroph_total)) %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_methanotroph = median(methanotroph_total, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

analysis_ratio <- tree_level_ratio %>%
  filter(!is.na(mcrA), !is.na(methanotroph_total), !is.na(log_ratio)) %>%
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

# Correlations for species-level
cor_area_mcra <- cor.test(log10(analysis_mcra$median_mcra + 1), 
                          analysis_mcra$median_flux)
cor_area_methanotroph <- cor.test(log10(analysis_methanotroph$median_methanotroph + 1),
                                  analysis_methanotroph$median_flux)
pearson_ratio <- cor.test(analysis_ratio$median_log_ratio, 
                          analysis_ratio$median_flux)

# ============================================================
# LEFT COLUMN: INDIVIDUAL TREE LEVEL - PSEUDOLOG FLUX
# ============================================================

# Panel A1: Individual - mcrA
tree_mcra_data <- tree_level_ratio %>% filter(!is.na(log_tree_mcra))
tree_lm_mcra <- lm(CH4_flux ~ species + log_tree_mcra, data = tree_mcra_data)
tree_mcra_r2 <- summary(tree_lm_mcra)$r.squared
tree_mcra_p <- summary(tree_lm_mcra)$coefficients["log_tree_mcra", "Pr(>|t|)"]

p_tree_mcra <- ggplot(tree_mcra_data, 
                      aes(x = log_tree_mcra, y = pseudolog10(CH4_flux), color = species)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", 
              fill = "gray80", alpha = 0.3, linewidth = 0.8) +
  geom_hline(yintercept = pseudolog10(0), linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           tree_mcra_r2,
                           tree_mcra_p),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  scale_y_continuous(
    breaks = pseudolog10(c(-0.01, 0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)),
    labels = c("-0.01", "0", "0.02", "0.05", "0.1", "0.2", "0.5", "1", "2")
  ) +
  labs(title = "A1) mcrA",
       x = expression("log"[10]*" mcrA"),
       y = expression("CH"[4]*" flux (pseudolog)")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

# Panel A2: Individual - pmoA
tree_pmoa_data <- tree_level_ratio %>% filter(!is.na(log_pmoa))
tree_lm_pmoa <- lm(CH4_flux ~ species + log_pmoa, data = tree_pmoa_data)
tree_pmoa_r2 <- summary(tree_lm_pmoa)$r.squared
tree_pmoa_p <- summary(tree_lm_pmoa)$coefficients["log_pmoa", "Pr(>|t|)"]

p_tree_pmoa <- ggplot(tree_pmoa_data,
                      aes(x = log_pmoa, y = pseudolog10(CH4_flux), color = species)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "gray80", alpha = 0.3, linewidth = 0.8) +
  geom_hline(yintercept = pseudolog10(0), linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           tree_pmoa_r2,
                           tree_pmoa_p),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  scale_y_continuous(
    breaks = pseudolog10(c(-0.01, 0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)),
    labels = c("-0.01", "0", "0.02", "0.05", "0.1", "0.2", "0.5", "1", "2")
  ) +
  labs(title = "A2) pmoA",
       x = expression("log"[10]*" pmoA"),
       y = "") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

# Panel A3: Individual - mmoX
tree_mmox_data <- tree_level_ratio %>% filter(!is.na(log_tree_mmox))
tree_lm_mmox <- lm(CH4_flux ~ species + log_tree_mmox, data = tree_mmox_data)
tree_mmox_r2 <- summary(tree_lm_mmox)$r.squared
tree_mmox_p <- summary(tree_lm_mmox)$coefficients["log_tree_mmox", "Pr(>|t|)"]

p_tree_mmox <- ggplot(tree_mmox_data,
                      aes(x = log_tree_mmox, y = pseudolog10(CH4_flux), color = species)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "gray80", alpha = 0.3, linewidth = 0.8) +
  geom_hline(yintercept = pseudolog10(0), linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           tree_mmox_r2,
                           tree_mmox_p),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  scale_y_continuous(
    breaks = pseudolog10(c(-0.01, 0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)),
    labels = c("-0.01", "0", "0.02", "0.05", "0.1", "0.2", "0.5", "1", "2")
  ) +
  labs(title = "A3) mmoX",
       x = expression("log"[10]*" mmoX"),
       y = expression("CH"[4]*" flux (pseudolog)")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

# Panel A4: Individual - Methanotrophs
tree_methanotroph_data <- tree_level_ratio %>% filter(!is.na(log_methanotroph))
tree_lm_methanotroph <- lm(CH4_flux ~ species + log_methanotroph, data = tree_methanotroph_data)
tree_methanotroph_r2 <- summary(tree_lm_methanotroph)$r.squared
tree_methanotroph_p <- summary(tree_lm_methanotroph)$coefficients["log_methanotroph", "Pr(>|t|)"]

p_tree_methanotroph <- ggplot(tree_methanotroph_data,
                              aes(x = log_methanotroph, y = pseudolog10(CH4_flux), color = species)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "gray80", alpha = 0.3, linewidth = 0.8) +
  geom_hline(yintercept = pseudolog10(0), linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           tree_methanotroph_r2,
                           tree_methanotroph_p),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  scale_y_continuous(
    breaks = pseudolog10(c(-0.01, 0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)),
    labels = c("-0.01", "0", "0.02", "0.05", "0.1", "0.2", "0.5", "1", "2")
  ) +
  labs(title = "A4) pmoA+mmoX",
       x = expression("log"[10]*" (pmoA+mmoX)"),
       y = "") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

# Panel A5: Individual - Ratio
tree_ratio_data <- tree_level_ratio %>% filter(!is.na(log_ratio))
tree_lm_ratio <- lm(CH4_flux ~ species + log_ratio, data = tree_ratio_data)
tree_ratio_r2 <- summary(tree_lm_ratio)$r.squared
tree_ratio_p <- summary(tree_lm_ratio)$coefficients["log_ratio", "Pr(>|t|)"]

p_tree_ratio <- ggplot(tree_ratio_data,
                       aes(x = log_ratio, y = pseudolog10(CH4_flux), color = species)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "gray80", alpha = 0.3, linewidth = 0.8) +
  geom_hline(yintercept = pseudolog10(0), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.5) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           tree_ratio_r2,
                           tree_ratio_p),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  scale_y_continuous(
    breaks = pseudolog10(c(-0.01, 0, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2)),
    labels = c("-0.01", "0", "0.02", "0.05", "0.1", "0.2", "0.5", "1", "2")
  ) +
  labs(title = "A5) Ratio",
       x = expression("log"[10]*" ratio"),
       y = expression("CH"[4]*" flux (pseudolog)")) +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

# ============================================================
# RIGHT COLUMN: SPECIES LEVEL - LINEAR FLUX (WITH ERROR BARS)
# ============================================================

# Panel B1: Species - mcrA (with x error bars)
analysis_mcra_full <- tree_genes_weighted %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm), 
            by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.na(mcrA)) %>%
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

cor_area_mcra <- cor.test(log10(analysis_mcra_full$median_mcra + 1), 
                          analysis_mcra_full$median_flux)

p_species_mcra <- ggplot(analysis_mcra_full,
                         aes(x = log10(median_mcra + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred",
              fill = "pink", alpha = 0.2, linewidth = 0.8) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.4, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(q25_mcra + 1), xmax = log10(q75_mcra + 1)),
                 height = 0.003, alpha = 0.4, color = "gray40") +
  geom_point(size = 2.5, alpha = 0.85, color = "darkred") +
  geom_text_repel(aes(label = species), size = 2, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           cor_area_mcra$estimate^2,
                           cor_area_mcra$p.value),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  labs(title = "B1) mcrA",
       x = expression("log"[10]*" median mcrA"),
       y = expression("Median CH"[4]*" flux")) +
  theme_classic(base_size = 9)

# Panel B2: Species - pmoA (with x error bars)
species_pmoa <- tree_genes_weighted %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm), by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.na(pmoA)) %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_pmoa = median(pmoA, na.rm = TRUE),
    q25_pmoa = quantile(pmoA, 0.25, na.rm = TRUE),
    q75_pmoa = quantile(pmoA, 0.75, na.rm = TRUE),
    median_flux = median(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_trees >= 5)

species_cor_pmoa <- cor.test(log10(species_pmoa$median_pmoa + 1),
                             species_pmoa$median_flux)

p_species_pmoa <- ggplot(species_pmoa,
                         aes(x = log10(median_pmoa + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue",
              fill = "lightblue", alpha = 0.2, linewidth = 0.8) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.4, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(q25_pmoa + 1), xmax = log10(q75_pmoa + 1)),
                 height = 0.003, alpha = 0.4, color = "gray40") +
  geom_point(size = 2.5, alpha = 0.85, color = "darkblue") +
  geom_text_repel(aes(label = species), size = 2, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           species_cor_pmoa$estimate^2,
                           species_cor_pmoa$p.value),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  labs(title = "B2) pmoA",
       x = expression("log"[10]*" median pmoA"),
       y = "") +
  theme_classic(base_size = 9)

# Panel B3: Species - mmoX (with x error bars)
species_mmox <- tree_genes_weighted %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm), by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.na(mmoX)) %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_mmox = median(mmoX, na.rm = TRUE),
    q25_mmox = quantile(mmoX, 0.25, na.rm = TRUE),
    q75_mmox = quantile(mmoX, 0.75, na.rm = TRUE),
    median_flux = median(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_trees >= 5)

species_cor_mmox <- cor.test(log10(species_mmox$median_mmox + 1),
                             species_mmox$median_flux)

p_species_mmox <- ggplot(species_mmox,
                         aes(x = log10(median_mmox + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkorange",
              fill = "lightyellow", alpha = 0.2, linewidth = 0.8) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.4, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(q25_mmox + 1), xmax = log10(q75_mmox + 1)),
                 height = 0.003, alpha = 0.4, color = "gray40") +
  geom_point(size = 2.5, alpha = 0.85, color = "darkorange") +
  geom_text_repel(aes(label = species), size = 2, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           species_cor_mmox$estimate^2,
                           species_cor_mmox$p.value),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  labs(title = "B3) mmoX",
       x = expression("log"[10]*" median mmoX"),
       y = "") +
  theme_classic(base_size = 9)

# Panel B4: Species - Methanotrophs (with x error bars)
analysis_methanotroph_full <- tree_genes_weighted %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm), 
            by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.na(methanotroph_total)) %>%
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

species_cor_methanotroph <- cor.test(log10(analysis_methanotroph_full$median_methanotroph + 1),
                                     analysis_methanotroph_full$median_flux)

p_species_methanotroph <- ggplot(analysis_methanotroph_full,
                                 aes(x = log10(median_methanotroph + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen",
              fill = "lightgreen", alpha = 0.2, linewidth = 0.8) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.4, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(q25_methanotroph + 1), xmax = log10(q75_methanotroph + 1)),
                 height = 0.003, alpha = 0.4, color = "gray40") +
  geom_point(size = 2.5, alpha = 0.85, color = "darkgreen") +
  geom_text_repel(aes(label = species), size = 2, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f",
                           species_cor_methanotroph$estimate^2,
                           species_cor_methanotroph$p.value),
           hjust = 1.1, vjust = 1.1, size = 2.5) +
  labs(title = "B4) pmoA+mmoX",
       x = expression("log"[10]*" median (pmoA+mmoX)"),
       y = "") +
  theme_classic(base_size = 9)

# Panel B5: Species - Ratio (already has x and y error bars)
p_species_ratio <- ggplot(analysis_ratio,
                          aes(x = median_log_ratio, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "purple",
              fill = "plum", alpha = 0.2, linewidth = 0.8) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.4, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(q25_ratio), xmax = log10(q75_ratio)),
                 height = 0.003, alpha = 0.4, color = "gray40") +
  geom_point(size = 2.5, alpha = 0.85, color = "purple") +
  geom_text_repel(aes(label = species), size = 2, fontface = "italic",
                  box.padding = 0.2, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.5) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R²=%.3f\np=%.3f **",
                           pearson_ratio$estimate^2,
                           pearson_ratio$p.value),
           hjust = 1.1, vjust = 1.1, size = 2.5, color = "purple") +
  labs(title = "B5) Ratio",
       x = expression("log"[10]*" ratio"),
       y = "") +
  theme_classic(base_size = 9)

# Update Panel B6 to use the new correlation values
species_comparison_data <- data.frame(
  Model = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  R2 = c(cor_area_mcra$estimate^2,
         species_cor_pmoa$estimate^2,
         species_cor_mmox$estimate^2,
         species_cor_methanotroph$estimate^2,
         pearson_ratio$estimate^2),
  P_value = c(cor_area_mcra$p.value,
              species_cor_pmoa$p.value,
              species_cor_mmox$p.value,
              species_cor_methanotroph$p.value,
              pearson_ratio$p.value)
) %>%
  mutate(
    Significant = P_value < 0.05,
    Model = factor(Model, levels = Model[order(R2)])
  )

# Calculate shared y-axis limit
max_r2 <- max(c(tree_comparison_data$R2, species_comparison_data$R2))
y_limit <- max_r2 * 1.2

p_species_comparison <- ggplot(species_comparison_data,
                               aes(x = Model, y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkgreen"),
                    labels = c("NS", "p < 0.05"),
                    name = "") +
  ylim(0, y_limit) +
  labs(title = "B6) Model Comparison",
       subtitle = sprintf("Correlations (n=%d species)", nrow(analysis_ratio)),
       x = "", y = "") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 7))

# ============================================================
# BOTTOM ROW: MODEL COMPARISONS (A6 and B6)
# ============================================================

# Panel A6: Individual level comparison
tree_comparison_data <- data.frame(
  Model = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  R2 = c(summary(tree_lm_mcra)$r.squared,
         summary(tree_lm_pmoa)$r.squared,
         summary(tree_lm_mmox)$r.squared,
         summary(tree_lm_methanotroph)$r.squared,
         summary(tree_lm_ratio)$r.squared),
  P_value = c(summary(tree_lm_mcra)$coefficients["log_tree_mcra", "Pr(>|t|)"],
              summary(tree_lm_pmoa)$coefficients["log_pmoa", "Pr(>|t|)"],
              summary(tree_lm_mmox)$coefficients["log_tree_mmox", "Pr(>|t|)"],
              summary(tree_lm_methanotroph)$coefficients["log_methanotroph", "Pr(>|t|)"],
              summary(tree_lm_ratio)$coefficients["log_ratio", "Pr(>|t|)"])
) %>%
  mutate(
    Significant = P_value < 0.05,
    Model = factor(Model, levels = Model[order(R2)])
  )

# Panel B6: Species level comparison
species_comparison_data <- data.frame(
  Model = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  R2 = c(cor_area_mcra$estimate^2,
         species_cor_pmoa$estimate^2,
         species_cor_mmox$estimate^2,
         species_cor_methanotroph$estimate^2,
         pearson_ratio$estimate^2),
  P_value = c(cor_area_mcra$p.value,
              species_cor_pmoa$p.value,
              species_cor_mmox$p.value,
              species_cor_methanotroph$p.value,
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
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "steelblue"),
                    labels = c("NS", "p < 0.05"),
                    name = "") +
  ylim(0, y_limit) +
  labs(title = "A6) Model Comparison",
       subtitle = sprintf("flux ~ species + gene (n=%d)", nrow(tree_level_ratio)),
       x = "", y = expression("R"^2)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 7))

p_species_comparison <- ggplot(species_comparison_data,
                               aes(x = Model, y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkgreen"),
                    labels = c("NS", "p < 0.05"),
                    name = "") +
  ylim(0, y_limit) +
  labs(title = "B6) Model Comparison",
       subtitle = sprintf("Correlations (n=%d species)", nrow(analysis_ratio)),
       x = "", y = "") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 7))

# ============================================================
# COMBINE LAYOUT (3x2 INDIVIDUAL | 3x2 SPECIES)
# ============================================================

# Left column: Individual tree panels (3 rows)
left_column <- p_tree_mcra / p_tree_pmoa / p_tree_mmox / 
  p_tree_methanotroph / p_tree_ratio / p_tree_comparison +
  plot_layout(ncol = 1)

# Right column: Species panels (3 rows)
right_column <- p_species_mcra / p_species_pmoa / p_species_mmox / 
  p_species_methanotroph / p_species_ratio / p_species_comparison +
  plot_layout(ncol = 1)

# Combine left and right columns
combined_layout <- left_column | right_column +
  plot_annotation(
    title = "Scale-Dependent Mechanisms: Production-Oxidation Balance Emerges at Species Level",
    subtitle = "Left: Individual trees (pseudolog flux) | Right: Species medians (linear flux)",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

ggsave("../../../outputs/figures/Figure_Scale_Dependent_3x2.pdf", combined_layout,
       width = 12, height = 16)
ggsave("../../../outputs/figures/Figure_Scale_Dependent_3x2.png", combined_layout,
       width = 12, height = 16, dpi = 300)

cat("\n============================================================\n")
cat("3x2 SIDE-BY-SIDE LAYOUT FIGURE CREATED\n")
cat("============================================================\n")
cat("\nFiles saved:\n")
cat("  - Figure_Scale_Dependent_3x2.pdf\n")
cat("  - Figure_Scale_Dependent_3x2.png\n")

cat("\n### SUMMARY TABLE ###\n")
summary_table <- rbind(
  tree_comparison_data %>% mutate(Scale = "Individual"),
  species_comparison_data %>% mutate(Scale = "Species")
) %>%
  dplyr::select(Scale, Model, R2, P_value, Significant) %>%
  arrange(Scale, desc(R2))

print(summary_table)

write.csv(summary_table, "../../../outputs/tables/Scale_Dependent_Summary.csv", row.names = FALSE)

cat("\n============================================================\n")