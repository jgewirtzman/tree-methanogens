# ==============================================================================
# Correlation Plots (Figure 3 partial)
# ==============================================================================
# Purpose: Correlation plots between gene abundances and flux with comprehensive
#   analysis. Sources 02_radial_cross_sections.R and 03_threshold_analysis.R.
#
# Pipeline stage: 4 — Publication Figures
# Run after: 00_harmonization/02_harmonize_all_data.R
#
# Inputs:
#   - methanogen_tree_flux_complete_dataset.csv (from data/processed/flux/)
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
# ==============================================================================

library(tidyverse)
library(scales)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(gridExtra)
library(cowplot)

# Load data
ymf2023 <- read.csv("../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv")
ymf2021 <- read.csv('../../data/processed/integrated/merged_tree_dataset_final.csv')
merged_final <- ymf2021  # Assuming ddPCR data is here

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

# ============================================================
# 1) Calculate TOTAL mcrA COPIES for each tree
# ============================================================

# Parse ddPCR data
prepare_long_mcra <- function(df) {
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
        location %in% c("Mineral","mineral") ~ "Mineral",
        location %in% c("Organic","organic") ~ "Organic",
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
        location == "Mineral" ~ "Mineral",
        location == "Organic" ~ "Organic",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(stringency == "loose", !is.na(sample_type)) %>%
    filter(gene == "mcrA", is_probe) %>%
    filter(sample_type %in% c("Heartwood", "Sapwood"))
}

# Function to calculate total copies using numerical integration
calculate_total_copies <- function(mcra_inner, mcra_outer, density_inner, density_outer, dbh_cm, dr = 0.01) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)  # End of inner zone
  r2 <- max(R - 5, r1)  # Start of outer zone
  
  # Create radial grid for integration
  r_edges <- seq(0, R, by = dr)
  if (tail(r_edges, 1) < R) r_edges <- c(r_edges, R)
  r_mid <- 0.5 * (r_edges[-1] + r_edges[-length(r_edges)])
  
  # Shell areas (2π cancels out, but we need it for actual area)
  shell_areas <- pi * (r_edges[-1]^2 - r_edges[-length(r_edges)]^2)  # cm²
  
  # C(r): copies/g with linear transition
  C_r <- ifelse(r_mid <= r1, mcra_inner,
                ifelse(r_mid >= r2, mcra_outer,
                       mcra_inner + (mcra_outer - mcra_inner) * (r_mid - r1) / max(r2 - r1, 1e-9)))
  
  # ρ(r): g/cm³ with linear transition
  # If densities missing, use typical wood density of 0.5 g/cm³
  if (!is.finite(density_inner) || !is.finite(density_outer)) {
    density_inner <- 0.5
    density_outer <- 0.5
  }
  
  rho_r <- ifelse(r_mid <= r1, density_inner,
                  ifelse(r_mid >= r2, density_outer,
                         density_inner + (density_outer - density_inner) * (r_mid - r1) / max(r2 - r1, 1e-9)))
  
  # Total copies per cm of stem height
  # For each shell: copies = C(copies/g) × ρ(g/cm³) × area(cm²) × 1cm height = copies
  total_copies <- sum(C_r * rho_r * shell_areas)
  
  return(total_copies)
}

# Vectorized version
calculate_total_copies_vec <- Vectorize(calculate_total_copies)

# Process mcrA data
long_data <- prepare_long_mcra(merged_final)

# Add species, densities, and calculate total copies per tree
mcra_with_species <- long_data %>%
  left_join(
    merged_final %>% dplyr::select(tree_id, species_id),
    by = "tree_id"
  ) %>%
  mutate(species = unname(species_mapping[species_id])) %>%
  left_join(
    ymf2021 %>% dplyr::select(
      tree_id, 
      dbh,
      density_inner = inner_density_final,
      density_outer = outer_density_final
    ),
    by = "tree_id"
  ) %>%
  filter(is.finite(dbh))

# Calculate inner/outer means and total copies per tree
tree_totals <- mcra_with_species %>%
  group_by(tree_id, species_id, species, dbh, density_inner, density_outer) %>%
  summarise(
    mcra_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
    mcra_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(mcra_inner), is.finite(mcra_outer), is.finite(dbh)) %>%
  mutate(
    mcra_total_copies = calculate_total_copies_vec(mcra_inner, mcra_outer, 
                                                   density_inner, density_outer, dbh),
    mcra_total_log10 = log10(mcra_total_copies + 1),
    species_label = ifelse(is.na(species) | species == "", species_id, species)
  )

# Calculate MEDIAN total copies by species
mcra_by_species <- tree_totals %>%
  group_by(species, species_id) %>%
  summarise(
    n_mcra = n(),
    median_mcra_total = median(mcra_total_copies, na.rm = TRUE),
    median_mcra_total_log = median(mcra_total_log10, na.rm = TRUE),
    mad_mcra_total = mad(mcra_total_copies, na.rm = TRUE),
    q25_mcra = quantile(mcra_total_copies, 0.25, na.rm = TRUE),
    q75_mcra = quantile(mcra_total_copies, 0.75, na.rm = TRUE),
    q25_mcra_log = quantile(mcra_total_log10, 0.25, na.rm = TRUE),
    q75_mcra_log = quantile(mcra_total_log10, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# ============================================================
# 2) Process FLUX DATA
# ============================================================

# Process 2023 flux data
flux_2023 <- ymf2023 %>%
  dplyr::select(
    species_id = Species.Code,
    CH4_flux = CH4_best.flux
  ) %>%
  mutate(
    species = species_mapping[species_id]
  ) %>%
  filter(!is.na(CH4_flux), !is.na(species))

# Process 2021 flux data
flux_2021 <- ymf2021 %>%
  dplyr::select(
    species_id,
    CH4_flux = CH4_best.flux_125cm
  ) %>%
  mutate(
    species = species_mapping[species_id]
  ) %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux), !is.na(species))

# Combine flux data
combined_flux <- bind_rows(flux_2023, flux_2021)

# Calculate MEDIAN flux by species
flux_by_species <- combined_flux %>%
  group_by(species, species_id) %>%
  summarise(
    n_flux = n(),
    median_flux = median(CH4_flux, na.rm = TRUE),
    mad_flux = mad(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# ============================================================
# 3) MERGE AND ANALYZE
# ============================================================

# Merge total mcrA with flux data
correlation_data <- inner_join(
  mcra_by_species,
  flux_by_species,
  by = c("species", "species_id")
) %>%
  filter(n_mcra >= 5, n_flux >= 5) %>%
  arrange(median_flux)

# Calculate correlations
cat("============================================================\n")
cat("TOTAL mcrA COPIES (per cm stem) vs FLUX CORRELATION\n")
cat("============================================================\n")

if(nrow(correlation_data) > 2) {
  # Pearson correlation on log-transformed total copies
  pearson_result <- cor.test(correlation_data$median_mcra_total_log, 
                             correlation_data$median_flux, 
                             method = "pearson")
  
  # Spearman correlation (rank-based, doesn't need log transform)
  spearman_result <- cor.test(correlation_data$median_mcra_total, 
                              correlation_data$median_flux, 
                              method = "spearman")
  
  # Kendall's tau
  kendall_result <- cor.test(correlation_data$median_mcra_total,
                             correlation_data$median_flux,
                             method = "kendall")
  
  cat(sprintf("Number of species: %d\n", nrow(correlation_data)))
  cat(sprintf("Pearson r = %.3f, p = %.4f (on log10-transformed)\n", 
              pearson_result$estimate, pearson_result$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n", 
              spearman_result$estimate, spearman_result$p.value))
  cat(sprintf("Kendall tau = %.3f, p = %.4f\n",
              kendall_result$estimate, kendall_result$p.value))
  
  cat("\nSpecies data:\n")
  print(correlation_data %>% 
          dplyr::select(species, n_mcra, n_flux, median_mcra_total, median_flux) %>%
          mutate(median_mcra_total = scientific(median_mcra_total, digits = 2)) %>%
          arrange(desc(median_mcra_total)))
} else {
  cat("Not enough species with n>=5 for both measurements\n")
  pearson_result <- list(estimate = NA, p.value = NA)
  spearman_result <- list(estimate = NA, p.value = NA)
  kendall_result <- list(estimate = NA, p.value = NA)
}

# ============================================================
# 4) CREATE PLOTS
# ============================================================

# Main correlation plot (log-scale x-axis)
p_main <- ggplot(correlation_data, aes(x = median_mcra_total_log, y = median_flux)) +
  # Add linear regression line with CI
  geom_smooth(method = "lm", se = TRUE, color = "blue", 
              fill = "lightblue", alpha = 0.3, size = 1.2) +
  # Add IQR error bars
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.015, alpha = 0.6, color = "gray40") +
  geom_errorbarh(aes(xmin = q25_mcra_log, xmax = q75_mcra_log),
                 height = 0.003, alpha = 0.6, color = "gray40") +
  # Add points
  geom_point(size = 4, alpha = 0.9, color = "darkblue") +
  # Add species labels
  geom_text_repel(aes(label = species), 
                  size = 3,
                  fontface = "italic",
                  alpha = 0.8,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.size = 0.3,
                  segment.alpha = 0.5,
                  max.overlaps = 20) +
  # Add zero line for flux
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7, size = 0.8) +
  # Add correlation text
  annotate("label", 
           x = Inf,
           y = Inf,
           label = sprintf("Pearson r = %.2f (p = %.3f)\nSpearman ρ = %.2f (p = %.3f)\nKendall τ = %.2f (p = %.3f)",
                           pearson_result$estimate, pearson_result$p.value,
                           spearman_result$estimate, spearman_result$p.value,
                           kendall_result$estimate, kendall_result$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.8,
           fill = "white", alpha = 0.9) +
  # X-axis with log scale labels
  scale_x_continuous(
    breaks = c(2, 3, 4, 5, 6, 7, 8),
    labels = c(expression(10^2), expression(10^3), 
               expression(10^4), expression(10^5), 
               expression(10^6), expression(10^7), expression(10^8))
  ) +
  labs(
    title = "Total mcrA Copies vs CH₄ Flux",
    subtitle = sprintf("n = %d species (≥5 observations each); error bars = IQR", nrow(correlation_data)),
    x = expression("Median total mcrA (copies per cm stem)"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(p_main)

# Pseudolog scale version for better visibility
p_pseudolog <- ggplot(correlation_data, aes(x = median_mcra_total, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", 
              fill = "lightblue", alpha = 0.3, size = 1.2) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.015, alpha = 0.6, color = "gray40") +
  geom_errorbarh(aes(xmin = q25_mcra, xmax = q75_mcra),
                 height = 0.003, alpha = 0.6, color = "gray40") +
  geom_point(size = 4, alpha = 0.9, color = "darkblue") +
  geom_text_repel(aes(label = species), 
                  size = 3,
                  fontface = "italic",
                  alpha = 0.8,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.size = 0.3,
                  segment.alpha = 0.5,
                  max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7, size = 0.8) +
  annotate("label", 
           x = Inf,
           y = Inf,
           label = sprintf("Pearson r = %.2f (p = %.3f)\nSpearman ρ = %.2f (p = %.3f)\nKendall τ = %.2f (p = %.3f)",
                           pearson_result$estimate, pearson_result$p.value,
                           spearman_result$estimate, spearman_result$p.value,
                           kendall_result$estimate, kendall_result$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.8,
           fill = "white", alpha = 0.9) +
  scale_x_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8),
    labels = scales::scientific
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 0.01),
    breaks = c(-0.1, -0.01, 0, 0.01, 0.1, 1),
    labels = c("-0.1", "-0.01", "0", "0.01", "0.1", "1")
  ) +
  labs(
    title = "Total mcrA Copies vs CH₄ Flux (Pseudolog Scale)",
    subtitle = sprintf("n = %d species (≥5 observations each); error bars = IQR", nrow(correlation_data)),
    x = expression("Median total mcrA (copies per cm stem)"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(p_pseudolog)

# Summary statistics
cat("\n============================================================\n")
cat("SUMMARY STATISTICS\n")
cat("============================================================\n")
cat(sprintf("Median total copies across species: %.2e copies/cm\n", 
            median(correlation_data$median_mcra_total)))
cat(sprintf("Range: %.2e to %.2e copies/cm\n", 
            min(correlation_data$median_mcra_total), 
            max(correlation_data$median_mcra_total)))

# Save plot
# ggsave("flux_mcra_total_copies_correlation.png", plot = p_main, 
#        width = 10, height = 8, dpi = 300)

# ============================================================
# SIMPLIFIED: Median Area-Weighted mcrA vs Median CH4 Flux
# ============================================================

library(tidyverse)
library(scales)
library(ggplot2)
library(ggrepel)

# Load data
ymf2023 <- read.csv("../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv")
ymf2021 <- read.csv('../../data/processed/integrated/merged_tree_dataset_final.csv')

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

# ============================================================
# PART 1: Calculate Area-Weighted mcrA for Each Tree
# ============================================================

# Function to parse ddPCR data
prepare_long_mcra <- function(df) {
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
      gene = if_else(gene == "mcra", "mcrA", gene),
      sample_type = case_when(
        location == "Inner" ~ "Heartwood",
        location == "Outer" ~ "Sapwood",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(stringency == "loose", 
           !is.na(sample_type),
           gene == "mcrA", 
           is_probe,
           sample_type %in% c("Heartwood", "Sapwood"))
}

# Area-weighted calculation function
area_weighted_mcra <- function(mcra_inner, mcra_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)  # End of inner zone (5cm from center)
  r2 <- max(R - 5, r1)  # Start of outer zone (5cm from bark)
  
  # Calculate weighting factor S
  S <- r1^2 + r1*r2 + r2^2
  
  # Area-weighted average
  mcra_weighted <- mcra_outer + (mcra_inner - mcra_outer) * (S / (3 * R^2))
  
  return(mcra_weighted)
}

# Process mcrA data
long_data <- prepare_long_mcra(ymf2021)

# Calculate area-weighted mcrA for each tree
tree_mcra <- long_data %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, species_id, dbh),
    by = "tree_id"
  ) %>%
  filter(is.finite(dbh)) %>%
  group_by(tree_id, species_id, dbh) %>%
  summarise(
    mcra_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
    mcra_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(mcra_inner), is.finite(mcra_outer)) %>%
  mutate(
    mcra_area_weighted = mapply(area_weighted_mcra, mcra_inner, mcra_outer, dbh),
    species = species_mapping[species_id]
  )

# Calculate MEDIAN area-weighted mcrA by species
mcra_by_species <- tree_mcra %>%
  group_by(species, species_id) %>%
  summarise(
    n_mcra = n(),
    median_mcra = median(mcra_area_weighted, na.rm = TRUE),
    median_mcra_log10 = log10(median_mcra + 1),
    q25_mcra = quantile(mcra_area_weighted, 0.25, na.rm = TRUE),
    q75_mcra = quantile(mcra_area_weighted, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# ============================================================
# PART 2: Process Flux Data
# ============================================================

# Combine 2023 and 2021 flux data
flux_all <- bind_rows(
  # 2023 data
  ymf2023 %>%
    dplyr::select(species_id = Species.Code, CH4_flux = CH4_best.flux) %>%
    mutate(species = species_mapping[species_id]),
  # 2021 data  
  ymf2021 %>%
    dplyr::select(species_id, CH4_flux = CH4_best.flux_125cm) %>%
    mutate(species = species_mapping[species_id])
) %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux), !is.na(species))

# Calculate MEDIAN flux by species
flux_by_species <- flux_all %>%
  group_by(species, species_id) %>%
  summarise(
    n_flux = n(),
    median_flux = median(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# ============================================================
# PART 3: Merge and Analyze
# ============================================================

# Merge data, requiring n≥5 for both measurements
analysis_data <- inner_join(
  mcra_by_species,
  flux_by_species,
  by = c("species", "species_id")
) %>%
  filter(n_mcra >= 5, n_flux >= 5)

# Calculate correlations
if(nrow(analysis_data) > 2) {
  pearson <- cor.test(analysis_data$median_mcra_log10, analysis_data$median_flux)
  spearman <- cor.test(analysis_data$median_mcra, analysis_data$median_flux, method = "spearman")
  kendall <- cor.test(analysis_data$median_mcra, analysis_data$median_flux, method = "kendall")
  
  cat("============================================================\n")
  cat("MEDIAN AREA-WEIGHTED mcrA vs MEDIAN FLUX CORRELATION\n")
  cat("============================================================\n")
  cat(sprintf("Number of species: %d\n", nrow(analysis_data)))
  cat(sprintf("Pearson r = %.3f, p = %.4f (log-transformed mcrA)\n", 
              pearson$estimate, pearson$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n", 
              spearman$estimate, spearman$p.value))
  cat(sprintf("Kendall tau = %.3f, p = %.4f\n",
              kendall$estimate, kendall$p.value))
  
  cat("\nSpecies data (sorted by flux):\n")
  print(analysis_data %>% 
          dplyr::select(species, n_mcra, n_flux, median_mcra, median_flux) %>%
          arrange(median_flux) %>%
          mutate(median_mcra = round(median_mcra, 1),
                 median_flux = round(median_flux, 4)))
}

# ============================================================
# PART 4: Create Plot
# ============================================================

# Main plot with log-transformed x-axis
p <- ggplot(analysis_data, aes(x = median_mcra_log10, y = median_flux)) +
  # Regression line (semi-transparent)
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", 
              fill = "lightblue", alpha = 0.2, size = 1) +
  # Error bars (IQR)
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(q25_mcra + 1), xmax = log10(q75_mcra + 1)),
                 height = 0.003, alpha = 0.5, color = "gray40") +
  # Points
  geom_point(size = 4, alpha = 0.85, color = "darkblue") +
  # Species labels (on top with better visibility)
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  # Correlation annotation
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nρ = %.2f, p = %.3f\nτ = %.2f, p = %.3f",
                           pearson$estimate, pearson$p.value,
                           spearman$estimate, spearman$p.value,
                           kendall$estimate, kendall$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  # Axis formatting
  scale_x_continuous(
    breaks = 0:5,
    labels = c(expression(10^0), expression(10^1), expression(10^2), 
               expression(10^3), expression(10^4), expression(10^5))
  ) +
  # Labels
  labs(
    #title = "Median Area-weighted mcrA vs Median CH₄ Flux by Species",
    #subtitle = sprintf("n = %d species (≥5 observations each); error bars show IQR", 
                       #nrow(analysis_data)),
    x = expression("Median area-weighted mcrA (copies g"^-1*")"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  # Theme
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", size = 0.3),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

print(p)

source('../../code/06_figures/02_radial_cross_sections.R')
source('../../code/06_figures/03_threshold_analysis.R')

combined_plot<-grid.arrange(p_main, p, nrow=2)

combined_plot <- grid.arrange(
  p_main, 
  p, p4+theme(legend.position="none"), 
  nrow = 2,
  ncol = 2,
  layout_matrix = rbind(c(1, 1),
                        c(2, 3))
)

combined_plot

# Save plot (uncomment to use)
# ggsave("median_areaweighted_mcra_vs_flux.png", plot = p, 
#        width = 10, height = 8, dpi = 300)


# ============================================================
# ADDITIONAL ANALYSIS: mcrA Copies per Unit Surface Area
# ============================================================

cat("\n\n============================================================\n")
cat("ADDITIONAL: mcrA COPIES PER SURFACE AREA ANALYSIS\n")
cat("============================================================\n")

# Function to calculate copies per unit surface area
calculate_copies_per_surface <- function(mcra_inner, mcra_outer, density_inner, density_outer, dbh_cm, dr = 0.01) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  
  # Create radial grid for integration
  r_edges <- seq(0, R, by = dr)
  if (tail(r_edges, 1) < R) r_edges <- c(r_edges, R)
  r_mid <- 0.5 * (r_edges[-1] + r_edges[-length(r_edges)])
  
  # Shell areas (cm²)
  shell_areas <- pi * (r_edges[-1]^2 - r_edges[-length(r_edges)]^2)
  
  # C(r): copies/g with linear transition
  C_r <- ifelse(r_mid <= r1, mcra_inner,
                ifelse(r_mid >= r2, mcra_outer,
                       mcra_inner + (mcra_outer - mcra_inner) * (r_mid - r1) / max(r2 - r1, 1e-9)))
  
  # ρ(r): g/cm³ with linear transition (default 0.5 if missing)
  if (!is.finite(density_inner) || !is.finite(density_outer)) {
    density_inner <- 0.5
    density_outer <- 0.5
  }
  
  rho_r <- ifelse(r_mid <= r1, density_inner,
                  ifelse(r_mid >= r2, density_outer,
                         density_inner + (density_outer - density_inner) * (r_mid - r1) / max(r2 - r1, 1e-9)))
  
  # Total copies in the cross-section (for 1 cm height)
  total_copies <- sum(C_r * rho_r * shell_areas)
  
  # Circumference (surface area for 1 cm height)
  circumference <- 2 * pi * R  # cm
  
  # Copies per cm of surface area
  copies_per_surface <- total_copies / circumference
  
  return(copies_per_surface)
}

# Get density data and calculate copies per surface area
tree_surface <- tree_mcra %>%
  left_join(
    ymf2021 %>% dplyr::select(
      tree_id,
      density_inner = inner_density_final,
      density_outer = outer_density_final
    ),
    by = "tree_id"
  ) %>%
  mutate(
    mcra_per_surface = mapply(calculate_copies_per_surface, 
                              mcra_inner, mcra_outer,
                              density_inner, density_outer, dbh),
    mcra_surface_log10 = log10(mcra_per_surface + 1)
  )

# Calculate MEDIAN copies per surface by species
surface_by_species <- tree_surface %>%
  group_by(species, species_id) %>%
  summarise(
    n_mcra = n(),
    median_surface = median(mcra_per_surface, na.rm = TRUE),
    median_surface_log10 = log10(median_surface + 1),
    q25_surface = quantile(mcra_per_surface, 0.25, na.rm = TRUE),
    q75_surface = quantile(mcra_per_surface, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# Merge with flux data
surface_analysis <- inner_join(
  surface_by_species,
  flux_by_species,
  by = c("species", "species_id")
) %>%
  filter(n_mcra >= 5, n_flux >= 5)

# Calculate correlations for surface area metric
if(nrow(surface_analysis) > 2) {
  pearson_surf <- cor.test(surface_analysis$median_surface_log10, surface_analysis$median_flux)
  spearman_surf <- cor.test(surface_analysis$median_surface, surface_analysis$median_flux, method = "spearman")
  kendall_surf <- cor.test(surface_analysis$median_surface, surface_analysis$median_flux, method = "kendall")
  
  cat(sprintf("Number of species: %d\n", nrow(surface_analysis)))
  cat(sprintf("Pearson r = %.3f, p = %.4f (log-transformed)\n", 
              pearson_surf$estimate, pearson_surf$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n", 
              spearman_surf$estimate, spearman_surf$p.value))
  cat(sprintf("Kendall tau = %.3f, p = %.4f\n",
              kendall_surf$estimate, kendall_surf$p.value))
  
  cat("\nSpecies data (sorted by flux):\n")
  print(surface_analysis %>% 
          dplyr::select(species, n_mcra, n_flux, median_surface, median_flux) %>%
          arrange(median_flux) %>%
          mutate(median_surface = scientific(median_surface, digits = 2),
                 median_flux = round(median_flux, 4)))
}

# Create plot for surface area analysis
p_surface <- ggplot(surface_analysis, aes(x = median_surface_log10, y = median_flux)) +
  # Regression line (semi-transparent)
  geom_smooth(method = "lm", se = TRUE, color = "lightsteelblue", 
              fill = "lightblue", alpha = 0.2, size = 1) +
  # Error bars (IQR)
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(q25_surface + 1), xmax = log10(q75_surface + 1)),
                 height = 0.003, alpha = 0.5, color = "gray40") +
  # Points
  geom_point(size = 4, alpha = 0.85, color = "darkblue") +
  # Species labels (on top with better visibility)
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  # Correlation annotation
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nρ = %.2f, p = %.3f\nτ = %.2f, p = %.3f",
                           pearson_surf$estimate, pearson_surf$p.value,
                           spearman_surf$estimate, spearman_surf$p.value,
                           kendall_surf$estimate, kendall_surf$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  # Axis formatting
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2), expression(10^3), 
               expression(10^4), expression(10^5), expression(10^6))
  ) +
  # Labels
  labs(
    #title = "Median mcrA Copies per Surface Area vs Median CH₄ Flux",
    #subtitle = sprintf("n = %d species (≥5 observations each); error bars show IQR", 
                       #nrow(surface_analysis)),
    x = expression("Median mcrA in tissue column beneath 1 cm"^2*" bark (total copies)"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  # Theme
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", size = 0.3),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

print(p_surface)

# Summary
cat("\n============================================================\n")
cat("COMPARISON OF METRICS\n")
cat("============================================================\n")
cat("Area-weighted (copies/g):\n")
cat(sprintf("  Pearson r = %.3f (p = %.4f)\n", pearson$estimate, pearson$p.value))
cat("\nPer surface area (copies/cm):\n")
cat(sprintf("  Pearson r = %.3f (p = %.4f)\n", pearson_surf$estimate, pearson_surf$p.value))
cat(sprintf("\nMedian copies per surface: %.2e copies/cm\n", 
            median(surface_analysis$median_surface)))
cat(sprintf("Range: %.2e to %.2e copies/cm\n", 
            min(surface_analysis$median_surface), 
            max(surface_analysis$median_surface)))


# ============================================================
# MEAN ± SE VERSIONS OF BOTH PLOTS
# ============================================================

cat("\n\n============================================================\n")
cat("MEAN ± SE ANALYSIS\n")
cat("============================================================\n")

# Calculate MEAN area-weighted mcrA by species with SE
mcra_by_species_mean <- tree_mcra %>%
  group_by(species, species_id) %>%
  summarise(
    n_mcra = n(),
    mean_mcra = mean(mcra_area_weighted, na.rm = TRUE),
    se_mcra = sd(mcra_area_weighted, na.rm = TRUE) / sqrt(n()),
    mean_mcra_log10 = log10(mean_mcra + 1),
    se_mcra_lower = mean_mcra - se_mcra,
    se_mcra_upper = mean_mcra + se_mcra,
    .groups = 'drop'
  )

# Calculate MEAN flux by species with SE
flux_by_species_mean <- flux_all %>%
  group_by(species, species_id) %>%
  summarise(
    n_flux = n(),
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    se_flux = sd(CH4_flux, na.rm = TRUE) / sqrt(n()),
    se_flux_lower = mean_flux - se_flux,
    se_flux_upper = mean_flux + se_flux,
    .groups = 'drop'
  )

# Merge data for mean analysis
analysis_data_mean <- inner_join(
  mcra_by_species_mean,
  flux_by_species_mean,
  by = c("species", "species_id")
) %>%
  filter(n_mcra >= 5, n_flux >= 5)

# Calculate correlations for mean data
if(nrow(analysis_data_mean) > 2) {
  pearson_mean <- cor.test(analysis_data_mean$mean_mcra_log10, analysis_data_mean$mean_flux)
  spearman_mean <- cor.test(analysis_data_mean$mean_mcra, analysis_data_mean$mean_flux, method = "spearman")
  kendall_mean <- cor.test(analysis_data_mean$mean_mcra, analysis_data_mean$mean_flux, method = "kendall")
  
  cat("MEAN AREA-WEIGHTED mcrA vs MEAN FLUX CORRELATION\n")
  cat("============================================================\n")
  cat(sprintf("Number of species: %d\n", nrow(analysis_data_mean)))
  cat(sprintf("Pearson r = %.3f, p = %.4f (log-transformed mcrA)\n", 
              pearson_mean$estimate, pearson_mean$p.value))
  cat(sprintf("R-squared = %.3f\n", pearson_mean$estimate^2))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n", 
              spearman_mean$estimate, spearman_mean$p.value))
  cat(sprintf("Kendall tau = %.3f, p = %.4f\n",
              kendall_mean$estimate, kendall_mean$p.value))
  
  cat("\nSpecies data (sorted by flux):\n")
  print(analysis_data_mean %>% 
          dplyr::select(species, n_mcra, n_flux, mean_mcra, mean_flux) %>%
          arrange(mean_flux) %>%
          mutate(mean_mcra = round(mean_mcra, 1),
                 mean_flux = round(mean_flux, 4)))
}

# PLOT 1: Mean Area-weighted mcrA vs Mean Flux
p_mean <- ggplot(analysis_data_mean, aes(x = mean_mcra_log10, y = mean_flux)) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "lightblue", 
              fill = "lightblue", alpha = 0.2, size = 1) +
  # Error bars (SE)
  geom_errorbar(aes(ymin = se_flux_lower, ymax = se_flux_upper),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(pmax(se_mcra_lower, 0) + 1), 
                     xmax = log10(se_mcra_upper + 1)),
                 height = 0.003, alpha = 0.5, color = "gray40") +
  # Points
  geom_point(size = 4, alpha = 0.85, color = "darkblue") +
  # Species labels
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  # Correlation annotation with R²
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR² = %.3f\nρ = %.2f, p = %.3f\nτ = %.2f, p = %.3f",
                           pearson_mean$estimate, pearson_mean$p.value,
                           pearson_mean$estimate^2,
                           spearman_mean$estimate, spearman_mean$p.value,
                           kendall_mean$estimate, kendall_mean$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  # Axis formatting
  scale_x_continuous(
    breaks = 0:5,
    labels = c(expression(10^0), expression(10^1), expression(10^2), 
               expression(10^3), expression(10^4), expression(10^5))
  ) +
  # Labels
  labs(
    x = expression("Mean area-weighted mcrA (copies g"^-1*")"),
    y = expression("Mean CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  # Theme
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", size = 0.3),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

print(p_mean)

# ============================================================
# MEAN ± SE for Surface Area Analysis
# ============================================================

cat("\n============================================================\n")
cat("MEAN mcrA COPIES PER SURFACE AREA ANALYSIS\n")
cat("============================================================\n")

# Calculate MEAN copies per surface by species with SE
surface_by_species_mean <- tree_surface %>%
  group_by(species, species_id) %>%
  summarise(
    n_mcra = n(),
    mean_surface = mean(mcra_per_surface, na.rm = TRUE),
    se_surface = sd(mcra_per_surface, na.rm = TRUE) / sqrt(n()),
    mean_surface_log10 = log10(mean_surface + 1),
    se_surface_lower = mean_surface - se_surface,
    se_surface_upper = mean_surface + se_surface,
    .groups = 'drop'
  )

# Merge with flux data
surface_analysis_mean <- inner_join(
  surface_by_species_mean,
  flux_by_species_mean,
  by = c("species", "species_id")
) %>%
  filter(n_mcra >= 5, n_flux >= 5)

# Calculate correlations for surface area metric with means
if(nrow(surface_analysis_mean) > 2) {
  pearson_surf_mean <- cor.test(surface_analysis_mean$mean_surface_log10, 
                                surface_analysis_mean$mean_flux)
  spearman_surf_mean <- cor.test(surface_analysis_mean$mean_surface, 
                                 surface_analysis_mean$mean_flux, method = "spearman")
  kendall_surf_mean <- cor.test(surface_analysis_mean$mean_surface, 
                                surface_analysis_mean$mean_flux, method = "kendall")
  
  cat(sprintf("Number of species: %d\n", nrow(surface_analysis_mean)))
  cat(sprintf("Pearson r = %.3f, p = %.4f (log-transformed)\n", 
              pearson_surf_mean$estimate, pearson_surf_mean$p.value))
  cat(sprintf("R-squared = %.3f\n", pearson_surf_mean$estimate^2))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n", 
              spearman_surf_mean$estimate, spearman_surf_mean$p.value))
  cat(sprintf("Kendall tau = %.3f, p = %.4f\n",
              kendall_surf_mean$estimate, kendall_surf_mean$p.value))
  
  cat("\nSpecies data (sorted by flux):\n")
  print(surface_analysis_mean %>% 
          dplyr::select(species, n_mcra, n_flux, mean_surface, mean_flux) %>%
          arrange(mean_flux) %>%
          mutate(mean_surface = scientific(mean_surface, digits = 2),
                 mean_flux = round(mean_flux, 4)))
}

# PLOT 2: Mean Surface Area Analysis
p_surface_mean <- ggplot(surface_analysis_mean, aes(x = mean_surface_log10, y = mean_flux)) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "blue", 
              fill = "lightblue", alpha = 0.2, size = 1) +
  # Error bars (SE)
  geom_errorbar(aes(ymin = se_flux_lower, ymax = se_flux_upper),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(pmax(se_surface_lower, 0) + 1), 
                     xmax = log10(se_surface_upper + 1)),
                 height = 0.003, alpha = 0.5, color = "gray40") +
  # Points
  geom_point(size = 4, alpha = 0.85, color = "darkblue") +
  # Species labels
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  # Correlation annotation with R²
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR² = %.3f\nρ = %.2f, p = %.3f\nτ = %.2f, p = %.3f",
                           pearson_surf_mean$estimate, pearson_surf_mean$p.value,
                           pearson_surf_mean$estimate^2,
                           spearman_surf_mean$estimate, spearman_surf_mean$p.value,
                           kendall_surf_mean$estimate, kendall_surf_mean$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  # Axis formatting
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2), expression(10^3), 
               expression(10^4), expression(10^5), expression(10^6))
  ) +
  # Labels
  labs(
    x = expression("Mean mcrA in tissue column beneath 1 cm"^2*" bark (total copies)"),
    y = expression("Mean CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  # Theme
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", size = 0.3),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

print(p_surface_mean)

# Summary comparison
cat("\n============================================================\n")
cat("SUMMARY: MEAN-BASED CORRELATIONS\n")
cat("============================================================\n")
cat("Area-weighted (copies/g):\n")
cat(sprintf("  Pearson r = %.3f (R² = %.3f), p = %.4f\n", 
            pearson_mean$estimate, pearson_mean$estimate^2, pearson_mean$p.value))
cat("\nPer surface area (copies/cm):\n")
cat(sprintf("  Pearson r = %.3f (R² = %.3f), p = %.4f\n", 
            pearson_surf_mean$estimate, pearson_surf_mean$estimate^2, pearson_surf_mean$p.value))
cat(sprintf("\nMean copies per surface: %.2e copies/cm\n", 
            mean(surface_analysis_mean$mean_surface)))
cat(sprintf("Range: %.2e to %.2e copies/cm\n", 
            min(surface_analysis_mean$mean_surface), 
            max(surface_analysis_mean$mean_surface)))

# Save plots (uncomment to use)
# ggsave("mean_areaweighted_mcra_vs_flux.png", plot = p_mean, 
#        width = 10, height = 8, dpi = 300)
# ggsave("mean_surface_mcra_vs_flux.png", plot = p_surface_mean, 
#        width = 10, height = 8, dpi = 300)


# ============================================================
# MEDIAN ± SE OF MEDIAN (BOOTSTRAP AND FORMULA APPROACHES)
# ============================================================

cat("\n\n============================================================\n")
cat("MEDIAN ± SE OF MEDIAN ANALYSIS\n")
cat("============================================================\n")

library(boot)  # For bootstrapping

# Function to calculate SE of median using bootstrap
bootstrap_median_se <- function(x, n_boot = 1000) {
  if(length(x) < 5) return(NA)
  
  boot_fun <- function(data, indices) {
    median(data[indices], na.rm = TRUE)
  }
  
  boot_result <- boot(x[!is.na(x)], boot_fun, R = n_boot)
  return(sd(boot_result$t))
}

# Function for formula-based SE of median (asymptotic approximation)
# SE(median) ≈ 1.253 * SD / sqrt(n) for normal distributions
# More robust: SE(median) = 1/(2*f(median)*sqrt(n)) where f is density
formula_median_se <- function(x) {
  x_clean <- x[!is.na(x)]
  n <- length(x_clean)
  if(n < 5) return(NA)
  
  # Using the 1.253 approximation (assumes approximately normal)
  return(1.253 * sd(x_clean) / sqrt(n))
}

# Calculate MEDIAN area-weighted mcrA with SE of median
mcra_by_species_median_se <- tree_mcra %>%
  group_by(species, species_id) %>%
  summarise(
    n_mcra = n(),
    median_mcra = median(mcra_area_weighted, na.rm = TRUE),
    median_mcra_log10 = log10(median_mcra + 1),
    # Bootstrap SE of median
    se_median_mcra = bootstrap_median_se(mcra_area_weighted),
    # Alternative: formula-based SE
    se_median_mcra_formula = formula_median_se(mcra_area_weighted),
    .groups = 'drop'
  ) %>%
  mutate(
    # Using bootstrap SE for error bars
    se_mcra_lower = pmax(0, median_mcra - se_median_mcra),
    se_mcra_upper = median_mcra + se_median_mcra
  )

# Calculate MEDIAN flux with SE of median
flux_by_species_median_se <- flux_all %>%
  group_by(species, species_id) %>%
  summarise(
    n_flux = n(),
    median_flux = median(CH4_flux, na.rm = TRUE),
    # Bootstrap SE of median
    se_median_flux = bootstrap_median_se(CH4_flux),
    # Alternative: formula-based SE
    se_median_flux_formula = formula_median_se(CH4_flux),
    .groups = 'drop'
  ) %>%
  mutate(
    se_flux_lower = median_flux - se_median_flux,
    se_flux_upper = median_flux + se_median_flux
  )

# Merge data
analysis_data_median_se <- inner_join(
  mcra_by_species_median_se,
  flux_by_species_median_se,
  by = c("species", "species_id")
) %>%
  filter(n_mcra >= 5, n_flux >= 5)

# Calculate correlations
if(nrow(analysis_data_median_se) > 2) {
  pearson_med_se <- cor.test(analysis_data_median_se$median_mcra_log10, 
                             analysis_data_median_se$median_flux)
  spearman_med_se <- cor.test(analysis_data_median_se$median_mcra, 
                              analysis_data_median_se$median_flux, method = "spearman")
  kendall_med_se <- cor.test(analysis_data_median_se$median_mcra, 
                             analysis_data_median_se$median_flux, method = "kendall")
  
  cat("MEDIAN AREA-WEIGHTED mcrA vs MEDIAN FLUX (with SE of median)\n")
  cat("============================================================\n")
  cat(sprintf("Number of species: %d\n", nrow(analysis_data_median_se)))
  cat(sprintf("Pearson r = %.3f, p = %.4f (log-transformed mcrA)\n", 
              pearson_med_se$estimate, pearson_med_se$p.value))
  cat(sprintf("R-squared = %.3f\n", pearson_med_se$estimate^2))
  
  cat("\nComparison of SE calculation methods:\n")
  print(analysis_data_median_se %>%
          dplyr::select(species, se_median_mcra, se_median_mcra_formula, 
                 se_median_flux, se_median_flux_formula) %>%
          mutate(across(where(is.numeric), ~round(., 4))))
}

# PLOT 1: Median Area-weighted with SE of median
p_median_se <- ggplot(analysis_data_median_se, aes(x = median_mcra_log10, y = median_flux)) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "blue", 
              fill = "lightblue", alpha = 0.2, size = 1) +
  # Error bars (SE of median)
  geom_errorbar(aes(ymin = se_flux_lower, ymax = se_flux_upper),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(se_mcra_lower + 1), 
                     xmax = log10(se_mcra_upper + 1)),
                 height = 0.003, alpha = 0.5, color = "gray40") +
  # Points
  geom_point(size = 4, alpha = 0.85, color = "darkblue") +
  # Species labels
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  # Annotation
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR² = %.3f\nρ = %.2f, p = %.3f\nτ = %.2f, p = %.3f",
                           pearson_med_se$estimate, pearson_med_se$p.value,
                           pearson_med_se$estimate^2,
                           spearman_med_se$estimate, spearman_med_se$p.value,
                           kendall_med_se$estimate, kendall_med_se$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  # Axis formatting
  scale_x_continuous(
    breaks = 0:5,
    labels = c(expression(10^0), expression(10^1), expression(10^2), 
               expression(10^3), expression(10^4), expression(10^5))
  ) +
  # Labels
  labs(
    x = expression("Median area-weighted mcrA (copies g"^-1*")"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  # Theme
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", size = 0.3),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

print(p_median_se)

# ============================================================
# Surface Area Analysis with SE of median
# ============================================================

cat("\n============================================================\n")
cat("MEDIAN mcrA COPIES PER SURFACE AREA (with SE of median)\n")
cat("============================================================\n")

# Calculate MEDIAN copies per surface with SE of median
surface_by_species_median_se <- tree_surface %>%
  group_by(species, species_id) %>%
  summarise(
    n_mcra = n(),
    median_surface = median(mcra_per_surface, na.rm = TRUE),
    median_surface_log10 = log10(median_surface + 1),
    # Bootstrap SE of median
    se_median_surface = bootstrap_median_se(mcra_per_surface),
    # Alternative: formula-based SE
    se_median_surface_formula = formula_median_se(mcra_per_surface),
    .groups = 'drop'
  ) %>%
  mutate(
    se_surface_lower = pmax(0, median_surface - se_median_surface),
    se_surface_upper = median_surface + se_median_surface
  )

# Merge with flux data
surface_analysis_median_se <- inner_join(
  surface_by_species_median_se,
  flux_by_species_median_se,
  by = c("species", "species_id")
) %>%
  filter(n_mcra >= 5, n_flux >= 5)

# Calculate correlations
if(nrow(surface_analysis_median_se) > 2) {
  pearson_surf_med_se <- cor.test(surface_analysis_median_se$median_surface_log10, 
                                  surface_analysis_median_se$median_flux)
  spearman_surf_med_se <- cor.test(surface_analysis_median_se$median_surface, 
                                   surface_analysis_median_se$median_flux, method = "spearman")
  kendall_surf_med_se <- cor.test(surface_analysis_median_se$median_surface, 
                                  surface_analysis_median_se$median_flux, method = "kendall")
  
  cat(sprintf("Number of species: %d\n", nrow(surface_analysis_median_se)))
  cat(sprintf("Pearson r = %.3f, p = %.4f (log-transformed)\n", 
              pearson_surf_med_se$estimate, pearson_surf_med_se$p.value))
  cat(sprintf("R-squared = %.3f\n", pearson_surf_med_se$estimate^2))
}

# PLOT 2: Median Surface Area with SE of median
p_surface_median_se <- ggplot(surface_analysis_median_se, aes(x = median_surface_log10, y = median_flux)) +
  # Regression line
  geom_smooth(method = "lm", se = TRUE, color = "blue", 
              fill = "lightblue", alpha = 0.2, size = 1) +
  # Error bars (SE of median)
  geom_errorbar(aes(ymin = se_flux_lower, ymax = se_flux_upper),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(se_surface_lower + 1), 
                     xmax = log10(se_surface_upper + 1)),
                 height = 0.003, alpha = 0.5, color = "gray40") +
  # Points
  geom_point(size = 4, alpha = 0.85, color = "darkblue") +
  # Species labels
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  # Annotation
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR² = %.3f\nρ = %.2f, p = %.3f\nτ = %.2f, p = %.3f",
                           pearson_surf_med_se$estimate, pearson_surf_med_se$p.value,
                           pearson_surf_med_se$estimate^2,
                           spearman_surf_med_se$estimate, spearman_surf_med_se$p.value,
                           kendall_surf_med_se$estimate, kendall_surf_med_se$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  # Axis formatting
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2), expression(10^3), 
               expression(10^4), expression(10^5), expression(10^6))
  ) +
  # Labels
  labs(
    x = expression("Median mcrA in tissue column beneath 1 cm"^2*" bark (total copies)"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  # Theme
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", size = 0.3),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

print(p_surface_median_se)

# Summary
cat("\n============================================================\n")
cat("SUMMARY: MEDIAN WITH SE OF MEDIAN\n")
cat("============================================================\n")
cat("Area-weighted (copies/g):\n")
cat(sprintf("  Pearson r = %.3f (R² = %.3f), p = %.4f\n", 
            pearson_med_se$estimate, pearson_med_se$estimate^2, pearson_med_se$p.value))
cat("\nPer surface area (copies/cm):\n")
cat(sprintf("  Pearson r = %.3f (R² = %.3f), p = %.4f\n", 
            pearson_surf_med_se$estimate, pearson_surf_med_se$estimate^2, pearson_surf_med_se$p.value))






























# ============================================================
# STATISTICS FOR MANUSCRIPT TEXT
# Add this section after the correlation analysis
# ============================================================

cat("\n============================================================\n")
cat("DETAILED STATISTICS FOR MANUSCRIPT\n")
cat("============================================================\n")

# 1. Area-weighted correlations at species level
mcra_by_species <- tree_weighted %>%
  group_by(species_id, species_label) %>%
  summarise(
    n_trees = n(),
    median_mcra = median(mcra_area_weighted, na.rm = TRUE),
    mean_mcra = mean(mcra_area_weighted, na.rm = TRUE),
    median_mcra_log = median(log10(mcra_area_weighted + 1), na.rm = TRUE),
    mean_mcra_log = mean(log10(mcra_area_weighted + 1), na.rm = TRUE),
    .groups = 'drop'
  )

# Get flux data by species
flux_by_species_full <- bind_rows(
  ymf2023 %>%
    dplyr::select(species_id = Species.Code, CH4_flux = CH4_best.flux),
  ymf2021 %>%
    dplyr::select(species_id, CH4_flux = CH4_best.flux_125cm)
) %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux)) %>%
  group_by(species_id) %>%
  summarise(
    n_flux = n(),
    median_flux = median(CH4_flux, na.rm = TRUE),
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    .groups = 'drop'
  )

# Merge for correlation analysis
correlation_data <- inner_join(
  mcra_by_species,
  flux_by_species_full,
  by = "species_id"
) %>%
  filter(n_trees >= 5, n_flux >= 5)

# Calculate linear model and R²
lm_median <- lm(median_flux ~ median_mcra_log, data = correlation_data)
r2_median <- summary(lm_median)$r.squared
p_median <- summary(lm_median)$coefficients[2,4]

lm_mean <- lm(mean_flux ~ mean_mcra_log, data = correlation_data)
r2_mean <- summary(lm_mean)$r.squared
p_mean <- summary(lm_mean)$coefficients[2,4]

cat("\n1. AREA-WEIGHTED mcrA vs FLUX CORRELATIONS:\n")
cat("--------------------------------------------\n")
cat(sprintf("Number of species with n≥5 for both: %d\n", nrow(correlation_data)))
cat(sprintf("\nUsing MEDIAN values:\n"))
cat(sprintf("  R² = %.3f, p = %.4f (log-transformed mcrA)\n", r2_median, p_median))
cat(sprintf("  Pearson r = %.3f\n", sqrt(r2_median)))
cat(sprintf("\nUsing MEAN values:\n"))
cat(sprintf("  R² = %.3f, p = %.4f (log-transformed mcrA)\n", r2_mean, p_mean))
cat(sprintf("  Pearson r = %.3f\n", sqrt(r2_mean)))

# Test robustness to transformations
cor_untransformed <- cor.test(correlation_data$median_mcra, correlation_data$median_flux)
cor_sqrt <- cor.test(sqrt(correlation_data$median_mcra), correlation_data$median_flux)
cor_log_both <- cor.test(log10(correlation_data$median_mcra + 1), 
                         log10(abs(correlation_data$median_flux) + 1))

cat("\n2. ROBUSTNESS TO TRANSFORMATIONS (median values):\n")
cat("--------------------------------------------------\n")
cat(sprintf("Untransformed: r = %.3f, p = %.4f\n", 
            cor_untransformed$estimate, cor_untransformed$p.value))
cat(sprintf("Square root mcrA: r = %.3f, p = %.4f\n", 
            cor_sqrt$estimate, cor_sqrt$p.value))
cat(sprintf("Log-log: r = %.3f, p = %.4f\n", 
            cor_log_both$estimate, cor_log_both$p.value))

# 3. Total mcrA copies per cm² bark surface calculation
# This would be the total copies in a 1 cm tall cylinder beneath 1 cm² of bark
calculate_total_per_cm2 <- function(mcra_inner, mcra_outer, density_inner, density_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  # Use default density if missing
  if (!is.finite(density_inner)) density_inner <- 0.5
  if (!is.finite(density_outer)) density_outer <- 0.5
  
  # Total copies per cm of stem height (from earlier calculation)
  # Then divide by circumference to get per cm² of bark surface
  circumference <- 2 * pi * R
  
  # This is simplified - you may want to use the full integration method
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  
  # Area-weighted average concentration
  S <- r1^2 + r1*r2 + r2^2
  mcra_weighted <- mcra_outer + (mcra_inner - mcra_outer) * (S / (3 * R^2))
  
  # Total mass in 1 cm height cylinder
  volume_cm3 <- pi * R^2  # cm³ for 1 cm height
  avg_density <- (density_inner + density_outer) / 2  # simplified
  mass_g <- volume_cm3 * avg_density
  
  # Total copies in cylinder
  total_copies <- mcra_weighted * mass_g
  
  # Per cm² of bark surface
  copies_per_cm2 <- total_copies / circumference
  
  return(copies_per_cm2)
}

# Add density data if available
tree_with_density <- tree_weighted %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, 
                       density_inner = inner_density_final,
                       density_outer = outer_density_final),
    by = "tree_id"
  ) %>%
  mutate(
    mcra_per_cm2 = mapply(calculate_total_per_cm2, 
                          mcra_inner, mcra_outer, 
                          density_inner, density_outer, dbh)
  )

# Calculate species-level statistics
species_per_cm2 <- tree_with_density %>%
  filter(!is.na(mcra_per_cm2)) %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    median_per_cm2 = median(mcra_per_cm2, na.rm = TRUE),
    median_per_cm2_log = median(log10(mcra_per_cm2 + 1), na.rm = TRUE),
    .groups = 'drop'
  )

# Merge with flux data
if(nrow(species_per_cm2) > 0) {
  correlation_per_cm2 <- inner_join(
    species_per_cm2,
    flux_by_species_full,
    by = "species_id"
  ) %>%
    filter(n >= 5, n_flux >= 5)
  
  if(nrow(correlation_per_cm2) > 2) {
    lm_per_cm2 <- lm(median_flux ~ median_per_cm2_log, data = correlation_per_cm2)
    r2_per_cm2 <- summary(lm_per_cm2)$r.squared
    p_per_cm2 <- summary(lm_per_cm2)$coefficients[2,4]
    
    cat("\n3. TOTAL mcrA PER CM² BARK SURFACE:\n")
    cat("------------------------------------\n")
    cat(sprintf("Number of species: %d\n", nrow(correlation_per_cm2)))
    cat(sprintf("R² = %.3f, p = %.4f\n", r2_per_cm2, p_per_cm2))
    cat(sprintf("Pearson r = %.3f\n", sqrt(r2_per_cm2)))
  }
}

# 4. Individual vs species-level comparison
individual_cors <- tree_weighted %>%
  left_join(
    bind_rows(
      ymf2023 %>% dplyr::select(species_id = Species.Code, CH4_flux = CH4_best.flux),
      ymf2021 %>% dplyr::select(tree_id, species_id, CH4_flux = CH4_best.flux_125cm)
    ),
    by = c("tree_id", "species_id")
  ) %>%
  filter(!is.na(CH4_flux))

if(nrow(individual_cors) > 10) {
  ind_cor <- cor.test(log10(individual_cors$mcra_area_weighted + 1), 
                      individual_cors$CH4_flux)
  lm_ind <- lm(CH4_flux ~ log10(mcra_area_weighted + 1), data = individual_cors)
  
  cat("\n4. INDIVIDUAL vs SPECIES-LEVEL CORRELATIONS:\n")
  cat("---------------------------------------------\n")
  cat(sprintf("Individual trees (n=%d):\n", nrow(individual_cors)))
  cat(sprintf("  R² = %.3f, r = %.3f, p = %.4f\n", 
              summary(lm_ind)$r.squared, ind_cor$estimate, ind_cor$p.value))
  cat(sprintf("Species-level aggregated (n=%d species):\n", nrow(correlation_data)))
  cat(sprintf("  R² = %.3f (%.1fx improvement)\n", 
              r2_median, r2_median/summary(lm_ind)$r.squared))
}

# 5. Summary paragraph values
cat("\n5. KEY VALUES FOR MANUSCRIPT TEXT:\n")
cat("-----------------------------------\n")
cat(sprintf("- Area-weighted R² (species median): %.2f (p < %.3f)\n", 
            r2_median, ifelse(p_median < 0.001, 0.001, p_median)))
cat(sprintf("- Area-weighted R² (species mean): %.2f (p < %.3f)\n", 
            r2_mean, ifelse(p_mean < 0.001, 0.001, p_mean)))
if(exists("r2_per_cm2")) {
  cat(sprintf("- Per cm² bark surface R²: %.2f (p < %.3f)\n", 
              r2_per_cm2, ifelse(p_per_cm2 < 0.001, 0.001, p_per_cm2)))
}
cat(sprintf("- Number of species in analysis: %d\n", nrow(correlation_data)))
cat(sprintf("- Sampling depth of cores: ~5 cm from pith and bark\n"))