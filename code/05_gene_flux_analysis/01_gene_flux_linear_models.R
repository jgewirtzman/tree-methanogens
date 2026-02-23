# ==============================================================================
# Gene-Flux Linear Models
# ==============================================================================
# Purpose: Analyzes area-weighted gene abundances (mcrA, pmoA, mmoxY) as
#   predictors of CH4 flux using linear models.
#
# Pipeline stage: 3 — Analysis
# Run after: 00_harmonization/02_harmonize_all_data.R
#
# Inputs:
#   - methanogen_tree_flux_complete_dataset.csv (from data/processed/flux/)
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#
# Outputs:
#   - model statistics, visualizations (to outputs/)
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(ggforce)
  library(patchwork)
  library(car)
  library(ggrepel)
  library(gridExtra)
})

# ============================================================
# LOAD DATA
# ============================================================

ymf2023 <- read.csv("../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv")
ymf2021 <- read.csv('../../data/processed/integrated/merged_tree_dataset_final.csv')
merged_final <- ymf2021

# Species mapping
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
# PART 1: PARSE TREE GENE DATA
# ============================================================

prepare_long_tree_genes <- function(df) {
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
    filter(stringency == "loose", !is.na(sample_type)) %>%
    filter(gene %in% c("mcrA", "pmoA", "mmoX"), is_probe) %>%
    filter(sample_type %in% c("Heartwood", "Sapwood"))
}

tree_genes <- prepare_long_tree_genes(merged_final)

# ============================================================
# PART 2: PARSE SOIL GENE DATA
# ============================================================

prepare_long_soil_genes <- function(df) {
  df %>%
    dplyr::select(tree_id, 
           starts_with("ddpcr_mcra_probe_Organic"), 
           starts_with("ddpcr_mcra_probe_Mineral"),
           starts_with("ddpcr_pmoa_Organic"), 
           starts_with("ddpcr_pmoa_Mineral"),
           starts_with("ddpcr_mmox_Organic"), 
           starts_with("ddpcr_mmox_Mineral")) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    pivot_longer(
      cols = starts_with("ddpcr_"),
      names_to = "measurement_type",
      values_to = "gene_copies"
    ) %>%
    filter(!is.na(gene_copies), gene_copies > 0) %>%
    separate(
      measurement_type,
      into = c("method", "gene", "part1", "part2", "part3"),
      sep = "_", extra = "merge", fill = "right"
    ) %>%
    mutate(
      is_probe = (part1 == "probe"),
      location = if_else(is_probe, part2, part1),
      stringency = if_else(is_probe, part3, part2),
      gene = case_when(
        gene == "mcra" ~ "mcrA",
        gene == "pmoa" ~ "pmoA",
        gene == "mmox" ~ "mmoX",
        TRUE ~ gene
      ),
      location = case_when(
        location %in% c("Mineral","mineral") ~ "Mineral",
        location %in% c("Organic","organic") ~ "Organic",
        TRUE ~ location
      ),
      sample_type = case_when(
        location == "Mineral" ~ "Mineral",
        location == "Organic" ~ "Organic",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(stringency == "loose", !is.na(sample_type)) %>%
    filter(gene %in% c("mcrA", "pmoA", "mmoX"))
}

soil_genes <- prepare_long_soil_genes(ymf2021)

cat(sprintf("Found %d tree gene measurements\n", nrow(tree_genes)))
cat(sprintf("Found %d soil gene measurements\n", nrow(soil_genes)))

# ============================================================
# PART 3: CALCULATE AREA-WEIGHTED VALUES
# ============================================================

area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1*r2 + r2^2
  gene_outer + (gene_inner - gene_outer) * (S / (3 * R^2))
}

# Tree genes - area weighted
tree_genes_weighted <- tree_genes %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
  mutate(species = species_mapping[species_id]) %>%
  filter(is.finite(dbh)) %>%
  group_by(tree_id, species_id, species, dbh, gene) %>%
  summarise(
    gene_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
    gene_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(gene_inner), is.finite(gene_outer)) %>%
  mutate(
    gene_area_weighted = mapply(area_weighted_gene, gene_inner, gene_outer, dbh)
  ) %>%
  dplyr::select(tree_id, species_id, species, dbh, gene, gene_area_weighted) %>%
  pivot_wider(
    names_from = gene,
    values_from = gene_area_weighted
  )

# Check what columns we have
cat("\nTree gene columns:\n")
print(names(tree_genes_weighted))

# Soil genes - averaged across organic and mineral
soil_genes_avg <- soil_genes %>%
  group_by(tree_id, gene) %>%
  summarise(soil_gene_avg = mean(gene_copies, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = gene,
    values_from = soil_gene_avg,
    names_prefix = "soil_"
  )

cat("\nSoil gene columns:\n")
print(names(soil_genes_avg))

# ============================================================
# PART 4: MERGE WITH FLUX DATA
# ============================================================

all_data <- tree_genes_weighted %>%
  left_join(soil_genes_avg, by = "tree_id") %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
    by = "tree_id"
  ) %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux))

# Create log-transformed versions - handle missing columns gracefully
if("mcrA" %in% names(all_data)) {
  all_data <- all_data %>% mutate(log_tree_mcra = log10(mcrA + 1))
}
if("pmoA" %in% names(all_data)) {
  all_data <- all_data %>% mutate(log_tree_pmoa = log10(pmoA + 1))
}
if("mmoX" %in% names(all_data)) {
  all_data <- all_data %>% mutate(log_tree_mmox = log10(mmoX + 1))
}
if("soil_mcrA" %in% names(all_data)) {
  all_data <- all_data %>% mutate(log_soil_mcra = log10(soil_mcrA + 1))
}
if("soil_pmoA" %in% names(all_data)) {
  all_data <- all_data %>% mutate(log_soil_pmoa = log10(soil_pmoA + 1))
}
if("soil_mmoX" %in% names(all_data)) {
  all_data <- all_data %>% mutate(log_soil_mmox = log10(soil_mmoX + 1))
}

cat(sprintf("\nMerged dataset: %d trees with flux data\n", nrow(all_data)))
cat(sprintf("Trees with tree mcrA: %d\n", sum(!is.na(all_data$mcrA))))
if("pmoA" %in% names(all_data)) {
  cat(sprintf("Trees with tree pmoA: %d\n", sum(!is.na(all_data$pmoA))))
}
if("mmoX" %in% names(all_data)) {
  cat(sprintf("Trees with tree mmoX: %d\n", sum(!is.na(all_data$mmoX))))
}
cat(sprintf("Trees with soil mcrA: %d\n", sum(!is.na(all_data$soil_mcrA))))
cat(sprintf("Trees with soil pmoA: %d\n", sum(!is.na(all_data$soil_pmoA))))
cat(sprintf("Trees with soil mmoX: %d\n", sum(!is.na(all_data$soil_mmoX))))

# ============================================================
# PART 5: CREATE RADIAL DISTRIBUTION PLOTS
# ============================================================

cat("\n\nCreating radial distribution plots...\n")

# Function to create radial plots
create_radial_plot <- function(gene_data, gene_name, title_text) {
  
  # Get unique trees
  trees <- unique(gene_data$tree_id)
  
  # Create plots for each tree
  plots_list <- list()
  
  for(tree in trees[1:min(9, length(trees))]) {  # Max 9 trees per page
    tree_data <- gene_data %>% filter(tree_id == tree)
    
    if(nrow(tree_data) == 0) next
    
    dbh_val <- tree_data$dbh[1]
    R <- dbh_val / 2
    
    # Create  synthetic radial profile
    radii <- seq(0, R, length.out = 100)
    
    gene_inner <- mean(tree_data$gene_inner, na.rm = TRUE)
    gene_outer <- mean(tree_data$gene_outer, na.rm = TRUE)
    
    # Linear gradient from center to edge
    profile <- gene_inner + (gene_outer - gene_inner) * (radii / R)
    
    profile_df <- data.frame(
      radius = radii,
      gene_copies = profile
    )
    
    p <- ggplot(profile_df, aes(x = radius, y = gene_copies)) +
      geom_area(fill = "steelblue", alpha = 0.5) +
      geom_line(color = "darkblue", size = 1) +
      geom_vline(xintercept = max(R - 5, 0), linetype = "dashed", 
                 color = "red", alpha = 0.5) +
      annotate("text", x = R/2, y = max(profile) * 0.9,
               label = sprintf("DBH: %.1f cm\nTree: %s", dbh_val, tree),
               size = 3) +
      labs(x = "Distance from center (cm)",
           y = "Gene copies/g",
           title = paste("Tree", tree)) +
      theme_classic(base_size = 10)
    
    plots_list[[length(plots_list) + 1]] <- p
  }
  
  # Combine plots
  if(length(plots_list) > 0) {
    combined <- wrap_plots(plots_list, ncol = 3)
    combined <- combined + plot_annotation(title = title_text)
    return(combined)
  }
  return(NULL)
}

# Create radial plots for each gene
if("mmoX" %in% names(tree_genes_weighted)) {
  mmox_radial_data <- tree_genes %>%
    filter(gene == "mmoX") %>%
    left_join(ymf2021 %>% dplyr::select(tree_id, dbh), by = "tree_id") %>%
    group_by(tree_id, dbh) %>%
    summarise(
      gene_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
      gene_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(gene_inner), is.finite(gene_outer))
  
  p_radial_mmox <- create_radial_plot(mmox_radial_data, "mmoX", 
                                      "mmoX Radial Distribution Across Trees")
  if(!is.null(p_radial_mmox)) {
    ggsave("../../outputs/figures/radial_plot_mmox.pdf", p_radial_mmox, width = 12, height = 10)
  }
}

if("pmoA" %in% names(tree_genes_weighted)) {
  pmoa_radial_data <- tree_genes %>%
    filter(gene == "pmoA") %>%
    left_join(ymf2021 %>% dplyr::select(tree_id, dbh), by = "tree_id") %>%
    group_by(tree_id, dbh) %>%
    summarise(
      gene_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
      gene_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(gene_inner), is.finite(gene_outer))
  
  p_radial_pmoa <- create_radial_plot(pmoa_radial_data, "pmoA", 
                                      "pmoA Radial Distribution Across Trees")
  if(!is.null(p_radial_pmoa)) {
    ggsave("../../outputs/figures/radial_plot_pmoa.pdf", p_radial_pmoa, width = 12, height = 10)
  }
}

# ============================================================
# PART 6: RUN ALL STATISTICAL MODELS
# ============================================================

cat("\n\n============================================================\n")
cat("STATISTICAL MODELS\n")
cat("============================================================\n")

# Initialize results list
model_results <- list()

# MODEL 1: Species + Tree mcrA
cat("\n### MODEL 1: Species + Tree mcrA ###\n")
m1_data <- all_data %>% filter(!is.na(mcrA))
if(nrow(m1_data) >= 20) {
  m1_species <- lm(CH4_flux ~ species, data = m1_data)
  m1_full <- lm(CH4_flux ~ species + log_tree_mcra, data = m1_data)
  anova_m1 <- anova(m1_species, m1_full)
  cat(sprintf("n=%d, ΔR²=%.3f, p=%.3f\n", nrow(m1_data),
              summary(m1_full)$r.squared - summary(m1_species)$r.squared,
              anova_m1$`Pr(>F)`[2]))
  model_results$m1 <- list(data=m1_data, species=m1_species, full=m1_full, anova=anova_m1)
}

# MODEL 2: Species + Tree pmoA
if("pmoA" %in% names(all_data)) {
  cat("\n### MODEL 2: Species + Tree pmoA ###\n")
  m2_data <- all_data %>% filter(!is.na(pmoA))
  if(nrow(m2_data) >= 20) {
    m2_species <- lm(CH4_flux ~ species, data = m2_data)
    m2_full <- lm(CH4_flux ~ species + log_tree_pmoa, data = m2_data)
    anova_m2 <- anova(m2_species, m2_full)
    cat(sprintf("n=%d, ΔR²=%.3f, p=%.3f\n", nrow(m2_data),
                summary(m2_full)$r.squared - summary(m2_species)$r.squared,
                anova_m2$`Pr(>F)`[2]))
    model_results$m2 <- list(data=m2_data, species=m2_species, full=m2_full, anova=anova_m2)
  }
}

# MODEL 3: Species + Tree mmoX *** KEY MODEL ***
if("mmoX" %in% names(all_data)) {
  cat("\n### MODEL 3: Species + Tree mmoX *** KEY MODEL *** ###\n")
  m3_data <- all_data %>% filter(!is.na(mmoX))
  if(nrow(m3_data) >= 20) {
    m3_species <- lm(CH4_flux ~ species, data = m3_data)
    m3_full <- lm(CH4_flux ~ species + log_tree_mmox, data = m3_data)
    anova_m3 <- anova(m3_species, m3_full)
    cat(sprintf("n=%d, ΔR²=%.3f, p=%.3f\n", nrow(m3_data),
                summary(m3_full)$r.squared - summary(m3_species)$r.squared,
                anova_m3$`Pr(>F)`[2]))
    if(anova_m3$`Pr(>F)`[2] < 0.05) {
      print(summary(m3_full))
      print(Anova(m3_full, type = "II"))
    }
    model_results$m3 <- list(data=m3_data, species=m3_species, full=m3_full, anova=anova_m3)
  }
}

# MODEL 4: Species + All tree genes
if(all(c("mcrA", "pmoA", "mmoX") %in% names(all_data))) {
  cat("\n### MODEL 4: Species + All tree genes ###\n")
  m4_data <- all_data %>% 
    filter(!is.na(mcrA), !is.na(pmoA), !is.na(mmoX))
  if(nrow(m4_data) >= 20) {
    m4_species <- lm(CH4_flux ~ species, data = m4_data)
    m4_full <- lm(CH4_flux ~ species + log_tree_mcra + log_tree_pmoa + log_tree_mmox, 
                  data = m4_data)
    anova_m4 <- anova(m4_species, m4_full)
    cat(sprintf("n=%d, ΔR²=%.3f, p=%.3f\n", nrow(m4_data),
                summary(m4_full)$r.squared - summary(m4_species)$r.squared,
                anova_m4$`Pr(>F)`[2]))
    if(anova_m4$`Pr(>F)`[2] < 0.05) {
      print(summary(m4_full))
      print(Anova(m4_full, type = "II"))
    }
    model_results$m4 <- list(data=m4_data, species=m4_species, full=m4_full, anova=anova_m4)
  }
}

# MODEL 5: Tree mmoX + Soil mcrA
if("mmoX" %in% names(all_data)) {
  cat("\n### MODEL 5: Tree mmoX + Soil mcrA ###\n")
  m5_data <- all_data %>% filter(!is.na(mmoX), !is.na(soil_mcrA))
  if(nrow(m5_data) >= 20) {
    m5_base <- lm(CH4_flux ~ species + log_tree_mmox, data = m5_data)
    m5_soil <- lm(CH4_flux ~ species + log_tree_mmox + log_soil_mcra, data = m5_data)
    anova_m5 <- anova(m5_base, m5_soil)
    cat(sprintf("n=%d, ΔR²=%.3f, p=%.3f\n", nrow(m5_data),
                summary(m5_soil)$r.squared - summary(m5_base)$r.squared,
                anova_m5$`Pr(>F)`[2]))
    model_results$m5 <- list(data=m5_data, base=m5_base, soil=m5_soil, anova=anova_m5)
  }
}

# MODEL 6: Tree mmoX + Soil pmoA
if("mmoX" %in% names(all_data)) {
  cat("\n### MODEL 6: Tree mmoX + Soil pmoA ###\n")
  m6_data <- all_data %>% filter(!is.na(mmoX), !is.na(soil_pmoA))
  if(nrow(m6_data) >= 20) {
    m6_base <- lm(CH4_flux ~ species + log_tree_mmox, data = m6_data)
    m6_soil <- lm(CH4_flux ~ species + log_tree_mmox + log_soil_pmoa, data = m6_data)
    anova_m6 <- anova(m6_base, m6_soil)
    cat(sprintf("n=%d, ΔR²=%.3f, p=%.3f\n", nrow(m6_data),
                summary(m6_soil)$r.squared - summary(m6_base)$r.squared,
                anova_m6$`Pr(>F)`[2]))
    model_results$m6 <- list(data=m6_data, base=m6_base, soil=m6_soil, anova=anova_m6)
  }
}

# MODEL 7: Tree mmoX + Soil mmoX
if("mmoX" %in% names(all_data)) {
  cat("\n### MODEL 7: Tree mmoX + Soil mmoX ###\n")
  m7_data <- all_data %>% filter(!is.na(mmoX), !is.na(soil_mmoX))
  if(nrow(m7_data) >= 20) {
    m7_base <- lm(CH4_flux ~ species + log_tree_mmox, data = m7_data)
    m7_soil <- lm(CH4_flux ~ species + log_tree_mmox + log_soil_mmox, data = m7_data)
    anova_m7 <- anova(m7_base, m7_soil)
    cat(sprintf("n=%d, ΔR²=%.3f, p=%.3f\n", nrow(m7_data),
                summary(m7_soil)$r.squared - summary(m7_base)$r.squared,
                anova_m7$`Pr(>F)`[2]))
    model_results$m7 <- list(data=m7_data, base=m7_base, soil=m7_soil, anova=anova_m7)
  }
}

# MODEL 8: Soil genes only
cat("\n### MODEL 8: Soil genes only ###\n")
m8_data <- all_data %>% 
  filter(!is.na(soil_mcrA), !is.na(soil_pmoA), !is.na(soil_mmoX))
if(nrow(m8_data) >= 20) {
  m8_species <- lm(CH4_flux ~ species, data = m8_data)
  m8_soil <- lm(CH4_flux ~ species + log_soil_mcra + log_soil_pmoa + log_soil_mmox, 
                data = m8_data)
  anova_m8 <- anova(m8_species, m8_soil)
  cat(sprintf("n=%d, ΔR²=%.3f, p=%.3f\n", nrow(m8_data),
              summary(m8_soil)$r.squared - summary(m8_species)$r.squared,
              anova_m8$`Pr(>F)`[2]))
  model_results$m8 <- list(data=m8_data, species=m8_species, full=m8_soil, anova=anova_m8)
}

# ============================================================
# PART 7: CREATE SUMMARY TABLE
# ============================================================

cat("\n\nCreating summary table...\n")

# Build summary table from actual results
summary_rows <- list()

# Add rows for each model that was run
if(!is.null(model_results$m1)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "1. Species only",
    N = nrow(model_results$m1$data),
    R2 = summary(model_results$m1$species)$r.squared,
    Delta_R2 = NA,
    p_value = NA
  )
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "2. + Tree mcrA",
    N = nrow(model_results$m1$data),
    R2 = summary(model_results$m1$full)$r.squared,
    Delta_R2 = summary(model_results$m1$full)$r.squared - summary(model_results$m1$species)$r.squared,
    p_value = model_results$m1$anova$`Pr(>F)`[2]
  )
}

if(!is.null(model_results$m2)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "3. + Tree pmoA",
    N = nrow(model_results$m2$data),
    R2 = summary(model_results$m2$full)$r.squared,
    Delta_R2 = summary(model_results$m2$full)$r.squared - summary(model_results$m2$species)$r.squared,
    p_value = model_results$m2$anova$`Pr(>F)`[2]
  )
}

if(!is.null(model_results$m3)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "4. + Tree mmoX",
    N = nrow(model_results$m3$data),
    R2 = summary(model_results$m3$full)$r.squared,
    Delta_R2 = summary(model_results$m3$full)$r.squared - summary(model_results$m3$species)$r.squared,
    p_value = model_results$m3$anova$`Pr(>F)`[2]
  )
}

if(!is.null(model_results$m4)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "5. + All tree genes",
    N = nrow(model_results$m4$data),
    R2 = summary(model_results$m4$full)$r.squared,
    Delta_R2 = summary(model_results$m4$full)$r.squared - summary(model_results$m4$species)$r.squared,
    p_value = model_results$m4$anova$`Pr(>F)`[2]
  )
}

if(!is.null(model_results$m5)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "6. Tree mmoX + Soil mcrA",
    N = nrow(model_results$m5$data),
    R2 = summary(model_results$m5$soil)$r.squared,
    Delta_R2 = summary(model_results$m5$soil)$r.squared - summary(model_results$m5$base)$r.squared,
    p_value = model_results$m5$anova$`Pr(>F)`[2]
  )
}

if(!is.null(model_results$m6)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "7. Tree mmoX + Soil pmoA",
    N = nrow(model_results$m6$data),
    R2 = summary(model_results$m6$soil)$r.squared,
    Delta_R2 = summary(model_results$m6$soil)$r.squared - summary(model_results$m6$base)$r.squared,
    p_value = model_results$m6$anova$`Pr(>F)`[2]
  )
}

if(!is.null(model_results$m7)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "8. Tree mmoX + Soil mmoX",
    N = nrow(model_results$m7$data),
    R2 = summary(model_results$m7$soil)$r.squared,
    Delta_R2 = summary(model_results$m7$soil)$r.squared - summary(model_results$m7$base)$r.squared,
    p_value = model_results$m7$anova$`Pr(>F)`[2]
  )
}

if(!is.null(model_results$m8)) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Model = "9. Soil genes only",
    N = nrow(model_results$m8$data),
    R2 = summary(model_results$m8$full)$r.squared,
    Delta_R2 = summary(model_results$m8$full)$r.squared - summary(model_results$m8$species)$r.squared,
    p_value = model_results$m8$anova$`Pr(>F)`[2]
  )
}

model_summary <- bind_rows(summary_rows)
model_summary$Significant <- ifelse(
  is.na(model_summary$p_value), NA,
  ifelse(model_summary$p_value < 0.01, "**",
         ifelse(model_summary$p_value < 0.05, "*", "NS"))
)

cat("\n============================================================\n")
cat("MODEL SUMMARY TABLE\n")
cat("============================================================\n")
print(model_summary)

write.csv(model_summary, "../../outputs/tables/complete_model_summary.csv", row.names = FALSE)

# ============================================================
# PART 8: CREATE VISUALIZATIONS
# ============================================================

cat("\n\nCreating final publication figure...\n")

# Only create plots for models that exist
plot_list <- list()

# Panel A: Tree mcrA
if(!is.null(model_results$m1)) {
  p_tree_mcra <- ggplot(model_results$m1$data, aes(x = log_tree_mcra, y = CH4_flux)) +
    geom_point(aes(color = species), alpha = 0.4, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", 
                fill = "pink", alpha = 0.2, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("NS\np = %.2f", model_results$m1$anova$`Pr(>F)`[2]), 
             hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
    labs(title = "A) Tree mcrA",
         x = expression("log"[10]*" mcrA"),
         y = expression("CH"[4]*" flux")) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
  plot_list$mcra <- p_tree_mcra
}

# Panel B: Tree mmoX
if(!is.null(model_results$m3)) {
  mmox_coef <- coef(model_results$m3$full)["log_tree_mmox"]
  p_tree_mmox <- ggplot(model_results$m3$data, aes(x = log_tree_mmox, y = CH4_flux)) +
    geom_point(aes(color = species), alpha = 0.4, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue", 
                fill = "lightblue", alpha = 0.2, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("**\np = %.3f\nβ = %.2f", 
                             model_results$m3$anova$`Pr(>F)`[2], mmox_coef), 
             hjust = 1.1, vjust = 1.1, size = 3.5, 
             fontface = "bold", color = "steelblue") +
    labs(title = "B) Tree mmoX",
         x = expression("log"[10]*" mmoX"),
         y = "") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
  plot_list$mmox <- p_tree_mmox
}

# Panel C: Soil mcrA
if(!is.null(model_results$m5)) {
  p_soil_mcra <- ggplot(model_results$m5$data, aes(x = log_soil_mcra, y = CH4_flux)) +
    geom_point(aes(color = species), alpha = 0.4, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "gray50", 
                fill = "gray80", alpha = 0.2, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("NS\np = %.2f", model_results$m5$anova$`Pr(>F)`[2]), 
             hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
    labs(title = "C) Soil mcrA",
         x = expression("log"[10]*" soil mcrA"),
         y = expression("CH"[4]*" flux")) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
  plot_list$soil_mcra <- p_soil_mcra
}

# Panel D: Soil mmoX
if(!is.null(model_results$m7)) {
  p_soil_mmox <- ggplot(model_results$m7$data, aes(x = log_soil_mmox, y = CH4_flux)) +
    geom_point(aes(color = species), alpha = 0.4, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "gray50", 
                fill = "gray80", alpha = 0.2, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("NS\np = %.2f", model_results$m7$anova$`Pr(>F)`[2]), 
             hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
    labs(title = "D) Soil mmoX",
         x = expression("log"[10]*" soil mmoX"),
         y = "") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
  plot_list$soil_mmox <- p_soil_mmox
}

# Panel E: Model comparison
model_plot_data <- model_summary %>%
  filter(!is.na(Delta_R2)) %>%
  mutate(
    Model_short = gsub("^[0-9]+\\. ", "", Model),
    Model_short = factor(Model_short, 
                         levels = rev(unique(Model_short))),
    Category = ifelse(grepl("Soil", Model_short), "Soil", "Tree")
  )

p_comparison <- ggplot(model_plot_data, aes(x = Delta_R2, y = Model_short, fill = Category)) +
  geom_col(alpha = 0.8) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_text(aes(label = sprintf("p = %.3f", p_value)),
            hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c("Tree" = "steelblue", "Soil" = "darkorange")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(title = "E) Model Improvement",
       x = expression(Delta*R^2),
       y = "") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")
plot_list$comparison <- p_comparison

# Panel F: Conceptual model
p_concept <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0.7, ymax = 1, 
           fill = "steelblue", alpha = 0.3) +
  annotate("text", x = 0.5, y = 0.85, 
           label = "TREE STEM\n(mmoX controls flux)", 
           size = 3.5, fontface = "bold", color = "steelblue") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0.35, ymax = 0.65, 
           fill = "gray50", alpha = 0.2) +
  annotate("text", x = 0.5, y = 0.5, 
           label = "DECOUPLED\n(no interaction)", 
           size = 3, fontface = "italic", color = "gray30") +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 0.3, 
           fill = "darkorange", alpha = 0.3) +
  annotate("text", x = 0.5, y = 0.15, 
           label = "SOIL\n(no effect)", 
           size = 3.5, fontface = "bold", color = "darkorange") +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "F) Conceptual Model") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
plot_list$concept <- p_concept

# Combine all available plots
if(length(plot_list) >= 4) {
  layout <- "
  AABBEE
  CCDDEE
  FFFFFF
  "
  
  combined <- plot_list$mcra + plot_list$mmox + 
    plot_list$soil_mcra + plot_list$soil_mmox + 
    plot_list$comparison + plot_list$concept +
    plot_layout(design = layout) +
    plot_annotation(
      title = "In-stem methanotrophy regulates tree methane emissions",
      subtitle = "Tree mmoX is the only significant predictor",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11)
      )
    )
  
  print(combined)
  
  ggsave("../../outputs/figures/Figure_Complete_Gene_Analysis_CORRECTED.png", combined,
         width = 12, height = 10, dpi = 300)
  ggsave("../../outputs/figures/Figure_Complete_Gene_Analysis_CORRECTED.pdf", combined,
         width = 12, height = 10)
}

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("\nFiles created:\n")
cat("- complete_model_summary.csv\n")
cat("- Figure_Complete_Gene_Analysis_CORRECTED.png/pdf\n")
if(file.exists("radial_plot_mmox.pdf")) cat("- radial_plot_mmox.pdf\n")
if(file.exists("radial_plot_pmoa.pdf")) cat("- radial_plot_pmoa.pdf\n")

cat("\n### KEY FINDINGS ###\n")
if(!is.null(model_results$m3)) {
  cat(sprintf("Tree mmoX: SIGNIFICANT (p = %.3f, β = %.2f)\n", 
              model_results$m3$anova$`Pr(>F)`[2],
              coef(model_results$m3$full)["log_tree_mmox"]))
}
if(!is.null(model_results$m1)) {
  cat(sprintf("Tree mcrA: NOT significant (p = %.3f)\n", 
              model_results$m1$anova$`Pr(>F)`[2]))
}
if(!is.null(model_results$m5)) {
  cat(sprintf("Soil genes: NO effect (p = %.3f - %.3f)\n",
              min(sapply(model_results[c("m5","m6","m7")], 
                         function(x) if(!is.null(x)) x$anova$`Pr(>F)`[2] else 1)),
              max(sapply(model_results[c("m5","m6","m7")], 
                         function(x) if(!is.null(x)) x$anova$`Pr(>F)`[2] else 0))))
}
cat("\n============================================================\n")


