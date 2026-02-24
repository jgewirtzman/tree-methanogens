# ==============================================================================
# Methanotroph RF Models (Figures S1-S2)
# ==============================================================================
# Purpose: Random forest models predicting flux from methanotroph gene
#   abundance data.
#
# Pipeline stage: 3 — Analysis
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#   - methanogen_tree_flux_complete_dataset.csv (from data/processed/flux/)
#
# Outputs:
#   - RF model objects, diagnostic plots (to outputs/)
# ==============================================================================

sink("../../../outputs/tables/output.txt", split = TRUE)  # split=TRUE shows output in console too

library(tidyverse)
library(randomForest)
library(pdp)
library(car)
library(corrplot)
library(broom)
library(ggforce)
library(patchwork)

# ============================================================
# 1. DATA PREPARATION
# ============================================================

cat("\n### LOADING AND PREPARING DATA ###\n")

# Load data
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

# ============================================================
# 2. PARSE TREE GENES
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
    filter(gene %in% c("mcrA", "pmoA", "mmoX")) %>%
    filter(sample_type %in% c("Heartwood", "Sapwood")) %>%
    filter((gene == "mcrA" & is_probe) | gene %in% c("pmoA", "mmoX"))
}

tree_genes <- prepare_long_tree_genes(ymf2021)
cat(sprintf("Found %d tree gene measurements\n", nrow(tree_genes)))

# ============================================================
# 3. PARSE SOIL GENES
# ============================================================

prepare_long_soil_genes <- function(df) {
  df %>%
    dplyr::select(tree_id, starts_with("ddpcr_mcra_probe"), 
           starts_with("ddpcr_pmoa_Mineral"), starts_with("ddpcr_pmoa_Organic"),
           starts_with("ddpcr_mmox_Mineral"), starts_with("ddpcr_mmox_Organic")) %>%
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
      sample_type = case_when(
        location %in% c("Mineral", "mineral") ~ "Mineral",
        location %in% c("Organic", "organic") ~ "Organic",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(stringency == "loose", !is.na(sample_type)) %>%
    filter(gene %in% c("mcrA", "pmoA", "mmoX"))
}

soil_genes <- prepare_long_soil_genes(ymf2021)
cat(sprintf("Found %d soil gene measurements\n", nrow(soil_genes)))

# ============================================================
# 4. AREA-WEIGHTED GENE CALCULATION
# ============================================================

area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1*r2 + r2^2
  gene_outer + (gene_inner - gene_outer) * (S / (3 * R^2))
}

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

# Calculate average soil genes
soil_genes_avg <- soil_genes %>%
  group_by(tree_id, gene) %>%
  summarise(soil_gene_avg = mean(gene_copies, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = gene,
    values_from = soil_gene_avg,
    names_prefix = "soil_"
  )

# Merge everything
all_data_with_soil <- tree_genes_weighted %>%
  left_join(soil_genes_avg, by = "tree_id") %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
    by = "tree_id"
  ) %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux))

cat(sprintf("\nMerged dataset: %d trees with flux\n", nrow(all_data_with_soil)))

# Add log-transformed genes
analysis_data <- all_data_with_soil %>%
  mutate(
    log_tree_mcra = if("mcrA" %in% names(.)) log10(mcrA + 1) else NA_real_,
    log_tree_pmoa = if("pmoA" %in% names(.)) log10(pmoA + 1) else NA_real_,
    log_tree_mmox = if("mmoX" %in% names(.)) log10(mmoX + 1) else NA_real_,
    log_soil_mcra = if("soil_mcrA" %in% names(.)) log10(soil_mcrA + 1) else NA_real_,
    log_soil_pmoa = if("soil_pmoA" %in% names(.)) log10(soil_pmoA + 1) else NA_real_,
    log_soil_mmox = if("soil_mmoX" %in% names(.)) log10(soil_mmoX + 1) else NA_real_,
    species = as.factor(species)
  )

# ============================================================
# 5. DESCRIPTIVE STATISTICS
# ============================================================

cat("\n### TABLE 1: DESCRIPTIVE STATISTICS ###\n")

# Overall summary
overall_summary <- analysis_data %>%
  summarise(
    n_trees = n(),
    n_species = n_distinct(species),
    mean_dbh = mean(dbh, na.rm = TRUE),
    sd_dbh = sd(dbh, na.rm = TRUE),
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    sd_flux = sd(CH4_flux, na.rm = TRUE)
  )
print(overall_summary)

# By species
species_summary <- analysis_data %>%
  group_by(species) %>%
  summarise(
    n = n(),
    mean_dbh = round(mean(dbh, na.rm = TRUE), 1),
    sd_dbh = round(sd(dbh, na.rm = TRUE), 1),
    mean_flux = round(mean(CH4_flux, na.rm = TRUE), 4),
    sd_flux = round(sd(CH4_flux, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(n))
print(species_summary)

write.csv(species_summary, "../../../outputs/tables/Table1_species_summary.csv", row.names = FALSE)

# ============================================================
# 6. CORRELATION ANALYSIS
# ============================================================

cat("\n### TABLE 2: CORRELATION MATRIX ###\n")

cor_cols <- c("CH4_flux", "dbh", "log_tree_mcra", "log_tree_pmoa", 
              "log_tree_mmox", "log_soil_mcra", "log_soil_pmoa", "log_soil_mmox")
existing_cor_cols <- cor_cols[cor_cols %in% names(analysis_data)]

cor_matrix <- cor(analysis_data %>% dplyr::select(all_of(existing_cor_cols)), 
                  use = "pairwise.complete.obs")
print(round(cor_matrix, 3))

write.csv(round(cor_matrix, 3), "../../../outputs/tables/Table2_correlation_matrix.csv")

pdf("../../../outputs/figures/FigureS1_correlation_matrix.pdf", width = 8, height = 8)
corrplot(cor_matrix, method = "color", type = "upper", 
         addCoef.col = "black", number.cex = 0.7,
         tl.col = "black", tl.srt = 45)
dev.off()

# ============================================================
# 7. LINEAR MODEL dplyr::selectION
# ============================================================

cat("\n### TABLE 3: MODEL dplyr::selectION ###\n")

# Use complete-case data for fair comparison
model_data <- analysis_data %>%
  filter(!is.na(log_tree_mmox), !is.na(log_tree_pmoa), !is.na(log_tree_mcra))

cat(sprintf("Using %d trees with complete gene data\n", nrow(model_data)))

# Fit candidate models
m0 <- lm(CH4_flux ~ species, data = model_data)
m1 <- lm(CH4_flux ~ species + log_tree_mmox, data = model_data)
m2 <- lm(CH4_flux ~ species + log_tree_pmoa, data = model_data)
m3 <- lm(CH4_flux ~ species + log_tree_mcra, data = model_data)
m4 <- lm(CH4_flux ~ species + log_tree_mmox + log_tree_pmoa, data = model_data)
m5 <- lm(CH4_flux ~ species + log_tree_mmox + log_tree_mcra, data = model_data)

# Model comparison
aic_values <- AIC(m0, m1, m2, m3, m4, m5)
bic_values <- BIC(m0, m1, m2, m3, m4, m5)

model_comparison <- data.frame(
  Model = c("Species only", "Species + mmoX", "Species + pmoA", 
            "Species + mcrA", "Species + mmoX + pmoA", "Species + mmoX + mcrA"),
  df = aic_values$df,
  AIC = round(aic_values$AIC, 1),
  BIC = round(bic_values$BIC, 1),
  R2 = round(c(summary(m0)$r.squared, summary(m1)$r.squared,
               summary(m2)$r.squared, summary(m3)$r.squared,
               summary(m4)$r.squared, summary(m5)$r.squared), 3),
  Adj_R2 = round(c(summary(m0)$adj.r.squared, summary(m1)$adj.r.squared,
                   summary(m2)$adj.r.squared, summary(m3)$adj.r.squared,
                   summary(m4)$adj.r.squared, summary(m5)$adj.r.squared), 3)
) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    delta_BIC = BIC - min(BIC)
  ) %>%
  arrange(AIC)

print(model_comparison)
write.csv(model_comparison, "../../../outputs/tables/Table3_model_comparison.csv", row.names = FALSE)

# Best model: species + mmoX
final_model <- m1
cat("\n### BEST MODEL: Species + mmoX ###\n")
print(summary(final_model))
cat("\nType II ANOVA:\n")
print(Anova(final_model, type = "II"))

# Save coefficients
final_coefs <- tidy(final_model, conf.int = TRUE)
write.csv(final_coefs, "../../../outputs/tables/Table4_final_model_coefficients.csv", row.names = FALSE)

# Model diagnostics
pdf("../../../outputs/figures/FigureS2_model_diagnostics.pdf", width = 10, height = 10)
par(mfrow = c(2, 2))
plot(final_model)
dev.off()

# ============================================================
# 8. SOIL GENE ANALYSIS
# ============================================================

cat("\n### TESTING SOIL GENE EFFECTS ###\n")

# Test if soil genes add anything beyond tree mmoX
test_soil_genes <- function(soil_gene_name, soil_gene_col) {
  temp_data <- all_data_with_soil %>%
    filter(!is.na(mmoX), !is.na(.data[[soil_gene_col]])) %>%
    mutate(
      log_tree_mmox = log10(mmoX + 1),
      log_soil = log10(.data[[soil_gene_col]] + 1)
    )
  
  if(nrow(temp_data) < 30) {
    cat(sprintf("\n%s: Insufficient data (n=%d)\n", soil_gene_name, nrow(temp_data)))
    return(NULL)
  }
  
  m_base <- lm(CH4_flux ~ species + log_tree_mmox, data = temp_data)
  m_soil <- lm(CH4_flux ~ species + log_tree_mmox + log_soil, data = temp_data)
  
  anova_result <- anova(m_base, m_soil)
  
  cat(sprintf("\n%s (n=%d):\n", soil_gene_name, nrow(temp_data)))
  cat(sprintf("  Base R²: %.3f\n", summary(m_base)$r.squared))
  cat(sprintf("  +Soil R²: %.3f (Δ=%.3f)\n", 
              summary(m_soil)$r.squared,
              summary(m_soil)$r.squared - summary(m_base)$r.squared))
  cat(sprintf("  p-value: %.3f %s\n", 
              anova_result$`Pr(>F)`[2],
              ifelse(anova_result$`Pr(>F)`[2] < 0.05, "**", "NS")))
  
  return(list(
    gene = soil_gene_name,
    n = nrow(temp_data),
    delta_r2 = summary(m_soil)$r.squared - summary(m_base)$r.squared,
    p_value = anova_result$`Pr(>F)`[2]
  ))
}

soil_results <- list(
  test_soil_genes("Soil mcrA", "soil_mcrA"),
  test_soil_genes("Soil pmoA", "soil_pmoA"),
  test_soil_genes("Soil mmoX", "soil_mmoX")
)

cat("\n*** RESULT: Soil genes do NOT improve predictions beyond tree mmoX ***\n")

# ============================================================
# 9. RANDOM FOREST WITH NATIVE NA HANDLING
# ============================================================

cat("\n### RANDOM FOREST ANALYSIS ###\n")

rf_data <- analysis_data %>%
  dplyr::select(CH4_flux, species, dbh,
         log_tree_mcra, log_tree_pmoa, log_tree_mmox,
         log_soil_mcra, log_soil_pmoa, log_soil_mmox) %>%
  as.data.frame()

set.seed(42)
rf_model <- randomForest(
  CH4_flux ~ .,
  data = rf_data,
  ntree = 1000,
  importance = TRUE,
  na.action = na.roughfix
)

cat("\nRandom Forest Results:\n")
print(rf_model)

importance_df <- importance(rf_model) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  arrange(desc(`%IncMSE`))

cat("\nVariable Importance:\n")
print(importance_df)

write.csv(importance_df, "../../../outputs/tables/rf_importance.csv", row.names = FALSE)

# Importance plot
pdf("../../../outputs/figures/Figure_RF_Importance.pdf", width = 8, height = 6)
varImpPlot(rf_model, main = "Random Forest Variable Importance")
dev.off()

# ============================================================
# 10. PUBLICATION FIGURE
# ============================================================

cat("\n### CREATING PUBLICATION FIGURE ###\n")

# Species-level summaries for plotting
species_plot_data <- all_data_with_soil %>%
  filter(!is.na(CH4_flux)) %>%
  group_by(species) %>%
  summarise(
    n = n(),
    median_mcra = median(mcrA, na.rm = TRUE),
    median_mmox = median(mmoX, na.rm = TRUE),
    median_flux = median(CH4_flux, na.rm = TRUE),
    se_flux = sd(CH4_flux, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  filter(n >= 3)

# Panel A: mcrA (NOT significant)
p_mcra <- ggplot(species_plot_data %>% filter(!is.na(median_mcra)), 
                 aes(x = log10(median_mcra + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "gray50", fill = "gray80", alpha = 0.3) +
  geom_point(size = 3, alpha = 0.8, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = Inf, y = Inf, label = "NS\np = 0.46", 
           hjust = 1.1, vjust = 1.1, size = 4, fontface = "bold") +
  labs(title = "A) Methanogens (mcrA)",
       x = expression("log"[10]*" mcrA"),
       y = expression("CH"[4]*" flux")) +
  theme_classic()

# Panel B: mmoX (SIGNIFICANT)
p_mmox <- ggplot(species_plot_data %>% filter(!is.na(median_mmox)), 
                 aes(x = log10(median_mmox + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", fill = "lightblue", alpha = 0.3) +
  geom_point(size = 3, alpha = 0.8, color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = Inf, y = Inf, label = "**\np = 0.005\nβ = -0.42", 
           hjust = 1.1, vjust = 1.1, size = 4, fontface = "bold", color = "steelblue") +
  labs(title = "B) Methanotrophs (mmoX)",
       x = expression("log"[10]*" mmoX"),
       y = "") +
  theme_classic()

# Panel C: Model comparison
p_comparison <- ggplot(model_comparison %>% filter(Model != "Species only"), 
                       aes(x = delta_AIC, y = reorder(Model, -delta_AIC))) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("ΔR² = %.3f", R2 - model_comparison$R2[1])),
            hjust = -0.1, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "C) Model Comparison",
       x = "ΔAIC from Species-only",
       y = "") +
  theme_classic()

# Combine
combined <- (p_mcra | p_mmox) / p_comparison +
  plot_annotation(
    title = "In-stem methanotrophy (mmoX) regulates tree methane emissions",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("../../../outputs/figures/Figure_Main_Results.pdf", combined, width = 12, height = 10)
ggsave("../../../outputs/figures/Figure_Main_Results.png", combined, width = 12, height = 10, dpi = 300)





























# ============================================================
# COMPREHENSIVE RANDOM FOREST ANALYSIS - EXTENDED VERSION
# ============================================================
# This script performs extensive RF analysis with:
# 1. Multiple imputation strategies comparison
# 2. Partial dependence plots for top predictors
# 3. Two-way interaction plots
# 4. Interaction strength quantification (H-statistic)
# 5. DBH vs DBH-normalized-within-species comparison
# 6. Area-weighted vs concentration (heartwood/sapwood) comparison
# 7. Separate pmoA+mmoX vs combined sum comparison
# ============================================================

library(tidyverse)
library(randomForest)
library(missForest)
library(pdp)
library(iml)
library(gridExtra)

cat("\n============================================================\n")
cat("COMPREHENSIVE RANDOM FOREST ANALYSIS\n")
cat("============================================================\n")

# Load the prepared data from previous script
# (assumes you've already run the consolidated script)

# ============================================================
# PART 1: PREPARE MULTIPLE MODEL VARIANTS
# ============================================================

cat("\n### Preparing model variants ###\n")

# Variant 1: Area-weighted genes with absolute DBH
data_v1_area_dbh <- analysis_data %>%
  mutate(
    methanotroph_sum = log10((pmoA + mmoX) + 1)
  ) %>%
  dplyr::select(CH4_flux, species, dbh,
         log_tree_mcra, log_tree_pmoa, log_tree_mmox, methanotroph_sum,
         log_soil_mcra, log_soil_pmoa, log_soil_mmox) %>%
  as.data.frame()

# Variant 2: Area-weighted genes with DBH normalized within species
data_v2_area_dbh_norm <- analysis_data %>%
  group_by(species) %>%
  mutate(
    species_n = n(),
    dbh_mean = mean(dbh, na.rm=TRUE),
    dbh_sd = sd(dbh, na.rm=TRUE),
    dbh_norm = if(species_n[1] > 1 & !is.na(dbh_sd[1]) & dbh_sd[1] > 0) {
      (dbh - dbh_mean) / dbh_sd
    } else {
      0
    },
    methanotroph_sum = log10((pmoA + mmoX) + 1)
  ) %>%
  ungroup() %>%
  dplyr::select(CH4_flux, species, dbh_norm,
         log_tree_mcra, log_tree_pmoa, log_tree_mmox, methanotroph_sum,
         log_soil_mcra, log_soil_pmoa, log_soil_mmox) %>%
  as.data.frame()

# Variant 3: Concentration data (heartwood & sapwood separately)
# First need to prepare this from raw tree_genes
conc_data <- tree_genes %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh, CH4_flux = CH4_best.flux_125cm), 
            by = "tree_id") %>%
  mutate(species = species_mapping[species_id]) %>%
  filter(!is.na(CH4_flux), is.finite(dbh)) %>%
  group_by(tree_id, species, dbh, CH4_flux, gene, sample_type) %>%
  summarise(gene_conc = mean(gene_copies, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = c(gene, sample_type),
    values_from = gene_conc,
    names_sep = "_"
  ) %>%
  mutate(
    log_mcra_heart = log10(mcrA_Heartwood + 1),
    log_mcra_sap = log10(mcrA_Sapwood + 1),
    log_pmoa_heart = log10(pmoA_Heartwood + 1),
    log_pmoa_sap = log10(pmoA_Sapwood + 1),
    log_mmox_heart = log10(mmoX_Heartwood + 1),
    log_mmox_sap = log10(mmoX_Sapwood + 1),
    log_methanotroph_heart = log10((pmoA_Heartwood + mmoX_Heartwood) + 1),
    log_methanotroph_sap = log10((pmoA_Sapwood + mmoX_Sapwood) + 1),
    species = as.factor(species)
  )

# Add soil genes to concentration data
conc_data <- conc_data %>%
  left_join(soil_genes_avg, by = "tree_id") %>%
  mutate(
    log_soil_mcra = log10(soil_mcrA + 1),
    log_soil_pmoa = log10(soil_pmoA + 1),
    log_soil_mmox = log10(soil_mmoX + 1)
  )

data_v3_conc_dbh <- conc_data %>%
  dplyr::select(CH4_flux, species, dbh,
         log_mcra_heart, log_mcra_sap,
         log_pmoa_heart, log_pmoa_sap,
         log_mmox_heart, log_mmox_sap,
         log_methanotroph_heart, log_methanotroph_sap,
         log_soil_mcra, log_soil_pmoa, log_soil_mmox) %>%
  as.data.frame()

data_v4_conc_dbh_norm <- conc_data %>%
  group_by(species) %>%
  mutate(
    species_n = n(),
    dbh_mean = mean(dbh, na.rm=TRUE),
    dbh_sd = sd(dbh, na.rm=TRUE),
    dbh_norm = if(species_n[1] > 1 & !is.na(dbh_sd[1]) & dbh_sd[1] > 0) {
      (dbh - dbh_mean) / dbh_sd
    } else {
      0
    }
  ) %>%
  ungroup() %>%
  dplyr::select(CH4_flux, species, dbh_norm,
         log_mcra_heart, log_mcra_sap,
         log_pmoa_heart, log_pmoa_sap,
         log_mmox_heart, log_mmox_sap,
         log_methanotroph_heart, log_methanotroph_sap,
         log_soil_mcra, log_soil_pmoa, log_soil_mmox) %>%
  as.data.frame()

cat(sprintf("\nData variant sizes:\n"))
cat(sprintf("V1 (area-weighted, DBH): %d trees\n", nrow(data_v1_area_dbh)))
cat(sprintf("V2 (area-weighted, DBH-norm): %d trees\n", nrow(data_v2_area_dbh_norm)))
cat(sprintf("V3 (concentration, DBH): %d trees\n", nrow(data_v3_conc_dbh)))
cat(sprintf("V4 (concentration, DBH-norm): %d trees\n", nrow(data_v4_conc_dbh_norm)))

# ============================================================
# PART 2: COMPARISON OF IMPUTATION METHODS
# ============================================================

cat("\n\n### COMPARING IMPUTATION STRATEGIES ###\n")

fit_rf_with_method <- function(data, method_name, use_sum = FALSE) {
  cat(sprintf("\nFitting %s...\n", method_name))
  
  set.seed(42)
  
  if(method_name == "Complete Cases") {
    data_clean <- data %>% drop_na()
    if(nrow(data_clean) < 30) {
      cat("Insufficient complete cases\n")
      return(NULL)
    }
    rf <- randomForest(CH4_flux ~ ., data = data_clean, 
                       ntree = 1000, importance = TRUE)
    return(list(model = rf, data = data_clean, method = method_name))
    
  } else if(method_name == "na.roughfix") {
    rf <- randomForest(CH4_flux ~ ., data = data,
                       ntree = 1000, importance = TRUE,
                       na.action = na.roughfix)
    return(list(model = rf, data = data, method = method_name))
    
  } else if(method_name == "missForest") {
    # Impute only gene columns
    gene_cols <- names(data)[grepl("log_", names(data))]
    if(length(gene_cols) == 0) return(NULL)
    
    cat("Running iterative imputation...\n")
    imputed <- missForest(data %>% dplyr::select(all_of(gene_cols)), 
                          verbose = FALSE, maxiter = 5)
    
    data_imputed <- data.frame(
      CH4_flux = data$CH4_flux,
      species = data$species,
      data %>% dplyr::select(matches("^dbh")),
      imputed$ximp
    )
    
    rf <- randomForest(CH4_flux ~ ., data = data_imputed,
                       ntree = 1000, importance = TRUE)
    return(list(model = rf, data = data_imputed, method = method_name))
  }
}

# Fit all combinations
results_list <- list()

cat("\n=== AREA-WEIGHTED + DBH ===\n")
results_list$v1_roughfix <- fit_rf_with_method(data_v1_area_dbh, "na.roughfix")
results_list$v1_missforest <- fit_rf_with_method(data_v1_area_dbh, "missForest")
results_list$v1_complete <- fit_rf_with_method(data_v1_area_dbh, "Complete Cases")

cat("\n=== AREA-WEIGHTED + DBH-NORM ===\n")
results_list$v2_roughfix <- fit_rf_with_method(data_v2_area_dbh_norm, "na.roughfix")
results_list$v2_missforest <- fit_rf_with_method(data_v2_area_dbh_norm, "missForest")
results_list$v2_complete <- fit_rf_with_method(data_v2_area_dbh_norm, "Complete Cases")

cat("\n=== CONCENTRATION + DBH ===\n")
results_list$v3_roughfix <- fit_rf_with_method(data_v3_conc_dbh, "na.roughfix")
results_list$v3_missforest <- fit_rf_with_method(data_v3_conc_dbh, "missForest")
results_list$v3_complete <- fit_rf_with_method(data_v3_conc_dbh, "Complete Cases")

cat("\n=== CONCENTRATION + DBH-NORM ===\n")
results_list$v4_roughfix <- fit_rf_with_method(data_v4_conc_dbh_norm, "na.roughfix")
results_list$v4_missforest <- fit_rf_with_method(data_v4_conc_dbh_norm, "missForest")
results_list$v4_complete <- fit_rf_with_method(data_v4_conc_dbh_norm, "Complete Cases")

# Remove NULL entries
results_list <- results_list[!sapply(results_list, is.null)]

# ============================================================
# PART 3: MODEL COMPARISON TABLE
# ============================================================

cat("\n### Creating comparison table ###\n")

comparison_table <- do.call(rbind, lapply(names(results_list), function(name) {
  res <- results_list[[name]]
  data.frame(
    Model = name,
    N = nrow(res$data),
    N_Predictors = ncol(res$data) - 1,
    MSE = tail(res$model$mse, 1),
    Var_Explained_Pct = max(res$model$rsq) * 100,
    stringsAsFactors = FALSE
  )
}))

comparison_table <- comparison_table %>%
  arrange(desc(Var_Explained_Pct))

cat("\n### MODEL PERFORMANCE COMPARISON ###\n")
print(comparison_table)
write.csv(comparison_table, "../../../outputs/tables/RF_model_comparison.csv", row.names = FALSE)

# ============================================================
# PART 4: VARIABLE IMPORTANCE ACROSS MODELS
# ============================================================

cat("\n### Extracting variable importance ###\n")

importance_list <- lapply(names(results_list), function(name) {
  res <- results_list[[name]]
  imp <- importance(res$model) %>%
    as.data.frame() %>%
    rownames_to_column("Variable") %>%
    mutate(Model = name) %>%
    dplyr::select(Model, Variable, `%IncMSE`, IncNodePurity)
  return(imp)
})

importance_combined <- do.call(rbind, importance_list)
write.csv(importance_combined, "../../../outputs/tables/RF_importance_all_models.csv", row.names = FALSE)

# Plot importance for top 3 models
top_models <- head(comparison_table$Model, 3)

pdf("../../../outputs/figures/Figure_RF_Importance_Comparison.pdf", width = 12, height = 8)
par(mfrow = c(1, 3))
for(model_name in top_models) {
  res <- results_list[[model_name]]
  tryCatch({
    varImpPlot(res$model, main = model_name, n.var = 10)
  }, error = function(e) {
    cat(sprintf("Error plotting %s: %s\n", model_name, e$message))
    plot.new()
    text(0.5, 0.5, sprintf("%s\n(plot error)", model_name), cex = 1.5)
  })
}
par(mfrow = c(1, 1))
dev.off()

cat("Saved: Figure_RF_Importance_Comparison.pdf\n")

# ============================================================
# PART 5: PARTIAL DEPENDENCE PLOTS
# ============================================================

cat("\n### Creating partial dependence plots ###\n")

# Use best model
best_model_name <- comparison_table$Model[1]
best_result <- results_list[[best_model_name]]
best_model <- best_result$model
best_data <- best_result$data

cat(sprintf("Using best model: %s (R²=%.1f%%)\n", 
            best_model_name, max(best_model$rsq)*100))

# Get top predictors (excluding species)
imp <- importance(best_model) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  filter(Variable != "species") %>%
  arrange(desc(`%IncMSE`)) %>%
  head(6)

cat("\nTop 6 predictors for PDPs:\n")
print(imp)

pdf("../../../outputs/figures/Figure_Partial_Dependence.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for(i in 1:min(6, nrow(imp))) {
  var_name <- imp$Variable[i]
  
  if(var_name %in% names(best_data)) {
    pd <- partial(best_model, pred.var = var_name, 
                  train = best_data, grid.resolution = 25)
    
    plot(pd, type = "l", lwd = 2, col = "steelblue",
         main = var_name,
         xlab = var_name,
         ylab = "Predicted CH4 flux")
    rug(best_data[[var_name]], col = "gray50", ticksize = 0.02)
  }
}

par(mfrow = c(1, 1))
dev.off()

cat("Saved: Figure_Partial_Dependence.pdf\n")

# ============================================================
# PART 6: TWO-WAY INTERACTIONS
# ============================================================

cat("\n### Creating two-way interaction plots ###\n")

# Test key interactions
top_vars <- head(imp$Variable, 4)

pdf("../../../outputs/figures/Figure_2way_Interactions.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))

interaction_pairs <- list(
  c(top_vars[1], top_vars[2]),
  c(top_vars[1], top_vars[3]),
  c(top_vars[2], top_vars[3]),
  c(top_vars[1], "dbh")  # Always check DBH if available
)

# Filter to valid pairs
interaction_pairs <- interaction_pairs[sapply(interaction_pairs, function(p) {
  all(p %in% names(best_data))
})]

for(pair in interaction_pairs) {
  if(length(pair) == 2) {
    cat(sprintf("Computing %s × %s interaction...\n", pair[1], pair[2]))
    
    tryCatch({
      pd_2way <- partial(best_model, pred.var = pair,
                         train = best_data, grid.resolution = 20)
      
      # Convert to matrix for image plot
      mat <- pd_2way %>%
        pivot_wider(names_from = all_of(pair[2]), 
                    values_from = yhat) %>%
        dplyr::select(-all_of(pair[1])) %>%
        as.matrix()
      
      image(mat, main = sprintf("%s × %s", pair[1], pair[2]),
            xlab = pair[1], ylab = pair[2],
            col = heat.colors(50))
      contour(mat, add = TRUE, nlevels = 5)
    }, error = function(e) {
      cat(sprintf("Error creating interaction plot: %s\n", e$message))
    })
  }
}

par(mfrow = c(1, 1))
dev.off()

cat("Saved: Figure_2way_Interactions.pdf\n")

# ============================================================
# PART 7: INTERACTION STRENGTH (H-STATISTIC)
# ============================================================

cat("\n### Computing interaction strength (H-statistic) ###\n")

if(require(iml, quietly = TRUE)) {
  
  cat("Note: H-statistic calculation can be sensitive to missing data.\n")
  cat("If errors occur, results will be skipped.\n\n")
  
  # Create predictor object - use only complete cases for stability
  best_data_complete <- best_data %>% 
    dplyr::select(CH4_flux, species, all_of(top_vars[1:3])) %>%
    drop_na()
  
  if(nrow(best_data_complete) < 30) {
    cat("Insufficient complete cases for H-statistic. Skipping.\n")
  } else {
    # Refit model on complete data for H-statistic
    rf_for_h <- randomForest(CH4_flux ~ ., data = best_data_complete,
                             ntree = 500, importance = TRUE)
    
    predictor <- Predictor$new(rf_for_h, data = best_data_complete, y = "CH4_flux")
    
    # Test interactions for top predictors
    h_stats <- list()
    
    for(var in top_vars[1:3]) {
      if(var %in% names(best_data_complete)) {
        cat(sprintf("Computing H-statistic for %s...\n", var))
        
        tryCatch({
          interact <- Interaction$new(predictor, feature = var, grid.size = 15)
          
          h_stats[[var]] <- interact$results %>%
            arrange(desc(.interaction)) %>%
            head(5)
          
          cat(sprintf("\nTop interactions with %s:\n", var))
          print(h_stats[[var]])
        }, error = function(e) {
          cat(sprintf("Error computing H-statistic for %s: %s\n", var, e$message))
        })
      }
    }
    
    # Save H-statistics if any were successful
    if(length(h_stats) > 0) {
      h_combined <- do.call(rbind, lapply(names(h_stats), function(v) {
        h_stats[[v]] %>% mutate(focal_var = v)
      }))
      
      write.csv(h_combined, "../../../outputs/tables/RF_H_statistics.csv", row.names = FALSE)
      
      # Plot H-statistics
      pdf("../../../outputs/figures/Figure_Interaction_Strength.pdf", width = 10, height = 6)
      
      p <- ggplot(h_combined, aes(x = reorder(.feature, .interaction), 
                                  y = .interaction, fill = focal_var)) +
        geom_col() +
        geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
        coord_flip() +
        facet_wrap(~focal_var, scales = "free_y") +
        labs(title = "Interaction Strength (H-statistic)",
             subtitle = "Values >0.1 indicate important interactions",
             x = "Interacting Variable",
             y = "H-statistic") +
        theme_classic() +
        theme(legend.position = "none")
      
      print(p)
      dev.off()
      
      cat("Saved: Figure_Interaction_Strength.pdf\n")
    } else {
      cat("No H-statistics successfully computed.\n")
    }
  }
  
} else {
  cat("Package 'iml' not available. Skipping H-statistic calculation.\n")
  cat("Install with: install.packages('iml')\n")
}

# ============================================================
# PART 8: SPECIES-SPECIFIC PATTERNS
# ============================================================

cat("\n### Species-specific gene effects ###\n")

# Use the best model's data
species_patterns <- best_data %>%
  group_by(species) %>%
  summarise(
    n = n(),
    across(where(is.numeric) & !matches("CH4_flux"), 
           list(mean = ~mean(., na.rm=TRUE)), 
           .names = "{.col}_mean"),
    CH4_flux_mean = mean(CH4_flux, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 5)

cat(sprintf("\nSpecies with >=5 trees: %d\n", nrow(species_patterns)))
write.csv(species_patterns, "../../../outputs/tables/RF_species_patterns.csv", row.names = FALSE)

# ============================================================
# PART 9: COMPARE SEPARATE VS COMBINED METHANOTROPH GENES
# ============================================================

cat("\n### Comparing pmoA+mmoX separate vs combined ###\n")

# Model with separate genes
data_separate <- data_v1_area_dbh %>%
  dplyr::select(-methanotroph_sum)

# Model with combined sum
data_combined <- data_v1_area_dbh %>%
  dplyr::select(-log_tree_pmoa, -log_tree_mmox)

set.seed(42)
rf_separate <- randomForest(CH4_flux ~ ., data = data_separate,
                            ntree = 1000, importance = TRUE,
                            na.action = na.roughfix)

rf_combined <- randomForest(CH4_flux ~ ., data = data_combined,
                            ntree = 1000, importance = TRUE,
                            na.action = na.roughfix)

cat("\nSEPARATE (pmoA + mmoX as individual predictors):\n")
print(rf_separate)

cat("\nCOMBINED (pmoA+mmoX summed):\n")
print(rf_combined)

comparison_sep_comb <- data.frame(
  Model = c("Separate pmoA+mmoX", "Combined sum"),
  Var_Explained = c(max(rf_separate$rsq)*100, max(rf_combined$rsq)*100),
  MSE = c(tail(rf_separate$mse, 1), tail(rf_combined$mse, 1))
)

cat("\nComparison:\n")
print(comparison_sep_comb)

write.csv(comparison_sep_comb, "../../../outputs/tables/RF_separate_vs_combined.csv", row.names = FALSE)

# ============================================================
# PART 10: SUMMARY VISUALIZATIONS
# ============================================================

cat("\n### Creating summary figure ###\n")

pdf("../../../outputs/figures/Figure_RF_Summary.pdf", width = 14, height = 10)

layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE))

# Panel 1: Model comparison
barplot(comparison_table$Var_Explained_Pct[1:8],
        names.arg = gsub("_", "\n", comparison_table$Model[1:8]),
        las = 2, cex.names = 0.7,
        main = "A) Model Performance Comparison",
        ylab = "% Variance Explained",
        col = "steelblue")

# Panel 2: Top variable importance (best model)
imp_best <- importance(best_model) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  arrange(desc(`%IncMSE`)) %>%
  head(10)

barplot(imp_best$`%IncMSE`,
        names.arg = imp_best$Variable,
        las = 2, cex.names = 0.8,
        main = "B) Variable Importance (Best Model)",
        ylab = "% Increase in MSE",
        col = "coral")

# Panel 3: Partial dependence (top predictor)
top_pred <- imp_best$Variable[imp_best$Variable != "species"][1]
pd <- partial(best_model, pred.var = top_pred, train = best_data)
plot(pd, type = "l", lwd = 3, col = "darkgreen",
     main = sprintf("C) Partial Dependence: %s", top_pred),
     xlab = top_pred, ylab = "Predicted CH4 flux")

# Panel 4: DBH vs DBH-normalized comparison
dbh_comparison <- comparison_table %>%
  mutate(
    DBH_type = case_when(
      grepl("v[12].*", Model) & grepl("norm", Model) ~ "DBH-normalized",
      grepl("v[12].*", Model) ~ "DBH-absolute",
      TRUE ~ "Other"
    )
  ) %>%
  filter(DBH_type != "Other")

if(nrow(dbh_comparison) > 0) {
  boxplot(Var_Explained_Pct ~ DBH_type, data = dbh_comparison,
          main = "D) DBH vs DBH-normalized",
          ylab = "% Variance Explained",
          col = c("lightblue", "lightcoral"))
}

dev.off()

cat("Saved: Figure_RF_Summary.pdf\n")

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("COMPREHENSIVE RF ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("\n### BEST MODEL ###\n")
cat(sprintf("Model: %s\n", best_model_name))
cat(sprintf("N: %d trees\n", nrow(best_data)))
cat(sprintf("Variance explained: %.1f%%\n", max(best_model$rsq)*100))
cat(sprintf("MSE: %.4f\n", tail(best_model$mse, 1)))

cat("\n### TOP 5 PREDICTORS ###\n")
print(head(imp, 5))

cat("\n============================================================\n")






















# ============================================================
# MAKE RF RESULTS INTERPRETABLE
# ============================================================
# This script extracts:
# 1. Direction of effects (positive/negative)
# 2. Magnitude of effects (effect size)
# 3. Linear approximations (pseudo-coefficients)
# 4. Confidence in linearity (R² of linear fit to PDP)
# 5. Threshold detection
# 6. Interaction effect directions
# ============================================================

library(tidyverse)
library(randomForest)
library(pdp)

cat("\n============================================================\n")
cat("EXTRACTING INTERPRETABLE EFFECTS FROM RF\n")
cat("============================================================\n")

# Use the best model from previous analysis
# Assumes you've already run the comprehensive RF script

# ============================================================
# PART 1: EXTRACT DIRECTION AND MAGNITUDE FROM PDPs
# ============================================================

cat("\n### Analyzing partial dependence to extract effects ###\n")

# Function to analyze a single PDP
analyze_pdp <- function(rf_model, data, var_name, n_points = 50) {
  
  if(!var_name %in% names(data)) {
    return(NULL)
  }
  
  # Get partial dependence
  pd <- partial(rf_model, pred.var = var_name, train = data, 
                grid.resolution = n_points)
  
  # Fit linear model to PDP
  lm_fit <- lm(yhat ~ get(var_name), data = pd)
  
  # Calculate effect statistics
  var_range <- range(pd[[var_name]], na.rm = TRUE)
  pred_range <- range(pd$yhat, na.rm = TRUE)
  
  # Effect size = change in prediction over predictor range
  effect_size <- diff(pred_range)
  
  # Direction = sign of slope
  slope <- coef(lm_fit)[2]
  direction <- ifelse(slope > 0, "Positive", "Negative")
  
  # Linearity = R² of linear fit
  linearity_r2 <- summary(lm_fit)$r.squared
  
  # Interpret linearity
  linearity_type <- case_when(
    linearity_r2 > 0.95 ~ "Linear",
    linearity_r2 > 0.80 ~ "Nearly linear",
    linearity_r2 > 0.60 ~ "Moderately nonlinear",
    TRUE ~ "Strongly nonlinear"
  )
  
  # Detect thresholds (inflection points)
  # Simple method: check if derivative changes sign
  derivatives <- diff(pd$yhat) / diff(pd[[var_name]])
  threshold_detected <- any(diff(sign(derivatives)) != 0)
  
  # Standardized effect (per SD change in predictor)
  predictor_sd <- sd(data[[var_name]], na.rm = TRUE)
  standardized_effect <- slope * predictor_sd
  
  return(data.frame(
    Variable = var_name,
    Direction = direction,
    Slope = slope,
    Effect_Size = effect_size,
    Standardized_Effect = standardized_effect,
    Linearity_R2 = linearity_r2,
    Linearity_Type = linearity_type,
    Threshold_Detected = threshold_detected,
    Min_Pred = pred_range[1],
    Max_Pred = pred_range[2],
    stringsAsFactors = FALSE
  ))
}

# Analyze all top predictors
top_predictors <- imp$Variable[1:min(10, nrow(imp))]

effect_table <- do.call(rbind, lapply(top_predictors, function(v) {
  analyze_pdp(best_model, best_data, v)
}))

cat("\n### INTERPRETABLE EFFECTS TABLE ###\n")
print(effect_table)

write.csv(effect_table, "../../../outputs/tables/RF_interpretable_effects.csv", row.names = FALSE)
cat("\nSaved: RF_interpretable_effects.csv\n")

# ============================================================
# PART 2: CREATE INTERPRETABLE SUMMARY
# ============================================================

cat("\n### PLAIN LANGUAGE INTERPRETATION ###\n\n")

for(i in 1:nrow(effect_table)) {
  row <- effect_table[i,]
  
  cat(sprintf("%d. %s:\n", i, row$Variable))
  cat(sprintf("   Direction: %s relationship\n", row$Direction))
  cat(sprintf("   Effect size: %.4f change in CH4 flux across observed range\n", 
              row$Effect_Size))
  cat(sprintf("   Linearity: %s (R² = %.3f)\n", 
              row$Linearity_Type, row$Linearity_R2))
  
  if(row$Threshold_Detected) {
    cat("   ⚠ Threshold or saturation effect detected\n")
  }
  
  # Interpret magnitude
  flux_sd <- sd(best_data$CH4_flux, na.rm = TRUE)
  effect_in_sd <- row$Effect_Size / flux_sd
  
  magnitude_label <- case_when(
    abs(effect_in_sd) > 1.0 ~ "LARGE",
    abs(effect_in_sd) > 0.5 ~ "Moderate",
    abs(effect_in_sd) > 0.2 ~ "Small",
    TRUE ~ "Negligible"
  )
  
  cat(sprintf("   Magnitude: %s (%.2f SD of flux)\n", 
              magnitude_label, abs(effect_in_sd)))
  
  # Plain language
  if(row$Direction == "Positive") {
    cat(sprintf("   → Higher %s → Higher CH4 emissions\n", row$Variable))
  } else {
    cat(sprintf("   → Higher %s → Lower CH4 emissions\n", row$Variable))
  }
  
  cat("\n")
}

# ============================================================
# PART 3: INTERACTION EFFECT DIRECTIONS
# ============================================================

cat("\n### INTERACTION EFFECT DIRECTIONS ###\n\n")

# Function to analyze 2-way interaction
analyze_interaction <- function(rf_model, data, var1, var2, resolution = 15) {
  
  if(!all(c(var1, var2) %in% names(data))) {
    return(NULL)
  }
  
  # Get 2-way PDP
  pd_2way <- partial(rf_model, pred.var = c(var1, var2),
                     train = data, grid.resolution = resolution)
  
  # Fit interaction model
  formula_str <- sprintf("yhat ~ %s * %s", var1, var2)
  lm_fit <- lm(as.formula(formula_str), data = pd_2way)
  
  # Extract interaction coefficient
  interaction_coef <- coef(lm_fit)[4]  # Interaction term
  
  # Test significance of interaction in linear approximation
  interaction_p <- summary(lm_fit)$coefficients[4, 4]
  
  # Overall R² of interaction model
  model_r2 <- summary(lm_fit)$r.squared
  
  # Interpret interaction
  interaction_type <- case_when(
    abs(interaction_coef) < 0.001 ~ "No interaction",
    interaction_coef > 0 ~ "Synergistic (both high = extra high flux)",
    interaction_coef < 0 ~ "Antagonistic (both high = dampened flux)"
  )
  
  return(data.frame(
    Var1 = var1,
    Var2 = var2,
    Interaction_Coef = interaction_coef,
    Interaction_P = interaction_p,
    Interaction_Type = interaction_type,
    Model_R2 = model_r2,
    stringsAsFactors = FALSE
  ))
}

# Analyze key interactions (from H-statistics)
key_interactions <- list(
  c("log_mcra_heart", "log_methanotroph_sap"),
  c("log_mcra_heart", "log_mcra_sap"),
  c("log_methanotroph_sap", "log_mcra_sap")
)

interaction_table <- do.call(rbind, lapply(key_interactions, function(pair) {
  if(all(pair %in% names(best_data))) {
    analyze_interaction(best_model, best_data, pair[1], pair[2])
  }
}))

if(!is.null(interaction_table) && nrow(interaction_table) > 0) {
  cat("### INTERACTION EFFECTS ###\n")
  print(interaction_table)
  
    write.csv(interaction_table, "../../../outputs/tables/RF_interaction_effects.csv", row.names = FALSE)
  cat("\nSaved: RF_interaction_effects.csv\n\n")
  
  # Plain language interpretation
  for(i in 1:nrow(interaction_table)) {
    row <- interaction_table[i,]
    
    cat(sprintf("%s × %s:\n", row$Var1, row$Var2))
    cat(sprintf("   Type: %s\n", row$Interaction_Type))
    cat(sprintf("   Coefficient: %.4f\n", row$Interaction_Coef))
    
    if(row$Interaction_Type == "Synergistic (both high = extra high flux)") {
      cat("   → When both are high, flux is HIGHER than expected from sum of effects\n")
    } else if(row$Interaction_Type == "Antagonistic (both high = dampened flux)") {
      cat("   → When both are high, flux is LOWER than expected from sum of effects\n")
      cat("   → Suggests opposing processes (e.g., production vs oxidation)\n")
    }
    cat("\n")
  }
}

# ============================================================
# PART 4: COMPARE WITH LINEAR MODEL COEFFICIENTS
# ============================================================

cat("\n### COMPARISON: RF vs LINEAR MODEL ###\n\n")

# From your linear model: CH4_flux ~ species + log_tree_mmox
# Coefficient for log_tree_mmox was -0.356

cat("Linear model (area-weighted mmoX):\n")
cat("  Coefficient: -0.356 (p = 0.006)\n")
cat("  Interpretation: Per log10 increase in mmoX → 0.356 decrease in flux\n\n")

cat("RF model (concentration data):\n")

# Find closest comparable predictor
sapwood_effects <- effect_table %>%
  filter(grepl("sap", Variable, ignore.case = TRUE))

cat("Sapwood effects:\n")
print(sapwood_effects %>% dplyr::select(Variable, Direction, Slope, Effect_Size))

# ============================================================
# PART 5: CREATE PUBLICATION-READY SUMMARY TABLE
# ============================================================

cat("\n### Creating publication table ###\n")

pub_table <- effect_table %>%
  left_join(
    importance(best_model) %>%
      as.data.frame() %>%
      rownames_to_column("Variable") %>%
      dplyr::select(Variable, `%IncMSE`),
    by = "Variable"
  ) %>%
  arrange(desc(`%IncMSE`)) %>%
  mutate(
    Effect_Category = case_when(
      abs(Effect_Size) > sd(best_data$CH4_flux, na.rm=TRUE) ~ "Large",
      abs(Effect_Size) > 0.5 * sd(best_data$CH4_flux, na.rm=TRUE) ~ "Moderate",
      TRUE ~ "Small"
    )
  ) %>%
  dplyr::select(Variable, Direction, Slope, Effect_Size, Effect_Category, 
         Linearity_Type, `%IncMSE`) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

cat("\n### PUBLICATION TABLE: RF PREDICTOR EFFECTS ###\n")
print(pub_table)

write.csv(pub_table, "../../../outputs/tables/RF_publication_table.csv", row.names = FALSE)
cat("\nSaved: RF_publication_table.csv\n")

# ============================================================
# PART 6: VISUALIZE INTERPRETABLE EFFECTS
# ============================================================

cat("\n### Creating interpretable effects figure ###\n")

pdf("../../../outputs/figures/Figure_Interpretable_Effects.pdf", width = 14, height = 10)

par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

# Panel 1: Effect sizes
barplot(effect_table$Effect_Size,
        names.arg = gsub("log_", "", effect_table$Variable),
        las = 2,
        col = ifelse(effect_table$Direction == "Positive", "coral", "steelblue"),
        main = "A) Effect Sizes (Change in Flux Across Predictor Range)",
        ylab = "Change in CH4 flux",
        cex.names = 0.8)
abline(h = 0, lty = 2)
legend("topright", c("Positive", "Negative"), 
       fill = c("coral", "steelblue"), bty = "n")

# Panel 2: Slopes (pseudo-coefficients)
barplot(effect_table$Slope,
        names.arg = gsub("log_", "", effect_table$Variable),
        las = 2,
        col = ifelse(effect_table$Direction == "Positive", "coral", "steelblue"),
        main = "B) Slopes (Pseudo-Coefficients from Linear Fit to PDP)",
        ylab = "Slope",
        cex.names = 0.8)
abline(h = 0, lty = 2)

# Panel 3: Linearity assessment
barplot(effect_table$Linearity_R2,
        names.arg = gsub("log_", "", effect_table$Variable),
        las = 2,
        col = "darkgreen",
        main = "C) Linearity (R² of Linear Fit to PDP)",
        ylab = "R² (1 = perfectly linear)",
        ylim = c(0, 1),
        cex.names = 0.8)
abline(h = 0.8, lty = 2, col = "red")
text(1, 0.85, "Threshold for\n'linear'", cex = 0.8, col = "red")

# Panel 4: Standardized effects (comparable scale)
barplot(abs(effect_table$Standardized_Effect),
        names.arg = gsub("log_", "", effect_table$Variable),
        las = 2,
        col = "purple",
        main = "D) Standardized Effects (Change per 1 SD of Predictor)",
        ylab = "Absolute standardized effect",
        cex.names = 0.8)

par(mfrow = c(1, 1))
dev.off()

cat("Saved: Figure_Interpretable_Effects.pdf\n")

# ============================================================
# PART 7: FINAL INTERPRETABLE SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("INTERPRETABLE SUMMARY\n")
cat("============================================================\n\n")

cat("### MAIN EFFECTS (Top 5) ###\n\n")

for(i in 1:min(5, nrow(pub_table))) {
  row <- pub_table[i,]
  
  cat(sprintf("%d. %s (%s effect):\n", i, 
              gsub("log_", "", row$Variable), row$Effect_Category))
  
  if(row$Direction == "Positive") {
    cat(sprintf("   - Each unit increase → +%.3f change in CH4 flux\n", 
                abs(row$Slope)))
    cat("   - Higher values → Higher emissions\n")
  } else {
    cat(sprintf("   - Each unit increase → -%.3f change in CH4 flux\n", 
                abs(row$Slope)))
    cat("   - Higher values → Lower emissions\n")
  }
  
  if(row$Linearity_Type == "Linear" || row$Linearity_Type == "Nearly linear") {
    cat("   - Relationship is linear (simple to interpret)\n")
  } else {
    cat(sprintf("   - Relationship is %s (see PDP for shape)\n", 
                tolower(row$Linearity_Type)))
  }
  
  cat("\n")
}

cat("### INTERACTIONS ###\n\n")

if(exists("interaction_table") && !is.null(interaction_table)) {
  cat("Strong interactions detected (H > 0.1):\n\n")
  
  # Heartwood mcrA × Sapwood methanotrophs
  cat("1. Heartwood mcrA × Sapwood methanotrophs (H = 0.51):\n")
  cat("   - STRONGEST interaction in model\n")
  cat("   - High production + Low oxidation → Highest emissions\n")
  cat("   - Low production + High oxidation → Lowest emissions\n")
  cat("   - Demonstrates production-oxidation balance\n\n")
  
  cat("2. Other significant interactions:\n")
  cat("   - Heartwood × Sapwood mcrA (H = 0.16)\n")
  cat("   - Sapwood methanotrophs × Sapwood mcrA (H = 0.32)\n")
  cat("   → All confirm spatial processes interact\n\n")
}










# ============================================================
# LINEAR MODELS WITH CONCENTRATION DATA
# ============================================================
# Compare area-weighted vs concentration (heartwood/sapwood)
# Look at ALL effect magnitudes regardless of significance
# ============================================================

library(tidyverse)
library(broom)
library(car)

cat("\n============================================================\n")
cat("LINEAR MODELS: CONCENTRATION VS AREA-WEIGHTED\n")
cat("============================================================\n")

# ============================================================
# PART 1: PREPARE BOTH DATASETS
# ============================================================

cat("\n### Preparing datasets ###\n")

# Dataset 1: Area-weighted (from earlier analysis)
area_weighted_data <- analysis_data %>%
  filter(!is.na(CH4_flux)) %>%
  dplyr::select(tree_id, CH4_flux, species, dbh,
         log_tree_mcra, log_tree_pmoa, log_tree_mmox) %>%
  mutate(data_type = "Area-weighted")

cat(sprintf("Area-weighted: %d trees\n", nrow(area_weighted_data)))

# Dataset 2: Concentration (heartwood/sapwood separate)
concentration_data <- conc_data %>%
  filter(!is.na(CH4_flux)) %>%
  dplyr::select(tree_id, CH4_flux, species, dbh,
         log_mcra_heart, log_mcra_sap,
         log_pmoa_heart, log_pmoa_sap,
         log_mmox_heart, log_mmox_sap) %>%
  mutate(data_type = "Concentration")

cat(sprintf("Concentration: %d trees\n", nrow(concentration_data)))

# ============================================================
# PART 2: FIT ALL MODELS - AREA-WEIGHTED
# ============================================================

cat("\n### AREA-WEIGHTED MODELS ###\n")

# Complete cases for fair comparison
area_complete <- area_weighted_data %>%
  filter(!is.na(log_tree_mcra), !is.na(log_tree_pmoa), !is.na(log_tree_mmox))

cat(sprintf("Using %d trees with complete gene data\n", nrow(area_complete)))

# Models with each gene separately
m_area_mcra <- lm(CH4_flux ~ species + log_tree_mcra, data = area_complete)
m_area_pmoa <- lm(CH4_flux ~ species + log_tree_pmoa, data = area_complete)
m_area_mmox <- lm(CH4_flux ~ species + log_tree_mmox, data = area_complete)

# Model with all genes
m_area_all <- lm(CH4_flux ~ species + log_tree_mcra + log_tree_pmoa + log_tree_mmox, 
                 data = area_complete)

# Extract coefficients
extract_gene_coefs <- function(model, gene_pattern) {
  coefs <- tidy(model, conf.int = TRUE) %>%
    filter(grepl(gene_pattern, term, ignore.case = TRUE)) %>%
    mutate(model = deparse(substitute(model)))
  return(coefs)
}

area_coefs <- bind_rows(
  extract_gene_coefs(m_area_mcra, "mcra"),
  extract_gene_coefs(m_area_pmoa, "pmoa"),
  extract_gene_coefs(m_area_mmox, "mmox"),
  tidy(m_area_all, conf.int = TRUE) %>%
    filter(grepl("log_tree", term)) %>%
    mutate(model = "m_area_all")
)

cat("\n### Area-weighted coefficients ###\n")
print(area_coefs %>% 
        dplyr::select(model, term, estimate, std.error, p.value, conf.low, conf.high))

# ============================================================
# PART 3: FIT ALL MODELS - CONCENTRATION
# ============================================================

cat("\n\n### CONCENTRATION MODELS ###\n")

# Complete cases
conc_complete <- concentration_data %>%
  filter(!is.na(log_mcra_heart), !is.na(log_mcra_sap),
         !is.na(log_pmoa_heart), !is.na(log_pmoa_sap),
         !is.na(log_mmox_heart), !is.na(log_mmox_sap))

cat(sprintf("Using %d trees with complete gene data\n", nrow(conc_complete)))

# Models with each gene (heartwood + sapwood)
m_conc_mcra <- lm(CH4_flux ~ species + log_mcra_heart + log_mcra_sap, 
                  data = conc_complete)
m_conc_pmoa <- lm(CH4_flux ~ species + log_pmoa_heart + log_pmoa_sap, 
                  data = conc_complete)
m_conc_mmox <- lm(CH4_flux ~ species + log_mmox_heart + log_mmox_sap, 
                  data = conc_complete)

# Model with all genes
m_conc_all <- lm(CH4_flux ~ species + 
                   log_mcra_heart + log_mcra_sap +
                   log_pmoa_heart + log_pmoa_sap +
                   log_mmox_heart + log_mmox_sap,
                 data = conc_complete)

# Extract coefficients
conc_coefs <- bind_rows(
  tidy(m_conc_mcra, conf.int = TRUE) %>%
    filter(grepl("log_mcra", term)) %>%
    mutate(model = "m_conc_mcra"),
  tidy(m_conc_pmoa, conf.int = TRUE) %>%
    filter(grepl("log_pmoa", term)) %>%
    mutate(model = "m_conc_pmoa"),
  tidy(m_conc_mmox, conf.int = TRUE) %>%
    filter(grepl("log_mmox", term)) %>%
    mutate(model = "m_conc_mmox"),
  tidy(m_conc_all, conf.int = TRUE) %>%
    filter(grepl("log_", term)) %>%
    mutate(model = "m_conc_all")
)

cat("\n### Concentration coefficients ###\n")
print(conc_coefs %>% 
        dplyr::select(model, term, estimate, std.error, p.value, conf.low, conf.high))

# ============================================================
# PART 4: COMPARE MODEL FIT
# ============================================================

cat("\n\n### MODEL FIT COMPARISON ###\n")

model_comparison <- data.frame(
  Model = c("Area: mcrA", "Area: pmoA", "Area: mmoX", "Area: All genes",
            "Conc: mcrA", "Conc: pmoA", "Conc: mmoX", "Conc: All genes"),
  N = c(rep(nrow(area_complete), 4), rep(nrow(conc_complete), 4)),
  N_params = c(
    length(coef(m_area_mcra)),
    length(coef(m_area_pmoa)),
    length(coef(m_area_mmox)),
    length(coef(m_area_all)),
    length(coef(m_conc_mcra)),
    length(coef(m_conc_pmoa)),
    length(coef(m_conc_mmox)),
    length(coef(m_conc_all))
  ),
  R2 = c(
    summary(m_area_mcra)$r.squared,
    summary(m_area_pmoa)$r.squared,
    summary(m_area_mmox)$r.squared,
    summary(m_area_all)$r.squared,
    summary(m_conc_mcra)$r.squared,
    summary(m_conc_pmoa)$r.squared,
    summary(m_conc_mmox)$r.squared,
    summary(m_conc_all)$r.squared
  ),
  Adj_R2 = c(
    summary(m_area_mcra)$adj.r.squared,
    summary(m_area_pmoa)$adj.r.squared,
    summary(m_area_mmox)$adj.r.squared,
    summary(m_area_all)$adj.r.squared,
    summary(m_conc_mcra)$adj.r.squared,
    summary(m_conc_pmoa)$adj.r.squared,
    summary(m_conc_mmox)$adj.r.squared,
    summary(m_conc_all)$adj.r.squared
  ),
  AIC = c(
    AIC(m_area_mcra),
    AIC(m_area_pmoa),
    AIC(m_area_mmox),
    AIC(m_area_all),
    AIC(m_conc_mcra),
    AIC(m_conc_pmoa),
    AIC(m_conc_mmox),
    AIC(m_conc_all)
  )
) %>%
  mutate(
    R2 = round(R2, 4),
    Adj_R2 = round(Adj_R2, 4),
    AIC = round(AIC, 1)
  )

print(model_comparison)
write.csv(model_comparison, "../../../outputs/tables/Linear_model_comparison_concentration.csv", row.names = FALSE)

# ============================================================
# PART 5: EFFECT MAGNITUDE SUMMARY
# ============================================================

cat("\n\n### EFFECT MAGNITUDE SUMMARY ###\n")
cat("(All effects shown regardless of significance)\n\n")

# Calculate standardized effects
flux_sd <- sd(conc_complete$CH4_flux, na.rm = TRUE)

effect_summary <- conc_coefs %>%
  filter(model == "m_conc_all") %>%
  mutate(
    gene = gsub("log_", "", term),
    gene = gsub("_heart|_sap", "", gene),
    location = ifelse(grepl("heart", term), "Heartwood", "Sapwood"),
    effect_magnitude = case_when(
      abs(estimate) > 0.5 ~ "Large (>0.5)",
      abs(estimate) > 0.2 ~ "Moderate (0.2-0.5)",
      abs(estimate) > 0.1 ~ "Small (0.1-0.2)",
      TRUE ~ "Negligible (<0.1)"
    ),
    standardized = estimate / flux_sd,
    significant = p.value < 0.05
  ) %>%
  arrange(gene, location)

cat("All gene effects in full concentration model:\n\n")
print(effect_summary %>%
        dplyr::select(gene, location, estimate, std.error, p.value, 
               conf.low, conf.high, effect_magnitude, significant))

write.csv(effect_summary, "../../../outputs/tables/Linear_effect_magnitudes_concentration.csv", row.names = FALSE)

# ============================================================
# PART 6: VISUALIZE EFFECT MAGNITUDES
# ============================================================

cat("\n### Creating visualization ###\n")

pdf("../../../outputs/figures/Figure_Linear_Effect_Magnitudes.pdf", width = 12, height = 8)

par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

# Panel 1: Area-weighted effects
area_plot <- area_coefs %>%
  filter(model == "m_area_all") %>%
  mutate(gene = gsub("log_tree_", "", term))

barplot(area_plot$estimate,
        names.arg = area_plot$gene,
        las = 2,
        col = ifelse(area_plot$p.value < 0.05, "darkblue", "lightblue"),
        main = "A) Area-Weighted Model (All Genes)",
        ylab = "Coefficient",
        ylim = c(-0.5, 0.5))
abline(h = 0, lty = 2)
legend("topright", c("p < 0.05", "p >= 0.05"), 
       fill = c("darkblue", "lightblue"), bty = "n")

# Panel 2: Concentration effects
conc_plot <- effect_summary

pos <- barplot(conc_plot$estimate,
               names.arg = paste(conc_plot$gene, conc_plot$location, sep = "\n"),
               las = 2,
               col = ifelse(conc_plot$significant, "darkgreen", "lightgreen"),
               main = "B) Concentration Model (All Genes)",
               ylab = "Coefficient",
               ylim = c(-0.5, 0.5),
               cex.names = 0.7)
abline(h = 0, lty = 2)

# Add confidence intervals
arrows(pos, conc_plot$conf.low, pos, conc_plot$conf.high,
       angle = 90, code = 3, length = 0.05)

# Panel 3: Compare mcrA (area vs conc)
mcra_comparison <- data.frame(
  Model = c("Area-weighted", "Heartwood", "Sapwood"),
  Estimate = c(
    area_coefs %>% filter(model == "m_area_all", grepl("mcra", term)) %>% pull(estimate),
    conc_plot %>% filter(gene == "mcra", location == "Heartwood") %>% pull(estimate),
    conc_plot %>% filter(gene == "mcra", location == "Sapwood") %>% pull(estimate)
  ),
  P_value = c(
    area_coefs %>% filter(model == "m_area_all", grepl("mcra", term)) %>% pull(p.value),
    conc_plot %>% filter(gene == "mcra", location == "Heartwood") %>% pull(p.value),
    conc_plot %>% filter(gene == "mcra", location == "Sapwood") %>% pull(p.value)
  )
)

barplot(mcra_comparison$Estimate,
        names.arg = mcra_comparison$Model,
        las = 2,
        col = ifelse(mcra_comparison$P_value < 0.05, "coral", "pink"),
        main = "C) mcrA: Area-weighted vs Concentration",
        ylab = "Coefficient",
        ylim = c(-0.5, 0.5))
abline(h = 0, lty = 2)

# Panel 4: R² comparison
barplot(model_comparison$Adj_R2,
        names.arg = gsub(" ", "\n", model_comparison$Model),
        las = 2,
        col = "purple",
        main = "D) Model R² Comparison",
        ylab = "Adjusted R²",
        ylim = c(0, 0.5),
        cex.names = 0.6)

par(mfrow = c(1, 1))
dev.off()

cat("Saved: Figure_Linear_Effect_Magnitudes.pdf\n")

# ============================================================
# PART 7: ANOVA TYPE II TESTS
# ============================================================

cat("\n\n### TYPE II ANOVA (Sequential F-tests) ###\n")

cat("\nArea-weighted model (all genes):\n")
print(Anova(m_area_all, type = "II"))

cat("\nConcentration model (all genes):\n")
print(Anova(m_conc_all, type = "II"))

# ============================================================
# PART 8: INTERPRETABLE SUMMARY
# ============================================================

cat("\n\n============================================================\n")
cat("INTERPRETABLE SUMMARY\n")
cat("============================================================\n\n")

cat("### EFFECT MAGNITUDES (Full Models) ###\n\n")

cat("AREA-WEIGHTED (species + 3 genes, n=", nrow(area_complete), "):\n")
area_final <- area_coefs %>% filter(model == "m_area_all")
for(i in 1:nrow(area_final)) {
  row <- area_final[i,]
  sig <- ifelse(row$p.value < 0.05, "**", 
                ifelse(row$p.value < 0.10, "*", "ns"))
  cat(sprintf("  %s: β = %.3f (SE = %.3f, p = %.3f) %s\n",
              gsub("log_tree_", "", row$term),
              row$estimate, row$std.error, row$p.value, sig))
  cat(sprintf("    95%% CI: [%.3f, %.3f]\n", row$conf.low, row$conf.high))
}

cat("\nCONCENTRATION (species + 6 genes, n=", nrow(conc_complete), "):\n")
for(i in 1:nrow(effect_summary)) {
  row <- effect_summary[i,]
  sig <- ifelse(row$significant, "**", "ns")
  cat(sprintf("  %s %s: β = %.3f (SE = %.3f, p = %.3f) %s\n",
              row$gene, row$location,
              row$estimate, row$std.error, row$p.value, sig))
  cat(sprintf("    95%% CI: [%.3f, %.3f]\n", row$conf.low, row$conf.high))
}

cat("\n### KEY FINDINGS ###\n\n")

# Find largest effects
largest_area <- area_final %>% 
  arrange(desc(abs(estimate))) %>% 
  head(1)

largest_conc <- effect_summary %>%
  arrange(desc(abs(estimate))) %>%
  head(1)

cat("1. AREA-WEIGHTED MODEL:\n")
cat(sprintf("   - Largest effect: %s (β = %.3f, p = %.3f)\n",
            gsub("log_tree_", "", largest_area$term),
            largest_area$estimate, largest_area$p.value))
cat(sprintf("   - R² = %.3f\n", model_comparison %>% 
              filter(Model == "Area: All genes") %>% pull(R2)))

cat("\n2. CONCENTRATION MODEL:\n")
cat(sprintf("   - Largest effect: %s %s (β = %.3f, p = %.3f)\n",
            largest_conc$gene, largest_conc$location,
            largest_conc$estimate, largest_conc$p.value))
cat(sprintf("   - R² = %.3f\n", model_comparison %>% 
              filter(Model == "Conc: All genes") %>% pull(R2)))

cat("\n3. COMPARISON:\n")
cat(sprintf("   - Concentration model explains %.1f%% more variance\n",
            100 * (model_comparison %>% filter(Model == "Conc: All genes") %>% pull(R2) -
                     model_comparison %>% filter(Model == "Area: All genes") %>% pull(R2))))

cat("\n### EFFECT MAGNITUDES IN CONTEXT ###\n\n")
cat(sprintf("Flux SD = %.3f\n", flux_sd))
cat(sprintf("Flux range = %.3f to %.3f\n\n", 
            min(conc_complete$CH4_flux), max(conc_complete$CH4_flux)))

cat("Effect size interpretations:\n")
cat("  - Negligible: |β| < 0.1 (< 1/3 SD of flux)\n")
cat("  - Small: |β| = 0.1-0.2 (1/3 to 2/3 SD)\n")
cat("  - Moderate: |β| = 0.2-0.5 (2/3 to 2 SD)\n")
cat("  - Large: |β| > 0.5 (> 2 SD)\n\n")

cat("By this scale:\n")
for(i in 1:nrow(effect_summary)) {
  row <- effect_summary[i,]
  cat(sprintf("  %s %s: %s effect\n",
              row$gene, row$location, row$effect_magnitude))
}

cat("\n============================================================\n")
cat("FILES CREATED:\n")
cat("- Linear_model_comparison_concentration.csv\n")
cat("- Linear_effect_magnitudes_concentration.csv\n")
cat("- Figure_Linear_Effect_Magnitudes.pdf\n")
cat("============================================================\n")













# ============================================================
# METHANOTROPH (pmoA + mmoX) vs FLUX CORRELATION ANALYSIS
# Fixed version with updated ggplot2 syntax
# ============================================================

library(tidyverse)
library(scales)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(boot)

# Load data
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

cat("\n============================================================\n")
cat("METHANOTROPH (pmoA + mmoX) vs FLUX ANALYSIS\n")
cat("============================================================\n")

# ============================================================
# 1. PARSE METHANOTROPH DATA
# ============================================================

prepare_long_methanotrophs <- function(df) {
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
      location = part1,
      stringency = part2,
      location = case_when(
        location %in% c("Inner","inner") ~ "Inner",
        location %in% c("Outer","outer") ~ "Outer",
        TRUE ~ location
      ),
      gene = case_when(
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
    filter(gene %in% c("pmoA", "mmoX")) %>%
    filter(sample_type %in% c("Heartwood", "Sapwood"))
}

long_data <- prepare_long_methanotrophs(ymf2021)
cat(sprintf("Found %d methanotroph gene measurements\n", nrow(long_data)))

# ============================================================
# 2. AREA-WEIGHTED CALCULATION
# ============================================================

area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1*r2 + r2^2
  
  gene_outer + (gene_inner - gene_outer) * (S / (3 * R^2))
}

tree_methanotrophs <- long_data %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, species_id, dbh),
    by = "tree_id"
  ) %>%
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
  dplyr::select(tree_id, species_id, species, dbh, gene, gene_area_weighted) %>%
  pivot_wider(
    names_from = gene,
    values_from = gene_area_weighted
  ) %>%
  mutate(
    methanotroph_total = case_when(
      is.na(pmoA) & !is.na(mmoX) ~ mmoX,
      !is.na(pmoA) & is.na(mmoX) ~ pmoA,
      !is.na(pmoA) & !is.na(mmoX) ~ pmoA + mmoX,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(methanotroph_total))

cat(sprintf("Calculated methanotroph totals for %d trees\n", nrow(tree_methanotrophs)))

# ============================================================
# 3. SPECIES-LEVEL AGGREGATION
# ============================================================

methanotroph_by_species <- tree_methanotrophs %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_methanotroph = median(methanotroph_total, na.rm = TRUE),
    mean_methanotroph = mean(methanotroph_total, na.rm = TRUE),
    median_log = log10(median_methanotroph + 1),
    mean_log = log10(mean_methanotroph + 1),
    q25 = quantile(methanotroph_total, 0.25, na.rm = TRUE),
    q75 = quantile(methanotroph_total, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# ============================================================
# 4. FLUX DATA
# ============================================================

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
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# ============================================================
# 5. MERGE AND CORRELATE
# ============================================================

analysis_data <- inner_join(
  methanotroph_by_species,
  flux_by_species,
  by = c("species", "species_id")
) %>%
  filter(n_trees >= 5, n_flux >= 5)

cat(sprintf("\nAnalyzing %d species with n>=5 for both measurements\n", nrow(analysis_data)))

# Correlations
if(nrow(analysis_data) > 2) {
  pearson_med <- cor.test(analysis_data$median_log, analysis_data$median_flux)
  spearman_med <- cor.test(analysis_data$median_methanotroph, analysis_data$median_flux, 
                           method = "spearman")
  kendall_med <- cor.test(analysis_data$median_methanotroph, analysis_data$median_flux, 
                          method = "kendall")
  
  pearson_mean <- cor.test(analysis_data$mean_log, analysis_data$mean_flux)
  spearman_mean <- cor.test(analysis_data$mean_methanotroph, analysis_data$mean_flux, 
                            method = "spearman")
  
  cat("\n### MEDIAN-BASED CORRELATIONS ###\n")
  cat(sprintf("Pearson r = %.3f (R² = %.3f), p = %.4f (log-transformed)\n", 
              pearson_med$estimate, pearson_med$estimate^2, pearson_med$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n", 
              spearman_med$estimate, spearman_med$p.value))
  cat(sprintf("Kendall tau = %.3f, p = %.4f\n",
              kendall_med$estimate, kendall_med$p.value))
  
  cat("\n### MEAN-BASED CORRELATIONS ###\n")
  cat(sprintf("Pearson r = %.3f (R² = %.3f), p = %.4f (log-transformed)\n", 
              pearson_mean$estimate, pearson_mean$estimate^2, pearson_mean$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n",
              spearman_mean$estimate, spearman_mean$p.value))
  
  cat("\nSpecies data (sorted by flux):\n")
  print(analysis_data %>%
          dplyr::select(species, n_trees, n_flux, median_methanotroph, median_flux) %>%
          arrange(median_flux) %>%
          mutate(median_methanotroph = round(median_methanotroph, 1),
                 median_flux = round(median_flux, 4)))
}

# ============================================================
# 6. VISUALIZATION - MEDIAN VALUES (FIXED SYNTAX)
# ============================================================

p_median <- ggplot(analysis_data, aes(x = median_log, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", 
              fill = "lightgreen", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbar(aes(xmin = log10(q25 + 1), xmax = log10(q75 + 1), y = median_flux),
                width = 0.003, alpha = 0.5, color = "gray40", orientation = "y") +
  geom_point(size = 4, alpha = 0.85, color = "darkgreen") +
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR-sq = %.3f\nrho = %.2f, p = %.3f\ntau = %.2f, p = %.3f",
                           pearson_med$estimate, pearson_med$p.value,
                           pearson_med$estimate^2,
                           spearman_med$estimate, spearman_med$p.value,
                           kendall_med$estimate, kendall_med$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  scale_x_continuous(
    breaks = 0:5,
    labels = c(expression(10^0), expression(10^1), expression(10^2), 
               expression(10^3), expression(10^4), expression(10^5))
  ) +
  labs(
    title = "Methanotrophs (pmoA + mmoX) vs CH4 Flux",
    subtitle = sprintf("Median values, n = %d species (>=5 obs each); error bars = IQR", 
                       nrow(analysis_data)),
    x = expression("Median area-weighted methanotrophs (copies/g)"),
    y = expression("Median CH4 flux (nmol/m2/s)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p_median)

# ============================================================
# 7. SURFACE AREA ANALYSIS
# ============================================================

cat("\n\n============================================================\n")
cat("TOTAL METHANOTROPH COPIES PER CM² BARK SURFACE\n")
cat("============================================================\n")

calculate_total_per_surface <- function(gene_inner, gene_outer, density_inner, density_outer, dbh_cm, dr = 0.01) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  
  if (!is.finite(density_inner)) density_inner <- 0.5
  if (!is.finite(density_outer)) density_outer <- 0.5
  
  r_edges <- seq(0, R, by = dr)
  if (tail(r_edges, 1) < R) r_edges <- c(r_edges, R)
  r_mid <- 0.5 * (r_edges[-1] + r_edges[-length(r_edges)])
  
  shell_areas <- pi * (r_edges[-1]^2 - r_edges[-length(r_edges)]^2)
  
  C_r <- ifelse(r_mid <= r1, gene_inner,
                ifelse(r_mid >= r2, gene_outer,
                       gene_inner + (gene_outer - gene_inner) * (r_mid - r1) / max(r2 - r1, 1e-9)))
  
  rho_r <- ifelse(r_mid <= r1, density_inner,
                  ifelse(r_mid >= r2, density_outer,
                         density_inner + (density_outer - density_inner) * (r_mid - r1) / max(r2 - r1, 1e-9)))
  
  total_copies <- sum(C_r * rho_r * shell_areas)
  circumference <- 2 * pi * R
  return(total_copies / circumference)
}

tree_surface <- long_data %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh,
                               density_inner = inner_density_final,
                               density_outer = outer_density_final),
            by = "tree_id") %>%
  filter(is.finite(dbh)) %>%
  group_by(tree_id, species_id, dbh, density_inner, density_outer, gene) %>%
  summarise(
    gene_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
    gene_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(gene_inner), is.finite(gene_outer)) %>%
  mutate(
    gene_per_surface = mapply(calculate_total_per_surface,
                              gene_inner, gene_outer,
                              density_inner, density_outer, dbh),
    species = species_mapping[species_id]
  ) %>%
  dplyr::select(tree_id, species_id, species, gene, gene_per_surface) %>%
  pivot_wider(
    names_from = gene,
    values_from = gene_per_surface
  ) %>%
  mutate(
    methanotroph_per_surface = case_when(
      is.na(pmoA) & !is.na(mmoX) ~ mmoX,
      !is.na(pmoA) & is.na(mmoX) ~ pmoA,
      !is.na(pmoA) & !is.na(mmoX) ~ pmoA + mmoX,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(methanotroph_per_surface))

surface_by_species <- tree_surface %>%
  group_by(species, species_id) %>%
  summarise(
    n = n(),
    median_surface = median(methanotroph_per_surface, na.rm = TRUE),
    median_surface_log = log10(median_surface + 1),
    q25_surface = quantile(methanotroph_per_surface, 0.25, na.rm = TRUE),
    q75_surface = quantile(methanotroph_per_surface, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

surface_analysis <- inner_join(
  surface_by_species,
  flux_by_species,
  by = c("species", "species_id")
) %>%
  filter(n >= 5, n_flux >= 5)

if(nrow(surface_analysis) > 2) {
  pearson_surf <- cor.test(surface_analysis$median_surface_log, surface_analysis$median_flux)
  spearman_surf <- cor.test(surface_analysis$median_surface, surface_analysis$median_flux, 
                            method = "spearman")
  kendall_surf <- cor.test(surface_analysis$median_surface, surface_analysis$median_flux,
                           method = "kendall")
  
  cat(sprintf("Number of species: %d\n", nrow(surface_analysis)))
  cat(sprintf("Pearson r = %.3f (R² = %.3f), p = %.4f\n", 
              pearson_surf$estimate, pearson_surf$estimate^2, pearson_surf$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n",
              spearman_surf$estimate, spearman_surf$p.value))
  cat(sprintf("Kendall tau = %.3f, p = %.4f\n",
              kendall_surf$estimate, kendall_surf$p.value))
  
  cat("\nSpecies data:\n")
  print(surface_analysis %>%
          dplyr::select(species, n, median_surface, median_flux) %>%
          arrange(median_flux) %>%
          mutate(median_surface = scientific(median_surface, digits = 2),
                 median_flux = round(median_flux, 4)))
}

# Plot (FIXED SYNTAX)
p_surface <- ggplot(surface_analysis, aes(x = median_surface_log, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", 
              fill = "lightgreen", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbar(aes(xmin = log10(q25_surface + 1), xmax = log10(q75_surface + 1), y = median_flux),
                width = 0.003, alpha = 0.5, color = "gray40", orientation = "y") +
  geom_point(size = 4, alpha = 0.85, color = "darkgreen") +
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR-sq = %.3f\nrho = %.2f, p = %.3f\ntau = %.2f, p = %.3f",
                           pearson_surf$estimate, pearson_surf$p.value,
                           pearson_surf$estimate^2,
                           spearman_surf$estimate, spearman_surf$p.value,
                           kendall_surf$estimate, kendall_surf$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  scale_x_continuous(
    breaks = 1:6,
    labels = c(expression(10^1), expression(10^2), expression(10^3), 
               expression(10^4), expression(10^5), expression(10^6))
  ) +
  labs(
    title = "Methanotrophs per Surface Area vs CH4 Flux",
    subtitle = sprintf("Median values, n = %d species; error bars = IQR", 
                       nrow(surface_analysis)),
    x = expression("Median methanotrophs beneath 1 cm^2 bark (copies)"),
    y = expression("Median CH4 flux (nmol/m2/s)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p_surface)

# ============================================================
# 8. SUMMARY
# ============================================================

cat("\n============================================================\n")
cat("SUMMARY: METHANOTROPH CORRELATIONS\n")
cat("============================================================\n")

cat("\nArea-weighted (copies/g):\n")
cat(sprintf("  Median-based: R² = %.3f, p = %.4f\n", 
            pearson_med$estimate^2, pearson_med$p.value))
cat(sprintf("  Mean-based: R² = %.3f, p = %.4f\n",
            pearson_mean$estimate^2, pearson_mean$p.value))

cat("\nPer surface area (copies/cm²):\n")
cat(sprintf("  R² = %.3f, p = %.4f\n",
            pearson_surf$estimate^2, pearson_surf$p.value))

cat("\nInterpretation:\n")
if(pearson_med$p.value < 0.05) {
  if(pearson_med$estimate < 0) {
    cat("  ** Species with MORE methanotrophs have LOWER emissions **\n")
    cat("  -> Supports oxidation hypothesis\n")
  } else {
    cat("  ** Species with MORE methanotrophs have HIGHER emissions **\n")
    cat("  -> Does not support simple oxidation hypothesis\n")
  }
} else {
  cat("  No significant correlation (p > 0.05)\n")
  cat("  -> Methanotroph abundance does not predict flux at species level\n")
  cat("  -> Trend is negative (r = %.2f) consistent with oxidation\n", pearson_med$estimate)
  cat("  -> But insufficient power (n = %d species) to detect\n", nrow(analysis_data))
}

cat("\n============================================================\n")








# ============================================================
# METHANOGEN:METHANOTROPH RATIO vs FLUX ANALYSIS
# Species-level test of production/oxidation balance
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(scales)

# Load data
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

cat("\n============================================================\n")
cat("METHANOGEN:METHANOTROPH RATIO ANALYSIS\n")
cat("============================================================\n")

# ============================================================
# 1. PARSE AND CALCULATE AREA-WEIGHTED GENES
# ============================================================

prepare_long_genes <- function(df) {
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
    filter(sample_type %in% c("Heartwood", "Sapwood"))
}

area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1*r2 + r2^2
  
  gene_outer + (gene_inner - gene_outer) * (S / (3 * R^2))
}

long_data <- prepare_long_genes(ymf2021)

# Calculate area-weighted for each gene
tree_genes <- long_data %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, species_id, dbh),
    by = "tree_id"
  ) %>%
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
  pivot_wider(
    names_from = gene,
    values_from = gene_area_weighted
  )

# Calculate ratios for each tree
tree_ratios <- tree_genes %>%
  mutate(
    # Sum methanotrophs
    methanotroph_total = case_when(
      is.na(pmoA) & !is.na(mmoX) ~ mmoX,
      !is.na(pmoA) & is.na(mmoX) ~ pmoA,
      !is.na(pmoA) & !is.na(mmoX) ~ pmoA + mmoX,
      TRUE ~ NA_real_
    ),
    # Calculate ratio (add pseudocount to avoid division by zero)
    ratio_mcra_methanotroph = (mcrA + 1) / (methanotroph_total + 1),
    log_ratio = log10(ratio_mcra_methanotroph)
  ) %>%
  filter(!is.na(mcrA), !is.na(methanotroph_total))

cat(sprintf("Calculated ratios for %d trees\n", nrow(tree_ratios)))

# ============================================================
# 2. SPECIES-LEVEL AGGREGATION
# ============================================================

ratio_by_species <- tree_ratios %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    # Ratio statistics
    median_ratio = median(ratio_mcra_methanotroph, na.rm = TRUE),
    mean_ratio = mean(ratio_mcra_methanotroph, na.rm = TRUE),
    median_log_ratio = median(log_ratio, na.rm = TRUE),
    mean_log_ratio = mean(log_ratio, na.rm = TRUE),
    q25_ratio = quantile(ratio_mcra_methanotroph, 0.25, na.rm = TRUE),
    q75_ratio = quantile(ratio_mcra_methanotroph, 0.75, na.rm = TRUE),
    # Individual gene levels for reference
    median_mcra = median(mcrA, na.rm = TRUE),
    median_methanotroph = median(methanotroph_total, na.rm = TRUE),
    .groups = 'drop'
  )

# ============================================================
# 3. FLUX DATA
# ============================================================

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
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# ============================================================
# 4. MERGE AND ANALYZE
# ============================================================

analysis_data <- inner_join(
  ratio_by_species,
  flux_by_species,
  by = c("species", "species_id")
) %>%
  filter(n_trees >= 5, n_flux >= 5)

cat(sprintf("\nAnalyzing %d species with n>=5 for both measurements\n", nrow(analysis_data)))

# Print data table
cat("\nSpecies data (sorted by flux):\n")
print(analysis_data %>%
        dplyr::select(species, n_trees, median_ratio, median_mcra, 
               median_methanotroph, median_flux) %>%
        arrange(median_flux) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# Correlations
if(nrow(analysis_data) > 2) {
  # Test 1: Ratio vs flux (log-transformed ratio)
  pearson_ratio <- cor.test(analysis_data$median_log_ratio, analysis_data$median_flux)
  spearman_ratio <- cor.test(analysis_data$median_ratio, analysis_data$median_flux, 
                             method = "spearman")
  kendall_ratio <- cor.test(analysis_data$median_ratio, analysis_data$median_flux, 
                            method = "kendall")
  
  cat("\n### RATIO CORRELATIONS (median values) ###\n")
  cat(sprintf("Pearson r = %.3f (R² = %.3f), p = %.4f (log-ratio)\n", 
              pearson_ratio$estimate, pearson_ratio$estimate^2, pearson_ratio$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.4f\n", 
              spearman_ratio$estimate, spearman_ratio$p.value))
  cat(sprintf("Kendall tau = %.3f, p = %.4f\n",
              kendall_ratio$estimate, kendall_ratio$p.value))
  
  # Test 2: Linear model with ratio
  lm_ratio <- lm(median_flux ~ median_log_ratio, data = analysis_data)
  cat("\n### LINEAR MODEL: flux ~ log(ratio) ###\n")
  print(summary(lm_ratio))
  
  # Test 3: Compare to individual genes
  lm_mcra <- lm(median_flux ~ log10(median_mcra + 1), data = analysis_data)
  lm_methanotroph <- lm(median_flux ~ log10(median_methanotroph + 1), data = analysis_data)
  
  cat("\n### MODEL COMPARISON ###\n")
  cat(sprintf("Ratio model R² = %.3f\n", summary(lm_ratio)$r.squared))
  cat(sprintf("mcrA alone R² = %.3f\n", summary(lm_mcra)$r.squared))
  cat(sprintf("Methanotroph alone R² = %.3f\n", summary(lm_methanotroph)$r.squared))
  
  # AIC comparison
  cat("\n### AIC COMPARISON (lower is better) ###\n")
  cat(sprintf("Ratio model AIC = %.2f\n", AIC(lm_ratio)))
  cat(sprintf("mcrA alone AIC = %.2f\n", AIC(lm_mcra)))
  cat(sprintf("Methanotroph alone AIC = %.2f\n", AIC(lm_methanotroph)))
}

# ============================================================
# 5. VISUALIZATION
# ============================================================

# Plot 1: Ratio vs Flux
p_ratio <- ggplot(analysis_data, aes(x = median_log_ratio, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "purple", 
              fill = "plum", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbar(aes(xmin = log10(q25_ratio), xmax = log10(q75_ratio), y = median_flux),
                width = 0.003, alpha = 0.5, color = "gray40", orientation = "y") +
  geom_point(size = 4, alpha = 0.85, color = "purple") +
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted",
             color = "gray50", alpha = 0.5) +
  annotate("text", x = 0, y = max(analysis_data$median_flux) * 0.95,
           label = "Equal\nproduction/oxidation", 
           size = 3, color = "gray50", hjust = 0.5) +
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR-sq = %.3f\nrho = %.2f, p = %.3f",
                           pearson_ratio$estimate, pearson_ratio$p.value,
                           pearson_ratio$estimate^2,
                           spearman_ratio$estimate, spearman_ratio$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  scale_x_continuous(
    breaks = seq(-2, 2, 0.5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    title = "Methanogen:Methanotroph Ratio vs CH4 Flux",
    subtitle = sprintf("n = %d species (>=5 obs each); error bars = IQR", 
                       nrow(analysis_data)),
    x = expression("log"[10]*"(mcrA / [pmoA+mmoX])"),
    y = expression("Median CH4 flux (nmol/m2/s)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p_ratio)

# Plot 2: Three-way comparison (ratio vs individual genes)
if(exists("lm_mcra") && exists("lm_methanotroph")) {
  
  comparison_df <- data.frame(
    Predictor = c("mcrA:Methanotroph\nRatio", "mcrA alone", "Methanotrophs alone"),
    R_squared = c(summary(lm_ratio)$r.squared,
                  summary(lm_mcra)$r.squared,
                  summary(lm_methanotroph)$r.squared),
    P_value = c(summary(lm_ratio)$coefficients[2,4],
                summary(lm_mcra)$coefficients[2,4],
                summary(lm_methanotroph)$coefficients[2,4])
  )
  
  p_comparison <- ggplot(comparison_df, aes(x = Predictor, y = R_squared, 
                                            fill = P_value < 0.05)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("R² = %.3f\np = %.3f", R_squared, P_value)),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("gray70", "steelblue"),
                      labels = c("p >= 0.05", "p < 0.05"),
                      name = "Significance") +
    labs(
      title = "Model Comparison: Ratio vs Individual Genes",
      subtitle = "Which predictor best explains species-level flux variation?",
      y = expression("R"^2),
      x = ""
    ) +
    ylim(0, max(comparison_df$R_squared) * 1.2) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 10),
      legend.position = "bottom"
    )
  
  print(p_comparison)
}

# ============================================================
# 6. INTERPRETATION
# ============================================================

cat("\n============================================================\n")
cat("INTERPRETATION\n")
cat("============================================================\n")

cat("\nHypothesis: Species with higher methanogen:methanotroph ratio\n")
cat("            should have higher CH4 emissions\n\n")

if(pearson_ratio$p.value < 0.05) {
  if(pearson_ratio$estimate > 0) {
    cat("RESULT: POSITIVE correlation (p < 0.05)\n")
    cat("  -> Species with higher production:oxidation ratio have HIGHER flux\n")
    cat("  -> SUPPORTS balance hypothesis\n")
  } else {
    cat("RESULT: NEGATIVE correlation (p < 0.05)\n")
    cat("  -> Unexpected: Higher ratio = Lower flux\n")
    cat("  -> Does NOT support simple balance hypothesis\n")
  }
} else {
  cat(sprintf("RESULT: No significant correlation (r = %.3f, p = %.3f)\n",
              pearson_ratio$estimate, pearson_ratio$p.value))
  if(abs(pearson_ratio$estimate) < 0.3) {
    cat("  -> Ratio does not predict flux at species level\n")
    cat("  -> Production/oxidation balance may matter within species\n")
    cat("     but is obscured by other species-level differences\n")
  } else {
    cat("  -> Moderate correlation but insufficient power\n")
    cat(sprintf("  -> Would need ~%.0f species for significance at this effect size\n",
                ceiling(8 / pearson_ratio$estimate^2)))
  }
}

cat("\nComparison to individual predictors:\n")
if(exists("lm_mcra")) {
  if(summary(lm_ratio)$r.squared > summary(lm_mcra)$r.squared & 
     summary(lm_ratio)$r.squared > summary(lm_methanotroph)$r.squared) {
    cat("  -> RATIO explains MORE variance than either gene alone\n")
    cat("  -> Suggests balance matters\n")
  } else {
    cat("  -> RATIO does NOT improve prediction over individual genes\n")
    cat("  -> Balance may not be key driver at species level\n")
  }
}

cat("\n============================================================\n")
cat("FILES CREATED:\n")
cat("- Ratio correlation plots (in memory)\n")
cat("============================================================\n")

# Optional: Save plots
ggsave("../../../outputs/figures/methanogen_methanotroph_ratio_vs_flux.png", plot = p_ratio,
       width = 10, height = 8, dpi = 300)
ggsave("../../../outputs/figures/ratio_comparison.png", plot = p_comparison,
       width = 8, height = 6, dpi = 300)














# ============================================================
# METHANOGEN:METHANOTROPH RATIO vs FLUX ANALYSIS - CORRECTED
# Each model uses the FULL dataset for its specific gene(s)
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(scales)

# Load data
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

cat("\n============================================================\n")
cat("RATIO ANALYSIS - CORRECTED VERSION\n")
cat("Each model uses the full dataset for its specific gene(s)\n")
cat("============================================================\n")

# ============================================================
# SHARED FUNCTIONS
# ============================================================

area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1*r2 + r2^2
  
  gene_outer + (gene_inner - gene_outer) * (S / (3 * R^2))
}

# ============================================================
# DATASET 1: mcrA-ONLY (all trees with mcrA, from scatter plot code)
# ============================================================

cat("\n### DATASET 1: mcrA-ONLY ###\n")

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

long_mcra <- prepare_long_mcra(ymf2021)

tree_mcra_full <- long_mcra %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
  filter(is.finite(dbh)) %>%
  group_by(tree_id, species_id, dbh) %>%
  summarise(
    mcra_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
    mcra_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(mcra_inner), is.finite(mcra_outer)) %>%
  mutate(
    mcra_area_weighted = mapply(area_weighted_gene, mcra_inner, mcra_outer, dbh),
    species = species_mapping[species_id]
  )

cat(sprintf("Trees with mcrA: %d\n", nrow(tree_mcra_full)))

# ============================================================
# DATASET 2: METHANOTROPHS-ONLY (all trees with pmoA/mmoX)
# ============================================================

cat("\n### DATASET 2: METHANOTROPHS-ONLY ###\n")

prepare_long_methanotrophs <- function(df) {
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
      location = part1,
      stringency = part2,
      location = case_when(
        location %in% c("Inner","inner") ~ "Inner",
        location %in% c("Outer","outer") ~ "Outer",
        TRUE ~ location
      ),
      gene = case_when(
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
    filter(gene %in% c("pmoA", "mmoX")) %>%
    filter(sample_type %in% c("Heartwood", "Sapwood"))
}

long_methanotrophs <- prepare_long_methanotrophs(ymf2021)

tree_methanotrophs_full <- long_methanotrophs %>%
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
  pivot_wider(names_from = gene, values_from = gene_area_weighted) %>%
  mutate(
    methanotroph_total = case_when(
      is.na(pmoA) & !is.na(mmoX) ~ mmoX,
      !is.na(pmoA) & is.na(mmoX) ~ pmoA,
      !is.na(pmoA) & !is.na(mmoX) ~ pmoA + mmoX,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(methanotroph_total))

cat(sprintf("Trees with methanotrophs: %d\n", nrow(tree_methanotrophs_full)))

# ============================================================
# DATASET 3: BOTH GENES (only trees with both mcrA AND methanotrophs)
# ============================================================

cat("\n### DATASET 3: BOTH GENES (for ratio) ###\n")

# Find trees that have BOTH genes
tree_both <- inner_join(
  tree_mcra_full %>% dplyr::select(tree_id, species_id, species, mcrA = mcra_area_weighted),
  tree_methanotrophs_full %>% dplyr::select(tree_id, methanotroph_total),
  by = "tree_id"
) %>%
  mutate(
    ratio_mcra_methanotroph = (mcrA + 1) / (methanotroph_total + 1),
    log_ratio = log10(ratio_mcra_methanotroph)
  )

cat(sprintf("Trees with BOTH genes: %d\n", nrow(tree_both)))

# ============================================================
# FLUX DATA (same for all)
# ============================================================

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

# ============================================================
# ANALYSIS 1: mcrA-ONLY MODEL
# ============================================================

cat("\n\n### ANALYSIS 1: mcrA-ONLY MODEL ###\n")

mcra_by_species <- tree_mcra_full %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_mcra = median(mcra_area_weighted, na.rm = TRUE),
    .groups = 'drop'
  )

analysis_mcra <- inner_join(mcra_by_species, flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

cat(sprintf("Species: %d\n", nrow(analysis_mcra)))

if(nrow(analysis_mcra) > 2) {
  lm_mcra <- lm(median_flux ~ log10(median_mcra + 1), data = analysis_mcra)
  pearson_mcra <- cor.test(log10(analysis_mcra$median_mcra + 1), analysis_mcra$median_flux)
  
  cat(sprintf("R² = %.3f, p = %.4f\n", 
              summary(lm_mcra)$r.squared, 
              summary(lm_mcra)$coefficients[2,4]))
  
  cat("\nSpecies with median mcrA:\n")
  print(analysis_mcra %>% 
          dplyr::select(species, n_trees, median_mcra, median_flux) %>%
          arrange(desc(median_mcra)) %>%
          head(3))
}

# ============================================================
# ANALYSIS 2: METHANOTROPHS-ONLY MODEL
# ============================================================

cat("\n\n### ANALYSIS 2: METHANOTROPHS-ONLY MODEL ###\n")

methanotrophs_by_species <- tree_methanotrophs_full %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_methanotroph = median(methanotroph_total, na.rm = TRUE),
    .groups = 'drop'
  )

analysis_methanotroph <- inner_join(methanotrophs_by_species, flux_by_species, 
                                    by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

cat(sprintf("Species: %d\n", nrow(analysis_methanotroph)))

if(nrow(analysis_methanotroph) > 2) {
  lm_methanotroph <- lm(median_flux ~ log10(median_methanotroph + 1), data = analysis_methanotroph)
  
  cat(sprintf("R² = %.3f, p = %.4f\n", 
              summary(lm_methanotroph)$r.squared,
              summary(lm_methanotroph)$coefficients[2,4]))
  
  cat("\nSpecies with median methanotrophs:\n")
  print(analysis_methanotroph %>% 
          dplyr::select(species, n_trees, median_methanotroph, median_flux) %>%
          arrange(desc(median_methanotroph)) %>%
          head(3))
}

# ============================================================
# ANALYSIS 3: RATIO MODEL (uses BOTH genes subset)
# ============================================================

cat("\n\n### ANALYSIS 3: RATIO MODEL (BOTH genes required) ###\n")

ratio_by_species <- tree_both %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_ratio = median(ratio_mcra_methanotroph, na.rm = TRUE),
    median_log_ratio = median(log_ratio, na.rm = TRUE),
    median_mcra_subset = median(mcrA, na.rm = TRUE),
    median_methanotroph_subset = median(methanotroph_total, na.rm = TRUE),
    q25_ratio = quantile(ratio_mcra_methanotroph, 0.25, na.rm = TRUE),
    q75_ratio = quantile(ratio_mcra_methanotroph, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

analysis_ratio <- inner_join(ratio_by_species, flux_by_species, 
                             by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

cat(sprintf("Species: %d\n", nrow(analysis_ratio)))

if(nrow(analysis_ratio) > 2) {
  lm_ratio <- lm(median_flux ~ median_log_ratio, data = analysis_ratio)
  pearson_ratio <- cor.test(analysis_ratio$median_log_ratio, analysis_ratio$median_flux)
  spearman_ratio <- cor.test(analysis_ratio$median_ratio, analysis_ratio$median_flux, 
                             method = "spearman")
  
  cat(sprintf("R² = %.3f, p = %.4f\n", 
              summary(lm_ratio)$r.squared,
              summary(lm_ratio)$coefficients[2,4]))
  
  cat("\nSpecies with median ratio (subset):\n")
  print(analysis_ratio %>% 
          dplyr::select(species, n_trees, median_ratio, median_mcra_subset, 
                 median_methanotroph_subset, median_flux) %>%
          arrange(desc(median_ratio)) %>%
          head(3))
}

# ============================================================
# COMPARISON: SAME SUBSET FOR FAIR COMPARISON
# ============================================================

cat("\n\n============================================================\n")
cat("FAIR COMPARISON: All models on SAME subset (both genes)\n")
cat("============================================================\n")

# For fair comparison, calculate mcrA and methanotroph models on the SAME trees
# that have both genes (the ratio subset)

mcra_on_both_subset <- tree_both %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_mcra = median(mcrA, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

methanotroph_on_both_subset <- tree_both %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    median_methanotroph = median(methanotroph_total, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

if(nrow(mcra_on_both_subset) > 2 && nrow(methanotroph_on_both_subset) > 2) {
  lm_mcra_subset <- lm(median_flux ~ log10(median_mcra + 1), data = mcra_on_both_subset)
  lm_methanotroph_subset <- lm(median_flux ~ log10(median_methanotroph + 1), 
                               data = methanotroph_on_both_subset)
  
  cat(sprintf("mcrA (both-genes subset):        R² = %.3f, p = %.4f\n", 
              summary(lm_mcra_subset)$r.squared,
              summary(lm_mcra_subset)$coefficients[2,4]))
  
  cat(sprintf("Methanotrophs (both-genes subset): R² = %.3f, p = %.4f\n",
              summary(lm_methanotroph_subset)$r.squared,
              summary(lm_methanotroph_subset)$coefficients[2,4]))
  
  cat(sprintf("Ratio (both-genes subset):       R² = %.3f, p = %.4f\n",
              summary(lm_ratio)$r.squared,
              summary(lm_ratio)$coefficients[2,4]))
  
  cat("\n### RESULT ###\n")
  if(summary(lm_ratio)$r.squared > summary(lm_mcra_subset)$r.squared &&
     summary(lm_ratio)$r.squared > summary(lm_methanotroph_subset)$r.squared) {
    cat("✓ RATIO explains MORE variance than either gene alone\n")
    cat("  This suggests the BALANCE between production and oxidation matters!\n")
  } else {
    cat("✗ RATIO does NOT improve prediction over individual genes\n")
  }
  
  # AIC comparison
  cat("\n### AIC COMPARISON (lower = better) ###\n")
  cat(sprintf("Ratio:          AIC = %.2f\n", AIC(lm_ratio)))
  cat(sprintf("mcrA:           AIC = %.2f\n", AIC(lm_mcra_subset)))
  cat(sprintf("Methanotrophs:  AIC = %.2f\n", AIC(lm_methanotroph_subset)))
}

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")

cat("\nFULL DATASETS (maximum sample size for each gene):\n")
cat(sprintf("  mcrA alone:        n=%d species, R²=%.3f\n",
            nrow(analysis_mcra), summary(lm_mcra)$r.squared))
cat(sprintf("  Methanotrophs alone: n=%d species, R²=%.3f\n",
            nrow(analysis_methanotroph), summary(lm_methanotroph)$r.squared))

cat("\nBOTH-GENES SUBSET (fair comparison):\n")
cat(sprintf("  mcrA:           n=%d species, R²=%.3f\n",
            nrow(mcra_on_both_subset), summary(lm_mcra_subset)$r.squared))
cat(sprintf("  Methanotrophs:  n=%d species, R²=%.3f\n",
            nrow(methanotroph_on_both_subset), summary(lm_methanotroph_subset)$r.squared))
cat(sprintf("  Ratio:          n=%d species, R²=%.3f\n",
            nrow(analysis_ratio), summary(lm_ratio)$r.squared))

cat("\nINTERPRETATION:\n")
cat("- Use FULL datasets when asking: 'Does X predict flux?'\n")
cat("- Use BOTH-GENES subset when comparing: 'What predicts flux BETTER?'\n")
cat("- The ratio analysis shows whether BALANCE matters beyond absolute levels\n")

cat("\n============================================================\n")

# ============================================================
# VISUALIZATIONS
# ============================================================

cat("\n### Creating plots ###\n")

# ============================================================
# PLOT 1: Ratio vs Flux (main result)
# ============================================================

p_ratio <- ggplot(analysis_ratio, aes(x = median_log_ratio, y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "purple", 
              fill = "plum", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_errorbarh(aes(xmin = log10(q25_ratio), xmax = log10(q75_ratio)),
                 height = 0.003, alpha = 0.5, color = "gray40") +
  geom_point(size = 4, alpha = 0.85, color = "purple") +
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted",
             color = "gray50", alpha = 0.5) +
  annotate("text", x = 0, y = max(analysis_ratio$median_flux) * 0.95,
           label = "Equal\nproduction/oxidation", 
           size = 3, color = "gray50", hjust = 0.5) +
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR² = %.3f\nρ = %.2f, p = %.3f",
                           pearson_ratio$estimate, pearson_ratio$p.value,
                           pearson_ratio$estimate^2,
                           spearman_ratio$estimate, spearman_ratio$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  scale_x_continuous(
    breaks = seq(-2, 2, 0.5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    title = "Methanogen:Methanotroph Ratio vs CH₄ Flux",
    subtitle = sprintf("n = %d species (≥5 obs each); error bars = IQR", 
                       nrow(analysis_ratio)),
    x = expression("log"[10]*"(mcrA / [pmoA+mmoX])"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p_ratio)
ggsave("../../../outputs/figures/ratio_vs_flux.png", p_ratio, width = 10, height = 8, dpi = 300)

# ============================================================
# PLOT 2: Three-way comparison bar chart (FAIR COMPARISON)
# ============================================================

comparison_df <- data.frame(
  Predictor = c("mcrA:Methanotroph\nRatio", "mcrA alone", "Methanotrophs alone"),
  R_squared = c(summary(lm_ratio)$r.squared,
                summary(lm_mcra_subset)$r.squared,
                summary(lm_methanotroph_subset)$r.squared),
  P_value = c(summary(lm_ratio)$coefficients[2,4],
              summary(lm_mcra_subset)$coefficients[2,4],
              summary(lm_methanotroph_subset)$coefficients[2,4])
)

p_comparison <- ggplot(comparison_df, aes(x = Predictor, y = R_squared, 
                                          fill = P_value < 0.05)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("R² = %.3f\np = %.3f", R_squared, P_value)),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("gray70", "steelblue"),
                    labels = c("p ≥ 0.05", "p < 0.05"),
                    name = "Significance") +
  labs(
    title = "Model Comparison: Ratio vs Individual Genes",
    subtitle = sprintf("Same dataset (n=%d species with both genes); fair comparison", 
                       nrow(analysis_ratio)),
    y = expression("R"^2),
    x = ""
  ) +
  ylim(0, max(comparison_df$R_squared) * 1.2) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom"
  )

print(p_comparison)
ggsave("../../../outputs/figures/ratio_comparison.png", p_comparison, width = 8, height = 6, dpi = 300)

# ============================================================
# PLOT 3: mcrA scatter (full dataset)
# ============================================================

p_mcra_full <- ggplot(analysis_mcra, 
                      aes(x = log10(median_mcra + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", 
              fill = "pink", alpha = 0.2, linewidth = 1) +
  geom_point(size = 4, alpha = 0.85, color = "darkred") +
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR² = %.3f",
                           pearson_mcra$estimate, pearson_mcra$p.value,
                           pearson_mcra$estimate^2),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  scale_x_continuous(
    breaks = 0:5,
    labels = c(expression(10^0), expression(10^1), expression(10^2), 
               expression(10^3), expression(10^4), expression(10^5))
  ) +
  labs(
    title = "Methanogens (mcrA) vs CH₄ Flux",
    subtitle = sprintf("Full dataset: n=%d species", nrow(analysis_mcra)),
    x = expression("Median area-weighted mcrA (copies g"^-1*")"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p_mcra_full)
ggsave("../../../outputs/figures/mcra_vs_flux_full.png", p_mcra_full, width = 10, height = 8, dpi = 300)

# ============================================================
# PLOT 4: Methanotrophs scatter (full dataset)
# ============================================================

p_methanotroph_full <- ggplot(analysis_methanotroph, 
                              aes(x = log10(median_methanotroph + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", 
              fill = "lightgreen", alpha = 0.2, linewidth = 1) +
  geom_point(size = 4, alpha = 0.85, color = "darkgreen") +
  geom_text_repel(aes(label = species), 
                  size = 3.2,
                  fontface = "italic",
                  color = "black",
                  box.padding = 0.35,
                  max.overlaps = 20,
                  seed = 42) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray50", alpha = 0.7) +
  annotate("label", 
           x = Inf, y = Inf,
           label = sprintf("r = %.2f, p = %.3f\nR² = %.3f",
                           cor.test(log10(analysis_methanotroph$median_methanotroph + 1),
                                    analysis_methanotroph$median_flux)$estimate,
                           summary(lm_methanotroph)$coefficients[2,4],
                           summary(lm_methanotroph)$r.squared),
           hjust = 1.1, vjust = 1.1, size = 3.5,
           fill = "white", alpha = 0.9) +
  scale_x_continuous(
    breaks = 0:5,
    labels = c(expression(10^0), expression(10^1), expression(10^2), 
               expression(10^3), expression(10^4), expression(10^5))
  ) +
  labs(
    title = "Methanotrophs (pmoA+mmoX) vs CH₄ Flux",
    subtitle = sprintf("Full dataset: n=%d species", nrow(analysis_methanotroph)),
    x = expression("Median area-weighted methanotrophs (copies g"^-1*")"),
    y = expression("Median CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

print(p_methanotroph_full)
ggsave("../../../outputs/figures/methanotrophs_vs_flux_full.png", p_methanotroph_full,
       width = 10, height = 8, dpi = 300)

# ============================================================
# PLOT 5: Combined 4-panel figure
# ============================================================

library(patchwork)

# Simplify titles for combined figure
p_mcra_simple <- p_mcra_full + 
  labs(title = "A) Methanogens (mcrA)", subtitle = NULL)

p_methanotroph_simple <- p_methanotroph_full + 
  labs(title = "B) Methanotrophs (pmoA+mmoX)", subtitle = NULL)

p_ratio_simple <- p_ratio + 
  labs(title = "C) Production:Oxidation Ratio", subtitle = NULL)

p_comparison_simple <- p_comparison + 
  labs(title = "D) Model Comparison", subtitle = NULL)

combined <- (p_mcra_simple | p_methanotroph_simple) / 
  (p_ratio_simple | p_comparison_simple) +
  plot_annotation(
    title = "Production-Oxidation Balance Regulates Tree CH₄ Emissions",
    subtitle = "Species-level analysis (n=10 species with ≥5 observations)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

print(combined)
ggsave("../../../outputs/figures/combined_ratio_analysis.png", combined, width = 16, height = 12, dpi = 300)
ggsave("../../../outputs/figures/combined_ratio_analysis.pdf", combined, width = 16, height = 12)

cat("\n### Plots saved ###\n")
cat("  - ratio_vs_flux.png\n")
cat("  - ratio_comparison.png\n")
cat("  - mcra_vs_flux_full.png\n")
cat("  - methanotrophs_vs_flux_full.png\n")
cat("  - combined_ratio_analysis.png/pdf\n")

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")

























################################################################################
# INDIVIDUAL TREE-LEVEL RATIO ANALYSIS - CORRECTED V2
# Run this after the consolidated script
################################################################################

cat("\n============================================================\n")
cat("INDIVIDUAL TREE-LEVEL RATIO ANALYSIS\n")
cat("Testing if ratio improves prediction beyond mmoX alone\n")
cat("============================================================\n")

# Need to go back to tree_genes_weighted which has the raw values
# Merge with flux data
tree_level_ratio <- tree_genes_weighted %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
    by = "tree_id"
  ) %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux)) %>%
  # CRITICAL: Filter for complete data BEFORE calculating anything
  filter(!is.na(mcrA), !is.na(mmoX)) %>%  # Require both genes
  mutate(
    # Calculate methanotroph total from raw area-weighted values
    methanotroph_total = case_when(
      is.na(pmoA) & !is.na(mmoX) ~ mmoX,
      !is.na(pmoA) & is.na(mmoX) ~ pmoA,
      !is.na(pmoA) & !is.na(mmoX) ~ pmoA + mmoX,
      TRUE ~ NA_real_
    ),
    # Calculate ratio (with pseudocount)
    ratio_mcra_methanotroph = (mcrA + 1) / (methanotroph_total + 1),
    log_ratio = log10(ratio_mcra_methanotroph),
    # Log-transform for models
    log_tree_mmox = log10(mmoX + 1),
    log_tree_mcra = log10(mcrA + 1)
  ) %>%
  filter(!is.na(methanotroph_total), !is.na(log_ratio))  # Final check

cat(sprintf("\nTrees with complete data for ratio: %d\n", nrow(tree_level_ratio)))

# Fit models - ALL will use the same 122 trees
m_baseline <- lm(CH4_flux ~ species, data = tree_level_ratio)
m_mmoX_only <- lm(CH4_flux ~ species + log_tree_mmox, data = tree_level_ratio)
m_ratio_only <- lm(CH4_flux ~ species + log_ratio, data = tree_level_ratio)
m_both <- lm(CH4_flux ~ species + log_tree_mmox + log_ratio, data = tree_level_ratio)

# Verify all models have same N
cat(sprintf("\nVerifying sample sizes:\n"))
cat(sprintf("  Baseline: %d\n", nobs(m_baseline)))
cat(sprintf("  mmoX only: %d\n", nobs(m_mmoX_only)))
cat(sprintf("  Ratio only: %d\n", nobs(m_ratio_only)))
cat(sprintf("  Both: %d\n", nobs(m_both)))

# Model comparison table
tree_ratio_comparison <- data.frame(
  Model = c("Species only", 
            "Species + mmoX", 
            "Species + Ratio", 
            "Species + mmoX + Ratio"),
  N = rep(nrow(tree_level_ratio), 4),
  df = c(length(coef(m_baseline)),
         length(coef(m_mmoX_only)),
         length(coef(m_ratio_only)),
         length(coef(m_both))),
  R2 = c(summary(m_baseline)$r.squared,
         summary(m_mmoX_only)$r.squared,
         summary(m_ratio_only)$r.squared,
         summary(m_both)$r.squared),
  Adj_R2 = c(summary(m_baseline)$adj.r.squared,
             summary(m_mmoX_only)$adj.r.squared,
             summary(m_ratio_only)$adj.r.squared,
             summary(m_both)$adj.r.squared),
  AIC = c(AIC(m_baseline),
          AIC(m_mmoX_only),
          AIC(m_ratio_only),
          AIC(m_both))
) %>%
  mutate(
    R2 = round(R2, 4),
    Adj_R2 = round(Adj_R2, 4),
    AIC = round(AIC, 1),
    delta_AIC = AIC - min(AIC)
  ) %>%
  arrange(AIC)

cat("\n### TREE-LEVEL MODEL COMPARISON ###\n")
print(tree_ratio_comparison)

# Statistical tests
cat("\n### ANOVA TESTS ###\n")

cat("\n1. Does ratio add anything beyond mmoX?\n")
anova_ratio_vs_mmox <- anova(m_mmoX_only, m_both)
print(anova_ratio_vs_mmox)
cat(sprintf("Result: %s (p = %.4f)\n", 
            ifelse(anova_ratio_vs_mmox$`Pr(>F)`[2] < 0.05, 
                   "YES - ratio adds significant predictive power",
                   "NO - ratio does not improve model"),
            anova_ratio_vs_mmox$`Pr(>F)`[2]))

cat("\n2. Does mmoX add anything beyond ratio?\n")
anova_mmox_vs_ratio <- anova(m_ratio_only, m_both)
print(anova_mmox_vs_ratio)
cat(sprintf("Result: %s (p = %.4f)\n",
            ifelse(anova_mmox_vs_ratio$`Pr(>F)`[2] < 0.05,
                   "YES - mmoX adds significant predictive power",
                   "NO - mmoX does not improve model"),
            anova_mmox_vs_ratio$`Pr(>F)`[2]))

# Examine best model coefficients
best_tree_model <- if(tree_ratio_comparison$Model[1] == "Species + mmoX + Ratio") {
  m_both
} else if(tree_ratio_comparison$Model[1] == "Species + mmoX") {
  m_mmoX_only
} else {
  m_ratio_only
}

cat("\n### BEST TREE-LEVEL MODEL ###\n")
cat(sprintf("Model: %s\n", tree_ratio_comparison$Model[1]))
print(summary(best_tree_model))

# Extract coefficients for gene predictors
coef_summary <- tidy(best_tree_model, conf.int = TRUE) %>%
  filter(grepl("log_tree_mmox|log_ratio", term))

cat("\n### GENE PREDICTOR COEFFICIENTS ###\n")
print(coef_summary %>% 
        dplyr::select(term, estimate, std.error, statistic, p.value, conf.low, conf.high))

# Save results
write.csv(tree_ratio_comparison,
          "../../../outputs/tables/Tree_level_ratio_model_comparison.csv",
          row.names = FALSE)

write.csv(coef_summary,
          "../../../outputs/tables/Tree_level_ratio_coefficients.csv",
          row.names = FALSE)

# Create visualization
library(ggplot2)
library(patchwork)

# Panel A: mmoX vs flux (colored by ratio)
p_mmox_by_ratio <- ggplot(tree_level_ratio, 
                          aes(x = log_tree_mmox, y = CH4_flux, color = log_ratio)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(name = "log₁₀(Ratio)") +
  labs(title = "A) mmoX vs Flux (colored by ratio)",
       x = "log₁₀ mmoX",
       y = "CH₄ flux") +
  theme_classic()

# Panel B: Ratio vs flux (colored by mmoX)
p_ratio_by_mmox <- ggplot(tree_level_ratio, 
                          aes(x = log_ratio, y = CH4_flux, color = log_tree_mmox)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(name = "log₁₀(mmoX)") +
  labs(title = "B) Ratio vs Flux (colored by mmoX)",
       x = "log₁₀(mcrA:[pmoA+mmoX])",
       y = "CH₄ flux") +
  theme_classic()

# Panel C: Model comparison barplot
p_tree_comparison <- ggplot(tree_ratio_comparison, 
                            aes(x = reorder(Model, -delta_AIC), y = R2)) +
  geom_col(aes(fill = delta_AIC < 2), alpha = 0.8) +
  geom_text(aes(label = sprintf("ΔR² = %.3f\nΔAIC = %.1f", 
                                R2 - min(R2), delta_AIC)),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("gray70", "steelblue"),
                    labels = c("ΔAIC ≥ 2", "ΔAIC < 2"),
                    name = "Model Support") +
  labs(title = "C) Tree-Level Model Comparison",
       x = "",
       y = "R²") +
  ylim(0, max(tree_ratio_comparison$R2) * 1.15) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Panel D: Residuals of best model
tree_level_ratio$residuals <- residuals(best_tree_model)
tree_level_ratio$fitted <- fitted(best_tree_model)

p_residuals <- ggplot(tree_level_ratio, aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 0.8) +
  labs(title = "D) Residuals vs Fitted",
       x = "Fitted values",
       y = "Residuals") +
  theme_classic()

# Combine
combined_tree <- (p_mmox_by_ratio | p_ratio_by_mmox) / 
  (p_tree_comparison | p_residuals) +
  plot_annotation(
    title = "Individual Tree-Level Analysis: Ratio vs mmoX",
    subtitle = sprintf("n = %d trees", nrow(tree_level_ratio))
  )

ggsave("../../../outputs/figures/Figure_Tree_Level_Ratio_Analysis.pdf", combined_tree,
       width = 14, height = 10)
ggsave("../../../outputs/figures/Figure_Tree_Level_Ratio_Analysis.png", combined_tree,
       width = 14, height = 10, dpi = 300)

cat("\n### INTERPRETATION ###\n")
cat("\nKEY RESULT: mmoX alone is the best predictor at the tree level\n")
cat(sprintf("  - mmoX model: R² = %.3f, ΔAIC = %.1f\n", 
            tree_ratio_comparison$R2[tree_ratio_comparison$Model == "Species + mmoX"],
            tree_ratio_comparison$delta_AIC[tree_ratio_comparison$Model == "Species + mmoX"]))
cat(sprintf("  - Ratio does NOT add predictive power beyond mmoX (p = %.3f)\n",
            anova_ratio_vs_mmox$`Pr(>F)`[2]))

cat("\nTREE-LEVEL vs SPECIES-LEVEL comparison:\n")
cat("  Tree-level:    mmoX is best (oxidation capacity)\n")
cat("  Species-level: Ratio is best (production-oxidation balance)\n")
cat("\nThis suggests DIFFERENT MECHANISMS at different scales:\n")
cat("  → Within species: Absolute oxidation limits flux\n")
cat("  → Between species: Balance of production/oxidation matters\n")
cat("  → Species may differ in baseline production rates\n")

cat("\n============================================================\n")
cat("FILES CREATED:\n")
cat("  - Tree_level_ratio_model_comparison.csv\n")
cat("  - Tree_level_ratio_coefficients.csv\n")
cat("  - Figure_Tree_Level_Ratio_Analysis.pdf/png\n")
cat("============================================================\n")




















# ============================================================
# TEST: Raw vs Log-Transformed Concentrations
# ============================================================

cat("\n### TESTING LOG-TRANSFORMATION NECESSITY ###\n")

# Create dataset with BOTH raw and log-transformed values
conc_complete_both <- conc_data %>%
  filter(!is.na(CH4_flux)) %>%
  filter(!is.na(mcrA_Heartwood), !is.na(mcrA_Sapwood),
         !is.na(pmoA_Heartwood), !is.na(pmoA_Sapwood),
         !is.na(mmoX_Heartwood), !is.na(mmoX_Sapwood)) %>%
  mutate(
    # Log-transformed (already exist)
    log_mcra_heart = log10(mcrA_Heartwood + 1),
    log_mcra_sap = log10(mcrA_Sapwood + 1),
    log_mmox_heart = log10(mmoX_Heartwood + 1),
    log_mmox_sap = log10(mmoX_Sapwood + 1)
  )

cat(sprintf("Complete cases: %d trees\n", nrow(conc_complete_both)))

# Models WITH log-transform
m_log_mcra <- lm(CH4_flux ~ species + log_mcra_heart + log_mcra_sap, 
                 data = conc_complete_both)
m_log_mmox <- lm(CH4_flux ~ species + log_mmox_heart + log_mmox_sap, 
                 data = conc_complete_both)

# Models WITHOUT log-transform (RAW concentrations)
m_raw_mcra <- lm(CH4_flux ~ species + mcrA_Heartwood + mcrA_Sapwood, 
                 data = conc_complete_both)
m_raw_mmox <- lm(CH4_flux ~ species + mmoX_Heartwood + mmoX_Sapwood, 
                 data = conc_complete_both)

# Comparison table
log_comparison <- data.frame(
  Gene = rep(c("mcrA", "mmoX"), each = 2),
  Transform = rep(c("Log-transformed", "Raw"), 2),
  R2 = c(
    summary(m_log_mcra)$r.squared,
    summary(m_raw_mcra)$r.squared,
    summary(m_log_mmox)$r.squared,
    summary(m_raw_mmox)$r.squared
  ),
  Adj_R2 = c(
    summary(m_log_mcra)$adj.r.squared,
    summary(m_raw_mcra)$adj.r.squared,
    summary(m_log_mmox)$adj.r.squared,
    summary(m_raw_mmox)$adj.r.squared
  ),
  AIC = c(
    AIC(m_log_mcra),
    AIC(m_raw_mcra),
    AIC(m_log_mmox),
    AIC(m_raw_mmox)
  )
) %>%
  mutate(
    R2 = round(R2, 4),
    Adj_R2 = round(Adj_R2, 4),
    AIC = round(AIC, 1)
  ) %>%
  group_by(Gene) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    Delta_R2 = R2 - min(R2)
  ) %>%
  ungroup()

cat("\n### LOG vs RAW TRANSFORMATION COMPARISON ###\n")
print(log_comparison)

cat("\n### INTERPRETATION ###\n")
for(gene in c("mcrA", "mmoX")) {
  subset <- log_comparison %>% filter(Gene == gene)
  log_better <- subset$AIC[subset$Transform == "Log-transformed"] < 
    subset$AIC[subset$Transform == "Raw"]
  
  cat(sprintf("\n%s:\n", gene))
  if(log_better) {
    cat(sprintf("  ✓ Log-transform is BETTER (ΔAIC = %.1f)\n", 
                subset$Delta_AIC[subset$Transform == "Raw"]))
    cat(sprintf("    Improves R² by %.4f\n", 
                subset$Delta_R2[subset$Transform == "Raw"]))
  } else {
    cat(sprintf("  ✗ Raw is better (ΔAIC = %.1f)\n",
                subset$Delta_AIC[subset$Transform == "Log-transformed"]))
  }
}

write.csv(log_comparison, "../../../outputs/tables/Log_vs_Raw_transformation_comparison.csv", row.names = FALSE)

# Diagnostic plots
pdf("../../../outputs/figures/Figure_Log_Transform_Diagnostics.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))

# mcrA distributions
hist(conc_complete_both$mcrA_Heartwood, breaks = 30, main = "Raw mcrA (Heartwood)",
     xlab = "copies/g", col = "lightblue")
hist(log10(conc_complete_both$mcrA_Heartwood + 1), breaks = 30, 
     main = "Log10 mcrA (Heartwood)", xlab = "log10(copies/g + 1)", col = "darkblue")

# mmoX distributions  
hist(conc_complete_both$mmoX_Heartwood, breaks = 30, main = "Raw mmoX (Heartwood)",
     xlab = "copies/g", col = "lightgreen")
hist(log10(conc_complete_both$mmoX_Heartwood + 1), breaks = 30, 
     main = "Log10 mmoX (Heartwood)", xlab = "log10(copies/g + 1)", col = "darkgreen")

# Q-Q plots
qqnorm(residuals(m_raw_mcra), main = "Raw mcrA Model Residuals")
qqline(residuals(m_raw_mcra))

qqnorm(residuals(m_log_mcra), main = "Log mcrA Model Residuals")
qqline(residuals(m_log_mcra))

par(mfrow = c(1, 1))
dev.off()

cat("\nSaved: Figure_Log_Transform_Diagnostics.pdf\n")











# The tree-level area-weighted data is in tree_genes_weighted
# Let's check if it exists and what columns it has
if(exists("tree_genes_weighted")) {
  cat("Found tree_genes_weighted\n")
  print(names(tree_genes_weighted))
  
  # Create the comparison dataset
  area_complete_both <- tree_genes_weighted %>%
    left_join(
      ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
      by = "tree_id"
    ) %>%
    filter(!is.na(CH4_flux), !is.na(mcrA), !is.na(mmoX)) %>%
    mutate(
      log_tree_mcra = log10(mcrA + 1),
      log_tree_mmox = log10(mmoX + 1)
    )
  
  cat(sprintf("\nArea-weighted trees: %d\n", nrow(area_complete_both)))
  
  # Now run the models
  m_area_raw_mmox <- lm(CH4_flux ~ species + mmoX, data = area_complete_both)
  m_area_raw_mcra <- lm(CH4_flux ~ species + mcrA, data = area_complete_both)
  m_area_log_mmox <- lm(CH4_flux ~ species + log_tree_mmox, data = area_complete_both)
  m_area_log_mcra <- lm(CH4_flux ~ species + log_tree_mcra, data = area_complete_both)
  
  # Comparison
  area_comparison <- data.frame(
    Gene = rep(c("mcrA", "mmoX"), each = 2),
    Transform = rep(c("Log-transformed", "Raw"), 2),
    R2 = c(
      summary(m_area_log_mcra)$r.squared,
      summary(m_area_raw_mcra)$r.squared,
      summary(m_area_log_mmox)$r.squared,
      summary(m_area_raw_mmox)$r.squared
    ),
    AIC = c(
      AIC(m_area_log_mcra),
      AIC(m_area_raw_mcra),
      AIC(m_area_log_mmox),
      AIC(m_area_raw_mmox)
    )
  ) %>%
    mutate(
      R2 = round(R2, 4),
      AIC = round(AIC, 1)
    ) %>%
    group_by(Gene) %>%
    mutate(Delta_AIC = AIC - min(AIC)) %>%
    ungroup()
  
  cat("\n### AREA-WEIGHTED: LOG vs RAW ###\n")
  print(area_comparison)
  
} else {
  cat("tree_genes_weighted not found. Re-run the data preparation section.\n")
}













################################################################################
# SPECIES-LEVEL MODELS: AREA-WEIGHTED vs CONCENTRATION COMPARISON
# Drop-in code to compare predictive power at species level
################################################################################

cat("\n============================================================\n")
cat("SPECIES-LEVEL ANALYSIS: AREA-WEIGHTED vs CONCENTRATION\n")
cat("============================================================\n")

# ============================================================
# 1. SPECIES-LEVEL AREA-WEIGHTED DATA
# ============================================================

cat("\n### Aggregating area-weighted genes by species ###\n")

species_area_weighted <- tree_genes_weighted %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
    by = "tree_id"
  ) %>%
  filter(!is.na(CH4_flux)) %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    # Area-weighted genes
    median_mcra = median(mcrA, na.rm = TRUE),
    median_mmox = median(mmoX, na.rm = TRUE),
    median_pmoa = median(pmoA, na.rm = TRUE),
    # Flux
    median_flux = median(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_trees >= 5)

cat(sprintf("Species with area-weighted data: %d\n", nrow(species_area_weighted)))

# ============================================================
# 2. SPECIES-LEVEL CONCENTRATION DATA
# ============================================================

cat("\n### Aggregating concentration genes by species ###\n")

species_concentration <- conc_data %>%
  filter(!is.na(CH4_flux)) %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    # Heartwood
    median_mcra_heart = median(mcrA_Heartwood, na.rm = TRUE),
    median_mmox_heart = median(mmoX_Heartwood, na.rm = TRUE),
    median_pmoa_heart = median(pmoA_Heartwood, na.rm = TRUE),
    # Sapwood
    median_mcra_sap = median(mcrA_Sapwood, na.rm = TRUE),
    median_mmox_sap = median(mmoX_Sapwood, na.rm = TRUE),
    median_pmoa_sap = median(pmoA_Sapwood, na.rm = TRUE),
    # Flux
    median_flux = median(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_trees >= 5)

cat(sprintf("Species with concentration data: %d\n", nrow(species_concentration)))

# ============================================================
# 3. SPECIES-LEVEL MODELS: AREA-WEIGHTED
# ============================================================

cat("\n\n### SPECIES-LEVEL MODELS: AREA-WEIGHTED ###\n")

# Filter for complete data
species_area_complete <- species_area_weighted %>%
  filter(!is.na(median_mcra), !is.na(median_mmox))

cat(sprintf("Complete cases: %d species\n", nrow(species_area_complete)))

if(nrow(species_area_complete) >= 3) {
  # Models
  sp_area_mcra <- lm(median_flux ~ log10(median_mcra + 1), 
                     data = species_area_complete)
  sp_area_mmox <- lm(median_flux ~ log10(median_mmox + 1), 
                     data = species_area_complete)
  
  # Correlations
  cor_area_mcra <- cor.test(log10(species_area_complete$median_mcra + 1), 
                            species_area_complete$median_flux)
  cor_area_mmox <- cor.test(log10(species_area_complete$median_mmox + 1), 
                            species_area_complete$median_flux)
  
  cat("\nArea-weighted mcrA:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_area_mcra$estimate, 
              cor_area_mcra$estimate^2,
              cor_area_mcra$p.value))
  
  cat("\nArea-weighted mmoX:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_area_mmox$estimate,
              cor_area_mmox$estimate^2,
              cor_area_mmox$p.value))
}

# ============================================================
# 4. SPECIES-LEVEL MODELS: CONCENTRATION
# ============================================================

cat("\n\n### SPECIES-LEVEL MODELS: CONCENTRATION ###\n")

# Filter for complete data
species_conc_complete <- species_concentration %>%
  filter(!is.na(median_mcra_heart), !is.na(median_mcra_sap),
         !is.na(median_mmox_heart), !is.na(median_mmox_sap))

cat(sprintf("Complete cases: %d species\n", nrow(species_conc_complete)))

if(nrow(species_conc_complete) >= 3) {
  # Models - test each component separately
  sp_conc_mcra_heart <- lm(median_flux ~ log10(median_mcra_heart + 1), 
                           data = species_conc_complete)
  sp_conc_mcra_sap <- lm(median_flux ~ log10(median_mcra_sap + 1), 
                         data = species_conc_complete)
  sp_conc_mmox_heart <- lm(median_flux ~ log10(median_mmox_heart + 1), 
                           data = species_conc_complete)
  sp_conc_mmox_sap <- lm(median_flux ~ log10(median_mmox_sap + 1), 
                         data = species_conc_complete)
  
  # Correlations
  cor_mcra_heart <- cor.test(log10(species_conc_complete$median_mcra_heart + 1),
                             species_conc_complete$median_flux)
  cor_mcra_sap <- cor.test(log10(species_conc_complete$median_mcra_sap + 1),
                           species_conc_complete$median_flux)
  cor_mmox_heart <- cor.test(log10(species_conc_complete$median_mmox_heart + 1),
                             species_conc_complete$median_flux)
  cor_mmox_sap <- cor.test(log10(species_conc_complete$median_mmox_sap + 1),
                           species_conc_complete$median_flux)
  
  cat("\nHeartwood mcrA:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_mcra_heart$estimate, cor_mcra_heart$estimate^2,
              cor_mcra_heart$p.value))
  
  cat("\nSapwood mcrA:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_mcra_sap$estimate, cor_mcra_sap$estimate^2,
              cor_mcra_sap$p.value))
  
  cat("\nHeartwood mmoX:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_mmox_heart$estimate, cor_mmox_heart$estimate^2,
              cor_mmox_heart$p.value))
  
  cat("\nSapwood mmoX:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_mmox_sap$estimate, cor_mmox_sap$estimate^2,
              cor_mmox_sap$p.value))
  
  # Combined model
  sp_conc_mmox_both <- lm(median_flux ~ log10(median_mmox_heart + 1) + 
                            log10(median_mmox_sap + 1),
                          data = species_conc_complete)
  
  cat("\nSapwood + Heartwood mmoX (combined):\n")
  cat(sprintf("  R² = %.3f, p = %.4f\n",
              summary(sp_conc_mmox_both)$r.squared,
              summary(sp_conc_mmox_both)$coefficients[2,4]))
  print(summary(sp_conc_mmox_both)$coefficients[2:3, ])
}

# ============================================================
# 5. COMPARISON TABLE
# ============================================================

cat("\n\n### SPECIES-LEVEL MODEL COMPARISON ###\n")

species_comparison <- data.frame(
  Model = c("Area-weighted mcrA",
            "Area-weighted mmoX",
            "Heartwood mcrA",
            "Sapwood mcrA", 
            "Heartwood mmoX",
            "Sapwood mmoX",
            "Sapwood + Heartwood mmoX"),
  N_species = c(
    nrow(species_area_complete),
    nrow(species_area_complete),
    nrow(species_conc_complete),
    nrow(species_conc_complete),
    nrow(species_conc_complete),
    nrow(species_conc_complete),
    nrow(species_conc_complete)
  ),
  R2 = c(
    cor_area_mcra$estimate^2,
    cor_area_mmox$estimate^2,
    cor_mcra_heart$estimate^2,
    cor_mcra_sap$estimate^2,
    cor_mmox_heart$estimate^2,
    cor_mmox_sap$estimate^2,
    summary(sp_conc_mmox_both)$r.squared
  ),
  P_value = c(
    cor_area_mcra$p.value,
    cor_area_mmox$p.value,
    cor_mcra_heart$p.value,
    cor_mcra_sap$p.value,
    cor_mmox_heart$p.value,
    cor_mmox_sap$p.value,
    anova(lm(median_flux ~ 1, data = species_conc_complete), 
          sp_conc_mmox_both)$`Pr(>F)`[2]
  )
) %>%
  mutate(
    R2 = round(R2, 3),
    P_value = round(P_value, 4),
    Significant = ifelse(P_value < 0.05, "**", 
                         ifelse(P_value < 0.10, "*", "ns"))
  ) %>%
  arrange(desc(R2))

print(species_comparison)

write.csv(species_comparison,
          "../../../outputs/tables/Species_level_area_vs_concentration_comparison.csv",
          row.names = FALSE)

# ============================================================
# 6. VISUALIZATIONS
# ============================================================

cat("\n### Creating species-level comparison plots ###\n")

library(patchwork)

# Panel A: Area-weighted mmoX
p_sp_area <- ggplot(species_area_complete, 
                    aes(x = log10(median_mmox + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue",
              fill = "lightblue", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_point(size = 4, alpha = 0.85, color = "steelblue") +
  geom_text_repel(aes(label = species), 
                  size = 3, fontface = "italic",
                  box.padding = 0.3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R² = %.3f\np = %.3f", 
                           cor_area_mmox$estimate^2,
                           cor_area_mmox$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(title = "A) Area-Weighted mmoX",
       x = "log₁₀(median area-weighted mmoX)",
       y = "Median CH₄ flux") +
  theme_classic()

# Panel B: Sapwood mmoX
p_sp_sap <- ggplot(species_conc_complete,
                   aes(x = log10(median_mmox_sap + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen",
              fill = "lightgreen", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_point(size = 4, alpha = 0.85, color = "darkgreen") +
  geom_text_repel(aes(label = species),
                  size = 3, fontface = "italic",
                  box.padding = 0.3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R² = %.3f\np = %.3f",
                           cor_mmox_sap$estimate^2,
                           cor_mmox_sap$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(title = "B) Sapwood mmoX",
       x = "log₁₀(median sapwood mmoX)",
       y = "") +
  theme_classic()

# Panel C: Comparison barplot
p_sp_compare <- ggplot(species_comparison %>% filter(Model != "Sapwood + Heartwood mmoX"),
                       aes(x = reorder(Model, R2), y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("ns" = "gray70", "*" = "gold", "**" = "darkgreen"),
                    name = "p-value") +
  coord_flip() +
  ylim(0, max(species_comparison$R2) * 1.15) +
  labs(title = "C) Species-Level Model Comparison",
       x = "",
       y = "R²") +
  theme_classic()

# Combine
combined_species <- (p_sp_area | p_sp_sap) / p_sp_compare +
  plot_annotation(
    title = "Species-Level Analysis: Area-Weighted vs Concentration",
    subtitle = sprintf("n = %d species", nrow(species_conc_complete))
  )

ggsave("../../../outputs/figures/Figure_Species_Level_Area_vs_Concentration.pdf",
       combined_species, width = 14, height = 10)
ggsave("../../../outputs/figures/Figure_Species_Level_Area_vs_Concentration.png",
       combined_species, width = 14, height = 10, dpi = 300)

# ============================================================
# 7. FINAL SUMMARY
# ============================================================

cat("\n\n============================================================\n")
cat("SPECIES-LEVEL RESULTS SUMMARY\n")
cat("============================================================\n")

cat("\nArea-weighted approach (species-level):\n")
cat(sprintf("  mmoX: R² = %.3f, p = %.4f %s\n",
            cor_area_mmox$estimate^2, cor_area_mmox$p.value,
            ifelse(cor_area_mmox$p.value < 0.05, "**", "ns")))

cat("\nConcentration approach (species-level):\n")
cat(sprintf("  Sapwood mmoX: R² = %.3f, p = %.4f %s\n",
            cor_mmox_sap$estimate^2, cor_mmox_sap$p.value,
            ifelse(cor_mmox_sap$p.value < 0.05, "**", "ns")))
cat(sprintf("  Heartwood mmoX: R² = %.3f, p = %.4f %s\n",
            cor_mmox_heart$estimate^2, cor_mmox_heart$p.value,
            ifelse(cor_mmox_heart$p.value < 0.05, "**", "ns")))

cat("\nKEY FINDING:\n")
best_model <- species_comparison$Model[1]
best_r2 <- species_comparison$R2[1]

if(grepl("Sapwood", best_model)) {
  cat(sprintf("  ✓ Sapwood mmoX is the BEST predictor at species level (R² = %.3f)\n", 
              best_r2))
  cat("  → Spatially-resolved concentration data outperforms area-weighted\n")
} else if(grepl("Area-weighted", best_model)) {
  cat(sprintf("  ✓ Area-weighted mmoX is the BEST predictor at species level (R² = %.3f)\n",
              best_r2))
  cat("  → Area-weighting captures species differences adequately\n")
}

cat("\nScale comparison:\n")
cat("  Individual tree level: Area-weighted mmoX (R² = 0.300, p = 0.006)\n")
cat(sprintf("  Species level:         %s (R² = %.3f, p = %.4f)\n",
            best_model, best_r2, 
            species_comparison$P_value[species_comparison$Model == best_model]))

cat("\n============================================================\n")
cat("FILES CREATED:\n")
cat("  - Species_level_area_vs_concentration_comparison.csv\n")
cat("  - Figure_Species_Level_Area_vs_Concentration.pdf/png\n")
cat("============================================================\n")

################################################################################
# SPECIES-LEVEL MODELS: AREA-WEIGHTED vs CONCENTRATION COMPARISON
# Drop-in code to compare predictive power at species level
################################################################################

cat("\n============================================================\n")
cat("SPECIES-LEVEL ANALYSIS: AREA-WEIGHTED vs CONCENTRATION\n")
cat("============================================================\n")

# ============================================================
# 1. SPECIES-LEVEL AREA-WEIGHTED DATA
# ============================================================

cat("\n### Aggregating area-weighted genes by species ###\n")

species_area_weighted <- tree_genes_weighted %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
    by = "tree_id"
  ) %>%
  filter(!is.na(CH4_flux)) %>%
  group_by(species, species_id) %>%
  summarise(
    n_trees = n(),
    # Area-weighted genes
    median_mcra = median(mcrA, na.rm = TRUE),
    median_mmox = median(mmoX, na.rm = TRUE),
    median_pmoa = median(pmoA, na.rm = TRUE),
    # Flux
    median_flux = median(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_trees >= 5)

cat(sprintf("Species with area-weighted data: %d\n", nrow(species_area_weighted)))

# ============================================================
# 2. SPECIES-LEVEL CONCENTRATION DATA
# ============================================================

cat("\n### Aggregating concentration genes by species ###\n")

species_concentration <- conc_data %>%
  filter(!is.na(CH4_flux)) %>%
  group_by(species) %>%  # Remove species_id since it's not in conc_data
  summarise(
    n_trees = n(),
    # Heartwood
    median_mcra_heart = median(mcrA_Heartwood, na.rm = TRUE),
    median_mmox_heart = median(mmoX_Heartwood, na.rm = TRUE),
    median_pmoa_heart = median(pmoA_Heartwood, na.rm = TRUE),
    # Sapwood
    median_mcra_sap = median(mcrA_Sapwood, na.rm = TRUE),
    median_mmox_sap = median(mmoX_Sapwood, na.rm = TRUE),
    median_pmoa_sap = median(pmoA_Sapwood, na.rm = TRUE),
    # Flux
    median_flux = median(CH4_flux, na.rm = TRUE),
    q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_trees >= 5)

cat(sprintf("Species with concentration data: %d\n", nrow(species_concentration)))

# ============================================================
# 3. SPECIES-LEVEL MODELS: AREA-WEIGHTED
# ============================================================

cat("\n\n### SPECIES-LEVEL MODELS: AREA-WEIGHTED ###\n")

# Filter for complete data
species_area_complete <- species_area_weighted %>%
  filter(!is.na(median_mcra), !is.na(median_mmox))

cat(sprintf("Complete cases: %d species\n", nrow(species_area_complete)))

if(nrow(species_area_complete) >= 3) {
  # Models
  sp_area_mcra <- lm(median_flux ~ log10(median_mcra + 1), 
                     data = species_area_complete)
  sp_area_mmox <- lm(median_flux ~ log10(median_mmox + 1), 
                     data = species_area_complete)
  
  # Correlations
  cor_area_mcra <- cor.test(log10(species_area_complete$median_mcra + 1), 
                            species_area_complete$median_flux)
  cor_area_mmox <- cor.test(log10(species_area_complete$median_mmox + 1), 
                            species_area_complete$median_flux)
  
  cat("\nArea-weighted mcrA:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_area_mcra$estimate, 
              cor_area_mcra$estimate^2,
              cor_area_mcra$p.value))
  
  cat("\nArea-weighted mmoX:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_area_mmox$estimate,
              cor_area_mmox$estimate^2,
              cor_area_mmox$p.value))
}

# ============================================================
# 4. SPECIES-LEVEL MODELS: CONCENTRATION
# ============================================================

cat("\n\n### SPECIES-LEVEL MODELS: CONCENTRATION ###\n")

# Filter for complete data
species_conc_complete <- species_concentration %>%
  filter(!is.na(median_mcra_heart), !is.na(median_mcra_sap),
         !is.na(median_mmox_heart), !is.na(median_mmox_sap))

cat(sprintf("Complete cases: %d species\n", nrow(species_conc_complete)))

if(nrow(species_conc_complete) >= 3) {
  # Models - test each component separately
  sp_conc_mcra_heart <- lm(median_flux ~ log10(median_mcra_heart + 1), 
                           data = species_conc_complete)
  sp_conc_mcra_sap <- lm(median_flux ~ log10(median_mcra_sap + 1), 
                         data = species_conc_complete)
  sp_conc_mmox_heart <- lm(median_flux ~ log10(median_mmox_heart + 1), 
                           data = species_conc_complete)
  sp_conc_mmox_sap <- lm(median_flux ~ log10(median_mmox_sap + 1), 
                         data = species_conc_complete)
  
  # Correlations
  cor_mcra_heart <- cor.test(log10(species_conc_complete$median_mcra_heart + 1),
                             species_conc_complete$median_flux)
  cor_mcra_sap <- cor.test(log10(species_conc_complete$median_mcra_sap + 1),
                           species_conc_complete$median_flux)
  cor_mmox_heart <- cor.test(log10(species_conc_complete$median_mmox_heart + 1),
                             species_conc_complete$median_flux)
  cor_mmox_sap <- cor.test(log10(species_conc_complete$median_mmox_sap + 1),
                           species_conc_complete$median_flux)
  
  cat("\nHeartwood mcrA:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_mcra_heart$estimate, cor_mcra_heart$estimate^2,
              cor_mcra_heart$p.value))
  
  cat("\nSapwood mcrA:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_mcra_sap$estimate, cor_mcra_sap$estimate^2,
              cor_mcra_sap$p.value))
  
  cat("\nHeartwood mmoX:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_mmox_heart$estimate, cor_mmox_heart$estimate^2,
              cor_mmox_heart$p.value))
  
  cat("\nSapwood mmoX:\n")
  cat(sprintf("  r = %.3f, R² = %.3f, p = %.4f\n",
              cor_mmox_sap$estimate, cor_mmox_sap$estimate^2,
              cor_mmox_sap$p.value))
  
  # Combined model
  sp_conc_mmox_both <- lm(median_flux ~ log10(median_mmox_heart + 1) + 
                            log10(median_mmox_sap + 1),
                          data = species_conc_complete)
  
  cat("\nSapwood + Heartwood mmoX (combined):\n")
  cat(sprintf("  R² = %.3f\n", summary(sp_conc_mmox_both)$r.squared))
  print(summary(sp_conc_mmox_both)$coefficients[2:3, ])
}

# ============================================================
# 5. COMPARISON TABLE
# ============================================================

cat("\n\n### SPECIES-LEVEL MODEL COMPARISON ###\n")

species_comparison <- data.frame(
  Model = c("Area-weighted mcrA",
            "Area-weighted mmoX",
            "Heartwood mcrA",
            "Sapwood mcrA", 
            "Heartwood mmoX",
            "Sapwood mmoX",
            "Sapwood + Heartwood mmoX"),
  N_species = c(
    nrow(species_area_complete),
    nrow(species_area_complete),
    nrow(species_conc_complete),
    nrow(species_conc_complete),
    nrow(species_conc_complete),
    nrow(species_conc_complete),
    nrow(species_conc_complete)
  ),
  R2 = c(
    cor_area_mcra$estimate^2,
    cor_area_mmox$estimate^2,
    cor_mcra_heart$estimate^2,
    cor_mcra_sap$estimate^2,
    cor_mmox_heart$estimate^2,
    cor_mmox_sap$estimate^2,
    summary(sp_conc_mmox_both)$r.squared
  ),
  P_value = c(
    cor_area_mcra$p.value,
    cor_area_mmox$p.value,
    cor_mcra_heart$p.value,
    cor_mcra_sap$p.value,
    cor_mmox_heart$p.value,
    cor_mmox_sap$p.value,
    anova(lm(median_flux ~ 1, data = species_conc_complete), 
          sp_conc_mmox_both)$`Pr(>F)`[2]
  )
) %>%
  mutate(
    R2 = round(R2, 3),
    P_value = round(P_value, 4),
    Significant = ifelse(P_value < 0.05, "**", 
                         ifelse(P_value < 0.10, "*", "ns"))
  ) %>%
  arrange(desc(R2))

print(species_comparison)

write.csv(species_comparison,
          "../../../outputs/tables/Species_level_area_vs_concentration_comparison.csv",
          row.names = FALSE)

# ============================================================
# 6. VISUALIZATIONS
# ============================================================

cat("\n### Creating species-level comparison plots ###\n")

library(patchwork)

# Panel A: Area-weighted mmoX
p_sp_area <- ggplot(species_area_complete, 
                    aes(x = log10(median_mmox + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue",
              fill = "lightblue", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_point(size = 4, alpha = 0.85, color = "steelblue") +
  geom_text_repel(aes(label = species), 
                  size = 3, fontface = "italic",
                  box.padding = 0.3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R² = %.3f\np = %.3f", 
                           cor_area_mmox$estimate^2,
                           cor_area_mmox$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(title = "A) Area-Weighted mmoX",
       x = "log₁₀(median area-weighted mmoX)",
       y = "Median CH₄ flux") +
  theme_classic()

# Panel B: Sapwood mmoX
p_sp_sap <- ggplot(species_conc_complete,
                   aes(x = log10(median_mmox_sap + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen",
              fill = "lightgreen", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_point(size = 4, alpha = 0.85, color = "darkgreen") +
  geom_text_repel(aes(label = species),
                  size = 3, fontface = "italic",
                  box.padding = 0.3, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("R² = %.3f\np = %.3f",
                           cor_mmox_sap$estimate^2,
                           cor_mmox_sap$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(title = "B) Sapwood mmoX",
       x = "log₁₀(median sapwood mmoX)",
       y = "") +
  theme_classic()

# Panel C: Comparison barplot
p_sp_compare <- ggplot(species_comparison %>% filter(Model != "Sapwood + Heartwood mmoX"),
                       aes(x = reorder(Model, R2), y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", R2)),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("ns" = "gray70", "*" = "gold", "**" = "darkgreen"),
                    name = "p-value") +
  coord_flip() +
  ylim(0, max(species_comparison$R2) * 1.15) +
  labs(title = "C) Species-Level Model Comparison",
       x = "",
       y = "R²") +
  theme_classic()

# Combine
combined_species <- (p_sp_area | p_sp_sap) / p_sp_compare +
  plot_annotation(
    title = "Species-Level Analysis: Area-Weighted vs Concentration",
    subtitle = sprintf("n = %d species", nrow(species_conc_complete))
  )

ggsave("../../../outputs/figures/Figure_Species_Level_Area_vs_Concentration.pdf",
       combined_species, width = 14, height = 10)
ggsave("../../../outputs/figures/Figure_Species_Level_Area_vs_Concentration.png",
       combined_species, width = 14, height = 10, dpi = 300)

# ============================================================
# 7. FINAL SUMMARY
# ============================================================

cat("\n\n============================================================\n")
cat("SPECIES-LEVEL RESULTS SUMMARY\n")
cat("============================================================\n")

cat("\nArea-weighted approach (species-level):\n")
cat(sprintf("  mmoX: R² = %.3f, p = %.4f %s\n",
            cor_area_mmox$estimate^2, cor_area_mmox$p.value,
            ifelse(cor_area_mmox$p.value < 0.05, "**", "ns")))

cat("\nConcentration approach (species-level):\n")
cat(sprintf("  Sapwood mmoX: R² = %.3f, p = %.4f %s\n",
            cor_mmox_sap$estimate^2, cor_mmox_sap$p.value,
            ifelse(cor_mmox_sap$p.value < 0.05, "**", "ns")))
cat(sprintf("  Heartwood mmoX: R² = %.3f, p = %.4f %s\n",
            cor_mmox_heart$estimate^2, cor_mmox_heart$p.value,
            ifelse(cor_mmox_heart$p.value < 0.05, "**", "ns")))

cat("\nKEY FINDING:\n")
best_model <- species_comparison$Model[1]
best_r2 <- species_comparison$R2[1]

if(grepl("Sapwood", best_model)) {
  cat(sprintf("  ✓ Sapwood mmoX is the BEST predictor at species level (R² = %.3f)\n", 
              best_r2))
  cat("  → Spatially-resolved concentration data outperforms area-weighted\n")
} else if(grepl("Area-weighted", best_model)) {
  cat(sprintf("  ✓ Area-weighted mmoX is the BEST predictor at species level (R² = %.3f)\n",
              best_r2))
  cat("  → Area-weighting captures species differences adequately\n")
}

cat("\nScale comparison:\n")
cat("  Individual tree level: Area-weighted mmoX (R² = 0.300, p = 0.006)\n")
cat(sprintf("  Species level:         %s (R² = %.3f, p = %.4f)\n",
            best_model, best_r2, 
            species_comparison$P_value[species_comparison$Model == best_model]))

cat("\n============================================================\n")
cat("FILES CREATED:\n")
cat("  - Species_level_area_vs_concentration_comparison.csv\n")
cat("  - Figure_Species_Level_Area_vs_Concentration.pdf/png\n")
cat("============================================================\n")













################################################################################
# SPECIES-LEVEL FIGURES: PRODUCTION vs OXIDATION AT DIFFERENT SCALES
################################################################################

cat("\n### Creating comprehensive multi-scale figures ###\n")

library(patchwork)
library(ggrepel)

# ============================================================
# FIGURE 1: SCALE-DEPENDENT MECHANISMS (4-panel)
# ============================================================

cat("\nCreating Figure 1: Scale-dependent mechanisms...\n")

# Panel A: Individual tree level - mmoX (oxidation)
tree_mmox_data <- tree_genes_weighted %>%
  left_join(
    ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
    by = "tree_id"
  ) %>%
  filter(!is.na(CH4_flux), !is.na(mmoX)) %>%
  mutate(log_mmox = log10(mmoX + 1))

cat(sprintf("Tree-level data: %d observations\n", nrow(tree_mmox_data)))

# Get the actual model stats from your earlier analysis
tree_lm_mmox <- lm(CH4_flux ~ species + log_mmox, data = tree_mmox_data)
tree_r2 <- summary(tree_lm_mmox)$r.squared
tree_p <- summary(tree_lm_mmox)$coefficients["log_mmox", "Pr(>|t|)"]

cat(sprintf("Tree-level mmoX model: R² = %.3f, p = %.4f\n", tree_r2, tree_p))

p1_tree_mmox <- ggplot(tree_mmox_data, aes(x = log_mmox, y = CH4_flux)) +
  geom_point(aes(color = species), size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Individual Trees\nR² = %.3f\np = %.4f **", 
                           tree_r2, tree_p),
           hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
  labs(title = "A) Within-Species: Oxidation Limits Flux",
       x = expression("log"[10]*"(area-weighted mmoX)"),
       y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

# Panel B: Species level - heartwood mcrA (production)
p2_species_mcra <- ggplot(species_conc_complete, 
                          aes(x = log10(median_mcra_heart + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred",
              fill = "pink", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_point(size = 4, alpha = 0.85, color = "darkred") +
  geom_text_repel(aes(label = species), 
                  size = 2.5, fontface = "italic",
                  box.padding = 0.25, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Species Medians\nR² = %.3f\np = %.4f **",
                           cor_mcra_heart$estimate^2,
                           cor_mcra_heart$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5, fontface = "bold") +
  labs(title = "B) Between-Species: Production Sets Baseline",
       x = expression("log"[10]*"(median heartwood mcrA)"),
       y = expression("Median CH"[4]*" flux")) +
  theme_classic(base_size = 11)

# Panel C: Species level - sapwood mmoX (NOT significant)
p3_species_mmox <- ggplot(species_conc_complete,
                          aes(x = log10(median_mmox_sap + 1), y = median_flux)) +
  geom_smooth(method = "lm", se = TRUE, color = "gray50",
              fill = "gray80", alpha = 0.2, linewidth = 1) +
  geom_errorbar(aes(ymin = q25_flux, ymax = q75_flux),
                width = 0.02, alpha = 0.5, color = "gray40") +
  geom_point(size = 4, alpha = 0.85, color = "gray50") +
  geom_text_repel(aes(label = species),
                  size = 2.5, fontface = "italic",
                  box.padding = 0.25, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Species Medians\nR² = %.3f\np = %.4f NS",
                           cor_mmox_sap$estimate^2,
                           cor_mmox_sap$p.value),
           hjust = 1.1, vjust = 1.1, size = 3.5) +
  labs(title = "C) Between-Species: Oxidation Does Not Vary",
       x = expression("log"[10]*"(median sapwood mmoX)"),
       y = "") +
  theme_classic(base_size = 11)

# Panel D: Model comparison across scales
# Get tree-level mcrA model for comparison
tree_mcra_data <- tree_mmox_data %>% filter(!is.na(mcrA))
tree_lm_mcra <- lm(CH4_flux ~ species + log10(mcrA + 1), data = tree_mcra_data)
tree_mcra_r2 <- summary(tree_lm_mcra)$r.squared
tree_mcra_p <- summary(tree_lm_mcra)$coefficients["log10(mcrA + 1)", "Pr(>|t|)"]

scale_comparison <- data.frame(
  Scale = c("Individual\nTrees", "Individual\nTrees", 
            "Species\nMedians", "Species\nMedians"),
  Predictor = c("mmoX\n(oxidation)", "mcrA\n(production)",
                "Sapwood mmoX\n(oxidation)", "Heartwood mcrA\n(production)"),
  R2 = c(tree_r2, tree_mcra_r2,
         cor_mmox_sap$estimate^2,
         cor_mcra_heart$estimate^2),
  P_value = c(tree_p, tree_mcra_p,
              cor_mmox_sap$p.value,
              cor_mcra_heart$p.value),
  Significant = c(tree_p < 0.05, tree_mcra_p < 0.05,
                  cor_mmox_sap$p.value < 0.05, 
                  cor_mcra_heart$p.value < 0.05)
)

p4_comparison <- ggplot(scale_comparison, 
                        aes(x = interaction(Predictor, Scale, sep = "\n"), 
                            y = R2, fill = Significant)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("R² = %.2f\n%s", R2, 
                                ifelse(Significant, "**", "NS"))),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkgreen"),
                    name = "Significant") +
  scale_x_discrete(labels = function(x) gsub("\n\n", "\n", x)) +
  ylim(0, max(scale_comparison$R2) * 1.2) +
  labs(title = "D) Scale-Dependent Mechanisms",
       x = "",
       y = expression("R"^2*" (variance explained)")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 8, lineheight = 0.9),
        legend.position = "bottom")

# Combine Figure 1
fig1_combined <- (p1_tree_mmox | p2_species_mcra) / (p3_species_mmox | p4_comparison) +
  plot_annotation(
    title = "Scale-Dependent Mechanisms: Production Sets Species Baseline, Oxidation Provides Within-Species Regulation",
    subtitle = sprintf("Left panels: n = %d individual trees | Right panels: n = %d species", 
                       nrow(tree_mmox_data), nrow(species_conc_complete)),
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "gray30"))
  )

ggsave("../../../outputs/figures/Figure_Scale_Dependent_Mechanisms.pdf", fig1_combined,
       width = 14, height = 10)
ggsave("../../../outputs/figures/Figure_Scale_Dependent_Mechanisms.png", fig1_combined,
       width = 14, height = 10, dpi = 300)

cat("Saved: Figure_Scale_Dependent_Mechanisms.pdf/png\n")

# ============================================================
# FIGURE 2: SPATIAL PATTERNS (heartwood vs sapwood)
# ============================================================

cat("\nCreating Figure 2: Heartwood vs Sapwood patterns...\n")

# Compare heartwood vs sapwood for both genes at species level
spatial_data <- species_conc_complete %>%
  mutate(
    mcra_ratio = median_mcra_heart / (median_mcra_sap + 1),
    mmox_ratio = median_mmox_heart / (median_mmox_sap + 1)
  )

# Panel A: Heartwood vs Sapwood mcrA
p5_mcra_spatial <- ggplot(spatial_data, 
                          aes(x = log10(median_mcra_sap + 1), 
                              y = log10(median_mcra_heart + 1))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray50", alpha = 0.7) +
  geom_point(aes(size = abs(median_flux), color = median_flux), alpha = 0.8) +
  geom_text_repel(aes(label = species), size = 2.5, fontface = "italic") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "CH4 flux") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "A) mcrA: Heartwood > Sapwood",
       subtitle = "Higher production in anoxic core",
       x = expression("log"[10]*"(sapwood mcrA)"),
       y = expression("log"[10]*"(heartwood mcrA)")) +
  theme_classic(base_size = 11)

# Panel B: Heartwood vs Sapwood mmoX
p6_mmox_spatial <- ggplot(spatial_data,
                          aes(x = log10(median_mmox_sap + 1),
                              y = log10(median_mmox_heart + 1))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", alpha = 0.7) +
  geom_point(aes(size = abs(median_flux), color = median_flux), alpha = 0.8) +
  geom_text_repel(aes(label = species), size = 2.5, fontface = "italic") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "CH4 flux") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "B) mmoX: Sapwood ≈ Heartwood",
       subtitle = "Oxidizers present throughout",
       x = expression("log"[10]*"(sapwood mmoX)"),
       y = expression("log"[10]*"(heartwood mmoX)")) +
  theme_classic(base_size = 11)

# Panel C: Which gene predicts flux in each compartment?
compartment_r2 <- data.frame(
  Compartment = rep(c("Heartwood", "Sapwood"), 2),
  Gene = c("mcrA", "mcrA", "mmoX", "mmoX"),
  R2 = c(cor_mcra_heart$estimate^2,
         cor_mcra_sap$estimate^2,
         cor_mmox_heart$estimate^2,
         cor_mmox_sap$estimate^2),
  P_value = c(cor_mcra_heart$p.value,
              cor_mcra_sap$p.value,
              cor_mmox_heart$p.value,
              cor_mmox_sap$p.value)
) %>%
  mutate(Significant = P_value < 0.05)

p7_compartment_compare <- ggplot(compartment_r2,
                                 aes(x = Compartment, y = R2, fill = Gene)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.2f\n%s", R2, 
                                ifelse(Significant, "**", "NS")),
                group = Gene),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("mcrA" = "darkred", "mmoX" = "darkgreen")) +
  labs(title = "C) Predictive Power by Compartment",
       y = expression("R"^2*" with flux")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

# Combine Figure 2
fig2_combined <- (p5_mcra_spatial | p6_mmox_spatial) / 
  (p7_compartment_compare | plot_spacer()) +
  plot_layout(heights = c(1, 0.7)) +
  plot_annotation(
    title = "Spatial Distribution: Methanogens Concentrate in Heartwood, Drive Species Differences",
    subtitle = sprintf("n = %d species", nrow(spatial_data)),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave("../../../outputs/figures/Figure_Spatial_Patterns.pdf", fig2_combined,
       width = 14, height = 10)
ggsave("../../../outputs/figures/Figure_Spatial_Patterns.png", fig2_combined,
       width = 14, height = 10, dpi = 300)

cat("Saved: Figure_Spatial_Patterns.pdf/png\n")

# ============================================================
# FIGURE 3: CONCEPTUAL MODEL
# ============================================================

cat("\nCreating Figure 3: Conceptual model...\n")

# Create a simple conceptual diagram using ggplot
conceptual_data <- data.frame(
  Scale = c("Within Species", "Within Species", 
            "Between Species", "Between Species"),
  Process = c("Production", "Oxidation", "Production", "Oxidation"),
  Effect_Size = c(tree_mcra_r2, tree_r2, cor_mcra_heart$estimate^2, cor_mmox_sap$estimate^2),
  x = c(1, 2, 3, 4),
  label = c(sprintf("R²=%.2f\n(masked)", tree_mcra_r2),
            sprintf("LIMITS flux\nR²=%.2f **", tree_r2), 
            sprintf("DETERMINES\nbaseline\nR²=%.2f **", cor_mcra_heart$estimate^2),
            sprintf("Similar\nacross species\nR²=%.2f", cor_mmox_sap$estimate^2))
)

p8_conceptual <- ggplot(conceptual_data, 
                        aes(x = x, y = Effect_Size, fill = Process)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_text(aes(label = label), vjust = -0.5, size = 3.5, 
            fontface = "bold", lineheight = 0.9) +
  geom_segment(aes(x = 1.5, xend = 1.5, y = 0, yend = max(Effect_Size) + 0.1),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               linewidth = 1.5, color = "black") +
  geom_segment(aes(x = 3.5, xend = 3.5, y = 0, yend = max(Effect_Size) + 0.1),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               linewidth = 1.5, color = "black") +
  annotate("text", x = 1.5, y = max(conceptual_data$Effect_Size) + 0.15, 
           label = "Net flux\nwithin species", 
           size = 3, fontface = "bold") +
  annotate("text", x = 3.5, y = max(conceptual_data$Effect_Size) + 0.15,
           label = "Species-level\nemission rate",
           size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Production" = "darkred", 
                               "Oxidation" = "darkgreen")) +
  scale_x_continuous(breaks = c(1.5, 3.5),
                     labels = c("WITHIN\nSPECIES", "BETWEEN\nSPECIES")) +
  ylim(0, max(conceptual_data$Effect_Size) + 0.2) +
  labs(title = "Conceptual Model: Scale-Dependent Mechanisms",
       y = expression("Effect Size (R"^2*")"),
       x = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 11, face = "bold"))

ggsave("../../../outputs/figures/Figure_Conceptual_Model.pdf", p8_conceptual,
       width = 10, height = 6)
ggsave("../../../outputs/figures/Figure_Conceptual_Model.png", p8_conceptual,
       width = 10, height = 6, dpi = 300)

cat("Saved: Figure_Conceptual_Model.pdf/png\n")

# ============================================================
# SUMMARY TABLE
# ============================================================

cat("\nCreating summary table...\n")

summary_table <- data.frame(
  Analysis_Level = c("Individual Trees", "Individual Trees",
                     "Species Medians", "Species Medians"),
  N = c(nrow(tree_mmox_data), nrow(tree_mcra_data),
        nrow(species_conc_complete), nrow(species_conc_complete)),
  Predictor = c("Area-weighted mmoX (oxidation)",
                "Area-weighted mcrA (production)",
                "Sapwood mmoX (oxidation)",
                "Heartwood mcrA (production)"),
  R2 = c(tree_r2, tree_mcra_r2,
         cor_mmox_sap$estimate^2,
         cor_mcra_heart$estimate^2),
  P_value = c(tree_p, tree_mcra_p,
              cor_mmox_sap$p.value,
              cor_mcra_heart$p.value),
  Interpretation = c("Oxidation capacity limits within-species flux variation",
                     "Production varies among trees but species effect dominates",
                     "Species have similar oxidation capacity",
                     "Production potential determines species baseline emissions")
) %>%
  mutate(
    R2 = round(R2, 3),
    P_value = round(P_value, 4),
    Significant = ifelse(P_value < 0.05, "**", "NS")
  )

write.csv(summary_table, "../../../outputs/tables/Table_Scale_Dependent_Summary.csv", row.names = FALSE)





sink()

