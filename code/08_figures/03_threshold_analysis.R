# ==============================================================================
# Threshold/Logistic Analysis (Figure S9)
# ==============================================================================
# Purpose: Threshold/logistic analysis for gene abundance effects on flux.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
# ==============================================================================

ymf2021 <- read.csv('data/processed/integrated/merged_tree_dataset_final.csv')


library(tidyverse)
library(ggplot2)
library(patchwork)
library(broom)

cat("============================================================\n")
cat("mcrA ANALYSIS: SPECIES AND DBH EFFECTS\n")
cat("FILTERED FOR SPECIES WITH n>5 OBSERVATIONS\n")
cat("============================================================\n\n")

# ============================================================
# STEP 1: Process ddPCR data to get area-weighted mcrA
# ============================================================

cat("STEP 1: Processing mcrA data from ymf2021\n")
cat("------------------------------------------\n")

# Parse ddPCR data
prepare_long_mcra <- function(df) {
  df %>%
    select(tree_id, starts_with("ddpcr_")) %>%
    select(where(~ !all(is.na(.)))) %>%
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

# Area-weighted calculation
area_weighted_mcra <- function(mcra_inner, mcra_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1*r2 + r2^2
  mcra_weighted <- mcra_outer + (mcra_inner - mcra_outer) * (S / (3 * R^2))
  
  return(mcra_weighted)
}

# Process the data
long_data <- prepare_long_mcra(ymf2021)

# Create initial analysis data
initial_data <- long_data %>%
  left_join(
    ymf2021 %>% select(tree_id, species_id, dbh),
    by = "tree_id"
  ) %>%
  filter(is.finite(dbh), !is.na(species_id)) %>%
  group_by(tree_id, species_id, dbh) %>%
  summarise(
    mcra_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
    mcra_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(mcra_inner), is.finite(mcra_outer))

# FILTER FOR SPECIES WITH n>5
species_counts <- initial_data %>%
  group_by(species_id) %>%
  summarise(n = n(), .groups = 'drop') %>%
  filter(n > 5)

cat(sprintf("Species with >5 observations: %d\n", nrow(species_counts)))
cat("Species included:\n")
print(species_counts)
cat("\n")

# Apply filter to analysis data
analysis_data <- initial_data %>%
  filter(species_id %in% species_counts$species_id) %>%
  mutate(
    mcra_weighted = mapply(area_weighted_mcra, mcra_inner, mcra_outer, dbh),
    log_mcra = log10(mcra_weighted + 1),
    elevated_1000 = mcra_weighted >= 1000,
    elevated_10000 = mcra_weighted >= 10000,
    dbh_scaled = scale(dbh)[,1]
  )

cat(sprintf("Total trees with complete data (n>5 filter): %d\n", nrow(analysis_data)))
cat(sprintf("Number of species: %d\n\n", n_distinct(analysis_data$species_id)))

# ============================================================
# STEP 2: Descriptive statistics
# ============================================================

cat("DESCRIPTIVE STATISTICS (n>5 filter)\n")
cat("-----------------------------------\n")

# Overall rates
cat("Overall elevation rates:\n")
cat(sprintf("  ≥1000 copies/g: %d trees (%.1f%%)\n", 
            sum(analysis_data$elevated_1000),
            100*mean(analysis_data$elevated_1000)))
cat(sprintf("  ≥10000 copies/g: %d trees (%.1f%%)\n\n", 
            sum(analysis_data$elevated_10000),
            100*mean(analysis_data$elevated_10000)))

# By species (already filtered for n>5)
species_summary <- analysis_data %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    median_mcra = round(median(mcra_weighted), 0),
    mean_dbh = round(mean(dbh), 1),
    pct_1000 = round(100*mean(elevated_1000), 1),
    pct_10000 = round(100*mean(elevated_10000), 1),
    .groups = 'drop'
  ) %>%
  arrange(desc(pct_1000))

cat("Summary by species (all have n>5):\n")
print(species_summary)
cat("\n")

# ============================================================
# STEP 3: Logistic regression - 1000 threshold
# ============================================================

cat("============================================================\n")
cat("LOGISTIC REGRESSION: 1000 copies/g threshold (n>5 filter)\n")
cat("============================================================\n\n")

# Check for complete separation
sep_check_1000 <- analysis_data %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    n_elevated = sum(elevated_1000),
    .groups = 'drop'
  ) %>%
  mutate(issue = ifelse(n_elevated == 0, "All non-elevated",
                        ifelse(n_elevated == n, "All elevated", "OK")))

problem_species_1000 <- filter(sep_check_1000, issue != "OK")
if(nrow(problem_species_1000) > 0) {
  cat("Species with complete separation:\n")
  print(problem_species_1000)
  cat("\n")
}

# Fit model
model_1000 <- glm(elevated_1000 ~ species_id + dbh_scaled,
                  data = analysis_data,
                  family = binomial(link = "logit"))

# Check convergence
if(!model_1000$converged) {
  cat("Warning: Model did not converge properly.\n\n")
}

# Model summary
summary_1000 <- summary(model_1000)
cat("Model summary:\n")
print(summary_1000)

# Extract key results
cat("\n\nKEY RESULTS:\n")
cat("------------\n")

# DBH effect
dbh_coef <- coef(summary_1000)["dbh_scaled",]
cat(sprintf("DBH effect: OR = %.3f, p = %.4f\n", 
            exp(dbh_coef["Estimate"]), dbh_coef["Pr(>|z|)"]))

# Species with highest odds (relative to reference)
species_coefs <- coef(summary_1000)[grep("species_id", rownames(coef(summary_1000))),]
if(nrow(as.matrix(species_coefs)) > 0) {
  species_or <- data.frame(
    species = gsub("species_id", "", rownames(as.matrix(species_coefs))),
    OR = exp(species_coefs[,"Estimate"]),
    p_value = species_coefs[,"Pr(>|z|)"]
  ) %>%
    arrange(desc(OR))
  
  cat("\nTop 5 species effects (OR relative to reference):\n")
  print(head(species_or, 5))
}

# Model fit
cat(sprintf("\nModel fit:\n"))
cat(sprintf("  AIC: %.1f\n", AIC(model_1000)))
cat(sprintf("  McFadden R²: %.3f\n", 1 - model_1000$deviance/model_1000$null.deviance))

# ============================================================
# STEP 4: Logistic regression - 10000 threshold
# ============================================================

cat("\n============================================================\n")
cat("LOGISTIC REGRESSION: 10000 copies/g threshold (n>5 filter)\n")
cat("============================================================\n\n")

if(sum(analysis_data$elevated_10000) >= 5) {
  
  # Check separation
  sep_check_10000 <- analysis_data %>%
    group_by(species_id) %>%
    summarise(
      n = n(),
      n_elevated = sum(elevated_10000),
      .groups = 'drop'
    ) %>%
    mutate(issue = ifelse(n_elevated == 0, "All non-elevated",
                          ifelse(n_elevated == n, "All elevated", "OK")))
  
  problem_species_10000 <- filter(sep_check_10000, issue != "OK")
  if(nrow(problem_species_10000) > 0) {
    cat("Species with complete separation:\n")
    print(problem_species_10000)
    cat("\n")
  }
  
  # Fit model
  model_10000 <- glm(elevated_10000 ~ species_id + dbh_scaled,
                     data = analysis_data,
                     family = binomial(link = "logit"))
  
  if(!model_10000$converged) {
    cat("Warning: Model did not converge properly.\n\n")
  }
  
  # Summary
  summary_10000 <- summary(model_10000)
  
  # Key results
  cat("KEY RESULTS:\n")
  cat("------------\n")
  
  # DBH effect
  dbh_coef_10000 <- coef(summary_10000)["dbh_scaled",]
  cat(sprintf("DBH effect: OR = %.3f, p = %.4f\n", 
              exp(dbh_coef_10000["Estimate"]), dbh_coef_10000["Pr(>|z|)"]))
  
  cat(sprintf("\nModel fit:\n"))
  cat(sprintf("  AIC: %.1f\n", AIC(model_10000)))
  cat(sprintf("  McFadden R²: %.3f\n", 1 - model_10000$deviance/model_10000$null.deviance))
  
} else {
  cat("Too few elevated cases for modeling.\n")
}

# ============================================================
# STEP 5: Visualization
# ============================================================

cat("\n============================================================\n")
cat("CREATING VISUALIZATIONS (n>5 filter)\n")
cat("============================================================\n\n")

# Plot 1: Species comparison (all species already have n>5)
species_plot_data <- analysis_data %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    pct_1000 = 100*mean(elevated_1000),
    pct_10000 = 100*mean(elevated_10000),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(pct_1000, pct_10000),
               names_to = "threshold",
               values_to = "percent") %>%
  mutate(threshold = factor(ifelse(threshold == "pct_1000", "≥1000", "≥10000"),
                            levels = c("≥1000", "≥10000")))

p1 <- ggplot(species_plot_data, 
             aes(x = reorder(species_id, percent), y = percent, fill = threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("≥1000" = "#2166AC", "≥10000" = "#B2182B")) +
  labs(
    title = "Elevated mcrA by Species",
    subtitle = "Species with n>5",
    x = "Species",
    y = "Percentage of trees (%)",
    fill = "Threshold\n(copies/g)"
  ) +
  theme_classic()

print(p1)

# Plot 2: DBH relationship
p2 <- ggplot(analysis_data, aes(x = dbh, y = log_mcra)) +
  geom_point(aes(color = species_id), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  geom_hline(yintercept = log10(1000), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = log10(10000), linetype = "dashed", color = "darkred", alpha = 0.7) +
  labs(
    title = "mcrA Concentration vs Tree Size",
    subtitle = "Red lines show 1000 and 10000 thresholds (species with n>5)",
    x = "DBH (cm)",
    y = "log10(mcrA copies/g)",
    color = "Species"
  ) +
  theme_classic() +
  theme(legend.position = "right")

print(p2)

cat("\nAnalysis complete.\n")
cat("============================================================\n")









# ============================================================
# Visualization of mcrA probability by species
# Multiple approaches for categorical predictor
# FILTERED FOR SPECIES WITH n>5 OBSERVATIONS
# ============================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

cat("============================================================\n")
cat("SPECIES-SPECIFIC PROBABILITY PLOTS (n>5 filter)\n")
cat("============================================================\n\n")

# ============================================================
# OPTION 1: Forest plot showing predicted probabilities
# ============================================================

cat("Creating forest plot of predicted probabilities...\n")

# Get predicted probabilities from the model for each species
# Using mean DBH for predictions
mean_dbh_scaled <- 0  # standardized mean

# Create prediction dataset (all species already have n>5)
species_list <- unique(analysis_data$species_id)
pred_data <- data.frame(
  species_id = species_list,
  dbh_scaled = mean_dbh_scaled
)

# Get predictions for 1000 threshold
pred_data$prob_1000 <- predict(model_1000, newdata = pred_data, type = "response")
pred_data$se_1000 <- predict(model_1000, newdata = pred_data, type = "response", se.fit = TRUE)$se.fit

# Calculate confidence intervals
pred_data$lower_1000 <- pred_data$prob_1000 - 1.96 * pred_data$se_1000
pred_data$upper_1000 <- pred_data$prob_1000 + 1.96 * pred_data$se_1000
pred_data$lower_1000[pred_data$lower_1000 < 0] <- 0
pred_data$upper_1000[pred_data$upper_1000 > 1] <- 1

# Add observed proportions (all species already have n>5)
obs_props <- analysis_data %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    obs_prop_1000 = mean(elevated_1000),
    .groups = 'drop'
  )

pred_data <- pred_data %>%
  left_join(obs_props, by = "species_id")

# Forest plot
p1 <- ggplot(pred_data, aes(x = prob_1000, y = reorder(species_id, prob_1000))) +
  # Predicted probabilities with CI
  geom_errorbarh(aes(xmin = lower_1000, xmax = upper_1000), 
                 height = 0.3, color = "gray50", alpha = 0.7) +
  geom_point(aes(size = n), color = "#1976D2", alpha = 0.7) +
  # Observed proportions
  geom_point(aes(x = obs_prop_1000), shape = 4, size = 3, color = "red") +
  # Reference lines
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = mean(analysis_data$elevated_1000), 
             linetype = "dotted", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), labels = percent) +
  scale_size_continuous(range = c(2, 6), name = "Sample size") +
  labs(
    title = "Predicted Probability of Elevated mcrA (≥1000 copies/g)",
    subtitle = "Blue circles = model predictions, Red X = observed proportions (species n>5)",
    x = "Probability",
    y = "Species"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(face = "italic"),
    panel.grid.major.x = element_line(color = "gray90")
  )

print(p1)

# ============================================================
# OPTION 2: Grouped bar plot with multiple thresholds
# ============================================================

cat("\nCreating grouped bar plot for multiple thresholds...\n")

# Calculate proportions for multiple thresholds (all species already have n>5)
threshold_props <- analysis_data %>%
  mutate(
    elevated_500 = mcra_weighted >= 500,
    elevated_2000 = mcra_weighted >= 2000,
    elevated_5000 = mcra_weighted >= 5000
  ) %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    `500` = mean(elevated_500),
    `1000` = mean(elevated_1000),
    `2000` = mean(elevated_2000),
    `5000` = mean(elevated_5000),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(`500`, `1000`, `2000`, `5000`),
               names_to = "threshold",
               values_to = "proportion") %>%
  mutate(threshold = factor(threshold, levels = c("500", "1000", "2000", "5000")))

# Create color palette similar to the heart rot paper
threshold_colors <- c("500" = "#4CAF50", "1000" = "#2196F3", 
                      "2000" = "#FF9800", "5000" = "#F44336")

p2 <- ggplot(threshold_props, 
             aes(x = reorder(species_id, proportion), y = proportion, fill = threshold)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), labels = percent) +
  scale_fill_manual(values = threshold_colors, name = "Threshold\n(copies/g)") +
  labs(
    title = "Probability of Elevated mcrA by Species",
    subtitle = "Multiple threshold comparison (species with n>5)",
    x = "Species",
    y = "Proportion of trees"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(face = "italic"),
    panel.grid.major.x = element_line(color = "gray90"),
    legend.position = "right"
  )

print(p2)

# ============================================================
# OPTION 3: Dot plot with DBH interaction
# ============================================================

cat("\nCreating species × DBH interaction plot...\n")

# Generate predictions across DBH range for each species
dbh_range <- seq(min(analysis_data$dbh), max(analysis_data$dbh), length.out = 50)
dbh_scaled_range <- (dbh_range - mean(analysis_data$dbh)) / sd(analysis_data$dbh)

# All species in analysis_data already have n>5, so we can use them all
# Or select top species for clarity
top_species <- analysis_data %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    prop = mean(elevated_1000),
    .groups = 'drop'
  ) %>%
  arrange(desc(prop)) %>%
  head(6) %>%
  pull(species_id)

# Create prediction grid
pred_grid <- expand.grid(
  species_id = top_species,
  dbh = dbh_range,
  dbh_scaled = dbh_scaled_range
)

# Get predictions
pred_grid$probability <- predict(model_1000, newdata = pred_grid, type = "response")

# Add actual data points
actual_points <- analysis_data %>%
  filter(species_id %in% top_species) %>%
  select(species_id, dbh, elevated_1000)

p3 <- ggplot() +
  # Prediction lines
  geom_line(data = pred_grid,
            aes(x = dbh, y = probability, color = species_id),
            size = 1.2, alpha = 0.8) +
  # Actual data points
  geom_point(data = actual_points,
             aes(x = dbh, y = as.numeric(elevated_1000), color = species_id),
             alpha = 0.4, size = 2) +
  facet_wrap(~ species_id, ncol = 3, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), labels = percent) +
  scale_color_brewer(palette = "Set1", guide = "none") +
  labs(
    title = "Probability of Elevated mcrA vs DBH by Species",
    subtitle = "Lines = model predictions, Points = actual observations (species n>5)",
    x = "DBH (cm)",
    y = "Probability of mcrA ≥1000 copies/g"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "italic"),
    strip.background = element_blank()
  )

print(p3)

# ============================================================
# OPTION 4: Heatmap showing all species × threshold combinations
# ============================================================

cat("\nCreating heatmap of species × threshold combinations...\n")

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

# Calculate proportions for all combinations (all species already have n>5)
heatmap_data <- analysis_data %>%
  mutate(
    `1000` = mcra_weighted >= 1000,
    `5000` = mcra_weighted >= 5000,
    `10000` = mcra_weighted >= 10000
  ) %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    across(c(`1000`, `5000`, `10000`), mean),
    .groups = 'drop'
  )

# Get species order based on 10000 threshold, then 5000, then 1000
species_order <- heatmap_data %>%
  arrange(`10000`, `5000`, `1000`) %>%
  pull(species_id)

# Reshape data and apply ordering
heatmap_data <- heatmap_data %>%
  select(-n) %>%
  pivot_longer(cols = -species_id,
               names_to = "threshold",
               values_to = "proportion") %>%
  mutate(
    threshold = factor(as.numeric(threshold)),
    species_name = species_mapping[species_id],
    species_name = factor(species_name, 
                          levels = species_mapping[species_order])
  )

p4 <- ggplot(heatmap_data, 
             aes(x = threshold, y = species_name, fill = proportion)) +
  geom_tile(color = "white", size = 0.1, alpha=0.75) +
  geom_text(aes(label = sprintf("%.0f%%", proportion * 100),
                color = proportion < 0.5), 
            size = 3) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
#  scale_fill_gradient(low = "white", high = "#560591", limits = c(0, 1), name = "Proportion\nelevated") +
  scale_fill_viridis_c(option="plasma",name = "Proportion\nelevated") +
  
   labs(
    #title = "Proportion of Trees with Elevated mcrA",
    #subtitle = "By species and threshold (species n>5)",
    x = "Threshold (copies/g)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "italic", size=8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p4)

# ============================================================
# Combine favorite plots
# ============================================================

cat("\nCombining plots...\n")

# Combine forest plot and grouped bar plot
combined <- p1 / p2
print(combined)

cat("\nAll plots created successfully!\n")
cat("Note: All analyses limited to species with >5 observations\n")
cat("============================================================\n")