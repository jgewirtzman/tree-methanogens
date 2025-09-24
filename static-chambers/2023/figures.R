library(ggplot2)
library(ggridges)
library(dplyr)
library(viridis)

# Assuming your dataframe is called ymf2023
# Using CH4_best.flux as the flux variable and vwc_mean for coloring

# Prepare the data - filter for species with n > 3
plot_data <- ymf2023 %>%
  filter(!is.na(CH4_best.flux) & !is.na(Species) & !is.na(vwc_mean)) %>%
  group_by(Species) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  # Calculate species order by mean flux for better visualization
  mutate(Species = reorder(Species, CH4_best.flux, mean, na.rm = TRUE))

# Calculate mean and SE for each species
species_stats <- plot_data %>%
  group_by(Species) %>%
  summarise(
    mean_flux = mean(CH4_best.flux, na.rm = TRUE),
    se_flux = sd(CH4_best.flux, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

# Create the ridgeline plot
p <- ggplot(plot_data, aes(x = CH4_best.flux, y = Species)) +
  # Add the density ridges with a single fill color
  geom_density_ridges(
    fill = "lightgray",
    alpha = 0.6,
    scale = 1.5,
    rel_min_height = 0.01,
    color = "gray30",
    size = 0.3
  ) +
  # Add jittered points colored by VWC
  geom_point(
    aes(color = vwc_mean),
    position = position_jitter(width = 0, height = 0.15),
    size = 1.5,
    alpha = 0.7
  ) +
  # Add mean point with error bars
  geom_pointrange(
    data = species_stats,
    aes(x = mean_flux, 
        xmin = mean_flux - se_flux, 
        xmax = mean_flux + se_flux,
        y = Species),
    color = "black",
    size = 0.5,
    linewidth = 0.7,
    fatten = 2,
    inherit.aes = FALSE
  ) +
  # Color scale for VWC
  scale_color_viridis(
    name = "Mean VWC (%)",
    #option = "plasma",
    limits = c(0, max(plot_data$vwc_mean, na.rm = TRUE))
  ) +
  # Add vertical line at zero to show emission vs uptake
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  # Use pseudo-log scale (symmetric log scale) for x-axis
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.01),
    breaks = c(-1, -0.1, -0.01, 0, 0.01, 0.1, 1),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  # Labels and theme
  labs(
    #title = "Distribution of CH₄ Flux by Tree Species",
    #subtitle = "Points colored by mean VWC; Black points show mean ± SE; Species with n > 3 only",
    x = expression(CH[4]~Flux~(nmol~m^{-2}~s^{-1})),
    y = "Species"
  ) +
  theme_ridges(grid = TRUE) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    panel.grid.minor.x = element_line(color = "gray95", size = 0.2)
  )

# Print the plot
print(p)

# Print summary statistics by species
cat("\nSummary of CH4 Flux by Species (n > 3):\n")
cat("========================================\n")
summary_stats <- plot_data %>%
  group_by(Species) %>%
  summarise(
    n = n(),
    mean_flux = mean(CH4_best.flux, na.rm = TRUE),
    se_flux = sd(CH4_best.flux, na.rm = TRUE) / sqrt(n()),
    median_flux = median(CH4_best.flux, na.rm = TRUE),
    sd_flux = sd(CH4_best.flux, na.rm = TRUE),
    min_flux = min(CH4_best.flux, na.rm = TRUE),
    max_flux = max(CH4_best.flux, na.rm = TRUE),
    mean_vwc = mean(vwc_mean, na.rm = TRUE),
    pct_uptake = sum(CH4_best.flux < 0, na.rm = TRUE) / n() * 100,
    pct_emission = sum(CH4_best.flux > 0, na.rm = TRUE) / n() * 100,
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_flux))

print(summary_stats)

# Print overall statistics
cat("\n\nOverall CH4 Flux Statistics (filtered data):\n")
cat("=============================================\n")
cat(sprintf("Number of species included: %d\n", n_distinct(plot_data$Species)))
cat(sprintf("Total observations: %d\n", nrow(plot_data)))
cat(sprintf("CH4 uptake (negative flux): %.1f%%\n", 
            sum(plot_data$CH4_best.flux < 0, na.rm = TRUE) / nrow(plot_data) * 100))
cat(sprintf("CH4 emission (positive flux): %.1f%%\n", 
            sum(plot_data$CH4_best.flux > 0, na.rm = TRUE) / nrow(plot_data) * 100))
cat(sprintf("Mean flux: %.4f nmol/m²/s\n", mean(plot_data$CH4_best.flux, na.rm = TRUE)))
cat(sprintf("Median flux: %.4f nmol/m²/s\n", median(plot_data$CH4_best.flux, na.rm = TRUE)))

# Optional: Save the plot
# ggsave("ch4_flux_ridgeplot.png", plot = p, width = 10, height = 8, dpi = 300)

# Create scatterplot of flux vs VWC faceted by species
# ========================================================

# Calculate correlation and sample size for each species
species_correlations <- plot_data %>%
  group_by(Species) %>%
  summarise(
    cor_value = cor(CH4_best.flux, vwc_mean, use = "complete.obs"),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    cor_label = paste0("r = ", sprintf("%.2f", cor_value), "\nn = ", n),
    # Position for text annotation
    x_pos = Inf,
    y_pos = Inf
  )

# Create the scatterplot
p2 <- ggplot(plot_data, aes(x = vwc_mean, y = CH4_best.flux)) +
  # Add points
  geom_point(
    aes(color = vwc_mean),
    alpha = 0.6,
    size = 2
  ) +
  # Add smoothing line with confidence interval
  geom_smooth(
    method = "loess",
    se = TRUE,
    color = "black",
    fill = "gray80",
    alpha = 0.3,
    size = 0.8
  ) +
  # Add horizontal line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  # Add correlation text
  geom_text(
    data = species_correlations,
    aes(x = x_pos, y = y_pos, label = cor_label),
    hjust = 1.1,
    vjust = 1.1,
    size = 3,
    color = "gray30"
  ) +
  # Facet by species
  facet_wrap(~ Species, scales = "free", ncol = 3) +
  # Color scale for VWC (same as ridgeplot)
  scale_color_viridis(
    name = "Mean VWC (%)",
    option = "plasma",
    limits = c(0, max(plot_data$vwc_mean, na.rm = TRUE))
  ) +
  # Use pseudo-log scale for y-axis
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.01),
    breaks = c(-1, -0.1, -0.01, 0, 0.01, 0.1, 1),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  # Labels and theme
  labs(
    title = "CH₄ Flux vs. Volumetric Water Content by Tree Species",
    subtitle = "LOESS smoothing with 95% CI; Species with n > 3 only",
    x = "Mean VWC (%)",
    y = expression(CH[4]~Flux~(nmol~m^{-2}~s^{-1}))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Print the scatterplot
print(p2)

# Create a combined plot with both visualizations
# ================================================
library(patchwork)

# Combine plots vertically
combined_plot <- p / p2 + 
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "CH₄ Flux Analysis by Tree Species and Soil Moisture",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Print combined plot
print(combined_plot)

# Optional: Save the plots
# ggsave("ch4_flux_ridgeplot.png", plot = p, width = 10, height = 8, dpi = 300)
# ggsave("ch4_flux_vs_vwc_scatter.png", plot = p2, width = 12, height = 10, dpi = 300)
# ggsave("ch4_flux_combined.png", plot = combined_plot, width = 12, height = 16, dpi = 300)


library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(broom)
library(car)
library(vegan)

# Prepare data for analysis
analysis_data <- ymf2023 %>%
  select(
    # Species
    Species,
    # Environmental drivers
    DBH..cm.,
    air_temp_C,
    stem_temp_C,
    Soil.Temp...C.,
    vwc_mean,
    # Response variable
    CH4_best.flux
  ) %>%
  # Remove rows with any NA values
  drop_na() %>%
  # Filter for species with n > 3
  group_by(Species) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  # Rename for cleaner analysis
  rename(
    DBH = DBH..cm.,
    Air_temp = air_temp_C,
    Stem_temp = stem_temp_C,
    Soil_temp = Soil.Temp...C.,
    VWC = vwc_mean,
    CH4_flux = CH4_best.flux
  )

# Scale environmental variables for comparable coefficients
env_vars <- c("DBH", "Air_temp", "Stem_temp", "Soil_temp", "VWC")
analysis_data[env_vars] <- scale(analysis_data[env_vars])

cat("========================================\n")
cat("VARIANCE PARTITIONING ANALYSIS\n")
cat("CH4 Flux: Environment vs Species Effects\n")
cat("========================================\n\n")

# 1. NULL MODEL (baseline)
model_null <- lm(CH4_flux ~ 1, data = analysis_data)

# 2. ENVIRONMENT ONLY MODEL
model_env <- lm(CH4_flux ~ DBH + Air_temp + Stem_temp + Soil_temp + VWC, 
                data = analysis_data)

# 3. SPECIES ONLY MODEL
model_species <- lm(CH4_flux ~ Species, data = analysis_data)

# 4. FULL MODEL (Environment + Species)
model_full <- lm(CH4_flux ~ DBH + Air_temp + Stem_temp + Soil_temp + VWC + Species, 
                 data = analysis_data)

# 5. INTERACTION MODEL (Environment * Species)
model_interaction <- lm(CH4_flux ~ (DBH + Air_temp + Stem_temp + Soil_temp + VWC) * Species, 
                        data = analysis_data)

# Calculate R-squared values
r2_env <- summary(model_env)$r.squared
r2_species <- summary(model_species)$r.squared
r2_full <- summary(model_full)$r.squared
r2_interaction <- summary(model_interaction)$r.squared

# Variance partitioning
cat("VARIANCE EXPLAINED (R²):\n")
cat("========================\n")
cat(sprintf("Environment only:        %.1f%%\n", r2_env * 100))
cat(sprintf("Species only:            %.1f%%\n", r2_species * 100))
cat(sprintf("Environment + Species:   %.1f%%\n", r2_full * 100))
cat(sprintf("With interactions:       %.1f%%\n", r2_interaction * 100))
cat("\n")

# Calculate unique and shared contributions
unique_env <- r2_full - r2_species
unique_species <- r2_full - r2_env
shared <- r2_env + r2_species - r2_full

cat("VARIANCE DECOMPOSITION:\n")
cat("=======================\n")
cat(sprintf("Unique to Environment:   %.1f%%\n", max(0, unique_env) * 100))
cat(sprintf("Unique to Species:       %.1f%%\n", max(0, unique_species) * 100))
cat(sprintf("Shared (confounded):     %.1f%%\n", max(0, shared) * 100))
cat(sprintf("Unexplained:             %.1f%%\n", (1 - r2_full) * 100))
cat("\n")

# ANOVA to test significance
cat("STATISTICAL SIGNIFICANCE:\n")
cat("=========================\n")
cat("\nModel Comparison (ANOVA):\n")
print(anova(model_null, model_env, model_species, model_full))

# Type III ANOVA for the full model (tests each factor controlling for others)
cat("\n\nType III ANOVA (Full Model):\n")
cat("------------------------------\n")
Anova_type3 <- Anova(model_full, type = "III")
print(Anova_type3)

# Calculate effect sizes (partial eta-squared)
SS_total <- sum(Anova_type3$`Sum Sq`)
Anova_type3$Partial_Eta2 <- Anova_type3$`Sum Sq` / 
  (Anova_type3$`Sum Sq` + Anova_type3$`Sum Sq`[length(Anova_type3$`Sum Sq`)])

cat("\n\nEFFECT SIZES (Partial η²):\n")
cat("===========================\n")
effect_sizes <- data.frame(
  Variable = rownames(Anova_type3)[-length(rownames(Anova_type3))],
  Partial_Eta2 = round(Anova_type3$Partial_Eta2[-length(Anova_type3$Partial_Eta2)], 3),
  Percent = round(Anova_type3$Partial_Eta2[-length(Anova_type3$Partial_Eta2)] * 100, 1)
)
print(effect_sizes)

# Environmental variable importance
cat("\n\nENVIRONMENTAL VARIABLE COEFFICIENTS:\n")
cat("=====================================\n")
env_coef <- summary(model_env)$coefficients
print(round(env_coef, 4))

# Create visualization of variance partitioning
variance_data <- data.frame(
  Component = c("Environment\nUnique", "Species\nUnique", "Shared\n(Confounded)", "Unexplained"),
  Variance = c(max(0, unique_env), max(0, unique_species), max(0, shared), 1 - r2_full) * 100
)

p_variance <- ggplot(variance_data, aes(x = Component, y = Variance, fill = Component)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Variance)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Environment\nUnique" = "#2E7D32",
                               "Species\nUnique" = "#1976D2", 
                               "Shared\n(Confounded)" = "#7B1FA2",
                               "Unexplained" = "gray70")) +
  labs(title = "Variance Partitioning of CH₄ Flux",
       subtitle = "Relative importance of environmental factors vs. species identity",
       y = "Variance Explained (%)",
       x = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 11),
        legend.position = "none") +
  ylim(0, max(variance_data$Variance) * 1.2)

print(p_variance)

# PCA visualization with environment vs species coloring
# Perform PCA on environmental variables
env_matrix <- as.matrix(analysis_data[env_vars])
pca_result <- prcomp(env_matrix, center = FALSE, scale. = FALSE)

# Create PCA plot colored by species and sized by flux
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$Species <- analysis_data$Species
pca_scores$CH4_flux <- analysis_data$CH4_flux

# Calculate centroids for each species
species_centroids <- pca_scores %>%
  group_by(Species) %>%
  summarise(PC1_mean = mean(PC1),
            PC2_mean = mean(PC2),
            mean_flux = mean(CH4_flux),
            .groups = 'drop')

var_explained <- summary(pca_result)$importance[2, 1:2] * 100

p_pca <- ggplot() +
  # Points colored by species
  geom_point(data = pca_scores,
             aes(x = PC1, y = PC2, color = Species, size = abs(CH4_flux)),
             alpha = 0.6) +
  # Add species centroids
  geom_point(data = species_centroids,
             aes(x = PC1_mean, y = PC2_mean, fill = Species),
             color = "black", shape = 21, size = 5) +
  geom_text(data = species_centroids,
            aes(x = PC1_mean, y = PC2_mean, label = substr(Species, 1, 3)),
            size = 3, fontface = "bold") +
  scale_size_continuous(name = expression("|CH"[4]*" Flux|"),
                        range = c(1, 5)) +
  labs(title = "Environmental Space Occupied by Different Species",
       subtitle = "PCA of environmental conditions; Point size = |CH₄ flux|",
       x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")) +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "right") +
  coord_equal()

print(p_pca)

# Calculate within vs between species variance
cat("\n\nWITHIN vs BETWEEN SPECIES VARIATION:\n")
cat("=====================================\n")

# Calculate ICC (Intraclass Correlation Coefficient)
library(lme4)
library(performance)

# Fit random effects model
model_random <- lmer(CH4_flux ~ (1|Species), data = analysis_data)
icc_species <- icc(model_random)

cat(sprintf("Intraclass Correlation (Species): %.3f\n", icc_species$ICC_adjusted))
cat(sprintf("  → %.1f%% of variance is between species\n", icc_species$ICC_adjusted * 100))
cat(sprintf("  → %.1f%% of variance is within species\n", (1 - icc_species$ICC_adjusted) * 100))

# Summary recommendations
cat("\n\n========================================\n")
cat("KEY FINDINGS:\n")
cat("========================================\n")
if (r2_species > r2_env) {
  cat("→ Species identity explains more variance than environment\n")
  cat("  Suggests: Strong phylogenetic/functional group control\n")
} else {
  cat("→ Environmental factors explain more variance than species\n")
  cat("  Suggests: CH₄ flux primarily driven by local conditions\n")
}

if (shared > 0.1) {
  cat("→ Substantial shared variance between species and environment\n")
  cat("  Suggests: Species occupy different environmental niches\n")
}

if (r2_interaction - r2_full > 0.1) {
  cat("→ Strong species × environment interactions\n")
  cat("  Suggests: Species respond differently to environmental gradients\n")
}

# Save plots
# ggsave("variance_partitioning.png", plot = p_variance, width = 8, height = 6, dpi = 300)
# ggsave("pca_species_environment.png", plot = p_pca, width = 10, height = 8, dpi = 300)








library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(ggforce)
library(broom)

# Assuming analysis_data and models are already loaded from previous script

# ============================================
# 1. VENN DIAGRAM STYLE VARIANCE PARTITIONING
# ============================================

# Create a more intuitive variance partitioning visualization
library(VennDiagram)
library(grid)

# Calculate variance components
r2_env <- 0.014  # From your output
r2_species <- 0.084
r2_full <- 0.089
r2_interaction <- 0.315

unique_env <- r2_full - r2_species  # 0.005
unique_species <- r2_full - r2_env  # 0.075
shared <- r2_env + r2_species - r2_full  # 0.009

# Create a custom Euler diagram
variance_components <- data.frame(
  category = c("Environment\nalone", "Species\nalone", "Shared", 
               "Interaction\neffects", "Unexplained"),
  value = c(unique_env, unique_species, shared, 
            r2_interaction - r2_full, 1 - r2_interaction) * 100,
  x = c(-1, 1, 0, 0, 3),
  y = c(0, 0, 0, -1.5, 0),
  color = c("#2E7D32", "#1976D2", "#7B1FA2", "#FF6F00", "gray70")
)

p_euler <- ggplot(variance_components) +
  # Main circles for Environment and Species
  geom_circle(aes(x0 = -0.5, y0 = 0, r = sqrt(r2_env * 10)), 
              fill = "#2E7D32", alpha = 0.3, color = "#2E7D32", size = 2) +
  geom_circle(aes(x0 = 0.5, y0 = 0, r = sqrt(r2_species * 10)), 
              fill = "#1976D2", alpha = 0.3, color = "#1976D2", size = 2) +
  # Labels
  geom_text(data = variance_components[1:3,],
            aes(x = x, y = y, label = paste0(category, "\n", round(value, 1), "%")),
            size = 4, fontface = "bold") +
  # Interaction box
  annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -2, ymax = -1,
           fill = "#FF6F00", alpha = 0.3) +
  annotate("text", x = 0, y = -1.5, 
           label = paste0("Species × Environment\nInteractions\n", 
                          round((r2_interaction - r2_full) * 100, 1), "%"),
           size = 4, fontface = "bold") +
  # Unexplained
  annotate("text", x = 3, y = 0,
           label = paste0("Unexplained\n", 
                          round((1 - r2_interaction) * 100, 1), "%"),
           size = 4, color = "gray50") +
  coord_equal() +
  xlim(-3, 4) + ylim(-2.5, 2) +
  theme_void() +
  labs(title = "Variance Components of CH₄ Flux",
       subtitle = "Unique and shared contributions + interaction effects") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5))

print(p_euler)

# ============================================
# 2. COEFFICIENT PLOT - DIRECTION OF EFFECTS
# ============================================

# Extract coefficients from the full model
model_full <- lm(CH4_flux ~ DBH + Air_temp + Stem_temp + Soil_temp + VWC + Species, 
                 data = analysis_data)

# Get coefficients with confidence intervals
coef_data <- tidy(model_full, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    category = ifelse(grepl("Species", term), "Species", "Environment"),
    term_clean = gsub("Species", "", term),
    significant = p.value < 0.05
  )

# Separate environment and species effects
env_coef <- coef_data %>% filter(category == "Environment")
species_coef <- coef_data %>% filter(category == "Species")

# Create coefficient plot
p_coef <- ggplot(coef_data, aes(x = estimate, y = reorder(term_clean, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = category),
                 height = 0.2, size = 0.8, alpha = 0.7) +
  geom_point(aes(color = category, shape = significant), size = 3) +
  scale_color_manual(values = c("Environment" = "#2E7D32", "Species" = "#1976D2")) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("Not significant", "Significant (p<0.05)")) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Direction and Magnitude of Effects on CH₄ Flux",
       subtitle = "Coefficient estimates with 95% CI from full model",
       x = "Standardized Effect Size",
       y = "",
       color = "Factor Type",
       shape = "Significance") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom")

print(p_coef)

# ============================================
# 3. SPECIES-SPECIFIC ENVIRONMENTAL RESPONSES
# ============================================

# Prepare data for species-specific slopes
species_list <- unique(analysis_data$Species)
env_vars <- c("VWC", "Air_temp", "Soil_temp")

# Calculate species-specific slopes for key environmental variables
slopes_data <- data.frame()

for (sp in species_list) {
  sp_data <- analysis_data %>% filter(Species == sp)
  for (env in env_vars) {
    if (nrow(sp_data) > 3) {
      formula_str <- paste("CH4_flux ~", env)
      model <- lm(as.formula(formula_str), data = sp_data)
      slope <- coef(model)[2]
      se <- summary(model)$coefficients[2, 2]
      
      slopes_data <- rbind(slopes_data, 
                           data.frame(Species = sp, 
                                      Environment = env,
                                      Slope = slope,
                                      SE = se,
                                      n = nrow(sp_data)))
    }
  }
}

# Create heatmap of species-specific responses
p_heatmap <- ggplot(slopes_data, aes(x = Environment, y = Species, fill = Slope)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Slope)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,
                       name = "Slope\n(effect on CH₄ flux)") +
  labs(title = "Species-Specific Environmental Responses",
       subtitle = "How each species responds to environmental gradients",
       x = "Environmental Variable",
       y = "Species") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_heatmap)

# ============================================
# 4. REACTION NORM PLOT
# ============================================

# Focus on VWC as the most important environmental variable
p_reaction <- ggplot(analysis_data, aes(x = VWC, y = CH4_flux, color = Species)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  facet_wrap(~ Species, scales = "free_y", ncol = 3) +
  labs(title = "Species-Specific Reaction Norms",
       subtitle = "CH₄ flux response to soil moisture (VWC) by species",
       x = "Volumetric Water Content (standardized)",
       y = expression(CH[4]~Flux~(nmol~m^{-2}~s^{-1}))) +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(face = "bold"))

print(p_reaction)

# ============================================
# 5. INTERACTION VISUALIZATION
# ============================================

# Create predicted values for interaction plot
vwc_range <- seq(min(analysis_data$VWC), max(analysis_data$VWC), length.out = 50)
predictions <- data.frame()

for (sp in species_list) {
  pred_data <- data.frame(
    VWC = vwc_range,
    DBH = 0,  # Set other variables to mean (0 after scaling)
    Air_temp = 0,
    Stem_temp = 0,
    Soil_temp = 0,
    Species = sp
  )
  
  preds <- predict(model_full, newdata = pred_data, se.fit = TRUE)
  pred_data$CH4_flux_pred <- preds$fit
  pred_data$SE <- preds$se.fit
  pred_data$CI_lower <- preds$fit - 1.96 * preds$se.fit
  pred_data$CI_upper <- preds$fit + 1.96 * preds$se.fit
  
  predictions <- rbind(predictions, pred_data)
}

p_interaction <- ggplot() +
  # Confidence ribbons
  geom_ribbon(data = predictions,
              aes(x = VWC, ymin = CI_lower, ymax = CI_upper, fill = Species),
              alpha = 0.2) +
  # Prediction lines
  geom_line(data = predictions,
            aes(x = VWC, y = CH4_flux_pred, color = Species),
            size = 1.2) +
  # Actual data points
  geom_point(data = analysis_data,
             aes(x = VWC, y = CH4_flux, color = Species),
             alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Species × Environment Interaction",
       subtitle = "Predicted CH₄ flux across moisture gradient with 95% CI",
       x = "Volumetric Water Content (standardized)",
       y = expression(CH[4]~Flux~(nmol~m^{-2}~s^{-1}))) +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "right")

print(p_interaction)

# ============================================
# 6. COMBINED SUMMARY FIGURE
# ============================================

# Combine key plots
combined_fig <- (p_euler | p_coef) / (p_heatmap | p_interaction) +
  plot_annotation(
    title = "Environmental vs Species Control of CH₄ Flux: Complete Analysis",
    subtitle = paste0("Key finding: Species × Environment interactions explain ", 
                      round((r2_interaction - r2_full) * 100, 1), 
                      "% of variance, suggesting species-specific responses"),
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

print(combined_fig)

# ============================================
# 7. SUMMARY STATISTICS TABLE
# ============================================

cat("\n\nSUMMARY: DIRECTION OF EFFECTS\n")
cat("================================\n\n")

# Environmental effects
cat("ENVIRONMENTAL EFFECTS (standardized coefficients):\n")
env_summary <- coef_data %>%
  filter(category == "Environment") %>%
  select(term_clean, estimate, p.value) %>%
  mutate(
    direction = ifelse(estimate > 0, "↑ Increases flux", "↓ Decreases flux"),
    significance = ifelse(p.value < 0.05, "***", ifelse(p.value < 0.1, "*", "ns"))
  )
print(env_summary)

cat("\n\nSPECIES EFFECTS (relative to reference species):\n")
species_summary <- coef_data %>%
  filter(category == "Species") %>%
  select(term_clean, estimate, p.value) %>%
  mutate(
    direction = ifelse(estimate > 0, "Higher flux", "Lower flux"),
    significance = ifelse(p.value < 0.05, "***", ifelse(p.value < 0.1, "*", "ns"))
  ) %>%
  arrange(desc(estimate))
print(species_summary)

# Save plots
# ggsave("variance_decomposition_enhanced.png", plot = p_euler, width = 10, height = 8, dpi = 300)
# ggsave("coefficient_plot.png", plot = p_coef, width = 10, height = 8, dpi = 300)
# ggsave("species_heatmap.png", plot = p_heatmap, width = 10, height = 8, dpi = 300)
# ggsave("interaction_plot.png", plot = p_interaction, width = 10, height = 8, dpi = 300)
# ggsave("combined_analysis.png", plot = combined_fig, width = 16, height = 12, dpi = 300)