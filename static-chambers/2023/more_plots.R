library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(broom)
library(patchwork)
library(viridis)
library(car)

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

# Prepare 2023 data
data_2023 <- ymf2023 %>%
  select(
    Species.Code,
    DBH = DBH..cm.,
    Air_temp = air_temp_C,
    Stem_temp = stem_temp_C,
    Soil_temp = Soil.Temp...C.,
    VWC = vwc_mean,
    CH4_flux = CH4_best.flux
  ) %>%
  mutate(
    Year = "2023",
    Species_Latin = species_mapping[Species.Code]
  ) %>%
  drop_na(CH4_flux)

# Prepare 2021 data - use only 125cm height
data_2021 <- ymf2021 %>%
  select(
    Species.Code = species_id,
    DBH = dbh,
    Air_temp = Temp_Air_125cm,
    Soil_temp = SoilTemp_mean,
    VWC = VWC_mean,
    CH4_flux = CH4_best.flux_125cm
  ) %>%
  mutate(
    Year = "2021",
    Species_Latin = species_mapping[Species.Code],
    Stem_temp = NA
  ) %>%
  drop_na(CH4_flux) %>%
  filter(!is.nan(CH4_flux))

# Combine datasets
combined_data <- bind_rows(data_2023, data_2021) %>%
  filter(!is.na(Species_Latin)) %>%
  group_by(Species_Latin) %>%
  filter(n() > 3) %>%
  ungroup() %>%
  mutate(Species_Latin = reorder(Species_Latin, CH4_flux, mean, na.rm = TRUE))

# Standardize environmental variables
env_vars <- c("DBH", "Air_temp", "Stem_temp", "Soil_temp", "VWC")
for(var in env_vars) {
  if(var %in% names(combined_data)) {
    combined_data[[paste0(var, "_std")]] <- scale(combined_data[[var]], center = TRUE, scale = TRUE)[,1]
  }
}

# Print data summary
cat("========================================\n")
cat("COMBINED DATASET SUMMARY\n")
cat("========================================\n")
cat(sprintf("Total observations: %d\n", nrow(combined_data)))
cat(sprintf("2021 observations: %d\n", sum(combined_data$Year == "2021")))
cat(sprintf("2023 observations: %d\n", sum(combined_data$Year == "2023")))
cat(sprintf("Number of species: %d\n", n_distinct(combined_data$Species_Latin)))
cat("\n")

# ========================================
# FIT ALL MODELS
# ========================================

# Simple models
model_env_only <- lm(CH4_flux ~ DBH_std + Air_temp_std + Soil_temp_std + VWC_std, 
                     data = combined_data)
model_species_only <- lm(CH4_flux ~ Species_Latin, data = combined_data)

# Full additive model
model_full <- lm(CH4_flux ~ DBH_std + Air_temp_std + Soil_temp_std + VWC_std + Species_Latin, 
                 data = combined_data)

# Interaction model
model_interaction <- lm(CH4_flux ~ (DBH_std + Air_temp_std + Soil_temp_std + VWC_std) * Species_Latin, 
                        data = combined_data)

# Get R-squared values
r2_env <- summary(model_env_only)$r.squared
r2_species <- summary(model_species_only)$r.squared
r2_full <- summary(model_full)$r.squared
r2_interaction <- summary(model_interaction)$r.squared

# Get reference species
ref_species <- levels(factor(combined_data$Species_Latin))[1]

# ========================================
# VARIANCE PARTITIONING - METHOD 2
# ========================================

variance_method2 <- data.frame(
  Component = factor(c("Environment", "Species", "Env × Species\nInteraction", "Unexplained"),
                     levels = c("Unexplained", "Env × Species\nInteraction", "Species", "Environment")),
  Variance = c(
    max(0, r2_full - r2_species) * 100,      # Unique environment
    max(0, r2_full - r2_env) * 100,          # Unique species  
    max(0, r2_interaction - r2_full) * 100,  # Interaction
    (1 - r2_interaction) * 100               # Unexplained
  )
)

# ========================================
# VARIANCE PARTITIONING - METHOD 3
# ========================================

anova_type3 <- Anova(model_full, type = "III")
ss_total <- sum(anova_type3$`Sum Sq`)
ss_env <- sum(anova_type3$`Sum Sq`[grepl("DBH|Air_temp|Soil_temp|VWC", rownames(anova_type3))])
ss_species <- anova_type3$`Sum Sq`[rownames(anova_type3) == "Species_Latin"]
ss_residual <- anova_type3$`Sum Sq`[rownames(anova_type3) == "Residuals"]

variance_method3 <- data.frame(
  Component = factor(c("Environment", "Species", "Unexplained"),
                     levels = c("Unexplained", "Species", "Environment")),
  Variance = c(
    ss_env / ss_total * 100,
    ss_species / ss_total * 100,
    ss_residual / ss_total * 100
  )
)

# ========================================
# CREATE PLOTS
# ========================================

# Extract coefficients for effects plot
all_coef <- tidy(model_full, conf.int = TRUE) %>%
  filter(term != "(Intercept)")

# Process environmental coefficients
env_coef <- all_coef %>%
  filter(grepl("_std", term)) %>%
  mutate(
    category = "Environment",
    term_clean = case_when(
      grepl("DBH", term) ~ "Tree diameter (DBH)",
      grepl("Air_temp", term) ~ "Air temperature",
      grepl("Soil_temp", term) ~ "Soil temperature",
      grepl("VWC", term) ~ "Soil moisture (VWC)",
      TRUE ~ term
    ),
    significant = p.value < 0.05
  )

# Process species coefficients
species_coef <- all_coef %>%
  filter(grepl("Species_Latin", term)) %>%
  mutate(
    category = "Species",
    term_clean = gsub("Species_Latin", "", term),
    significant = p.value < 0.05
  )

# Combine all coefficients
all_effects <- bind_rows(env_coef, species_coef)

# Calculate species statistics for ridgeline
species_stats <- combined_data %>%
  group_by(Species_Latin) %>%
  summarise(
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    se_flux = sd(CH4_flux, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

# 1. CREATE RIDGELINE PLOT
p_ridgeline <- ggplot(combined_data, aes(x = CH4_flux, y = Species_Latin)) +
  geom_density_ridges(
    fill = "lightgray",
    alpha = 0.6,
    scale = 1.5,
    rel_min_height = 0.01,
    color = "gray30",
    linewidth = 0.3
  ) +
  geom_point(
    aes(color = VWC),
    position = position_jitter(width = 0, height = 0.15),
    size = 1.5,
    alpha = 0.5
  ) +
  geom_pointrange(
    data = species_stats,
    aes(x = mean_flux, 
        xmin = mean_flux - se_flux, 
        xmax = mean_flux + se_flux,
        y = Species_Latin),
    color = "black",
    size = 0.75,
    linewidth = 0.75,
    fatten = 2,
    inherit.aes = FALSE
  ) +
  scale_color_viridis(
    name = "VWC (%)",
    limits = c(0, max(combined_data$VWC, na.rm = TRUE)),
    breaks = c(0, 25, 50),
    guide = guide_colorbar(
      position = "bottom",
      barwidth = 10,
      barheight = 0.5,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.01),
    breaks = c(-0.1, 0, 0.1, 1),
    labels = c("-0.1", "0", "0.1", "1")
  ) +
  labs(
    x = expression(CH[4]~Flux~(nmol~m^{-2}~s^{-1})),
    y = ""
  ) +
  theme_ridges(grid = TRUE) +
  theme(
    axis.title.x = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    panel.grid.minor.x = element_line(color = "gray95", size = 0.2),
    plot.margin = margin(5, 2, 5, 5, "pt")
  )

# 2. CREATE VARIANCE PLOT - METHOD 2 (with y-axis on right)
p_var_method2 <- ggplot(variance_method2, aes(x = Component, y = Variance, fill = Component)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.6) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Environment" = "#2E7D32",
                               "Species" = "#1976D2",
                               "Env × Species\nInteraction" = "#FBC02D",
                               "Unexplained" = "#9A8C98"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 100),
                     expand = c(0, 0),
                     position = "right") +
  labs(x = "",
       y = "Variance (%)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y.right = element_text(size = 9),
    plot.margin = margin(2, 2, 2, 2, "pt")
  )

# Skip Method 3 plot creation since we're only using Method 2

# 4. CREATE EFFECTS PLOT
p_effects <- ggplot(all_effects, aes(x = estimate, y = reorder(term_clean, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = category),
                 height = 0, size = 0.6, alpha = 0.8) +
  geom_point(aes(color = category, shape = significant), size = 2.5) +
  scale_color_manual(values = c("Environment" = "#2E7D32", 
                                "Species" = "#1976D2"),
                     name = "Effect type") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("p ≥ 0.05", "p < 0.05"),
                     name = "Significance") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Standardized effect size (nmol/m²/s)",
       y = "") +
  theme_bw() +
  theme(
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.1, "cm"),
    plot.margin = margin(5, 5, 5, 2, "pt")
  )

# ========================================
# COMBINE ALL PLOTS
# ========================================

combined_plot <- p_ridgeline | (p_var_method2 / p_effects) +
  plot_layout(widths = c(1.2, 1), heights = c(1, 2))

print(combined_plot)

# ========================================
# PRINT SUMMARY STATISTICS
# ========================================

cat("\n\nVARIANCE PARTITIONING RESULTS\n")
cat("========================================\n")
cat("\nMETHOD 2 (With Interactions):\n")
print(variance_method2)
cat(sprintf("\nTotal R² with interactions: %.1f%%\n", r2_interaction * 100))

cat("\nMETHOD 3 (Type III SS):\n")
print(variance_method3)

cat("\n\nMODEL R-SQUARED VALUES\n")
cat("========================================\n")
cat(sprintf("Environment only: %.1f%%\n", r2_env * 100))
cat(sprintf("Species only: %.1f%%\n", r2_species * 100))
cat(sprintf("Full additive model: %.1f%%\n", r2_full * 100))
cat(sprintf("With interactions: %.1f%%\n", r2_interaction * 100))

# Save plot if desired
# ggsave("ch4_flux_complete_analysis.png", plot = combined_plot, width = 16, height = 10, dpi = 300)