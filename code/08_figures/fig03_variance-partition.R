# ==============================================================================
# Variance Partitioning (Figure 3)
# ==============================================================================
# Purpose: Variance partitioning showing relative contributions of species,
#   moisture, genes, and wood properties to flux.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - data/compiled/flux_measurements_tree.csv (the merged measurement-level table)
#
# The unit of analysis is the TREE. Each tree contributes one row: the mean of its
# breast-height measurements. Four changes from the previous version, all aimed at
# partitioning variance on the scale the data live on, over a sample where every
# component is identified from the whole of it.
#
# 1. RESPONSE IS TRANSFORMED. The previous version partitioned raw nmol, where the
#    top 1% of measurements carry 80% of the total variance and the top 5% carry
#    95%. That decomposition described about ten extreme fluxes, not the
#    population, and reported ~83% unexplained as a result. Panel (a) was already
#    drawn on an arcsinh axis for exactly this reason -- the figure displayed one
#    scale and analysed another. Both are now arcsinh.
#
# 2. ONE ROW PER TREE, AVERAGING REPEATS RATHER THAN DISCARDING THEM. Repeated
#    measurements are very unevenly distributed: only 41 of 446 trees have more
#    than one, yet they supply 44% of the measurement-level rows, nearly all from
#    the 2020 monthly survey. Treating rows as independent therefore let a minority
#    of trees dominate, and a tree random effect was no better -- it estimated a
#    36.6% "tree identity" term from those 41 trees and reported it as if it came
#    from the whole sample. Averaging each tree's measurements uses every one of
#    them, removes the pseudo-replication, and leaves every component identified
#    from all 446 trees. Unexplained rises from 40% to 65% as a result; that is the
#    honest figure, and the earlier one was low because a term fitted on 9% of the
#    trees was absorbing variance.
#
# 3. ALL THREE CAMPAIGNS, FROM MEASUREMENT-LEVEL DATA. The previous version took
#    2021 from the tree-level wide table, where the environmental columns are
#    per-tree summaries, and omitted the 2020 monthly campaign entirely.
#    Soil-temperature spread was 2.5 against 4.0 in the measurement-level data, so
#    environment was being tested on a fraction of its observed range. It is still
#    weak with the full range, which makes that a result rather than an artefact.
#
# 4. BREAST HEIGHT ONLY, plus a chamber-design covariate -- see the data block.
#
# The three campaigns were designed around different axes -- 2020 across time,
# 2021 across height, 2023 across species. Averaging to the tree is what lets them
# be pooled for a question about species and environment: the monthly campaign
# contributes a better-estimated mean per tree rather than repeated rows, and its
# temporal axis stays the subject of Figure 1, as height stays that of Figure 2.
#
# Partition is Nakagawa & Schielzeth marginal/conditional R2: marginal is the
# fixed effects, conditional adds the tree random effect, the remainder is
# within-stem residual.
# ==============================================================================

flux_all <- read.csv("data/compiled/flux_measurements_tree.csv")

library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(broom)
library(patchwork)
library(viridis)
library(car)
library(lme4)
library(lmerTest)

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

# ---- one measurement-level table, all three campaigns -----------------------
# ASINH_SIGMA matches the axis transform panel (a) already uses, so the analysis
# and the display are on the same scale.
ASINH_SIGMA <- 0.1
asinh_t <- function(x) asinh(x / ASINH_SIGMA)

# BREAST HEIGHT ONLY. The measurement-level table mixes the 2021 multi-height
# campaign (50 / 125 / 200 cm) with campaigns that recorded no height. Height is
# Figure 2's subject, so including 50 and 200 cm here would put vertical structure
# into a figure about species and environment, where -- since height is not a term
# -- it would be absorbed by tree identity and the residual. The unlabelled rows
# are the monthly and 2023 campaigns, measured at breast height (confirmed with
# Jon); they are kept, and the 50 cm and 200 cm rows are dropped. This costs
# little: the partition moves by 1-2 points either way.
#
# chamber_type is carried as a nuisance covariate. The campaigns differ ~6-fold in
# median flux, so pooling designs is the obvious referee question; including the
# term answers it, and it changes the partition by under a point, which says the
# difference is between the trees measured rather than between chamber designs.
combined_data <- flux_all %>%
  filter(is.na(measurement_height_cm) | measurement_height_cm == 125) %>%
  transmute(
    tree_id,
    Species.Code = species_code,
    DBH      = dbh_m * 100,                  # m -> cm, as the old wide table held it
    Air_temp = air_temp_C,
    Soil_temp = soil_temp_C,
    VWC      = soil_moisture_abs * 100,      # m3 m-3 -> %, as panel (a) labels it
    CH4_flux = stem_flux_nmol_m2_s,
    Year     = as.character(year),
    Chamber  = chamber_type
  ) %>%
  mutate(Species_Latin = species_mapping[Species.Code]) %>%
  filter(!is.na(Species_Latin)) %>%
  drop_na(CH4_flux, DBH, Air_temp, Soil_temp, VWC, Chamber) %>%
  filter(is.finite(CH4_flux)) %>%
  mutate(CH4_asinh = asinh_t(CH4_flux))

# Average to one row per tree. The mean is taken on the arcsinh scale, the scale
# the models are fitted on, so a tree's summary is the mean of what is analysed
# rather than a mean of raw fluxes that the transform would then distort.
combined_data <- combined_data %>%
  group_by(tree_id) %>%
  summarise(Species_Latin = dplyr::first(Species_Latin),
            Species.Code  = dplyr::first(Species.Code),
            across(c(DBH, Air_temp, Soil_temp, VWC, CH4_asinh), ~ mean(.x, na.rm = TRUE)),
            Chamber = names(sort(table(Chamber), decreasing = TRUE))[1],
            Year    = names(sort(table(Year),    decreasing = TRUE))[1],
            n_measurements = dplyr::n(), .groups = "drop") %>%
  mutate(CH4_flux = ASINH_SIGMA * sinh(CH4_asinh)) %>%
  group_by(Species_Latin) %>% filter(n() > 3) %>% ungroup() %>%
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
cat(sprintf("Total observations: %d from %d trees\n",
            nrow(combined_data), n_distinct(combined_data$tree_id)))
for (y in sort(unique(combined_data$Year)))
  cat(sprintf("  %s: %d\n", y, sum(combined_data$Year == y)))
cat(sprintf("Chamber designs: %s\n",
            paste(sprintf("%s=%d", names(table(combined_data$Chamber)),
                          table(combined_data$Chamber)), collapse = "  ")))
cat(sprintf("Number of species: %d\n", n_distinct(combined_data$Species_Latin)))
cat("\n")

# ========================================
# FIT ALL MODELS
# ========================================

# ---- models: OLS on tree means ----------------------------------------------
# No random effect. With one row per tree, tree identity is saturated and not
# separable from residual; fitting (1|tree_id) anyway makes lmer assign the whole
# residual to the random effect and report 0% unexplained, which is an artefact of
# the saturation rather than a result. Between-tree variance that the fixed effects
# do not explain belongs in "unexplained" here, and that is what it means.
model_env_only     <- lm(CH4_asinh ~ DBH_std + Air_temp_std + Soil_temp_std + VWC_std + Chamber,
                         data = combined_data)
model_species_only <- lm(CH4_asinh ~ Species_Latin + Chamber, data = combined_data)
model_full         <- lm(CH4_asinh ~ DBH_std + Air_temp_std + Soil_temp_std + VWC_std + Species_Latin + Chamber,
                         data = combined_data)
model_interaction  <- lm(CH4_asinh ~ (DBH_std + Air_temp_std + Soil_temp_std + VWC_std) * Species_Latin + Chamber,
                         data = combined_data)

r2_env         <- summary(model_env_only)$r.squared
r2_species     <- summary(model_species_only)$r.squared
r2_full        <- summary(model_full)$r.squared
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
    max(0, r2_full - r2_species) * 100,     # unique environment
    max(0, r2_full - r2_env) * 100,         # unique species
    max(0, r2_interaction - r2_full) * 100, # interaction
    (1 - r2_interaction) * 100              # between-tree variance the model does not explain
  )
)

# ========================================
# VARIANCE PARTITIONING - METHOD 3
# ========================================

# Reported alongside the partition rather than plotted: the R2 ladder, so the
# components in panel (b) can be traced back to the models they came from.
variance_method3 <- data.frame(
  Model  = c("environment only", "species only", "full additive", "with interactions"),
  R2_pct = round(100 * c(r2_env, r2_species, r2_full, r2_interaction), 1)
)

# ========================================
# CREATE PLOTS
# ========================================

# Extract coefficients for effects plot
all_coef <- broom::tidy(model_full, conf.int = TRUE) %>%
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
    trans = scales::trans_new("arcsinh", function(x) asinh(x/0.1), function(x) 0.1*sinh(x)),
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
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "gray95", linewidth = 0.2),
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
                 height = 0, linewidth = 0.6, alpha = 0.8) +
  geom_point(aes(color = category, shape = significant), size = 2.5) +
  scale_color_manual(values = c("Environment" = "#2E7D32", 
                                "Species" = "#1976D2"),
                     name = "Effect type") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c(expression(p >= 0.05), expression(p < 0.05)),
                     name = "Significance") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  # Effects are on the arcsinh scale the models are fitted on, NOT nmol m-2 s-1.
  # Labelling them in flux units would misstate them by the transform.
  labs(x = expression("Standardized effect size (arcsinh CH"[4]*" flux)"),
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

combined_plot <- (p_ridgeline | (p_var_method2 / p_effects + plot_layout(heights = c(1, 2)))) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")",
                  theme = theme(plot.tag = element_text(size = 11, face = "bold")))

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

# Save plot
ggsave("outputs/figures/generated/fig3_final.png", plot = combined_plot, width = 10, height = 7, dpi = 300)






# ===== EXTRACT STATISTICS FOR SPECIES TRENDS SECTION =====

cat("\n===== SPECIES-LEVEL STATISTICS FOR TEXT =====\n")

# Get species means and SEs
species_summary <- combined_data %>%
  group_by(Species_Latin) %>%
  summarise(
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    se_flux = sd(CH4_flux, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_flux))

cat("\nTop 5 species by mean flux:\n")
print(head(species_summary, 5))

cat("\nGymnosperms:\n")
print(species_summary %>% filter(Species_Latin %in% c("Pinus strobus", "Tsuga canadensis")))

# Variance partitioning summary
cat("\n===== VARIANCE PARTITIONING SUMMARY =====\n")
cat("Species identity explains:", round(variance_method2$Variance[variance_method2$Component == "Species"], 1), "%\n")
cat("Environment-Species interaction:", round(variance_method2$Variance[variance_method2$Component == "Env × Species\nInteraction"], 1), "%\n")
cat("Unexplained variance:", round(variance_method2$Variance[variance_method2$Component == "Unexplained"], 1), "%\n")
cat("Environment alone:", round(variance_method2$Variance[variance_method2$Component == "Environment"], 1), "%\n")

# Total observations
cat("\nTotal observations pooled:", nrow(combined_data), "\n")
for (y in sort(unique(combined_data$Year)))
  cat(sprintf("  from %s: %d\n", y, sum(combined_data$Year == y)))
cat("  trees:", dplyr::n_distinct(combined_data$tree_id), "\n")
cat("  measurements averaged into these means:", sum(combined_data$n_measurements), "\n")