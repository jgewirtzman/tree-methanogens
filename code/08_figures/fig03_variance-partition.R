# ==============================================================================
# fig03_variance-partition.R -- Figure 3: species and environment in stem CH4 flux
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) flux by species (ridgelines; points coloured by soil VWC; mean +/- SE),
#       arcsinh axis (code/lib/flux_transform.R)
#   (b) variance partition of tree-mean flux: unique environment, unique species,
#       environment x species interaction, unexplained (from the R2 of four OLS fits)
#   (c) standardized coefficients of the full additive model
# Input
#   outputs/data/flux_measurements_tree.csv  every stem deployment of the three
#                                            campaigns (in_rf_training flags the model's
#                                            subset; all rows are used here)
# Output
#   outputs/figures/generated/fig3_final.png; the partition and R2 ladder to the console
#
# Rules
# - Unit of analysis: the TREE, one row per tree, the mean of its deployments on the
#   arcsinh scale (flux_asinh()), the scale the models are fitted on and panel (a)
#   draws. Repeat measurements are averaged rather than treated as independent rows,
#   so the few trees of the monthly survey do not dominate.
# - Breast height only: deployments with no height (monthly and 2023 surveys) and the
#   125 cm deployments of the 2021 multi-height survey. Height is Figure 2's subject.
# - Growing season only (May-September), so tree means of the monthly survey are
#   comparable with the single-visit campaigns.
# - Chamber design is a nuisance covariate in every model.
# - Species with more than three trees.
# - OLS on tree means, no tree random effect: with one row per tree it is not
#   separable from the residual. Between-tree variance the fixed effects do not
#   explain is "unexplained".
# ==============================================================================

library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(broom)
library(patchwork)
library(viridis)
source("code/lib/tables_read.R")
source("code/lib/flux_transform.R")
source("code/lib/geometry.R")   # GYMNOSPERMS

flux_all <- read.csv("outputs/data/flux_measurements_tree.csv")
species_mapping <- SPECIES_NAMES   # data/tables/species.csv

GROWING_SEASON <- 5:9

# Every filter below is counted and printed, so a drop is visible.
cat(sprintf("Figure 3 input: %d deployments (%d outside the model's training set)\n",
            nrow(flux_all), sum(!flux_all$in_rf_training)))
.step <- function(d, label) { cat(sprintf("  %-46s %5d rows\n", label, nrow(d))); d }
combined_data <- flux_all %>%
  .step("all deployments") %>%
  filter(is.na(measurement_height_cm) | measurement_height_cm == 125) %>% .step("breast height") %>%
  filter(month %in% GROWING_SEASON) %>% .step("growing season (May-Sep)") %>%
  transmute(
    tree_id,
    Species.Code = species_code,
    DBH      = dbh_m * 100,                  # m -> cm
    Air_temp = air_temp_C,                   # measured at the deployment (flux_stem air_temp_C)
    Soil_temp = soil_temp_C,
    VWC      = soil_moisture_abs * 100,      # m3 m-3 -> %, as panel (a) labels it
    CH4_flux = stem_flux_nmol_m2_s,
    Year     = as.character(year),
    Chamber  = chamber_type
  ) %>%
  mutate(Species_Latin = species_mapping[Species.Code]) %>%
  filter(!is.na(Species_Latin)) %>%
  .step("species identified (species.csv)") %>%
  drop_na(CH4_flux, DBH, Air_temp, Soil_temp, VWC, Chamber) %>% .step("complete covariates") %>%
  filter(is.finite(CH4_flux)) %>%
  mutate(CH4_asinh = flux_asinh(CH4_flux))

# Average to one row per tree, on the arcsinh scale. Chamber and year: the tree's
# most frequent value.
combined_data <- combined_data %>%
  group_by(tree_id) %>%
  summarise(Species_Latin = dplyr::first(Species_Latin),
            Species.Code  = dplyr::first(Species.Code),
            across(c(DBH, Air_temp, Soil_temp, VWC, CH4_asinh), ~ mean(.x, na.rm = TRUE)),
            Chamber = names(sort(table(Chamber), decreasing = TRUE))[1],
            Year    = names(sort(table(Year),    decreasing = TRUE))[1],
            n_measurements = dplyr::n(), .groups = "drop") %>%
  mutate(CH4_flux = FLUX_COFACTOR * sinh(CH4_asinh)) %>%
  group_by(Species_Latin) %>% filter(n() > 3) %>% ungroup() %>%
  mutate(Species_Latin = reorder(Species_Latin, CH4_flux, mean, na.rm = TRUE))

# Standardize environmental variables
for (var in c("DBH", "Air_temp", "Soil_temp", "VWC"))
  combined_data[[paste0(var, "_std")]] <- scale(combined_data[[var]], center = TRUE, scale = TRUE)[,1]

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

# ---- models: OLS on tree means (no tree random effect; see header) -----------
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

# ========================================
# VARIANCE PARTITIONING
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

# The R2 ladder, printed (not plotted) so the components in panel (b) can be traced
# back to the models they came from.
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
    name = "Soil VWC (%)",
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
    trans = flux_asinh_trans(),
    breaks = c(-0.1, 0, 0.1, 1),
    labels = c("\u22120.1", "0", "0.1", "1")
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = ""
  ) +
  theme_ridges(grid = TRUE) +
  theme(
    axis.title.x = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 9), legend.text = element_text(size = 8),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "gray95", linewidth = 0.2),
    plot.margin = margin(5, 2, 5, 5, "pt")
  )

# 2. CREATE VARIANCE PLOT (y-axis on the right)
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

# 3. CREATE EFFECTS PLOT
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
  scale_y_discrete(labels = function(x) parse(text = ifelse(x %in% all_effects$term_clean[all_effects$category == "Species"],
                                                            paste0("italic('", x, "')"), paste0("'", x, "'")))) +
  scale_x_continuous(labels = function(b) sub("-", "\u2212", format(b, drop0trailing = TRUE, trim = TRUE))) +
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
                  theme = theme(plot.tag = element_text(size = 11, face = "bold"))) &
  theme(plot.tag = element_text(size = 11, face = "bold"))

# ========================================
# PRINT SUMMARY STATISTICS
# ========================================

cat("\n\nVARIANCE PARTITIONING RESULTS\n")
cat("========================================\n")
cat("\nPartition (with interactions):\n")
print(variance_method2)
cat(sprintf("\nTotal R² with interactions: %.1f%%\n", r2_interaction * 100))

cat("\nR2 ladder:\n")
print(variance_method3)

cat("\n\nMODEL R-SQUARED VALUES\n")
cat("========================================\n")
cat(sprintf("Environment only: %.1f%%\n", r2_env * 100))
cat(sprintf("Species only: %.1f%%\n", r2_species * 100))
cat(sprintf("Full additive model: %.1f%%\n", r2_full * 100))
cat(sprintf("With interactions: %.1f%%\n", r2_interaction * 100))

# Save plot
ggsave("outputs/figures/generated/fig3_final.png", plot = combined_plot, width = 8.2, height = 6.2, dpi = 300, bg = "white")   # near print size

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
print(species_summary %>% filter(Species_Latin %in% GYMNOSPERMS))

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