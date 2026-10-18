# ==============================================================================
# fig02_height-flux.R -- Figure 2: stem CH4 flux with height, by species
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) flux at 0.5, 1.25 and 2 m per species: half-boxplots and deployments, one
#       shared arcsinh axis (code/lib/flux_transform.R)
#   (b) height effect per species: lmer(flux ~ height + (1 | tree)), slope per m
#       with 95% CI; p from the normal approximation of t
#   (c) species means of soil VWC and soil mcrA (organic, mineral), z-scored
# The 2021 multi-height deployments at 50, 125 and 200 cm are used; Kalmia latifolia
# (measured at 75 cm and the root crown, no height_cm for the latter) is left out.
# Inputs (data/tables/, via code/lib/tables_read.R)
#   flux_stem.csv (2021_multiheight), trees.csv, species.csv
#   survey_2021_trees(soil = TRUE): per-tree soil VWC (code/lib/soil_fill.R) and soil
#     mcrA (ddpcr.csv + samples.csv)
# Output
#   outputs/figures/generated/fig2_final.png; height-effect statistics to the console
# ==============================================================================
library(ggplot2)
library(dplyr)
library(lme4)
library(tidyr)
library(scales)
library(patchwork)
library(gghalves)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

source("code/lib/tables_read.R")
source("code/lib/flux_transform.R")
# 2021 multi-height deployments; species and landscape position (plot) from trees.csv
final_dataset <- flux_stem_campaign("2021_multiheight") %>%
  transmute(UniqueID = measurement_id, tree_id, species = species_code,
            plot = trees_table()$landscape_position[match(tree_id, trees_table()$tree_id)],
            measurement_height = height_cm, CH4_best.flux = CH4_flux_nmol_m2_s)

merged_data <- survey_2021_trees(soil = TRUE)   # 2021 survey per tree (data/tables/)


# ==============================================================================
# 2. PREPARE DATA
# ==============================================================================

# Publication colors -- red for positive height effect, blue for negative
colors_pub <- c(
  "Significant Negative" = "#4575B4",
  "Significant Positive" = "#D73027",
  "Non-significant" = "grey60"
)

shapes_pub <- c(
  "Significant Negative" = 16,
  "Significant Positive" = 16,
  "Non-significant" = 1
)

boxplot_color <- "#756BB1"

# Deployments with a height (the root crown has none) and flux; Kalmia left out
plot_data_top <- final_dataset %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux) & !is.na(height_numeric)) %>%
  mutate(
    height_m = height_numeric / 100
  ) %>%
  filter(species != "KALA") %>%
  mutate(
    tree_unique = paste(plot, tree_id, sep = "_"),
    species_latin = species_name(species)
  )

# Calculate tree counts per species
tree_counts <- plot_data_top %>%
  group_by(species, species_latin) %>%
  summarise(n_trees = n_distinct(tree_unique), .groups = 'drop') %>%
  mutate(
    species_label_with_n = paste0(species_latin, "\n(", n_trees, ")"),
    species_label_no_n = species_latin
  )

# Add labels to plot data
plot_data_top <- plot_data_top %>%
  left_join(tree_counts %>% dplyr::select(species, species_label_with_n, species_label_no_n),
            by = "species") %>%
  mutate(species_label = factor(species_label_with_n, levels = sort(unique(species_label_with_n))))

# ==============================================================================
# 3. MIXED-EFFECTS HEIGHT ANALYSIS
# ==============================================================================

analysis_data <- plot_data_top

# Species with fewer than 3 trees measured at more than one height are not tested
test_species_height_effect <- function(species_name, data) {
  species_data <- data %>% filter(species == species_name)

  n_obs <- nrow(species_data)
  n_trees <- length(unique(species_data$tree_unique))
  n_trees_multi_height <- species_data %>%
    group_by(tree_unique) %>%
    summarise(n_heights = n_distinct(height_numeric)) %>%
    filter(n_heights > 1) %>%
    nrow()

  if(n_trees_multi_height < 3) {
    return(data.frame(
      species = species_name, n_obs = n_obs, n_trees = n_trees,
      n_trees_multi_height = n_trees_multi_height,
      height_coef = NA, height_se = NA, height_t = NA, height_p = NA,
      model_type = "insufficient_data"
    ))
  }

  model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = species_data)
  coef_table <- summary(model)$coefficients
  data.frame(
    species = species_name, n_obs = n_obs, n_trees = n_trees,
    n_trees_multi_height = n_trees_multi_height,
    height_coef = coef_table["height_numeric", "Estimate"],
    height_se = coef_table["height_numeric", "Std. Error"],
    height_t = coef_table["height_numeric", "t value"],
    height_p = 2 * (1 - pnorm(abs(coef_table["height_numeric", "t value"]))),
    model_type = "mixed_effects"
  )
}

species_list <- sort(unique(plot_data_top$species))
species_results <- do.call(rbind, lapply(species_list, test_species_height_effect, analysis_data))

# ==============================================================================
# 4. PANEL A -- species half-boxplots, one faceted plot
# ==============================================================================
# One shared arcsinh flux axis (cofactor 0.1, as Figs 1 and 3); dashed line = zero
# net flux. Axis text only on the outer panels, species names as strips.
trans_symlog_fig2 <- flux_asinh_trans()
shared_brks_fig2  <- c(0, 0.1, 1)   # 3 well-separated ticks (no label crowding near zero)
shared_lims_fig2  <- range(c(plot_data_top$CH4_best.flux, shared_brks_fig2), na.rm = TRUE)
shared_lab_fig2   <- function(x) sub("-", "\u2212", format(x, drop0trailing = TRUE, trim = TRUE, scientific = FALSE))

lab_df <- plot_data_top %>% distinct(species, species_label)
top_df <- plot_data_top %>%
  mutate(species_label = factor(species_label, levels = lab_df$species_label[match(species_list, lab_df$species)]))
p_top_gg <- ggplot(top_df, aes(x = height_m, y = CH4_best.flux)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5, alpha = 0.8) +
  geom_half_boxplot(aes(group = height_m), alpha = 0.6, outlier.shape = NA,
                    fill = boxplot_color, color = "gray30", side = "r", linewidth = 0.3) +
  geom_jitter(alpha = 0.4, size = 0.9, color = "gray20", position = position_jitter(width = 0.1, height = 0)) +
  facet_wrap(~ species_label, ncol = 5) +
  coord_flip() +
  scale_x_continuous(breaks = c(0.5, 1.25, 2), labels = c("0.5", "1.25", "2"), expand = expansion(mult = c(0.13, 0.13))) +
  scale_y_continuous(trans = trans_symlog_fig2, breaks = shared_brks_fig2, limits = shared_lims_fig2,
                     labels = shared_lab_fig2, expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Height (m)", y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*"),  arcsinh scale")) +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(size = 7, face = "italic", margin = margin(b = 1, t = 1), lineheight = 0.9),
        axis.text = element_text(size = 7), axis.title = element_text(size = 9),
        panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray40", fill = NA, linewidth = 0.4),
        panel.spacing.x = unit(0.25, "lines"), panel.spacing.y = unit(0.35, "lines"),
        plot.margin = margin(t = 0, r = 2, b = 0, l = 2))
p_top_gg <- wrap_elements(full = p_top_gg)   # not aligned to (b)/(c), so it can use the full width

# ==============================================================================
# 5. PANEL B — Height effect coefficients (middle)
# ==============================================================================

plot_data_middle <- species_results %>%
  filter(!is.na(height_coef)) %>%
  left_join(tree_counts %>% dplyr::select(species, species_label_with_n, species_label_no_n),
            by = "species") %>%
  mutate(
    species_label = factor(species_label_with_n, levels = sort(unique(species_label_with_n))),
    significant = ifelse(is.na(height_p), FALSE, height_p < 0.05),
    negative_trend = ifelse(is.na(height_coef), FALSE, height_coef < 0),
    color_group = case_when(
      significant & negative_trend ~ "Significant Negative",
      significant & !negative_trend ~ "Significant Positive",
      TRUE ~ "Non-significant"
    )
  ) %>%
  arrange(height_coef)

p_middle <- ggplot(plot_data_middle,
                   aes(x = reorder(species_label, -height_coef),
                       y = height_coef * 100,
                       color = color_group, shape = color_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = height_coef * 100 - 1.96 * height_se * 100,
                    ymax = height_coef * 100 + 1.96 * height_se * 100),
                width = 0.3, linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = colors_pub, name = "") +
  scale_shape_manual(values = shapes_pub, name = "") +
  scale_y_continuous(labels = function(b) sub("-", "\u2212", format(b, drop0trailing = TRUE, trim = TRUE))) +
  labs(
    x = "",
    y = expression(atop("Height effect",
                        "(nmol m"^-2*" s"^-1*" per m)"))
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 8, lineheight = 1.1),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    legend.margin = margin(l = -5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.3, color = "gray85"),
    plot.margin = margin(t = 2, r = 5, b = 0, l = 2)
  )

# ==============================================================================
# 6. PANEL C — Soil conditions heatmap (bottom)
# ==============================================================================

# Trees in the analysis: flux_stem and the per-tree table share the canonical tree_id.
trees_in_analysis <- unique(plot_data_top$tree_id)
no_soil_row <- setdiff(trees_in_analysis, merged_data$tree_id)
cat(sprintf("Panel C: %d trees in the height analysis; %d without a row in the per-tree table%s\n",
            length(trees_in_analysis), length(no_soil_row),
            if (length(no_soil_row)) paste0(": ", paste(no_soil_row, collapse = ", ")) else ""))

# Calculate species-level soil/mcrA means
soil_mcra_data <- merged_data %>%
  filter(tree_id %in% trees_in_analysis) %>%
  mutate(
    species = species_id,
    mcra_organic = ddpcr_mcra_probe_Organic_loose,
    mcra_mineral = ddpcr_mcra_probe_Mineral_loose,
    organic_depth = OrganicDepth_mean,
    mineral_depth = MineralDepth_mean,
    mcra_weighted = case_when(
      !is.na(mcra_organic) & !is.na(mcra_mineral) &
        !is.na(organic_depth) & !is.na(mineral_depth) &
        is.finite(mcra_organic) & is.finite(mcra_mineral) &
        is.finite(organic_depth) & is.finite(mineral_depth) ~
        (mcra_organic * organic_depth + mcra_mineral * mineral_depth) /
        (organic_depth + mineral_depth),
      !is.na(mcra_organic) & is.finite(mcra_organic) &
        (is.na(mcra_mineral) | !is.finite(mcra_mineral)) ~ mcra_organic,
      !is.na(mcra_mineral) & is.finite(mcra_mineral) &
        (is.na(mcra_organic) | !is.finite(mcra_organic)) ~ mcra_mineral,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(species) %>%
  summarise(
    mean_VWC = mean(VWC_mean, na.rm = TRUE),
    mean_mcra_organic = mean(mcra_organic, na.rm = TRUE),
    mean_mcra_mineral = mean(mcra_mineral, na.rm = TRUE),
    mean_mcra_weighted = mean(mcra_weighted, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  left_join(tree_counts %>% dplyr::select(species, species_label_with_n, species_label_no_n),
            by = "species") %>%
  filter(!is.na(species_label_with_n))

# Species order from height effect panel (highest -> lowest slope, matches panels a/b)
species_order <- plot_data_middle %>% arrange(desc(height_coef)) %>% pull(species_label)
species_order_no_n <- plot_data_middle %>% arrange(desc(height_coef)) %>% pull(species_label_no_n)

# Z-score normalize for heatmap, using log-transformed mcrA
heatmap_data <- soil_mcra_data %>%
  filter(species %in% plot_data_middle$species) %>%
  mutate(
    log_mcra_organic = log10(mean_mcra_organic + 1),
    log_mcra_mineral = log10(mean_mcra_mineral + 1),
    log_mcra_weighted = log10(mean_mcra_weighted + 1)
  ) %>%
  dplyr::select(species_label_no_n, mean_VWC, log_mcra_organic, log_mcra_mineral) %>%
  # the depth-weighted mean of the two soil layers was dropped: it adds nothing to them
  pivot_longer(cols = c(mean_VWC, log_mcra_organic, log_mcra_mineral),
               names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(z_score = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    variable = case_when(
      variable == "mean_VWC" ~ "VWC",
      variable == "log_mcra_organic" ~ "log mcrA (Organic)",
      variable == "log_mcra_mineral" ~ "log mcrA (Mineral)",
      variable == "log_mcra_weighted" ~ "log mcrA (Weighted)",
      TRUE ~ variable
    ),
    variable = factor(variable, levels = c("VWC", "log mcrA (Organic)", "log mcrA (Mineral)"))
  ) %>%
  filter(species_label_no_n %in% species_order_no_n) %>%
  mutate(
    species_label_no_n = factor(species_label_no_n, levels = species_order_no_n),
    z_score = case_when(
      is.infinite(z_score) ~ NA_real_,
      is.nan(z_score) ~ NA_real_,
      TRUE ~ z_score
    )
  ) %>%
  filter(!is.na(species_label_no_n))

p_bottom <- ggplot(heatmap_data, aes(x = species_label_no_n, y = variable, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_y_discrete(labels = c("VWC" = "Soil VWC",
                              "log mcrA (Organic)"  = expression("Soil"~italic(mcrA)*", organic"),
                              "log mcrA (Mineral)"  = expression("Soil"~italic(mcrA)*", mineral"))) +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "#4575B4",
                       midpoint = 0, na.value = "grey90",
                       name = "Z-score") +
  labs(x = "", y = "") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(size = 7.5, angle = 35, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 7.5, hjust = 1, lineheight = 0.9),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.margin = margin(l = -5),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 3, l = 5)
  )

# ==============================================================================
# 7. COMBINE PANELS AND SAVE
# ==============================================================================

# Panel labels via patchwork annotation
combined_plot <- (p_top_gg / p_middle / p_bottom) +
  plot_layout(heights = c(4.2, 1.2, 1.2)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")",
                  theme = theme(plot.tag = element_text(size = 11, face = "bold"))) &
  theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave("outputs/figures/generated/fig2_final.png", plot = combined_plot, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
cat("Wrote outputs/figures/generated/fig2_final.png\n")


# ==============================================================================
# 8. STATISTICS — HEIGHT EFFECTS
# ==============================================================================

cat("\n===== HEIGHT EFFECT STATISTICS =====\n")

# Betula alleghaniensis detail
beal_data <- analysis_data %>% filter(species == "BEAL")
beal_heights <- beal_data %>%
  group_by(height_numeric) %>%
  summarise(
    mean_flux = mean(CH4_best.flux, na.rm = TRUE),
    se_flux = sd(CH4_best.flux, na.rm = TRUE) / sqrt(n()),
    n = n(), .groups = 'drop'
  )
cat("\nBetula alleghaniensis fluxes by height:\n")
print(beal_heights)

beal_stats <- species_results %>% filter(species == "BEAL")
cat("\nBetula alleghaniensis height effect:\n")
cat("  Coefficient:", round(beal_stats$height_coef, 6), "\n")
cat("  P-value:", round(beal_stats$height_p, 4), "\n")

# Summary of significant effects
n_sig <- sum(species_results$height_p < 0.05, na.rm = TRUE)
n_tested <- sum(!is.na(species_results$height_coef))
cat("\nSpecies with significant height effects:", n_sig, "out of", n_tested, "tested\n")

sig_species <- species_results %>%
  filter(height_p < 0.05) %>%
  select(species, height_coef, height_p)
cat("\nSpecies with significant effects:\n")
print(sig_species)

# High-flux species without significant height effects
high_flux_species <- plot_data_top %>%
  group_by(species) %>%
  summarise(mean_flux = mean(CH4_best.flux, na.rm = TRUE),
            max_flux = max(CH4_best.flux, na.rm = TRUE), .groups = 'drop') %>%
  filter(mean_flux > 0.05) %>%
  left_join(species_results %>% select(species, height_p), by = "species") %>%
  filter(is.na(height_p) | height_p >= 0.05)
cat("\nSpecies with high fluxes but no significant height effect:\n")
print(high_flux_species)

cat("\n===== FIGURE 2 SUMMARY STATS =====\n")
cat("Total flux measurements in height analysis:", nrow(plot_data_top), "\n")
cat("Total individual trees:", length(unique(plot_data_top$tree_unique)), "\n")
cat("Total species analyzed:", length(unique(plot_data_top$species)), "\n")

# Height coefficients per meter
cat("\n===== HEIGHT COEFFICIENTS IN NMOL M-2 S-1 PER METER =====\n")
species_results_converted <- species_results %>%
  mutate(height_coef_per_m = height_coef * 100, height_se_per_m = height_se * 100) %>%
  filter(height_p < 0.05)
cat("\nSignificant height effects (nmol m-2 s-1 per meter of height):\n")
for(i in 1:nrow(species_results_converted)) {
  cat(species_results_converted$species[i], ": ",
      round(species_results_converted$height_coef_per_m[i], 3),
      " (p = ", round(species_results_converted$height_p[i], 4), ")\n", sep = "")
}


# ==============================================================================
# 9. ROBUSTNESS TESTING — LEAVE-ONE-TREE-OUT
# ==============================================================================

cat("\n===== TESTING ROBUSTNESS OF HEIGHT EFFECTS =====\n")

test_robustness <- function(species_name, data) {
  species_data <- data %>% filter(species == species_name)
  unique_trees <- unique(species_data$tree_unique)

  if(length(unique_trees) < 4) {
    return(data.frame(species = species_name, n_trees = length(unique_trees),
                      robust = FALSE, note = "Too few trees for robustness test"))
  }

  original_model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = species_data)
  original_coef <- summary(original_model)$coefficients["height_numeric", "Estimate"]

  jackknife_results <- data.frame()
  for(tree in unique_trees) {
    subset_data <- species_data %>% filter(tree_unique != tree)
    if(nrow(subset_data) > 6) {
      model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = subset_data)
      coef_val <- summary(model)$coefficients["height_numeric", "Estimate"]
      p_val <- 2 * (1 - pnorm(abs(summary(model)$coefficients["height_numeric", "t value"])))
      jackknife_results <- rbind(jackknife_results,
                                 data.frame(excluded_tree = tree, coef = coef_val, p = p_val))
    }
  }

  if(nrow(jackknife_results) > 0) {
    prop_sig <- sum(jackknife_results$p < 0.05) / nrow(jackknife_results)
    same_sign <- all(sign(jackknife_results$coef) == sign(original_coef))
    coef_range <- range(jackknife_results$coef)

    data.frame(
      species = species_name, n_trees = length(unique_trees),
      original_coef = original_coef,
      original_p = 2 * (1 - pnorm(abs(summary(original_model)$coefficients["height_numeric", "t value"]))),
      prop_still_sig = prop_sig, all_same_sign = same_sign,
      coef_min = coef_range[1], coef_max = coef_range[2],
      robust = prop_sig > 0.75 & same_sign
    )
  } else {
    data.frame(species = species_name, n_trees = length(unique_trees),
               robust = FALSE, note = "Jackknife failed")
  }
}

sig_species_list <- c("BEAL", "TSCA", "ACSA")   # typed: the species named in the text
robustness_results <- do.call(rbind, lapply(sig_species_list, test_robustness, analysis_data))
cat("\nRobustness of significant height effects:\n")
print(robustness_results)


# ==============================================================================
# 10. SOIL CONDITION RELATIONSHIPS
# ==============================================================================

cat("\n===== SOIL CONDITION RELATIONSHIPS =====\n")

soil_by_species <- merged_data %>%
  filter(species_id %in% unique(analysis_data$species)) %>%
  group_by(species_id) %>%
  summarise(
    mean_VWC = mean(VWC_mean, na.rm = TRUE),
    se_VWC = sd(VWC_mean, na.rm = TRUE) / sqrt(n()),
    mean_mcra_mineral = mean(ddpcr_mcra_probe_Mineral_loose, na.rm = TRUE),
    mean_mcra_organic = mean(ddpcr_mcra_probe_Organic_loose, na.rm = TRUE),
    n_trees = n(), .groups = 'drop'
  ) %>%
  rename(species = species_id)

species_env <- species_results %>%
  left_join(soil_by_species, by = "species") %>%
  filter(!is.na(height_coef) & !is.na(mean_VWC))

cat("\nCorrelation between height effect and environmental variables:\n")
vwc_cor <- cor.test(species_env$height_coef, species_env$mean_VWC)
cat("  Height coefficient vs VWC: r =", round(vwc_cor$estimate, 3),
    ", p =", round(vwc_cor$p.value, 4), "\n")

species_env$log_mcra_mineral <- log10(species_env$mean_mcra_mineral + 1)
mcra_cor <- cor.test(species_env$height_coef, species_env$log_mcra_mineral)
cat("  Height coefficient vs log(mcrA+1): r =", round(mcra_cor$estimate, 3),
    ", p =", round(mcra_cor$p.value, 4), "\n")

cat("\nEnvironmental conditions for species with significant height effects:\n")
env_sig <- species_env %>%
  filter(species %in% sig_species_list) %>%
  select(species, height_coef, height_p, mean_VWC, mean_mcra_mineral, mean_mcra_organic) %>%
  arrange(height_coef)
print(env_sig)


# ==============================================================================
# 11. WITHIN-SPECIES VARIATION
# ==============================================================================

cat("\n===== WITHIN-SPECIES VARIATION IN HEIGHT RESPONSE =====\n")

tree_slopes <- analysis_data %>%
  group_by(species, tree_unique) %>%
  filter(n() >= 2) %>%
  summarise(
    n_heights = n(),
    slope = ifelse(length(unique(height_numeric)) >= 2,
                   coef(lm(CH4_best.flux ~ height_numeric))[2], NA),
    .groups = 'drop'
  ) %>%
  filter(!is.na(slope))

species_slope_variation <- tree_slopes %>%
  group_by(species) %>%
  summarise(
    n_trees_with_slopes = n(),
    mean_slope = mean(slope), sd_slope = sd(slope),
    cv_slope = sd(slope) / abs(mean(slope)),
    prop_negative = sum(slope < 0) / n(), .groups = 'drop'
  ) %>%
  filter(n_trees_with_slopes >= 3)

cat("\nVariation in individual tree slopes within species:\n")
print(species_slope_variation %>% arrange(desc(n_trees_with_slopes)))

beal_slopes <- tree_slopes %>% filter(species == "BEAL")
cat("\nBEAL individual tree slopes:\n")
cat("  Total trees with slopes:", nrow(beal_slopes), "\n")
cat("  Trees with negative slopes:", sum(beal_slopes$slope < 0), "\n")
cat("  Proportion negative:", round(sum(beal_slopes$slope < 0) / nrow(beal_slopes), 2), "\n")


# ==============================================================================
# 12. ADDITIONAL CORRELATIONS
# ==============================================================================

cat("\n===== SOIL VWC vs MCRA CORRELATION =====\n")
soil_correlation <- merged_data %>%
  filter(!is.na(VWC_mean) & !is.na(ddpcr_mcra_probe_Mineral_loose)) %>%
  mutate(log_mcra = log10(ddpcr_mcra_probe_Mineral_loose + 1))

vwc_mcra_cor <- cor.test(soil_correlation$VWC_mean, soil_correlation$log_mcra)
cat("VWC vs log(mcrA+1): r =", round(vwc_mcra_cor$estimate, 3),
    ", p =", round(vwc_mcra_cor$p.value, 4), "\n")

cat("\n===== ACSA HEIGHT PATTERN ANALYSIS =====\n")
acsa_data <- analysis_data %>% filter(species == "ACSA")
acsa_200cm <- acsa_data %>% filter(height_numeric == 200) %>% arrange(desc(CH4_best.flux))
cat("ACSA fluxes at 200cm height:\n")
print(acsa_200cm %>% select(tree_unique, CH4_best.flux))

outlier_tree <- acsa_200cm$tree_unique[1]
cat("\nOutlier tree:", outlier_tree, "with flux:", round(acsa_200cm$CH4_best.flux[1], 3), "\n")

acsa_no_outlier <- acsa_data %>% filter(tree_unique != outlier_tree)
if(nrow(acsa_no_outlier) > 10) {
  model_no_outlier <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique),
                           data = acsa_no_outlier)
  coef_no_outlier <- summary(model_no_outlier)$coefficients["height_numeric", "Estimate"]
  p_no_outlier <- 2 * (1 - pnorm(abs(summary(model_no_outlier)$coefficients["height_numeric", "t value"])))

  cat("\nACSA with outlier removed:\n")
  cat("  Height coefficient:", round(coef_no_outlier, 6), "\n")
  cat("  P-value:", round(p_no_outlier, 4), "\n")
}
