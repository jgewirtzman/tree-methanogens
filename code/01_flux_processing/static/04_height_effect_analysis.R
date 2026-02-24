# ==============================================================================
# Height Effect Analysis (Figure 2)
# ==============================================================================
# Purpose: Analyzes height-dependent CH4 flux patterns across species with
#   statistical tests and visualization.
#
# Pipeline stage: 03 Analysis
# Run after: 00_harmonization/02_harmonize_all_data.R
#
# Inputs:
#   - goflux_auxfile.csv (from data/processed/flux/) — has measurement_height
#   - CH4_best_flux_lgr_results.csv (from data/processed/flux/) — has CH4_best.flux
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#   - tree_id_comprehensive_mapping.csv (from data/processed/tree_data/)
#
# Outputs:
#   - outputs/figures/main/fig2_height_dependent_flux.png
# ==============================================================================

library(ggplot2)
library(dplyr)
library(purrr)
library(lme4)
library(gridExtra)
library(grid)
library(tidyr)
library(scales)
library(patchwork)
library(gghalves)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

aux <- read.csv('data/processed/flux/goflux_auxfile.csv')
ch4_flux <- read.csv('data/processed/flux/CH4_best_flux_lgr_results.csv')
final_dataset <- merge(ch4_flux, aux[, c("UniqueID", "measurement_height", "tree_id",
                                          "species", "plot")],
                       by = "UniqueID", all.x = TRUE)
names(final_dataset)[names(final_dataset) == "best.flux"] <- "CH4_best.flux"

merged_data <- read.csv('data/processed/integrated/merged_tree_dataset_final.csv')

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

# ==============================================================================
# 2. PREPARE DATA
# ==============================================================================

# Publication colors — red for positive height effect, blue for negative
colors_pub <- c(
  "Significant Negative" = "#4575B4",
  "Significant Positive" = "#D73027",
  "Non-significant Negative" = "#74ADD1",
  "Non-significant Positive" = "#F46D43"
)

shapes_pub <- c(
  "Significant Negative" = 16,
  "Significant Positive" = 16,
  "Non-significant Negative" = 1,
  "Non-significant Positive" = 1
)

boxplot_color <- "#756BB1"

# Prepare data for top panel (boxplots)
plot_data_top <- final_dataset %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux) & !is.na(height_numeric)) %>%
  mutate(
    species = ifelse(species == "TSLA", "TSCA", species),
    height_m = height_numeric / 100
  ) %>%
  filter(species != "KALA") %>%
  mutate(
    tree_unique = paste(plot, tree_id, sep = "_"),
    species_latin = species_mapping[species]
  )

# Calculate tree counts per species
tree_counts <- plot_data_top %>%
  group_by(species, species_latin) %>%
  summarise(n_trees = n_distinct(tree_unique), .groups = 'drop') %>%
  mutate(
    species_label_with_n = paste0(species_latin, "\n(n=", n_trees, ")"),
    species_label_no_n = species_latin
  )

# Add labels to plot data
plot_data_top <- plot_data_top %>%
  left_join(tree_counts %>% dplyr::select(species, species_label_with_n, species_label_no_n),
            by = "species") %>%
  mutate(species_label = factor(species_label_with_n, levels = sort(unique(species_label_with_n))))

# Calculate axis breaks for each species panel
species_breaks <- plot_data_top %>%
  group_by(species, species_label) %>%
  summarise(
    data_min = min(CH4_best.flux, na.rm = TRUE),
    data_max = max(CH4_best.flux, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    breaks_list = case_when(
      species == "SAAL" ~ map2(data_min, data_max, ~ c(.x, .y)),
      TRUE ~ map2(data_min, data_max, ~ c(.x, (.x + .y)/2, .y))
    )
  )

breaks_lookup <- setNames(species_breaks$breaks_list, species_breaks$species_label)

# ==============================================================================
# 3. MIXED-EFFECTS HEIGHT ANALYSIS
# ==============================================================================

analysis_data <- plot_data_top %>%
  mutate(tree_unique = paste(plot, tree_id, sep = "_"))

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

  tryCatch({
    model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = species_data)
    coef_table <- summary(model)$coefficients
    height_row <- which(rownames(coef_table) == "height_numeric")

    data.frame(
      species = species_name, n_obs = n_obs, n_trees = n_trees,
      n_trees_multi_height = n_trees_multi_height,
      height_coef = coef_table[height_row, "Estimate"],
      height_se = coef_table[height_row, "Std. Error"],
      height_t = coef_table[height_row, "t value"],
      height_p = 2 * (1 - pnorm(abs(coef_table[height_row, "t value"]))),
      model_type = "mixed_effects"
    )
  }, error = function(e) {
    tryCatch({
      model <- lm(CH4_best.flux ~ height_numeric, data = species_data)
      coef_table <- summary(model)$coefficients
      height_row <- which(rownames(coef_table) == "height_numeric")

      data.frame(
        species = species_name, n_obs = n_obs, n_trees = n_trees,
        n_trees_multi_height = n_trees_multi_height,
        height_coef = coef_table[height_row, "Estimate"],
        height_se = coef_table[height_row, "Std. Error"],
        height_t = coef_table[height_row, "t value"],
        height_p = coef_table[height_row, "Pr(>|t|)"],
        model_type = "linear_model"
      )
    }, error = function(e2) {
      data.frame(
        species = species_name, n_obs = n_obs, n_trees = n_trees,
        n_trees_multi_height = n_trees_multi_height,
        height_coef = NA, height_se = NA, height_t = NA, height_p = NA,
        model_type = "failed"
      )
    })
  })
}

species_list <- sort(unique(plot_data_top$species))
species_results <- do.call(rbind, lapply(species_list, test_species_height_effect, analysis_data))

# ==============================================================================
# 4. PANEL A — Species boxplots (top)
# ==============================================================================

create_species_plot_half <- function(species_name, species_data, breaks_data) {
  current_breaks <- breaks_data$breaks_list[[which(breaks_data$species == species_name)]]
  current_label <- breaks_data$species_label[breaks_data$species == species_name]

  ggplot(species_data, aes(x = factor(height_m), y = CH4_best.flux)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40",
               linewidth = 0.6, alpha = 0.8) +
    geom_half_boxplot(alpha = 0.6, outlier.shape = NA, fill = boxplot_color,
                      color = "gray30", side = "r") +
    geom_jitter(alpha = 0.4, size = 0.7, color = "gray20",
                position = position_jitter(width = 0.15, height = 0)) +
    coord_flip() +
    scale_y_continuous(
      breaks = current_breaks,
      expand = expansion(mult = c(0.05, 0.12)),
      labels = function(x) {
        ifelse(abs(x) < 0.001, "0",
               ifelse(abs(x) < 0.01, sprintf("%.3f", x),
                      ifelse(abs(x) < 0.1, sprintf("%.2f", x),
                             sprintf("%.1f", x))))
      }
    ) +
    labs(x = "", y = "") +
    ggtitle(current_label) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 7, face = "italic", margin = margin(b = 1),
                                hjust = 0.5),
      axis.text = element_text(size = 7),
      axis.text.x = element_text(margin = margin(t = 1)),
      axis.text.y = element_text(margin = margin(r = 1)),
      panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray40", fill = NA, linewidth = 0.4),
      plot.margin = margin(t = 0, r = 1, b = 0, l = 1)
    )
}

# Build individual species plots
plot_list_half <- list()
for(sp in species_list) {
  sp_data <- plot_data_top %>% filter(species == sp)
  if(nrow(sp_data) > 0) {
    plot_list_half[[sp]] <- create_species_plot_half(sp, sp_data, species_breaks)
  }
}

# Arrange into grid with shared axis labels
# Convert plots to grobs with clipping off so axis text isn't cut
plot_grobs_half <- lapply(plot_list_half, function(p) {
  g <- ggplotGrob(p)
  g$layout$clip <- "off"
  g
})

p_top_final <- arrangeGrob(
  grobs = plot_grobs_half,
  ncol = 5,
  left = textGrob("Height (m)", rot = 90, vjust = 1, gp = gpar(fontsize = 10)),
  bottom = textGrob(
    expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")"),
    vjust = -0.3, gp = gpar(fontsize = 10)
  ),
  padding = unit(0.2, "lines")
)

p_top_gg <- wrap_elements(full = p_top_final)

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
      !significant & negative_trend ~ "Non-significant Negative",
      TRUE ~ "Non-significant Positive"
    )
  ) %>%
  arrange(height_coef)

p_middle <- ggplot(plot_data_middle,
                   aes(x = reorder(species_label, height_coef),
                       y = height_coef * 100,
                       color = color_group, shape = color_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = height_coef * 100 - 1.96 * height_se * 100,
                    ymax = height_coef * 100 + 1.96 * height_se * 100),
                width = 0.3, linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = colors_pub, name = "") +
  scale_shape_manual(values = shapes_pub, name = "") +
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

# Load tree ID mapping for cross-referencing datasets
mapping <- read.csv("data/processed/tree_data/tree_id_comprehensive_mapping.csv")
mapping <- mapping %>% mutate(across(everything(), ~ifelse(.x == "NA", NA, .x)))

create_comprehensive_lookup <- function(mapping) {
  name_columns <- names(mapping)[grepl("^name_in_|^variant_", names(mapping))]
  name_columns <- c(name_columns, "primary_id")

  lookup_list <- list()
  for (col in name_columns) {
    if (col %in% names(mapping)) {
      temp_lookup <- mapping %>%
        filter(!is.na(.data[[col]])) %>%
        dplyr::select(Tree_ID_normalized, original_name = all_of(col)) %>%
        mutate(original_name = as.character(original_name))
      lookup_list[[col]] <- temp_lookup
    }
  }

  bind_rows(lookup_list) %>%
    distinct() %>%
    mutate(
      original_name_lower = tolower(trimws(original_name)),
      Tree_ID_normalized = tolower(trimws(Tree_ID_normalized))
    )
}

tree_lookup <- create_comprehensive_lookup(mapping)

standardize_tree_id <- function(tree_ids) {
  clean_ids <- tolower(trimws(as.character(tree_ids)))
  standardized <- tree_lookup$Tree_ID_normalized[match(clean_ids, tree_lookup$original_name_lower)]
  ifelse(is.na(standardized), clean_ids, standardized)
}

# Get tree IDs in the analysis and standardize them
trees_in_analysis_raw <- unique(plot_data_top$tree_id)
trees_in_analysis_standardized <- standardize_tree_id(trees_in_analysis_raw)

# Calculate species-level soil/mcrA means
soil_mcra_data <- merged_data %>%
  filter(tree_id %in% trees_in_analysis_standardized) %>%
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

# Species order from height effect panel
species_order <- plot_data_middle %>% arrange(height_coef) %>% pull(species_label)
species_order_no_n <- plot_data_middle %>% arrange(height_coef) %>% pull(species_label_no_n)

# Z-score normalize for heatmap, using log-transformed mcrA
heatmap_data <- soil_mcra_data %>%
  filter(species %in% plot_data_middle$species) %>%
  mutate(
    log_mcra_organic = log10(mean_mcra_organic + 1),
    log_mcra_mineral = log10(mean_mcra_mineral + 1),
    log_mcra_weighted = log10(mean_mcra_weighted + 1)
  ) %>%
  dplyr::select(species_label_no_n, mean_VWC, log_mcra_organic,
                log_mcra_mineral, log_mcra_weighted) %>%
  pivot_longer(cols = c(mean_VWC, log_mcra_organic, log_mcra_mineral, log_mcra_weighted),
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
    variable = factor(variable, levels = c("VWC", "log mcrA (Organic)",
                                           "log mcrA (Mineral)", "log mcrA (Weighted)"))
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
  plot_layout(heights = c(5, 1.2, 1.2)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")",
                  theme = theme(plot.tag = element_text(size = 11, face = "bold")))

print(combined_plot)

ggsave("outputs/figures/main/fig2_height_dependent_flux.png",
       plot = combined_plot,
       width = 7, height = 7.5,
       units = "in", dpi = 300)


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
      tryCatch({
        model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = subset_data)
        coef_val <- summary(model)$coefficients["height_numeric", "Estimate"]
        p_val <- 2 * (1 - pnorm(abs(summary(model)$coefficients["height_numeric", "t value"])))
        jackknife_results <- rbind(jackknife_results,
                                   data.frame(excluded_tree = tree, coef = coef_val, p = p_val))
      }, error = function(e) {})
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

sig_species_list <- c("BEAL", "TSCA", "ACSA")
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
