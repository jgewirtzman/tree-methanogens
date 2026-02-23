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
#   - methanogen_tree_flux_complete_dataset.csv (from data/processed/flux/)
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#   - tree_id_comprehensive_mapping.csv (from data/processed/tree_data/)
#
# Outputs:
#   - combined_flux_plot.png
# ==============================================================================

library(ggplot2)
library(dplyr)
library(lme4)
library(gridExtra)
library(grid)
library(tidyr)
library(scales)
library(patchwork)
library(gghalves)

# Load the merged dataset
final_dataset <- read.csv('../../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv')
merged_data <- read.csv('../../../data/processed/integrated/merged_tree_dataset_final.csv')

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

# Define publication colors and shapes - red for positive, blue for negative
colors_pub <- c(
  "Significant Negative" = "#4575B4",      # Blue for negative
  "Significant Positive" = "#D73027",      # Red for positive  
  "Non-significant Negative" = "#74ADD1",  # Light blue for negative
  "Non-significant Positive" = "#F46D43"   # Light red for positive
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
  left_join(tree_counts %>% dplyr::select(species, species_label_with_n, species_label_no_n), by = "species") %>%
  mutate(species_label = factor(species_label_with_n, levels = sort(unique(species_label_with_n))))

# Calculate actual data breaks for each species
# Calculate actual data breaks for each species
species_breaks <- plot_data_top %>%
  group_by(species, species_label) %>%
  summarise(
    data_min = min(CH4_best.flux, na.rm = TRUE),
    data_max = max(CH4_best.flux, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    data_mid = (data_min + data_max) / 2,
    # Special handling for Sassafras (SAAL) - only min and max, no middle value
    breaks_list = case_when(
      species == "SAAL" ~ map2(data_min, data_max, ~ c(.x, .y)),  # Only min and max for Sassafras
      TRUE ~ map2(data_min, data_max, ~ c(.x, (.x + .y)/2, .y))   # Min, mid, max for all others
    )
  )

breaks_lookup <- setNames(species_breaks$breaks_list, species_breaks$species_label)

# Analysis for height effects
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
      species = species_name,
      n_obs = n_obs,
      n_trees = n_trees,
      n_trees_multi_height = n_trees_multi_height,
      height_coef = NA,
      height_se = NA,
      height_t = NA,
      height_p = NA,
      model_type = "insufficient_data"
    ))
  }
  
  tryCatch({
    model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = species_data)
    coef_table <- summary(model)$coefficients
    height_row <- which(rownames(coef_table) == "height_numeric")
    
    return(data.frame(
      species = species_name,
      n_obs = n_obs,
      n_trees = n_trees,
      n_trees_multi_height = n_trees_multi_height,
      height_coef = coef_table[height_row, "Estimate"],
      height_se = coef_table[height_row, "Std. Error"],
      height_t = coef_table[height_row, "t value"],
      height_p = 2 * (1 - pnorm(abs(coef_table[height_row, "t value"]))),
      model_type = "mixed_effects"
    ))
  }, error = function(e) {
    tryCatch({
      model <- lm(CH4_best.flux ~ height_numeric, data = species_data)
      coef_table <- summary(model)$coefficients
      height_row <- which(rownames(coef_table) == "height_numeric")
      
      return(data.frame(
        species = species_name,
        n_obs = n_obs,
        n_trees = n_trees,
        n_trees_multi_height = n_trees_multi_height,
        height_coef = coef_table[height_row, "Estimate"],
        height_se = coef_table[height_row, "Std. Error"],
        height_t = coef_table[height_row, "t value"],
        height_p = coef_table[height_row, "Pr(>|t|)"],
        model_type = "linear_model"
      ))
    }, error = function(e2) {
      return(data.frame(
        species = species_name,
        n_obs = n_obs,
        n_trees = n_trees,
        n_trees_multi_height = n_trees_multi_height,
        height_coef = NA,
        height_se = NA,
        height_t = NA,
        height_p = NA,
        model_type = "failed"
      ))
    })
  })
}

species_list <- sort(unique(plot_data_top$species))
species_results <- do.call(rbind, lapply(species_list, test_species_height_effect, analysis_data))

# Middle panel: Height effect coefficients
plot_data_middle <- species_results %>%
  filter(!is.na(height_coef)) %>%
  left_join(tree_counts %>% dplyr::select(species, species_label_with_n, species_label_no_n), by = "species") %>%
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

p_middle <- ggplot(plot_data_middle, aes(x = reorder(species_label, height_coef), y = height_coef*100, color = color_group, shape = color_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = height_coef*100 - 1.96*height_se*100, 
                    ymax = height_coef*100 + 1.96*height_se*100), 
                width = 0.3, alpha = 0.8) +
  scale_color_manual(values = colors_pub, name = "") +
  scale_shape_manual(values = shapes_pub, name = "") +
  labs(
    x = "",
    y = expression("Height effect")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels since they overlap with bottom plot
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(t = 2, r = 5, b = 2, l = 2)  # Reduced margins
  )

# =============================================================================
# Soil conditions and mcrA data for heatmap
# =============================================================================

# Load the mapping file and create lookup (from harmonization script)
mapping <- read.csv("../../../data/processed/tree_data/tree_id_comprehensive_mapping.csv")
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
  
  comprehensive_lookup <- bind_rows(lookup_list) %>%
    distinct() %>%
    mutate(
      original_name_lower = tolower(trimws(original_name)),
      Tree_ID_normalized = tolower(trimws(Tree_ID_normalized))
    )
  return(comprehensive_lookup)
}

tree_lookup <- create_comprehensive_lookup(mapping)

standardize_tree_id <- function(tree_ids) {
  clean_ids <- tolower(trimws(as.character(tree_ids)))
  standardized <- tree_lookup$Tree_ID_normalized[match(clean_ids, tree_lookup$original_name_lower)]
  result <- ifelse(is.na(standardized), clean_ids, standardized)
  return(result)
}

# Get tree IDs that are in the height effect analysis and standardize them
trees_in_analysis_raw <- unique(plot_data_top$tree_id)
trees_in_analysis_standardized <- standardize_tree_id(trees_in_analysis_raw)

# Prepare soil/mcrA data for species in the height effect plot
soil_mcra_data <- merged_data %>%
  filter(tree_id %in% trees_in_analysis_standardized) %>%
  # Map species_id to the species codes used in height analysis
  mutate(
    species = species_id  # Use species_id directly since it's already in the right format
  ) %>%
  # Calculate depth-weighted mcrA with better handling of missing values
  mutate(
    mcra_organic = ddpcr_mcra_probe_Organic_loose,
    mcra_mineral = ddpcr_mcra_probe_Mineral_loose,
    organic_depth = OrganicDepth_mean,
    mineral_depth = MineralDepth_mean,
    mcra_weighted = case_when(
      # Both available and depths available
      !is.na(mcra_organic) & !is.na(mcra_mineral) & !is.na(organic_depth) & !is.na(mineral_depth) & 
        is.finite(mcra_organic) & is.finite(mcra_mineral) & is.finite(organic_depth) & is.finite(mineral_depth) ~ 
        (mcra_organic * organic_depth + mcra_mineral * mineral_depth) / (organic_depth + mineral_depth),
      # Only organic available
      !is.na(mcra_organic) & is.finite(mcra_organic) & (is.na(mcra_mineral) | !is.finite(mcra_mineral)) ~ mcra_organic,
      # Only mineral available  
      !is.na(mcra_mineral) & is.finite(mcra_mineral) & (is.na(mcra_organic) | !is.finite(mcra_organic)) ~ mcra_mineral,
      TRUE ~ NA_real_
    )
  ) %>%
  # Calculate species means
  group_by(species) %>%
  summarise(
    mean_VWC = mean(VWC_mean, na.rm = TRUE),
    mean_mcra_organic = mean(mcra_organic, na.rm = TRUE),
    mean_mcra_mineral = mean(mcra_mineral, na.rm = TRUE),
    mean_mcra_weighted = mean(mcra_weighted, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Add species labels and filter to species in height analysis
  left_join(tree_counts %>% dplyr::select(species, species_label_with_n, species_label_no_n), by = "species") %>%
  filter(!is.na(species_label_with_n))

# Get the species order from the height effect plot
species_order <- plot_data_middle %>% 
  arrange(height_coef) %>% 
  pull(species_label)

# Get the species order without n= for heatmap
species_order_no_n <- plot_data_middle %>% 
  arrange(height_coef) %>% 
  pull(species_label_no_n)

# Prepare heatmap data with z-score normalization using log-transformed mcrA values
heatmap_data <- soil_mcra_data %>%
  # Filter to only species that are in the height effect analysis
  filter(species %in% plot_data_middle$species) %>%
  # Apply log transformation to mcrA values before selecting columns
  mutate(
    # Add 1 to avoid log(0), then log transform
    log_mcra_organic = log10(mean_mcra_organic + 1),
    log_mcra_mineral = log10(mean_mcra_mineral + 1),
    log_mcra_weighted = log10(mean_mcra_weighted + 1)
  ) %>%
  # Select the variables for the heatmap (VWC remains untransformed, mcrA values are now log-transformed)
  dplyr::select(species_label_no_n, mean_VWC, log_mcra_organic, log_mcra_mineral, log_mcra_weighted) %>%
  # Convert to long format
  pivot_longer(cols = c(mean_VWC, log_mcra_organic, log_mcra_mineral, log_mcra_weighted), 
               names_to = "variable", values_to = "value") %>%
  # Calculate z-scores for each variable (now on log-transformed mcrA values)
  group_by(variable) %>%
  mutate(
    z_score = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Clean variable names
  mutate(
    variable = case_when(
      variable == "mean_VWC" ~ "VWC",
      variable == "log_mcra_organic" ~ "log mcrA (Organic)",
      variable == "log_mcra_mineral" ~ "log mcrA (Mineral)", 
      variable == "log_mcra_weighted" ~ "log mcrA (Weighted)",
      TRUE ~ variable
    ),
    variable = factor(variable, levels = c("VWC", "log mcrA (Organic)", "log mcrA (Mineral)", "log mcrA (Weighted)"))
  ) %>%
  # Apply species ordering from height effect plot - only include species present in both datasets
  filter(species_label_no_n %in% species_order_no_n) %>%
  mutate(species_label_no_n = factor(species_label_no_n, levels = species_order_no_n)) %>%
  # Keep infinite and NaN z-scores as missing
  mutate(z_score = case_when(
    is.infinite(z_score) ~ NA_real_,
    is.nan(z_score) ~ NA_real_,
    TRUE ~ z_score
  )) %>%
  # Remove any rows where species_label is NA
  filter(!is.na(species_label_no_n))

# Create a single unified heatmap using z-scores
p_bottom <- ggplot(heatmap_data, aes(x = species_label_no_n, y = variable, fill = z_score)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "#4575B4", 
                       midpoint = 0, na.value = "grey90",
                       name = "Z-score") +
  labs(
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 35, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 8, hjust = 1),
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 2, r = 5, b = 5, l = 50),  # Reduced top margin
    strip.text = element_blank(),
    panel.spacing = unit(0, "lines")
  )

# SOLUTION 1: Use patchwork instead of grid.arrange for much tighter control
library(patchwork)

# Create the plots with minimal margins
create_species_plot_half <- function(species_name, species_data, breaks_data) {
  current_breaks <- breaks_data$breaks_list[[which(breaks_data$species == species_name)]]
  current_label <- breaks_data$species_label[breaks_data$species == species_name]
  
  ggplot(species_data, aes(x = factor(height_m), y = CH4_best.flux)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.7) +
    geom_half_boxplot(alpha = 0.6, outlier.shape = NA, fill = boxplot_color, color = "gray30", side = "r") +
    geom_jitter(alpha = 0.4, size = 0.8, color = "gray20", 
                position = position_jitter(width = 0.15, height = 0)) +
    coord_flip() +
    scale_y_continuous(
      breaks = current_breaks,
      labels = function(x) {
        ifelse(abs(x) < 0.001, "0", 
               ifelse(abs(x) < 0.01, sprintf("%.3f", x),
                      ifelse(abs(x) < 0.1, sprintf("%.2f", x), 
                             sprintf("%.1f", x))))
      }
    ) +
    labs(x = "", y = "") +
    ggtitle(current_label) +
    theme_void() +  # Start with void theme
    theme(
      plot.title = element_text(size = 8, face = "italic", margin = margin(b = 2)),
      axis.text = element_text(size = 8),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      plot.margin = margin(0, 0, 5, 0),  # Zero margins
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
}

# Recreate the plot list
plot_list_half <- list()
for(sp in species_list) {
  sp_data <- plot_data_top %>% filter(species == sp)
  if(nrow(sp_data) > 0) {
    plot_list_half[[sp]] <- create_species_plot_half(sp, sp_data, species_breaks)
  }
}

# SOLUTION 1A: Using patchwork with wrap_plots
p_top_patchwork <- wrap_plots(plot_list_half, ncol = 5) + 
  plot_annotation(
    title = NULL,
    theme = theme(plot.margin = margin(0, 0, 2, 0))
  )


# Choose your preferred solution and replace p_top_final with one of:
p_top_final <- p_top_patchwork
p_top_final


# Updated create_species_plot_half function with all spacing reduced
create_species_plot_half <- function(species_name, species_data, breaks_data) {
  current_breaks <- breaks_data$breaks_list[[which(breaks_data$species == species_name)]]
  current_label <- breaks_data$species_label[breaks_data$species == species_name]
  
  ggplot(species_data, aes(x = factor(height_m), y = CH4_best.flux)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", size=1, alpha = 0.7) +
    geom_half_boxplot(alpha = 0.6, outlier.shape = NA, fill = boxplot_color, color = "gray30", side = "r") +
    geom_jitter(alpha = 0.4, size = 0.8, color = "gray20", 
                position = position_jitter(width = 0.15, height = 0)) +
    coord_flip() +
    scale_y_continuous(
      breaks = current_breaks,
      labels = function(x) {
        ifelse(abs(x) < 0.001, "0", 
               ifelse(abs(x) < 0.01, sprintf("%.3f", x),
                      ifelse(abs(x) < 0.1, sprintf("%.2f", x), 
                             sprintf("%.1f", x))))
      }
    ) +
    labs(x = "", y = "") +
    ggtitle(current_label) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8, face = "italic", margin = margin(b = 1)), # Reduced title margin
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      # Minimal plot margins
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
      # Add border around the plot
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      # Remove panel spacing
      panel.spacing = unit(0, "lines"),
      # Reduce axis margins
      axis.title = element_text(margin = margin(0, 2, 0, 0)),
      axis.text.x = element_text(margin = margin(t = 1, b = 1)),
      axis.text.y = element_text(margin = margin(l = 1, r = 1))
    )
}

# Updated grid arrangement with all spacing controls
library(grid)
library(gridExtra)

p_top_final <- arrangeGrob(
  grobs = plot_list_half, 
  ncol = 5,
  left = textGrob(
    "Height (m)",
    rot = 90, vjust = 1,
    gp = gpar(fontsize = 10)
  ),
  bottom = textGrob(
    expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")"),
    vjust = -0.5,
    gp = gpar(fontsize = 10)
  ),
  # No padding between plots
  padding = unit(1.25, "lines"),
  # Control spacing between plots
  widths = rep(1, 5),  # Equal width for all columns
  heights = rep(1, ceiling(length(plot_list_half)/5)),  # Equal height for all rows
  # Remove margins around the entire grid
  vp = viewport(layout.pos.row = 1, layout.pos.col = 1)
)

# Alternative approach using grid.arrange with additional spacing controls
# p_top_final <- grid.arrange(
#   grobs = plot_list_half, 
#   ncol = 5,
#   left = textGrob(
#     "Height (m)",
#     rot = 90, vjust = 1,
#     gp = gpar(fontsize = 10)
#   ),
#   bottom = textGrob(
#     expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")"),
#     vjust = -0.5,
#     gp = gpar(fontsize = 10)
#   ),
#   padding = unit(0, "lines"),
#   # Additional parameters for tight spacing
#   respect = FALSE,  # Don't respect aspect ratios
#   clip = FALSE      # Allow elements to extend beyond boundaries
# )

# Convert grid.arrange object to a ggplot-compatible object for patchwork
p_top_gg <- wrap_elements(full = p_top_final)

# Top gets 2/3; middle and bottom split the remaining 1/3 in a 1:0.8 ratio
combined_plot <- p_top_gg / p_middle / p_bottom +
  plot_layout(heights = c(5, 1, 1), axes = "collect_x")

# Display the combined plot
print(combined_plot)

# Save the plot as 8x7 inches
ggsave("../../../outputs/figures/combined_flux_plot.png",
       plot = combined_plot,
       width = 8,
       height = 8,
       units = "in",
       dpi = 300)

# Alternative: save as PDF
# ggsave("../../../outputs/figures/combined_flux_plot.pdf",
#        plot = combined_plot,
#        width = 8,
#        height = 7,
#        units = "in")

# Alternative: save as TIFF for publication
# ggsave("../../../outputs/figures/combined_flux_plot.tiff",
#        plot = combined_plot,
#        width = 8,
#        height = 7,
#        units = "in",
#        dpi = 300,
#        compression = "lzw")












# ===== CORRECTED STATISTICS EXTRACTION =====

cat("\n===== HEIGHT EFFECT STATISTICS =====\n")

# Find Betula alleghaniensis specifically
beal_data <- analysis_data %>% 
  filter(species == "BEAL")

# Get flux values at different heights for BEAL
beal_heights <- beal_data %>%
  group_by(height_numeric) %>%
  summarise(
    mean_flux = mean(CH4_best.flux, na.rm = TRUE),
    se_flux = sd(CH4_best.flux, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )
cat("\nBetula alleghaniensis fluxes by height:\n")
print(beal_heights)

# Get the BEAL model statistics
beal_stats <- species_results %>% filter(species == "BEAL")
cat("\nBetula alleghaniensis height effect:\n")
cat("  Coefficient:", round(beal_stats$height_coef, 6), "\n")
cat("  P-value:", round(beal_stats$height_p, 4), "\n")

# Count species with significant effects (using p < 0.05 threshold)
n_sig <- sum(species_results$height_p < 0.05, na.rm = TRUE)
n_tested <- sum(!is.na(species_results$height_coef))
cat("\nSpecies with significant height effects:", n_sig, "out of", n_tested, "tested\n")

# Check which species have significant effects
sig_species <- species_results %>% 
  filter(height_p < 0.05) %>%
  select(species, height_coef, height_p)
cat("\nSpecies with significant effects:\n")
print(sig_species)

# List species with non-significant effects but high mean fluxes
high_flux_species <- plot_data_top %>%
  group_by(species) %>%
  summarise(
    mean_flux = mean(CH4_best.flux, na.rm = TRUE),
    max_flux = max(CH4_best.flux, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(mean_flux > 0.05) %>%  # >0.05 nmol m⁻² s⁻¹
  left_join(species_results %>% select(species, height_p), by = "species") %>%
  filter(is.na(height_p) | height_p >= 0.05)

cat("\nSpecies with high fluxes but no significant height effect:\n")
print(high_flux_species)

# Count total observations
cat("\n===== FIGURE 2 SUMMARY STATS =====\n")
cat("Total flux measurements in height analysis:", nrow(plot_data_top), "\n")
cat("Total individual trees:", length(unique(plot_data_top$tree_unique)), "\n")
cat("Total species analyzed:", length(unique(plot_data_top$species)), "\n")




# Convert height coefficients from per cm to per meter
cat("\n===== HEIGHT COEFFICIENTS IN NMOL M-2 S-1 PER METER =====\n")

species_results_converted <- species_results %>%
  mutate(
    height_coef_per_m = height_coef * 100,  # Convert from per cm to per meter
    height_se_per_m = height_se * 100
  ) %>%
  filter(height_p < 0.05)  # Show only significant ones

cat("\nSignificant height effects (nmol m⁻² s⁻¹ per meter of height):\n")
for(i in 1:nrow(species_results_converted)) {
  cat(species_results_converted$species[i], ": ", 
      round(species_results_converted$height_coef_per_m[i], 3), 
      " (p = ", round(species_results_converted$height_p[i], 4), ")\n", sep="")
}












# ===== ROBUSTNESS TESTING FOR HEIGHT EFFECTS =====

cat("\n===== TESTING ROBUSTNESS OF HEIGHT EFFECTS =====\n")

# For each species with significant effect, test influence of individual trees
test_robustness <- function(species_name, data) {
  species_data <- data %>% filter(species == species_name)
  unique_trees <- unique(species_data$tree_unique)
  
  if(length(unique_trees) < 4) {
    return(data.frame(
      species = species_name,
      n_trees = length(unique_trees),
      robust = FALSE,
      note = "Too few trees for robustness test"
    ))
  }
  
  # Get original model result
  original_model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = species_data)
  original_coef <- summary(original_model)$coefficients["height_numeric", "Estimate"]
  original_p <- 2 * (1 - pnorm(abs(summary(original_model)$coefficients["height_numeric", "t value"])))
  
  # Leave-one-tree-out analysis
  jackknife_results <- data.frame()
  for(tree in unique_trees) {
    subset_data <- species_data %>% filter(tree_unique != tree)
    if(nrow(subset_data) > 6) {  # Need enough data for model
      tryCatch({
        model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = subset_data)
        coef <- summary(model)$coefficients["height_numeric", "Estimate"]
        p_val <- 2 * (1 - pnorm(abs(summary(model)$coefficients["height_numeric", "t value"])))
        jackknife_results <- rbind(jackknife_results, 
                                   data.frame(excluded_tree = tree, coef = coef, p = p_val))
      }, error = function(e) {})
    }
  }
  
  if(nrow(jackknife_results) > 0) {
    # Check if effect remains significant after removing any single tree
    prop_sig <- sum(jackknife_results$p < 0.05) / nrow(jackknife_results)
    same_sign <- all(sign(jackknife_results$coef) == sign(original_coef))
    coef_range <- range(jackknife_results$coef)
    
    return(data.frame(
      species = species_name,
      n_trees = length(unique_trees),
      original_coef = original_coef,
      original_p = original_p,
      prop_still_sig = prop_sig,
      all_same_sign = same_sign,
      coef_min = coef_range[1],
      coef_max = coef_range[2],
      robust = prop_sig > 0.75 & same_sign
    ))
  } else {
    return(data.frame(
      species = species_name,
      n_trees = length(unique_trees),
      robust = FALSE,
      note = "Jackknife failed"
    ))
  }
}

# Test significant species
sig_species_list <- c("BEAL", "TSCA", "ACSA")
robustness_results <- do.call(rbind, lapply(sig_species_list, test_robustness, analysis_data))

cat("\nRobustness of significant height effects:\n")
print(robustness_results)

# ===== RELATIONSHIP WITH SOIL CONDITIONS =====

cat("\n===== SOIL CONDITION RELATIONSHIPS =====\n")

# Get soil data for species with height data
soil_by_species <- merged_data %>%
  filter(species_id %in% unique(analysis_data$species)) %>%
  group_by(species_id) %>%
  summarise(
    mean_VWC = mean(VWC_mean, na.rm = TRUE),
    se_VWC = sd(VWC_mean, na.rm = TRUE) / sqrt(n()),
    mean_mcra_mineral = mean(ddpcr_mcra_probe_Mineral_loose, na.rm = TRUE),
    mean_mcra_organic = mean(ddpcr_mcra_probe_Organic_loose, na.rm = TRUE),
    n_trees = n(),
    .groups = 'drop'
  ) %>%
  rename(species = species_id)

# Merge with height effects
species_env <- species_results %>%
  left_join(soil_by_species, by = "species") %>%
  filter(!is.na(height_coef) & !is.na(mean_VWC))

# Test correlations
cat("\nCorrelation between height effect and environmental variables:\n")

# VWC correlation
vwc_cor <- cor.test(species_env$height_coef, species_env$mean_VWC)
cat("  Height coefficient vs VWC: r =", round(vwc_cor$estimate, 3), 
    ", p =", round(vwc_cor$p.value, 4), "\n")

# Log mcrA correlation (add 1 to avoid log(0))
species_env$log_mcra_mineral <- log10(species_env$mean_mcra_mineral + 1)
mcra_cor <- cor.test(species_env$height_coef, species_env$log_mcra_mineral)
cat("  Height coefficient vs log(mcrA+1): r =", round(mcra_cor$estimate, 3), 
    ", p =", round(mcra_cor$p.value, 4), "\n")

# Show environmental conditions for species with significant effects
cat("\nEnvironmental conditions for species with significant height effects:\n")
env_sig <- species_env %>%
  filter(species %in% sig_species_list) %>%
  select(species, height_coef, height_p, mean_VWC, mean_mcra_mineral, mean_mcra_organic) %>%
  arrange(height_coef)
print(env_sig)

# ===== WITHIN-SPECIES VARIATION ANALYSIS =====

cat("\n===== WITHIN-SPECIES VARIATION IN HEIGHT RESPONSE =====\n")

# For species with enough trees, calculate individual tree slopes
tree_slopes <- analysis_data %>%
  group_by(species, tree_unique) %>%
  filter(n() >= 2) %>%  # Need at least 2 heights
  summarise(
    n_heights = n(),
    slope = ifelse(length(unique(height_numeric)) >= 2,
                   coef(lm(CH4_best.flux ~ height_numeric))[2],
                   NA),
    .groups = 'drop'
  ) %>%
  filter(!is.na(slope))

# Summarize by species
species_slope_variation <- tree_slopes %>%
  group_by(species) %>%
  summarise(
    n_trees_with_slopes = n(),
    mean_slope = mean(slope),
    sd_slope = sd(slope),
    cv_slope = sd(slope) / abs(mean(slope)),
    prop_negative = sum(slope < 0) / n(),
    .groups = 'drop'
  ) %>%
  filter(n_trees_with_slopes >= 3)

cat("\nVariation in individual tree slopes within species:\n")
print(species_slope_variation %>% arrange(desc(n_trees_with_slopes)))

# Check if BEAL trees consistently show negative slopes
beal_slopes <- tree_slopes %>% filter(species == "BEAL")
cat("\nBEAL individual tree slopes:\n")
cat("  Total trees with slopes:", nrow(beal_slopes), "\n")
cat("  Trees with negative slopes:", sum(beal_slopes$slope < 0), "\n")
cat("  Proportion negative:", round(sum(beal_slopes$slope < 0) / nrow(beal_slopes), 2), "\n")








# ===== CHECK SOIL VWC vs MCRA CORRELATION =====
cat("\n===== SOIL VWC vs MCRA CORRELATION =====\n")

soil_correlation <- merged_data %>%
  filter(!is.na(VWC_mean) & !is.na(ddpcr_mcra_probe_Mineral_loose)) %>%
  mutate(log_mcra = log10(ddpcr_mcra_probe_Mineral_loose + 1))

vwc_mcra_cor <- cor.test(soil_correlation$VWC_mean, soil_correlation$log_mcra)
cat("VWC vs log(mcrA+1): r =", round(vwc_mcra_cor$estimate, 3), 
    ", p =", round(vwc_mcra_cor$p.value, 4), "\n")

# ===== CHECK ACSA OUTLIER =====
cat("\n===== ACSA HEIGHT PATTERN ANALYSIS =====\n")

acsa_data <- analysis_data %>% 
  filter(species == "ACSA")

# Find the tree with highest flux at 2m
acsa_200cm <- acsa_data %>% 
  filter(height_numeric == 200) %>%
  arrange(desc(CH4_best.flux))

cat("ACSA fluxes at 200cm height:\n")
print(acsa_200cm %>% select(tree_unique, CH4_best.flux))

# Get the outlier tree
outlier_tree <- acsa_200cm$tree_unique[1]
cat("\nOutlier tree:", outlier_tree, "with flux:", round(acsa_200cm$CH4_best.flux[1], 3), "\n")

# Test ACSA model without the outlier
acsa_no_outlier <- acsa_data %>% 
  filter(tree_unique != outlier_tree)

if(nrow(acsa_no_outlier) > 10) {
  model_no_outlier <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), 
                           data = acsa_no_outlier)
  coef_no_outlier <- summary(model_no_outlier)$coefficients["height_numeric", "Estimate"]
  p_no_outlier <- 2 * (1 - pnorm(abs(summary(model_no_outlier)$coefficients["height_numeric", "t value"])))
  
  cat("\nACSA with outlier removed:\n")
  cat("  Height coefficient:", round(coef_no_outlier, 6), "\n")
  cat("  P-value:", round(p_no_outlier, 4), "\n")
}