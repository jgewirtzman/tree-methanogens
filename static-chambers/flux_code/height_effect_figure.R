library(ggplot2)
library(dplyr)
library(lme4)
library(gridExtra)
library(grid)

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

# Define publication colors and shapes - solid and hollow points
colors_pub <- c(
  "Significant Negative" = "#D73027",      # Deep red
  "Significant Positive" = "#4575B4",      # Deep blue  
  "Non-significant Negative" = "#F46D43",  # Light red-orange
  "Non-significant Positive" = "#74ADD1"   # Light blue
)

shapes_pub <- c(
  "Significant Negative" = 16,    # Solid circle
  "Significant Positive" = 16,    # Solid circle
  "Non-significant Negative" = 1, # Hollow circle
  "Non-significant Positive" = 1  # Hollow circle
)

# Color for boxplots - single subtle color
boxplot_color <- "#756BB1"  # Muted purple

# Prepare data for top panel (boxplots)
plot_data_top <- final_dataset %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux) & !is.na(height_numeric)) %>%
  # Fix TSLA -> TSCA typo, remove Kalmia, convert height to meters, and map to Latin names
  mutate(
    species = ifelse(species == "TSLA", "TSCA", species),
    height_m = height_numeric / 100  # Convert cm to meters
  ) %>%
  filter(species != "KALA") %>%  # Remove Kalmia latifolia
  mutate(
    tree_unique = paste(plot, tree_id, sep = "_"),
    species_latin = species_mapping[species]
  )

# Calculate tree counts per species
tree_counts <- plot_data_top %>%
  group_by(species, species_latin) %>%
  summarise(n_trees = n_distinct(tree_unique), .groups = 'drop') %>%
  mutate(species_label = paste0(species_latin, " (n=", n_trees, ")"))

# Add labels to plot data
plot_data_top <- plot_data_top %>%
  left_join(tree_counts %>% select(species, species_label), by = "species") %>%
  mutate(species_label = factor(species_label, levels = sort(unique(species_label))))

# Top panel: Boxplots with linear scale and free scales
p_top <- ggplot(plot_data_top, aes(x = factor(height_m), y = CH4_best.flux)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.7) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, fill = boxplot_color, color = "gray30") +
  geom_jitter(alpha = 0.4, width = 0.2, size = 0.8, color = "gray20") +
  facet_wrap(~ species_label, scales = "free", ncol = 3) +
  coord_flip() +
  scale_y_continuous(
    breaks = function(x) {
      max_val <- max(x, na.rm = TRUE)
      min_val <- min(x, na.rm = TRUE)
      
      # Handle edge cases
      if(is.na(max_val) || is.na(min_val)) return(c(0))
      if(max_val == min_val) return(c(min_val))
      
      # Always use actual data range for better tick placement
      # Create three evenly spaced breaks across the actual data range
      mid_val <- (min_val + max_val) / 2
      result <- c(min_val, mid_val, max_val)
      
      return(result)
    },
    labels = function(x) {
      ifelse(abs(x) < 0.001, "0", 
             ifelse(abs(x) < 0.01, sprintf("%.3f", x),
                    ifelse(abs(x) < 0.1, sprintf("%.2f", x), 
                           sprintf("%.1f", x))))
    }
  ) +
  labs(
    x = "Height (m)",
    y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8, face = "italic"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

# Prepare data for bottom panel (coefficients)
# Re-run the species analysis with corrected names and exclude Kalmia
analysis_data <- plot_data_top %>%
  mutate(tree_unique = paste(plot, tree_id, sep = "_"))

# Function to test height effect for each species
test_species_height_effect <- function(species_name, data) {
  species_data <- data %>% filter(species == species_name)
  
  # Check if we have enough data
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
  
  # Try mixed-effects model first
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
    # Fall back to simple linear model if mixed model fails
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

# Test each species
species_list <- unique(analysis_data$species)
species_results <- do.call(rbind, lapply(species_list, test_species_height_effect, analysis_data))

# Prepare data for bottom panel
plot_data_bottom <- species_results %>%
  filter(!is.na(height_coef)) %>%
  left_join(tree_counts %>% select(species, species_label), by = "species") %>%
  mutate(
    species_label = factor(species_label, levels = sort(unique(species_label))),
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

# Bottom panel: Coefficient plot with solid/hollow points
p_bottom <- ggplot(plot_data_bottom, aes(x = reorder(species_label, height_coef), y = height_coef, color = color_group, shape = color_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = height_coef - 1.96*height_se, 
                    ymax = height_coef + 1.96*height_se), 
                width = 0.3, alpha = 0.8) +
  scale_color_manual(values = colors_pub, name = "") +
  scale_shape_manual(values = shapes_pub, name = "") +
  labs(
    x = "",
    y = expression("Height effect")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8, face = "italic"),
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 50)
  )

# Combine panels with better proportions
combined_plot <- grid.arrange(p_top, p_bottom, ncol = 1, heights = c(2.5, 1))

# Display the combined plot
print(combined_plot)