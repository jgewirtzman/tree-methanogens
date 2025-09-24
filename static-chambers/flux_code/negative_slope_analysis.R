library(ggplot2)
library(dplyr)
library(lme4)      # for mixed-effects models
library(broom.mixed) # for tidy model output

final_dataset<-read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/methanogen_tree_flux_complete_dataset.csv')

# Prepare the data
analysis_data <- final_dataset %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux) & !is.na(height_numeric)) %>%
  # Create a unique tree identifier
  mutate(tree_unique = paste(plot, tree_id, sep = "_"))

# Check data structure
print("Data summary:")
print(paste("Total observations:", nrow(analysis_data)))
print("Trees with multiple height measurements:")
tree_heights <- analysis_data %>%
  group_by(tree_unique) %>%
  summarise(
    n_heights = n_distinct(height_numeric),
    heights = paste(sort(unique(height_numeric)), collapse = ", "),
    .groups = 'drop'
  ) %>%
  filter(n_heights > 1)
print(paste("Number of trees with multiple heights:", nrow(tree_heights)))
print(head(tree_heights))

# Mixed-effects model: flux ~ height with random intercept for each tree
model_simple <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), 
                     data = analysis_data)

# Mixed-effects model with species as fixed effect
model_species <- lmer(CH4_best.flux ~ height_numeric + species + (1|tree_unique), 
                      data = analysis_data)

# Mixed-effects model with random slopes (if enough data)
if(nrow(tree_heights) > 10) {
  model_slopes <- lmer(CH4_best.flux ~ height_numeric + species + (height_numeric|tree_unique), 
                       data = analysis_data)
  print("Model with random slopes fitted")
} else {
  model_slopes <- NULL
  print("Not enough trees with multiple heights for random slopes model")
}

# Print model summaries
print("=== Simple Model: flux ~ height + (1|tree) ===")
print(summary(model_simple))

print("=== Model with Species: flux ~ height + species + (1|tree) ===")
print(summary(model_species))

# Extract and display key results
# Get coefficients and calculate p-values manually
simple_coef <- summary(model_simple)$coefficients
species_coef <- summary(model_species)$coefficients

print("=== HEIGHT EFFECT RESULTS ===")
print("Simple model:")
height_row_simple <- which(rownames(simple_coef) == "height_numeric")
height_estimate_simple <- simple_coef[height_row_simple, "Estimate"]
height_se_simple <- simple_coef[height_row_simple, "Std. Error"]
height_t_simple <- simple_coef[height_row_simple, "t value"]
# Calculate p-value using normal approximation (common for large samples)
height_p_simple <- 2 * (1 - pnorm(abs(height_t_simple)))

print(paste("Height coefficient:", round(height_estimate_simple, 6)))
print(paste("Standard Error:", round(height_se_simple, 6)))
print(paste("t-value:", round(height_t_simple, 3)))
print(paste("P-value (approx):", round(height_p_simple, 4)))

print("Species model:")
height_row_species <- which(rownames(species_coef) == "height_numeric")
height_estimate_species <- species_coef[height_row_species, "Estimate"]
height_se_species <- species_coef[height_row_species, "Std. Error"]
height_t_species <- species_coef[height_row_species, "t value"]
height_p_species <- 2 * (1 - pnorm(abs(height_t_species)))

print(paste("Height coefficient:", round(height_estimate_species, 6)))
print(paste("Standard Error:", round(height_se_species, 6)))
print(paste("t-value:", round(height_t_species, 3)))
print(paste("P-value (approx):", round(height_p_species, 4)))

# Create visualization showing individual tree trends
viz_data <- analysis_data %>%
  filter(tree_unique %in% tree_heights$tree_unique)

# Plot individual tree trajectories
p1 <- ggplot(viz_data, aes(x = height_numeric, y = CH4_best.flux)) +
  geom_line(aes(group = tree_unique), alpha = 0.3, color = "gray50") +
  geom_point(aes(color = species), size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red", size = 1.2) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    x = "Height (cm)",
    y = "CH₄ Flux",
    title = "CH₄ Flux vs Height: Individual Tree Trajectories",
    subtitle = "Gray lines = individual trees, Red line = overall trend"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p1)

# Alternative: Show within-tree changes
if(nrow(tree_heights) > 5) {
  # Calculate within-tree height effects
  within_tree_effects <- viz_data %>%
    group_by(tree_unique, species) %>%
    filter(n() >= 2) %>%  # Need at least 2 heights
    do({
      if(nrow(.) >= 2 & length(unique(.$height_numeric)) >= 2) {
        mod <- lm(CH4_best.flux ~ height_numeric, data = .)
        data.frame(
          slope = coef(mod)[2],
          p_value = summary(mod)$coefficients[2, 4],
          n_obs = nrow(.)
        )
      } else {
        data.frame(slope = NA, p_value = NA, n_obs = nrow(.))
      }
    }) %>%
    ungroup() %>%
    filter(!is.na(slope))
  
  print("=== WITHIN-TREE ANALYSIS ===")
  print(paste("Trees with calculable slopes:", nrow(within_tree_effects)))
  print("Summary of individual tree slopes:")
  print(summary(within_tree_effects$slope))
  
  negative_slopes <- sum(within_tree_effects$slope < 0, na.rm = TRUE)
  total_slopes <- sum(!is.na(within_tree_effects$slope))
  print(paste("Trees with negative slopes:", negative_slopes, "out of", total_slopes))
  print(paste("Proportion with negative slopes:", round(negative_slopes/total_slopes, 3)))
  
  # Test if slopes are significantly different from zero
  slope_test <- t.test(within_tree_effects$slope)
  print("T-test of individual tree slopes vs 0:")
  print(slope_test)
  
  # Plot distribution of individual tree slopes
  p2 <- ggplot(within_tree_effects, aes(x = slope)) +
    geom_histogram(bins = 15, alpha = 0.7, fill = "steelblue") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = mean(within_tree_effects$slope, na.rm = TRUE), 
               color = "orange", size = 1) +
    labs(
      x = "Within-Tree Slope (CH₄ flux change per cm height)",
      y = "Number of Trees",
      title = "Distribution of Within-Tree Height Effects",
      subtitle = "Red line = 0, Orange line = mean"
    ) +
    theme_minimal()
  
  print(p2)
}

# Summary interpretation
print("=== INTERPRETATION ===")
if(height_p_species < 0.05) {
  if(height_estimate_species < 0) {
    print("SIGNIFICANT DECREASE in CH4 flux with height")
  } else {
    print("SIGNIFICANT INCREASE in CH4 flux with height")
  }
} else {
  print("NO SIGNIFICANT relationship between CH4 flux and height")
}

print(paste("Effect size: For every 1 cm increase in height, CH4 flux changes by", 
            round(height_estimate_species, 6), "units"))

# Calculate effect over the full height range (50 to 200 cm)
height_range_effect <- height_estimate_species * (200 - 50)
print(paste("Effect over full height range (50-200 cm):", round(height_range_effect, 4), "units"))

# Compare models
print("=== MODEL COMPARISON ===")
print(paste("Simple model AIC:", round(AIC(model_simple), 1)))
print(paste("Species model AIC:", round(AIC(model_species), 1)))
print("Lower AIC indicates better model fit")


















library(ggplot2)
library(dplyr)
library(lme4)

# Prepare the data
analysis_data <- final_dataset %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux) & !is.na(height_numeric)) %>%
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

# Get all species
species_list <- unique(analysis_data$species)
print(paste("Testing", length(species_list), "species"))

# Test each species
species_results <- do.call(rbind, lapply(species_list, test_species_height_effect, analysis_data))

# Clean up results and add significance indicators
species_results <- species_results %>%
  mutate(
    significant = ifelse(is.na(height_p), FALSE, height_p < 0.05),
    negative_trend = ifelse(is.na(height_coef), FALSE, height_coef < 0),
    sig_negative = significant & negative_trend,
    height_coef_rounded = round(height_coef, 6),
    height_p_rounded = round(height_p, 4)
  ) %>%
  arrange(height_p)

# Display results
print("=== SPECIES-SPECIFIC HEIGHT EFFECTS ===")
print("Summary table:")
print(species_results %>% 
        select(species, n_obs, n_trees_multi_height, height_coef_rounded, height_p_rounded, significant, negative_trend, sig_negative, model_type))

# Species with significant negative trends
sig_negative_species <- species_results %>%
  filter(sig_negative == TRUE) %>%
  arrange(height_p)

print("\n=== SPECIES WITH SIGNIFICANT NEGATIVE TRENDS ===")
if(nrow(sig_negative_species) > 0) {
  print(sig_negative_species %>% 
          select(species, height_coef_rounded, height_p_rounded, n_trees_multi_height))
} else {
  print("No species show significant negative trends individually")
}

# Species with any significant trends (positive or negative)
sig_any_species <- species_results %>%
  filter(significant == TRUE) %>%
  arrange(height_p)

print("\n=== ALL SPECIES WITH SIGNIFICANT HEIGHT EFFECTS ===")
if(nrow(sig_any_species) > 0) {
  print(sig_any_species %>% 
          select(species, height_coef_rounded, height_p_rounded, negative_trend, n_trees_multi_height))
} else {
  print("No species show significant height effects individually")
}

# Summary statistics
print("\n=== SUMMARY ===")
print(paste("Total species tested:", nrow(species_results)))
print(paste("Species with sufficient data:", sum(species_results$model_type != "insufficient_data")))
print(paste("Species with significant effects:", sum(species_results$significant, na.rm = TRUE)))
print(paste("Species with significant NEGATIVE effects:", sum(species_results$sig_negative, na.rm = TRUE)))
print(paste("Species with negative (but not necessarily significant) trends:", sum(species_results$negative_trend, na.rm = TRUE)))

# Create a summary plot
plot_data <- species_results %>%
  filter(!is.na(height_coef)) %>%
  mutate(
    species = reorder(species, height_coef),
    color_group = case_when(
      sig_negative ~ "Significant Negative",
      significant ~ "Significant Positive", 
      negative_trend ~ "Non-significant Negative",
      TRUE ~ "Non-significant Positive"
    )
  )

p <- ggplot(plot_data, aes(x = species, y = height_coef, color = color_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = height_coef - 1.96*height_se, 
                    ymax = height_coef + 1.96*height_se), 
                width = 0.2) +
  scale_color_manual(values = c(
    "Significant Negative" = "red",
    "Significant Positive" = "blue",
    "Non-significant Negative" = "pink",
    "Non-significant Positive" = "lightblue"
  )) +
  labs(
    x = "Species",
    y = "Height Effect Coefficient",
    title = "Height Effect on CH₄ Flux by Species",
    subtitle = "Error bars show 95% confidence intervals",
    color = "Effect Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)