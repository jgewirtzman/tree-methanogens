# ==============================================================================
# Chamber Comparison Utility
# ==============================================================================
# Purpose: Compares rigid vs semirigid chamber measurements for quality control.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)
library(scales)

# Load the datasets
methanogen_data <- read.csv('../../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv', stringsAsFactors = FALSE)
semirigid_data <- read.csv('../../../data/processed/flux/semirigid_tree_final_complete_dataset.csv', stringsAsFactors = FALSE)
species_mapping <- read.csv('../../../data/raw/lgr/semirigid_2020-2021/Yale Myers Methane Project/spatial_data/YM_trees_measured.csv.csv', stringsAsFactors = FALSE)

# Data preparation for methanogen dataset
methanogen_clean <- methanogen_data %>%
  # Remove rows with missing species or flux data
  filter(!is.na(species) & species != "" & species != "NA") %>%
  # Convert flux columns to numeric (handle "NA" strings)
  mutate(
    CO2_best.flux = as.numeric(ifelse(CO2_best.flux == "NA", NA, CO2_best.flux)),
    CH4_best.flux = as.numeric(ifelse(CH4_best.flux == "NA", NA, CH4_best.flux))
  ) %>%
  # Remove rows with missing flux data
  filter(!is.na(CO2_best.flux) | !is.na(CH4_best.flux)) %>%
  # Add method identifier
  mutate(method = "Methanogen Chamber")

# Data preparation for semirigid dataset
semirigid_clean <- semirigid_data %>%
  # Convert Plot.Tag to character to match Label type
  mutate(Plot.Tag = as.character(Plot.Tag)) %>%
  # Map Plot Tag to actual species using the species mapping file
  left_join(species_mapping, by = c("Plot.Tag" = "Label")) %>%
  # Filter out rows where species mapping failed
  filter(!is.na(Species)) %>%
  # Rename the Species column to species for consistency
  rename(species = Species) %>%
  # Use the .x version of CH4 flux (assuming this is the primary measurement)
  mutate(
    CH4_best.flux = CH4_best.flux.x,
    # Add site information for additional grouping if needed
    site_species = paste0(species, "_", Site)
  ) %>%
  # Remove rows with missing flux data
  filter(!is.na(CO2_best.flux) | !is.na(CH4_best.flux)) %>%
  # Add method identifier
  mutate(method = "Semirigid Chamber") %>%
  # dplyr::select relevant columns to match methanogen dataset structure
  dplyr::select(UniqueID, species, site_species, Site, CO2_best.flux, CH4_best.flux, method)

# Combine datasets for comparison
combined_data <- bind_rows(
  methanogen_clean %>% dplyr::select(UniqueID, species, CO2_best.flux, CH4_best.flux, method),
  semirigid_clean %>% dplyr::select(UniqueID, species, CO2_best.flux, CH4_best.flux, method)
)

# Create long format for easier plotting
combined_long <- combined_data %>%
  pivot_longer(
    cols = c(CO2_best.flux, CH4_best.flux),
    names_to = "gas_type",
    values_to = "flux_value"
  ) %>%
  filter(!is.na(flux_value)) %>%
  mutate(
    gas_type = case_when(
      gas_type == "CO2_best.flux" ~ "CO2",
      gas_type == "CH4_best.flux" ~ "CH4"
    )
  )

# Function to create species distribution plots
create_flux_distribution_plot <- function(data, gas, method_name, top_n = 10) {
  # Filter data for specific gas and method
  plot_data <- data %>%
    filter(gas_type == gas & method == method_name)
  
  # Get top species by sample count
  top_species <- plot_data %>%
    count(species) %>%
    top_n(top_n, n) %>%
    pull(species)
  
  # Filter to top species
  plot_data <- plot_data %>%
    filter(species %in% top_species) %>%
    mutate(species = factor(species, levels = names(sort(table(species), decreasing = TRUE))))
  
  # Create the plot
  ggplot(plot_data, aes(x = species, y = flux_value, fill = species)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(alpha = 0.4, width = 0.2, size = 0.8) +
    scale_fill_viridis_d(name = "Species/Plot") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    labs(
      title = paste(gas, "Flux Distribution -", method_name),
      subtitle = paste("Top", min(top_n, length(top_species)), "species/plots by sample count"),
      x = "Species/Plot ID",
      y = paste(gas, "Flux (μmol/m²/s)")
    )
}

# Create individual plots
co2_methanogen <- create_flux_distribution_plot(combined_long, "CO2", "Methanogen Chamber")
ch4_methanogen <- create_flux_distribution_plot(combined_long, "CH4", "Methanogen Chamber")
co2_semirigid <- create_flux_distribution_plot(combined_long, "CO2", "Semirigid Chamber")
ch4_semirigid <- create_flux_distribution_plot(combined_long, "CH4", "Semirigid Chamber")

# Display plots in a grid
grid.arrange(co2_methanogen, ch4_methanogen, co2_semirigid, ch4_semirigid, ncol = 2)

# Summary statistics by method and gas type
summary_stats <- combined_long %>%
  group_by(method, gas_type) %>%
  summarise(
    n_observations = n(),
    n_species = n_distinct(species),
    mean_flux = mean(flux_value, na.rm = TRUE),
    median_flux = median(flux_value, na.rm = TRUE),
    sd_flux = sd(flux_value, na.rm = TRUE),
    min_flux = min(flux_value, na.rm = TRUE),
    max_flux = max(flux_value, na.rm = TRUE),
    .groups = 'drop'
  )

print("Summary Statistics by Method and Gas Type:")
print(summary_stats)

# Method comparison plots
comparison_plot <- function(gas) {
  plot_data <- combined_long %>%
    filter(gas_type == gas)
  
  ggplot(plot_data, aes(x = method, y = flux_value, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_jitter(alpha = 0.3, width = 0.2, size = 0.8) +
    scale_fill_manual(values = c("Methanogen Chamber" = "#440154", "Semirigid Chamber" = "#31688e")) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = paste(gas, "Flux: Method Comparison"),
      x = "Measurement Method",
      y = paste(gas, "Flux (μmol/m²/s)")
    )
}

# Create method comparison plots
co2_comparison <- comparison_plot("CO2")
ch4_comparison <- comparison_plot("CH4")

# Display comparison plots
grid.arrange(co2_comparison, ch4_comparison, ncol = 2)

# Statistical tests to compare methods
cat("\n=== Statistical Comparisons Between Methods ===\n")

# CO2 comparison
co2_methanogen_vals <- combined_long %>% 
  filter(gas_type == "CO2" & method == "Methanogen Chamber") %>% 
  pull(flux_value)
co2_semirigid_vals <- combined_long %>% 
  filter(gas_type == "CO2" & method == "Semirigid Chamber") %>% 
  pull(flux_value)

if(length(co2_methanogen_vals) > 0 & length(co2_semirigid_vals) > 0) {
  cat("\nCO2 Flux Comparison (Wilcoxon Rank Sum Test):\n")
  co2_test <- wilcox.test(co2_methanogen_vals, co2_semirigid_vals)
  print(co2_test)
}

# CH4 comparison
ch4_methanogen_vals <- combined_long %>% 
  filter(gas_type == "CH4" & method == "Methanogen Chamber") %>% 
  pull(flux_value)
ch4_semirigid_vals <- combined_long %>% 
  filter(gas_type == "CH4" & method == "Semirigid Chamber") %>% 
  pull(flux_value)

if(length(ch4_methanogen_vals) > 0 & length(ch4_semirigid_vals) > 0) {
  cat("\nCH4 Flux Comparison (Wilcoxon Rank Sum Test):\n")
  ch4_test <- wilcox.test(ch4_methanogen_vals, ch4_semirigid_vals)
  print(ch4_test)
}

# Species-level analysis for methanogen data
cat("\n=== Top Species in Methanogen Dataset ===\n")
methanogen_species_summary <- methanogen_clean %>%
  group_by(species) %>%
  summarise(
    n_observations = n(),
    co2_mean = mean(CO2_best.flux, na.rm = TRUE),
    co2_median = median(CO2_best.flux, na.rm = TRUE),
    ch4_mean = mean(CH4_best.flux, na.rm = TRUE),
    ch4_median = median(CH4_best.flux, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_observations))

print(head(methanogen_species_summary, 10))

# Species comparison between methods
cat("\n=== Species Overlap Between Methods ===\n")
methanogen_species <- unique(methanogen_clean$species)
semirigid_species <- unique(semirigid_clean$species)
common_species <- intersect(methanogen_species, semirigid_species)

cat("Species in methanogen dataset:", length(methanogen_species), "\n")
cat("Species in semirigid dataset:", length(semirigid_species), "\n")
cat("Common species:", length(common_species), "\n")
cat("Common species list:", paste(common_species, collapse = ", "), "\n")

# Analysis by site for semirigid data
cat("\n=== Semirigid Data by Site ===\n")
semirigid_site_summary <- semirigid_clean %>%
  group_by(Site, species) %>%
  summarise(
    n_observations = n(),
    co2_mean = mean(CO2_best.flux, na.rm = TRUE),
    ch4_mean = mean(CH4_best.flux, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Site, desc(n_observations))

print(semirigid_site_summary)

# Create site-specific plots for semirigid data
if(nrow(semirigid_clean) > 0) {
  semirigid_site_plot <- semirigid_clean %>%
    pivot_longer(
      cols = c(CO2_best.flux, CH4_best.flux),
      names_to = "gas_type",
      values_to = "flux_value"
    ) %>%
    filter(!is.na(flux_value)) %>%
    mutate(
      gas_type = case_when(
        gas_type == "CO2_best.flux" ~ "CO2",
        gas_type == "CH4_best.flux" ~ "CH4"
      )
    )
  
  site_flux_plot <- ggplot(semirigid_site_plot, aes(x = species, y = flux_value, fill = Site)) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
    facet_wrap(~gas_type, scales = "free_y") +
    scale_fill_viridis_d(name = "Site") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      title = "Semirigid Chamber: Flux by Species and Site",
      x = "Species",
      y = "Flux (μmol/m²/s)"
    )
  
  print(site_flux_plot)
}

# Create a correlation plot if both gases have data
correlation_data <- combined_data %>%
  filter(!is.na(CO2_best.flux) & !is.na(CH4_best.flux))

if(nrow(correlation_data) > 10) {
  correlation_plot <- ggplot(correlation_data, aes(x = CO2_best.flux, y = CH4_best.flux, color = method)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.3) +
    scale_color_manual(values = c("Methanogen Chamber" = "#440154", "Semirigid Chamber" = "#31688e")) +
    theme_minimal() +
    labs(
      title = "CO2 vs CH4 Flux Correlation",
      x = "CO2 Flux (μmol/m²/s)",
      y = "CH4 Flux (μmol/m²/s)",
      color = "Method"
    ) +
    theme(legend.position = "bottom")
  
  print(correlation_plot)
}

# Direct species comparison for common species
if(length(common_species) > 0) {
  cat("\n=== Direct Species Comparison Between Methods ===\n")
  
  common_species_data <- combined_long %>%
    filter(species %in% common_species)
  
  if(nrow(common_species_data) > 0) {
    species_comparison_plot <- ggplot(common_species_data, 
                                      aes(x = species, y = flux_value, fill = method)) +
      geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
      facet_wrap(~gas_type, scales = "free_y") +
      scale_fill_manual(values = c("Methanogen Chamber" = "#440154", "Semirigid Chamber" = "#31688e")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      ) +
      labs(
        title = "Method Comparison for Common Species",
        x = "Species",
        y = "Flux (μmol/m²/s)",
        fill = "Method"
      )
    
    print(species_comparison_plot)
    
    # Statistical comparison for each common species
    for(sp in common_species) {
      for(gas in c("CO2", "CH4")) {
        meth_data <- common_species_data %>%
          filter(species == sp & gas_type == gas & method == "Methanogen Chamber") %>%
          pull(flux_value)
        semi_data <- common_species_data %>%
          filter(species == sp & gas_type == gas & method == "Semirigid Chamber") %>%
          pull(flux_value)
        
        if(length(meth_data) >= 3 & length(semi_data) >= 3) {
          test_result <- wilcox.test(meth_data, semi_data)
          cat(sprintf("%s %s flux - %s: p-value = %.4f\n", 
                      sp, gas, 
                      ifelse(test_result$p.value < 0.05, "Significant", "Not significant"),
                      test_result$p.value))
        }
      }
    }
  }
}

cat("\nAnalysis complete! Check the plots and statistical summaries above.\n")




# Statistical Analysis: Chamber Effect Controlling for Confounding Variables
# Testing whether chamber type matters when controlling for species, height, season, and site

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(broom)
library(broom.mixed)
library(lubridate)
library(ggeffects)

# Load the datasets
methanogen_data <- read.csv('../../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv', stringsAsFactors = FALSE)
semirigid_data <- read.csv('../../../data/processed/flux/semirigid_tree_final_complete_dataset.csv', stringsAsFactors = FALSE)
species_mapping <- read.csv('../../../data/raw/lgr/semirigid_2020-2021/Yale Myers Methane Project/spatial_data/YM_trees_measured.csv.csv', stringsAsFactors = FALSE)

# Prepare methanogen data with detailed covariates
methanogen_prep <- methanogen_data %>%
  filter(!is.na(species) & species != "" & species != "NA") %>%
  mutate(
    CO2_best.flux = as.numeric(ifelse(CO2_best.flux == "NA", NA, CO2_best.flux)),
    CH4_best.flux = as.numeric(ifelse(CH4_best.flux == "NA", NA, CH4_best.flux)),
    date = as.Date(start.time),
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring", 
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall"
    ),
    height_cm = as.numeric(measurement_height),
    height_category = case_when(
      height_cm < 100 ~ "Low",
      height_cm >= 100 & height_cm < 150 ~ "Mid", 
      height_cm >= 150 ~ "High"
    ),
    site_type = "Upland",  # All methanogen measurements are upland
    chamber_type = "Methanogen",
    method = "Methanogen Chamber"
  ) %>%
  filter(!is.na(CO2_best.flux) | !is.na(CH4_best.flux))

# Prepare semirigid data with mapping
semirigid_prep <- semirigid_data %>%
  mutate(Plot.Tag = as.character(Plot.Tag)) %>%
  left_join(species_mapping, by = c("Plot.Tag" = "Label")) %>%
  filter(!is.na(Species)) %>%
  rename(species = Species) %>%
  mutate(
    date = as.Date(Date),
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer", 
      month %in% c(9, 10, 11) ~ "Fall"
    ),
    CO2_best.flux = as.numeric(CO2_best.flux),
    CH4_best.flux = as.numeric(CH4_best.flux.x),
    height_cm = NA,  # Semirigid used only middle height
    height_category = "Mid",  # All semirigid at middle height
    site_type = Site,
    chamber_type = "Semirigid",
    method = "Semirigid Chamber"
  ) %>%
  filter(!is.na(CO2_best.flux) | !is.na(CH4_best.flux)) %>%
  dplyr::select(UniqueID, species, site_type, season, month, height_category, height_cm, 
         chamber_type, method, CO2_best.flux, CH4_best.flux, date)

# Combine datasets
combined_analysis <- bind_rows(
  methanogen_prep %>% dplyr::select(UniqueID, species, site_type, season, month, height_category, 
                             height_cm, chamber_type, method, CO2_best.flux, CH4_best.flux, date),
  semirigid_prep
)

# Create long format for mixed models
combined_long <- combined_analysis %>%
  pivot_longer(
    cols = c(CO2_best.flux, CH4_best.flux),
    names_to = "gas_type", 
    values_to = "flux_value"
  ) %>%
  filter(!is.na(flux_value)) %>%
  mutate(
    gas_type = gsub("_best.flux", "", gas_type),
    # Log transform flux for better model fit (add small constant for any zeros)
    log_flux = log(flux_value + 0.001),
    # Center continuous predictors
    month_centered = scale(month, center = TRUE, scale = FALSE)[,1]
  )

cat("=== DATA OVERVIEW ===\n")
cat("Total observations:", nrow(combined_long), "\n")

# Summary by chamber type and confounders
summary_table <- combined_long %>%
  group_by(chamber_type, gas_type, species, season, site_type, height_category) %>%
  summarise(
    n_obs = n(),
    mean_flux = mean(flux_value, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_obs >= 3)

cat("Data distribution across factors:\n")
print(table(combined_long$chamber_type, combined_long$season))
cat("\n")
print(table(combined_long$chamber_type, combined_long$site_type))
cat("\n")
print(table(combined_long$chamber_type, combined_long$height_category))

# ===== APPROACH 1: MATCHED SUBSET ANALYSIS =====
cat("\n\n=== APPROACH 1: MATCHED SUBSET ANALYSIS ===\n")
cat("Analyzing only summer, upland, mid-height measurements to control confounders\n")

matched_data <- combined_long %>%
  filter(
    season == "Summer",
    site_type == "Upland", 
    height_category == "Mid"
  )

cat("Matched subset size:", nrow(matched_data), "observations\n")
cat("Species in matched subset:", length(unique(matched_data$species)), "\n")

if(nrow(matched_data) > 20) {
  
  # Run mixed effects models for matched data
  for(gas in c("CO2", "CH4")) {
    cat(sprintf("\n--- %s Analysis (Matched Subset) ---\n", gas))
    
    gas_data <- matched_data %>% filter(gas_type == gas)
    
    if(nrow(gas_data) > 10 && length(unique(gas_data$species)) > 2) {
      
      # Model with species as random effect
      model_matched <- lmer(log_flux ~ chamber_type + (1|species), 
                            data = gas_data, 
                            REML = FALSE)
      
      cat("Mixed Effects Model Results:\n")
      print(summary(model_matched))
      
      # Test chamber effect
      model_null <- lmer(log_flux ~ 1 + (1|species), data = gas_data, REML = FALSE)
      lr_test <- anova(model_null, model_matched)
      cat("\nLikelihood Ratio Test for Chamber Effect:\n")
      print(lr_test)
      
      # Effect sizes
      emmeans_result <- emmeans(model_matched, "chamber_type")
      cat("\nEstimated Marginal Means:\n")
      print(emmeans_result)
      
      contrast_result <- contrast(emmeans_result, "pairwise")
      cat("\nPairwise Contrast:\n")
      print(contrast_result)
    }
  }
}

# ===== APPROACH 2: FULL MODEL WITH INTERACTIONS =====
cat("\n\n=== APPROACH 2: FULL MODEL WITH ALL DATA ===\n")
cat("Including all data with chamber × covariate interactions\n")

# Check which species appear in both chamber types
species_both <- combined_long %>%
  group_by(species) %>%
  summarise(n_chambers = n_distinct(chamber_type)) %>%
  filter(n_chambers == 2) %>%
  pull(species)

cat("Species measured by both chamber types:", length(species_both), "\n")
cat("Species list:", paste(species_both, collapse = ", "), "\n")

# Full model analysis
for(gas in c("CO2", "CH4")) {
  cat(sprintf("\n--- %s Analysis (Full Model) ---\n", gas))
  
  gas_data <- combined_long %>% 
    filter(gas_type == gas) %>%
    # Only include species measured by both methods for fair comparison
    filter(species %in% species_both)
  
  if(nrow(gas_data) > 30 && length(unique(gas_data$species)) > 2) {
    
    cat("Sample size:", nrow(gas_data), "\n")
    cat("Species:", length(unique(gas_data$species)), "\n")
    
    # Full model with interactions
    tryCatch({
      model_full <- lmer(log_flux ~ chamber_type * season + 
                           chamber_type * site_type +
                           chamber_type * height_category +
                           month_centered +
                           (1|species), 
                         data = gas_data,
                         REML = FALSE)
      
      cat("Full Model Summary:\n")
      print(summary(model_full))
      
      # Type II tests for main effects and interactions
      cat("\nType II Tests:\n")
      print(Anova(model_full, type = "II"))
      
    }, error = function(e) {
      cat("Full model failed, trying simplified model...\n")
      
      # Simplified model
      model_simple <- lmer(log_flux ~ chamber_type + season + site_type + 
                             height_category + month_centered + (1|species),
                           data = gas_data, 
                           REML = FALSE)
      
      cat("Simplified Model Summary:\n")
      print(summary(model_simple))
      
      cat("\nType II Tests:\n")
      print(Anova(model_simple, type = "II"))
      
      # Test chamber effect specifically
      model_no_chamber <- lmer(log_flux ~ season + site_type + height_category + 
                                 month_centered + (1|species),
                               data = gas_data, REML = FALSE)
      
      lr_test <- anova(model_no_chamber, model_simple)
      cat("\nLikelihood Ratio Test for Chamber Effect:\n")
      print(lr_test)
    })
  }
}

# ===== APPROACH 3: SPECIES-SPECIFIC ANALYSIS =====
cat("\n\n=== APPROACH 3: SPECIES-SPECIFIC ANALYSIS ===\n")
cat("Testing chamber effects within each species that has both chamber types\n")

species_results <- data.frame()

for(sp in species_both) {
  for(gas in c("CO2", "CH4")) {
    species_data <- combined_long %>%
      filter(species == sp & gas_type == gas) %>%
      filter(n_distinct(chamber_type) == 2)  # Must have both chamber types
    
    if(nrow(species_data) >= 6) {  # Minimum sample size
      
      # Simple t-test or Mann-Whitney U test
      meth_vals <- species_data %>% filter(chamber_type == "Methanogen") %>% pull(flux_value)
      semi_vals <- species_data %>% filter(chamber_type == "Semirigid") %>% pull(flux_value)
      
      # Use Mann-Whitney U test (non-parametric)
      test_result <- wilcox.test(meth_vals, semi_vals)
      
      # Calculate effect size (median difference)
      median_diff <- median(meth_vals, na.rm = TRUE) - median(semi_vals, na.rm = TRUE)
      
      species_results <- rbind(species_results, data.frame(
        species = sp,
        gas_type = gas,
        n_methanogen = length(meth_vals),
        n_semirigid = length(semi_vals), 
        p_value = test_result$p.value,
        median_methanogen = median(meth_vals, na.rm = TRUE),
        median_semirigid = median(semi_vals, na.rm = TRUE),
        median_difference = median_diff,
        method = "Mann-Whitney U"
      ))
    }
  }
}

# Multiple testing correction
if(nrow(species_results) > 0) {
  species_results$p_adjusted <- p.adjust(species_results$p_value, method = "fdr")
  species_results$significant <- species_results$p_adjusted < 0.05
  
  cat("Species-specific chamber comparisons:\n")
  print(species_results %>% arrange(p_value))
  
  cat("\nNumber of significant comparisons (FDR < 0.05):", sum(species_results$significant), "\n")
}

# ===== VISUALIZATION OF CONFOUNDING =====
cat("\n\n=== CONFOUNDING VISUALIZATION ===\n")

# Plot showing the confounding structure
confound_plot <- combined_long %>%
  group_by(chamber_type, season, site_type, height_category, gas_type) %>%
  summarise(
    mean_flux = mean(flux_value, na.rm = TRUE),
    n_obs = n(),
    .groups = 'drop'
  ) %>%
  filter(n_obs >= 3) %>%
  ggplot(aes(x = paste(season, site_type, height_category, sep = "-"), 
             y = mean_flux, 
             fill = chamber_type)) +
  geom_col(position = "dodge", alpha = 0.7) +
  facet_wrap(~gas_type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Mean Flux by Chamber Type and Experimental Conditions",
    subtitle = "Shows confounding between chamber type and season/site/height",
    x = "Season - Site - Height",
    y = "Mean Flux",
    fill = "Chamber Type"
  )

print(confound_plot)

# Summary of analysis approaches
cat("\n\n=== SUMMARY ===\n")
cat("Three approaches were used to test chamber effects:\n")
cat("1. MATCHED SUBSET: Compare chambers using only summer/upland/mid-height data\n")
cat("2. FULL MODEL: Include all data with chamber × covariate interactions\n") 
cat("3. SPECIES-SPECIFIC: Test chamber effects within each species separately\n")
cat("\nConfounding factors:\n")
cat("- Season: Methanogen = summer only, Semirigid = all seasons\n")
cat("- Site: Methanogen = upland only, Semirigid = gradient\n") 
cat("- Height: Methanogen = 3 heights, Semirigid = 1 height (mid)\n")
cat("\nRecommendation: Focus on matched subset and species-specific analyses\n")
cat("for the most reliable chamber comparison.\n")




























# Matched Subset Overlap Analysis
# Only summer/upland/mid-height measurements to control confounders

library(dplyr)

# First, need to recreate the matched subset from the ridge plot analysis
# This filters to summer, upland, mid-height to control for confounding

# Create matched subset data using the same approach as statistical analysis
# Use the combined_analysis dataset that already has all covariates

# Recreate the combined_analysis dataset with all covariates
# (This mirrors the statistical analysis code exactly)

# Prepare rigid data with detailed covariates  
rigid_prep <- rigid_data %>%
  filter(!is.na(species) & species != "" & species != "NA") %>%
  mutate(
    CO2_best.flux = as.numeric(ifelse(CO2_best.flux == "NA", NA, CO2_best.flux)),
    CH4_best.flux = as.numeric(ifelse(CH4_best.flux == "NA", NA, CH4_best.flux)),
    date = as.Date(start.time),
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring", 
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall"
    ),
    height_cm = as.numeric(measurement_height),
    height_category = case_when(
      height_cm < 100 ~ "Low",
      height_cm >= 100 & height_cm < 150 ~ "Mid", 
      height_cm >= 150 ~ "High"
    ),
    site_type = "Upland",  # All rigid measurements are upland
    chamber_type = "Rigid",
    method = "Rigid Chamber"
  ) %>%
  filter(!is.na(CO2_best.flux) | !is.na(CH4_best.flux))

# Prepare semirigid data with mapping
semirigid_prep <- semirigid_data %>%
  mutate(Plot.Tag = as.character(Plot.Tag)) %>%
  left_join(species_mapping, by = c("Plot.Tag" = "Label")) %>%
  filter(!is.na(Species)) %>%
  rename(species = Species) %>%
  mutate(
    date = as.Date(Date),
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer", 
      month %in% c(9, 10, 11) ~ "Fall"
    ),
    CO2_best.flux = as.numeric(CO2_best.flux),
    CH4_best.flux = as.numeric(CH4_best.flux.x),
    height_cm = NA,  # Semirigid used only middle height
    height_category = "Mid",  # All semirigid at middle height
    site_type = Site,
    chamber_type = "Semirigid",
    method = "Semirigid Chamber"
  ) %>%
  filter(!is.na(CO2_best.flux) | !is.na(CH4_best.flux)) %>%
  dplyr::select(UniqueID, species, site_type, season, month, height_category, height_cm, 
         chamber_type, method, CO2_best.flux, CH4_best.flux, date)

# Combine datasets (exactly like statistical analysis)
combined_analysis_fixed <- bind_rows(
  rigid_prep %>% dplyr::select(UniqueID, species, site_type, season, month, height_category, 
                        height_cm, chamber_type, method, CO2_best.flux, CH4_best.flux, date),
  semirigid_prep
)

# Create long format
combined_long_fixed <- combined_analysis_fixed %>%
  pivot_longer(
    cols = c(CO2_best.flux, CH4_best.flux),
    names_to = "gas_type", 
    values_to = "flux_value"
  ) %>%
  filter(!is.na(flux_value)) %>%
  mutate(
    gas_type = gsub("_best.flux", "", gas_type)
  )

# NOW create matched subset: summer + upland + mid-height only
matched_subset <- combined_long_fixed %>%
  filter(
    season == "Summer",
    site_type == "Upland", 
    height_category == "Mid"
  )

cat("=== MATCHED SUBSET ANALYSIS ===\n")
cat("Filtered to: Summer + Upland + Mid-height measurements only\n")
cat("Sample size:", nrow(matched_subset), "observations\n\n")

# Calculate overlap statistics for matched subset
matched_overlap_stats <- matched_subset %>%
  group_by(species, gas_type) %>%
  # Only include species-gas combinations with both chamber types
  filter(n_distinct(method) == 2 & n() >= 6) %>%  # Minimum 6 total observations
  summarise(
    # Sample sizes
    n_rigid = sum(method == "Rigid Chamber"),
    n_semirigid = sum(method == "Semirigid Chamber"),
    
    # Medians
    rigid_median = median(flux_value[method == "Rigid Chamber"], na.rm = TRUE),
    semirigid_median = median(flux_value[method == "Semirigid Chamber"], na.rm = TRUE),
    
    # Quartiles for IQR overlap calculation
    rigid_q25 = quantile(flux_value[method == "Rigid Chamber"], 0.25, na.rm = TRUE),
    rigid_q75 = quantile(flux_value[method == "Rigid Chamber"], 0.75, na.rm = TRUE),
    semirigid_q25 = quantile(flux_value[method == "Semirigid Chamber"], 0.25, na.rm = TRUE),
    semirigid_q75 = quantile(flux_value[method == "Semirigid Chamber"], 0.75, na.rm = TRUE),
    
    # Min/Max for range overlap
    rigid_min = min(flux_value[method == "Rigid Chamber"], na.rm = TRUE),
    rigid_max = max(flux_value[method == "Rigid Chamber"], na.rm = TRUE),
    semirigid_min = min(flux_value[method == "Semirigid Chamber"], na.rm = TRUE),
    semirigid_max = max(flux_value[method == "Semirigid Chamber"], na.rm = TRUE),
    
    .groups = 'drop'
  ) %>%
  mutate(
    # Calculate IQR overlap
    iqr_overlap_start = pmax(rigid_q25, semirigid_q25),
    iqr_overlap_end = pmin(rigid_q75, semirigid_q75),
    iqr_overlap_width = pmax(0, iqr_overlap_end - iqr_overlap_start),
    
    # Calculate individual IQR widths
    rigid_iqr = rigid_q75 - rigid_q25,
    semirigid_iqr = semirigid_q75 - semirigid_q25,
    
    # Overlap as percentage of average IQR
    iqr_overlap_percent = (iqr_overlap_width / ((rigid_iqr + semirigid_iqr) / 2)) * 100,
    
    # Range overlap
    range_overlap_start = pmax(rigid_min, semirigid_min),
    range_overlap_end = pmin(rigid_max, semirigid_max),
    range_overlap_width = pmax(0, range_overlap_end - range_overlap_start),
    
    rigid_range = rigid_max - rigid_min,
    semirigid_range = semirigid_max - semirigid_min,
    range_overlap_percent = (range_overlap_width / ((rigid_range + semirigid_range) / 2)) * 100,
    
    # Fold difference (rigid / semirigid)
    median_fold_difference = rigid_median / semirigid_median,
    
    # Overlap categories
    overlap_category = case_when(
      iqr_overlap_percent > 50 ~ "High overlap",
      iqr_overlap_percent > 25 ~ "Moderate overlap",
      iqr_overlap_percent > 0 ~ "Low overlap",
      TRUE ~ "No overlap"
    )
  )

cat("Individual Species-Gas Combinations (Matched Subset):\n")
print(matched_overlap_stats %>% 
        dplyr::select(species, gas_type, n_rigid, n_semirigid, median_fold_difference, 
               iqr_overlap_percent, range_overlap_percent, overlap_category) %>%
        arrange(gas_type, desc(iqr_overlap_percent)))

# Summary statistics for matched subset
cat("\n=== MATCHED SUBSET SUMMARY STATISTICS ===\n")

matched_summary_by_gas <- matched_overlap_stats %>%
  group_by(gas_type) %>%
  summarise(
    n_species = n(),
    mean_fold_diff = mean(median_fold_difference, na.rm = TRUE),
    median_fold_diff = median(median_fold_difference, na.rm = TRUE),
    mean_iqr_overlap = mean(iqr_overlap_percent, na.rm = TRUE),
    median_iqr_overlap = median(iqr_overlap_percent, na.rm = TRUE),
    mean_range_overlap = mean(range_overlap_percent, na.rm = TRUE),
    high_overlap_count = sum(overlap_category == "High overlap"),
    moderate_overlap_count = sum(overlap_category == "Moderate overlap"),
    low_overlap_count = sum(overlap_category == "Low overlap"),
    no_overlap_count = sum(overlap_category == "No overlap"),
    .groups = 'drop'
  )

cat("Matched Subset Summary by Gas Type:\n")
print(matched_summary_by_gas)

matched_overall_summary <- matched_overlap_stats %>%
  summarise(
    total_comparisons = n(),
    mean_fold_diff = mean(median_fold_difference, na.rm = TRUE),
    mean_iqr_overlap = mean(iqr_overlap_percent, na.rm = TRUE),
    mean_range_overlap = mean(range_overlap_percent, na.rm = TRUE),
    percent_high_overlap = sum(overlap_category == "High overlap") / n() * 100,
    percent_moderate_plus = sum(overlap_category %in% c("High overlap", "Moderate overlap")) / n() * 100
  )

cat("\nMatched Subset Overall Summary:\n")
print(matched_overall_summary)

# Species ranking correlations for matched subset
cat("\n=== MATCHED SUBSET SPECIES RANKING CORRELATIONS ===\n")

matched_species_ranks <- matched_subset %>%
  group_by(species, gas_type, method) %>%
  summarise(median_flux = median(flux_value, na.rm = TRUE), .groups = 'drop') %>%
  # Only include species with both methods
  group_by(species, gas_type) %>%
  filter(n() == 2) %>%
  pivot_wider(names_from = method, values_from = median_flux) %>%
  group_by(gas_type) %>%
  summarise(
    spearman_cor = cor(`Rigid Chamber`, `Semirigid Chamber`, method = "spearman", use = "complete.obs"),
    pearson_cor = cor(`Rigid Chamber`, `Semirigid Chamber`, method = "pearson", use = "complete.obs"),
    n_species = n(),
    .groups = 'drop'
  )

cat("Matched Subset Species Median Flux Correlations:\n")
print(matched_species_ranks)

# Key numbers for manuscript (matched subset)
cat("\n=== MATCHED SUBSET KEY NUMBERS FOR MANUSCRIPT ===\n")
cat(sprintf("Rigid chambers measured %.1fx higher CO2 fluxes (median, matched subset)\n", 
            matched_summary_by_gas$median_fold_diff[matched_summary_by_gas$gas_type == "CO2"]))
cat(sprintf("Rigid chambers measured %.1fx higher CH4 fluxes (median, matched subset)\n", 
            matched_summary_by_gas$median_fold_diff[matched_summary_by_gas$gas_type == "CH4"]))
cat(sprintf("Mean IQR overlap (matched subset): %.1f%%\n", matched_overall_summary$mean_iqr_overlap))
cat(sprintf("Mean range overlap (matched subset): %.1f%%\n", matched_overall_summary$mean_range_overlap))
cat(sprintf("Percentage with moderate+ overlap (matched subset): %.1f%%\n", matched_overall_summary$percent_moderate_plus))
cat(sprintf("CO2 species ranking correlation (matched): r = %.2f\n", matched_species_ranks$spearman_cor[matched_species_ranks$gas_type == "CO2"]))
cat(sprintf("CH4 species ranking correlation (matched): r = %.2f\n", matched_species_ranks$spearman_cor[matched_species_ranks$gas_type == "CH4"]))

# Compare with full dataset
cat("\n=== COMPARISON: FULL vs MATCHED SUBSET ===\n")
cat("FULL DATASET:\n")
cat(sprintf("  CO2 fold difference: %.1fx\n", summary_by_gas$median_fold_diff[summary_by_gas$gas_type == "CO2"]))
cat(sprintf("  CH4 fold difference: %.1fx\n", summary_by_gas$median_fold_diff[summary_by_gas$gas_type == "CH4"]))
cat(sprintf("  Mean IQR overlap: %.1f%%\n", overall_summary$mean_iqr_overlap))
cat(sprintf("  CO2 correlation: r = %.2f\n", species_ranks$spearman_cor[species_ranks$gas_type == "CO2"]))
cat(sprintf("  CH4 correlation: r = %.2f\n", species_ranks$spearman_cor[species_ranks$gas_type == "CH4"]))

cat("\nMATCHED SUBSET:\n")
cat(sprintf("  CO2 fold difference: %.1fx\n", matched_summary_by_gas$median_fold_diff[matched_summary_by_gas$gas_type == "CO2"]))
cat(sprintf("  CH4 fold difference: %.1fx\n", matched_summary_by_gas$median_fold_diff[matched_summary_by_gas$gas_type == "CH4"]))
cat(sprintf("  Mean IQR overlap: %.1f%%\n", matched_overall_summary$mean_iqr_overlap))
cat(sprintf("  CO2 correlation: r = %.2f\n", matched_species_ranks$spearman_cor[matched_species_ranks$gas_type == "CO2"]))
cat(sprintf("  CH4 correlation: r = %.2f\n", matched_species_ranks$spearman_cor[matched_species_ranks$gas_type == "CH4"]))

cat("\nConclusion: Use MATCHED SUBSET numbers for manuscript as they control for confounding.\n")


# Final Ridge Plot with dplyr::selective Densities from Highlighted Points Only
# Comparing rigid vs semirigid chamber methods

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(scales)
library(lubridate)

# Load the datasets with your file paths
rigid_data <- read.csv('../../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv', stringsAsFactors = FALSE)
semirigid_data <- read.csv('../../../data/processed/flux/semirigid_tree_final_complete_dataset.csv', stringsAsFactors = FALSE)
species_mapping <- read.csv('../../../data/raw/lgr/semirigid_2020-2021/Yale Myers Methane Project/spatial_data/YM_trees_measured.csv.csv', stringsAsFactors = FALSE)

# Data preparation for rigid dataset
rigid_clean <- rigid_data %>%
  filter(!is.na(species) & species != "" & species != "NA") %>%
  mutate(
    CO2_best.flux = as.numeric(ifelse(CO2_best.flux == "NA", NA, CO2_best.flux)),
    CH4_best.flux = as.numeric(ifelse(CH4_best.flux == "NA", NA, CH4_best.flux)),
    date = as.Date(start.time),
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring", 
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall"
    ),
    color_category = ifelse(measurement_height == 125, "rigid_highlight", "rigid_other")
  ) %>%
  filter(!is.na(CO2_best.flux) | !is.na(CH4_best.flux)) %>%
  mutate(method = "Rigid Chamber")

# Data preparation for semirigid dataset
semirigid_clean <- semirigid_data %>%
  mutate(Plot.Tag = as.character(Plot.Tag)) %>%
  left_join(species_mapping, by = c("Plot.Tag" = "Label")) %>%
  filter(!is.na(Species)) %>%
  rename(species = Species) %>%
  mutate(
    date = as.Date(Date),
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer", 
      month %in% c(9, 10, 11) ~ "Fall"
    ),
    CH4_best.flux = CH4_best.flux.x,
    color_category = ifelse(
      month %in% c(6, 7, 8) & Site %in% c("Upland", "Intermediate"), 
      "semirigid_highlight", 
      "semirigid_other"
    )
  ) %>%
  filter(!is.na(CO2_best.flux) | !is.na(CH4_best.flux)) %>%
  mutate(method = "Semirigid Chamber") %>%
  dplyr::select(UniqueID, species, Site, CO2_best.flux, CH4_best.flux, method, color_category)

# Combine datasets
combined_data <- bind_rows(
  rigid_clean %>% dplyr::select(UniqueID, species, CO2_best.flux, CH4_best.flux, method, color_category),
  semirigid_clean %>% dplyr::select(UniqueID, species, CO2_best.flux, CH4_best.flux, method, color_category)
)

# Create long format for plotting
combined_long <- combined_data %>%
  pivot_longer(
    cols = c(CO2_best.flux, CH4_best.flux),
    names_to = "gas_type",
    values_to = "flux_value"
  ) %>%
  filter(!is.na(flux_value)) %>%
  mutate(
    gas_type = case_when(
      gas_type == "CO2_best.flux" ~ "CO2",
      gas_type == "CH4_best.flux" ~ "CH4"
    ),
    gas_with_units = case_when(
      gas_type == "CO2" ~ "CO2 (μmol/m²/s)",
      gas_type == "CH4" ~ "CH4 (nmol/m²/s)"
    )
  )

# Get all species with n>5 individuals for each gas type
species_with_sufficient_data <- combined_long %>%
  group_by(species, gas_type) %>%
  summarise(n_obs = n(), .groups = 'drop') %>%
  filter(n_obs > 5) %>%
  pull(species) %>%
  unique()

cat("Species with >5 observations:", paste(species_with_sufficient_data, collapse = ", "), "\n")

# STANDARD PLOT: Ridges from all data, dplyr::selective point colors
cat("=== STANDARD PLOT ===\n")

# Prepare data for individual species (both gases) - filter by n>5
species_data <- combined_long %>%
  filter(species %in% species_with_sufficient_data) %>%
  mutate(group_label = species)

# Prepare pooled data (both gases)
pooled_data <- combined_long %>%
  mutate(group_label = "ALL SPECIES POOLED")

# Combine data
combined_plot_data <- bind_rows(species_data, pooled_data) %>%
  mutate(
    group_label = factor(group_label, 
                         levels = c("ALL SPECIES POOLED", sort(species_with_sufficient_data)))
  )

# Point color palette
point_color_palette <- c(
  "rigid_highlight" = "#440154",  # Purple (matching ridge)
  "rigid_other" = "#CCCCCC",      # Grey
  "semirigid_highlight" = "#31688e",   # Teal (matching ridge)
  "semirigid_other" = "#CCCCCC"        # Grey
)

standard_plot <- ggplot(combined_plot_data, aes(x = flux_value, y = group_label)) +
  # Ridge densities from all data
  geom_density_ridges(
    aes(fill = method, color = method),
    alpha = 0.4,
    scale = 0.8,
    size = 0.8,
    adjust = 2
  ) +
  # Points with dplyr::selective coloring
  geom_point(
    aes(color = color_category),
    position = position_jitter(width = 0, height = 0.12),
    size = 0.8,
    alpha = 0.7
  ) +
  facet_wrap(~gas_with_units, scales = "free_x", ncol = 2) +
  scale_x_continuous(trans = "pseudo_log") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey60", alpha = 0.7) +
  scale_fill_manual(
    values = c("Rigid Chamber" = "#440154", "Semirigid Chamber" = "#31688e"),
    name = "Method"
  ) +
  scale_color_manual(
    values = point_color_palette,
    name = "Sample Type",
    labels = c(
      "rigid_highlight" = "Rigid (125cm)",
      "rigid_other" = "Rigid (other)",
      "semirigid_highlight" = "Semirigid\n(Summer/Upland/Int.)",
      "semirigid_other" = "Semirigid (other)"
    )
  ) +
  theme_ridges(grid = TRUE) +
  theme(
    legend.position = "bottom",
    axis.title = NULL,
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_line(color = "grey90", size = 0.3),
    legend.box = "horizontal"
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.7, size = 1), nrow = 1),
    color = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)
  )

print(standard_plot)

# dplyr::selectIVE RIDGES PLOT: Ridges only from highlighted points
cat("\n=== dplyr::selectIVE RIDGES PLOT ===\n")
cat("Ridges computed only from highlighted samples, all points shown\n")

# Create datasets for ridge densities (highlighted points only) and points (all points)
ridge_data <- combined_plot_data %>%
  filter(color_category %in% c("rigid_highlight", "semirigid_highlight"))

all_points_data <- combined_plot_data

dplyr::selective_plot <- ggplot() +
  # Ridge densities from highlighted data only
  geom_density_ridges(
    data = ridge_data,
    aes(x = flux_value, y = group_label, fill = method, color = method),
    alpha = 0.4,
    scale = 0.8,
    size = 0.8,
    adjust = 2
  ) +
  # All points with dplyr::selective coloring
  geom_point(
    data = all_points_data,
    aes(x = flux_value, y = group_label, color = color_category),
    position = position_jitter(width = 0, height = 0.12),
    size = 0.8,
    alpha = 0.7
  ) +
  facet_wrap(~gas_with_units, scales = "free_x", ncol = 2) +
  scale_x_continuous(trans = "pseudo_log") +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey60", alpha = 0.7) +
  scale_fill_manual(
    values = c("Rigid Chamber" = "#440154", "Semirigid Chamber" = "#31688e"),
    name = "Method"
  ) +
  scale_color_manual(
    values = point_color_palette,
    name = "Sample Type",
    labels = c(
      "rigid_highlight" = "Rigid (125cm)",
      "rigid_other" = "Rigid (other)",
      "semirigid_highlight" = "Semirigid\n(Summer/Upland/Int.)",
      "semirigid_other" = "Semirigid (other)"
    )
  ) +
  theme_ridges(grid = TRUE) +
  theme(
    legend.position = "bottom",
    axis.title = NULL,
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_line(color = "grey90", size = 0.3),
    legend.box = "horizontal"
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 0.7, size = 1), nrow = 1),
    color = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)
  )

print(dplyr::selective_plot)

# Summary
cat("\nSummary:\n")
cat("Standard plot: Ridges from all data, dplyr::selective point colors\n")
cat("dplyr::selective plot: Ridges only from highlighted samples (Rigid 125cm + Semirigid Summer/Upland/Int.), all points shown\n")
cat("Species included:", paste(species_with_sufficient_data, collapse = ", "), "\n")