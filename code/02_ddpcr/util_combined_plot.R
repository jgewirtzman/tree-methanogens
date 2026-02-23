# ==============================================================================
# Combined Gene Abundance Plot (Figure 4)
# ==============================================================================
# Purpose: Combined visualization of multiple gene types (barplot + scatter).
#
# Dependencies:
#   - 04_species_barplots.R must be sourced first (defines plotting functions
#     and species_mapping)
#
# Inputs:
#   - data/processed/integrated/merged_tree_dataset_final.csv
#
# Outputs:
#   - outputs/figures/main/fig4_methanogen_methanotroph_abundance.png
# ==============================================================================

library(tidyverse)
library(cowplot)

# Load data if not already in environment
if (!exists("merged_final")) {
  merged_final <- read_csv("data/processed/integrated/merged_tree_dataset_final.csv",
                           show_col_types = FALSE)
  cat("Loaded merged_final:", nrow(merged_final), "rows\n")
}

# Load plotting functions and species_mapping if not already defined
if (!exists("species_mapping") || !exists("create_mcra_barplot_by_species")) {
  source("code/02_ddpcr/04_species_barplots.R")
  cat("Sourced 04_species_barplots.R\n")
}

# Load scatter plot function if not already defined
if (!exists("create_gene_scatter_ggside_transformed_probe_mcra")) {
  source("code/02_ddpcr/util_ridge_plots.R")
  cat("Sourced util_ridge_plots.R\n")
}

# Generate your plots
result <- create_mcra_barplot_by_species(merged_final, species_mapping)
barplot <- result$plot

scatterplot <- create_gene_scatter_ggside_transformed_probe_mcra(merged_final)


# Option 1: Smaller text and better spacing
barplot_improved <- result$plot +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),  # Smaller text
    panel.spacing = unit(1, "lines")  # More space between facets
  )



# Then use the improved version:
combined_plot <- plot_grid(barplot_improved, scatterplot, ncol = 2,
                           labels = c("(a)", "(b)"),
                           label_size = 11,
                           label_fontface = "bold")

print(combined_plot)

# Save combined figure
ggsave("outputs/figures/main/fig4_methanogen_methanotroph_abundance.png",
       combined_plot, width = 12, height = 5.3, dpi = 300)



# ===== EXTRACT STATISTICS FOR METHANOGEN QUANTIFICATION SECTION =====

# Get summary statistics for each compartment
compartment_summary <- merged_final %>%
  select(tree_id, species_id, 
         ddpcr_mcra_probe_Inner_loose, ddpcr_mcra_probe_Outer_loose,
         ddpcr_mcra_probe_Mineral_loose, ddpcr_mcra_probe_Organic_loose,
         ddpcr_pmoa_Inner_loose, ddpcr_pmoa_Outer_loose,
         ddpcr_pmoa_Mineral_loose, ddpcr_pmoa_Organic_loose,
         ddpcr_mmox_Inner_loose, ddpcr_mmox_Outer_loose,
         ddpcr_mmox_Mineral_loose, ddpcr_mmox_Organic_loose) %>%
  pivot_longer(cols = starts_with("ddpcr"), names_to = "measurement", values_to = "copies") %>%
  filter(!is.na(copies), copies > 0) %>%
  separate(measurement, into = c("method", "gene", "probe", "location", "stringency")) %>%
  mutate(
    compartment = case_when(
      location == "Inner" ~ "Heartwood",
      location == "Outer" ~ "Sapwood", 
      location == "Mineral" ~ "Mineral soil",
      location == "Organic" ~ "Organic soil"
    )
  ) %>%
  group_by(gene, compartment) %>%
  summarise(
    n = n(),
    mean_copies = mean(copies, na.rm = TRUE),
    sd_copies = sd(copies, na.rm = TRUE),
    median_copies = median(copies, na.rm = TRUE),
    max_copies = max(copies, na.rm = TRUE),
    pct_above_500 = sum(copies >= 500) / n() * 100,
    .groups = 'drop'
  )

cat("\n===== METHANOGEN (mcrA) STATISTICS BY COMPARTMENT =====\n")
print(compartment_summary %>% filter(gene == "mcra"))

cat("\n===== METHANOTROPH STATISTICS BY COMPARTMENT =====\n")
print(compartment_summary %>% filter(gene %in% c("pmoa", "mmox")))

# Count total trees sampled
n_trees <- merged_final %>%
  filter(!is.na(ddpcr_mcra_probe_Inner_loose) | !is.na(ddpcr_mcra_probe_Outer_loose)) %>%
  nrow()
cat("\nTotal trees with mcrA measurements:", n_trees, "\n")

# Get species count
n_species <- merged_final %>%
  filter(!is.na(ddpcr_mcra_probe_Inner_loose) | !is.na(ddpcr_mcra_probe_Outer_loose)) %>%
  pull(species_id) %>%
  unique() %>%
  length()
cat("Number of species:", n_species, "\n")







# ===== FIXED CODE FOR METHANOGEN QUANTIFICATION STATISTICS =====

# Fix the methanotroph compartment issue by processing them separately
compartment_summary_mcra <- merged_final %>%
  select(tree_id, species_id, 
         Inner = ddpcr_mcra_probe_Inner_loose, 
         Outer = ddpcr_mcra_probe_Outer_loose,
         Mineral = ddpcr_mcra_probe_Mineral_loose, 
         Organic = ddpcr_mcra_probe_Organic_loose) %>%
  pivot_longer(cols = c(Inner, Outer, Mineral, Organic), 
               names_to = "location", values_to = "copies") %>%
  filter(!is.na(copies)) %>%
  mutate(
    compartment = case_when(
      location == "Inner" ~ "Heartwood",
      location == "Outer" ~ "Sapwood", 
      location == "Mineral" ~ "Mineral soil",
      location == "Organic" ~ "Organic soil"
    ),
    # Add 1 to handle zeros for log transformation
    log_copies = log10(copies + 1)
  ) %>%
  group_by(compartment) %>%
  summarise(
    n_total = n(),
    n_positive = sum(copies > 0),
    mean_log = mean(log_copies[copies > 0], na.rm = TRUE),
    sd_log = sd(log_copies[copies > 0], na.rm = TRUE),
    median_log = median(log_copies[copies > 0], na.rm = TRUE),
    max_log = max(log_copies[copies > 0], na.rm = TRUE),
    pct_above_500 = sum(copies >= 500) / n_total * 100,
    pct_detectable = n_positive / n_total * 100,
    .groups = 'drop'
  )

cat("\n===== METHANOGEN (mcrA) LOG-SCALE STATISTICS =====\n")
print(compartment_summary_mcra)

# Process methanotrophs separately with correct compartment assignment
methanotroph_summary <- merged_final %>%
  select(tree_id, 
         # pmoA
         pmoA_Inner = ddpcr_pmoa_Inner_loose,
         pmoA_Outer = ddpcr_pmoa_Outer_loose,
         pmoA_Mineral = ddpcr_pmoa_Mineral_loose,
         pmoA_Organic = ddpcr_pmoa_Organic_loose,
         # mmoX  
         mmoX_Inner = ddpcr_mmox_Inner_loose,
         mmoX_Outer = ddpcr_mmox_Outer_loose,
         mmoX_Mineral = ddpcr_mmox_Mineral_loose,
         mmoX_Organic = ddpcr_mmox_Organic_loose) %>%
  pivot_longer(cols = -tree_id, names_to = "measurement", values_to = "copies") %>%
  filter(!is.na(copies), copies > 0) %>%
  separate(measurement, into = c("gene", "location"), sep = "_") %>%
  mutate(
    compartment = case_when(
      location == "Inner" ~ "Heartwood",
      location == "Outer" ~ "Sapwood",
      location == "Mineral" ~ "Mineral soil", 
      location == "Organic" ~ "Organic soil"
    ),
    log_copies = log10(copies + 1)
  ) %>%
  group_by(compartment, gene) %>%
  summarise(
    n = n(),
    mean_log = mean(log_copies, na.rm = TRUE),
    sd_log = sd(log_copies, na.rm = TRUE),
    median_log = median(log_copies, na.rm = TRUE),
    max_log = max(log_copies, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n===== METHANOTROPH LOG-SCALE STATISTICS =====\n")
print(methanotroph_summary)

# Combined methanotroph abundance (sum pmoA + mmoX)
combined_methanotroph <- merged_final %>%
  mutate(
    methanotroph_Inner = ddpcr_pmoa_Inner_loose + ddpcr_mmox_Inner_loose,
    methanotroph_Outer = ddpcr_pmoa_Outer_loose + ddpcr_mmox_Outer_loose,
    methanotroph_Mineral = ddpcr_pmoa_Mineral_loose + ddpcr_mmox_Mineral_loose,
    methanotroph_Organic = ddpcr_pmoa_Organic_loose + ddpcr_mmox_Organic_loose
  ) %>%
  select(tree_id, starts_with("methanotroph_")) %>%
  pivot_longer(cols = starts_with("methanotroph"), names_to = "location", values_to = "copies") %>%
  filter(!is.na(copies), copies > 0) %>%
  mutate(
    compartment = case_when(
      grepl("Inner", location) ~ "Heartwood",
      grepl("Outer", location) ~ "Sapwood",
      grepl("Mineral", location) ~ "Mineral soil",
      grepl("Organic", location) ~ "Organic soil"
    ),
    log_copies = log10(copies)
  ) %>%
  group_by(compartment) %>%
  summarise(
    n = n(),
    mean_log = mean(log_copies, na.rm = TRUE),
    median_log = median(log_copies, na.rm = TRUE),
    max_log = max(log_copies, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n===== COMBINED METHANOTROPH (pmoA + mmoX) LOG-SCALE =====\n")
print(combined_methanotroph)

# Print in readable format
cat("\n===== SUMMARY IN 10^X FORMAT =====\n")
for(i in 1:nrow(compartment_summary_mcra)) {
  cat(sprintf("%s: mean = 10^%.1f, median = 10^%.1f, max = 10^%.1f (n=%d, %.0f%% detectable)\n",
              compartment_summary_mcra$compartment[i],
              compartment_summary_mcra$mean_log[i],
              compartment_summary_mcra$median_log[i],
              compartment_summary_mcra$max_log[i],
              compartment_summary_mcra$n_total[i],
              compartment_summary_mcra$pct_detectable[i]))
}





# Find species with highest mcrA abundance
species_mcra <- merged_final %>%
  select(species_id, 
         Inner = ddpcr_mcra_probe_Inner_loose,
         Outer = ddpcr_mcra_probe_Outer_loose) %>%
  pivot_longer(cols = c(Inner, Outer), names_to = "location", values_to = "copies") %>%
  filter(!is.na(copies), copies > 0) %>%
  mutate(log_copies = log10(copies + 1)) %>%
  group_by(species_id, location) %>%
  summarise(
    n = n(),
    mean_log = mean(log_copies, na.rm = TRUE),
    max_log = max(log_copies, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(max_log))

cat("\n===== SPECIES WITH HIGHEST mcrA ABUNDANCE =====\n")
print(head(species_mcra, 10))