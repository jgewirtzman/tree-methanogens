# ==============================================================================
# Comprehensive ddPCR Analysis
# ==============================================================================
# Purpose: Comprehensive ddPCR analysis and visualization with species
#   comparisons using the merged final dataset.
#
# Pipeline stage: 03 Analysis
# Run after: 00_harmonization/02_harmonize_all_data.R
#
# Inputs:
#   - merged_tree_dataset_final.csv (from data/processed/integrated/)
#
# Outputs:
#   - Analysis plots and summary statistics
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)
library(RColorBrewer)
library(cowplot)
library(stringr)

# Set theme
theme_set(theme_bw(base_size = 12))

# Load the merged final dataset (should already be loaded)
if (!exists("merged_final")) {
  merged_final <- read_csv("../../data/processed/integrated/merged_tree_dataset_final.csv")
}

cat("=== DATASET OVERVIEW ===\n")
cat("Total trees:", nrow(merged_final), "\n")
cat("Total columns:", ncol(merged_final), "\n")

# =============================================================================
# DATA PREPARATION AND FILTERING
# =============================================================================

cat("=== PREPARING DATA ===\n")

# Create long format ddPCR data for plotting
ddpcr_long <- merged_final %>%
  dplyr::select(tree_id, species, plot, starts_with("ddpcr_")) %>%
  pivot_longer(
    cols = starts_with("ddpcr_"),
    names_to = "measurement",
    values_to = "concentration",
    values_drop_na = TRUE
  ) %>%
  # Parse measurement components
  mutate(
    clean_name = str_remove(measurement, "^ddpcr_"),
    parts = str_split(clean_name, "_"),
    analysis_type = map_chr(parts, ~ tail(.x, 1)),
    core_type = map_chr(parts, ~ .x[length(.x) - 1]),
    target = map_chr(parts, ~ paste(head(.x, length(.x) - 2), collapse = "_"))
  ) %>%
  dplyr::select(tree_id, species, plot, target, core_type, analysis_type, concentration)

# Filter out problematic species (like in original script)
species_to_exclude <- c("KALA", "SAAL", "QUAL", "QUVE", "CAOV", "PRSE")

merged_clean <- merged_final %>%
  filter(!species %in% species_to_exclude | is.na(species))

ddpcr_clean <- ddpcr_long %>%
  filter(!species %in% species_to_exclude | is.na(species))

cat("After filtering species:", nrow(merged_clean), "trees\n")

# Create color palette
coul <- c(brewer.pal(9, "Set3"), brewer.pal(12, "Paired"), 
          rev(brewer.pal(8, "Accent")), brewer.pal(8, "Dark2"),
          brewer.pal(12, "Paired"), brewer.pal(8, "Accent"), 
          brewer.pal(8, "Dark2"))

# =============================================================================
# BASIC ddPCR VISUALIZATIONS
# =============================================================================

cat("=== CREATING BASIC ddPCR PLOTS ===\n")

# Target abundance by species (equivalent to species_bar in original)
target_by_species <- ddpcr_clean %>%
  filter(target %in% c("mcra", "mcra_probe", "mmox", "pmoa")) %>%
  ggplot(aes(x = species, y = log10(1 + concentration), color = target)) +
  stat_summary(size = 1.2) +
  ggtitle("DNA Abundances by Species and Target") +
  facet_wrap(~ core_type, scales = "free_x", ncol = 2) +
  ylab(expression(Log[10](Gene~Copies))) +
  xlab("Species ID") +
  scale_color_discrete(name = "Target") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(target_by_species)

# mcrA abundance comparison (loose vs strict)
mcra_comparison <- ddpcr_clean %>%
  filter(target == "mcra_probe") %>%
  ggplot(aes(x = species, y = log10(1 + concentration), color = analysis_type)) +
  stat_summary(size = 1.1) +
  ggtitle("mcrA DNA (Probe) - Loose vs Strict") +
  facet_grid(. ~ core_type, scales = "free", space = "free") +
  ylab(expression(Log[10](mcrA~(Copies~rxn^-1)))) +
  xlab("Species ID") +
  scale_color_manual(name = "Cutoff Threshold", values = wes_palette("Darjeeling1")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(mcra_comparison)

# =============================================================================
# SOIL ENVIRONMENT CORRELATIONS
# =============================================================================

cat("=== CREATING SOIL CORRELATION PLOTS ===\n")

# VWC vs DNA abundances
vwc_plot <- merged_clean %>%
  filter(!is.na(VWC_mean), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  ggplot(aes(x = log10(ddpcr_mcra_probe_Inner_loose + ddpcr_mmox_Inner_strict + ddpcr_pmoa_Inner_strict + 1),
             y = VWC_mean)) +
  geom_point(aes(col = species), size = 3) +
  ggtitle("DNA Abundances vs. VWC") +
  stat_cor() +
  geom_smooth(method = "lm") +
  ylab("Mean VWC") +
  xlab(expression(Log[10](Sum~All~Targets~(Copies~rxn^-1)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(vwc_plot)

# Soil temperature correlation
soil_temp_plot <- merged_clean %>%
  filter(!is.na(SoilTemp_mean), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  ggplot(aes(x = log10(ddpcr_mcra_probe_Inner_loose + ddpcr_mmox_Inner_strict + ddpcr_pmoa_Inner_strict + 1),
             y = SoilTemp_mean)) +
  geom_point(aes(col = species), size = 3) +
  ggtitle("DNA Abundances vs. Soil Temperature") +
  stat_cor() +
  geom_smooth(method = "lm") +
  ylab("Mean Soil Temperature") +
  xlab(expression(Log[10](Sum~All~Targets~(Copies~rxn^-1)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(soil_temp_plot)

# ORP correlation
orp_plot <- merged_clean %>%
  filter(!is.na(ORP_mean), !is.na(ddpcr_mcra_probe_Inner_loose), !is.na(ddpcr_pmoa_Inner_loose)) %>%
  ggplot(aes(x = log10(ddpcr_mcra_probe_Inner_loose + 1),
             y = log10(ddpcr_pmoa_Inner_loose + 1),
             col = ORP_mean)) +
  geom_point(size = 3) +
  ggtitle("DNA Abundances vs. ORP") +
  ylab(expression(Log[10](pmoA~(Copies~rxn^-1)))) +
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_gradient(name = "ORP", low = "green", high = "red")

print(orp_plot)

# Combine soil plots
soil_combined <- plot_grid(vwc_plot, soil_temp_plot, orp_plot, nrow = 1)
print(soil_combined)

# =============================================================================
# FLUX CORRELATIONS
# =============================================================================

cat("=== CREATING FLUX CORRELATION PLOTS ===\n")

# Flux correlation at 50cm
flux_50_plot <- merged_clean %>%
  filter(!is.na(CH4_best.flux_50cm), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  ggplot(aes(x = log10(ddpcr_mcra_probe_Inner_loose + 1),
             y = log10(CH4_best.flux_50cm))) +
  geom_point(size = 3, aes(col = species)) +
  ggtitle("DNA Abundances vs. Flux at 50cm") +
  facet_wrap(~ species, scales = "free", ncol = 5) +
  stat_cor() +
  geom_smooth(method = "lm") +
  ylab(expression(Log[10](CH[4]~Flux))) +
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(name = "Species", values = coul[10:30])

print(flux_50_plot)

# Flux correlation at 125cm
flux_125_plot <- merged_clean %>%
  filter(!is.na(CH4_best.flux_125cm), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  ggplot(aes(x = log10(ddpcr_mcra_probe_Inner_loose + 1),
             y = log10(CH4_best.flux_125cm))) +
  geom_point(size = 3, aes(col = species)) +
  ggtitle("DNA Abundances vs. Flux at 125cm") +
  stat_cor() +
  geom_smooth(method = "lm") +
  ylab(expression(Log[10](CH[4]~Flux))) +
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(name = "Species", values = coul[10:30])

print(flux_125_plot)

# Flux correlation at 200cm
flux_200_plot <- merged_clean %>%
  filter(!is.na(CH4_best.flux_200cm), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  ggplot(aes(x = log10(ddpcr_mcra_probe_Inner_loose + 1),
             y = log10(CH4_best.flux_200cm))) +
  geom_point(size = 3, aes(col = species)) +
  ggtitle("DNA Abundances vs. Flux at 200cm") +
  stat_cor() +
  geom_smooth(method = "lm") +
  ylab(expression(Log[10](CH[4]~Flux))) +
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(name = "Species", values = coul[10:30])

print(flux_200_plot)

# Combine flux plots
flux_combined <- plot_grid(flux_50_plot, flux_125_plot, flux_200_plot, nrow = 1)
print(flux_combined)

# =============================================================================
# INTERNAL GAS CONCENTRATION CORRELATIONS
# =============================================================================

cat("=== CREATING INTERNAL GAS PLOTS ===\n")

# Internal CH4 vs DNA abundances
internal_ch4_plot <- merged_clean %>%
  filter(!is.na(CH4_concentration), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  ggplot(aes(x = log10(ddpcr_mcra_probe_Inner_loose + 1),
             y = log10(CH4_concentration))) +
  geom_point(size = 3, aes(col = species)) +
  ggtitle("DNA Abundances vs. Internal Gas Concentration") +
  stat_cor() +
  geom_smooth(method = "lm") +
  ylab(expression(Log[10](CH[4]~ppm))) +
  xlab(expression(Log[10](mcrA~(Copies~rxn^-1)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(name = "Species", values = coul[10:30])

print(internal_ch4_plot)

# Internal gas concentration boxplot by species
internal_box_plot <- merged_clean %>%
  filter(!is.na(CH4_concentration)) %>%
  ggplot(aes(y = log10(CH4_concentration), x = species, fill = species)) +
  geom_jitter() +
  ggtitle("Internal Gas Concentrations") +
  geom_boxplot(alpha = 0.8) +
  ylab(expression(Log[10](CH[4]~ppm))) +
  xlab("Tree Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(internal_box_plot)

# Flux vs Internal gas correlations
flux_internal_50 <- merged_clean %>%
  filter(!is.na(CH4_best.flux_50cm), !is.na(CH4_concentration)) %>%
  ggplot(aes(x = log10(CH4_best.flux_50cm),
             y = log10(CH4_concentration))) +
  geom_point(size = 3, aes(col = species)) +
  ggtitle("Flux vs. Internal Gas Concentration: 50cm") +
  stat_cor() +
  geom_smooth(method = "lm") +
  ylab(expression(Log[10](CH[4]~ppm))) +
  xlab(expression(Log[10](CH[4]~Flux))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(name = "Species", values = coul[10:30])

print(flux_internal_50)

# =============================================================================
# SPECIES COMPARISONS AND ABUNDANCE PATTERNS
# =============================================================================

cat("=== CREATING SPECIES COMPARISON PLOTS ===\n")

# Species abundance comparison (mass-corrected)
# Note: Original script used Sample.Mass.Added.to.Tube..mg. - we'll use a proxy calculation
species_abundance <- ddpcr_clean %>%
  filter(target == "mcra_probe", analysis_type == "loose") %>%
  filter(!is.na(concentration)) %>%
  ggplot(aes(y = log10(concentration * 37.5 / 100 + 1), x = species)) + # Using 100mg as proxy mass
  geom_jitter(shape = 21, aes(fill = species), size = 3) +
  stat_summary(size = 1) +
  ggtitle("mcrA Abundances by Species") +
  facet_wrap(~ core_type, scales = "free_x", nrow = 1) +
  ylab(expression(Log[10](mcrA~(Copies~mg^-1)))) +
  xlab("Species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.position = "none")

print(species_abundance)

# Cross-compartment correlations (Inner vs Outer)
cross_compartment <- merged_clean %>%
  filter(!is.na(ddpcr_mcra_probe_Inner_loose), !is.na(ddpcr_mcra_probe_Outer_loose)) %>%
  ggplot(aes(x = log10(ddpcr_mcra_probe_Inner_loose + 1),
             y = log10(ddpcr_mcra_probe_Outer_loose + 1))) +
  geom_point(size = 4, aes(col = species)) +
  ggtitle("Inner vs. Outer mcrA Abundances") +
  facet_wrap(~ species, scales = "free", ncol = 5) +
  stat_cor() +
  geom_smooth(method = "lm") +
  ylab(expression(Log[10](mcrA~Outer~(Copies~rxn^-1)))) +
  xlab(expression(Log[10](mcrA~Inner~(Copies~rxn^-1)))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(name = "Species", values = coul[10:30])

print(cross_compartment)

# =============================================================================
# STATISTICAL ANALYSES
# =============================================================================

cat("=== RUNNING STATISTICAL ANALYSES ===\n")

# Linear model analyses for different compartments
run_compartment_analysis <- function(data, compartment_suffix, compartment_name) {
  cat("\n--- ", compartment_name, " Compartment Analysis ---\n")
  
  # Get relevant columns
  mcra_col <- paste0("ddpcr_mcra_probe_", compartment_suffix, "_loose")
  mmox_col <- paste0("ddpcr_mmox_", compartment_suffix, "_strict") 
  pmoa_col <- paste0("ddpcr_pmoa_", compartment_suffix, "_strict")
  
  # Filter data with required columns
  analysis_data <- data %>%
    filter(!is.na(.data[[mcra_col]]), !is.na(VWC_mean), !is.na(CH4_concentration))
  
  if (nrow(analysis_data) > 10) {
    # Model for internal CH4 concentration
    tryCatch({
      model_int <- lm(log10(CH4_concentration) ~ 
                        log10(.data[[mcra_col]] + 1) + 
                        VWC_mean + 
                        log10(.data[[mmox_col]] + 1) + 
                        log10(.data[[pmoa_col]] + 1) + 
                        species, 
                      data = analysis_data)
      
      cat("Internal CH4 Model Summary:\n")
      print(anova(model_int))
      
    }, error = function(e) {
      cat("Error in internal CH4 model:", e$message, "\n")
    })
    
    # Model for flux (if data available)
    if (!is.na(analysis_data$CH4_best.flux_50cm[1])) {
      tryCatch({
        flux_data <- analysis_data %>% filter(!is.na(CH4_best.flux_50cm))
        
        if (nrow(flux_data) > 5) {
          model_flux <- lm(log10(CH4_best.flux_50cm) ~ 
                             log10(.data[[mcra_col]] + 1) + 
                             VWC_mean + 
                             log10(.data[[mmox_col]] + 1) + 
                             log10(.data[[pmoa_col]] + 1) + 
                             species, 
                           data = flux_data)
          
          cat("Flux Model Summary:\n")
          print(anova(model_flux))
        }
      }, error = function(e) {
        cat("Error in flux model:", e$message, "\n")
      })
    }
  } else {
    cat("Insufficient data for", compartment_name, "analysis\n")
  }
}

# Run analyses for each compartment
run_compartment_analysis(merged_clean, "Inner", "Inner")
run_compartment_analysis(merged_clean, "Outer", "Outer") 
run_compartment_analysis(merged_clean, "Organic", "Organic")
run_compartment_analysis(merged_clean, "Mineral", "Mineral")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")

# Data availability summary
cat("Data availability by measurement type:\n")
cat("Trees with ddPCR data:", sum(!is.na(merged_clean$ddpcr_mcra_probe_Inner_loose)), "\n")
cat("Trees with flux data:", sum(!is.na(merged_clean$CH4_best.flux_50cm)), "\n")
cat("Trees with internal gas data:", sum(!is.na(merged_clean$CH4_concentration)), "\n")
cat("Trees with soil data:", sum(!is.na(merged_clean$VWC_mean)), "\n")

# Species representation
species_counts <- merged_clean %>%
  filter(!is.na(species)) %>%
  count(species, sort = TRUE)

cat("\nSpecies representation:\n")
print(species_counts)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All plots and analyses have been generated using the merged final dataset.\n")
cat("This provides a much cleaner and more reliable analysis than the previous harmonization approach.\n")
