# ==============================================================================
# Chamber Dimension QC Plots
# ==============================================================================
# Purpose: Quality control visualization of chamber dimensions and consistency
#   across samples.
# ==============================================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)

# Read the combined results
#data <- read_excel("combined_corrected_results.xlsx", sheet = "all_data")
data <- results

# Clean and prepare data
plot_data <- data %>%
  filter(!is.na(Sample) & Sample != "") %>%
  mutate(
    stem_diameter = as.numeric(`D stem (cm)`),
    surface_area = as.numeric(Sc_cm2_corrected),
    volume = as.numeric(Vc_cm3_corrected),
    sample_id = as.factor(Sample)
  ) %>%
  filter(!is.na(stem_diameter) & !is.na(surface_area) & !is.na(volume)) %>%
  select(sample_id, stem_diameter, surface_area, volume)

# Create summary statistics by sample
sample_summary <- plot_data %>%
  group_by(sample_id) %>%
  summarise(
    n_measurements = n(),
    mean_diameter = mean(stem_diameter),
    mean_surface_area = mean(surface_area),
    mean_volume = mean(volume),
    sd_diameter = sd(stem_diameter),
    sd_surface_area = sd(surface_area),
    sd_volume = sd(volume),
    .groups = 'drop'
  ) %>%
  arrange(mean_diameter)

# Print summary statistics
cat("Summary by Sample ID:\n")
cat("Total samples:", nrow(sample_summary), "\n")
cat("Total measurements:", nrow(plot_data), "\n")
cat("Measurements per sample range:", min(sample_summary$n_measurements), "to", max(sample_summary$n_measurements), "\n\n")

print(head(sample_summary, 10))

# 1. Stem Diameter Distribution by Sample
p1 <- ggplot(plot_data, aes(x = reorder(sample_id, stem_diameter), y = stem_diameter, fill = sample_id)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title = "Stem Diameter Distribution by Tree Sample",
    x = "Sample ID (ordered by diameter)",
    y = "Stem Diameter (cm)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold")
  )

# 2. Surface Area Distribution by Sample
p2 <- ggplot(plot_data, aes(x = reorder(sample_id, surface_area), y = surface_area, fill = sample_id)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title = "Corrected Surface Area Distribution by Tree Sample",
    x = "Sample ID (ordered by surface area)",
    y = "Surface Area (cm²)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold")
  )

# 3. Volume Distribution by Sample
p3 <- ggplot(plot_data, aes(x = reorder(sample_id, volume), y = volume, fill = sample_id)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title = "Corrected Volume Distribution by Tree Sample",
    x = "Sample ID (ordered by volume)",
    y = "Volume (cm³)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold")
  )

# 4. Combined scatter plot showing relationships
p4 <- ggplot(plot_data, aes(x = stem_diameter, y = surface_area, color = sample_id)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_color_viridis_d(name = "Sample ID") +
  labs(
    title = "Surface Area vs Stem Diameter by Sample",
    x = "Stem Diameter (cm)",
    y = "Corrected Surface Area (cm²)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")  # Too many samples for legend

# 5. Volume vs Surface Area relationship
p5 <- ggplot(plot_data, aes(x = surface_area, y = volume, color = sample_id)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_color_viridis_d(name = "Sample ID") +
  labs(
    title = "Volume vs Surface Area by Sample",
    x = "Corrected Surface Area (cm²)",
    y = "Corrected Volume (cm³)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 6. Coefficient of variation plot (shows consistency within each sample)
cv_data <- plot_data %>%
  group_by(sample_id) %>%
  summarise(
    cv_diameter = sd(stem_diameter) / mean(stem_diameter) * 100,
    cv_surface_area = sd(surface_area) / mean(surface_area) * 100,
    cv_volume = sd(volume) / mean(volume) * 100,
    n_measurements = n(),
    .groups = 'drop'
  ) %>%
  filter(n_measurements > 1) %>%  # Only samples with multiple measurements
  pivot_longer(cols = starts_with("cv_"), names_to = "metric", values_to = "cv") %>%
  mutate(metric = case_when(
    metric == "cv_diameter" ~ "Diameter",
    metric == "cv_surface_area" ~ "Surface Area", 
    metric == "cv_volume" ~ "Volume"
  ))

p6 <- ggplot(cv_data, aes(x = reorder(sample_id, cv), y = cv, fill = metric)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_viridis_d(name = "Metric") +
  labs(
    title = "Coefficient of Variation by Sample (Multiple Measurements Only)",
    x = "Sample ID",
    y = "Coefficient of Variation (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold")
  )

# Display plots
cat("\nGenerating plots...\n")

# Print individual plots
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

# Create a combined overview plot
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, 
                              top = "Tree Chamber Geometry Distributions by Sample ID")

# Save plots
ggsave("../../../outputs/figures/stem_diameter_by_sample.png", p1, width = 12, height = 8, dpi = 300)
ggsave("../../../outputs/figures/surface_area_by_sample.png", p2, width = 12, height = 8, dpi = 300)
ggsave("../../../outputs/figures/volume_by_sample.png", p3, width = 12, height = 8, dpi = 300)
ggsave("../../../outputs/figures/diameter_vs_surface_area.png", p4, width = 10, height = 8, dpi = 300)
ggsave("../../../outputs/figures/volume_vs_surface_area.png", p5, width = 10, height = 8, dpi = 300)
ggsave("../../../outputs/figures/coefficient_variation.png", p6, width = 12, height = 8, dpi = 300)

cat("Plots saved as PNG files in working directory\n")

# Create summary table for export
write.csv(sample_summary, "sample_summary_statistics.csv", row.names = FALSE)
cat("Summary statistics saved as 'sample_summary_statistics.csv'\n")

# Additional analysis: Identify outliers
outlier_analysis <- plot_data %>%
  group_by(sample_id) %>%
  mutate(
    diameter_z = abs(scale(stem_diameter)),
    surface_z = abs(scale(surface_area)),
    volume_z = abs(scale(volume))
  ) %>%
  filter(diameter_z > 2 | surface_z > 2 | volume_z > 2) %>%
  arrange(desc(pmax(diameter_z, surface_z, volume_z)))

if (nrow(outlier_analysis) > 0) {
  cat("\nPotential outliers (>2 standard deviations):\n")
  print(outlier_analysis %>% 
          select(sample_id, stem_diameter, surface_area, volume, diameter_z, surface_z, volume_z))
} else {
  cat("\nNo significant outliers detected.\n")
}

cat("\nAnalysis complete!")






# Chamber-Stem Combination Consistency Check
# Verify that identical physical setups produce consistent calculations

library(readxl)
library(dplyr)
library(ggplot2)
library(knitr)

# Read the data
data <- read_excel("combined_corrected_results.xlsx", sheet = "all_data")

# Create unique combination identifier
consistency_data <- data %>%
  filter(!is.na(Sample) & Sample != "") %>%
  mutate(
    stem_diameter = as.numeric(`D stem (cm)`),
    chamber_length = as.numeric(`Chamber length (cm)`),
    chamber_height = as.numeric(`Chamber height (cm)`),
    chamber_thickness = as.numeric(`T (sleeve thickness, cm)`),
    chamber_used = `Chamber used`,
    surface_area = as.numeric(Sc_cm2_corrected),
    volume = as.numeric(Vc_cm3_corrected),
    k_factor = as.numeric(K_eff),
    sample_id = Sample
  ) %>%
  filter(!is.na(stem_diameter) & !is.na(surface_area) & !is.na(volume)) %>%
  # Create combination key: Chamber + Stem + Physical dimensions
  mutate(
    combo_key = paste(
      chamber_used, 
      sample_id,
      round(stem_diameter, 1),
      round(chamber_length, 1), 
      round(chamber_height, 1),
      round(chamber_thickness, 1),
      sep = "_"
    )
  ) %>%
  select(combo_key, chamber_used, sample_id, stem_diameter, chamber_length, 
         chamber_height, chamber_thickness, surface_area, volume, k_factor)

# Find combinations with multiple measurements
combo_analysis <- consistency_data %>%
  group_by(combo_key) %>%
  summarise(
    n_measurements = n(),
    chamber_used = first(chamber_used),
    sample_id = first(sample_id),
    stem_diameter = first(stem_diameter),
    chamber_length = first(chamber_length),
    chamber_height = first(chamber_height),
    chamber_thickness = first(chamber_thickness),
    
    # Surface area consistency
    mean_surface_area = mean(surface_area),
    sd_surface_area = sd(surface_area),
    cv_surface_area = ifelse(n_measurements > 1, sd(surface_area) / mean(surface_area) * 100, 0),
    min_surface_area = min(surface_area),
    max_surface_area = max(surface_area),
    range_surface_area = max(surface_area) - min(surface_area),
    
    # Volume consistency  
    mean_volume = mean(volume),
    sd_volume = sd(volume),
    cv_volume = ifelse(n_measurements > 1, sd(volume) / mean(volume) * 100, 0),
    min_volume = min(volume),
    max_volume = max(volume),
    range_volume = max(volume) - min(volume),
    
    # K factor consistency
    mean_k_factor = mean(k_factor),
    sd_k_factor = sd(k_factor),
    cv_k_factor = ifelse(n_measurements > 1, sd(k_factor) / mean(k_factor) * 100, 0),
    
    .groups = 'drop'
  ) %>%
  arrange(desc(n_measurements), desc(cv_surface_area))

# Summary statistics
cat("CHAMBER-STEM COMBINATION CONSISTENCY ANALYSIS\n")
cat("=" , rep("=", 50), "\n", sep = "")
cat("Total unique combinations:", nrow(combo_analysis), "\n")
cat("Combinations with multiple measurements:", sum(combo_analysis$n_measurements > 1), "\n")
cat("Total measurements:", sum(combo_analysis$n_measurements), "\n")

# Focus on combinations with multiple measurements
multiple_measurements <- combo_analysis %>%
  filter(n_measurements > 1) %>%
  arrange(desc(cv_surface_area))

if (nrow(multiple_measurements) > 0) {
  cat("\nCOMBINATIONS WITH MULTIPLE MEASUREMENTS:\n")
  cat("Number of repeated combinations:", nrow(multiple_measurements), "\n")
  cat("Range of measurements per combination:", min(multiple_measurements$n_measurements), 
      "to", max(multiple_measurements$n_measurements), "\n\n")
  
  # Identify problematic combinations (high coefficient of variation)
  problematic <- multiple_measurements %>%
    filter(cv_surface_area > 1 | cv_volume > 1)  # >1% CV suggests inconsistency
  
  cat("CONSISTENCY ASSESSMENT:\n")
  if (nrow(problematic) > 0) {
    cat("⚠️  INCONSISTENT combinations (CV > 1%):", nrow(problematic), "\n")
    cat("✅ CONSISTENT combinations (CV ≤ 1%):", nrow(multiple_measurements) - nrow(problematic), "\n\n")
    
    cat("MOST INCONSISTENT COMBINATIONS:\n")
    print(problematic %>%
            select(chamber_used, sample_id, stem_diameter, n_measurements, 
                   cv_surface_area, cv_volume, range_surface_area, range_volume) %>%
            head(10))
    
  } else {
    cat("✅ ALL combinations are consistent (CV ≤ 1%)\n")
  }
  
  # Detailed table of all repeated combinations
  cat("\n\nDETAILED CONSISTENCY TABLE:\n")
  detailed_table <- multiple_measurements %>%
    mutate(
      surface_consistency = case_when(
        cv_surface_area <= 0.1 ~ "Excellent (≤0.1%)",
        cv_surface_area <= 1.0 ~ "Good (≤1%)", 
        cv_surface_area <= 5.0 ~ "Moderate (≤5%)",
        TRUE ~ "Poor (>5%)"
      ),
      volume_consistency = case_when(
        cv_volume <= 0.1 ~ "Excellent (≤0.1%)",
        cv_volume <= 1.0 ~ "Good (≤1%)",
        cv_volume <= 5.0 ~ "Moderate (≤5%)", 
        TRUE ~ "Poor (>5%)"
      )
    ) %>%
    select(chamber_used, sample_id, stem_diameter, n_measurements,
           mean_surface_area, cv_surface_area, surface_consistency,
           mean_volume, cv_volume, volume_consistency)
  
  print(detailed_table)
  
  # Visualization of consistency
  cat("\nGenerating consistency plots...\n")
  
  # Plot 1: CV distribution
  p1 <- ggplot(multiple_measurements, aes(x = cv_surface_area)) +
    geom_histogram(bins = 20, fill = "skyblue", alpha = 0.7, color = "black") +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = "Distribution of Surface Area Coefficient of Variation",
      subtitle = "Red line at 1% - values above suggest inconsistency", 
      x = "Coefficient of Variation (%)",
      y = "Number of Combinations"
    ) +
    theme_minimal()
  
  # Plot 2: Surface area range vs mean
  p2 <- ggplot(multiple_measurements, aes(x = mean_surface_area, y = range_surface_area)) +
    geom_point(aes(size = n_measurements, color = cv_surface_area), alpha = 0.7) +
    scale_color_gradient(low = "green", high = "red", name = "CV (%)") +
    scale_size_continuous(name = "N Measurements") +
    labs(
      title = "Surface Area Consistency Check",
      x = "Mean Surface Area (cm²)",
      y = "Range of Surface Areas (cm²)",
      subtitle = "Points near y=0 indicate consistent calculations"
    ) +
    theme_minimal()
  
  # Plot 3: Volume consistency
  p3 <- ggplot(multiple_measurements, aes(x = mean_volume, y = range_volume)) +
    geom_point(aes(size = n_measurements, color = cv_volume), alpha = 0.7) +
    scale_color_gradient(low = "green", high = "red", name = "CV (%)") +
    scale_size_continuous(name = "N Measurements") +
    labs(
      title = "Volume Consistency Check", 
      x = "Mean Volume (cm³)",
      y = "Range of Volumes (cm³)",
      subtitle = "Points near y=0 indicate consistent calculations"
    ) +
    theme_minimal()
  
  print(p1)
  print(p2) 
  print(p3)
  
  # Save outputs
  write.csv(detailed_table, "chamber_stem_consistency_check.csv", row.names = FALSE)
  ggsave("../../../outputs/figures/surface_area_cv_distribution.png", p1, width = 10, height = 6, dpi = 300)
  ggsave("../../../outputs/figures/surface_area_consistency_plot.png", p2, width = 10, height = 6, dpi = 300)
  ggsave("../../../outputs/figures/volume_consistency_plot.png", p3, width = 10, height = 6, dpi = 300)
  
} else {
  cat("\n⚠️  NO REPEATED COMBINATIONS FOUND\n")
  cat("Each chamber-stem combination appears only once in the dataset.\n")
  cat("This could mean:\n")
  cat("- Each measurement was on a unique setup\n") 
  cat("- Chamber/sample IDs are not being recorded consistently\n")
  cat("- Physical dimensions vary slightly between 'repeat' measurements\n")
}

# Show some example raw data to help diagnose
cat("\n\nSAMPLE OF RAW COMBINATION KEYS:\n")
sample_combos <- consistency_data %>%
  select(combo_key, chamber_used, sample_id, stem_diameter, chamber_length, surface_area, volume) %>%
  head(20)
print(sample_combos)

# Check for near-duplicates (allowing small rounding differences)
cat("\n\nCHECKING FOR NEAR-DUPLICATE COMBINATIONS:\n")
near_duplicate_check <- consistency_data %>%
  mutate(
    stem_diameter_rounded = round(stem_diameter, 0),  # Round to nearest cm
    chamber_length_rounded = round(chamber_length, 0),
    combo_key_loose = paste(chamber_used, sample_id, stem_diameter_rounded, 
                            chamber_length_rounded, sep = "_")
  ) %>%
  group_by(combo_key_loose) %>%
  summarise(
    n_measurements = n(),
    chamber_used = first(chamber_used),
    sample_id = first(sample_id),
    mean_surface_area = mean(surface_area),
    cv_surface_area = ifelse(n() > 1, sd(surface_area) / mean(surface_area) * 100, 0),
    .groups = 'drop'
  ) %>%
  filter(n_measurements > 1) %>%
  arrange(desc(cv_surface_area))

if (nrow(near_duplicate_check) > 0) {
  cat("Found", nrow(near_duplicate_check), "near-duplicate combinations (allowing ±0.5cm rounding):\n")
  print(near_duplicate_check)
} else {
  cat("No near-duplicate combinations found even with loose matching.\n")
}

cat("\nConsistency check complete!\n")
