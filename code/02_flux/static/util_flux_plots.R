# ==============================================================================
# Static Flux Exploratory Plots
# ==============================================================================
# Purpose: General exploratory visualization of static flux patterns.
# ==============================================================================

library(ggplot2)
library(dplyr)

final_dataset<-read.csv('../../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv')

# Prepare the data
plot_data <- final_dataset %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux) & !is.na(height_numeric))

# Simple boxplot with jittered points
ggplot(plot_data, aes(x = factor(height_numeric), y = CH4_best.flux)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Remove outliers from boxplot since we'll show all points
  geom_jitter(alpha = 0.5, width = 0.2, size = 1) +  # Add jittered points
  facet_wrap(~ species, scales = "free") +
  coord_flip() +  # Flip coordinates
  labs(
    x = "Measurement Height (cm)",
    y = "CH₄ Flux (best model)",
    title = "CH₄ Fluxes by Tree Height",
    subtitle = "Faceted by Species Code"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )