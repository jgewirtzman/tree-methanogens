# ==============================================================================
# Quadratic Models for mcrA-Flux Relationships
# ==============================================================================
# Purpose: Quadratic models for mcrA-flux relationships.
#
# Pipeline stage: 3 — Analysis
# ==============================================================================

# Multivariate Quadratic Regression
mv_model <- lm(cbind(flux_z, ch4_z, mcra_abundance_z) ~ heights_m + I(heights_m^2), data = df)
summary(mv_model)


library(ggplot2)

ggplot(df, aes(x = heights_m, y = flux_z)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
  ggtitle("Flux_z vs. Height with Quadratic Fit")

# Repeat for ch4_z and mcra_abundance_z
# Melt the data frame for ggplot2
library(reshape2)
melted_data <- melt(your_data_frame, id.vars = "heights_m", 
                    measure.vars = c("flux_z", "ch4_z", "mcra_abundance_z"))

ggplot(melted_data, aes(x = heights_m, y = value, color = variable)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
  ggtitle("All Variables vs. Height with Quadratic Fits")



library(VGAM)

# Load necessary library
library(VGAM)

str(normalized_df)

# Fit a VGAM with multiple response variables modeled as functions of height and tissue
vgam_model <- vgam(cbind(flux_z, ch4_z, mcra_abundance_z) ~ s(height_cm) + tissue, family = gaussianff, data = normalized_df)

# Display summary of the VGAM model
summary(vgam_model)

plot(vgam_model)






# Load necessary libraries for plotting
library(ggplot2)
library(dplyr)

# Generate a new data frame for prediction over a range of heights and tissues
height_range <- seq(min(normalized_df$height_cm), max(normalized_df$height_cm), length.out = 100)
tissue_levels <- unique(normalized_df$tissue)
predict_df <- expand.grid(height_cm = height_range, tissue = tissue_levels)

# Get predictions from the VGAM model
predictions <- predict(vgam_model, newdata = predict_df, type = "response")

# Add predictions to the data frame
predict_df$flux_pred <- predictions[, 1]   # Predictions for flux_z
predict_df$ch4_pred <- predictions[, 2]    # Predictions for ch4_z
predict_df$mcra_pred <- predictions[, 3]   # Predictions for mcra_abundance_z

# Plot smooth terms for each response variable
# Flux
ggplot(predict_df, aes(x = height_cm, y = flux_pred, color = tissue)) +
  geom_line(size = 1) +
  labs(title = "Predicted Flux (Z-score) by Height and Tissue", x = "Height (cm)", y = "Predicted Flux (Z-score)") +
  theme_minimal()

# CH4
ggplot(predict_df, aes(x = height_cm, y = ch4_pred, color = tissue)) +
  geom_line(size = 1) +
  labs(title = "Predicted CH4 (Z-score) by Height and Tissue", x = "Height (cm)", y = "Predicted CH4 (Z-score)") +
  theme_minimal()

# mcrA abundance
ggplot(predict_df, aes(x = height_cm, y = mcra_pred, color = tissue)) +
  geom_line(size = 1) +
  labs(title = "Predicted mcrA Abundance (Z-score) by Height and Tissue", x = "Height (cm)", y = "Predicted mcrA Abundance (Z-score)") +
  theme_minimal()





# Fit VGAM as before, if not already done
vgam_model <- vgam(cbind(flux_z, ch4_z, mcra_abundance_z) ~ s(height_cm) + tissue, family = gaussianff, data = normalized_df)

# Extract residuals for each response
residuals_df <- as.data.frame(resid(vgam_model))

# Use spline correlogram to test for spatial autocorrelation in residuals
library(ncf)

# Spline correlogram for flux residuals
spline_corr_flux <- spline.correlog(x = normalized_df$height_cm, y = normalized_df$height_cm, z = residuals_df[, 1])
plot(spline_corr_flux, main = "Spatial Autocorrelation in Flux Residuals")

# Spline correlogram for CH4 residuals
spline_corr_ch4 <- spline.correlog(x = normalized_df$height_cm, y = normalized_df$height_cm, z = residuals_df[, 2])
plot(spline_corr_ch4, main = "Spatial Autocorrelation in CH4 Residuals")

# Spline correlogram for mcrA abundance residuals
spline_corr_mcra <- spline.correlog(x = normalized_df$height_cm, y = normalized_df$height_cm, z = residuals_df[, 3])
plot(spline_corr_mcra, main = "Spatial Autocorrelation in mcrA Abundance Residuals")

