# ==============================================================================
# Methanogen MANOVA
# ==============================================================================
# Purpose: Multivariate ANOVA testing methanogen abundance differences across
#   species and compartments.
#
# Pipeline stage: 3 — Analysis
# ==============================================================================

# Create the data frame with heights, flux_z, ch4_z, mcra_abundance_z, and tissue
df <- data.frame(
  heights_m = c(0.5, 0.5, 1.25, 1.25, 2.0, 2.0, 4.0, 4.0, 6.0, 6.0, 8.0, 8.0, 10.0, 10.0),
  flux_z = c(-0.68449456, -0.68449456, -0.36725165, -0.36725165, -0.37819106, -0.37819106, 
             2.31290399, 2.31290399, -0.44382752, -0.44382752, 0.01562773, 0.01562773, 
             -0.45476693, -0.45476693),
  ch4_z = c(-0.8733952, -0.8733952, -0.8638322, -0.8638322, 1.7960488, 1.7960488, 
            0.9423530, 0.9423530, 0.2353051, 0.2353051, -0.7875693, -0.7875693, 
            -0.4489101, -0.4489101),
  mcra_abundance_z = c(0.004228406, -0.793673912, -0.434127534, -0.793673912, 
                       -0.436476476, -0.793673912, 1.715163386, 0.670376174, 
                       1.798366697, 1.488047812, -0.793673912, -0.793673912, 
                       -0.793673912, -0.043534992),
  tissue = rep(c("Heartwood", "Sapwood"), 7)  # Alternating tissue types
)

# Add binned heights (categorical variable) as "Low", "Mid", "High"
df$group <- cut(df$heights_m, breaks = c(0, 2, 6, Inf), labels = c("Low", "Mid", "High"))

# Add quadratic term for heights
df$heights_m2 <- df$heights_m^2

# Display the data frame
df


# Prepare the response matrix
response_matrix <- as.matrix(df[, c("flux_z", "ch4_z", "mcra_abundance_z")])

# MANOVA with binned heights
manova_binned <- manova(response_matrix ~ group, data = df)

# Summary of MANOVA using multiple test statistics
summary(manova_binned, test = "Wilks")
summary(manova_binned, test = "Pillai")
summary(manova_binned, test = "Hotelling-Lawley")
summary(manova_binned, test = "Roy")

# Add the quadratic term
df$heights_m2 <- df$heights_m^2

# MANOVA with continuous height and height²
manova_continuous <- manova(response_matrix ~ heights_m + heights_m2 + tissue, data = df)

# Summary of MANOVA using multiple test statistics
summary(manova_continuous, test = "Wilks")
summary(manova_continuous, test = "Pillai")
summary(manova_continuous, test = "Hotelling-Lawley")
summary(manova_continuous, test = "Roy")


# Load necessary libraries
library(heplots)

# Perform MANOVA on the data with group and tissue as predictors
manova_model <- manova(cbind(flux_z, ch4_z, mcra_abundance_z) ~ group + tissue, data = df)

# Plot Canonical Variates to visualize the MANOVA results
plot(manova_model, type = "canonical", main = "Canonical Variate Plot for MANOVA")


# Install candisc if not already installed
if (!requireNamespace("candisc", quietly = TRUE)) {
  install.packages("candisc")
}
library(candisc)

# Run MANOVA with the group and tissue factors
manova_model <- manova(cbind(flux_z, ch4_z, mcra_abundance_z) ~ group + tissue, data = df)

# Perform Canonical Discriminant Analysis
cva <- candisc(manova_model)

# Plot the Canonical Variate Analysis results
plot(cva, main = "Canonical Variate Plot for MANOVA by Group and Tissue")
