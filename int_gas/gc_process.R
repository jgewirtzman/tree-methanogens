# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(broom)
library(readxl)
library(gridExtra)
library(viridis)

# Load the data files
raw_data <- read_excel('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/int_gas/Internal Concentration.xlsx', sheet = "Raw Data")
standards_conc <- read_excel('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/int_gas/Internal Concentration.xlsx', sheet = "Standard Concentrations")

# Handle redo samples - add this after loading raw_data but before creating GC_data

# First, remove the original SH1 and WA104 samples (keep only the redo versions)
raw_data <- raw_data %>%
  filter(!(`Tree ID` == "SH1" & !grepl("Redo|redo", `Tree ID`, ignore.case = TRUE))) %>%
  filter(!(`Tree ID` == "WA104" & !grepl("Redo|redo", `Tree ID`, ignore.case = TRUE)))

# Then rename the redo samples to remove the "_Redo" suffix
raw_data <- raw_data %>%
  mutate(`Tree ID` = case_when(
    `Tree ID` == "SH1_Redo" ~ "SH1",
    `Tree ID` == "WA104_Redo" ~ "WA104", 
    TRUE ~ `Tree ID`
  ))

# Clean column names for easier processing
colnames(raw_data) <- make.names(colnames(raw_data))

# Create a comprehensive dataset by merging FID and ECD data
GC_data <- raw_data

# Handle N2 samples - set O2 area to NA for N2 samples (since N2 has no O2)
GC_data$O2.Area[which(GC_data$Species.ID == "N2")] <- NA

# Convert relevant columns to numeric to avoid coercion issues
GC_data$CH4.Area <- as.numeric(GC_data$CH4.Area)
GC_data$CO2.Area <- as.numeric(GC_data$CO2.Area)
GC_data$O2.Area <- as.numeric(GC_data$O2.Area)
GC_data$N2O.Area <- as.numeric(GC_data$N2O.Area)

# Extract standards data
ghg_standard_names <- c("N2", "Outdoor Air 1", "Outdoor Air 2", "Outdoor Air 3", 
                        "Outdoor Air 4", "Outdoor Air 5", "Outdoor Air 6", "Outdoor Air 7",
                        "SB1", "SB2", "SB3", "SB4", "SB5", "SB6")

o2_standard_names <- c("N2", "Oxygen Standard 1", "Oxygen Standard 2", "Oxygen Standard 3", 
                       "Oxygen Standard 4", "Oxygen Standard 5", "SB1", "SB2", "SB3", "SB4", "SB5", "SB6")

# Create standards datasets by filtering and merging with concentration data
standards_conc_clean <- standards_conc %>%
  select(Sample, `[CO2] (ppm)`, `[CH4] (ppm)`, `[N2O] (ppm)`, `[O2] (ppm)`) %>%
  rename(
    CO2_ppm = `[CO2] (ppm)`,
    CH4_ppm = `[CH4] (ppm)`,
    N2O_ppm = `[N2O] (ppm)`,
    O2_ppm = `[O2] (ppm)`
  )

# Filter GC data for standards
GHG_standards_data <- GC_data %>%
  filter(Species.ID %in% ghg_standard_names) %>%
  left_join(standards_conc_clean, by = c("Species.ID" = "Sample"))

O2_standards_data <- GC_data %>%
  filter(Species.ID %in% o2_standard_names) %>%
  left_join(standards_conc_clean, by = c("Species.ID" = "Sample"))

# Convert concentration columns to numeric
GHG_standards_data$CO2_ppm <- as.numeric(GHG_standards_data$CO2_ppm)
GHG_standards_data$CH4_ppm <- as.numeric(GHG_standards_data$CH4_ppm)
GHG_standards_data$N2O_ppm <- as.numeric(GHG_standards_data$N2O_ppm)
O2_standards_data$O2_ppm <- as.numeric(O2_standards_data$O2_ppm)

# Filter out any NAs or missing values for standards
GHG_standards_data <- GHG_standards_data %>%
  filter(!is.na(CH4.Area) & !is.na(CH4_ppm)) %>%
  filter(!is.na(CO2.Area) & !is.na(CO2_ppm)) %>%
  filter(!is.na(N2O.Area) & !is.na(N2O_ppm))

O2_standards_data <- O2_standards_data %>%
  filter(!is.na(O2.Area) & !is.na(O2_ppm))

# Create linear and polynomial calibration curves (using all standards data)
# Linear curves
CO2_linear <- lm(CO2_ppm ~ CO2.Area, data = GHG_standards_data)
CH4_linear <- lm(CH4_ppm ~ CH4.Area, data = GHG_standards_data)
N2O_linear <- lm(N2O_ppm ~ N2O.Area, data = GHG_standards_data)
O2_linear <- lm(O2_ppm ~ O2.Area, data = O2_standards_data)

# Polynomial curves (2nd degree)
CO2_poly <- lm(CO2_ppm ~ poly(CO2.Area, 2), data = GHG_standards_data)
CH4_poly <- lm(CH4_ppm ~ poly(CH4.Area, 2), data = GHG_standards_data)
N2O_poly <- lm(N2O_ppm ~ poly(N2O.Area, 2), data = GHG_standards_data)
O2_poly <- lm(O2_ppm ~ poly(O2.Area, 2), data = O2_standards_data)

# Function to apply concentrations using the better fitting curve
apply_best_curve <- function(areas, linear_curve, poly_curve) {
  linear_r2 <- summary(linear_curve)$r.squared
  poly_r2 <- summary(poly_curve)$r.squared
  
  if (poly_r2 > linear_r2) {
    return(predict(poly_curve, newdata = data.frame(areas)))
  } else {
    return(predict(linear_curve, newdata = data.frame(areas)))
  }
}

# Apply concentration calculations using the best fitting curves
GC_data <- GC_data %>%
  mutate(
    CO2_concentration = case_when(
      is.na(CO2.Area) ~ NA_real_,
      TRUE ~ if (summary(CO2_poly)$r.squared > summary(CO2_linear)$r.squared) {
        predict(CO2_poly, newdata = data.frame(CO2.Area = CO2.Area))
      } else {
        predict(CO2_linear, newdata = data.frame(CO2.Area = CO2.Area))
      }
    ),
    
    CH4_concentration = case_when(
      is.na(CH4.Area) ~ NA_real_,
      TRUE ~ if (summary(CH4_poly)$r.squared > summary(CH4_linear)$r.squared) {
        predict(CH4_poly, newdata = data.frame(CH4.Area = CH4.Area))
      } else {
        predict(CH4_linear, newdata = data.frame(CH4.Area = CH4.Area))
      }
    ),
    
    N2O_concentration = case_when(
      is.na(N2O.Area) ~ NA_real_,
      TRUE ~ if (summary(N2O_poly)$r.squared > summary(N2O_linear)$r.squared) {
        predict(N2O_poly, newdata = data.frame(N2O.Area = N2O.Area))
      } else {
        predict(N2O_linear, newdata = data.frame(N2O.Area = N2O.Area))
      }
    ),
    
    O2_concentration = case_when(
      is.na(O2.Area) ~ NA_real_,
      TRUE ~ if (summary(O2_poly)$r.squared > summary(O2_linear)$r.squared) {
        predict(O2_poly, newdata = data.frame(O2.Area = O2.Area))
      } else {
        predict(O2_linear, newdata = data.frame(O2.Area = O2.Area))
      }
    )
  )

# Function to plot linear vs polynomial comparison
plot_curve_comparison <- function(standards_data, linear_curve, poly_curve, analyte, x_col, y_col) {
  
  linear_summary <- summary(linear_curve)
  poly_summary <- summary(poly_curve)
  
  linear_r2 <- round(linear_summary$r.squared, 4)
  poly_r2 <- round(poly_summary$r.squared, 4)
  
  # Generate smooth prediction lines
  area_range <- seq(min(standards_data[[x_col]], na.rm = TRUE), 
                    max(standards_data[[x_col]], na.rm = TRUE), 
                    length.out = 100)
  
  linear_pred <- predict(linear_curve, newdata = setNames(data.frame(area_range), x_col))
  poly_pred <- predict(poly_curve, newdata = setNames(data.frame(area_range), x_col))
  
  pred_data <- data.frame(
    area = area_range,
    linear = linear_pred,
    poly = poly_pred
  )
  
  ggplot(standards_data, aes_string(x = x_col, y = y_col)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line(data = pred_data, aes(x = area, y = linear), color = "blue", size = 1) +
    geom_line(data = pred_data, aes(x = area, y = poly), color = "red", size = 1) +
    annotate("text", x = Inf, y = Inf, 
             label = paste("Linear R² =", linear_r2, "\nPolynomial R² =", poly_r2), 
             hjust = 1.1, vjust = 1.1, size = 4, 
             color = ifelse(poly_r2 > linear_r2, "red", "blue")) +
    labs(title = paste(analyte, "Calibration: Linear vs Polynomial"),
         x = "Peak Area",
         y = "Concentration (ppm)",
         subtitle = paste("Best fit:", ifelse(poly_r2 > linear_r2, "Polynomial", "Linear"))) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "bottom"
    )
}

# Create calibration comparison plots
CO2_comparison <- plot_curve_comparison(GHG_standards_data, CO2_linear, CO2_poly, "CO2", "CO2.Area", "CO2_ppm")
CH4_comparison <- plot_curve_comparison(GHG_standards_data, CH4_linear, CH4_poly, "CH4", "CH4.Area", "CH4_ppm")
N2O_comparison <- plot_curve_comparison(GHG_standards_data, N2O_linear, N2O_poly, "N2O", "N2O.Area", "N2O_ppm")
O2_comparison <- plot_curve_comparison(O2_standards_data, O2_linear, O2_poly, "O2", "O2.Area", "O2_ppm")

# Display calibration plots
print(CO2_comparison)
print(CH4_comparison)
print(N2O_comparison)
print(O2_comparison)

# Print curve statistics
cat("=== Calibration Curve Comparison ===\n")
cat("\nCO2:\n")
cat("Linear R²:", round(summary(CO2_linear)$r.squared, 4), "\n")
cat("Polynomial R²:", round(summary(CO2_poly)$r.squared, 4), "\n")
cat("Best fit:", ifelse(summary(CO2_poly)$r.squared > summary(CO2_linear)$r.squared, "Polynomial", "Linear"), "\n")

cat("\nCH4:\n")
cat("Linear R²:", round(summary(CH4_linear)$r.squared, 4), "\n")
cat("Polynomial R²:", round(summary(CH4_poly)$r.squared, 4), "\n")
cat("Best fit:", ifelse(summary(CH4_poly)$r.squared > summary(CH4_linear)$r.squared, "Polynomial", "Linear"), "\n")

cat("\nN2O:\n")
cat("Linear R²:", round(summary(N2O_linear)$r.squared, 4), "\n")
cat("Polynomial R²:", round(summary(N2O_poly)$r.squared, 4), "\n")
cat("Best fit:", ifelse(summary(N2O_poly)$r.squared > summary(N2O_linear)$r.squared, "Polynomial", "Linear"), "\n")

cat("\nO2:\n")
cat("Linear R²:", round(summary(O2_linear)$r.squared, 4), "\n")
cat("Polynomial R²:", round(summary(O2_poly)$r.squared, 4), "\n")
cat("Best fit:", ifelse(summary(O2_poly)$r.squared > summary(O2_linear)$r.squared, "Polynomial", "Linear"), "\n")

# Check for unrealistic concentrations (>100% = 1,000,000 ppm)
cat("\n=== Checking for Unrealistic Concentrations ===\n")
high_co2 <- sum(GC_data$CO2_concentration > 1000000, na.rm = TRUE)
high_ch4 <- sum(GC_data$CH4_concentration > 1000000, na.rm = TRUE)
high_n2o <- sum(GC_data$N2O_concentration > 1000000, na.rm = TRUE)
high_o2 <- sum(GC_data$O2_concentration > 1000000, na.rm = TRUE)

cat("Samples with CO2 > 1,000,000 ppm (>100%):", high_co2, "\n")
cat("Samples with CH4 > 1,000,000 ppm (>100%):", high_ch4, "\n")
cat("Samples with N2O > 1,000,000 ppm (>100%):", high_n2o, "\n")
cat("Samples with O2 > 1,000,000 ppm (>100%):", high_o2, "\n")

# Show peak area ranges for standards vs samples
cat("\n=== Peak Area Ranges ===\n")
cat("CO2 Standards - Min Area:", min(GHG_standards_data$CO2.Area, na.rm = TRUE), 
    "Max Area:", max(GHG_standards_data$CO2.Area, na.rm = TRUE), "\n")
cat("CH4 Standards - Min Area:", min(GHG_standards_data$CH4.Area, na.rm = TRUE), 
    "Max Area:", max(GHG_standards_data$CH4.Area, na.rm = TRUE), "\n")
cat("N2O Standards - Min Area:", min(GHG_standards_data$N2O.Area, na.rm = TRUE), 
    "Max Area:", max(GHG_standards_data$N2O.Area, na.rm = TRUE), "\n")
cat("O2 Standards - Min Area:", min(O2_standards_data$O2.Area, na.rm = TRUE), 
    "Max Area:", max(O2_standards_data$O2.Area, na.rm = TRUE), "\n")

cat("\nAll Samples - Peak Area Ranges:\n")
cat("CO2 - Min Area:", min(GC_data$CO2.Area, na.rm = TRUE), 
    "Max Area:", max(GC_data$CO2.Area, na.rm = TRUE), "\n")
cat("CH4 - Min Area:", min(GC_data$CH4.Area, na.rm = TRUE), 
    "Max Area:", max(GC_data$CH4.Area, na.rm = TRUE), "\n")
cat("N2O - Min Area:", min(GC_data$N2O.Area, na.rm = TRUE), 
    "Max Area:", max(GC_data$N2O.Area, na.rm = TRUE), "\n")
cat("O2 - Min Area:", min(GC_data$O2.Area, na.rm = TRUE), 
    "Max Area:", max(GC_data$O2.Area, na.rm = TRUE), "\n")

# Set negative concentrations to 0 (negative concentrations are not physically possible)
GC_data <- GC_data %>%
  mutate(
    CO2_concentration = pmax(0, CO2_concentration, na.rm = TRUE),
    CH4_concentration = pmax(0, CH4_concentration, na.rm = TRUE),
    N2O_concentration = pmax(0, N2O_concentration, na.rm = TRUE),
    O2_concentration = pmax(0, O2_concentration, na.rm = TRUE)
  )

# Filter for samples only (exclude standards, blanks, check standards, etc.)
excluded_patterns <- c("Lab Air Blank", "N2", "Outdoor Air 1", "Outdoor Air 2", "Outdoor Air 3", 
                       "Outdoor Air 4", "Outdoor Air 5", "Outdoor Air 6", "Outdoor Air 7",
                       "Oxygen Standard", "BLANK", "Check Standard", "Ambient", "SB1", "SB2", 
                       "SB3", "SB4", "SB5", "SB6")

sample_data <- GC_data %>%
  filter(!Species.ID %in% excluded_patterns) %>%
  filter(!grepl("Std|Standard|Blank|blank|Ambient|SB|Outdoor Air|Oxygen", Species.ID, ignore.case = TRUE)) %>%
  filter(!is.na(CO2_concentration) | !is.na(CH4_concentration) | !is.na(N2O_concentration) | !is.na(O2_concentration))

# Create distribution plots for each gas (with log scale, adding 0.1 to handle zeros)
co2_dist <- ggplot(sample_data, aes(x = CO2_concentration + 0.1)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "CO2 Concentration Distribution",
       x = "CO2 (ppm) - Log Scale", y = "Frequency") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

ch4_dist <- ggplot(sample_data, aes(x = CH4_concentration + 0.1)) +
  geom_histogram(bins = 30, fill = "forestgreen", alpha = 0.7, color = "black") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "CH4 Concentration Distribution",
       x = "CH4 (ppm) - Log Scale", y = "Frequency") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

n2o_dist <- ggplot(sample_data, aes(x = N2O_concentration + 0.1)) +
  geom_histogram(bins = 30, fill = "orange", alpha = 0.7, color = "black") +
  scale_x_log10() +
  labs(title = "N2O Concentration Distribution",
       x = "N2O (ppm) - Log Scale", y = "Frequency") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

o2_dist <- ggplot(sample_data, aes(x = O2_concentration + 0.1)) +
  geom_histogram(bins = 30, fill = "red", alpha = 0.7, color = "black") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "O2 Concentration Distribution",
       x = "O2 (ppm) - Log Scale", y = "Frequency") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# Display distribution plots
print(co2_dist)
print(ch4_dist)
print(n2o_dist)
print(o2_dist)

# Create GHG vs O2 plots
co2_vs_o2 <- ggplot(sample_data, aes(x = O2_concentration, y = CO2_concentration)) +
  geom_point(aes(color = Species.ID), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_viridis_d(name = "Species") +
  labs(title = "CO2 vs O2 Concentration",
       x = "O2 (ppm)", y = "CO2 (ppm)") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "bottom"
  )

ch4_vs_o2 <- ggplot(sample_data, aes(x = O2_concentration, y = CH4_concentration)) +
  geom_point(aes(color = Species.ID), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_viridis_d(name = "Species") +
  labs(title = "CH4 vs O2 Concentration",
       x = "O2 (ppm)", y = "CH4 (ppm)") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "bottom"
  )

n2o_vs_o2 <- ggplot(sample_data, aes(x = O2_concentration, y = N2O_concentration)) +
  geom_point(aes(color = Species.ID), alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_viridis_d(name = "Species") +
  labs(title = "N2O vs O2 Concentration",
       x = "O2 (ppm)", y = "N2O (ppm)") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "bottom"
  )

# Display GHG vs O2 plots
print(co2_vs_o2)
print(ch4_vs_o2)
print(n2o_vs_o2)

# Export processed data
write.csv(GC_data, "processed_GC_data_internal_conc.csv", row.names = FALSE)
write.csv(sample_data, "sample_data_only.csv", row.names = FALSE)

# Summary statistics for samples only
sample_summary <- sample_data %>%
  group_by(Species.ID) %>%
  summarise(
    count = n(),
    mean_CO2 = round(mean(CO2_concentration, na.rm = TRUE), 2),
    mean_CH4 = round(mean(CH4_concentration, na.rm = TRUE), 2),
    mean_N2O = round(mean(N2O_concentration, na.rm = TRUE), 3),
    mean_O2 = round(mean(O2_concentration, na.rm = TRUE), 0),
    .groups = 'drop'
  ) %>%
  arrange(Species.ID)

cat("\n=== Sample Summary (Excluding Standards/Blanks) ===\n")
print(sample_summary)