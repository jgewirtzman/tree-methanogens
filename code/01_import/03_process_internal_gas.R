# ==============================================================================
# Process Internal Gas Concentrations
# ==============================================================================
# Purpose: Processes gas chromatograph data for internal tree gas concentrations
#   (CH4, CO2, N2O, O2) from extracted stem samples.
#
# Calibration (revised 2026-09-30; ported from tree-gas-traits
#   code/00a_process_ymf_gc.R). The previous version fitted one unweighted
#   quadratic through all standards (0-60,340 ppm CH4) and set negative
#   predictions to 0. That fit is dominated by the high standards, so
#   near-ambient samples predicted < 0 and 92/157 trees (and the lab-air
#   blanks) became CH4 = 0 ppm. Now:
#   - CH4: weighted linear fit through N2 + SB1-SB5; above SB5, linear
#     interpolation SB5 -> SB6 (SB6 is off-trend: response factor 6.0 vs
#     4.1-4.5 for SB1-SB5)
#   - CO2, N2O: weighted quadratics
#   - O2: weighted quadratic; NOT rescaled to the field "Ambient YYMMDD" vials
#     (they read ~19.6% O2 against ~20.7% for fresh outdoor air, a likely
#     storage effect; reported as QC only)
#   - no clamp and no substitution: detection limits are computed and values
#     are flagged (below LOD, near ambient), never overwritten. No tree sample
#     falls below the CH4 LOD.
#
# Pipeline stage: 01 Tree Data Processing
# Run after: None
#
# Inputs:
#   - Internal Concentration.xlsx (from data/raw/internal_gas/)
#
# Outputs:
#   - sample_data_only.csv
#   - processed_GC_data_internal_conc.csv
#   - internal_gas_calibration_check.csv (back-predicted standards)
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(broom)
library(readxl)
library(gridExtra)
library(viridis)

# Load the data files
raw_data <- read_excel('../../data/raw/internal_gas/Internal Concentration.xlsx', sheet = "Raw Data")
standards_conc <- read_excel('../../data/raw/internal_gas/Internal Concentration.xlsx', sheet = "Standard Concentrations")

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

# SB6 has no certified O2; the outdoor-air series anchors the ambient end
o2_standard_names <- c("N2", "Oxygen Standard 1", "Oxygen Standard 2", "Oxygen Standard 3", 
                       "Oxygen Standard 4", "Oxygen Standard 5",
                       "Outdoor Air 1", "Outdoor Air 2", "Outdoor Air 3", "Outdoor Air 4",
                       "Outdoor Air 5", "Outdoor Air 6", "Outdoor Air 7",
                       "SB1", "SB2", "SB3", "SB4", "SB5")

# Create standards datasets by filtering and merging with concentration data
standards_conc_clean <- standards_conc %>%
  select(Sample, `[CO2] (ppm)`, `[CH4] (ppm)`, `[N2O] (ppm)`, `[O2] (ppm)`) %>%
  rename(
    CO2_ppm = `[CO2] (ppm)`,
    CH4_ppm = `[CH4] (ppm)`,
    N2O_ppm = `[N2O] (ppm)`,
    O2_ppm = `[O2] (ppm)`
  ) %>%
  mutate(across(ends_with("_ppm"), as.numeric))

standards_data <- GC_data %>%
  filter(Species.ID %in% union(ghg_standard_names, o2_standard_names)) %>%
  left_join(standards_conc_clean, by = c("Species.ID" = "Sample"))

# ==============================================================================
# Calibration
# ==============================================================================
# Standards are pooled over the four run days (SB areas vary < 5% between
# days). Weights 1/(conc + c0)^2 make the fits minimise *relative* error, so
# ambient-level standards carry as much weight as the high ones.
#
# Standard choice (from response factors, ppm per area unit):
#  - CH4: SB1-SB5 give 4.1-4.5 (linear, ~zero intercept); SB6 gives 6.0, so a
#    single curve through SB6 biases everything below it by ~20%. The
#    outdoor-air dilutions assume a nominal 1.8 ppm stock, but undiluted
#    outdoor air reads ~2.3 ppm against the SB tanks, so they are not used for
#    CH4 or CO2. -> linear through N2 + SB1-SB5; above SB5, linear
#    interpolation SB5 -> SB6 (flagged).
#  - CO2: response factor rises steadily SB2 -> SB6 (2.6 -> 4.2), i.e. a
#    genuinely curved response -> weighted quadratic through N2 + SB1-SB6.
#  - N2O: response factor 0.0035-0.0041 for all standards; outdoor-air series
#    kept -> weighted quadratic.
#  - O2: N2 + O2 standards + outdoor-air series + SB1-SB5 -> weighted quadratic.

c0 <- c(CH4 = 0.3, CO2 = 50, N2O = 0.05, O2 = 20000)
sb_names <- paste0("SB", 1:6)

cal_data <- function(gas, std_names) {
  standards_data %>%
    filter(Species.ID %in% std_names) %>%
    transmute(Species.ID, area = .data[[paste0(gas, ".Area")]], conc = .data[[paste0(gas, "_ppm")]]) %>%
    filter(!is.na(area), !is.na(conc))
}
fit_weighted_quadratic <- function(d, gas) {
  lm(conc ~ area + I(area^2), data = d, weights = 1 / (conc + c0[[gas]])^2)
}

CH4_cal_data <- cal_data("CH4", c("N2", sb_names))
CH4_linear <- lm(conc ~ area, data = filter(CH4_cal_data, Species.ID != "SB6"),
                 weights = 1 / (conc + c0[["CH4"]])^2)
CH4_top <- CH4_cal_data %>%
  filter(Species.ID %in% c("SB5", "SB6")) %>%
  group_by(Species.ID) %>%
  summarise(area = mean(area), conc = mean(conc), .groups = "drop") %>%
  arrange(area)

cal_curves <- list(
  CO2 = fit_weighted_quadratic(cal_data("CO2", c("N2", sb_names)), "CO2"),
  N2O = fit_weighted_quadratic(cal_data("N2O", ghg_standard_names), "N2O"),
  O2  = fit_weighted_quadratic(cal_data("O2",  o2_standard_names), "O2")
)

predict_gas <- function(gas, area) {
  if (gas == "CH4") {
    lin <- as.numeric(predict(CH4_linear, newdata = data.frame(area = area)))
    # SB5 -> SB6 line, held at the SB6 value beyond its area (no extrapolation
    # past the highest standard). One tree, RO8, lies ~17% above the SB6 area:
    # it is reported as 60,340 ppm and flagged CH4_above_SB6 (a lower bound).
    # Matches tree-gas-traits code/00a_process_ymf_gc.R exactly.
    top <- approx(CH4_top$area, CH4_top$conc, xout = area, rule = 2)$y
    # above the SB5 area, interpolate towards SB6 (never below the linear fit)
    return(if_else(!is.na(area) & area > CH4_top$area[1], pmax(lin, top), lin))
  }
  as.numeric(predict(cal_curves[[gas]], newdata = data.frame(area = area)))
}

# Detection limit: 3 x SD of replicate predictions of the lowest certified
# standard (SB1; one injection per run day)
SB1_data <- standards_data %>% filter(Species.ID == "SB1")
lod <- sapply(c(CH4 = "CH4", CO2 = "CO2", N2O = "N2O"), function(gas) {
  3 * sd(predict_gas(gas, SB1_data[[paste0(gas, ".Area")]]), na.rm = TRUE)
})
SB5_CH4_ppm <- CH4_top$conc[CH4_top$Species.ID == "SB5"]

# Apply calibration to all injections (uncensored, un-normalised)
GC_data <- GC_data %>%
  mutate(
    CO2_concentration = predict_gas("CO2", CO2.Area),
    CH4_concentration = predict_gas("CH4", CH4.Area),
    N2O_concentration = predict_gas("N2O", N2O.Area),
    O2_concentration  = predict_gas("O2",  O2.Area)
  )

# Back-prediction of standards, including the independent check standard
# (an SB2-like blend: 3.008 ppm CH4, 426.2 ppm CO2, 0.379 ppm N2O, 240,100 ppm O2)
calibration_check <- GC_data %>%
  filter(Species.ID %in% c(ghg_standard_names, o2_standard_names, "Check Standard")) %>%
  select(Species.ID, FID.Date, CO2_concentration, CH4_concentration,
         N2O_concentration, O2_concentration) %>%
  pivot_longer(ends_with("_concentration"), names_to = "gas", values_to = "predicted") %>%
  mutate(gas = sub("_concentration", "", gas)) %>%
  left_join(standards_conc_clean %>%
              bind_rows(tibble(Sample = "Check Standard", CH4_ppm = 3.008, CO2_ppm = 426.2,
                               N2O_ppm = 0.379, O2_ppm = 240100)) %>%
              pivot_longer(ends_with("_ppm"), names_to = "gas", values_to = "certified") %>%
              mutate(gas = sub("_ppm", "", gas)),
            by = c("Species.ID" = "Sample", "gas")) %>%
  filter(!is.na(certified), !is.na(predicted),
         !(gas == "O2" & !Species.ID %in% c(o2_standard_names, "Check Standard"))) %>%
  mutate(rel_error = (predicted - certified) / certified)

cat("=== Calibration ===\n")
cat("Detection limits (ppm):\n"); print(round(lod, 3))
cat("\nBack-predicted standards (median relative error by level):\n")
calibration_check %>%
  filter(certified > 0) %>%
  group_by(gas, Species.ID, certified) %>%
  summarise(median_pred = median(predicted), rel_err = round(median(rel_error), 3),
            .groups = "drop") %>%
  arrange(gas, certified) %>%
  print(n = 80)

write.csv(calibration_check, "../../data/processed/internal_gas/internal_gas_calibration_check.csv",
          row.names = FALSE)

# Calibration plot: back-prediction error of each standard
calibration_plot <- calibration_check %>%
  filter(certified > 0) %>%
  ggplot(aes(certified, rel_error * 100, colour = Species.ID == "Check Standard")) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red"),
                      labels = c("Standards", "Check standard"), name = NULL) +
  facet_wrap(~ gas, scales = "free") +
  labs(x = "Certified concentration (ppm)", y = "Back-prediction error (%)",
       title = "Internal gas calibration") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "bottom")
print(calibration_plot)

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

# ==============================================================================
# Field-ambient vials: QC only
# ==============================================================================
# "Ambient YYMMDD" vials were filled in the field on each sampling day and
# stored with the tree samples until the GC run (~14 months). Their O2 reads
# ~19.6% against ~20.7% for fresh outdoor air in the same runs. We report this
# as QC and do not rescale O2 to them (decision shared with tree-gas-traits).

ambient_vials <- GC_data %>% filter(grepl("^Ambient", Species.ID))
ambient_O2 <- median(ambient_vials$O2_concentration, na.rm = TRUE)
cat(sprintf("\nField ambient vials (n=%d): O2 %.0f ppm, CO2 %.0f ppm, CH4 %.2f ppm\n",
            nrow(ambient_vials), ambient_O2,
            median(ambient_vials$CO2_concentration, na.rm = TRUE),
            median(ambient_vials$CH4_concentration, na.rm = TRUE)))

# ==============================================================================
# Detection flags (replace the old zero clamp)
# ==============================================================================
# Values are kept as calibrated; nothing is clamped or substituted. Flags mark
# below-LOD and near-ambient (< 3 ppm CH4) samples for sensitivity runs.
# Standards and blanks are processed identically so they can be inspected.

GC_data <- GC_data %>%
  mutate(
    CH4_below_lod    = CH4_concentration < lod[["CH4"]],
    CO2_below_lod    = CO2_concentration < lod[["CO2"]],
    N2O_below_lod    = N2O_concentration < lod[["N2O"]],
    CH4_near_ambient = CH4_concentration < 3,
    CH4_above_SB5    = CH4_concentration > SB5_CH4_ppm,
    CH4_above_SB6    = CH4.Area > CH4_top$area[2]
  )

cat("\nLab-air blanks (median ppm): CH4",
    round(median(GC_data$CH4_concentration[GC_data$Species.ID == "Lab Air Blank"], na.rm = TRUE), 2),
    "| lab CH4.ppm column",
    round(median(as.numeric(GC_data$CH4.ppm[GC_data$Species.ID == "Lab Air Blank"]), na.rm = TRUE), 2), "\n")

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
write.csv(GC_data, "../../data/processed/internal_gas/processed_GC_data_internal_conc.csv", row.names = FALSE)
write.csv(sample_data, "../../data/processed/internal_gas/sample_data_only.csv", row.names = FALSE)

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
cat("\nTree samples:", nrow(sample_data),
    "| CH4 < LOD:", sum(sample_data$CH4_below_lod, na.rm = TRUE),
    "| CH4 < 3 ppm:", sum(sample_data$CH4_near_ambient, na.rm = TRUE),
    "| CH4 above SB5 (", SB5_CH4_ppm, "ppm; SB5->SB6 interpolation):",
    sum(sample_data$CH4_above_SB5, na.rm = TRUE),
    "| above SB6 (held at SB6, lower bound):", sum(sample_data$CH4_above_SB6, na.rm = TRUE), "\n")
