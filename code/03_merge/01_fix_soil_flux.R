# ==============================================================================
# Fix Soil Flux Calculations
# ==============================================================================
# Purpose: Recalculates CH4 fluxes for soil data with corrected chamber
#   geometry parameters.
#
# Pipeline stage: 02 Integration
# Run after: 01_flux_processing/semirigid/ pipeline
#
# Inputs:
#   - lgr_manual_identification_results_soil.csv
#   - Flux result CSVs (from semirigid processing)
#
# Outputs (data/processed/flux/):
#   - CH4_flux_lgr_results_soil.csv, CH4_best_flux_lgr_results_soil.csv (read by the
#     screens, detection and the soil flux table)
#   - semirigid_tree_final_complete_dataset_soil_CORRECTED.csv
# (Removed 2026-10-17: a later diagnostic section that wrote six files nothing reads --
#  CH4_best_flux_MANUAL_soil, FINAL_soil_dataset_CH4_corrected, CH4_best_flux_PROPER_soil,
#  FINAL_soil_flux_dataset_complete / _simplified / _summary_by_treatment.)
# ==============================================================================

library(goFlux)
source("../lib/flux_windows.R")   # maybe_flux_plot(), maybe_flux2pdf()
library(dplyr)
library(readr)

source("../lib/chamber_constants.R")  # ANALYZER_VOLUME_CM3

cat("=== FIXING CH4 FLUX CALCULATIONS FOR SOIL DATA ===\n\n")

# Define correct chamber geometry (same as your main script)
CORRECT_CHAMBER_SURFACE_AREA_CM2 <- 507.7
CORRECT_CHAMBER_VOLUME_L <- 7.53
correct_tubing_volume_cm3 <- pi * (1/16)^2 * 12 * 12 * 16.387
correct_system_volume_cm3 <- ANALYZER_VOLUME_CM3
correct_total_system_volume_cm3 <- CORRECT_CHAMBER_VOLUME_L * 1000 + 
  correct_tubing_volume_cm3 + correct_system_volume_cm3
CORRECT_VTOT_L <- correct_total_system_volume_cm3 / 1000

cat("Using chamber volume:", CORRECT_CHAMBER_VOLUME_L, "L\n")
cat("Using total system volume:", round(CORRECT_VTOT_L, 3), "L\n\n")

# Load the manual identification data with correct volume
manID.lgr3 <- read_csv(window_path("lgr_manual_identification_results_soil.csv"), 
                       show_col_types = FALSE)

# Ensure correct volume is set; chamber temperature from the current auxfile (station
# air in true time), not the value saved with the windows
soil_aux <- read_csv("../../data/processed/flux/auxfile_goFlux_soilflux_with_weather.csv", show_col_types = FALSE)
manID.lgr3 <- manID.lgr3 %>%
  mutate(
    Vtot = CORRECT_VTOT_L,
    Area = CORRECT_CHAMBER_SURFACE_AREA_CM2,
    Tcham = soil_aux$Tcham[match(UniqueID, soil_aux$UniqueID)]
  )
stopifnot(!anyNA(manID.lgr3$Tcham))

cat("Loaded manual ID data with", nrow(unique(manID.lgr3[c("UniqueID")])), "unique measurements\n")

# Verify CH4 data exists and looks correct
ch4_check <- manID.lgr3 %>%
  summarise(
    has_CH4 = "CH4dry_ppb" %in% names(.),
    n_CH4_data = sum(!is.na(CH4dry_ppb)),
    CH4_mean = mean(CH4dry_ppb, na.rm = TRUE),
    CH4_min = min(CH4dry_ppb, na.rm = TRUE),
    CH4_max = max(CH4dry_ppb, na.rm = TRUE)
  )

cat("\nCH4 data verification:\n")
cat("- CH4dry_ppb column exists:", ch4_check$has_CH4, "\n")
cat("- Data points with CH4:", ch4_check$n_CH4_data, "\n")
cat("- CH4 concentration mean:", round(ch4_check$CH4_mean, 1), "ppb\n")
cat("- CH4 concentration range:", round(ch4_check$CH4_min, 1), "-", 
    round(ch4_check$CH4_max, 1), "ppb\n")

# Also check CO2 for comparison
co2_check <- manID.lgr3 %>%
  summarise(
    CO2_mean = mean(CO2dry_ppm, na.rm = TRUE),
    CO2_min = min(CO2dry_ppm, na.rm = TRUE),
    CO2_max = max(CO2dry_ppm, na.rm = TRUE)
  )

cat("\nCO2 data for comparison:\n")
cat("- CO2 concentration mean:", round(co2_check$CO2_mean, 1), "ppm\n")
cat("- CO2 concentration range:", round(co2_check$CO2_min, 1), "-", 
    round(co2_check$CO2_max, 1), "ppm\n\n")

# =============================================================================
# RECALCULATE CH4 FLUXES WITH CORRECT GASTYPE
# =============================================================================

cat("=== CALCULATING CH4 FLUXES (USING CH4dry_ppb) ===\n")

# Calculate CH4 fluxes - CRITICAL: Use CH4dry_ppb, not CO2dry_ppm!
CH4_flux_lgr3 <- goFlux(
  dataframe = manID.lgr3,
  gastype = "CH4dry_ppb",  # THIS IS THE CRITICAL FIX!
  H2O_col = "H2O_ppm",
  warn.length = 60
)

cat("CH4 flux calculation complete:", nrow(CH4_flux_lgr3), "measurements\n")

# Check the results
ch4_flux_summary <- CH4_flux_lgr3 %>%
  summarise(
    n = n(),
    LM_flux_mean = mean(LM.flux, na.rm = TRUE),
    LM_flux_median = median(LM.flux, na.rm = TRUE),
    LM_flux_min = min(LM.flux, na.rm = TRUE),
    LM_flux_max = max(LM.flux, na.rm = TRUE),
    HM_flux_mean = mean(HM.flux, na.rm = TRUE)
  )

cat("\nCH4 Flux Results Summary:\n")
cat("- Linear model flux mean:", round(ch4_flux_summary$LM_flux_mean, 4), "nmol/m²/s\n")
cat("- Linear model flux median:", round(ch4_flux_summary$LM_flux_median, 4), "nmol/m²/s\n")
cat("- Linear model flux range:", round(ch4_flux_summary$LM_flux_min, 4), "to", 
    round(ch4_flux_summary$LM_flux_max, 4), "nmol/m²/s\n")

# For comparison, load and check CO2 fluxes
CO2_flux_lgr3 <- read_csv("../../data/processed/flux/CO2_flux_lgr_results_soil.csv", show_col_types = FALSE)
co2_flux_mean <- mean(CO2_flux_lgr3$LM.flux, na.rm = TRUE)

cat("\nComparison with CO2:\n")
cat("- CO2 flux mean:", round(co2_flux_mean, 2), "nmol/m²/s\n")
cat("- CH4 flux mean:", round(ch4_flux_summary$LM_flux_mean, 4), "nmol/m²/s\n")
cat("- Ratio CO2/CH4:", round(co2_flux_mean / abs(ch4_flux_summary$LM_flux_mean), 1), "\n")

if(abs(ch4_flux_summary$LM_flux_mean - co2_flux_mean) < 0.1) {
  cat("\n❌ WARNING: CH4 and CO2 fluxes are still identical!\n")
  cat("   Check that CH4dry_ppb column has different values from CO2dry_ppm\n")
} else {
  cat("\n✓ Success! CH4 and CO2 fluxes are now different.\n")
}

# Save the corrected CH4 flux results
write_csv(CH4_flux_lgr3, "../../data/processed/flux/CH4_flux_lgr_results_soil.csv")
cat("\nSaved corrected CH4 flux results to: data/processed/flux/CH4_flux_lgr_results_soil.csv\n")

# =============================================================================
# RUN BEST.FLUX ANALYSIS ON CORRECTED CH4 DATA
# =============================================================================

cat("\n=== RUNNING BEST.FLUX ON CORRECTED CH4 DATA ===\n")

CH4_best_lgr3 <- best.flux(
  flux.result = CH4_flux_lgr3,
  criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", 
               "MDF", "nb.obs", "intercept", "p-value"),
  intercept.lim = NULL,
  g.limit = 2,
  p.val = 0.05,
  k.ratio = 1,
  warn.length = 60
)

# Save corrected best flux results
write_csv(CH4_best_lgr3, "../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")
cat("Saved corrected CH4 best flux results\n")

# =============================================================================
# FINAL STATISTICS COMPARISON
# =============================================================================

cat("\n=== FINAL CORRECTED FLUX STATISTICS ===\n")

# Load CO2 best flux for comparison
CO2_best_lgr3 <- read_csv("../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv", 
                          show_col_types = FALSE)

# CO2 statistics
co2_stats <- CO2_best_lgr3 %>%
  summarise(
    gas = "CO2",
    n = n(),
    mean_flux = mean(flux.term, na.rm = TRUE),
    median_flux = median(flux.term, na.rm = TRUE),
    sd_flux = sd(flux.term, na.rm = TRUE),
    min_flux = min(flux.term, na.rm = TRUE),
    max_flux = max(flux.term, na.rm = TRUE)
  )

# CH4 statistics
ch4_stats <- CH4_best_lgr3 %>%
  summarise(
    gas = "CH4",
    n = n(),
    mean_flux = mean(flux.term, na.rm = TRUE),
    median_flux = median(flux.term, na.rm = TRUE),
    sd_flux = sd(flux.term, na.rm = TRUE),
    min_flux = min(flux.term, na.rm = TRUE),
    max_flux = max(flux.term, na.rm = TRUE)
  )

cat("\nCO2 Flux (nmol/m²/s):\n")
cat("- Mean:", round(co2_stats$mean_flux, 2), "\n")
cat("- Median:", round(co2_stats$median_flux, 2), "\n")
cat("- Range:", round(co2_stats$min_flux, 2), "to", round(co2_stats$max_flux, 2), "\n")

cat("\nCH4 Flux (nmol/m²/s):\n")
cat("- Mean:", round(ch4_stats$mean_flux, 4), "\n")
cat("- Median:", round(ch4_stats$median_flux, 4), "\n")
cat("- Range:", round(ch4_stats$min_flux, 4), "to", round(ch4_stats$max_flux, 4), "\n")

# =============================================================================
# UPDATE FINAL DATASET WITH CORRECTED CH4 VALUES
# =============================================================================

cat("\n=== UPDATING FINAL DATASET ===\n")

# Load original data
original_data <- read_csv('../../data/processed/flux/auxfile_goFlux_soilflux_with_weather_formatted.csv',
                          show_col_types = FALSE)

# Update volume
original_data <- original_data %>%
  mutate(
    Vtot = CORRECT_VTOT_L,
    Area = CORRECT_CHAMBER_SURFACE_AREA_CM2
  )

# Add CO2_ prefix to CO2 results
co2_results <- CO2_best_lgr3 %>%
  rename_with(~ paste0("CO2_", .), -UniqueID)

# Add CH4_ prefix to corrected CH4 results
ch4_results <- CH4_best_lgr3 %>%
  rename_with(~ paste0("CH4_", .), -UniqueID)

# Merge everything together
final_dataset <- original_data %>%
  left_join(co2_results, by = "UniqueID") %>%
  left_join(ch4_results, by = "UniqueID")

# Save the corrected final dataset
write_csv(final_dataset, "../../data/processed/flux/semirigid_tree_final_complete_dataset_soil_CORRECTED.csv")

cat("Updated final dataset saved with", nrow(final_dataset), "rows\n")
cat("File: data/processed/flux/semirigid_tree_final_complete_dataset_soil_CORRECTED.csv\n")

# =============================================================================
# CREATE CORRECTED PLOTS
# =============================================================================

cat("\n=== CREATING CORRECTED PLOTS ===\n")

# Fix any missing model values
CH4_best_lgr3 <- CH4_best_lgr3 %>%
  mutate(
    model = case_when(
      is.na(model) ~ "LM",
      TRUE ~ model
    )
  )

# Create CH4 plots with corrected data
CH4_plots_lgr3 <- maybe_flux_plot(
  flux.results = CH4_best_lgr3,
  dataframe = manID.lgr3,
  gastype = "CH4dry_ppb",  # Correct gastype for plots
  shoulder = 30,
  plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
  quality.check = TRUE,
  best.model = TRUE,
  p.val.disp = "round"
)

# Save corrected CH4 plots
maybe_flux2pdf(
  plot.list = CH4_plots_lgr3,
  outfile = "../../outputs/figures/CH4_flux_plots_soil_CORRECTED.pdf",
  width = 11.6,
  height = 8.2
)

cat("Corrected CH4 plots saved to: outputs/figures/CH4_flux_plots_soil_CORRECTED.pdf\n")

cat("\n=== CH4 FLUX CORRECTION COMPLETE ===\n")
cat("All CH4 fluxes have been recalculated using CH4dry_ppb data.\n")
cat("Files created/updated:\n")
cat("- data/processed/flux/CH4_flux_lgr_results_soil.csv\n")
cat("- data/processed/flux/CH4_best_flux_lgr_results_soil.csv\n")
cat("- data/processed/flux/semirigid_tree_final_complete_dataset_soil_CORRECTED.csv\n")
cat("- outputs/figures/CH4_flux_plots_soil_CORRECTED.pdf\n")
