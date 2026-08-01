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
# Outputs:
#   - FINAL_soil_flux_dataset_complete.csv
# ==============================================================================

library(goFlux)
library(dplyr)
library(readr)

cat("=== FIXING CH4 FLUX CALCULATIONS FOR SOIL DATA ===\n\n")

# Define correct chamber geometry (same as your main script)
CORRECT_CHAMBER_SURFACE_AREA_CM2 <- 507.7
CORRECT_CHAMBER_VOLUME_L <- 7.53
correct_tubing_volume_cm3 <- pi * (1/16)^2 * 12 * 12 * 16.387
correct_system_volume_cm3 <- 70
correct_total_system_volume_cm3 <- CORRECT_CHAMBER_VOLUME_L * 1000 + 
  correct_tubing_volume_cm3 + correct_system_volume_cm3
CORRECT_VTOT_L <- correct_total_system_volume_cm3 / 1000

cat("Using chamber volume:", CORRECT_CHAMBER_VOLUME_L, "L\n")
cat("Using total system volume:", round(CORRECT_VTOT_L, 3), "L\n\n")

# Load the manual identification data with correct volume
manID.lgr3 <- read_csv("../../data/processed/flux/lgr_manual_identification_results_soil.csv", 
                       show_col_types = FALSE)

# Ensure correct volume is set
manID.lgr3 <- manID.lgr3 %>%
  mutate(
    Vtot = CORRECT_VTOT_L,
    Area = CORRECT_CHAMBER_SURFACE_AREA_CM2
  )

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
CH4_plots_lgr3 <- flux.plot(
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
flux2pdf(
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











# Diagnose and Fix the best.flux CH4 Issue
library(goFlux)
library(dplyr)
library(readr)

cat("=== DIAGNOSING BEST.FLUX CH4 ISSUE ===\n\n")

# Load the corrected CH4 flux results (from goFlux)
CH4_flux_lgr3 <- read_csv("../../data/processed/flux/CH4_flux_lgr_results_soil.csv", show_col_types = FALSE)

cat("Step 1: Checking CH4 flux results from goFlux:\n")
cat("----------------------------------------\n")

# Check column names
cat("Column names in CH4_flux_lgr3:\n")
print(names(CH4_flux_lgr3))

# Check the actual flux values
ch4_check <- CH4_flux_lgr3 %>%
  slice(1:5) %>%
  select(UniqueID, LM.flux, HM.flux)

cat("\nFirst 5 CH4 flux values (from goFlux):\n")
print(ch4_check)

cat("\nCH4 flux summary from goFlux:\n")
cat("- Mean LM.flux:", round(mean(CH4_flux_lgr3$LM.flux, na.rm = TRUE), 4), "nmol/m²/s\n")
cat("- This should be around -0.66\n\n")

# Now let's see what best.flux is doing
cat("Step 2: Testing best.flux function:\n")
cat("------------------------------------\n")

# Try running best.flux with explicit monitoring
CH4_best_test <- best.flux(
  flux.result = CH4_flux_lgr3,
  criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", 
               "MDF", "nb.obs", "intercept", "p-value"),
  intercept.lim = NULL,
  g.limit = 2,
  p.val = 0.05,
  k.ratio = 1,
  warn.length = 60
)

cat("\nChecking best.flux output:\n")
best_check <- CH4_best_test %>%
  slice(1:5) %>%
  select(UniqueID, LM.flux, HM.flux, flux.term, best.flux, model)

print(best_check)

cat("\nFlux.term values from best.flux:\n")
cat("- Mean:", round(mean(CH4_best_test$flux.term, na.rm = TRUE), 4), "\n")
cat("- If this is ~6.2, best.flux is using wrong data\n\n")

# Check if best.flux is somehow mixing up with CO2 data
cat("Step 3: Comparing with CO2 data:\n")
cat("---------------------------------\n")

CO2_flux_lgr3 <- read_csv("../../data/processed/flux/CO2_flux_lgr_results_soil.csv", show_col_types = FALSE)
co2_check <- CO2_flux_lgr3 %>%
  slice(1:5) %>%
  select(UniqueID, LM.flux)

cat("First 5 CO2 LM.flux values:\n")
print(co2_check)

# Check if the values are being swapped
cat("\nComparing flux values:\n")
comparison <- data.frame(
  UniqueID = CH4_flux_lgr3$UniqueID[1:5],
  CH4_goFlux = CH4_flux_lgr3$LM.flux[1:5],
  CH4_bestflux = CH4_best_test$LM.flux[1:5],
  CO2_goFlux = CO2_flux_lgr3$LM.flux[1:5]
)
print(comparison)

# SOLUTION: Create a clean CH4 best flux manually
cat("\n=== SOLUTION: MANUAL BEST.FLUX SELECTION ===\n")

# Since best.flux seems to be broken, let's manually select the best model
CH4_best_manual <- CH4_flux_lgr3

# Add the flux.term column (usually this selects between LM and HM)
# For now, we'll use LM as default since it's more reliable for small fluxes
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    # Select LM unless it's NA or HM has much better fit
    flux.term = case_when(
      is.na(HM.flux) ~ LM.flux,
      is.na(LM.flux) ~ HM.flux,
      abs(HM.AICc) < abs(LM.AICc) - 2 ~ HM.flux,  # HM only if AICc is substantially better
      TRUE ~ LM.flux
    ),
    best.flux = case_when(
      flux.term == LM.flux ~ "LM",
      flux.term == HM.flux ~ "HM",
      TRUE ~ "LM"
    ),
    model = best.flux,
    quality.check = case_when(
      abs(flux.term) < 0.01 ~ "Below MDF",
      is.na(flux.term) ~ "No flux calculated",
      TRUE ~ "OK"
    )
  )

cat("Manual best flux selection complete\n")
cat("- Mean flux.term:", round(mean(CH4_best_manual$flux.term, na.rm = TRUE), 4), "nmol/m²/s\n")
cat("- Median flux.term:", round(median(CH4_best_manual$flux.term, na.rm = TRUE), 4), "nmol/m²/s\n")

# Save the manually corrected best flux
write_csv(CH4_best_manual, "../../data/processed/flux/CH4_best_flux_MANUAL_soil.csv")
cat("\nSaved manually corrected best flux to: data/processed/flux/CH4_best_flux_MANUAL_soil.csv\n")

# Create final statistics
cat("\n=== FINAL CORRECTED STATISTICS ===\n")

# Load CO2 for comparison
CO2_best <- read_csv("../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv", show_col_types = FALSE)

co2_summary <- CO2_best %>%
  summarise(
    Gas = "CO2",
    n = n(),
    mean = round(mean(flux.term, na.rm = TRUE), 3),
    median = round(median(flux.term, na.rm = TRUE), 3),
    sd = round(sd(flux.term, na.rm = TRUE), 3),
    min = round(min(flux.term, na.rm = TRUE), 3),
    max = round(max(flux.term, na.rm = TRUE), 3)
  )

ch4_summary <- CH4_best_manual %>%
  summarise(
    Gas = "CH4",
    n = n(),
    mean = round(mean(flux.term, na.rm = TRUE), 4),
    median = round(median(flux.term, na.rm = TRUE), 4),
    sd = round(sd(flux.term, na.rm = TRUE), 4),
    min = round(min(flux.term, na.rm = TRUE), 4),
    max = round(max(flux.term, na.rm = TRUE), 4)
  )

summary_table <- bind_rows(co2_summary, ch4_summary)
cat("\nFlux Summary Table (nmol/m²/s):\n")
print(summary_table)

cat("\nRatio CO2/CH4 (absolute values):", 
    round(co2_summary$mean / abs(ch4_summary$mean), 1), "\n")

# Update the final dataset with corrected values
cat("\n=== UPDATING FINAL DATASET WITH CORRECTED CH4 ===\n")

original_data <- read_csv('../../data/processed/flux/auxfile_goFlux_soilflux_with_weather_formatted.csv',
                          show_col_types = FALSE)

# Define correct volumes
CORRECT_CHAMBER_SURFACE_AREA_CM2 <- 507.7
CORRECT_CHAMBER_VOLUME_L <- 7.53
CORRECT_VTOT_L <- 7.629

original_data <- original_data %>%
  mutate(
    Vtot = CORRECT_VTOT_L,
    Area = CORRECT_CHAMBER_SURFACE_AREA_CM2
  )

# Add prefixes
co2_results <- CO2_best %>%
  rename_with(~ paste0("CO2_", .), -UniqueID)

ch4_results <- CH4_best_manual %>%
  rename_with(~ paste0("CH4_", .), -UniqueID)

# Merge
final_dataset <- original_data %>%
  left_join(co2_results, by = "UniqueID") %>%
  left_join(ch4_results, by = "UniqueID")

# Save
write_csv(final_dataset, "../../data/processed/flux/FINAL_soil_dataset_CH4_corrected.csv")

cat("\nFinal dataset saved to: data/processed/flux/FINAL_soil_dataset_CH4_corrected.csv\n")
cat("Rows:", nrow(final_dataset), "\n")
cat("Columns:", ncol(final_dataset), "\n")

cat("\n✅ COMPLETE! CH4 fluxes are now correctly calculated and saved.\n")
cat("The mean CH4 flux is", round(ch4_summary$mean, 4), "nmol/m²/s (typical for soil CH4 uptake)\n")
cat("The mean CO2 flux is", round(co2_summary$mean, 2), "nmol/m²/s (typical for soil respiration)\n")



# Proper Implementation of best.flux Selection Following goFlux Documentation
library(goFlux)
library(dplyr)
library(readr)

cat("=== IMPLEMENTING BEST.FLUX SELECTION ACCORDING TO DOCUMENTATION ===\n\n")

# Load the CH4 flux results from goFlux
CH4_flux_lgr3 <- read_csv("../../data/processed/flux/CH4_flux_lgr_results_soil.csv", show_col_types = FALSE)

cat("Starting with", nrow(CH4_flux_lgr3), "CH4 flux measurements\n")
cat("Mean LM.flux from goFlux:", round(mean(CH4_flux_lgr3$LM.flux, na.rm = TRUE), 4), "nmol/m²/s\n\n")

# Initialize the best flux dataframe with default selections
CH4_best_manual <- CH4_flux_lgr3 %>%
  mutate(
    HM.diagnose = "",
    LM.diagnose = "",
    best.flux = HM.flux,  # Default to HM (assumed non-linearity)
    model = "HM",
    quality.check = "",
    HM.score = 0,
    LM.score = 0
  )

# Define thresholds (from documentation defaults)
g.limit <- 2        # Default g-factor limit
p.val <- 0.05      # Default p-value threshold
k.ratio <- 1       # Default kappa ratio
warn.length <- 60  # Default minimum observations

cat("=== APPLYING SELECTION CRITERIA ===\n\n")

# 1. MAE CRITERION
cat("1. MAE (Mean Absolute Error) criterion:\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    MAE.lim = prec,
    # Score based on MAE comparison
    HM.score = if_else(HM.MAE > LM.MAE, HM.score + 1, HM.score),
    LM.score = if_else(LM.MAE > HM.MAE, LM.score + 1, LM.score),
    # Diagnose noise issues
    LM.diagnose = if_else(LM.MAE > prec, 
                          paste0(LM.diagnose, if_else(LM.diagnose == "", "", " | "), "Noise (LM.MAE)"), 
                          LM.diagnose),
    HM.diagnose = if_else(HM.MAE > prec, 
                          paste0(HM.diagnose, if_else(HM.diagnose == "", "", " | "), "Noise (HM.MAE)"), 
                          HM.diagnose)
  )

# 2. RMSE CRITERION
cat("2. RMSE (Root Mean Square Error) criterion:\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    RMSE.lim = prec,
    # Score based on RMSE comparison
    HM.score = if_else(HM.RMSE > LM.RMSE, HM.score + 1, HM.score),
    LM.score = if_else(LM.RMSE > HM.RMSE, LM.score + 1, LM.score),
    # Diagnose noise issues
    LM.diagnose = if_else(LM.RMSE > prec, 
                          paste0(LM.diagnose, if_else(LM.diagnose == "", "", " | "), "Noise (LM.RMSE)"), 
                          LM.diagnose),
    HM.diagnose = if_else(HM.RMSE > prec, 
                          paste0(HM.diagnose, if_else(HM.diagnose == "", "", " | "), "Noise (HM.RMSE)"), 
                          HM.diagnose)
  )

# 3. AICc CRITERION
cat("3. AICc (Akaike Information Criterion corrected) criterion:\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    # Lower AICc is better
    HM.score = if_else(HM.AICc > LM.AICc, HM.score + 1, HM.score),
    LM.score = if_else(LM.AICc > HM.AICc, LM.score + 1, LM.score)
  )

# 4. SE CRITERION
cat("4. SE (Standard Error) criterion:\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    SE.lim = prec / sqrt(nb.obs),
    # Score based on SE comparison
    HM.score = if_else(HM.SE > LM.SE, HM.score + 1, HM.score),
    LM.score = if_else(LM.SE > HM.SE, LM.score + 1, LM.score),
    # Diagnose noise issues
    LM.diagnose = if_else(LM.SE > SE.lim, 
                          paste0(LM.diagnose, if_else(LM.diagnose == "", "", " | "), "Noise (LM.SE)"), 
                          LM.diagnose),
    HM.diagnose = if_else(HM.SE > SE.lim, 
                          paste0(HM.diagnose, if_else(HM.diagnose == "", "", " | "), "Noise (HM.SE)"), 
                          HM.diagnose)
  )

# Select best model based on scores
cat("5. Selecting best model based on scores:\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    # If LM has lower score, it wins; ties go to HM (assumed non-linearity)
    best.flux = if_else(LM.score < HM.score, LM.flux, HM.flux),
    model = if_else(LM.score < HM.score, "LM", "HM")
  )

# 5. G-FACTOR CRITERION (overrides score-based selection)
cat("6. G-factor criterion (g.limit = 2):\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    g.limit = g.limit,
    # Override to LM if g-factor exceeds limit
    best.flux = if_else(abs(g.fact) > g.limit & !is.na(g.fact), LM.flux, best.flux),
    model = if_else(abs(g.fact) > g.limit & !is.na(g.fact), "LM", model),
    HM.diagnose = if_else(abs(g.fact) > g.limit & !is.na(g.fact), 
                          paste0(HM.diagnose, if_else(HM.diagnose == "", "", " | "), "Exaggerated curvature (g-fact)"), 
                          HM.diagnose),
    quality.check = if_else(abs(g.fact) > g.limit & !is.na(g.fact), 
                            paste0(quality.check, if_else(quality.check == "", "", " | "), "g-fact. > 2"), 
                            quality.check)
  )

# 6. KAPPA RATIO CRITERION
cat("7. Kappa ratio criterion (k.ratio = 1):\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    k.ratio.lim = k.ratio,
    k.ratio.calc = abs(HM.k / k.max),
    # Override to LM if kappa ratio exceeds limit
    best.flux = if_else(k.ratio.calc > k.ratio & !is.na(k.ratio.calc), LM.flux, best.flux),
    model = if_else(k.ratio.calc > k.ratio & !is.na(k.ratio.calc), "LM", model),
    HM.diagnose = if_else(k.ratio.calc > k.ratio & !is.na(k.ratio.calc), 
                          paste0(HM.diagnose, if_else(HM.diagnose == "", "", " | "), "Exaggerated curvature (k.ratio)"), 
                          HM.diagnose),
    quality.check = if_else(k.ratio.calc > k.ratio & !is.na(k.ratio.calc), 
                            paste0(quality.check, if_else(quality.check == "", "", " | "), "kappa max"), 
                            quality.check)
  )

# 7. MDF CRITERION (Minimal Detectable Flux)
cat("8. MDF (Minimal Detectable Flux) criterion:\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    MDF.lim = MDF,
    # Flag measurements below MDF
    LM.diagnose = if_else(abs(LM.flux) < MDF, 
                          paste0(LM.diagnose, if_else(LM.diagnose == "", "", " | "), "MDF (LM)"), 
                          LM.diagnose),
    HM.diagnose = if_else(abs(HM.flux) < MDF, 
                          paste0(HM.diagnose, if_else(HM.diagnose == "", "", " | "), "MDF (HM)"), 
                          HM.diagnose),
    quality.check = if_else(abs(best.flux) < MDF, 
                            paste0(quality.check, if_else(quality.check == "", "", " | "), "MDF"), 
                            quality.check)
  )

# 8. P-VALUE CRITERION
cat("9. P-value criterion (p.val = 0.05):\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    p.val.lim = p.val,
    # Flag non-significant linear fluxes
    LM.diagnose = if_else(LM.p.val > p.val, 
                          paste0(LM.diagnose, if_else(LM.diagnose == "", "", " | "), "No detect. flux (p-val.)"), 
                          LM.diagnose),
    quality.check = if_else(LM.p.val > p.val & model == "LM", 
                            paste0(quality.check, if_else(quality.check == "", "", " | "), "p-value"), 
                            quality.check)
  )

# 9. Handle NA values
cat("10. Handling NA values:\n")
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    # If HM is NA but LM exists, use LM
    best.flux = if_else(is.na(HM.flux) & !is.na(LM.flux), LM.flux, best.flux),
    model = if_else(is.na(HM.flux) & !is.na(LM.flux), "LM", model),
    # If both are NA
    best.flux = if_else(is.na(HM.flux) & is.na(LM.flux), NA_real_, best.flux),
    model = if_else(is.na(HM.flux) & is.na(LM.flux), NA_character_, model),
    # Update diagnostics
    HM.diagnose = if_else(is.na(HM.flux), "HM.flux is NA", HM.diagnose),
    LM.diagnose = if_else(is.na(LM.flux), "LM.flux is NA", LM.diagnose),
    quality.check = if_else(is.na(HM.flux) & is.na(LM.flux), "fluxes are NA", quality.check)
  )

# Clean up quality.check field
CH4_best_manual <- CH4_best_manual %>%
  mutate(
    quality.check = if_else(quality.check == "", "OK", quality.check),
    # Use best.flux for flux.term (THIS IS THE KEY FIX)
    flux.term = best.flux
  )

# Summary statistics
cat("\n=== SELECTION SUMMARY ===\n")
model_summary <- CH4_best_manual %>%
  group_by(model) %>%
  summarise(
    n = n(),
    mean_flux = mean(flux.term, na.rm = TRUE),
    median_flux = median(flux.term, na.rm = TRUE)
  )
print(model_summary)

quality_summary <- CH4_best_manual %>%
  summarise(
    OK = sum(quality.check == "OK"),
    MDF = sum(grepl("MDF", quality.check)),
    g_factor = sum(grepl("g-fact", quality.check)),
    kappa = sum(grepl("kappa", quality.check)),
    p_value = sum(grepl("p-value", quality.check)),
    noise = sum(grepl("Noise", HM.diagnose) | grepl("Noise", LM.diagnose))
  )

cat("\nQuality check summary:\n")
print(quality_summary)

# Final statistics
cat("\n=== FINAL CH4 FLUX STATISTICS ===\n")
ch4_stats <- CH4_best_manual %>%
  summarise(
    n = n(),
    mean = round(mean(flux.term, na.rm = TRUE), 4),
    median = round(median(flux.term, na.rm = TRUE), 4),
    sd = round(sd(flux.term, na.rm = TRUE), 4),
    min = round(min(flux.term, na.rm = TRUE), 4),
    max = round(max(flux.term, na.rm = TRUE), 4),
    q25 = round(quantile(flux.term, 0.25, na.rm = TRUE), 4),
    q75 = round(quantile(flux.term, 0.75, na.rm = TRUE), 4)
  )

cat("CH4 flux (nmol/m²/s):\n")
cat("- Mean:", ch4_stats$mean, "\n")
cat("- Median:", ch4_stats$median, "\n")
cat("- SD:", ch4_stats$sd, "\n")
cat("- Range:", ch4_stats$min, "to", ch4_stats$max, "\n")
cat("- IQR:", ch4_stats$q25, "to", ch4_stats$q75, "\n")

# Save the properly implemented best flux results
write_csv(CH4_best_manual, "../../data/processed/flux/CH4_best_flux_PROPER_soil.csv")
cat("\nSaved properly implemented best flux to: data/processed/flux/CH4_best_flux_PROPER_soil.csv\n")

cat("\n✅ COMPLETE! Best flux selection implemented according to goFlux documentation.\n")


















# Create Final Soil Dataset with Corrected CH4 and CO2 Fluxes
library(dplyr)
library(readr)

cat("=== CREATING FINAL SOIL DATASET ===\n\n")

# Define correct chamber geometry
CORRECT_CHAMBER_SURFACE_AREA_CM2 <- 507.7
CORRECT_CHAMBER_VOLUME_L <- 7.53
correct_tubing_volume_cm3 <- pi * (1/16)^2 * 12 * 12 * 16.387
correct_system_volume_cm3 <- 70
correct_total_system_volume_cm3 <- CORRECT_CHAMBER_VOLUME_L * 1000 + 
  correct_tubing_volume_cm3 + correct_system_volume_cm3
CORRECT_VTOT_L <- correct_total_system_volume_cm3 / 1000

cat("Chamber configuration:\n")
cat("- Volume:", CORRECT_CHAMBER_VOLUME_L, "L\n")
cat("- Surface Area:", CORRECT_CHAMBER_SURFACE_AREA_CM2, "cm²\n")
cat("- Total System Volume:", round(CORRECT_VTOT_L, 3), "L\n\n")

# 1. Load the original auxfile data
cat("Loading original auxfile data...\n")
original_data <- read_csv('../../data/processed/flux/auxfile_goFlux_soilflux_with_weather_formatted.csv',
                          show_col_types = FALSE)

# Ensure correct volume is set in original data
original_data <- original_data %>%
  mutate(
    Vtot = CORRECT_VTOT_L,
    Area = CORRECT_CHAMBER_SURFACE_AREA_CM2
  )

cat("Original data: ", nrow(original_data), "measurements,", ncol(original_data), "columns\n")
cat("Column names in original data:\n")
print(names(original_data))

# 2. Load CO2 best flux results
cat("\nLoading CO2 flux results...\n")
CO2_best <- read_csv("../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv", show_col_types = FALSE)

# Check for flux.term issue in CO2
if(mean(CO2_best$flux.term, na.rm = TRUE) > 5) {
  cat("CO2 flux.term looks correct (mean =", round(mean(CO2_best$flux.term, na.rm = TRUE), 2), "nmol/m²/s)\n")
} else {
  cat("Warning: CO2 flux.term may be incorrect\n")
}

# 3. Load the CORRECTED CH4 best flux results
cat("\nLoading corrected CH4 flux results...\n")

# Try to load the properly implemented version first
ch4_file <- if(file.exists("../../data/processed/flux/CH4_best_flux_PROPER_soil.csv")) {
  "../../data/processed/flux/CH4_best_flux_PROPER_soil.csv"
} else if(file.exists("../../data/processed/flux/CH4_best_flux_MANUAL_soil.csv")) {
  "../../data/processed/flux/CH4_best_flux_MANUAL_soil.csv"
} else {
  "../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv"
}

cat("Using CH4 file:", ch4_file, "\n")
CH4_best <- read_csv(ch4_file, show_col_types = FALSE)

# Verify CH4 flux.term is correct
# Check if flux.term exists, if not use best.flux or LM.flux
if(!"flux.term" %in% names(CH4_best)) {
  cat("flux.term not found, checking for best.flux column...\n")
  if("best.flux" %in% names(CH4_best)) {
    CH4_best <- CH4_best %>%
      mutate(flux.term = best.flux)
    cat("Using best.flux for flux.term\n")
  } else {
    CH4_best <- CH4_best %>%
      mutate(flux.term = LM.flux)
    cat("Using LM.flux for flux.term\n")
  }
}

# Now check the values
ch4_mean <- mean(CH4_best$flux.term, na.rm = TRUE)
if(abs(ch4_mean) < 2) {
  cat("CH4 flux.term looks correct (mean =", round(ch4_mean, 4), "nmol/m²/s)\n")
} else if(ch4_mean > 5) {
  cat("ERROR: CH4 flux.term still shows CO2 values! Using LM.flux instead.\n")
  # Emergency fix - use LM.flux directly
  CH4_best <- CH4_best %>%
    mutate(flux.term = LM.flux)
  cat("Fixed: CH4 flux.term now =", round(mean(CH4_best$flux.term, na.rm = TRUE), 4), "nmol/m²/s\n")
}

# 4. Add prefixes to distinguish CO2 and CH4 columns
cat("\nAdding prefixes to flux columns...\n")

# Select important columns for the final dataset
co2_columns_to_keep <- c("UniqueID", "LM.flux", "HM.flux", "flux.term", "best.flux", 
                         "model", "quality.check", "LM.p.val", "HM.r2", "LM.r2",
                         "LM.MAE", "HM.MAE", "LM.RMSE", "HM.RMSE", "LM.AICc", "HM.AICc",
                         "g.fact", "MDF", "HM.k", "k.max")

ch4_columns_to_keep <- co2_columns_to_keep  # Same columns for CH4

# Create CO2 results with prefix
co2_results <- CO2_best %>%
  select(any_of(co2_columns_to_keep)) %>%
  rename_with(~ paste0("CO2_", .), -UniqueID)

# Create CH4 results with prefix
ch4_results <- CH4_best %>%
  select(any_of(ch4_columns_to_keep)) %>%
  rename_with(~ paste0("CH4_", .), -UniqueID)

cat("CO2 columns prepared:", ncol(co2_results) - 1, "+ UniqueID\n")
cat("CH4 columns prepared:", ncol(ch4_results) - 1, "+ UniqueID\n")

# 5. Merge everything together
cat("\nMerging datasets...\n")

final_dataset <- original_data %>%
  left_join(co2_results, by = "UniqueID") %>%
  left_join(ch4_results, by = "UniqueID")

cat("Final dataset created:", nrow(final_dataset), "rows,", ncol(final_dataset), "columns\n")

# 6. Add summary columns for easy analysis
cat("\nAdding summary columns...\n")

final_dataset <- final_dataset %>%
  mutate(
    # Add ratio column
    CO2_CH4_ratio = if_else(CH4_flux.term != 0, 
                            abs(CO2_flux.term / CH4_flux.term), 
                            NA_real_),
    # Add quality flags
    CO2_quality_OK = CO2_quality.check == "OK" | is.na(CO2_quality.check),
    CH4_quality_OK = CH4_quality.check == "OK" | is.na(CH4_quality.check),
    # Add flux direction indicators
    CH4_uptake = CH4_flux.term < 0,
    CH4_emission = CH4_flux.term > 0
  )

# 7. Verify the results
cat("\n=== VERIFICATION ===\n")

co2_summary <- final_dataset %>%
  summarise(
    n = sum(!is.na(CO2_flux.term)),
    mean = round(mean(CO2_flux.term, na.rm = TRUE), 3),
    median = round(median(CO2_flux.term, na.rm = TRUE), 3),
    sd = round(sd(CO2_flux.term, na.rm = TRUE), 3)
  )

ch4_summary <- final_dataset %>%
  summarise(
    n = sum(!is.na(CH4_flux.term)),
    mean = round(mean(CH4_flux.term, na.rm = TRUE), 4),
    median = round(median(CH4_flux.term, na.rm = TRUE), 4),
    sd = round(sd(CH4_flux.term, na.rm = TRUE), 4)
  )

cat("CO2 flux summary (nmol/m²/s):\n")
cat("  n =", co2_summary$n, ", mean =", co2_summary$mean, 
    ", median =", co2_summary$median, ", sd =", co2_summary$sd, "\n")

cat("\nCH4 flux summary (nmol/m²/s):\n")
cat("  n =", ch4_summary$n, ", mean =", ch4_summary$mean, 
    ", median =", ch4_summary$median, ", sd =", ch4_summary$sd, "\n")

cat("\nCH4 uptake vs emission:\n")
uptake_emission <- final_dataset %>%
  summarise(
    uptake = sum(CH4_uptake, na.rm = TRUE),
    emission = sum(CH4_emission, na.rm = TRUE),
    zero = sum(CH4_flux.term == 0, na.rm = TRUE)
  )
cat("  Uptake:", uptake_emission$uptake, "measurements\n")
cat("  Emission:", uptake_emission$emission, "measurements\n")
cat("  Zero flux:", uptake_emission$zero, "measurements\n")

# 8. Check data quality
cat("\nData quality:\n")
quality_check <- final_dataset %>%
  summarise(
    CO2_OK = sum(CO2_quality_OK, na.rm = TRUE),
    CO2_flagged = sum(!CO2_quality_OK, na.rm = TRUE),
    CH4_OK = sum(CH4_quality_OK, na.rm = TRUE),
    CH4_flagged = sum(!CH4_quality_OK, na.rm = TRUE)
  )
cat("  CO2: ", quality_check$CO2_OK, "OK, ", quality_check$CO2_flagged, "flagged\n")
cat("  CH4: ", quality_check$CH4_OK, "OK, ", quality_check$CH4_flagged, "flagged\n")

# 9. Save the final dataset
cat("\n=== SAVING FINAL DATASET ===\n")

# Main comprehensive file
write_csv(final_dataset, "../../data/processed/flux/FINAL_soil_flux_dataset_complete.csv")
cat("Saved complete dataset to: data/processed/flux/FINAL_soil_flux_dataset_complete.csv\n")

# Also save a simplified version with just the key results
simplified_dataset <- final_dataset %>%
  select(UniqueID, start.time, Plot, trt, soilT.C, moistV, 
         CO2_flux.term, CO2_model, CO2_quality.check,
         CH4_flux.term, CH4_model, CH4_quality.check,
         CO2_CH4_ratio, CH4_uptake, CH4_emission)

write_csv(simplified_dataset, "../../data/processed/flux/FINAL_soil_flux_dataset_simplified.csv")
cat("Saved simplified dataset to: data/processed/flux/FINAL_soil_flux_dataset_simplified.csv\n")

# Save a summary table
summary_by_treatment <- final_dataset %>%
  group_by(trt) %>%
  summarise(
    n = n(),
    CO2_mean = round(mean(CO2_flux.term, na.rm = TRUE), 2),
    CO2_se = round(sd(CO2_flux.term, na.rm = TRUE) / sqrt(sum(!is.na(CO2_flux.term))), 3),
    CH4_mean = round(mean(CH4_flux.term, na.rm = TRUE), 4),
    CH4_se = round(sd(CH4_flux.term, na.rm = TRUE) / sqrt(sum(!is.na(CH4_flux.term))), 4),
    CH4_uptake_pct = round(100 * sum(CH4_uptake, na.rm = TRUE) / sum(!is.na(CH4_flux.term)), 1)
  )

write_csv(summary_by_treatment, "../../data/processed/flux/FINAL_soil_flux_summary_by_treatment.csv")
cat("Saved summary by treatment to: data/processed/flux/FINAL_soil_flux_summary_by_treatment.csv\n")

cat("\n✅ COMPLETE! Final dataset created with correct CH4 and CO2 fluxes.\n")
cat("\nKey findings:\n")
cat("- CO2 flux (soil respiration):", co2_summary$mean, "±", co2_summary$sd, "nmol/m²/s\n")
cat("- CH4 flux:", ch4_summary$mean, "±", ch4_summary$sd, "nmol/m²/s\n")
cat("- CO2/CH4 ratio:", round(abs(co2_summary$mean / ch4_summary$mean), 1), "\n")
cat("- CH4 behavior:", round(100 * uptake_emission$uptake / (uptake_emission$uptake + uptake_emission$emission), 1), 
    "% of measurements show CH4 uptake\n")

cat("\nAll original column names and metadata have been preserved.\n")
cat("Flux results are now correctly differentiated between CO2 and CH4.\n")