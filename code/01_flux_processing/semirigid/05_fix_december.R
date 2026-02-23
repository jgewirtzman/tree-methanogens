# ==============================================================================
# Fix December 2020 Flux Data
# ==============================================================================
# Purpose: Handles December 2020 data processing separately due to date
#   handling issues in the main pipeline.
#
# Pipeline stage: 01 Flux Processing (semirigid)
# Run after: 04_goflux_trees.R / 04_goflux_soils.R
#
# Inputs:
#   - Raw LGR3 December data files
#
# Outputs:
#   - Updates existing flux result files with December data
# ==============================================================================

library(goFlux)
library(dplyr)
library(readr)
library(lubridate)

# =============================================================================
# STEP 1: LOAD AND PREPARE LGR3 DATA (DECEMBER ONLY)
# =============================================================================

cat("=== STEP 1: LOADING LGR3 DATA (DECEMBER ONLY) ===\n")

# Set data path
data_path <- '../../../data/raw/lgr/semirigid_2020-2021'

# Check for zip files that need extraction
zip_files <- list.files(data_path, recursive = TRUE, pattern = "\\.zip$", full.names = TRUE)
cat("Found", length(zip_files), "zip files\n")

# Extract zip files to temporary directory
temp_extract_dir <- tempfile("lgr3_extract_")
dir.create(temp_extract_dir, recursive = TRUE)

cat("Extracting zip files...\n")
for(zip_file in zip_files) {
  tryCatch({
    unzip(zip_file, exdir = temp_extract_dir, overwrite = TRUE)
    cat("Extracted:", basename(zip_file), "\n")
  }, error = function(e) {
    cat("Error with", basename(zip_file), ":", e$message, "\n")
  })
}

# Combine existing and extracted files
existing_data_files <- list.files(data_path, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE)
existing_data_files <- existing_data_files[file.size(existing_data_files) > 0]

new_txt_files <- list.files(temp_extract_dir, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE)

# Create complete data directory
complete_lgr3_dir <- tempfile("lgr3_complete_")
dir.create(complete_lgr3_dir, recursive = TRUE)

# Copy all files
file.copy(existing_data_files, complete_lgr3_dir)
file.copy(new_txt_files, complete_lgr3_dir)

cat("Total files prepared:", length(list.files(complete_lgr3_dir, pattern = "\\.txt$")), "\n")

# Import complete LGR3 dataset
cat("Importing LGR3 data...\n")
lgr3_data_complete <- import2RData(
  path = complete_lgr3_dir,
  instrument = "UGGA",
  date.format = "mdy",
  timezone = "UTC",
  keep_all = FALSE,
  prec = c(0.35, 0.9, 200),
  merge = TRUE
)

# Clean up temporary directories
unlink(temp_extract_dir, recursive = TRUE)
unlink(complete_lgr3_dir, recursive = TRUE)

cat("LGR3 data loaded:", dim(lgr3_data_complete), "\n")

# =============================================================================
# STEP 2: FILTER FOR DECEMBER DATES ONLY
# =============================================================================

cat("\n=== STEP 2: FILTERING FOR DECEMBER DATES ===\n")

# Load auxfile to identify December measurements
lgr3_auxfile_full <- read_csv("auxfile_goFlux_soilflux_with_weather.csv") %>%
  mutate(start.time = as.POSIXct(start.time, tz = "UTC"),
         month = month(start.time))

# Filter auxfile for December only
lgr3_auxfile_december <- lgr3_auxfile_full %>%
  filter(month == 12)

cat("Total measurements in auxfile:", nrow(lgr3_auxfile_full), "\n")
cat("December measurements:", nrow(lgr3_auxfile_december), "\n")

# Create observation windows for December only
cat("Creating observation windows for December measurements...\n")
ow.lgr3.december <- obs.win(
  inputfile = lgr3_data_complete,
  auxfile = lgr3_auxfile_december,
  gastype = "CO2dry_ppm",
  obs.length = 300,
  shoulder = 300
)

# Check window quality
window_sizes <- sapply(ow.lgr3.december, nrow)
good_windows <- which(window_sizes >= 30)

# Assign names
unique_ids <- sapply(ow.lgr3.december, function(x) unique(x$UniqueID)[1])
names(ow.lgr3.december) <- unique_ids

cat("December observation windows:", length(ow.lgr3.december), "\n")
cat("Good December windows (>=30 obs):", length(good_windows), "\n")

# =============================================================================
# STEP 3: MANUAL IDENTIFICATION FOR DECEMBER ONLY
# =============================================================================

cat("\n=== STEP 3: MANUAL IDENTIFICATION (DECEMBER ONLY) ===\n")

# Process in batches of 20
batch_size <- 20
total_measurements <- length(good_windows)
num_batches <- ceiling(total_measurements / batch_size)

cat("Processing", total_measurements, "December measurements in", num_batches, "batches\n")

manID_batches <- list()

for(batch_num in 1:num_batches) {
  start_idx <- (batch_num - 1) * batch_size + 1
  end_idx <- min(batch_num * batch_size, total_measurements)
  batch_windows <- good_windows[start_idx:end_idx]
  
  cat("\n=== Processing December Batch", batch_num, "of", num_batches, "===\n")
  cat("Measurements", start_idx, "to", end_idx, "\n")
  
  # Process this batch
  manID_batch <- click.peak2(
    ow.lgr3.december,
    seq = batch_windows,
    gastype = "CO2dry_ppm",
    sleep = 3,
    plot.lim = c(200, 5000),
    save.plots = paste0("../../../data/processed/flux/lgr3_december_batch_", batch_num, "_plots")
  )
  
  # Store the batch result
  manID_batches[[batch_num]] <- manID_batch
  
  # Force close graphics devices
  while(dev.cur() > 1) dev.off()
  
  cat("December Batch", batch_num, "complete. Processed", nrow(manID_batch), "measurements\n")
  cat("Press Enter to continue to next batch...\n")
  readline()
}

# Combine all batches
manID.lgr3.december <- do.call(rbind, manID_batches)

cat("\nDecember manual identification complete! Total:", nrow(manID.lgr3.december), "measurements\n")

# Save December manual identification results
write_csv(manID.lgr3.december, "../../../data/processed/flux/lgr_manual_identification_results_december_soil.csv")

# =============================================================================
# STEP 4: FLUX CALCULATIONS FOR DECEMBER
# =============================================================================

cat("\n=== STEP 4: DECEMBER FLUX CALCULATIONS ===\n")

# Calculate CO2 fluxes for December
cat("Calculating December CO2 fluxes...\n")
CO2_flux_lgr3_december <- goFlux(
  dataframe = manID.lgr3.december,
  gastype = "CO2dry_ppm",
  H2O_col = "H2O_ppm",
  warn.length = 60
)

# Calculate CH4 fluxes if available
if("CH4dry_ppb" %in% names(manID.lgr3.december)) {
  cat("Calculating December CH4 fluxes...\n")
  CH4_flux_lgr3_december <- goFlux(
    dataframe = manID.lgr3.december,
    gastype = "CH4dry_ppb",
    H2O_col = "H2O_ppm",
    warn.length = 60
  )
} else {
  cat("CH4dry_ppb column not found - skipping CH4 flux calculation\n")
}

cat("December CO2 flux calculation complete:", nrow(CO2_flux_lgr3_december), "measurements\n")

# =============================================================================
# STEP 5: BEST FLUX ANALYSIS FOR DECEMBER
# =============================================================================

cat("\n=== STEP 5: DECEMBER BEST FLUX ANALYSIS ===\n")

# Run best.flux on December CO2 results
CO2_best_lgr3_december <- best.flux(
  flux.result = CO2_flux_lgr3_december,
  criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", "MDF", "nb.obs", "intercept", "p-value"),
  intercept.lim = NULL,
  g.limit = 2,
  p.val = 0.05,
  k.ratio = 1,
  warn.length = 60
)

# Run best.flux on December CH4 results if available
if(exists("CH4_flux_lgr3_december")) {
  CH4_best_lgr3_december <- best.flux(
    flux.result = CH4_flux_lgr3_december,
    criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", "MDF", "nb.obs", "intercept", "p-value"),
    intercept.lim = NULL,
    g.limit = 2,
    p.val = 0.05,
    k.ratio = 1,
    warn.length = 60
  )
}

# Quality summary for December
quality_summary_december <- CO2_best_lgr3_december %>%
  summarise(
    clean_measurements = sum(quality.check == "OK", na.rm = TRUE),
    flagged_measurements = sum(quality.check != "OK", na.rm = TRUE),
    below_MDF = sum(grepl("MDF", quality.check), na.rm = TRUE),
    high_g_factor = sum(grepl("g.factor", quality.check), na.rm = TRUE),
    poor_fit = sum(grepl("MAE|RMSE|SE", quality.check), na.rm = TRUE)
  )

cat("December Quality Assessment:\n")
cat("Clean measurements:", quality_summary_december$clean_measurements, "\n")
cat("Flagged measurements:", quality_summary_december$flagged_measurements, "\n")

# =============================================================================
# STEP 6: CREATE DECEMBER PLOTS
# =============================================================================

cat("\n=== STEP 6: CREATING DECEMBER PLOTS ===\n")

# Create December CO2 flux plots
CO2_plots_lgr3_december <- flux.plot(
  flux.results = CO2_best_lgr3_december,
  dataframe = manID.lgr3.december,
  gastype = "CO2dry_ppm",
  shoulder = 30,
  plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
  quality.check = TRUE,
  best.model = TRUE,
  p.val.disp = "round"
)

# Create December CH4 plots if available
if(exists("CH4_best_lgr3_december")) {
  CH4_plots_lgr3_december <- flux.plot(
    flux.results = CH4_best_lgr3_december,
    dataframe = manID.lgr3.december,
    gastype = "CH4dry_ppb",
    shoulder = 30,
    plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
    plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
    quality.check = TRUE,
    best.model = TRUE,
    p.val.disp = "round"
  )
  
  # Combine plots
  all_plots_lgr3_december <- c(CO2_plots_lgr3_december, CH4_plots_lgr3_december)
} else {
  all_plots_lgr3_december <- CO2_plots_lgr3_december
}

# Save December plots to PDF
flux2pdf(
  plot.list = all_plots_lgr3_december,
  outfile = "../../../data/processed/flux/LGR3_flux_plots_december_soil.pdf",
  width = 11.6,
  height = 8.2
)
cat("December plots saved successfully\n")

# =============================================================================
# STEP 7: UPDATE ALL ORIGINAL OUTPUT FILES WITH DECEMBER DATA
# =============================================================================

cat("\n=== STEP 7: UPDATING ALL ORIGINAL OUTPUT FILES WITH DECEMBER DATA ===\n")

# Get UniqueIDs for December 14th measurements (new data)
december_14_unique_ids <- lgr3_auxfile_december$UniqueID

# Get UniqueIDs for December 8th measurements (to remove)
# These will be any UniqueID containing "20201208"
december_8_pattern <- "20201208"

# =============================================================================
# UPDATE MANUAL IDENTIFICATION RESULTS
# =============================================================================

cat("Updating manual identification results...\n")

# Load existing manual ID results
if(file.exists("../../../data/processed/flux/lgr_manual_identification_results_soil.csv")) {
  existing_manID <- read_csv("../../../data/processed/flux/lgr_manual_identification_results_soil.csv") %>%
    mutate(DATE = as.character(DATE))  # Convert DATE to character to match new data
  
  # Remove existing December 8th data
  non_december_8_manID <- existing_manID %>%
    filter(!grepl(december_8_pattern, UniqueID))
  
  # Convert new December data DATE column to character to match
  manID.lgr3.december_fixed <- manID.lgr3.december %>%
    mutate(DATE = as.character(DATE))
  
  # Combine with new December 14th data
  updated_manID <- bind_rows(non_december_8_manID, manID.lgr3.december_fixed)
  
  # Save updated manual ID results
  write_csv(updated_manID, "../../../data/processed/flux/lgr_manual_identification_results_soil.csv")
  
  cat("Manual ID results updated:", nrow(updated_manID), "total rows\n")
  cat("December 8th rows removed:", nrow(existing_manID) - nrow(non_december_8_manID), "\n")
  cat("December 14th rows added:", nrow(manID.lgr3.december_fixed), "\n")
} else {
  # If original doesn't exist, just save December data
  write_csv(manID.lgr3.december, "../../../data/processed/flux/lgr_manual_identification_results_soil.csv")
  cat("Created new manual ID results file\n")
}

# =============================================================================
# UPDATE CO2 FLUX RESULTS
# =============================================================================

cat("Updating CO2 flux results...\n")

# Load existing CO2 flux results
if(file.exists("../../../data/processed/flux/CO2_flux_lgr_results_soil.csv")) {
  existing_CO2_flux <- read_csv("../../../data/processed/flux/CO2_flux_lgr_results_soil.csv")
  
  # Remove existing December 8th data
  non_december_8_CO2_flux <- existing_CO2_flux %>%
    filter(!grepl(december_8_pattern, UniqueID))
  
  # Combine with new December 14th data
  updated_CO2_flux <- bind_rows(non_december_8_CO2_flux, CO2_flux_lgr3_december)
  
  # Save updated CO2 flux results
  write_csv(updated_CO2_flux, "../../../data/processed/flux/CO2_flux_lgr_results_soil.csv")
  
  cat("CO2 flux results updated:", nrow(updated_CO2_flux), "total rows\n")
  cat("December 8th rows removed:", nrow(existing_CO2_flux) - nrow(non_december_8_CO2_flux), "\n")
  cat("December 14th rows added:", nrow(CO2_flux_lgr3_december), "\n")
} else {
  # If original doesn't exist, just save December data
  write_csv(CO2_flux_lgr3_december, "../../../data/processed/flux/CO2_flux_lgr_results_soil.csv")
  cat("Created new CO2 flux results file\n")
}

# =============================================================================
# UPDATE CO2 BEST FLUX RESULTS
# =============================================================================

cat("Updating CO2 best flux results...\n")

# Load existing CO2 best flux results
if(file.exists("../../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv")) {
  existing_CO2_best <- read_csv("../../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv")
  
  # Remove existing December 8th data
  non_december_8_CO2_best <- existing_CO2_best %>%
    filter(!grepl(december_8_pattern, UniqueID))
  
  # Combine with new December 14th data
  updated_CO2_best <- bind_rows(non_december_8_CO2_best, CO2_best_lgr3_december)
  
  # Save updated CO2 best flux results
  write_csv(updated_CO2_best, "../../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv")
  
  cat("CO2 best flux results updated:", nrow(updated_CO2_best), "total rows\n")
  cat("December 8th rows removed:", nrow(existing_CO2_best) - nrow(non_december_8_CO2_best), "\n")
  cat("December 14th rows added:", nrow(CO2_best_lgr3_december), "\n")
} else {
  # If original doesn't exist, just save December data
  write_csv(CO2_best_lgr3_december, "../../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv")
  cat("Created new CO2 best flux results file\n")
}

# =============================================================================
# UPDATE CH4 FLUX RESULTS (IF AVAILABLE)
# =============================================================================

if(exists("CH4_flux_lgr3_december")) {
  cat("Updating CH4 flux results...\n")
  
  # Load existing CH4 flux results
  if(file.exists("../../../data/processed/flux/CH4_flux_lgr_results_soil.csv")) {
    existing_CH4_flux <- read_csv("../../../data/processed/flux/CH4_flux_lgr_results_soil.csv")
    
    # Remove existing December 8th data
    non_december_8_CH4_flux <- existing_CH4_flux %>%
      filter(!grepl(december_8_pattern, UniqueID))
    
    # Combine with new December 14th data
    updated_CH4_flux <- bind_rows(non_december_8_CH4_flux, CH4_flux_lgr3_december)
    
    # Save updated CH4 flux results
    write_csv(updated_CH4_flux, "../../../data/processed/flux/CH4_flux_lgr_results_soil.csv")
    
    cat("CH4 flux results updated:", nrow(updated_CH4_flux), "total rows\n")
    cat("December 8th rows removed:", nrow(existing_CH4_flux) - nrow(non_december_8_CH4_flux), "\n")
    cat("December 14th rows added:", nrow(CH4_flux_lgr3_december), "\n")
  } else {
    # If original doesn't exist, just save December data
    write_csv(CH4_flux_lgr3_december, "../../../data/processed/flux/CH4_flux_lgr_results_soil.csv")
    cat("Created new CH4 flux results file\n")
  }
}

# =============================================================================
# UPDATE CH4 BEST FLUX RESULTS (IF AVAILABLE)
# =============================================================================

if(exists("CH4_best_lgr3_december")) {
  cat("Updating CH4 best flux results...\n")
  
  # Load existing CH4 best flux results
  if(file.exists("../../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")) {
    existing_CH4_best <- read_csv("../../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")
    
    # Remove existing December data
    non_december_CH4_best <- existing_CH4_best %>%
      filter(!UniqueID %in% december_unique_ids)
    
    # Combine with new December data
    updated_CH4_best <- bind_rows(non_december_CH4_best, CH4_best_lgr3_december)
    
    # Save updated CH4 best flux results
    write_csv(updated_CH4_best, "../../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")
    
    cat("CH4 best flux results updated:", nrow(updated_CH4_best), "total rows\n")
  } else {
    # If original doesn't exist, just save December data
    write_csv(CH4_best_lgr3_december, "../../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")
    cat("Created new CH4 best flux results file\n")
  }
}

# =============================================================================
# UPDATE FINAL COMPLETE DATASET
# =============================================================================

cat("Updating final complete dataset...\n")

# Load the existing complete dataset
existing_dataset <- read_csv("../../../data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv")

# Prepare December CO2 results with CO2_ prefix
co2_results_december <- CO2_best_lgr3_december %>%
  rename_with(~ paste0("CO2_", .), -UniqueID)

# Prepare December CH4 results if they exist
if(exists("CH4_best_lgr3_december")) {
  ch4_results_december <- CH4_best_lgr3_december %>%
    rename_with(~ paste0("CH4_", .), -UniqueID)
}

# Remove existing December 8th data from the dataset
cat("Removing existing December 8th data from final dataset...\n")
non_december_8_dataset <- existing_dataset %>%
  filter(!grepl(december_8_pattern, UniqueID))

cat("Rows before December 8th removal:", nrow(existing_dataset), "\n")
cat("Rows after December 8th removal:", nrow(non_december_8_dataset), "\n")
cat("December 8th rows removed:", nrow(existing_dataset) - nrow(non_december_8_dataset), "\n")

# Prepare new December 14th data to add back
december_base_data <- lgr3_auxfile_december

# Add December 14th flux results
december_final_data <- december_base_data %>%
  left_join(co2_results_december, by = "UniqueID")

# Add December 14th CH4 results if they exist
if(exists("ch4_results_december")) {
  december_final_data <- december_final_data %>%
    left_join(ch4_results_december, by = "UniqueID")
}

# Combine non-December-8th data with updated December 14th data
updated_final_dataset <- bind_rows(non_december_8_dataset, december_final_data)

cat("Final dataset rows:", nrow(updated_final_dataset), "\n")
cat("December 14th measurements added:", nrow(december_final_data), "\n")

# Save the updated final dataset
write_csv(updated_final_dataset, "../../../data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv")

# =============================================================================
# UPDATE PLOTS (COMBINE EXISTING NON-DECEMBER WITH NEW DECEMBER PLOTS)
# =============================================================================

cat("Creating updated complete plots with new December data...\n")

# Create plots for all updated data (you may want to regenerate all plots to be safe)
# Load the updated complete manual ID data for plotting
updated_manID_complete <- read_csv("../../../data/processed/flux/lgr_manual_identification_results_soil.csv")
updated_CO2_best_complete <- read_csv("../../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv")

# Create complete CO2 plots with updated data
CO2_plots_complete_updated <- flux.plot(
  flux.results = updated_CO2_best_complete,
  dataframe = updated_manID_complete,
  gastype = "CO2dry_ppm",
  shoulder = 30,
  plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
  quality.check = TRUE,
  best.model = TRUE,
  p.val.disp = "round"
)

# Create complete CH4 plots if CH4 data exists
if(file.exists("../../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")) {
  updated_CH4_best_complete <- read_csv("../../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")
  
  CH4_plots_complete_updated <- flux.plot(
    flux.results = updated_CH4_best_complete,
    dataframe = updated_manID_complete,
    gastype = "CH4dry_ppb",
    shoulder = 30,
    plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
    plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
    quality.check = TRUE,
    best.model = TRUE,
    p.val.disp = "round"
  )
  
  all_plots_complete_updated <- c(CO2_plots_complete_updated, CH4_plots_complete_updated)
} else {
  all_plots_complete_updated <- CO2_plots_complete_updated
}

# Save updated complete plots (overwriting the original)
flux2pdf(
  plot.list = all_plots_complete_updated,
  outfile = "../../../data/processed/flux/LGR3_flux_plots_complete_soil.pdf",
  width = 11.6,
  height = 8.2
)

cat("Updated complete plots saved successfully\n")

# Also save December-specific results for reference
write_csv(CO2_best_lgr3_december, "../../../data/processed/flux/CO2_best_flux_december_lgr_results_soil.csv")
if(exists("CH4_best_lgr3_december")) {
  write_csv(CH4_best_lgr3_december, "../../../data/processed/flux/CH4_best_flux_december_lgr_results_soil.csv")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== DECEMBER LGR3 PROCESSING COMPLETE ===\n")
cat("ALL ORIGINAL OUTPUT FILES UPDATED with new December data\n")
cat("\nFiles updated:\n")
cat("- data/processed/flux/lgr_manual_identification_results_soil.csv (UPDATED)\n")
cat("- data/processed/flux/CO2_flux_lgr_results_soil.csv (UPDATED)\n")
cat("- data/processed/flux/CO2_best_flux_lgr_results_soil.csv (UPDATED)\n")
if(file.exists("../../../data/processed/flux/CH4_flux_lgr_results_soil.csv")) {
  cat("- data/processed/flux/CH4_flux_lgr_results_soil.csv (UPDATED)\n")
}
if(file.exists("../../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")) {
  cat("- data/processed/flux/CH4_best_flux_lgr_results_soil.csv (UPDATED)\n")
}
cat("- data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv (UPDATED)\n")
cat("- data/processed/flux/LGR3_flux_plots_complete_soil.pdf (UPDATED with new December plots)\n")

cat("\nAdditional December-specific files created:\n")
cat("- data/processed/flux/CO2_best_flux_december_lgr_results_soil.csv (December only)\n")
if(exists("CH4_best_lgr3_december")) {
  cat("- data/processed/flux/CH4_best_flux_december_lgr_results_soil.csv (December only)\n")
}
cat("- data/processed/flux/LGR3_flux_plots_december_soil.pdf (December plots only)\n")

cat("\nDecember summary statistics:\n")
cat("Total December measurements processed:", nrow(december_final_data), "\n")
cat("Clean December measurements:", quality_summary_december$clean_measurements, "\n")
cat("Flagged December measurements:", quality_summary_december$flagged_measurements, "\n")

cat("\n*** ALL ORIGINAL OUTPUT FILES NOW CONTAIN UPDATED DECEMBER DATA ***\n")
cat("December 8th data removed and December 14th data added to ALL output files.\n")
cat("December-only LGR3 flux processing workflow complete!\n")