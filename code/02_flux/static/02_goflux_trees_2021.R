# ==============================================================================
# Static Chamber Tree Flux Calculation - 2021 (goFlux)
# ==============================================================================
# Purpose: Main static chamber tree flux calculation for 2021 data using goFlux.
#
# Pipeline stage: 01 Flux Processing (static)
# Run after: 01_prep_auxfile.R
#
# Inputs:
#   - goflux_auxfile.csv (from 01_prep_auxfile.R)
#   - Raw LGR3 files (from data/raw/lgr/static_2021/)
#
# Outputs:
#   - methanogen_tree_flux_complete_dataset.csv (KEY OUTPUT)
#   - Flux result CSVs
# ==============================================================================

library(goFlux)
library(dplyr)
library(readr)

# =============================================================================
# STEP 1: LOAD AND PREPARE LGR3 DATA
# =============================================================================

cat("=== STEP 1: LOADING LGR3 DATA ===\n")

# Set data path
data_path <- '../../../data/raw/lgr/static_2021'

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
# STEP 2: PREPARE AUXFILE AND CREATE OBSERVATION WINDOWS
# =============================================================================

cat("\n=== STEP 2: PREPARING AUXFILE ===\n")

# Check timezone and load auxfile
lgr3_data_tz <- attr(lgr3_data_complete$POSIX.time, "tzone")
if(is.null(lgr3_data_tz)) lgr3_data_tz <- ""
cat("LGR3 data timezone:", ifelse(lgr3_data_tz == "", "UTC (default)", lgr3_data_tz), "\n")

# Load LGR3 auxfile
lgr3_auxfile <- read_csv('../../../data/processed/flux/goflux_auxfile.csv') %>%
  mutate(start.time = as.POSIXct(start.time, tz = ifelse(lgr3_data_tz == "", "UTC", lgr3_data_tz)))

cat("Auxfile loaded:", nrow(lgr3_auxfile), "measurements\n")

# Create observation windows
cat("Creating observation windows...\n")
ow.lgr3.complete <- obs.win(
  inputfile = lgr3_data_complete,
  auxfile = lgr3_auxfile,
  gastype = "CO2dry_ppm",
  obs.length = 300,
  shoulder = 1200
)

# Check window quality
window_sizes <- sapply(ow.lgr3.complete, nrow)
good_windows <- which(window_sizes >= 30)

# Assign names
unique_ids <- sapply(ow.lgr3.complete, function(x) unique(x$UniqueID)[1])
names(ow.lgr3.complete) <- unique_ids

cat("Total observation windows:", length(ow.lgr3.complete), "\n")
cat("Good windows (>=60 obs):", length(good_windows), "\n")

# =============================================================================
# STEP 3: MANUAL IDENTIFICATION IN BATCHES
# =============================================================================

cat("\n=== STEP 3: MANUAL IDENTIFICATION ===\n")

# Process in batches of 20
batch_size <- 20
total_measurements <- length(good_windows)
num_batches <- ceiling(total_measurements / batch_size)

cat("Processing", total_measurements, "measurements in", num_batches, "batches\n")

manID_batches <- list()

for(batch_num in 1:num_batches) {
  start_idx <- (batch_num - 1) * batch_size + 1
  end_idx <- min(batch_num * batch_size, total_measurements)
  batch_windows <- good_windows[start_idx:end_idx]
  
  cat("\n=== Processing Batch", batch_num, "of", num_batches, "===\n")
  cat("Measurements", start_idx, "to", end_idx, "\n")
  
  # Process this batch
  manID_batch <- click.peak2(
    ow.lgr3.complete,
    seq = batch_windows,
    gastype = "CO2dry_ppm",
    sleep = 3,
    plot.lim = c(200, 5000),
    save.plots = paste0("../../../data/processed/flux/lgr3_batch_", batch_num, "_plots")
  )
  
  # Store the batch result
  manID_batches[[batch_num]] <- manID_batch
  
  # Force close graphics devices
  while(dev.cur() > 1) dev.off()
  
  cat("Batch", batch_num, "complete. Processed", nrow(manID_batch), "measurements\n")
  cat("Press Enter to continue to next batch...\n")
  readline()
}

# 
# # First, close any stuck graphics devices
# while(dev.cur() > 1) dev.off()
# 
# # The error occurred at batch 11, measurements 201-220
# # You successfully completed batches 1-10, so let's continue from where it failed
# 
# # Check which measurements were already processed
# completed_batches <- 10  # You completed batches 1-10
# completed_measurements <- completed_batches * batch_size
# 
# # Continue from batch 11, but skip the problematic measurement
# # The error occurred around measurement 94 (Mar_23_38_BL60_stem)
# 
# # Let's identify and skip problematic measurements
# problem_windows <- which(window_sizes < 10)  # Very small windows
# cat("Problematic windows with <10 observations:", length(problem_windows), "\n")
# 
# # Filter out problematic windows from good_windows
# good_windows_filtered <- good_windows[!good_windows %in% problem_windows]
# 
# # Continue processing from where you left off
# remaining_start <- completed_measurements + 1
# remaining_measurements <- length(good_windows_filtered) - completed_measurements
# 
# if(remaining_measurements > 0) {
#   remaining_batches <- ceiling(remaining_measurements / batch_size)
#   
#   cat("Continuing from measurement", remaining_start, "\n")
#   cat("Remaining measurements:", remaining_measurements, "in", remaining_batches, "batches\n")
#   
#   for(batch_num in (completed_batches + 1):(completed_batches + remaining_batches)) {
#     start_idx <- (batch_num - completed_batches - 1) * batch_size + remaining_start
#     end_idx <- min(start_idx + batch_size - 1, length(good_windows_filtered))
#     batch_windows <- good_windows_filtered[start_idx:end_idx]
#     
#     cat("\n=== Processing Batch", batch_num, "===\n")
#     cat("Measurements", start_idx, "to", end_idx, "\n")
#     
#     tryCatch({
#       manID_batch <- click.peak2(
#         ow.lgr3.complete,
#         seq = batch_windows,
#         gastype = "CO2dry_ppm",
#         sleep = 3,
#         plot.lim = c(200, 1200),
#         save.plots = paste0("../../../data/processed/flux/lgr3_batch_", batch_num, "_plots")
#       )
#       
#       manID_batches[[batch_num]] <- manID_batch
#       cat("Batch", batch_num, "complete. Processed", nrow(manID_batch), "measurements\n")
#       
#     }, error = function(e) {
#       cat("Error in batch", batch_num, ":", e$message, "\n")
#       cat("Skipping problematic measurement and continuing...\n")
#     })
#     
#     while(dev.cur() > 1) dev.off()
#     cat("Press Enter to continue...\n")
#     readline()
#   }
# }



# Combine all successful batches
successful_batches <- manID_batches[!sapply(manID_batches, is.null)]
manID.lgr3 <- do.call(rbind, successful_batches)

cat("Processing recovered! Total measurements:", nrow(manID.lgr3), "\n")


# Combine all batches
manID.lgr3 <- do.call(rbind, manID_batches)

cat("\nManual identification complete! Total:", nrow(manID.lgr3), "measurements\n")

# Save manual identification results
write_csv(manID.lgr3, "../../../data/processed/flux/lgr_manual_identification_results.csv")

# =============================================================================
# STEP 4: FLUX CALCULATIONS
# =============================================================================

cat("\n=== STEP 4: FLUX CALCULATIONS ===\n")

# Calculate CO2 fluxes
cat("Calculating CO2 fluxes...\n")
CO2_flux_lgr3 <- goFlux(
  dataframe = manID.lgr3,
  gastype = "CO2dry_ppm",
  H2O_col = "H2O_ppm",
  warn.length = 60
)

# Calculate CH4 fluxes if available
if("CH4dry_ppb" %in% names(manID.lgr3)) {
  cat("Calculating CH4 fluxes...\n")
  CH4_flux_lgr3 <- goFlux(
    dataframe = manID.lgr3,
    gastype = "CH4dry_ppb",
    H2O_col = "H2O_ppm",
    warn.length = 60
  )
} else {
  cat("CH4dry_ppb column not found - skipping CH4 flux calculation\n")
}

# Save flux results
write_csv(CO2_flux_lgr3, "../../../data/processed/flux/CO2_flux_lgr_results.csv")
if(exists("CH4_flux_lgr3")) {
  write_csv(CH4_flux_lgr3, "../../../data/processed/flux/CH4_flux_lgr_results.csv")
}

cat("CO2 flux calculation complete:", nrow(CO2_flux_lgr3), "measurements\n")

# =============================================================================
# STEP 5: BEST FLUX ANALYSIS
# =============================================================================

cat("\n=== STEP 5: BEST FLUX ANALYSIS ===\n")

# Run best.flux on CO2 results
CO2_best_lgr3 <- best.flux(
  flux.result = CO2_flux_lgr3,
  criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", "MDF", "nb.obs", "intercept", "p-value"),
  intercept.lim = NULL,
  g.limit = 2,
  p.val = 0.05,
  k.ratio = 1,
  warn.length = 60
)

# Run best.flux on CH4 results if available
if(exists("CH4_flux_lgr3")) {
  CH4_best_lgr3 <- best.flux(
    flux.result = CH4_flux_lgr3,
    criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", "MDF", "nb.obs", "intercept", "p-value"),
    intercept.lim = NULL,
    g.limit = 2,
    p.val = 0.05,
    k.ratio = 1,
    warn.length = 60
  )
}

# Quality summary
quality_summary <- CO2_best_lgr3 %>%
  summarise(
    clean_measurements = sum(quality.check == "OK", na.rm = TRUE),
    flagged_measurements = sum(quality.check != "OK", na.rm = TRUE),
    below_MDF = sum(grepl("MDF", quality.check), na.rm = TRUE),
    high_g_factor = sum(grepl("g.factor", quality.check), na.rm = TRUE),
    poor_fit = sum(grepl("MAE|RMSE|SE", quality.check), na.rm = TRUE)
  )

cat("Quality Assessment:\n")
cat("Clean measurements:", quality_summary$clean_measurements, "\n")
cat("Flagged measurements:", quality_summary$flagged_measurements, "\n")

# Save best flux results
write_csv(CO2_best_lgr3, "../../../data/processed/flux/CO2_best_flux_lgr_results.csv")
if(exists("CH4_best_lgr3")) {
  write_csv(CH4_best_lgr3, "../../../data/processed/flux/CH4_best_flux_lgr_results.csv")
}

# =============================================================================
# STEP 6: CREATE PLOTS
# =============================================================================
# 
# cat("\n=== STEP 6: CREATING PLOTS ===\n")
# 
# # Create CO2 flux plots
# CO2_plots_lgr3 <- flux.plot(
#   flux.results = CO2_best_lgr3,
#   dataframe = manID.lgr3,
#   gastype = "CO2dry_ppm",
#   shoulder = 30,
#   plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
#   plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
#   quality.check = TRUE,
#   best.model = TRUE,
#   p.val.disp = "round"
# )
# 
# # Create CH4 plots if available
# if(exists("CH4_best_lgr3")) {
#   CH4_plots_lgr3 <- flux.plot(
#     flux.results = CH4_best_lgr3,
#     dataframe = manID.lgr3,
#     gastype = "CH4dry_ppb",
#     shoulder = 30,
#     plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
#     plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
#     quality.check = TRUE,
#     best.model = TRUE,
#     p.val.disp = "round"
#   )
#   


# Check for missing values in the model column
cat("Checking for missing model values...\n")
missing_models <- sum(is.na(CO2_best_lgr3$model))
cat("Number of missing model values:", missing_models, "\n")

# Display rows with missing models
if(missing_models > 0) {
  cat("Rows with missing model values:\n")
  print(CO2_best_lgr3[is.na(CO2_best_lgr3$model), c("UniqueID", "model", "LM.flux", "HM.flux")])
}

# Option 1: Remove rows with missing model values
CO2_best_lgr3_clean <- CO2_best_lgr3[!is.na(CO2_best_lgr3$model), ]
cat("Original rows:", nrow(CO2_best_lgr3), "-> Clean rows:", nrow(CO2_best_lgr3_clean), "\n")

# Option 2: Alternative approach - set best.model = FALSE to avoid the issue
cat("\n=== CREATING PLOTS WITH CLEAN DATA ===\n")

# Method 1: Use cleaned data with best.model = TRUE
CO2_plots_lgr3_method1 <- flux.plot(
  flux.results = CO2_best_lgr3_clean,
  dataframe = manID.lgr3,
  gastype = "CO2dry_ppm",
  shoulder = 30,
  plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
  quality.check = TRUE,
  best.model = TRUE,
  p.val.disp = "round"
)

cat("Method 1 completed successfully!\n")

# Method 2: Use original data but disable best.model highlighting
CO2_plots_lgr3_method2 <- flux.plot(
  flux.results = CO2_best_lgr3,
  dataframe = manID.lgr3,
  gastype = "CO2dry_ppm",
  shoulder = 30,
  plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
  quality.check = TRUE,
  best.model = FALSE,  # Disable best model highlighting
  p.val.disp = "round"
)

cat("Method 2 completed successfully!\n")

# Create CH4 plots if data exists
if(exists("CH4_best_lgr3")) {
  cat("\n=== CREATING CH4 PLOTS ===\n")
  
  # Check CH4 data for missing models too
  missing_ch4_models <- sum(is.na(CH4_best_lgr3$model))
  cat("CH4 missing model values:", missing_ch4_models, "\n")
  
  if(missing_ch4_models > 0) {
    CH4_best_lgr3_clean <- CH4_best_lgr3[!is.na(CH4_best_lgr3$model), ]
  } else {
    CH4_best_lgr3_clean <- CH4_best_lgr3
  }
  
  CH4_plots_lgr3 <- flux.plot(
    flux.results = CH4_best_lgr3_clean,
    dataframe = manID.lgr3,
    gastype = "CH4dry_ppb",
    shoulder = 30,
    plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
    plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
    quality.check = TRUE,
    best.model = ifelse(missing_ch4_models > 0, FALSE, TRUE),
    p.val.disp = "round"
  )
  
  # Combine all plots
  all_plots_lgr3 <- c(CO2_plots_lgr3_method1, CH4_plots_lgr3)
  cat("Combined CO2 and CH4 plots created!\n")
} else {
  all_plots_lgr3 <- CO2_plots_lgr3_method1
  cat("Only CO2 plots created (no CH4 data found)\n")
}

# Display first few plots
cat("\n=== DISPLAYING PLOTS ===\n")
cat("Total number of plots created:", length(all_plots_lgr3), "\n")

# Show first 3 plots as examples
if(length(all_plots_lgr3) >= 3) {
  print(all_plots_lgr3[[1]])
  print(all_plots_lgr3[[2]]) 
  print(all_plots_lgr3[[3]])
} else {
  # Show all available plots
  for(i in 1:length(all_plots_lgr3)) {
    print(all_plots_lgr3[[i]])
  }
}

cat("Plot creation completed successfully!\n")

#   # Combine plots
#   all_plots_lgr3 <- c(CO2_plots_lgr3, CH4_plots_lgr3)
# } else {
#   all_plots_lgr3 <- CO2_plots_lgr3
# }
# 
# Save plots to PDF
flux2pdf(
  plot.list = all_plots_lgr3,
  outfile = "../../../data/processed/flux/LGR3_flux_plots_complete_ymf_methanogens.pdf",
  width = 11.6,
  height = 8.2
)
cat("Plots saved successfully\n")

# 
# 
# # Check for missing values in the model column
# cat("Checking for missing model values...\n")
# missing_models <- sum(is.na(CO2_best_lgr3$model))
# cat("Missing model values:", missing_models, "\n")
# 
# # Check what model values exist
# table(CO2_best_lgr3$model, useNA = "always")
# 
# # Fix the missing model values
# CO2_best_lgr3_fixed <- CO2_best_lgr3 %>%
#   mutate(
#     model = case_when(
#       is.na(model) ~ "LM",  # Default to linear model for missing values
#       TRUE ~ model
#     )
#   )
# 
# # Verify the fix
# cat("After fix - Missing model values:", sum(is.na(CO2_best_lgr3_fixed$model)), "\n")
# 
# # Now try creating plots with the fixed data
# CO2_plots_lgr3 <- flux.plot(
#   flux.results = CO2_best_lgr3_fixed,
#   dataframe = manID.lgr3,
#   gastype = "CO2dry_ppm",
#   shoulder = 30,
#   plot.legend = c("MAE", "AICc"),  # Simplified legend to avoid potential issues
#   plot.display = c("MDF", "prec"),  # Simplified display
#   quality.check = TRUE,
#   best.model = TRUE,
#   p.val.disp = "round"
# )
# 
# # If CH4 exists, apply the same fix
# if(exists("CH4_best_lgr3")) {
#   CH4_best_lgr3_fixed <- CH4_best_lgr3 %>%
#     mutate(
#       model = case_when(
#         is.na(model) ~ "LM",
#         TRUE ~ model
#       )
#     )
#   
#   CH4_plots_lgr3 <- flux.plot(
#     flux.results = CH4_best_lgr3_fixed,
#     dataframe = manID.lgr3,
#     gastype = "CH4dry_ppb",
#     shoulder = 30,
#     plot.legend = c("MAE", "AICc"),
#     plot.display = c("MDF", "prec"),
#     quality.check = TRUE,
#     best.model = TRUE,
#     p.val.disp = "round"
#   )
#   
#   all_plots_lgr3 <- c(CO2_plots_lgr3, CH4_plots_lgr3)
# } else {
#   all_plots_lgr3 <- CO2_plots_lgr3
# }
# 
# # Continue with saving plots
# flux2pdf(
#   plot.list = all_plots_lgr3,
#   outfile = "../../../data/processed/flux/LGR3_flux_plots_complete.pdf",
#   width = 11.6,
#   height = 8.2
# )
# 
# cat("Plots created successfully!\n")

# =============================================================================
# STEP 7: CREATE FINAL DATASET
# =============================================================================

cat("\n=== STEP 7: CREATING FINAL DATASET ===\n")

# Load original tree data
original_data <- read_csv('../../../data/processed/flux/goflux_auxfile.csv')

# Filter for LGR3 only
lgr3_original <- original_data #%>% filter(analyzer_id == "LGR3")

# Add CO2_ prefix to all flux columns except UniqueID
co2_results <- CO2_best_lgr3 %>%
  rename_with(~ paste0("CO2_", .), -UniqueID)

# Add CH4 results if they exist
if(exists("CH4_best_lgr3")) {
  ch4_results <- CH4_best_lgr3 %>%
    rename_with(~ paste0("CH4_", .), -UniqueID)
}

# Merge everything together
final_dataset <- lgr3_original %>%
  left_join(co2_results, by = "UniqueID")

# Add CH4 if it exists
if(exists("ch4_results")) {
  final_dataset <- final_dataset %>%
    left_join(ch4_results, by = "UniqueID")
}

# Save the final dataset
write_csv(final_dataset, "../../../data/processed/flux/methanogen_tree_flux_complete_dataset.csv")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== LGR3 PROCESSING COMPLETE ===\n")
cat("Final dataset created with", nrow(final_dataset), "rows and", ncol(final_dataset), "columns\n")
cat("\nFiles created:\n")
cat("- data/processed/flux/LGR3_final_complete_dataset.csv (complete final dataset)\n")
cat("- data/processed/flux/CO2_best_flux_lgr3_results.csv (CO2 flux results)\n")
if(exists("CH4_best_lgr3")) {
  cat("- data/processed/flux/CH4_best_flux_lgr3_results.csv (CH4 flux results)\n")
}
cat("- data/processed/flux/LGR3_flux_plots_complete.pdf (all plots)\n")
cat("- Multiple batch plot files for manual review\n")

cat("\nSummary statistics:\n")
cat("Total measurements processed:", nrow(final_dataset), "\n")
cat("Clean measurements:", quality_summary$clean_measurements, "\n")
cat("Flagged measurements:", quality_summary$flagged_measurements, "\n")

cat("\nLGR3 flux processing workflow complete!\n")
