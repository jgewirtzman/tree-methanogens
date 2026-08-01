# ==============================================================================
# Semirigid Soil Flux Calculation (goFlux)
# ==============================================================================
# Purpose: Main semirigid soil flux calculation using goFlux package on LGR3
#   data.
#
# Pipeline stage: 01 Flux Processing (semirigid)
# Run after: 03_prep_soil_auxfile.R
#
# Inputs:
#   - Raw LGR3 soil data files
#   - auxfile_goFlux_soilflux_with_weather.csv (from 03_prep_soil_auxfile.R)
#
# Outputs:
#   - semirigid_tree_final_complete_dataset_soil.csv
#   - Soil flux result CSVs
# ==============================================================================

library(goFlux)
library(dplyr)
library(readr)

# =============================================================================
# DEFINE CORRECT CHAMBER GEOMETRY UPFRONT
# =============================================================================

cat("=== CHAMBER GEOMETRY CONFIGURATION ===\n")

# UPDATE THESE VALUES WITH YOUR CORRECT CHAMBER SPECIFICATIONS
CORRECT_CHAMBER_SURFACE_AREA_CM2 <- 507.7  # Based on 25.43 cm interior diameter
CORRECT_CHAMBER_VOLUME_L <- 7.53  # UPDATE THIS TO YOUR CORRECT VALUE!

# Calculate system volumes
correct_tubing_volume_cm3 <- pi * (1/16)^2 * 12 * 12 * 16.387  # ≈ 18.6 cm³
correct_system_volume_cm3 <- 70  # Analyzer volume
correct_total_system_volume_cm3 <- CORRECT_CHAMBER_VOLUME_L * 1000 + 
  correct_tubing_volume_cm3 + 
  correct_system_volume_cm3
CORRECT_VTOT_L <- correct_total_system_volume_cm3 / 1000  # Convert to L

cat("Chamber Geometry:\n")
cat("- Surface Area:", CORRECT_CHAMBER_SURFACE_AREA_CM2, "cm²\n")
cat("- Chamber Volume:", CORRECT_CHAMBER_VOLUME_L, "L\n")
cat("- Total System Volume (Vtot):", round(CORRECT_VTOT_L, 3), "L\n\n")

# =============================================================================
# STEP 1: LOAD AND PREPARE LGR3 DATA
# =============================================================================

cat("=== STEP 1: LOADING LGR3 DATA ===\n")

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
# STEP 2: PREPARE AUXFILE AND CREATE OBSERVATION WINDOWS
# =============================================================================

cat("\n=== STEP 2: PREPARING AUXFILE ===\n")

# Check timezone
lgr3_data_tz <- attr(lgr3_data_complete$POSIX.time, "tzone")
if(is.null(lgr3_data_tz)) lgr3_data_tz <- ""
cat("LGR3 data timezone:", ifelse(lgr3_data_tz == "", "UTC (default)", lgr3_data_tz), "\n")

# Load auxfile WITH CORRECT VOLUME
lgr3_auxfile <- read_csv('../../../data/processed/flux/auxfile_goFlux_soilflux_with_weather.csv', show_col_types = FALSE) %>%
  mutate(
    start.time = as.POSIXct(start.time, tz = ifelse(lgr3_data_tz == "", "UTC", lgr3_data_tz)),
    # UPDATE VOLUME TO CORRECT VALUE
    Vtot = CORRECT_VTOT_L,
    Area = CORRECT_CHAMBER_SURFACE_AREA_CM2
  )

cat("Auxfile loaded:", nrow(lgr3_auxfile), "measurements\n")
cat("Using correct Vtot:", round(CORRECT_VTOT_L, 3), "L\n")

# Create observation windows
cat("Creating observation windows...\n")
ow.lgr3.complete <- obs.win(
  inputfile = lgr3_data_complete,
  auxfile = lgr3_auxfile,
  gastype = "CO2dry_ppm",
  obs.length = 300,
  shoulder = 300
)

# Check window quality
window_sizes <- sapply(ow.lgr3.complete, nrow)
good_windows <- which(window_sizes >= 30)

# Assign names
unique_ids <- sapply(ow.lgr3.complete, function(x) unique(x$UniqueID)[1])
names(ow.lgr3.complete) <- unique_ids

cat("Total observation windows:", length(ow.lgr3.complete), "\n")
cat("Good windows (>=30 obs):", length(good_windows), "\n")

# =============================================================================
# STEP 3: MANUAL IDENTIFICATION
# =============================================================================

cat("\n=== STEP 3: MANUAL IDENTIFICATION ===\n")

# Check if manual ID already exists
if(file.exists("../../../data/processed/flux/lgr_manual_identification_results_soil.csv")) {
  cat("Found existing manual identification file!\n")
  cat("Do you want to:\n")
  cat("1. Use existing manual identification (skip click.peak2)\n")
  cat("2. Redo manual identification (run click.peak2 again)\n")
  cat("Enter 1 or 2: ")
  
  choice <- readline()
  
  if(choice == "1") {
    cat("Loading existing manual identification...\n")
    manID.lgr3 <- read_csv("../../../data/processed/flux/lgr_manual_identification_results_soil.csv", 
                           show_col_types = FALSE)
    
    # Update volume in loaded data
    manID.lgr3 <- manID.lgr3 %>%
      mutate(
        Vtot = CORRECT_VTOT_L,
        Area = CORRECT_CHAMBER_SURFACE_AREA_CM2
      )
    
    cat("Loaded", nrow(manID.lgr3), "measurements with corrected volume\n")
  } else {
    cat("Proceeding with new manual identification...\n")
    # Run manual identification in batches
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
      
      if(batch_num < num_batches) {
        cat("Press Enter to continue to next batch...\n")
        readline()
      }
    }
    
    # Combine all batches
    successful_batches <- manID_batches[!sapply(manID_batches, is.null)]
    manID.lgr3 <- do.call(rbind, successful_batches)
    
    cat("\nManual identification complete! Total:", nrow(manID.lgr3), "measurements\n")
  }
} else {
  # No existing file, must run click.peak2
  cat("No existing manual identification found. Running click.peak2...\n")
  
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
    
    if(batch_num < num_batches) {
      cat("Press Enter to continue to next batch...\n")
      readline()
    }
  }
  
  # Combine all batches
  successful_batches <- manID_batches[!sapply(manID_batches, is.null)]
  manID.lgr3 <- do.call(rbind, successful_batches)
  
  cat("\nManual identification complete! Total:", nrow(manID.lgr3), "measurements\n")
}

# Save manual identification results (with correct volume)
write_csv(manID.lgr3, "../../../data/processed/flux/lgr_manual_identification_results_soil.csv")

# =============================================================================
# STEP 4: FLUX CALCULATIONS
# =============================================================================

cat("\n=== STEP 4: FLUX CALCULATIONS ===\n")

# Calculate CO2 fluxes
cat("Calculating CO2 fluxes with Vtot =", round(CORRECT_VTOT_L, 3), "L...\n")
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
  cat("CH4 flux calculation complete:", nrow(CH4_flux_lgr3), "measurements\n")
} else {
  cat("CH4dry_ppb column not found - skipping CH4 flux calculation\n")
}

# Save flux results
write_csv(CO2_flux_lgr3, "../../../data/processed/flux/CO2_flux_lgr_results_soil.csv")
if(exists("CH4_flux_lgr3")) {
  write_csv(CH4_flux_lgr3, "../../../data/processed/flux/CH4_flux_lgr_results_soil.csv")
}

cat("CO2 flux calculation complete:", nrow(CO2_flux_lgr3), "measurements\n\n")

# =============================================================================
# STEP 5: BEST FLUX ANALYSIS
# =============================================================================

cat("=== STEP 5: BEST FLUX ANALYSIS ===\n")

# Run best.flux on CO2 results
CO2_best_lgr3 <- best.flux(
  flux.result = CO2_flux_lgr3,
  criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", 
               "MDF", "nb.obs", "intercept", "p-value"),
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
    criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", 
                 "MDF", "nb.obs", "intercept", "p-value"),
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
cat("- Clean measurements:", quality_summary$clean_measurements, "\n")
cat("- Flagged measurements:", quality_summary$flagged_measurements, "\n")
cat("  - Below MDF:", quality_summary$below_MDF, "\n")
cat("  - High g-factor:", quality_summary$high_g_factor, "\n")  
cat("  - Poor fit:", quality_summary$poor_fit, "\n\n")

# Save best flux results
write_csv(CO2_best_lgr3, "../../../data/processed/flux/CO2_best_flux_lgr_results_soil.csv")
if(exists("CH4_best_lgr3")) {
  write_csv(CH4_best_lgr3, "../../../data/processed/flux/CH4_best_flux_lgr_results_soil.csv")
}

# =============================================================================
# STEP 6: CREATE PLOTS
# =============================================================================

cat("=== STEP 6: CREATING PLOTS ===\n")

# Fix any missing model values before plotting
CO2_best_lgr3 <- CO2_best_lgr3 %>%
  mutate(
    model = case_when(
      is.na(model) ~ "LM",
      TRUE ~ model
    )
  )

# Create CO2 flux plots
CO2_plots_lgr3 <- flux.plot(
  flux.results = CO2_best_lgr3,
  dataframe = manID.lgr3,
  gastype = "CO2dry_ppm",
  shoulder = 30,
  plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
  quality.check = TRUE,
  best.model = TRUE,
  p.val.disp = "round"
)

# Create CH4 plots if available
if(exists("CH4_best_lgr3")) {
  CH4_best_lgr3 <- CH4_best_lgr3 %>%
    mutate(
      model = case_when(
        is.na(model) ~ "LM",
        TRUE ~ model
      )
    )
  
  CH4_plots_lgr3 <- flux.plot(
    flux.results = CH4_best_lgr3,
    dataframe = manID.lgr3,
    gastype = "CH4dry_ppb",
    shoulder = 30,
    plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
    plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
    quality.check = TRUE,
    best.model = TRUE,
    p.val.disp = "round"
  )
  
  all_plots_lgr3 <- c(CO2_plots_lgr3, CH4_plots_lgr3)
} else {
  all_plots_lgr3 <- CO2_plots_lgr3
}

# Save plots to PDF
flux2pdf(
  plot.list = all_plots_lgr3,
  outfile = "../../../data/processed/flux/LGR3_flux_plots_complete_soil.pdf",
  width = 11.6,
  height = 8.2
)
cat("Plots saved successfully\n\n")



library(gridExtra)
library(grid)

# Custom function to replace flux2pdf
save_flux_plots <- function(plot_list, outfile, width = 11.6, height = 8.2) {
  
  pdf(outfile, width = width, height = height)
  
  # Process plots one by one
  for(i in 1:length(plot_list)) {
    
    # Try to extract UniqueID from the plot if available
    # This might not work with the new structure, so we'll use index as backup
    plot_id <- paste("Plot", i, "of", length(plot_list))
    
    # Create title and footnote
    title <- textGrob(plot_id, gp = gpar(fontsize = 16))
    footnote <- textGrob(paste("page", i, "of", length(plot_list)), 
                         gp = gpar(fontface = 3, fontsize = 12), hjust = 1)
    
    # Arrange plot with title and footnote
    grid.arrange(plot_list[[i]], 
                 top = title, 
                 bottom = footnote)
  }
  
  dev.off()
  cat("PDF saved to:", outfile, "\n")
}

# Use the custom function
save_flux_plots(
  plot_list = all_plots_lgr3,
  outfile = "../../../outputs/figures/LGR3_flux_plots_complete_soil.pdf",
  width = 11.6,
  height = 8.2
)

# =============================================================================
# STEP 7: CREATE FINAL DATASET
# =============================================================================

cat("=== STEP 7: CREATING FINAL DATASET ===\n")

# Load original auxfile data
original_data <- read_csv('../../../data/processed/flux/auxfile_goFlux_soilflux_with_weather_formatted.csv',
                          show_col_types = FALSE)

# Update volume in auxfile for consistency
original_data <- original_data %>%
  mutate(
    Vtot = CORRECT_VTOT_L,
    Area = CORRECT_CHAMBER_SURFACE_AREA_CM2
  )

# Add CO2_ prefix to all flux columns except UniqueID
co2_results <- CO2_best_lgr3 %>%
  rename_with(~ paste0("CO2_", .), -UniqueID)

# Add CH4 results if they exist
if(exists("CH4_best_lgr3")) {
  ch4_results <- CH4_best_lgr3 %>%
    rename_with(~ paste0("CH4_", .), -UniqueID)
}

# Merge everything together
final_dataset <- original_data %>%
  left_join(co2_results, by = "UniqueID")

# Add CH4 if it exists
if(exists("ch4_results")) {
  final_dataset <- final_dataset %>%
    left_join(ch4_results, by = "UniqueID")
}

# Save the final dataset
write_csv(final_dataset, "../../../data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv")

# Also update the auxfiles with correct volume
write_csv(original_data %>% select(-contains("CO2_"), -contains("CH4_")), 
          "auxfile_goFlux_soilflux_with_weather_formatted.csv")
write.table(original_data %>% select(-contains("CO2_"), -contains("CH4_"), -start.time_formatted), 
            "auxfile_goFlux_soilflux_with_weather.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
# STEP 8: FLUX SUMMARY STATISTICS
# =============================================================================

cat("=== STEP 8: FLUX SUMMARY STATISTICS ===\n")

# CO2 flux statistics
co2_stats <- CO2_best_lgr3 %>%
  summarise(
    n = n(),
    mean_flux = mean(flux.term, na.rm = TRUE),
    median_flux = median(flux.term, na.rm = TRUE),
    sd_flux = sd(flux.term, na.rm = TRUE),
    min_flux = min(flux.term, na.rm = TRUE),
    max_flux = max(flux.term, na.rm = TRUE),
    q25 = quantile(flux.term, 0.25, na.rm = TRUE),
    q75 = quantile(flux.term, 0.75, na.rm = TRUE)
  )

cat("\nCO2 Flux Statistics (nmol/m²/s):\n")
cat("- Mean:", round(co2_stats$mean_flux, 2), "\n")
cat("- Median:", round(co2_stats$median_flux, 2), "\n")
cat("- SD:", round(co2_stats$sd_flux, 2), "\n")
cat("- Range:", round(co2_stats$min_flux, 2), "to", round(co2_stats$max_flux, 2), "\n")
cat("- IQR:", round(co2_stats$q25, 2), "to", round(co2_stats$q75, 2), "\n")

# CH4 flux statistics if available
if(exists("CH4_best_lgr3")) {
  ch4_stats <- CH4_best_lgr3 %>%
    summarise(
      n = n(),
      mean_flux = mean(flux.term, na.rm = TRUE),
      median_flux = median(flux.term, na.rm = TRUE),
      sd_flux = sd(flux.term, na.rm = TRUE),
      min_flux = min(flux.term, na.rm = TRUE),
      max_flux = max(flux.term, na.rm = TRUE),
      q25 = quantile(flux.term, 0.25, na.rm = TRUE),
      q75 = quantile(flux.term, 0.75, na.rm = TRUE)
    )
  
  cat("\nCH4 Flux Statistics (nmol/m²/s):\n")
  cat("- Mean:", round(ch4_stats$mean_flux, 2), "\n")
  cat("- Median:", round(ch4_stats$median_flux, 2), "\n")
  cat("- SD:", round(ch4_stats$sd_flux, 2), "\n")
  cat("- Range:", round(ch4_stats$min_flux, 2), "to", round(ch4_stats$max_flux, 2), "\n")
  cat("- IQR:", round(ch4_stats$q25, 2), "to", round(ch4_stats$q75, 2), "\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== LGR3 PROCESSING COMPLETE ===\n")
cat("Chamber specifications applied:\n")
cat("- Surface Area:", CORRECT_CHAMBER_SURFACE_AREA_CM2, "cm²\n")
cat("- Total System Volume:", round(CORRECT_VTOT_L, 3), "L\n\n")

cat("Final dataset created with", nrow(final_dataset), "rows and", ncol(final_dataset), "columns\n")

cat("\nAll files saved with correct volume (", round(CORRECT_VTOT_L, 3), "L):\n")
cat("- data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv\n")
cat("- data/processed/flux/CO2_best_flux_lgr_results_soil.csv\n")
if(exists("CH4_best_lgr3")) {
  cat("- data/processed/flux/CH4_best_flux_lgr_results_soil.csv\n")
}
cat("- data/processed/flux/LGR3_flux_plots_complete_soil.pdf\n")
cat("- auxfile_goFlux_soilflux_with_weather.txt\n")
cat("- auxfile_goFlux_soilflux_with_weather_formatted.csv\n")
cat("- data/processed/flux/lgr_manual_identification_results_soil.csv\n")

cat("\nSummary:\n")
cat("- Total measurements processed:", nrow(final_dataset), "\n")
cat("- Clean measurements:", quality_summary$clean_measurements, "\n")
cat("- Flagged measurements:", quality_summary$flagged_measurements, "\n")

cat("\n✓ All flux values calculated with correct volume from the start\n")
cat("✓ All files saved with consistent, correct chamber specifications\n")
