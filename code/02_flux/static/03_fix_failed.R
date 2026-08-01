# ==============================================================================
# Reprocess Failed Flux Measurements
# ==============================================================================
# Purpose: Reprocesses failed or problematic flux measurements with manual
#   corrections.
#
# Pipeline stage: 01 Flux Processing (static)
# Run after: 02_goflux_trees_2021.R
#
# Inputs:
#   - lgr_manual_identification_results.csv
#
# Outputs:
#   - lgr_manual_identification_results_final.csv
# ==============================================================================

cat("=== REPROCESSING SPECIFIC FAILED MEASUREMENTS ===\n")

# Define the specific measurements that need reprocessing
failed_measurements <- c(
  "20210803_Lowland_YB12_200_3_182300",
  "20210803_Lowland_YB11_200_3_172500", 
  "20210803_Lowland_YB11_125_3_171900",
  "20210803_Lowland_YB12_50_3_181900",
  "20210806_Lowland_AB 6prime_125_6_102700"
)

cat("Failed measurements to reprocess:\n")
print(failed_measurements)

# Find the corresponding window indices in ow.lgr3.complete
failed_window_indices <- which(names(ow.lgr3.complete) %in% failed_measurements)

cat("\nFound window indices:", failed_window_indices, "\n")
cat("Number of windows found:", length(failed_window_indices), "\n")

# Check if all measurements were found
if(length(failed_window_indices) != length(failed_measurements)) {
  missing <- failed_measurements[!failed_measurements %in% names(ow.lgr3.complete)]
  cat("WARNING: Could not find windows for:", missing, "\n")
}

if(length(failed_window_indices) > 0) {
  cat("\nReprocessing", length(failed_window_indices), "failed measurements...\n")
  
  # Reprocess just the failed measurements
  manID_reprocessed <- click.peak2(
    ow.lgr3.complete,
    seq = failed_window_indices,
    gastype = "CO2dry_ppm",
    sleep = 3,
    plot.lim = c(200, 5000),  # Adjust if needed
    save.plots = "../../../data/processed/flux/lgr3_reprocessed_failed_plots"
  )
  
  # Close graphics devices
  while(dev.cur() > 1) dev.off()
  
  cat("Reprocessing complete!\n")
  cat("New results:\n")
  print(manID_reprocessed)
  
  # Create summary of reprocessed results
  reprocessed_summary <- manID_reprocessed %>%
    group_by(UniqueID) %>%
    summarise(
      nb.obs = n(),
      start.time = first(start.time_corr),
      end.time = last(end.time_corr),
      .groups = "drop"
    )
  
  cat("\nSummary of reprocessed measurements:\n")
  print(reprocessed_summary)
  
  # Now update your existing batch results
  cat("\nUpdating batch results...\n")
  
  # Remove the old failed measurements from existing batches
  for(i in 1:length(manID_batches)) {
    if(!is.null(manID_batches[[i]]) && nrow(manID_batches[[i]]) > 0) {
      # Remove failed measurements from this batch
      manID_batches[[i]] <- manID_batches[[i]] %>%
        filter(!UniqueID %in% failed_measurements)
    }
  }
  
  # Add the reprocessed data as a new batch
  manID_batches[[length(manID_batches) + 1]] <- manID_reprocessed
  
  # Recreate the complete manID.lgr3 object
  successful_batches <- manID_batches[!sapply(manID_batches, is.null)]
  successful_batches <- successful_batches[sapply(successful_batches, nrow) > 0]
  
  manID.lgr3 <- do.call(rbind, successful_batches)
  
  cat("\nUpdated manID.lgr3 object:\n")
  cat("Total rows:", nrow(manID.lgr3), "\n")
  
  # Create final summary statistics
  final_summary <- manID.lgr3 %>%
    group_by(UniqueID) %>%
    summarise(
      nb.obs = n(),
      start.time = first(start.time_corr),
      end.time = last(end.time_corr),
      .groups = "drop"
    )
  
  cat("\nFinal quality check:\n")
  quality_stats <- final_summary %>%
    summarise(
      total_measurements = n(),
      mean_obs = mean(nb.obs),
      min_obs = min(nb.obs),
      max_obs = max(nb.obs),
      measurements_under_60 = sum(nb.obs < 60),
      measurements_under_30 = sum(nb.obs < 30),
      measurements_under_20 = sum(nb.obs < 20)
    )
  
  print(quality_stats)
  
  # Save updated results
  write_csv(manID.lgr3, "../../../data/processed/flux/lgr_manual_identification_results_final.csv")
  write_csv(final_summary, "../../../data/processed/flux/lgr_manual_identification_summary_final.csv")
  
  cat("\n=== REPROCESSING COMPLETE ===\n")
  cat("✓ Updated manID.lgr3 object with corrected measurements\n")
  cat("✓ Files saved:\n")
  cat("  - flux_code/lgr_manual_identification_results_final.csv\n")
  cat("  - flux_code/lgr_manual_identification_summary_final.csv\n")
  cat("✓ Ready to continue with flux calculations!\n")
  
} else {
  cat("No window indices found for the specified measurements.\n")
  cat("Please check that the UniqueID values are correct.\n")
}

