# ==============================================================================
# Interactive Peak Identification Utility
# ==============================================================================
# Purpose: Interactive utility for manual peak identification in flux time
#   series curves using click.peak2.
# ==============================================================================

library(goFlux)
library(dplyr)
library(purrr)

# Check the combined dataset
cat("\nDataset summary:\n")
cat("Dimensions:", nrow(lgr_data), "rows,", ncol(lgr_data), "columns\n")
cat("Columns:", paste(names(lgr_data), collapse = ", "), "\n")

# Check for gas measurement columns
gas_columns <- names(lgr_data)[grepl("dry_pp|CH4|CO2", names(lgr_data))]
cat("Gas measurement columns:", paste(gas_columns, collapse = ", "), "\n")

# Check time range
if("POSIX.time" %in% names(lgr_data)) {
  time_range <- range(lgr_data$POSIX.time, na.rm = TRUE)
  cat("UGGA data time range:", as.character(time_range), "\n")
} else {
  stop("No POSIX.time column found in UGGA data!")
}

# Load auxfile (use the basic version for now)
cat("\n=== LOADING AUXFILE ===\n")
auxfile <- read.csv("../../../data/processed/flux/auxfile_goFlux_with_weather_formatted_datetime.csv")

# Convert start.time to POSIXct with proper handling
auxfile$start.time <- as.POSIXct(auxfile$start.time_formatted, tz = "UTC")

cat("Auxfile loaded with", nrow(auxfile), "measurements\n")
cat("Auxfile time range:", as.character(range(auxfile$start.time, na.rm = TRUE)), "\n")

# Check for time overlap
ugga_start <- min(lgr_data$POSIX.time, na.rm = TRUE)
ugga_end <- max(lgr_data$POSIX.time, na.rm = TRUE)
aux_start <- min(auxfile$start.time, na.rm = TRUE)
aux_end <- max(auxfile$start.time, na.rm = TRUE)

overlap_start <- max(ugga_start, aux_start)
overlap_end <- min(ugga_end, aux_end)

if(overlap_start < overlap_end) {
  cat("Time overlap confirmed from", as.character(overlap_start), "to", as.character(overlap_end), "\n")
} else {
  stop("No time overlap between UGGA data and auxfile! Check your date/time formats.")
}

# Check required columns
required_cols <- c("UniqueID", "start.time", "Area", "Vtot", "Tcham", "Pcham")
missing_cols <- required_cols[!required_cols %in% names(auxfile)]
if(length(missing_cols) > 0) {
  stop("Missing required columns in auxfile: ", paste(missing_cols, collapse = ", "))
} else {
  cat("All required auxfile columns present\n")
}

cat("\n=== CREATING OBSERVATION WINDOWS ===\n")

# Create observation windows for methane measurements
ow.data <- obs.win(
  inputfile = lgr_data,           # Combined UGGA data
  auxfile = auxfile,               # Your auxfile
  gastype = "CH4dry_ppb",          # Methane measurements
  obs.length = 600,                # 5 minutes chamber closure
  shoulder = 300                   # 2 minutes buffer before/after
)

cat("Created observation windows for", length(ow.data), "measurements\n")

# Check if observation windows contain data
if(length(ow.data) > 0) {
  first_window <- ow.data[[1]]
  cat("First window contains", nrow(first_window), "data points\n")
  
  if(nrow(first_window) > 0) {
    ch4_range <- range(first_window$CH4dry_ppb, na.rm = TRUE)
    cat("CH4 range in first window:", round(ch4_range, 1), "ppb\n")
    
    if(any(is.finite(ch4_range))) {
      cat("Observation windows look good!\n")
    } else {
      stop("CH4 data in observation windows is not finite. Check time matching.")
    }
  } else {
    stop("Observation windows are empty. Check time overlap between auxfile and UGGA data.")
  }
} else {
  stop("No observation windows created. Check auxfile format and time overlap.")
}

# Batch size recommendation
if(length(ow.data) > 20) {
  cat("You have", length(ow.data), "measurements.\n")
  cat("Recommended: Process in batches of 10-20 to avoid mistakes.\n")
  batch_size <- 10
} else {
  batch_size <- length(ow.data)
}

cat("\n=== SETTING UP GRAPHICS DEVICE ===\n")

# Save default graphics device
default.device <- getOption("device")

# Set appropriate graphics device for interactive clicking
if(.Platform$OS.type == "windows") {
  options(device = "windows")
  cat("Graphics device set to 'windows'\n")
} else if(Sys.info()["sysname"] == "Darwin") {  # macOS
  options(device = "quartz") 
  cat("Graphics device set to 'quartz'\n")
} else {  # Linux/Unix
  options(device = "X11")
  cat("Graphics device set to 'X11'\n")
}

cat("\n=== INTERACTIVE MANUAL IDENTIFICATION ===\n")
cat("Starting click.peak2 for methane measurements.\n")
cat("For each plot that appears:\n")
cat("1. Click on the START point of the methane change\n")
cat("2. Click on the END point of the measurement\n")
cat("3. Plot will close and move to next measurement\n")
cat("\nProcessing first", batch_size, "measurements...\n")

# Run click.peak2 for first batch
manID.data <- click.peak2(
  ow.list = ow.data,
  gastype = "CH4dry_ppb",          # Methane measurements
  sleep = 3,                       # 3 seconds between plots
  plot.lim = c(1800, 2500),        # Y-axis limits for methane (adjust if needed)
  seq = seq(1, batch_size),        # First batch
  warn.length = 60,                # Warning for short measurements
  save.plots = "../../../outputs/figures/methane_manual_ID_batch1"  # Save plots as PDF
)

cat("First batch complete! Processed", nrow(manID.data), "measurements\n")

# Process additional batches if needed
if(length(ow.data) > batch_size) {
  remaining <- length(ow.data) - batch_size
  cat("You have", remaining, "more measurements to process.\n")
  cat("Continue with next batch? (y/n): ")
  
  continue_response <- readline()
  
  if(tolower(continue_response) == "y") {
    # Process second batch
    batch2_end <- min(batch_size * 2, length(ow.data))
    
    manID.data.2 <- click.peak2(
      ow.list = ow.data,
      gastype = "CH4dry_ppb",
      seq = seq(batch_size + 1, batch2_end),
      save.plots = "../../../outputs/figures/methane_manual_ID_batch2"
    )
    
    # Combine batches
    manID.data <- rbind(manID.data, manID.data.2)
    cat("Combined batches. Total processed:", nrow(manID.data), "measurements\n")
    
    # Continue with additional batches if needed
    if(batch2_end < length(ow.data)) {
      cat("Continue with remaining", length(ow.data) - batch2_end, "measurements? (y/n): ")
      continue_more <- readline()
      
      if(tolower(continue_more) == "y") {
        cat("Process remaining measurements in batches of 10-20 using:\n")
        cat("manID.data.3 <- click.peak2(ow.data, gastype='CH4dry_ppb', seq=seq(", batch2_end + 1, ",", length(ow.data), "))\n")
        cat("manID.data <- rbind(manID.data, manID.data.3)\n")
      }
    }
  }
}

# Restore graphics device
options(device = default.device)
cat("Graphics device restored\n")

cat("\n=== SAVING RESULTS ===\n")

# Save manual identification results
write.csv(manID.data, "methane_manual_identification.csv", row.names = FALSE)
cat("Results saved as: methane_manual_identification.csv\n")

# Quality check summary
cat("\n=== QUALITY CHECK ===\n")
cat("Total measurements processed:", nrow(manID.data), "\n")

if("flag" %in% names(manID.data)) {
  flagged <- sum(grepl("flag", manID.data$flag, ignore.case = TRUE), na.rm = TRUE)
  cat("Flagged measurements:", flagged, "\n")
}

if("obs.length_corr" %in% names(manID.data)) {
  avg_length <- round(mean(manID.data$obs.length_corr, na.rm = TRUE), 1)
  cat("Average measurement length:", avg_length, "seconds\n")
}

cat("\n=== NEXT STEPS ===\n")
cat("1. Review saved plots: methane_manual_ID_batch*.pdf\n")
cat("2. Check results in: methane_manual_identification.csv\n")
cat("3. Run flux calculations:\n")
cat("   CH4_flux <- goFlux(manID.data, 'CH4dry_ppb')\n")
cat("   write.csv(CH4_flux, 'methane_flux_results.csv', row.names = FALSE)\n")

cat("\nManual identification workflow complete!\n")