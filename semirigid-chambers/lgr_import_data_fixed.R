# Load required packages
library(goFlux)
library(dplyr)
library(purrr)

# Set the path to your LGR files
file_path <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/LGR files"

# Create a directory for unzipped files
unzipped_path <- file.path(dirname(file_path), "LGR files unzipped")
dir.create(unzipped_path, showWarnings = FALSE, recursive = TRUE)

# Function to unzip all zip files while maintaining folder structure
unzip_all_files <- function(base_path, output_path) {
  # Find all zip files recursively
  zip_files <- list.files(base_path, pattern = "\\.zip$", 
                          full.names = TRUE, recursive = TRUE)
  
  cat("Found", length(zip_files), "zip files to extract\n")
  
  for (zip_file in zip_files) {
    # Get relative path to maintain folder structure
    rel_path <- gsub(paste0("^", base_path, "/"), "", zip_file)
    rel_dir <- dirname(rel_path)
    
    # Create corresponding directory in output path
    extract_dir <- file.path(output_path, rel_dir)
    dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Extract the zip file
    tryCatch({
      unzip(zip_file, exdir = extract_dir, overwrite = TRUE)
      cat("Extracted:", basename(zip_file), "to", extract_dir, "\n")
    }, error = function(e) {
      cat("Error extracting", basename(zip_file), ":", e$message, "\n")
    })
  }
}

cat("=== EXTRACTING ZIP FILES ===\n")
# Unzip all files
cat("Starting to unzip files...\n")
unzip_all_files(file_path, unzipped_path)
cat("Unzipping complete!\n\n")

cat("=== COLLECTING AND DEDUPLICATING FILES ===\n")
# Get all txt files after extraction
txt_files <- list.files(unzipped_path, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
cat("Found", length(txt_files), ".txt files after extraction\n")

# Filter and clean files
clean_txt_files <- txt_files[!grepl("/.*\\.txt/", txt_files)]  # Remove nested txt files
non_empty_files <- clean_txt_files[file.size(clean_txt_files) > 0]  # Remove empty files

cat("Clean .txt files:", length(clean_txt_files), "\n") 
cat("Non-empty files:", length(non_empty_files), "\n")

# Initialize lgr_data as NULL
lgr_data <- NULL

if (length(non_empty_files) > 0) {
  cat("Sample file paths:\n")
  cat(head(non_empty_files, 5), sep = "\n")
  
  # DEDUPLICATION: Remove duplicate files by basename and size
  cat("\n=== DEDUPLICATING FILES ===\n")
  
  # Create a data frame with file info for deduplication
  file_info <- data.frame(
    full_path = non_empty_files,
    basename = basename(non_empty_files),
    size = file.size(non_empty_files),
    stringsAsFactors = FALSE
  )
  
  # Group by basename and size, keep only the first occurrence of each unique file
  unique_files <- file_info %>%
    group_by(basename, size) %>%
    slice(1) %>%  # Keep only the first occurrence
    ungroup() %>%
    pull(full_path)
  
  cat("Files before deduplication:", length(non_empty_files), "\n")
  cat("Files after deduplication:", length(unique_files), "\n")
  cat("Duplicates removed:", length(non_empty_files) - length(unique_files), "\n")
  
  # Show which files were duplicated (optional diagnostic)
  duplicated_files <- file_info %>%
    group_by(basename, size) %>%
    filter(n() > 1) %>%
    arrange(basename, full_path)
  
  if(nrow(duplicated_files) > 0) {
    cat("\nDuplicated files found:\n")
    print(duplicated_files[,c("basename", "full_path")], row.names = FALSE)
  }
  
  # Focus on flux measurement files (f-files)
  f_files <- unique_files[grepl("_f[0-9]+\\.txt$", unique_files)]
  cat("Flux measurement files (f-files):", length(f_files), "\n")
  
  # Let's examine the first file to understand its structure
  cat("\nChecking first file structure...\n")
  first_file <- unique_files[1]
  cat("First file:", basename(first_file), "\n")
  
  # Read first few lines to understand date format
  first_lines <- readLines(first_file, n = 10)
  cat("First few lines of the file:\n")
  cat(first_lines[1:5], sep = "\n")
  
  cat("\n=== TESTING IMPORT FORMAT ===\n")
  # Test different date formats to find the correct one
  test_formats <- c("dmy", "ymd", "mdy")
  successful_format <- NULL
  
  for (format in test_formats) {
    cat("Testing date format:", format, "\n")
    tryCatch({
      test_import <- import.UGGA(
        inputfile = first_file,
        date.format = format,
        timezone = "UTC",
        prec = c(0.35, 0.9, 200),
        save = FALSE
      )
      cat("SUCCESS: Single file import worked with", format, "format!\n")
      successful_format <- format
      break
    }, error = function(e) {
      cat("FAILED:", format, "format -", e$message, "\n")
    })
  }
  
  # Create flat directory for import (always use this method since it's more reliable)
  cat("\n=== CREATING FLAT DIRECTORY FOR IMPORT ===\n")
  all_flat_dir <- file.path(dirname(file_path), "All_LGR_files_flat")
  dir.create(all_flat_dir, showWarnings = FALSE)
  
  cat("Copying all", length(unique_files), "unique files to flat directory...\n")
  
  # Copy all files with unique names to avoid conflicts
  for (i in seq_along(unique_files)) {
    original_file <- unique_files[i]
    # Create unique filename using date and original name
    date_folder <- basename(dirname(original_file))
    new_name <- paste0(date_folder, "_", basename(original_file))
    
    file.copy(original_file, file.path(all_flat_dir, new_name))
    
    # Progress indicator
    if (i %% 100 == 0) {
      cat("Copied", i, "of", length(unique_files), "files...\n")
    }
  }
  
  cat("\n=== IMPORTING DATA ===\n")
  # Try importing from flat directory with different formats
  if (is.null(successful_format)) {
    # If single file test failed, try all formats for batch import
    formats_to_try <- test_formats
  } else {
    # If single file test worked, try that format first, then others
    formats_to_try <- c(successful_format, test_formats[test_formats != successful_format])
  }
  
  for (format in formats_to_try) {
    cat("Trying batch import with", format, "format...\n")
    tryCatch({
      lgr_data <- import2RData(
        path = all_flat_dir,
        instrument = "UGGA",
        date.format = format,
        timezone = "UTC",
        prec = c(0.35, 0.9, 200),
        keep_all = FALSE,
        merge = TRUE
      )
      cat("SUCCESS: Batch import worked with", format, "format!\n")
      break
    }, error = function(e) {
      cat("FAILED: Batch import with", format, "format:", e$message, "\n")
      lgr_data <- NULL
    })
  }
  
  # If batch import failed, try to recover data from RData files
  if (is.null(lgr_data)) {
    cat("\nBatch import failed. Checking for RData files...\n")
    
    # Check if RData files were created despite the error
    rdata_path <- file.path(all_flat_dir, "RData")
    if (dir.exists(rdata_path)) {
      rdata_files <- list.files(rdata_path, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)
      cat("Found", length(rdata_files), "RData files in", rdata_path, "\n")
      
      if (length(rdata_files) > 0) {
        cat("Loading RData files manually...\n")
        all_data <- list()
        
        for (i in seq_along(rdata_files)) {
          tryCatch({
            # Create a new environment to load into
            temp_env <- new.env()
            load(rdata_files[i], envir = temp_env)
            
            # Get all objects from the environment
            obj_names <- ls(envir = temp_env)
            if (length(obj_names) > 0) {
              # Take the first data frame-like object
              for (obj_name in obj_names) {
                obj <- get(obj_name, envir = temp_env)
                if (is.data.frame(obj) && nrow(obj) > 0) {
                  all_data[[length(all_data) + 1]] <- obj
                  cat("Loaded file", i, "(", basename(rdata_files[i]), ") with", nrow(obj), "rows\n")
                  break
                }
              }
            }
          }, error = function(e) {
            cat("Error loading", basename(rdata_files[i]), ":", e$message, "\n")
          })
        }
        
        if (length(all_data) > 0) {
          # Combine all loaded data
          cat("Combining", length(all_data), "datasets...\n")
          lgr_data <- do.call(rbind, all_data)
          cat("Successfully combined datasets into one with", nrow(lgr_data), "rows\n")
        }
      }
    } else {
      cat("No RData directory found at", rdata_path, "\n")
    }
  }
  
  # Final check if we have data
  if (is.null(lgr_data)) {
    cat("\n=== IMPORT FAILED ===\n")
    cat("ERROR: No data was successfully imported.\n")
    cat("Please check:\n")
    cat("1. File formats and date formats\n")
    cat("2. Whether RData files were created in", file.path(all_flat_dir, "RData"), "\n")
    cat("3. Individual file imports using import.UGGA() function\n")
    stop("Data import failed")
  }
  
  cat("\n=== FINAL DATASET DEDUPLICATION ===\n")
  # Additional deduplication check on the final dataset
  initial_rows <- nrow(lgr_data)
  cat("Initial dataset rows:", initial_rows, "\n")
  
  # Remove duplicates based on timestamp and measurement values
  lgr_data_clean <- lgr_data %>%
    distinct(POSIX.time, CO2dry_ppm, CH4dry_ppb, H2O_ppm, .keep_all = TRUE)
  
  final_rows <- nrow(lgr_data_clean)
  cat("Rows after final deduplication:", final_rows, "\n")
  cat("Final duplicates removed:", initial_rows - final_rows, "\n")
  
  # Use the cleaned dataset
  lgr_data <- lgr_data_clean
  
  cat("\n=== IMPORT COMPLETE ===\n")
  cat("Final dataset dimensions:", dim(lgr_data), "\n")
  
  # Check column names and date columns
  cat("Column names:", paste(colnames(lgr_data), collapse = ", "), "\n")
  
  if ("POSIX.Date" %in% colnames(lgr_data)) {
    date_range <- range(lgr_data$POSIX.Date, na.rm = TRUE)
    cat("Date range:", date_range[1], "to", date_range[2], "\n")
    cat("Unique measurement days:", length(unique(as.Date(lgr_data$POSIX.Date))), "\n")
  } else if ("DATE" %in% colnames(lgr_data)) {
    unique_dates <- unique(as.Date(lgr_data$DATE))
    unique_dates <- sort(unique_dates[!is.na(unique_dates)])
    if (length(unique_dates) > 0) {
      cat("Date range:", min(unique_dates), "to", max(unique_dates), "\n")
      cat("Unique measurement days:", length(unique_dates), "\n")
    }
  }
  
  # Check for any remaining duplicates
  remaining_dupes <- lgr_data %>%
    group_by(POSIX.time, CO2dry_ppm, CH4dry_ppb) %>%
    filter(n() > 1) %>%
    nrow()
  
  cat("Remaining duplicate rows:", remaining_dupes, "\n")
  
  # Show first few rows
  cat("\nFirst few rows of final dataset:\n")
  print(head(lgr_data))
  
  cat("\n=== SAVE CLEAN DATASET ===\n")
  # Save the deduplicated dataset
  output_filename <- "lgr_clean_deduplicated_dataset.csv"
  write.csv(lgr_data, output_filename, row.names = FALSE)
  cat("Clean dataset saved to:", output_filename, "\n")
  
  # Data quality summary
  if ("DATE" %in% colnames(lgr_data)) {
    daily_counts <- table(as.Date(lgr_data$DATE))
    cat("Total measurements:", nrow(lgr_data), "\n")
    cat("Average measurements per day:", round(mean(daily_counts), 1), "\n")
    cat("Range of measurements per day:", min(daily_counts), "to", max(daily_counts), "\n")
  }
  
} else {
  cat("No .txt files found. Let's check what files were extracted:\n")
  all_files <- list.files(unzipped_path, recursive = TRUE, full.names = TRUE)
  cat("Total files found:", length(all_files), "\n")
  if (length(all_files) > 0) {
    cat("Sample files:\n")
    cat(head(all_files, 10), sep = "\n")
  }
}

cat("\n=== CLEANUP ===\n")
# Optional: Clean up temporary directories
cat("Cleaning up temporary directories...\n")
unlink(unzipped_path, recursive = TRUE)
if (exists("all_flat_dir") && dir.exists(all_flat_dir)) {
  unlink(all_flat_dir, recursive = TRUE)
}

cat("Process complete! Your LGR MGGA files have been successfully imported and deduplicated.\n")