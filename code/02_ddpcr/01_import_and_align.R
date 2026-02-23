# ==============================================================================
# Import and Align ddPCR Data
# ==============================================================================
# Purpose: Imports raw ddPCR data files, aligns with sample metadata, and
#   calculates gene abundances (mcrA, pmoA, mmoxY).
#
# Pipeline stage: 01 Molecular Processing
# Run after: None (first script in ddPCR pipeline)
#
# Inputs:
#   - Raw ddPCR CSVs (from data/raw/ddpcr/)
#   - Guide/metadata files
#
# Outputs:
#   - processed_ddpcr_data.csv
# ==============================================================================

library(plyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(wesanderson)
theme_set(theme_bw(base_size = 12))

# Define file paths
base_path <- '../../data/raw/ddpcr/'
raw_data_path <- file.path(base_path, 'Raw Data')

# ============================================================================
# 1. IMPORT AND PROCESS METADATA
# ============================================================================

import_metadata <- function(base_path) {
  # Import guide CSV files
  guide_files <- c('inner_guide.csv', 'outer_guide.csv', 'mineral_guide.csv', 'organic_guide.csv')
  
  metadata_list <- lapply(guide_files, function(file) {
    file_path <- file.path(base_path, file)
    if (file.exists(file_path)) {
      read.csv(file_path, skip = 1, header = TRUE, stringsAsFactors = TRUE)
    } else {
      warning(paste("File not found:", file_path))
      return(NULL)
    }
  })
  
  # Combine all metadata
  metadata <- do.call(rbind, metadata_list[!sapply(metadata_list, is.null)])
  
  # Extract sample information from Inner.Core.Sample.ID
  metadata$species <- substr(metadata$Inner.Core.Sample.ID, 1, 4)
  
  # Extract core type and material information
  core_info <- sub("\nL.*", "", metadata$Inner.Core.Sample.ID)
  core_info <- sub(".*\n", "", core_info)
  metadata$material <- sub("_.*", "", core_info)
  metadata$core_type <- sub(".*_", "", core_info)
  
  # Extract lab ID and create plate identifier
  metadata$lab_id <- sub(".*: ", "", metadata$Inner.Core.Sample.ID)
  metadata$plate_identifier <- paste(
    sub("Plate ", "", metadata$Extraction.Plate.ID), 
    metadata$Plate.Position, 
    sep = ""
  )
  
  return(metadata)
}

# ============================================================================
# 2. IMPORT AND PROCESS DDPCR DATA
# ============================================================================

import_ddpcr_data <- function(raw_data_path) {
  # Get list of CSV files, excluding system files
  data_files <- list.files(raw_data_path, 
                           pattern = "\\.csv$", 
                           full.names = TRUE)
  
  # Filter out system files
  data_files <- data_files[!grepl("Icon|\\._|DS_Store|\\.tmp", basename(data_files))]
  
  cat("Found", length(data_files), "data files to process\n")
  
  # Initialize empty dataframe
  ddpcr_data <- data.frame()
  failed_files <- c()
  
  # Process each file
  for (i in seq_along(data_files)) {
    file_path <- data_files[i]
    file_name <- basename(file_path)
    
    tryCatch({
      cat("Processing file", i, "of", length(data_files), ":", file_name, "\n")
      
      # Read the file
      data_file <- read.csv(file_path, stringsAsFactors = FALSE)
      
      # Determine analysis type from filename
      if (str_detect(file_name, "loose")) {
        data_file$analysis_type <- "loose"
      } else if (str_detect(file_name, "strict")) {
        data_file$analysis_type <- "strict"
      } else {
        data_file$analysis_type <- "none"
      }
      
      # Standardize to first 23 columns plus analysis_type
      n_cols <- min(23, ncol(data_file) - 1)  # -1 for analysis_type column
      data_file <- data_file[, c(1:n_cols, ncol(data_file))]
      
      # Append to main dataset
      ddpcr_data <- rbind(ddpcr_data, data_file)
      
    }, error = function(e) {
      warning(paste("Failed to import:", file_name, "- Error:", e$message))
      failed_files <<- c(failed_files, file_path)
    })
  }
  
  if (length(failed_files) > 0) {
    cat("Failed to import", length(failed_files), "files:\n")
    cat(paste("-", basename(failed_files), collapse = "\n"), "\n")
  }
  
  return(ddpcr_data)
}

# ============================================================================
# 3. CLEAN AND STANDARDIZE DDPCR DATA
# ============================================================================

clean_ddpcr_data <- function(ddpcr_data) {
  # Standardize well format (e.g., "A01" -> "A1")
  if ("Well" %in% names(ddpcr_data)) {
    ddpcr_data$Well <- paste0(
      substr(ddpcr_data$Well, 1, 1),
      as.numeric(substr(ddpcr_data$Well, 2, 4))
    )
  }
  
  # Clean sample description formatting
  if ("Sample.description.1" %in% names(ddpcr_data)) {
    ddpcr_data$Sample.description.1 <- gsub("_", " ", ddpcr_data$Sample.description.1)
  }
  
  # Create plate identifier for merging
  if ("Sample.description.1" %in% names(ddpcr_data) && "Well" %in% names(ddpcr_data)) {
    ddpcr_data$plate_identifier <- paste0(
      sub("Plate ", "", ddpcr_data$Sample.description.1),
      ddpcr_data$Well
    )
  }
  
  # Standardize target names
  if ("Target" %in% names(ddpcr_data)) {
    ddpcr_data$Target <- tolower(ddpcr_data$Target)
  }
  
  return(ddpcr_data)
}

# ============================================================================
# 4. MERGE DATA AND APPLY QUALITY FILTERS
# ============================================================================

merge_and_filter <- function(metadata, ddpcr_data) {
  # Merge metadata with ddPCR data
  merged_data <- merge(metadata, ddpcr_data, by = "plate_identifier", all.x = FALSE, all.y = FALSE)
  
  cat("Merged data dimensions:", nrow(merged_data), "rows x", ncol(merged_data), "columns\n")
  
  # Remove standards and blank wells
  if ("Sample.description.2" %in% names(merged_data)) {
    merged_data <- merged_data[!grepl("std", merged_data$Sample.description.2, ignore.case = TRUE), ]
    cat("After removing standards:", nrow(merged_data), "rows\n")
  }
  
  if ("Blank.Well." %in% names(merged_data)) {
    merged_data <- merged_data[!grepl("yes", merged_data$Blank.Well., ignore.case = TRUE), ]
    cat("After removing blanks:", nrow(merged_data), "rows\n")
  }
  
  # Quality filter: Remove samples with < 5000 accepted droplets
  if ("Accepted.Droplets" %in% names(merged_data)) {
    low_quality <- merged_data$Accepted.Droplets < 5000
    cat("Found", sum(low_quality, na.rm = TRUE), "samples with < 5000 droplets\n")
    merged_data <- merged_data[!low_quality, ]
    cat("After quality filtering:", nrow(merged_data), "rows\n")
  }
  
  return(merged_data)
}

# ============================================================================
# 5. MAIN PROCESSING FUNCTION
# ============================================================================

process_ddpcr_data <- function(base_path, raw_data_path) {
  cat("=== ddPCR Data Processing Pipeline ===\n\n")
  
  # Step 1: Import metadata
  cat("Step 1: Importing metadata...\n")
  metadata <- import_metadata(base_path)
  cat("Imported metadata for", nrow(metadata), "samples\n\n")
  
  # Step 2: Import ddPCR data
  cat("Step 2: Importing ddPCR data...\n")
  ddpcr_data <- import_ddpcr_data(raw_data_path)
  cat("Imported", nrow(ddpcr_data), "ddPCR measurements\n\n")
  
  # Step 3: Clean and standardize data
  cat("Step 3: Cleaning and standardizing data...\n")
  ddpcr_data <- clean_ddpcr_data(ddpcr_data)
  cat("Data cleaning complete\n\n")
  
  # Step 4: Merge and filter
  cat("Step 4: Merging data and applying quality filters...\n")
  final_data <- merge_and_filter(metadata, ddpcr_data)
  cat("Final dataset:", nrow(final_data), "rows x", ncol(final_data), "columns\n\n")
  
  # Summary statistics
  cat("=== Processing Summary ===\n")
  if ("Target" %in% names(final_data)) {
    cat("Targets analyzed:", paste(unique(final_data$Target), collapse = ", "), "\n")
  }
  if ("analysis_type" %in% names(final_data)) {
    cat("Analysis types:", paste(unique(final_data$analysis_type), collapse = ", "), "\n")
  }
  if ("species" %in% names(final_data)) {
    cat("Species represented:", length(unique(final_data$species)), "\n")
  }
  
  return(final_data)
}

# ============================================================================
# 6. RUN THE PIPELINE
# ============================================================================

# Execute the processing pipeline
final_data <- process_ddpcr_data(base_path, raw_data_path)

# Optional: Save the processed data
 write.csv(final_data, "../../data/processed/molecular/processed_ddpcr_data.csv", row.names = FALSE)
 