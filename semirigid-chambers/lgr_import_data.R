#if (!require("devtools", quietly = TRUE)) install.packages("devtools")
#try(detach("package:goFlux", unload = TRUE), silent = TRUE)
#devtools::install_github("Qepanna/goFlux")
library(goFlux)
library(dplyr)
#library(purrr)
library(readxl)
library(openxlsx)

# Set your data path
data_path <- '/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/LGR files'

cat("=== EXTRACTING ZIP FILES ===\n")

# Step 1: Extract all remaining zip files
zip_files <- list.files(data_path, recursive = TRUE, pattern = "\\.zip$", full.names = TRUE)
cat("Found", length(zip_files), "zip files to extract\n")

# Extract each zip file to its parent directory
for(zip_file in zip_files) {
  # Extract to the same directory as the zip file
  extract_dir <- dirname(zip_file)
  tryCatch({
    unzip(zip_file, exdir = extract_dir, overwrite = TRUE)
    cat("Extracted:", basename(zip_file), "\n")
  }, error = function(e) {
    cat("Error extracting", basename(zip_file), ":", e$message, "\n")
  })
}

cat("\n=== COLLECTING ALL DATA FILES ===\n")

# Step 2: Get all txt files after extraction
txt_files <- list.files(data_path, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE)

# Step 3: Filter for data files (non-empty, excluding files in nested .txt folders)
clean_txt_files <- txt_files[!grepl("/.*\\.txt/", txt_files)]
data_files <- clean_txt_files[file.size(clean_txt_files) > 0]

cat("Total .txt files after extraction:", length(txt_files), "\n")
cat("Clean .txt files:", length(clean_txt_files), "\n") 
cat("Non-empty data files:", length(data_files), "\n")

# Step 4: Verify we're getting the "f" files (flux measurements)
f_files <- data_files[grepl("_f[0-9]+\\.txt$", data_files)]
cat("Flux measurement files (f-files):", length(f_files), "\n")

# Step 5: Create temporary directory with only data files
temp_data_path <- tempfile(pattern = "lgr_complete_")
dir.create(temp_data_path, recursive = TRUE)

# Copy all data files
file.copy(data_files, temp_data_path)
cat("Copied", length(data_files), "files to temporary directory\n")

cat("\n=== IMPORTING WITH GOFLUX ===\n")

# Step 6: Import using goFlux
lgr_data <- import2RData(
  path = temp_data_path,
  instrument = "UGGA",
  date.format = "mdy",        # MM/DD/YYYY format from LGR
  timezone = "UTC",
  keep_all = FALSE,           
  prec = c(0.35, 0.9, 200),  # MGGA GLA131 precision: CO2, CH4, H2O
  merge = TRUE               
)

# Step 7: Clean up
unlink(temp_data_path, recursive = TRUE)

# Results summary
cat("\n=== IMPORT COMPLETE ===\n")
cat("Final dataset dimensions:", dim(lgr_data), "\n")
cat("Date range:", range(lgr_data$POSIX.Date, na.rm = TRUE), "\n")
cat("Unique measurement days:", length(unique(as.Date(lgr_data$POSIX.Date))), "\n")

# Show first few rows
head(lgr_data)

# Optional: Save the complete dataset
# write.csv(lgr_data, "lgr_gla131_complete_dataset.csv", row.names = FALSE)


# Check date coverage using the DATE column
cat("=== DATE COVERAGE ANALYSIS ===\n\n")

# Convert DATE column to proper Date format
unique_dates <- unique(as.Date(lgr_data$DATE))
unique_dates <- sort(unique_dates[!is.na(unique_dates)])

cat("Total unique measurement days:", length(unique_dates), "\n")
cat("Date range:", min(unique_dates), "to", max(unique_dates), "\n")
cat("Time span:", as.numeric(max(unique_dates) - min(unique_dates)), "days\n\n")

# Show all unique dates
cat("All measurement dates:\n")
print(unique_dates)

# Check data points per day
cat("\n=== DATA DENSITY ===\n")
daily_counts <- table(as.Date(lgr_data$DATE))
cat("Total measurements:", nrow(lgr_data), "\n")
cat("Average measurements per day:", round(mean(daily_counts), 1), "\n")
cat("Range of measurements per day:", min(daily_counts), "to", max(daily_counts), "\n")

# Show measurements per day
cat("\nMeasurements per day:\n")
print(daily_counts)

