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
      unzip(zip_file, exdir = extract_dir)
      cat("Extracted:", basename(zip_file), "to", extract_dir, "\n")
    }, error = function(e) {
      cat("Error extracting", basename(zip_file), ":", e$message, "\n")
    })
  }
}

# Unzip all files
cat("Starting to unzip files...\n")
unzip_all_files(file_path, unzipped_path)
cat("Unzipping complete!\n\n")

# First, let's check what files were actually extracted
cat("Checking extracted files...\n")
txt_files <- list.files(unzipped_path, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
cat("Found", length(txt_files), ".txt files\n")

if (length(txt_files) > 0) {
  cat("Sample file paths:\n")
  cat(head(txt_files, 5), sep = "\n")
  
  # Let's examine the first file to understand its structure
  cat("\nChecking first file structure...\n")
  first_file <- txt_files[1]
  cat("First file:", basename(first_file), "\n")
  
  # Read first few lines to understand date format
  first_lines <- readLines(first_file, n = 10)
  cat("First few lines of the file:\n")
  cat(first_lines[1:5], sep = "\n")
  
  # Since goFlux is having trouble finding files, let's try importing individual files
  # or copy files to a simpler directory structure
  
  # Option 1: Try importing a single file first to test
  cat("\nTrying to import a single file to test...\n")
  tryCatch({
    # Test with single file import
    test_import <- import.UGGA(
      inputfile = first_file,
      date.format = "dmy",  # Start with dmy
      timezone = "UTC",
      prec = c(0.35, 0.9, 200),
      save = FALSE
    )
    cat("Single file import successful with dmy format!\n")
    cat("Columns in imported data:", colnames(test_import), "\n")
    
    # If single file works, try the batch import again
    cat("Single file worked, trying batch import again...\n")
    import2RData(
      path = unzipped_path,
      instrument = "UGGA",
      date.format = "dmy",
      timezone = "UTC",
      prec = c(0.35, 0.9, 200),
      keep_all = FALSE,
      merge = FALSE
    )
    cat("Batch import successful!\n")
    
  }, error = function(e1) {
    cat("dmy format failed, trying ymd...\n")
    tryCatch({
      test_import <- import.UGGA(
        inputfile = first_file,
        date.format = "ymd",
        timezone = "UTC",
        prec = c(0.35, 0.9, 200),
        save = FALSE
      )
      cat("Single file import successful with ymd format!\n")
      
      # Try batch import with ymd
      import2RData(
        path = unzipped_path,
        instrument = "UGGA",
        date.format = "ymd",
        timezone = "UTC",
        prec = c(0.35, 0.9, 200),
        keep_all = FALSE,
        merge = FALSE
      )
      cat("Batch import successful!\n")
      
    }, error = function(e2) {
      cat("Both date formats failed for single file. Error:", e2$message, "\n")
      
      # Option 2: Create a flatter directory structure
      cat("Trying to create a flatter directory structure...\n")
      flat_dir <- file.path(dirname(file_path), "LGR_files_flat")
      dir.create(flat_dir, showWarnings = FALSE)
      
      # Copy a few files to test
      test_files <- head(txt_files, 10)
      for (i in seq_along(test_files)) {
        file.copy(test_files[i], 
                  file.path(flat_dir, paste0("test_", i, ".txt")))
      }
      
      cat("Created flat directory with", length(test_files), "test files\n")
      cat("Try running import2RData on:", flat_dir, "\n")
      
      # Try importing from flat directory
      tryCatch({
        import2RData(
          path = flat_dir,
          instrument = "UGGA",
          date.format = "dmy",
          timezone = "UTC",
          prec = c(0.35, 0.9, 200),
          keep_all = FALSE,
          merge = FALSE
        )
        cat("Flat directory import successful!\n")
      }, error = function(e3) {
        cat("Flat directory import also failed:", e3$message, "\n")
        cat("This suggests there might be an issue with the file format or content.\n")
        cat("You may need to check the actual content and format of the .txt files.\n")
      })
    })
  })
} else {
  cat("No .txt files found. Let's check what files were extracted:\n")
  all_files <- list.files(unzipped_path, recursive = TRUE, full.names = TRUE)
  cat("Total files found:", length(all_files), "\n")
  if (length(all_files) > 0) {
    cat("Sample files:\n")
    cat(head(all_files, 10), sep = "\n")
  }
}

# SUCCESS! The flat directory approach worked. 
# Now let's import ALL your files using this method.

# Create a flat directory for all files
all_flat_dir <- file.path(dirname(file_path), "All_LGR_files_flat")
dir.create(all_flat_dir, showWarnings = FALSE)

cat("Copying all", length(txt_files), "files to flat directory...\n")
cat("This may take a few minutes...\n")

# Copy all files with unique names to avoid conflicts
for (i in seq_along(txt_files)) {
  original_file <- txt_files[i]
  # Create unique filename using date and original name
  date_folder <- basename(dirname(original_file))
  new_name <- paste0(date_folder, "_", basename(original_file))
  
  file.copy(original_file, file.path(all_flat_dir, new_name))
  
  # Progress indicator
  if (i %% 100 == 0) {
    cat("Copied", i, "of", length(txt_files), "files...\n")
  }
}

cat("All files copied! Now importing...\n")

# Import all files from flat directory
import2RData(
  path = all_flat_dir,
  instrument = "UGGA",
  date.format = "dmy",  # We confirmed this works
  timezone = "UTC",
  prec = c(0.35, 0.9, 200),  # MGGA precision values
  keep_all = FALSE,
  merge = TRUE  # Set to TRUE if you want one combined dataframe
)

cat("Import complete! All files have been processed.\n")
cat("RData files are saved in the 'RData' folder in your working directory.\n")

# Optional: Clean up temporary directories
cat("Cleaning up temporary directories...\n")
unlink(unzipped_path, recursive = TRUE)
unlink(all_flat_dir, recursive = TRUE)

cat("Process complete! Your LGR MGGA files have been successfully imported.\n")