# Tree DBH Consensus List Generator with Comprehensive Variant Tracking
# Uses datasheets as authoritative source, fills missing values from other sources
# Handles inconsistent tree naming with manual corrections and comprehensive variant tracking

library(dplyr)
library(readr)
library(openxlsx)
library(stringr)
library(RecordLinkage)  # For fuzzy matching
library(tidyr)

# Functions for tree ID normalization and matching (improved to preserve all meaningful characters)
normalize_tree_id <- function(id) {
  # Handle NA and empty values
  result <- ifelse(is.na(id) | id == "", "", as.character(id))
  
  # Convert to lowercase
  result <- tolower(result)
  # Replace various apostrophe formats with standard apostrophe
  result <- str_replace_all(result, "'|'|'|prime", "'")
  # Normalize whitespace to single spaces
  result <- str_replace_all(result, "\\s+", " ")
  # Remove spaces around slashes and apostrophes for consistency
  result <- str_replace_all(result, "\\s*([/'`])\\s*", "\\1")
  # Remove ONLY non-alphanumeric punctuation, preserving letters, numbers, apostrophes, and slashes
  result <- str_replace_all(result, "[^a-z0-9/'\\s]", "")
  # Trim whitespace
  result <- str_trim(result)
  
  # Apply manual corrections based on domain knowledge
  result <- apply_manual_corrections(result)
  
  return(result)
}

# Function to apply manual corrections based on domain knowledge
apply_manual_corrections <- function(normalized_id) {
  # Create mapping of corrections
  corrections <- list(
    # ABT 14P variations
    "abt14p" = "abt 14p",
    
    # Bowser variations  
    "bowser dead red oak" = "bowser",
    
    # CUL8R variations
    "ybcul8r" = "cul8r",
    
    # H7/HM7 variations
    "hm7" = "h7",
    
    # Typo corrections
    "ob501" = "pb501",
    "wb403" = "wp403",
    
    # SF1/SA1 variations
    "sa1" = "sf1",
    
    # Berenice variations
    "wpberenice" = "wpbernice",
    
    # Confirmed formatting matches from user review
    "wpbernice" = "wp bernice",  # Standardize to spaced version
    "ab3/5" = "ab 3/5",         # Standardize to spaced version  
    "ab4/6" = "ab 4/6",         # Standardize to spaced version
    
    # Additional confirmed matches from latest review
    "abt14p" = "abt 14p",       # This was already included but ensuring consistency
    "ybcul8r" = "cul8r"         # This was already included but ensuring consistency
  )
  
  # Apply corrections
  for (i in seq_along(normalized_id)) {
    if (normalized_id[i] %in% names(corrections)) {
      normalized_id[i] <- corrections[[normalized_id[i]]]
    }
  }
  
  return(normalized_id)
}

# Function to find tree ID column with flexible naming
find_tree_id_column <- function(df) {
  # Common variations of tree ID column names
  possible_names <- c("Tree_ID", "tree_id", "Tree ID", "tree id", "TreeID", "treeid", 
                      "Tree.ID", "tree.id", "TREE_ID", "TREE.ID", "lab_id", "Lab_ID")
  
  # Find which one exists in the dataframe
  found_col <- intersect(possible_names, names(df))
  
  if (length(found_col) > 0) {
    return(found_col[1])  # Return the first match
  } else {
    # If no standard column found, return NULL and we'll handle it manually
    return(NULL)
  }
}

# Function to parse ddPCR tree IDs from Inner.Core.Sample.ID (vectorized)
parse_ddpcr_tree_id <- function(sample_id) {
  # Handle vectors properly
  result <- character(length(sample_id))
  
  for (i in seq_along(sample_id)) {
    if (is.na(sample_id[i]) || sample_id[i] == "") {
      result[i] <- ""
    } else {
      # Split on "Wood" and take the first part
      # Example: "ACRU_RM1Wood_OuterLab ID: 1a" becomes "ACRU_RM1"
      parts <- str_split(sample_id[i], "Wood")[[1]]
      if (length(parts) > 0) {
        tree_part <- parts[1]
        # Extract just the tree ID part (everything after the last underscore for species_treeID format)
        tree_id_match <- str_extract(tree_part, "[A-Za-z]+_([A-Za-z0-9/']+)$")
        if (!is.na(tree_id_match)) {
          # Extract just the tree ID (after the underscore)
          tree_id <- str_extract(tree_id_match, "(?<=_)[A-Za-z0-9/']+$")
          result[i] <- ifelse(is.na(tree_id), "", tree_id)
        } else {
          result[i] <- ""
        }
      } else {
        result[i] <- ""
      }
    }
  }
  
  return(result)
}

# Improved similarity function that respects ALL character differences
calculate_similarity_strict <- function(s1, s2) {
  if (is.na(s1) || is.na(s2) || s1 == "" || s2 == "") return(0)
  
  # Exact match gets highest score
  if (s1 == s2) return(1.0)
  
  # Extract the base letters (including x) and numbers separately
  extract_parts <- function(str) {
    # Extract all letter sequences (including x, which is meaningful)
    letters <- str_extract_all(str, "[a-z]+")[[1]]
    # Extract number parts (including those with slashes/apostrophes)
    numbers <- str_extract_all(str, "[0-9]+[/']*[0-9]*")[[1]]
    return(list(letters = paste(letters, collapse = ""), 
                numbers = paste(numbers, collapse = "")))
  }
  
  parts1 <- extract_parts(s1)
  parts2 <- extract_parts(s2)
  
  # If the number parts are different, these are different trees
  if (parts1$numbers != parts2$numbers && 
      parts1$numbers != "" && parts2$numbers != "") {
    return(0)  # Different numbered trees
  }
  
  # If letter parts are different (including x vs no x), these are different trees
  if (parts1$letters != parts2$letters && 
      parts1$letters != "" && parts2$letters != "") {
    # Only allow very minor differences (like spacing/punctuation)
    letter_sim <- jarowinkler(parts1$letters, parts2$letters)
    if (letter_sim < 0.95) return(0)  # Different letter combinations (including x differences)
  }
  
  # For cases where it's just formatting differences (spacing, punctuation), calculate similarity
  jw_sim <- jarowinkler(s1, s2)
  lev_sim <- 1 - (levenshteinDist(s1, s2) / max(nchar(s1), nchar(s2)))
  
  return(max(jw_sim, lev_sim))
}

# Read all data sources
cat("Reading data files...\n")

# 1. Read the authoritative datasheets file
datasheets <- read_csv("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/static-chambers/field/Gewirtzman_datasheets_fixing_dates.csv")

# 2. Read the Trees_DBH.csv file
trees_dbh <- read_csv("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/dbh/Trees_DBH.csv")

# 3. Read the Excel labels file
labels_excel <- read.xlsx("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/dbh/labels.xlsx", sheet = "Sheet1")

# 4. Read soil tree level means (no DBH, just for Tree IDs)
soil_data <- read_csv("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/static-chambers/field/soil_tree_level_means.csv")

# 5. Read flux data TEMP (has DBH data)
flux_temp <- read_csv("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/static-chambers/flux/flux_data_TEMP.csv")

# 6. Read ddPCR data (no DBH, just for Tree IDs) - with encoding fix
cat("Reading ddPCR data with encoding handling...\n")
tryCatch({
  ddpcr_data <- read_csv("/Users/jongewirtzman/Google Drive/Research/YMF Tree Microbiomes & Methane/Tree Methane Lab/Data/ddPCR Data/Data Files/Collated Data/ddPCR_tree_full_data.csv", 
                         locale = locale(encoding = "latin1"))
}, error = function(e) {
  cat("Error with latin1 encoding, trying UTF-8...\n")
  tryCatch({
    ddpcr_data <<- read_csv("/Users/jongewirtzman/Google Drive/Research/YMF Tree Microbiomes & Methane/Tree Methane Lab/Data/ddPCR Data/Data Files/Collated Data/ddPCR_tree_full_data.csv", 
                            locale = locale(encoding = "UTF-8"))
  }, error = function(e2) {
    cat("Error with UTF-8, trying cp1252...\n")
    ddpcr_data <<- read_csv("/Users/jongewirtzman/Google Drive/Research/YMF Tree Microbiomes & Methane/Tree Methane Lab/Data/ddPCR Data/Data Files/Collated Data/ddPCR_tree_full_data.csv", 
                            locale = locale(encoding = "cp1252"))
  })
})

# 7. Read flux data with dimensions (has DBH data)
flux_dimensions <- read_csv("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/static-chambers/flux_code/flux_data_with_dimensions_weather_and_final_airtemp.csv")

# Clean and standardize data with enhanced ID matching
cat("Cleaning and standardizing data...\n")

# Function to prepare ALL data from files with DBH columns (keeps records even without DBH)
prepare_data_with_dbh_comprehensive <- function(df, source_name, tree_id_col = NULL) {
  # Auto-detect tree ID column if not specified
  if (is.null(tree_id_col)) {
    tree_id_col <- find_tree_id_column(df)
    if (is.null(tree_id_col)) {
      stop(paste("Could not find tree ID column in", source_name, "file"))
    }
  }
  
  # Auto-detect DBH column (could be "DBH", "dbh", "Dbh")
  dbh_cols <- c("DBH", "dbh", "Dbh", "DBH_cm", "dbh_cm")
  dbh_col <- intersect(dbh_cols, names(df))
  if (length(dbh_col) == 0) {
    stop(paste("Could not find DBH column in", source_name, "file"))
  }
  dbh_col <- dbh_col[1]
  
  # Auto-detect Plot column
  plot_cols <- c("Plot", "plot", "PLOT", "Site", "site", "Location", "location")
  plot_col <- intersect(plot_cols, names(df))
  if (length(plot_col) == 0) {
    df$Plot <- "Unknown"
    plot_col <- "Plot"
  } else {
    plot_col <- plot_col[1]
  }
  
  cat("  ", source_name, "- using columns:", tree_id_col, "(Tree ID),", dbh_col, "(DBH),", plot_col, "(Plot)\n")
  
  df %>%
    select(all_of(c(plot_col, tree_id_col, dbh_col))) %>%
    rename(Plot = all_of(plot_col), Tree_ID_original = all_of(tree_id_col), DBH = all_of(dbh_col)) %>%
    filter(!is.na(Tree_ID_original), Tree_ID_original != "") %>%  # Only filter out missing Tree IDs
    mutate(
      DBH = as.numeric(DBH),  # Convert DBH but don't filter out NAs yet
      Tree_ID_normalized = normalize_tree_id(Tree_ID_original),
      source = source_name
    ) %>%
    filter(Tree_ID_normalized != "") %>%  # Only filter out trees that couldn't be normalized
    distinct()
}

# Function for DBH-containing files that separates trees with/without DBH
prepare_data_with_optional_dbh <- function(df, source_name, tree_id_col = NULL) {
  all_data <- prepare_data_with_dbh_comprehensive(df, source_name, tree_id_col)
  
  # Split into DBH and non-DBH records
  with_dbh <- all_data %>% filter(!is.na(DBH))
  without_dbh <- all_data %>% filter(is.na(DBH)) %>% select(-DBH)
  
  return(list(with_dbh = with_dbh, without_dbh = without_dbh))
}

# Special function for ddPCR data parsing
prepare_ddpcr_data <- function(df, source_name) {
  cat("  ", source_name, "- parsing Tree IDs from Inner.Core.Sample.ID\n")
  
  df %>%
    filter(!is.na(Inner.Core.Sample.ID), Inner.Core.Sample.ID != "") %>%
    mutate(
      Tree_ID_original = parse_ddpcr_tree_id(Inner.Core.Sample.ID),
      Plot = "Unknown",  # ddPCR data doesn't have plot info
      Tree_ID_normalized = normalize_tree_id(Tree_ID_original),
      source = source_name
    ) %>%
    filter(Tree_ID_normalized != "") %>%
    select(Plot, Tree_ID_original, Tree_ID_normalized, source) %>%
    distinct()
}

# Prepare data for files WITHOUT DBH (just for collecting Tree IDs)
prepare_data_no_dbh_auto <- function(df, source_name, tree_id_col = NULL) {
  # Auto-detect tree ID column if not specified
  if (is.null(tree_id_col)) {
    tree_id_col <- find_tree_id_column(df)
    if (is.null(tree_id_col)) {
      stop(paste("Could not find tree ID column in", source_name, "file"))
    }
  }
  
  # Auto-detect Plot column
  plot_cols <- c("Plot", "plot", "PLOT", "Site", "site", "Location", "location")
  plot_col <- intersect(plot_cols, names(df))
  if (length(plot_col) == 0) {
    df$Plot <- "Unknown"
    plot_col <- "Plot"
  } else {
    plot_col <- plot_col[1]
  }
  
  cat("  ", source_name, "- using columns:", tree_id_col, "(Tree ID),", plot_col, "(Plot)\n")
  
  df %>%
    select(all_of(c(plot_col, tree_id_col))) %>%
    rename(Plot = all_of(plot_col), Tree_ID_original = all_of(tree_id_col)) %>%
    filter(!is.na(Tree_ID_original), Tree_ID_original != "") %>%
    mutate(
      Tree_ID_normalized = normalize_tree_id(Tree_ID_original),
      source = source_name
    ) %>%
    filter(Tree_ID_normalized != "") %>%
    distinct()
}

# Process files with comprehensive approach
cat("Processing files with DBH data (comprehensive approach):\n")
datasheets_data <- prepare_data_with_optional_dbh(datasheets, "datasheets")
trees_dbh_data <- prepare_data_with_optional_dbh(trees_dbh, "trees_dbh")
labels_data <- prepare_data_with_optional_dbh(labels_excel, "labels")
flux_temp_data <- prepare_data_with_optional_dbh(flux_temp, "flux_temp")
flux_dimensions_data <- prepare_data_with_optional_dbh(flux_dimensions, "flux_dimensions")

cat("\nProcessing files without DBH data:\n")
soil_prep <- prepare_data_no_dbh_auto(soil_data, "soil_data")
ddpcr_prep <- prepare_ddpcr_data(ddpcr_data, "ddpcr_data")

# For files with multiple measurements per tree, take the most recent DBH
datasheets_clean <- datasheets_data$with_dbh %>%
  group_by(Tree_ID_normalized) %>%
  summarise(
    Tree_ID_original = last(Tree_ID_original),
    Plot = first(Plot),
    DBH = last(DBH),
    source = "datasheets",
    .groups = "drop"
  )

flux_temp_clean <- flux_temp_data$with_dbh %>%
  group_by(Tree_ID_normalized) %>%
  summarise(
    Tree_ID_original = last(Tree_ID_original),
    Plot = first(Plot),
    DBH = last(DBH),
    source = "flux_temp",
    .groups = "drop"
  )

flux_dimensions_clean <- flux_dimensions_data$with_dbh %>%
  group_by(Tree_ID_normalized) %>%
  summarise(
    Tree_ID_original = last(Tree_ID_original),
    Plot = first(Plot),
    DBH = last(DBH),
    source = "flux_dimensions",
    .groups = "drop"
  )

# Clean single-measurement files
trees_dbh_clean <- trees_dbh_data$with_dbh
labels_clean <- labels_data$with_dbh

cat("\nData summary after comprehensive cleaning:\n")
cat("Datasheets (authoritative):", nrow(datasheets_clean), "trees with DBH\n")
cat("Trees_DBH file:", nrow(trees_dbh_clean), "trees with DBH\n")
cat("Labels Excel:", nrow(labels_clean), "trees with DBH\n")
cat("Flux TEMP file:", nrow(flux_temp_clean), "trees with DBH\n")
cat("Flux dimensions file:", nrow(flux_dimensions_clean), "trees with DBH\n")
cat("Soil data file:", nrow(soil_prep), "unique tree IDs (no DBH)\n")
cat("ddPCR data file:", nrow(ddpcr_prep), "unique tree IDs (no DBH)\n")

# Additional trees without DBH from DBH-containing files
additional_trees <- bind_rows(
  datasheets_data$without_dbh,
  trees_dbh_data$without_dbh,
  labels_data$without_dbh,
  flux_temp_data$without_dbh,
  flux_dimensions_data$without_dbh
)

if (nrow(additional_trees) > 0) {
  cat("Additional trees found in DBH files but without DBH values:", nrow(additional_trees), "\n")
}

# Show sample of parsed ddPCR tree IDs
if (nrow(ddpcr_prep) > 0) {
  cat("\nSample ddPCR parsed tree IDs:\n")
  print(head(ddpcr_prep %>% select(Tree_ID_original), 10))
}

# ===========================================
# COMPREHENSIVE DIAGNOSTIC SECTION
# ===========================================

cat("\n=== COMPREHENSIVE TREE ID DIAGNOSTIC ===\n")

# Function to extract all unique tree IDs from each raw file
extract_all_tree_ids <- function(df, file_name, tree_id_col) {
  cat("\n", file_name, ":\n")
  
  if (tree_id_col %in% names(df)) {
    all_ids <- df[[tree_id_col]]
    # Remove NAs and empty strings
    clean_ids <- all_ids[!is.na(all_ids) & all_ids != ""]
    unique_ids <- unique(clean_ids)
    
    cat("  Total records:", nrow(df), "\n")
    cat("  Records with Tree ID:", length(clean_ids), "\n")
    cat("  Unique Tree IDs:", length(unique_ids), "\n")
    
    # Show first 15 unique IDs
    cat("  Sample IDs:", paste(head(unique_ids, 15), collapse = ", "), "\n")
    if (length(unique_ids) > 15) {
      cat("  ... and", length(unique_ids) - 15, "more\n")
    }
    
    return(data.frame(
      file = file_name,
      total_records = nrow(df),
      records_with_id = length(clean_ids),
      unique_ids = length(unique_ids),
      sample_ids = paste(head(unique_ids, 5), collapse = ", "),
      stringsAsFactors = FALSE
    ))
  } else {
    cat("  ERROR: Column", tree_id_col, "not found!\n")
    cat("  Available columns:", paste(names(df), collapse = ", "), "\n")
    return(NULL)
  }
}

# Extract from all raw files
cat("Raw file analysis:\n")
diagnostic_summary <- bind_rows(
  extract_all_tree_ids(datasheets, "datasheets", "Tree_ID"),
  extract_all_tree_ids(trees_dbh, "trees_dbh", "Tree ID"),
  extract_all_tree_ids(labels_excel, "labels_excel", "Tree.ID"),
  extract_all_tree_ids(flux_temp, "flux_temp", "Tree_ID"),
  extract_all_tree_ids(flux_dimensions, "flux_dimensions", "tree_id"),
  extract_all_tree_ids(soil_data, "soil_data", "Tree_ID")
)

# Special handling for ddPCR
cat("\nddPCR_data:\n")
if ("Inner.Core.Sample.ID" %in% names(ddpcr_data)) {
  sample_ids <- ddpcr_data$Inner.Core.Sample.ID
  clean_sample_ids <- sample_ids[!is.na(sample_ids) & sample_ids != ""]
  parsed_ids <- parse_ddpcr_tree_id(clean_sample_ids)
  clean_parsed <- parsed_ids[parsed_ids != ""]
  
  cat("  Total records:", nrow(ddpcr_data), "\n")
  cat("  Records with Sample ID:", length(clean_sample_ids), "\n")
  cat("  Successfully parsed Tree IDs:", length(clean_parsed), "\n")
  cat("  Unique parsed Tree IDs:", length(unique(clean_parsed)), "\n")
  cat("  Sample parsed IDs:", paste(head(unique(clean_parsed), 10), collapse = ", "), "\n")
  
  diagnostic_summary <- bind_rows(diagnostic_summary, data.frame(
    file = "ddpcr_data",
    total_records = nrow(ddpcr_data),
    records_with_id = length(clean_sample_ids),
    unique_ids = length(unique(clean_parsed)),
    sample_ids = paste(head(unique(clean_parsed), 5), collapse = ", "),
    stringsAsFactors = FALSE
  ))
}

cat("\n=== RAW FILE SUMMARY ===\n")
print(diagnostic_summary)

# Now check what's happening in our processed data
cat("\n=== PROCESSED DATA ANALYSIS ===\n")

processed_summary <- data.frame(
  source = c("datasheets", "trees_dbh", "labels", "flux_temp", "flux_dimensions", "soil_data", "ddpcr_data"),
  with_dbh = c(
    nrow(datasheets_clean),
    nrow(trees_dbh_clean), 
    nrow(labels_clean),
    nrow(flux_temp_clean),
    nrow(flux_dimensions_clean),
    0, # soil has no DBH
    0  # ddPCR has no DBH
  ),
  without_dbh = c(
    if(exists("datasheets_data")) nrow(datasheets_data$without_dbh) else 0,
    if(exists("trees_dbh_data")) nrow(trees_dbh_data$without_dbh) else 0,
    if(exists("labels_data")) nrow(labels_data$without_dbh) else 0,
    if(exists("flux_temp_data")) nrow(flux_temp_data$without_dbh) else 0,
    if(exists("flux_dimensions_data")) nrow(flux_dimensions_data$without_dbh) else 0,
    nrow(soil_prep),
    nrow(ddpcr_prep)
  ),
  stringsAsFactors = FALSE
)

processed_summary$total_processed <- processed_summary$with_dbh + processed_summary$without_dbh

cat("Processed data summary:\n")
print(processed_summary)

cat("\n=== POTENTIAL ISSUES ===\n")
for (i in 1:nrow(diagnostic_summary)) {
  raw_count <- diagnostic_summary$unique_ids[i]
  file_name <- diagnostic_summary$file[i]
  
  if (file_name %in% processed_summary$source) {
    processed_count <- processed_summary$total_processed[processed_summary$source == file_name]
    
    if (raw_count != processed_count) {
      cat("⚠️ ", file_name, ": Raw =", raw_count, "unique IDs, Processed =", processed_count, "trees\n")
      cat("   Difference:", raw_count - processed_count, "trees lost\n")
    } else {
      cat("✅ ", file_name, ": All", raw_count, "trees captured\n")
    }
  }
}

# Create master list of all unique tree IDs (normalized) including ALL sources
all_tree_ids <- bind_rows(
  datasheets_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot, source),
  trees_dbh_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot, source),
  labels_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot, source),
  flux_temp_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot, source),
  flux_dimensions_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot, source),
  soil_prep %>% select(Tree_ID_normalized, Tree_ID_original, Plot, source),
  ddpcr_prep %>% select(Tree_ID_normalized, Tree_ID_original, Plot, source),
  # Include additional trees found in DBH files but without DBH values
  additional_trees %>% select(Tree_ID_normalized, Tree_ID_original, Plot, source)
) %>%
  distinct()

cat("\nTotal unique normalized tree IDs:", length(unique(all_tree_ids$Tree_ID_normalized)), "\n")
cat("Total tree ID records (including all variants):", nrow(all_tree_ids), "\n")

# Show breakdown by source
cat("\nTree ID records by source:\n")
source_breakdown <- all_tree_ids %>%
  group_by(source) %>%
  summarise(
    unique_trees = n_distinct(Tree_ID_normalized),
    total_records = n(),
    .groups = "drop"
  )
print(source_breakdown)

# Find potential matches using improved fuzzy matching
cat("Performing enhanced fuzzy matching for similar tree IDs...\n")

# Get all unique normalized IDs
unique_ids <- unique(all_tree_ids$Tree_ID_normalized)
similarity_threshold <- 0.90  # Higher threshold since we're being more strict

# Create a matching matrix with improved logic
matches <- data.frame()
for (i in 1:length(unique_ids)) {
  for (j in (i+1):length(unique_ids)) {
    if (j > length(unique_ids)) break
    
    id1 <- unique_ids[i]
    id2 <- unique_ids[j]
    similarity <- calculate_similarity_strict(id1, id2)
    
    if (similarity >= similarity_threshold) {
      matches <- bind_rows(matches, data.frame(
        id1 = id1,
        id2 = id2,
        similarity = similarity
      ))
    }
  }
}

if (nrow(matches) > 0) {
  cat("Found", nrow(matches), "potential formatting matches (same tree, different formatting):\n")
  print(matches %>% arrange(desc(similarity)))
} else {
  cat("No formatting matches found above threshold", similarity_threshold, "\n")
}

# Create a comprehensive tree ID mapping with separate columns for each variant
cat("\nCreating comprehensive tree ID mapping with separate variant columns...\n")

# First, create the basic mapping
tree_id_mapping_basic <- all_tree_ids %>%
  group_by(Tree_ID_normalized) %>%
  summarise(
    primary_id = first(Tree_ID_original),
    all_variants = list(unique(Tree_ID_original)),
    sources = paste(unique(source), collapse = ", "),
    plots = paste(unique(Plot), collapse = ", "),
    name_variants = n(),
    .groups = "drop"
  )

# Create separate columns for each variant name
create_variant_columns <- function(mapping_df) {
  # Find the maximum number of variants any tree has
  max_variants <- max(sapply(mapping_df$all_variants, length))
  
  # Create base dataframe
  result_df <- data.frame(
    Tree_ID_normalized = mapping_df$Tree_ID_normalized,
    primary_id = mapping_df$primary_id,
    sources = mapping_df$sources,
    plots = mapping_df$plots,
    name_variants = mapping_df$name_variants,
    stringsAsFactors = FALSE
  )
  
  # Add columns for each variant
  for (i in 1:max_variants) {
    col_name <- paste0("variant_", i)
    result_df[[col_name]] <- sapply(mapping_df$all_variants, function(x) {
      if (length(x) >= i) x[i] else NA
    })
  }
  
  # Create source-specific variant columns to track which sources contributed which names
  source_list <- c("datasheets", "trees_dbh", "labels", "flux_temp", "flux_dimensions", "soil_data", "ddpcr_data")
  
  for (source in source_list) {
    col_name <- paste0("name_in_", source)
    result_df[[col_name]] <- sapply(seq_len(nrow(mapping_df)), function(i) {
      # Get variants from this specific source
      source_variants <- all_tree_ids %>%
        filter(Tree_ID_normalized == mapping_df$Tree_ID_normalized[i], 
               source == !!source) %>%
        pull(Tree_ID_original) %>%
        unique()
      
      if (length(source_variants) > 0) {
        paste(source_variants, collapse = " | ")
      } else {
        NA
      }
    })
  }
  
  return(result_df)
}

tree_id_mapping <- create_variant_columns(tree_id_mapping_basic)

# Show summary
cat("Maximum variants for any single tree:", sum(grepl("^variant_", names(tree_id_mapping))), "\n")
cat("Trees with multiple name variants:", sum(tree_id_mapping$name_variants > 1), "\n")

# Show example of variant tracking
cat("\nExample of variant tracking (trees with most variants):\n")
example_variants <- tree_id_mapping %>%
  filter(name_variants > 3) %>%
  arrange(desc(name_variants)) %>%
  head(3) %>%
  select(primary_id, starts_with("variant_"), name_variants)
print(example_variants)

# Create consensus list with proper matching and comprehensive variant tracking
cat("\nCreating consensus list with enhanced matching and complete variant tracking...\n")

# Step 1: Start with datasheets as authoritative source
consensus_base <- datasheets_clean %>%
  left_join(tree_id_mapping, by = c("Tree_ID_normalized" = "Tree_ID_normalized")) %>%
  select(Tree_ID_normalized, Tree_ID_original, Plot, DBH, source, 
         starts_with("variant_"), starts_with("name_in_"), name_variants)

# Step 2: Find trees missing from datasheets using normalized IDs
all_normalized_trees <- bind_rows(
  trees_dbh_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot),
  labels_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot),
  flux_temp_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot),
  flux_dimensions_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot),
  datasheets_clean %>% select(Tree_ID_normalized, Tree_ID_original, Plot),
  soil_prep %>% select(Tree_ID_normalized, Tree_ID_original, Plot),
  ddpcr_prep %>% select(Tree_ID_normalized, Tree_ID_original, Plot)
) %>%
  distinct(Tree_ID_normalized, .keep_all = TRUE)

missing_from_datasheets <- all_normalized_trees %>%
  filter(!Tree_ID_normalized %in% consensus_base$Tree_ID_normalized)

cat("Trees missing from datasheets (after normalization):", nrow(missing_from_datasheets), "\n")

# Step 3: Fill missing trees from Trees_DBH file first (priority 1)
missing_filled_from_trees <- missing_from_datasheets %>%
  inner_join(trees_dbh_clean, by = "Tree_ID_normalized") %>%
  left_join(tree_id_mapping, by = c("Tree_ID_normalized" = "Tree_ID_normalized")) %>%
  select(Tree_ID_normalized, Tree_ID_original = Tree_ID_original.y, 
         Plot = Plot.y, DBH, source, starts_with("variant_"), starts_with("name_in_"), name_variants)

# Step 4: Fill remaining missing trees from Labels Excel (priority 2)
still_missing_1 <- missing_from_datasheets %>%
  filter(!Tree_ID_normalized %in% missing_filled_from_trees$Tree_ID_normalized)

missing_filled_from_labels <- still_missing_1 %>%
  inner_join(labels_clean, by = "Tree_ID_normalized") %>%
  left_join(tree_id_mapping, by = c("Tree_ID_normalized" = "Tree_ID_normalized")) %>%
  select(Tree_ID_normalized, Tree_ID_original = Tree_ID_original.y, 
         Plot = Plot.y, DBH, source, starts_with("variant_"), starts_with("name_in_"), name_variants)

# Step 5: Fill from flux_temp file (priority 3)
still_missing_2 <- still_missing_1 %>%
  filter(!Tree_ID_normalized %in% missing_filled_from_labels$Tree_ID_normalized)

missing_filled_from_flux_temp <- still_missing_2 %>%
  inner_join(flux_temp_clean, by = "Tree_ID_normalized") %>%
  left_join(tree_id_mapping, by = c("Tree_ID_normalized" = "Tree_ID_normalized")) %>%
  select(Tree_ID_normalized, Tree_ID_original = Tree_ID_original.y, 
         Plot = Plot.y, DBH, source, starts_with("variant_"), starts_with("name_in_"), name_variants)

# Step 6: Fill from flux_dimensions file (priority 4)
still_missing_3 <- still_missing_2 %>%
  filter(!Tree_ID_normalized %in% missing_filled_from_flux_temp$Tree_ID_normalized)

missing_filled_from_flux_dims <- still_missing_3 %>%
  inner_join(flux_dimensions_clean, by = "Tree_ID_normalized") %>%
  left_join(tree_id_mapping, by = c("Tree_ID_normalized" = "Tree_ID_normalized")) %>%
  select(Tree_ID_normalized, Tree_ID_original = Tree_ID_original.y, 
         Plot = Plot.y, DBH, source, starts_with("variant_"), starts_with("name_in_"), name_variants)

# Step 7: Combine all sources into final consensus list
final_consensus <- bind_rows(
  consensus_base,
  missing_filled_from_trees,
  missing_filled_from_labels,
  missing_filled_from_flux_temp,
  missing_filled_from_flux_dims
) %>%
  arrange(Plot, Tree_ID_original) %>%
  mutate(
    consensus_tree_id = Tree_ID_original,
    normalized_id = Tree_ID_normalized
  ) %>%
  select(consensus_tree_id, normalized_id, Plot, DBH, source, 
         starts_with("variant_"), starts_with("name_in_"), name_variants)

# Step 8: Create a separate table for trees with multiple name variants
trees_with_variants <- final_consensus %>%
  filter(name_variants > 1) %>%
  arrange(desc(name_variants))

# Step 9: Identify trees still missing DBH data
trees_still_missing <- all_normalized_trees %>%
  filter(!Tree_ID_normalized %in% final_consensus$normalized_id) %>%
  left_join(tree_id_mapping, by = c("Tree_ID_normalized" = "Tree_ID_normalized")) %>%
  mutate(
    consensus_tree_id = Tree_ID_original,
    DBH = NA,
    source = "missing"
  ) %>%
  select(consensus_tree_id, normalized_id = Tree_ID_normalized, Plot, DBH, source, 
         starts_with("variant_"), starts_with("name_in_"), name_variants)

# Final summary with enhanced reporting
cat("\n=== ENHANCED CONSENSUS LIST SUMMARY ===\n")
cat("Total trees with DBH data:", nrow(final_consensus), "\n")
cat("From datasheets (authoritative):", sum(final_consensus$source == "datasheets"), "\n")
cat("From Trees_DBH.csv:", sum(final_consensus$source == "trees_dbh"), "\n")
cat("From Labels Excel:", sum(final_consensus$source == "labels"), "\n")
cat("From Flux TEMP:", sum(final_consensus$source == "flux_temp"), "\n")
cat("From Flux dimensions:", sum(final_consensus$source == "flux_dimensions"), "\n")
cat("Trees still missing DBH:", nrow(trees_still_missing), "\n")
cat("Trees with multiple name variants:", nrow(trees_with_variants), "\n")

# Display source breakdown by plot
cat("\nSource breakdown by plot:\n")
print(final_consensus %>%
        group_by(Plot, source) %>%
        summarise(count = n(), .groups = "drop") %>%
        pivot_wider(names_from = source, values_from = count, values_fill = 0))

# Show trees with the most name variants
if (nrow(trees_with_variants) > 0) {
  cat("\nTrees with multiple name variants (top 10):\n")
  print(trees_with_variants %>%
          select(consensus_tree_id, variant_1, variant_2, name_variants, source) %>%
          head(10))
}

# Check for potential data quality issues
cat("\n=== DATA QUALITY CHECKS ===\n")

# Check for very similar but unmatched IDs (potential typos)
unmatched_similar <- data.frame()
all_ids <- unique(c(final_consensus$normalized_id, trees_still_missing$normalized_id))
for (i in 1:(length(all_ids)-1)) {
  for (j in (i+1):length(all_ids)) {
    sim <- calculate_similarity_strict(all_ids[i], all_ids[j])
    # Look for moderate similarity that might indicate typos
    # but exclude cases where it's clearly just different numbers
    if (sim > 0.7 && sim < similarity_threshold) {
      unmatched_similar <- bind_rows(unmatched_similar, 
                                     data.frame(id1 = all_ids[i], id2 = all_ids[j], similarity = sim))
    }
  }
}

if (nrow(unmatched_similar) > 0) {
  cat("Potentially similar IDs that might be typos (review manually):\n")
  print(unmatched_similar %>% arrange(desc(similarity)) %>% head(15))
} else {
  cat("No potentially similar unmatched IDs found (good - suggests clean numbering).\n")
}

# Additional check: Look for trees that might have formatting issues
cat("\n=== FORMATTING ANALYSIS ===\n")
formatting_variants <- tree_id_mapping %>%
  filter(name_variants > 1) %>%
  arrange(desc(name_variants)) %>%
  head(20)

if (nrow(formatting_variants) > 0) {
  cat("Trees with the most name variants (likely formatting differences):\n")
  print(formatting_variants %>% select(primary_id, variant_1, variant_2, name_variants, sources))
}

# Save the enhanced results
cat("\nSaving comprehensive results...\n")

# Main consensus list with all variant information
write_csv(final_consensus, "tree_dbh_consensus_comprehensive.csv")

# Complete list including missing trees
complete_list <- bind_rows(final_consensus, trees_still_missing) %>%
  arrange(Plot, consensus_tree_id)
write_csv(complete_list, "tree_dbh_complete_comprehensive.csv")

# Tree ID mapping table for reference with all variants
write_csv(tree_id_mapping, "tree_id_comprehensive_mapping.csv")

# Trees with multiple variants (for manual review)
if (nrow(trees_with_variants) > 0) {
  write_csv(trees_with_variants, "trees_multiple_variants_comprehensive.csv")
}

# Fuzzy matches found (for manual review)
if (nrow(matches) > 0) {
  write_csv(matches, "potential_formatting_matches.csv")
}

# Moderately similar unmatched IDs (for manual review)
if (nrow(unmatched_similar) > 0) {
  write_csv(unmatched_similar, "potentially_similar_unmatched_ids.csv")
}

# Summary statistics
summary_stats <- final_consensus %>%
  filter(!is.na(DBH)) %>%
  group_by(Plot) %>%
  summarise(
    tree_count = n(),
    mean_dbh = round(mean(DBH, na.rm = TRUE), 2),
    median_dbh = round(median(DBH, na.rm = TRUE), 2),
    min_dbh = min(DBH, na.rm = TRUE),
    max_dbh = max(DBH, na.rm = TRUE),
    trees_with_variants = sum(name_variants > 1, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_stats, "dbh_summary_comprehensive.csv")

cat("\nComprehensive files created:\n")
cat("- tree_dbh_consensus_comprehensive.csv: Final consensus with complete variant tracking\n")
cat("- tree_dbh_complete_comprehensive.csv: Complete list including missing trees\n")  
cat("- tree_id_comprehensive_mapping.csv: Master mapping of all variants by source\n")
cat("- trees_multiple_variants_comprehensive.csv: Trees with multiple variants\n")
cat("- potential_formatting_matches.csv: Formatting matches for review\n")
cat("- potentially_similar_unmatched_ids.csv: Similar IDs for review\n")
cat("- dbh_summary_comprehensive.csv: Enhanced summary statistics\n")

cat("\n=== SAMPLE OF FINAL CONSENSUS LIST WITH ALL VARIANTS ===\n")
print(head(final_consensus %>% select(consensus_tree_id, normalized_id, Plot, DBH, source, 
                                      starts_with("variant_"), name_variants), 5))

cat("\n=== SAMPLE OF SOURCE-SPECIFIC NAME TRACKING ===\n")
print(head(final_consensus %>% select(consensus_tree_id, starts_with("name_in_")), 3))

if(nrow(trees_still_missing) > 0) {
  cat("\n=== TREES STILL MISSING DBH ===\n")
  print(head(trees_still_missing %>% select(consensus_tree_id, normalized_id, Plot, starts_with("variant_")), 10))
  cat("... and", max(0, nrow(trees_still_missing) - 10), "more\n")
}

cat("\n=== COMPREHENSIVE VARIANT TRACKING COMPLETED ===\n")
cat("The final table now includes:\n")
cat("- variant_1, variant_2, etc: All name variations found across all sources\n")
cat("- name_in_datasheets, name_in_trees_dbh, etc: Specific names used in each source\n")
cat("- This comprehensive mapping will help you match data across files in future analysis\n")
cat("- Manual corrections applied: ABT 14P, bowser, cul8r, h7, typo fixes, etc.\n")
cat("- Perfect for future data mapping and quality control\n")

cat("\nAnalysis completed successfully!\n")