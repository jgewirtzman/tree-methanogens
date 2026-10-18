# ==============================================================================
# 05_ddpcr_import_and_align.R -- import the ddPCR runs and align them to sample metadata
# ------------------------------------------------------------------------------
# Stage A; workdir: repository root.
#
# Sample metadata come from the four plate guides (inner and outer wood cores, mineral
# and organic soil); each well's sample is identified by extraction plate + position.
# Every QuantaSoft export in Raw Data/ (one per plate x target x threshold setting;
# "loose" or "strict" in the file name) is read, its first 23 columns kept, and joined
# to the metadata on plate + well. Standards, blank wells and wells with < 5000
# accepted droplets are removed. Concentrations are copies per uL of reaction as
# exported; conversion to copies per gram is done downstream (code/lib/ddpcr_constants.R).
#
# Inputs:  data/raw/ddpcr/Data Files/{inner,outer,mineral,organic}_guide.csv
#          data/raw/ddpcr/Data Files/Raw Data/*.csv
# Outputs: data/processed/molecular/processed_ddpcr_data.csv
# ==============================================================================

library(stringr)

base_path <- "data/raw/ddpcr/Data Files/"
raw_data_path <- file.path(base_path, "Raw Data")

# ---- sample metadata (plate guides) -------------------------------------------------
guide_files <- c("inner_guide.csv", "outer_guide.csv", "mineral_guide.csv", "organic_guide.csv")
metadata <- do.call(rbind, lapply(guide_files, function(f)
  read.csv(file.path(base_path, f), skip = 1, header = TRUE, stringsAsFactors = TRUE)))

# Inner.Core.Sample.ID holds species code, material_core type and lab ID
metadata$species <- substr(metadata$Inner.Core.Sample.ID, 1, 4)
core_info <- sub("\nL.*", "", metadata$Inner.Core.Sample.ID)
core_info <- sub(".*\n", "", core_info)
metadata$material <- sub("_.*", "", core_info)
metadata$core_type <- sub(".*_", "", core_info)
metadata$lab_id <- sub(".*: ", "", metadata$Inner.Core.Sample.ID)
metadata$plate_identifier <- paste(sub("Plate ", "", metadata$Extraction.Plate.ID),
                                   metadata$Plate.Position, sep = "")
cat("Metadata:", nrow(metadata), "samples\n")

# ---- ddPCR exports --------------------------------------------------------------------
data_files <- list.files(raw_data_path, pattern = "\\.csv$", full.names = TRUE)
data_files <- data_files[!grepl("Icon|\\._|DS_Store|\\.tmp", basename(data_files))]
ddpcr_data <- data.frame()
for (file_path in data_files) {
  data_file <- read.csv(file_path, stringsAsFactors = FALSE)
  file_name <- basename(file_path)
  data_file$analysis_type <- if (str_detect(file_name, "loose")) "loose" else
                             if (str_detect(file_name, "strict")) "strict" else "none"
  n_cols <- min(23, ncol(data_file) - 1)          # first 23 columns + analysis_type
  data_file <- data_file[, c(1:n_cols, ncol(data_file))]
  ddpcr_data <- rbind(ddpcr_data, data_file)
}
cat("ddPCR exports:", length(data_files), "files,", nrow(ddpcr_data), "wells\n")

# well "A01" -> "A1"; plate + well identifier as in the guides; lower-case targets
ddpcr_data$Well <- paste0(substr(ddpcr_data$Well, 1, 1), as.numeric(substr(ddpcr_data$Well, 2, 4)))
ddpcr_data$Sample.description.1 <- gsub("_", " ", ddpcr_data$Sample.description.1)
ddpcr_data$plate_identifier <- paste0(sub("Plate ", "", ddpcr_data$Sample.description.1), ddpcr_data$Well)
ddpcr_data$Target <- tolower(ddpcr_data$Target)

# ---- join and filter ----------------------------------------------------------------
final_data <- merge(metadata, ddpcr_data, by = "plate_identifier", all.x = FALSE, all.y = FALSE)
cat("Joined:", nrow(final_data), "rows x", ncol(final_data), "columns\n")
final_data <- final_data[!grepl("std", final_data$Sample.description.2, ignore.case = TRUE), ]
cat("After removing standards:", nrow(final_data), "rows\n")
final_data <- final_data[!grepl("yes", final_data$Blank.Well., ignore.case = TRUE), ]
cat("After removing blanks:", nrow(final_data), "rows\n")
low_quality <- final_data$Accepted.Droplets < 5000
cat("Wells with < 5000 accepted droplets:", sum(low_quality, na.rm = TRUE), "\n")
final_data <- final_data[!low_quality, ]
cat("Final:", nrow(final_data), "rows x", ncol(final_data), "columns | targets:",
    paste(unique(final_data$Target), collapse = ", "), "| species:", length(unique(final_data$species)), "\n")

write.csv(final_data, "data/processed/molecular/processed_ddpcr_data.csv", row.names = FALSE)
cat("Wrote data/processed/molecular/processed_ddpcr_data.csv\n")
