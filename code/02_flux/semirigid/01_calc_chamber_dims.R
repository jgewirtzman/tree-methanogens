# ==============================================================================
# 01_calc_chamber_dims.R -- semi-rigid stem chambers: area and volume from the chamber sheets
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/02_flux/semirigid/).
#
# For every row of the monthly chamber sheets (Tree_fluxes_<month>.xlsx, one sheet per
# chamber), the chamber's surface area on the stem and its headspace volume, by
# Siegenthaler et al. (2016): an annulus sector of the stem (diameter D, sleeve
# thickness T, chamber length L and height H; sector fraction K = L / (pi (D + 2T)),
# capped at 1 where the chamber wraps the whole stem) minus the measured wedges. The
# sheet's own Sc / Vc are kept beside the corrected values for comparison.
#
# The deployment geometry used for the fluxes is computed from the same sheets by
# code/lib/stem_chamber_geometry.R (02_join_flux_geometry.R); these workbooks are the
# per-row record of the correction.
#
# Inputs:  data/raw/flux_chamber_corrections/Chamber vol corr/Tree_fluxes_{june,july,
#          august,september,october,december}.xlsx
# Outputs: data/processed/flux/chamber_geometry/Tree_fluxes_<month>_corrected.xlsx
#          data/processed/flux/chamber_geometry/combined_corrected_results.xlsx
# ==============================================================================

source("../../lib/project_library.R")   # the project package library (renv.lock), however this script is started
library(readxl)
library(writexl)
library(dplyr)
library(stringr)
library(readr)
library(purrr)

input_dir  <- "../../../data/raw/flux_chamber_corrections/Chamber vol corr"
output_dir <- "../../../data/processed/flux/chamber_geometry"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

files <- c(
  "Tree_fluxes_june.xlsx",
  "Tree_fluxes_july.xlsx", 
  "Tree_fluxes_august.xlsx",
  "Tree_fluxes_september.xlsx",
  "Tree_fluxes_october.xlsx",
  "Tree_fluxes_december.xlsx"
)

num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))

# first column whose name matches one of the patterns, tried in order (the sheets'
# headers differ slightly between months)
find_col <- function(df, patterns) {
  for (pattern in patterns) {
    matches <- str_detect(names(df), regex(pattern, ignore_case = TRUE))
    if (any(matches)) {
      return(names(df)[which(matches)[1]])
    }
  }
  return(NULL)
}

# Area and volume per sheet row (see header)
calculate_chamber_geometry <- function(df) {
  
  # Column name patterns
  col_patterns <- list(
    Dstem   = c("^D\\s*stem", "Dstem", "D_stem", "stem_diameter"),
    T       = c("^T\\s*\\(sleeve\\s*thickness", "^T\\s*\\(", "thickness", "sleeve_thickness"),
    L       = c("^Chamber\\s*length", "^L\\b", "length", "chamber_length"),
    H       = c("^Chamber\\s*height", "^H\\b", "height", "chamber_height"),
    Swedges = c("^Swedges", "S_wedges", "wedge.*surface"),
    Vwedges = c("Two\\s*V\\s*wedges", "V_wedges", "wedge.*volume"),
    Sc_orig = c("^Sc\\s*\\(cm2\\)", "^Sc\\b", "S_chamber", "surface_area"),
    Vc_orig = c("^Vc\\b", "V_chamber", "chamber_volume")
  )
  
  # Find columns
  cols <- map(col_patterns, ~find_col(df, .x))
  
  required <- c("Dstem", "T", "L", "H", "Swedges", "Vwedges", "Sc_orig", "Vc_orig")
  missing <- required[sapply(required, function(x) is.null(cols[[x]]))]
  
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(names(df), collapse = ", "))
  }
  
  # Extract measured values
  Dstem_cm <- num(df[[cols$Dstem]])
  T_cm     <- num(df[[cols$T]])
  L_cm     <- num(df[[cols$L]])
  H_cm     <- num(df[[cols$H]])
  
  # Validate inputs
  valid_inputs <- is.finite(Dstem_cm) & is.finite(T_cm) & is.finite(L_cm) & is.finite(H_cm) &
    (Dstem_cm > 0) & (T_cm > 0) & (L_cm > 0) & (H_cm > 0)
  
  # measured wedge area and volume from the sheet (a blank cell counts as no wedge)
  Swedges_cm2 <- num(df[[cols$Swedges]])
  Vwedges_cm3 <- num(df[[cols$Vwedges]])
  wedge_source <- "measured"
  Swedges_cm2[!is.finite(Swedges_cm2)] <- 0
  Vwedges_cm3[!is.finite(Vwedges_cm3)] <- 0
  
  # Calculate basic geometry
  Dext_cm <- Dstem_cm + 2 * T_cm
  circumference_cm <- pi * Dext_cm
  
  # Calculate sector fraction with wraparound correction
  K_raw <- L_cm / circumference_cm
  K_eff <- ifelse(is.finite(K_raw), pmin(K_raw, 1.0), NA_real_)
  wraparound_flag <- is.finite(K_raw) & (K_raw > 1.0)
  
  # Calculate volumes - annulus sector minus wedges
  annulus_area_cm2 <- (pi/4) * (Dext_cm^2 - Dstem_cm^2)
  Vc_cm3_raw <- K_eff * annulus_area_cm2 * H_cm - Vwedges_cm3
  
  # Calculate surface areas - stem surface sector minus wedges  
  stem_surface_full_cm2 <- pi * Dstem_cm * H_cm
  Sc_cm2_raw <- K_eff * stem_surface_full_cm2 - Swedges_cm2
  
  # Ensure non-negative results
  Vc_cm3 <- ifelse(is.finite(Vc_cm3_raw), pmax(Vc_cm3_raw, 0), NA_real_)
  Sc_cm2 <- ifelse(is.finite(Sc_cm2_raw), pmax(Sc_cm2_raw, 0), NA_real_)
  
  # Validation metrics
  denom_surface <- (stem_surface_full_cm2 - Swedges_cm2)
  surface_utilization_pct <- ifelse(is.finite(denom_surface) & (denom_surface > 0),
                                    (Sc_cm2 / denom_surface) * 100, NA_real_)
  volume_per_area_ratio <- ifelse(is.finite(Sc_cm2) & (Sc_cm2 > 0),
                                  Vc_cm3 / Sc_cm2, NA_real_)
  
  # Create results dataframe
  results <- tibble(
    row_number = seq_len(nrow(df)),
    valid_inputs = valid_inputs,
    Dstem_cm = Dstem_cm,
    T_cm = T_cm,
    L_cm = L_cm,
    H_cm = H_cm,
    Swedges_cm2 = Swedges_cm2,
    Vwedges_cm3 = Vwedges_cm3,
    wedge_source = wedge_source,
    Dext_cm = Dext_cm,
    circumference_cm = circumference_cm,
    K_raw = K_raw,
    K_eff = K_eff,
    wraparound_occurred = wraparound_flag,
    Sc_cm2_corrected = Sc_cm2,
    Vc_cm3_corrected = Vc_cm3,
    Sc_m2_corrected = Sc_cm2 / 10000,
    max_possible_Sc_cm2 = stem_surface_full_cm2 - Swedges_cm2,
    surface_utilization_pct = surface_utilization_pct,
    volume_per_area_ratio = volume_per_area_ratio
  )
  
  # the sheet's own values, for comparison
  Sc_orig <- num(df[[cols$Sc_orig]])
  Vc_orig <- num(df[[cols$Vc_orig]])
  results %>%
    mutate(
      Sc_cm2_original = Sc_orig,
      Sc_difference = Sc_cm2_corrected - Sc_orig,
      Sc_rel_error_pct = ifelse(is.finite(Sc_orig) & (Sc_orig != 0),
                                100 * (Sc_orig - Sc_cm2_corrected) / Sc_cm2_corrected,
                                NA_real_),
      Vc_cm3_original = Vc_orig,
      Vc_difference = Vc_cm3_corrected - Vc_orig,
      Vc_rel_error_pct = ifelse(is.finite(Vc_orig) & (Vc_orig != 0),
                                100 * (Vc_orig - Vc_cm3_corrected) / Vc_cm3_corrected,
                                NA_real_)
    )
}

generate_summary <- function(all_results) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n", sep = "")
  cat("TREE CHAMBER GEOMETRY ANALYSIS SUMMARY\n")
  cat(paste(rep("=", 60), collapse=""), "\n", sep = "")
  
  total_samples   <- nrow(all_results)
  valid_samples   <- sum(all_results$valid_inputs, na.rm = TRUE)
  wraparound_cases<- sum(all_results$wraparound_occurred, na.rm = TRUE)
  wrap_pct        <- if (valid_samples > 0) round(wraparound_cases / valid_samples * 100, 1) else NA
  
  cat("Sample counts:\n")
  cat("  Total samples:", total_samples, "\n")
  cat("  Valid samples:", valid_samples, "\n")
  cat("  Invalid/missing data:", total_samples - valid_samples, "\n")
  cat("  Wraparound cases (K > 1):", wraparound_cases, "\n")
  cat("  Wraparound percentage:", ifelse(is.na(wrap_pct), "NA", paste0(wrap_pct, "%")), "\n\n")
  
  cat("Sector fraction (K) statistics:\n")
  cat("  Raw K - Mean:", round(mean(all_results$K_raw, na.rm = TRUE), 3), "\n")
  cat("  Raw K - Range:",
      round(min(all_results$K_raw, na.rm = TRUE), 3), "to",
      round(max(all_results$K_raw, na.rm = TRUE), 3), "\n")
  cat("  Effective K - Mean:", round(mean(all_results$K_eff, na.rm = TRUE), 3), "\n\n")
  
  cat("Corrected surface area (cm²):\n")
  cat("  Mean:", round(mean(all_results$Sc_cm2_corrected, na.rm = TRUE), 1), "\n")
  cat("  Range:",
      round(min(all_results$Sc_cm2_corrected, na.rm = TRUE), 1), "to",
      round(max(all_results$Sc_cm2_corrected, na.rm = TRUE), 1), "\n\n")
  
  cat("Corrected chamber volume (cm³):\n")
  cat("  Mean:", round(mean(all_results$Vc_cm3_corrected, na.rm = TRUE), 1), "\n")
  cat("  Range:",
      round(min(all_results$Vc_cm3_corrected, na.rm = TRUE), 1), "to",
      round(max(all_results$Vc_cm3_corrected, na.rm = TRUE), 1), "\n\n")
  
  cat("Surface area, sheet vs corrected (relative difference):\n")
  cat("  Mean:", round(mean(all_results$Sc_rel_error_pct, na.rm = TRUE), 1), "%\n")
  cat("  Range:",
      round(min(all_results$Sc_rel_error_pct, na.rm = TRUE), 1), "% to",
      round(max(all_results$Sc_rel_error_pct, na.rm = TRUE), 1), "%\n\n")
  cat("Volume, sheet vs corrected (relative difference):\n")
  cat("  Mean:", round(mean(all_results$Vc_rel_error_pct, na.rm = TRUE), 1), "%\n")
  cat("  Range:",
      round(min(all_results$Vc_rel_error_pct, na.rm = TRUE), 1), "% to",
      round(max(all_results$Vc_rel_error_pct, na.rm = TRUE), 1), "%\n\n")
  
  # Validation warnings
  negative_volumes     <- sum(all_results$Vc_cm3_corrected < 0, na.rm = TRUE)
  negative_surfaces    <- sum(all_results$Sc_cm2_corrected < 0, na.rm = TRUE)
  over_utilization     <- sum(all_results$surface_utilization_pct > 100, na.rm = TRUE)
  
  if (negative_volumes > 0 || negative_surfaces > 0 || over_utilization > 0) {
    cat("Validation warnings:\n")
    if (negative_volumes  > 0) cat("  ", negative_volumes,  "samples with negative volumes\n")
    if (negative_surfaces > 0) cat("  ", negative_surfaces, "samples with negative surface areas\n")
    if (over_utilization  > 0) cat("  ", over_utilization,  "samples with >100% surface utilization\n")
  }
  cat(paste(rep("=", 60), collapse=""), "\n", sep = "")
}

# One workbook: every non-empty sheet with its corrected geometry appended
process_file <- function(file_path) {
  file_name <- basename(file_path)
  cat("\nProcessing:", file_name, "\n")
  sheet_names <- excel_sheets(file_path)
  all_sheets <- map(sheet_names, function(sheet) {
    df <- read_excel(file_path, sheet = sheet, col_types = "text")
    df <- df[!apply(is.na(df) | df == "", 1, all), , drop = FALSE]   # drop empty rows
    if (nrow(df) == 0) {
      message("  (empty sheet skipped): ", sheet)
      return(NULL)
    }
    cbind(df, calculate_chamber_geometry(df))
  })
  names(all_sheets) <- sheet_names
  all_sheets <- compact(all_sheets)
  stopifnot(length(all_sheets) > 0)
  output_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(file_name), "_corrected.xlsx"))
  write_xlsx(all_sheets, output_file)
  cat("  Written:", output_file, "\n")
  result <- bind_rows(all_sheets, .id = "sheet")
  result$source_file <- file_name
  result
}

cat("Chamber geometry (Siegenthaler et al. 2016)\n")
all_results <- map(file.path(input_dir, files), process_file) %>%
  bind_rows(.id = "file_index")
generate_summary(all_results)
combined_output <- file.path(output_dir, "combined_corrected_results.xlsx")
write_xlsx(list(all_data = all_results), combined_output)
cat("\nCombined results written:", combined_output, "\n")
