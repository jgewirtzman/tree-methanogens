# ==============================================================================
# isotope_samples.R -- the one definition of the whole-tree delta13C sample set
# ------------------------------------------------------------------------------
# Applied once, by code/tables/33_isotopes.R, which reads the raw Picarro runs and flags
# the set in data/tables/isotopes.csv (in_whole_tree_set; tree_id via tree_tags,
# species via trees.csv). Analyses read that table (code/lib/tables_read.R), never this
# file.
#
# Rules: whole-tree single-borehole samples from every Picarro results run;
# exclude atmosphere and room air, incubations (V*), the paired heartwood/sapwood
# tissue set (names ending H/S), calibration standards; CH4 >= 1.5 ppm (reliable
# delta13C). No delta13C window is applied: extreme values are handled
# by robust (MM) regression in the source estimate, medians and rank statistics
# elsewhere, with a Tukey 3 x IQR exclusion reported as a sensitivity check
# (12_isotopes-canonical.R). Where a "_Redo" exists, the redo replaces the original -- the same
# rule 03_process_internal_gas.R applies to the GC data. Species is not attached here
# (it comes from trees.csv through the table).
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(purrr); library(stringr) })

ISO_STANDARDS <- c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c","SA1")
ISO_CH4_FLOOR <- 1.5

picarro_runs <- function(dir = "data/raw/internal_gas/picarro") {
  f <- list.files(dir, pattern = "_results.csv$", full.names = TRUE)
  map_dfr(f[!grepl("merged", f)], ~ read_csv(.x, show_col_types = FALSE))
}

isotope_whole_tree_samples <- function(raw = picarro_runs()) {
  d <- raw %>%
    transmute(SampleName,
              d13CH4 = HR_Delta_iCH4_Raw_mean, ch4_ppm = HR_12CH4_dry_mean,
              d13CO2 = Delta_Raw_iCO2_mean,     co2_ppm = `12CO2_mean`) %>%
    filter(!str_detect(SampleName, regex("^Amb|room\\s*air", ignore_case = TRUE)),
           !str_detect(SampleName, "^V"), !str_detect(SampleName, "[HS]$"),
           !SampleName %in% ISO_STANDARDS)
  redone <- sub("_Redo$", "", d$SampleName[grepl("_Redo$", d$SampleName)])
  d <- d %>% filter(!SampleName %in% redone) %>%                 # original replaced by its redo
    mutate(tree_id = sub("_Redo$", "", SampleName)) %>%
    filter(!is.na(d13CH4), ch4_ppm >= ISO_CH4_FLOOR)
  stopifnot(!anyDuplicated(d$tree_id))
  d
}

