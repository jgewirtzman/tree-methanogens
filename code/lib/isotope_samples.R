# ==============================================================================
# isotope_samples.R -- the one definition of the whole-tree delta13C sample set
# ------------------------------------------------------------------------------
# Used by 12_isotopes-canonical.R, figS17_isotope-sources.R,
# 11a_isotope_d13ch4_single.R and 20_manuscript_statistics.R. Until 2026-10-01 each
# selected samples its own way, and the paper quoted two sets that disagreed
# (n = 125, median -63.7 vs n = 130, median -63.0):
#   - three scripts kept only samples whose name matched a species in the ddPCR
#     metadata, dropping three real hemlocks (H6, H7, H304) and both redo runs;
#   - the canonical script kept a room-air vial ("ROOM AIR H8_NA": its name does
#     not start with "Amb") and counted WA104 twice (original and redo both pass
#     the 1.5 ppm floor).
#
# Rules: whole-tree single-borehole samples from every Picarro results run;
# exclude atmosphere and room air, incubations (V*), the paired heartwood/sapwood
# tissue set (names ending H/S), calibration standards; CH4 >= 1.5 ppm (reliable
# delta13C). No delta13C window is applied (2026-10-02): extreme values are handled
# by robust (MM) regression in the source estimate, medians and rank statistics
# elsewhere, with a Tukey 3 x IQR exclusion reported as a sensitivity check
# (12_isotopes-canonical.R). Where a "_Redo" exists, the redo replaces the original -- the same
# rule 03_process_internal_gas.R applies to the GC data. Species comes from the
# ddPCR metadata, else from the internal-gas table (same tree IDs).
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
  sp_dd <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv", stringsAsFactors = FALSE) %>%
    distinct(seq_id, species) %>% filter(!duplicated(seq_id))
  sp_gc <- read.csv("data/processed/internal_gas/sample_data_only.csv", stringsAsFactors = FALSE) %>%
    distinct(Tree.ID, Species.ID) %>% filter(!duplicated(Tree.ID))
  d %>% left_join(sp_dd, by = c("tree_id" = "seq_id")) %>%
    left_join(sp_gc, by = c("tree_id" = "Tree.ID")) %>%
    mutate(species = coalesce(species, Species.ID)) %>% select(-Species.ID)
}

# Atmosphere vials (for reference lines): "Amb..." with CH4 <= 5 ppm (higher = contaminated)
isotope_atmosphere_samples <- function(raw = picarro_runs()) {
  raw %>% transmute(SampleName, d13CH4 = HR_Delta_iCH4_Raw_mean, ch4_ppm = HR_12CH4_dry_mean) %>%
    filter(str_detect(SampleName, "^Amb"), !is.na(d13CH4), ch4_ppm >= ISO_CH4_FLOOR, ch4_ppm <= 5)
}
