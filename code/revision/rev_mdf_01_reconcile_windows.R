#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_01_reconcile_windows.R      STEP 1 of the empirical-MDF rebuild
# ------------------------------------------------------------------------------
# Build ONE master list of measurement windows covering every flux in the dataset.
#
# WHY THIS IS STEP 1
#   Empirical MDF needs, for each measurement, the raw 1 Hz trace over its closure
#   window plus the chamber geometry. Both exist, but the window definitions are
#   spread across five auxfiles written at different times by different scripts --
#   the original per-campaign processing plus the untagged/dead-snag rescue. Nothing
#   downstream can run until they are reconciled, so this script does only that and
#   reports coverage honestly.
#
# THE FIVE SOURCES
#   auxfile_goFlux_with_weather.csv       semirigid tree stems
#   auxfile_goFlux_soilflux_with_weather.csv  semirigid soil
#   goflux_auxfile.csv                    2021 height campaign (carries measurement_height)
#   ymf2023_goflux_auxfile.csv            2023 static chambers
#   outputs/revision/untagged_auxfile.csv rescued untagged / dead-snag stems
#
# All five carry the goFlux-required columns (start.time, Area, Vtot, Tcham, Pcham),
# so no window has to be reconstructed -- only harmonised and deduplicated.
#
# RAW TRACES
#   data/raw/lgr/semirigid_2020-2021   1970 files
#   data/raw/lgr/static_2021            132 files
#   data/raw/lgr/static_2023            199 files
#
# Output: outputs/revision/mdf_master_windows.csv
#         outputs/revision/mdf_window_coverage.txt
# ==============================================================================

suppressMessages({library(dplyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

pick <- function(d, ...) {           # first matching column name, else NA
  for (n in c(...)) if (n %in% names(d)) return(d[[n]])
  rep(NA, nrow(d))
}
load_aux <- function(path, campaign, raw_dir) {
  if (!file.exists(path)) { cat(sprintf("  MISSING: %s\n", path)); return(NULL) }
  d <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
  out <- data.frame(
    campaign   = campaign,
    UniqueID   = as.character(pick(d, "UniqueID", "unique_id")),
    start.time = as.character(pick(d, "start.time")),
    obs.length = suppressWarnings(as.numeric(pick(d, "obs.length"))),
    Area       = suppressWarnings(as.numeric(pick(d, "Area"))),
    Vtot       = suppressWarnings(as.numeric(pick(d, "Vtot"))),
    Tcham      = suppressWarnings(as.numeric(pick(d, "Tcham"))),
    Pcham      = suppressWarnings(as.numeric(pick(d, "Pcham"))),
    height_cm  = suppressWarnings(as.numeric(pick(d, "measurement_height", "measurement_height_cm"))),
    file_name  = as.character(pick(d, "File name", "file_name", "File.name")),
    raw_dir    = raw_dir,
    src        = basename(path), stringsAsFactors=FALSE)
  cat(sprintf("  %-44s n=%4d  start.time %3.0f%%  geom %3.0f%%  obs.length %3.0f%%\n",
      basename(path), nrow(out),
      100*mean(!is.na(out$start.time) & nzchar(out$start.time)),
      100*mean(!is.na(out$Area) & !is.na(out$Vtot)),
      100*mean(!is.na(out$obs.length))))
  out
}

rule("LOADING THE FIVE AUXFILES")
A <- bind_rows(
  load_aux("data/processed/flux/auxfile_goFlux_with_weather.csv",
           "semirigid_tree", "data/raw/lgr/semirigid_2020-2021"),
  load_aux("data/processed/flux/auxfile_goFlux_soilflux_with_weather.csv",
           "semirigid_soil", "data/raw/lgr/semirigid_2020-2021"),
  load_aux("data/processed/flux/goflux_auxfile.csv",
           "height_2021", "data/raw/lgr/static_2021"),
  load_aux("data/processed/flux/ymf2023_goflux_auxfile.csv",
           "static_2023", "data/raw/lgr/static_2023"),
  load_aux("outputs/revision/untagged_auxfile.csv",
           "untagged_rescue", "data/raw/lgr/semirigid_2020-2021"))

rule("HARMONISED MASTER LIST")
cat(sprintf("\n  total window definitions: %d\n\n", nrow(A)))
print(A %>% group_by(campaign) %>%
  summarise(n=n(), has_start=sum(!is.na(start.time) & nzchar(start.time)),
            has_geom=sum(!is.na(Area) & !is.na(Vtot)),
            has_obslen=sum(!is.na(obs.length)),
            has_height=sum(!is.na(height_cm)), .groups="drop") %>% as.data.frame(),
  row.names=FALSE)

dup <- A$UniqueID[duplicated(A$UniqueID) & !is.na(A$UniqueID)]
cat(sprintf("\n  duplicate UniqueIDs across sources: %d %s\n", length(unique(dup)),
    if (length(dup)) paste0("(", paste(head(unique(dup),4), collapse=", "), " ...)") else ""))

# ---- coverage against what actually needs MDF --------------------------------
rule("COVERAGE AGAINST THE ANALYSED DATA")
load("outputs/models/TRAINING_DATA.RData")
n_tree <- nrow(tree_train_complete); n_soil <- nrow(soil_train_complete)
cat(sprintf("
  tree training rows %d
  soil training rows %d
  master windows     %d
", n_tree, n_soil, nrow(A)))

have_mdf <- read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE)$UniqueID
have_mdf <- c(have_mdf, read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv",
                                 stringsAsFactors=FALSE)$UniqueID)
A$has_existing_mdf <- A$UniqueID %in% have_mdf
cat(sprintf("
  windows with an EXISTING goFlux MDF (manufacturer precision) : %d (%.0f%%)
  windows needing a fresh goFlux run                           : %d (%.0f%%)
", sum(A$has_existing_mdf), 100*mean(A$has_existing_mdf),
   sum(!A$has_existing_mdf), 100*mean(!A$has_existing_mdf)))
cat("\n  by campaign:\n\n")
print(A %>% group_by(campaign) %>%
  summarise(n=n(), existing_mdf=sum(has_existing_mdf), need_run=sum(!has_existing_mdf),
            .groups="drop") %>% as.data.frame(), row.names=FALSE)

# ---- raw trace availability ---------------------------------------------------
rule("RAW TRACE AVAILABILITY")
for (rd in unique(A$raw_dir)) {
  n <- length(list.files(rd, recursive=TRUE, pattern="\\.(txt|zip)$"))
  cat(sprintf("  %-38s %5d trace files   (%d windows point here)\n", rd, n, sum(A$raw_dir==rd)))
}

write.csv(A, file.path(outdir,"mdf_master_windows.csv"), row.names=FALSE)
cat(sprintf("\n  Written: outputs/revision/mdf_master_windows.csv  (%d windows)\n", nrow(A)))

rule("WHAT STEP 2 NEEDS")
cat("
  1. Match each window to its raw 1 Hz trace file and extract the closure segment.
     Windows carry start.time and obs.length; the semirigid traces are zipped and
     will need decompressing on the fly.
  2. Per-measurement Allan deviation at tau = 1 s, computed on the FIRST DIFFERENCE
     of the dry mole fraction so the flux trend does not inflate the noise estimate.
  3. Christiansen MDF = (SD_allan * 3 * t_alpha) / t * flux.term, at 90/95/99%.
  4. Re-emit C0 for every measurement, which closes the 61% gap that limited the
     task #16 contamination screen to the soil and height campaigns only.

  Any window that cannot be matched to a trace must be reported, not silently
  dropped -- an unmatched measurement is one we cannot assess, which is different
  from one that passed.
")
