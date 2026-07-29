#!/usr/bin/env Rscript
# ==============================================================================
# Revision pipeline. Run from the repo root:  Rscript code/revision/run_all.R
# Reads data/ (the Zenodo drop-in) + code/ (git); writes outputs/revision/ and
# outputs/revision/figures/{main,SI}/.
#
# ORDER IS EXPLICIT, NOT A GLOB. This script used to collect work by filename
# pattern -- ^rev_stat_, ^rev_tbl_, ^rev_fig -- which silently skipped 44 of the
# ~90 scripts here, including all three that produce Figure 9's inputs. On a
# clean checkout it reached rev_fig09_budget.R with no canonical_budget.csv, no
# scaling_full_grid.csv and no tree_flux_predictions.csv, and only appeared to
# work because those files were sitting in an untracked directory from earlier
# manual runs. The core chain below is ordered by dependency and runs first;
# anything matching the old patterns and not already named runs afterwards, and
# whatever is still never reached is reported at the end.
# Exploratory scripts (code/revision/exploratory/) are NOT run.
# ==============================================================================
run <- function(f, fatal = FALSE) {
  cat("\n>>>", basename(f), "\n")
  st <- tryCatch(system2("Rscript", f, stdout = "", stderr = ""),
                 warning = function(w) 1L, error = function(e) 1L)
  if (!identical(st, 0L)) {
    cat("   [non-zero exit]\n")
    if (fatal) stop("required step failed: ", basename(f), call. = FALSE)
  }
  invisible(st)
}

# --- run marker ---------------------------------------------------------------
# rev_00_assemble_figures.R compares every file it assembles against this marker's
# timestamp, so a generator that failed cannot slip its previous output into the
# manuscript set unnoticed. Written before anything else runs.
dir.create("outputs/revision", showWarnings = FALSE, recursive = TRUE)
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
           "outputs/revision/.pipeline_run_started")

# --- 0) original-pipeline figures the assembler copies through ---------------
# Two fail here by design (fig5, S12 methanome); both are replaced by revision
# versions, so the assembler still finds them.
cat(">>> generate_all_figures.R (original pipeline -- unchanged figures)\n")
run("code/generate_all_figures.R")

# --- 1) CORE CHAIN: each step consumes the previous one's output -------------
# inventory -> per-stem tree flux -> soil surface -> budget -> scaling grid.
# Fatal, because everything downstream reads what these write.
# The DRIVER builders come first. They used to sit in SUPPORT, which runs AFTER
# this block, so the chain only worked because their outputs were already on disk
# from a previous run -- a clean checkout would have stopped at the first
# prediction script. Both prediction scripts now read the same two climatologies
# and the same moisture surface, so those have to be built before either runs.
CORE <- c(
  "code/revision/rev_rf_grouped_cv.R",             # -> rf_grouped_cv.csv (budget reads it)
  "code/revision/rev_inventory_build.R",           # raw -> inventory_stems.csv
  "code/revision/rev_wb_reference_et.R",           # -> water balance (climatology input)
  "code/revision/rev_moisture_climatology.R",      # -> moisture_climatology_monthly.csv
  "code/revision/rev_soil_temp_climatology.R",     # -> soil_temp_climatology_monthly.csv
  "code/revision/rev_moisture_surface.R",          # -> moisture_surface_grid.csv
  "code/revision/rev_predict_tree_flux_current.R", # -> tree_flux_predictions.csv, tree_monthly_stand.csv
  "code/revision/rev_predict_soil_surface.R",      # -> soil_surface_{monthly,annual}.csv
  "code/revision/rev_budget_canonical.R",          # -> canonical_{budget,monthly}.csv
  "code/revision/rev_scaling_full_grid.R")
# fig_scaling_profiles / heatmap read the grid exports and run in the figure block         # -> scaling_full_grid.csv
cat(sprintf("\n== CORE CHAIN (%d steps, dependency-ordered) ==\n", length(CORE)))
for (f in CORE) run(f, fatal = TRUE)

# --- 2) supporting analyses (produce CSV/TXT that figures and prose cite) ----
SUPPORT <- c(
  "code/revision/rev_moisture_elevation_check.R",
  "code/revision/rev_moisture_interpolation.R",
  "code/revision/rev_qc_c0_screen.R",
  "code/revision/rev_mdf_FINAL_precision_and_detection.R",
  "code/revision/rev_surface_area_model.R",
  # Promoted out of exploratory/ 2026-07-29. Figure 6 panel (b) reads
  # outputs/revision/FAPROTAX_all_functions_HW_SW.csv, and this is its ONLY
  # producer -- but it lived in exploratory/, which the glob below never reaches
  # (non-recursive, by design). So rev_fig06_hydrogenotrophy.R aborted on every
  # run since 2026-07-23, and because rev_00_assemble_figures.R tests only
  # file.exists() and never mtime, the assembler copied a stale PNG and reported
  # success. A load-bearing producer must not sit in a directory documented as
  # not-run. It reads raw data only, so it has no ordering constraint beyond
  # preceding the figure block.
  "code/revision/rev_faprotax_dump_HW_SW.R",
  "code/revision/rev_wai_bottomup_and_rf_interactions.R",
  "code/revision/rev_area_distribution_scenarios.R",
  "code/revision/rev_height_form_crossvalidation.R",
  "code/revision/rev_rf_model_diagnostics.R",
  "code/revision/rev_rf_species_fallback_loso.R",
  "code/revision/rev_rf_species_pooling.R",
  "code/revision/rev_rf_species_bias_audit.R",
  "code/revision/rev_fig_model_findings.R",
  "code/revision/rev_model_family_comparison.R",
  "code/revision/rev_rf_height_extrapolation.R",
  "code/revision/rev_scaling_assumptions_audit.R")
SUPPORT <- SUPPORT[file.exists(SUPPORT)]
cat(sprintf("\n== SUPPORTING ANALYSES (%d) ==\n", length(SUPPORT)))
for (f in SUPPORT) run(f)

# --- 3) everything else matching the historical patterns ---------------------
rest <- sort(list.files("code/revision", "^rev_(stat|tbl|fig).*\\.R$", full.names = TRUE))
rest <- setdiff(rest, c(CORE, SUPPORT))
cat(sprintf("\n== STATS, TABLES AND FIGURES (%d) ==\n", length(rest)))
for (f in rest) run(f)

# --- 4) assemble -------------------------------------------------------------
cat("\n>>> assembling numbered manuscript figures\n")
source("code/revision/rev_00_assemble_figures.R")

# --- 5) report anything never reached ----------------------------------------
SOURCED <- c("code/revision/rev_geometry.R",        # sourced by others, not run alone
             "code/revision/rev_prep_species_data.R")
allR  <- list.files("code/revision", "\\.R$", full.names = TRUE)
never <- setdiff(allR, c(CORE, SUPPORT, rest, SOURCED,
                         "code/revision/run_all.R",
                         "code/revision/rev_00_assemble_figures.R"))
if (length(never)) {
  cat(sprintf("\n[note] %d script(s) in code/revision/ were not run by this pipeline:\n",
              length(never)))
  cat(paste("   -", basename(never), collapse = "\n"), "\n")
  cat("   (these are superseded or one-off; add to SUPPORT above if they should run)\n")
}
cat("\nDone. Numbered figures in outputs/revision/figures/{main,SI}/ (see MANIFEST.md).\n")
