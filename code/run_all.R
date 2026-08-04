#!/usr/bin/env Rscript
# ==============================================================================
# Revision pipeline. Run from the repo root:  Rscript code/run_all.R
# Reads data/ (the Zenodo drop-in) + code/ (git); writes outputs/ and
# outputs/figures/{main,SI}/.
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
LOGDIR <- "outputs/logs"
dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)

# Every script's console output is captured to outputs/logs/<name>.txt.
# Nine scripts used to DECLARE a .txt output in their headers that nothing ever wrote:
# they only cat() to the console, and system2(stdout = "") let it fall through to the
# terminal. The .txt files on disk had been produced by hand-redirection on one day in
# July and could never be refreshed, so four referee-facing audits were frozen against
# a superseded model while appearing to be pipeline products. Capturing centrally fixes
# the whole class at once and gives every step a transcript.
run <- function(f, fatal = FALSE) {
  cat("\n>>>", basename(f), "\n")
  logf <- file.path(LOGDIR, sub("\\.R$", ".txt", basename(f)))
  st <- tryCatch(system2("Rscript", f, stdout = logf, stderr = logf),
                 warning = function(w) 1L, error = function(e) 1L)
  if (file.exists(logf)) cat(readLines(logf, warn = FALSE), sep = "\n")
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
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
           "outputs/.pipeline_run_started")

# --- model prerequisite -------------------------------------------------------
# This pipeline SCORES and CONSUMES the locked forests; it does not fit them. The only
# producer of RF_MODELS.RData / TRAINING_DATA.RData is code/05_model/02_rf_models.R,
# which is deliberately not run here because it retrains (~10 min) and the model is
# frozen. Five fatal CORE steps load those files, so on a clean checkout the run would
# otherwise die deep in the chain with an opaque error. Fail early and say why.
local({
  need <- c("outputs/models/RF_MODELS.RData", "outputs/models/TRAINING_DATA.RData")
  miss <- need[!file.exists(need)]
  if (length(miss))
    stop("missing locked model file(s):\n  ", paste(miss, collapse = "\n  "),
         "\n\nBuild them first (from code/05_model/):\n",
         "  Rscript 01_load_and_prep_data.R && Rscript 02_rf_models.R\n",
         call. = FALSE)
  src <- "code/05_model/02_rf_models.R"
  if (file.exists(src) && file.mtime(src) > min(file.mtime(need)))
    cat("[note] 02_rf_models.R is NEWER than the locked model files;",
        "the model may need rebuilding.\n")
})

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
  "code/05_model/rev_rf_grouped_cv.R",             # -> rf_grouped_cv.csv (budget reads it)
  "code/01_import/rev_inventory_build.R",           # raw -> inventory_stems.csv
  "code/04_drivers/rev_wb_reference_et.R",           # -> water balance (climatology input)
  "code/04_drivers/rev_moisture_climatology.R",      # -> moisture_climatology_monthly.csv
  "code/04_drivers/rev_soil_temp_climatology.R",     # -> soil_temp_climatology_monthly.csv
  "code/04_drivers/rev_moisture_surface.R",          # -> moisture_surface_grid.csv
  "code/06_upscale/rev_predict_tree_flux_current.R", # -> tree_flux_predictions.csv, tree_monthly_stand.csv
  "code/06_upscale/rev_predict_soil_surface.R",      # -> soil_surface_{monthly,annual}.csv
  "code/06_upscale/rev_budget_canonical.R",          # -> canonical_{budget,monthly}.csv
  # MUST precede the grid: rev_scaling_full_grid.R now stop()s if wai_bottomup.csv is
  # absent, and outputs/ is gitignored, so with this in SUPPORT (which runs AFTER the
  # fatal CORE block) a clean checkout aborted at the last CORE step. That is the exact
  # failure this file's header describes for the driver builders.
  "code/06_upscale/rev_wai_bottomup_and_rf_interactions.R",  # -> wai_bottomup.csv
  "code/06_upscale/rev_scaling_full_grid.R")
# fig_scaling_profiles / heatmap read the grid exports and run in the figure block         # -> scaling_full_grid.csv
cat(sprintf("\n== CORE CHAIN (%d steps, dependency-ordered) ==\n", length(CORE)))
for (f in CORE) run(f, fatal = TRUE)

# --- 2) supporting analyses (produce CSV/TXT that figures and prose cite) ----
SUPPORT <- c(
  # Produces data/processed/environmental/soil_env_by_collar.csv, which
  # code/05_model/01_load_and_prep_data.R reads to give each soil collar and each
  # monthly tree its OWN measured temperature and moisture rather than a plot-level
  # constant. It was never in this pipeline, so that dependency was real but unwired
  # and the CSV survived only from a manual run on 2026-07-25.
  "code/03_merge/rev_compile_soil_env.R",
  "code/04_drivers/rev_moisture_elevation_check.R",
  "code/04_drivers/rev_moisture_interpolation.R",
  "code/02_flux/rev_qc_c0_screen.R",
  "code/02_flux/rev_mdf_FINAL_precision_and_detection.R",
  "code/06_upscale/rev_surface_area_model.R",
  # Promoted out of exploratory/ 2026-07-29. Figure 6 panel (b) reads
  # outputs/data/FAPROTAX_all_functions_HW_SW.csv, and this is its ONLY
  # producer -- but it lived in exploratory/, which the glob below never reaches
  # (non-recursive, by design). So rev_fig06_hydrogenotrophy.R aborted on every
  # run since 2026-07-23, and because rev_00_assemble_figures.R tests only
  # file.exists() and never mtime, the assembler copied a stale PNG and reported
  # success. A load-bearing producer must not sit in a directory documented as
  # not-run. It reads raw data only, so it has no ordering constraint beyond
  # preceding the figure block.
  "code/07_molecular/rev_faprotax_dump_HW_SW.R",
  "code/06_upscale/rev_area_distribution_scenarios.R",
  "code/05_model/rev_height_form_crossvalidation.R",
  "code/05_model/rev_rf_model_diagnostics.R",
  "code/05_model/rev_rf_species_fallback_loso.R",
  "code/05_model/rev_rf_species_pooling.R",
  "code/05_model/rev_rf_species_bias_audit.R",
  "code/08_figures/rev_figS21_rf-model-summary.R",   # absorbed rev_fig_model_findings.R
  "code/05_model/rev_model_family_comparison.R",
  "code/05_model/rev_rf_height_extrapolation.R",
  "code/06_upscale/rev_scaling_assumptions_audit.R",
  # Referee-facing evidence produced in the 2026-07-30 pass. Both were written but
  # never wired in, which is the same defect this file exists to prevent.
  "code/05_model/rev_rf_predictor_selection_current.R",  # -> rf_predictor_selection_current.csv
  "code/05_model/rev_rf_calibration_sensitivity.R")      # -> rf_calibration_sensitivity.csv
SUPPORT <- SUPPORT[file.exists(SUPPORT)]
cat(sprintf("\n== SUPPORTING ANALYSES (%d) ==\n", length(SUPPORT)))
for (f in SUPPORT) run(f)

# --- 3) stats, tables and figures --------------------------------------------
# NAMED, NOT GLOBBED. This was list.files("code/revision", "^rev_(stat|tbl|fig)")
# until the reorg. A glob that selects on a filename prefix silently empties when
# the prefix changes: retiring "rev_" -- an explicit goal of the reorg -- would
# have dropped all 41 scripts below while this file still exited 0. The gate would
# not have caught it either, because rev_check_consistency.R asserts agreement
# BETWEEN outputs, not that they were regenerated, so the stale CSVs already on
# disk keep all 28 invariants passing. Select on meaning, never on spelling.
# This list was generated from the glob it replaces and verified set-identical.
rest <- c(
  "code/08_figures/rev_fig01_temporal-flux.R",
  "code/08_figures/rev_fig02_height-flux.R",
  "code/08_figures/rev_fig02a_axis-support.R",
  "code/08_figures/rev_fig03_variance-partition.R",
  "code/08_figures/rev_fig04_gene-abundance.R",
  "code/08_figures/rev_fig05_methane-cycling.R",
  "code/08_figures/rev_fig06_hydrogenotrophy.R",
  "code/08_figures/rev_fig07_decay-methanogenesis.R",
  "code/08_figures/rev_fig07b_copies-per-g.R",
  "code/08_figures/rev_fig07c_flux-unit.R",
  "code/08_figures/rev_fig09_budget.R",
  "code/08_figures/rev_figS02_height-slope-moisture.R",
  "code/08_figures/rev_figS04_pmoa-mmox-coupling.R",
  "code/08_figures/rev_figS11_scale-dependent.R",
  "code/08_figures/rev_figS12_isotope-sources.R",
  "code/08_figures/rev_figS15_black-oak-methanome.R",
  "code/08_figures/rev_figS17_plant-traits.R",
  "code/08_figures/rev_figS19_mcra-probe-validation.R",
  "code/08_figures/rev_figS20_stem-deterioration.R",
  "code/08_figures/rev_figSI_detection.R",
  "code/08_figures/rev_figS_black-oak-cross-sections.R",
  "code/08_figures/rev_figS_rf-calibration.R",
  "code/08_figures/rev_fig_height_curves.R",
  "code/08_figures/rev_fig_scaling_diagnostics.R",
  "code/08_figures/rev_fig_scaling_heatmap.R",
  "code/08_figures/rev_fig_scaling_profiles.R",
  "code/09_tables_stats/rev_stat_campaign_counts.R",
  "code/09_tables_stats/rev_stat_copies-per-gram.R",
  "code/09_tables_stats/rev_stat_dbh_by_species_campaign.R",
  "code/09_tables_stats/rev_stat_faprotax-caveats.R",
  "code/09_tables_stats/rev_stat_isotopes-canonical.R",
  "code/09_tables_stats/rev_stat_known-putative-table.R",
  "code/09_tables_stats/rev_stat_multigene-models.R",
  "code/09_tables_stats/rev_stat_pmoa-mmox-robustness.R",
  "code/09_tables_stats/rev_stat_s1-rf-soil-arcsinh.R",
  "code/09_tables_stats/rev_stat_s1s2-arcsinh.R",
  "code/09_tables_stats/rev_stat_species-aggregation-rma.R",
  "code/09_tables_stats/rev_stat_tree-distribution.R",
  "code/09_tables_stats/rev_stat_tree_flux_merged.R",
  "code/09_tables_stats/rev_stat_variance-partition.R",
  "code/09_tables_stats/rev_tbl_ddpcr-16s-concordance.R"
)
rest <- setdiff(rest, c(CORE, SUPPORT))
local({
  missing <- rest[!file.exists(rest)]
  if (length(missing))
    stop("named script(s) not found -- did a move miss this list?\n  ",
         paste(missing, collapse = "\n  "), call. = FALSE)
})
cat(sprintf("\n== STATS, TABLES AND FIGURES (%d) ==\n", length(rest)))
for (f in rest) run(f)

# --- 3b) regenerate the parameter record from the canonical outputs ----------
# scaling_parameters.md section 0 claims "nothing here is typed by hand". This is what
# makes that true; without it the claim drifted through four rounds of changes.
run("code/09_tables_stats/rev_write_parameter_record.R")

# --- 4) assemble -------------------------------------------------------------
cat("\n>>> assembling numbered manuscript figures\n")
source("code/08_figures/rev_00_assemble_figures.R")

# --- 5) report anything never reached ----------------------------------------
SOURCED <- c("code/lib/rev_geometry.R",        # sourced by others, not run alone
             "code/lib/rev_species_levels.R",  # ditto -- the species->level mapping
             "code/lib/rev_prep_species_data.R")
# Scans the whole live tree, not one directory. This used to read
# list.files("code/revision", ...); the 2026-08 reorg emptied that directory of
# scripts, so the reachability check silently matched nothing and reported a
# clean bill on every run -- the exact failure it exists to catch, applied to
# itself. archive/ is excluded by design; lib/ is covered by SOURCED below.
allR <- setdiff(
  list.files("code", "\\.R$", recursive = TRUE, full.names = TRUE),
  list.files("code/archive", "\\.R$", recursive = TRUE, full.names = TRUE))
never <- setdiff(allR, c(CORE, SUPPORT, rest, SOURCED,
                         "code/run_all.R",
                         "code/08_figures/rev_00_assemble_figures.R",
                         # Run standalone at step 3b, not via CORE/SUPPORT/rest, so it
                         # was absent here and got reported as "never run" on every
                         # pass -- immediately after this script had just run it. The
                         # 2026-07-31 reorg audit nearly archived it on that evidence;
                         # it is the only producer of scaling_parameters.md section 0.
                         "code/09_tables_stats/rev_write_parameter_record.R",
                         # The gate. Deliberately not part of the pipeline (it checks
                         # agreement BETWEEN outputs, so it runs after, by hand), but
                         # it is not dead either.
                         "code/check_consistency.R"))
if (length(never)) {
  cat(sprintf("\n[note] %d script(s) in code/revision/ were not run by this pipeline:\n",
              length(never)))
  cat(paste("   -", basename(never), collapse = "\n"), "\n")
  cat("   (these are superseded or one-off; add to SUPPORT above if they should run)\n")
}
cat("\nDone. Numbered figures in outputs/figures/{main,SI}/ (see MANIFEST.md).\n")
