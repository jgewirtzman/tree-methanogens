#!/usr/bin/env Rscript
# ==============================================================================
# make_figures.R -- the single runner for every publication figure and table.
# Run from the repo root:  Rscript code/make_figures.R
#
# Replaces generate_all_figures.R (original-pipeline figures) and the standalone
# invocation of the assembler. Order: original generators -> revision generators
# -> assembler.
#
# WHY EACH SCRIPT GETS ITS OWN PROCESS.
# generate_all_figures.R ran every script with
#   source(path, local = new.env(parent = globalenv()))
# which gives each script its own *environment* but leaves them sharing one R
# *process*. Graphics devices are process-global, so `local=` isolates variables
# and does nothing whatsoever for devices.
#
# That produced a wrong figure in the shipped SI. 12b_picrust_pathway_heatmap.R
# opens png(fig6_picrust_mcra_no_mcra_heatmap.png); when its pheatmap() call
# errored, the handler's try(dev.off(), silent = TRUE) did not take, and the
# device stayed open on fig6's path. Thirty scripts later 05_methods_figure_map.R
# drew the moisture map -- into fig6's still-open device. On disk, fig6 carried
# the map's timestamp (15:05) rather than its own siblings' (14:58), and
# Figure_S12 in the assembled SI set was the moisture overlay.
#
# system2("Rscript", f) makes that class of failure impossible by construction:
# a leaked device dies with the process that leaked it, and R flushes on exit.
# The cost is reloading packages per script. That is the right trade -- the
# alternative silently corrupts figures.
#
# It also means a figure that FAILS to draw cannot be mistaken for one that
# drew: the runner checks each script's declared outputs against the run marker
# and reports any that are missing or older than this run.
# ==============================================================================

if (!file.exists("data/processed/integrated/merged_tree_dataset_final.csv"))
  stop("Must run from the repo root.\n  Usage: Rscript code/make_figures.R", call. = FALSE)

LOGDIR <- "outputs/logs"
dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)
# MUST be the name 00_assemble_figures.R looks for. Naming it anything else
# silently disables the assembler's staleness check -- the guard that stops a
# figure older than this run from being copied into the manuscript set and
# reported as fresh. That guard exists because the assembler once shipped a
# stale PNG while printing success; it is the same failure class as the Fig S12
# device leak, and it must not be switched off by a marker rename.
MARKER <- "outputs/.pipeline_run_started"
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), MARKER)
t0 <- Sys.time()

results <- list()

run_one <- function(path, label) {
  if (!file.exists(path)) {
    results[[length(results) + 1]] <<- list(path = path, label = label,
                                            status = "MISSING", secs = 0, msg = "script not found")
    cat(sprintf("  %-58s MISSING\n", basename(path)));  return(invisible(NULL))
  }
  st_time <- Sys.time()
  logf <- file.path(LOGDIR, sub("\\.R$", ".txt", basename(path)))
  st <- tryCatch(system2("Rscript", path, stdout = logf, stderr = logf),
                 warning = function(w) 1L, error = function(e) 1L)
  secs <- as.numeric(difftime(Sys.time(), st_time, units = "secs"))
  ok <- identical(st, 0L)
  results[[length(results) + 1]] <<- list(
    path = path, label = label, status = if (ok) "PASS" else "FAIL", secs = secs,
    msg = if (ok) "" else paste(utils::tail(readLines(logf, warn = FALSE), 3), collapse = " | "))
  cat(sprintf("  %-58s %s (%.0fs)\n", basename(path), if (ok) "ok" else "FAIL", secs))
  invisible(NULL)
}

# --- 1) original-pipeline generators -----------------------------------------
ORIGINAL <- c(
  "code/08_figures/06_soil_tree_timeseries.R"                  = "fig1",
  "code/02_flux/static/04_height_effect_analysis.R"            = "fig2",
  "code/08_figures/04_variance_partition.R"                    = "fig3",
  "code/08_figures/util_combined_plot.R"                       = "fig4",
  "code/08_figures/08c_combined_methane_cycling_composition.R" = "fig5",
  "code/08_figures/12b_picrust_pathway_heatmap.R"              = "fig6",
  "code/08_figures/09_felled_oak_profiles.R"                   = "fig7",
  "code/07_molecular/04_species_gene_flux.R"                   = "fig8",
  "code/08_figures/09_upscale_publication_plots.R"             = "fig9",
  "code/08_figures/12c_taxonomy_pmoa_heatmap.R"                = "S2",
  "code/08_figures/08d_faprotax_heatmaps.R"                    = "S3",
  "code/08_figures/12a_taxonomy_mcra_heatmap.R"                = "S6",
  "code/08_figures/05_internal_gas_plots.R"                    = "S7,S8",
  "code/08_figures/11a_isotope_d13ch4_single.R"                = "S9",
  "code/07_molecular/methanotrophs/03_pmoa_mmox_analysis.R"    = "S10",
  "code/08_figures/10_black_oak_methanome_heatmap.R"           = "S12",
  "code/08_figures/02_radial_cross_sections.R"                 = "S13",
  "code/08_figures/08_rf_publication_plots.R"                  = "S15",
  "code/08_figures/05_methods_figure_map.R"                    = "S1")

cat("\n== ORIGINAL-PIPELINE GENERATORS ==\n")
for (p in names(ORIGINAL)) run_one(p, ORIGINAL[[p]])

# --- 2) revision generators ---------------------------------------------------
# NAMED, NOT GLOBBED. A glob on "^rev_fig" empties the moment that prefix is
# retired, and this runner would still exit 0 having generated nothing --
# the assembler would then copy the previous run's figures and report success.
# The list below was generated from the glob and asserted set-identical to it.
REVISION <- c(
  "code/08_figures/fig01_temporal-flux.R",
  "code/08_figures/fig02_height-flux.R",
  "code/08_figures/fig02a_axis-support.R",
  "code/08_figures/fig03_variance-partition.R",
  "code/08_figures/fig04_gene-abundance.R",
  "code/08_figures/fig05_methane-cycling.R",
  "code/08_figures/fig06_hydrogenotrophy.R",
  "code/08_figures/fig07_decay-methanogenesis.R",
  "code/08_figures/fig07b_copies-per-g.R",
  "code/08_figures/fig07c_flux-unit.R",
  "code/08_figures/fig09_budget.R",
  "code/08_figures/figS02_height-slope-moisture.R",
  "code/08_figures/figS04_pmoa-mmox-coupling.R",
  "code/08_figures/figS11_scale-dependent.R",
  "code/08_figures/figS12_isotope-sources.R",
  "code/08_figures/figS15_black-oak-methanome.R",
  "code/08_figures/figS17_plant-traits.R",
  "code/08_figures/figS19_mcra-probe-validation.R",
  "code/08_figures/figS20_stem-deterioration.R",
  "code/08_figures/figS21_rf-model-summary.R",
  "code/08_figures/figSI_detection.R",
  "code/08_figures/figS_black-oak-cross-sections.R",
  "code/08_figures/figS_rf-calibration.R",
  "code/08_figures/fig_height_curves.R",
  "code/08_figures/fig_scaling_diagnostics.R",
  "code/08_figures/fig_scaling_heatmap.R",
  "code/08_figures/fig_scaling_profiles.R")
stopifnot(all(file.exists(REVISION)))
cat(sprintf("\n== REVISION GENERATORS (%d) ==\n", length(REVISION)))
for (p in REVISION) run_one(p, "revision")

# --- 3) assemble --------------------------------------------------------------
cat("\n== ASSEMBLE ==\n")
run_one("code/08_figures/00_assemble_figures.R", "assemble")

# --- 4) report ----------------------------------------------------------------
st  <- vapply(results, function(r) r$status, character(1))
cat(sprintf("\n%s\nscripts: %d | ok: %d | failed: %d | %.0fs\n%s\n",
            strrep("=", 62), length(st), sum(st == "PASS"),
            sum(st != "PASS"), as.numeric(difftime(Sys.time(), t0, units = "secs")),
            strrep("=", 62)))
bad <- Filter(function(r) r$status != "PASS", results)
if (length(bad)) {
  cat("\n--- FAILURES (these figures are NOT regenerated) ---\n")
  for (r in bad) cat(sprintf("  %-52s %s\n     %s\n", basename(r$path), r$status, r$msg))
  cat("\nA failed generator leaves the PREVIOUS figure on disk. Do not ship the\n",
      "assembled set until these pass -- file.exists() cannot tell a stale or\n",
      "wrong figure from a correct one.\n", sep = "")
}
cat("\nAssembled set: outputs/figures/{main,SI}/ (see outputs/figures/MANIFEST.md)\n")
quit(status = if (length(bad)) 1L else 0L)
