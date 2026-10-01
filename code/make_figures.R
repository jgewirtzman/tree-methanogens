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
# That produced a wrong figure in the shipped SI. figS11-13_picrust-heatmaps.R
# opens png(fig6_picrust_mcra_no_mcra_heatmap.png); when its pheatmap() call
# errored, the handler's try(dev.off(), silent = TRUE) did not take, and the
# device stayed open on fig6's path. Thirty scripts later figS01_moisture-overlay.R
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
# MUST be the name zz_assemble_figures.R looks for. Naming it anything else
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
# 04_variance_partition.R is deliberately absent. It is the pre-revision Figure 3
# -- raw-scale OLS on measurement rows, no growing-season restriction -- and the
# assembler takes Figure 3 from fig03_variance-partition.R instead. Left in the
# runner it regenerated outputs/figures/original/main/fig3_variance_partitioning.png
# on every pass: a file that looks authoritative, is named like Figure 3, and
# reports 82.9% unexplained against the current 65.3%. Archived.
source("code/lib/figure_scripts.R")
ORIGINAL <- FIGURE_SCRIPTS_ORIGINAL

cat("\n== ORIGINAL-PIPELINE GENERATORS ==\n")
for (p in names(ORIGINAL)) run_one(p, ORIGINAL[[p]])

# --- 2) revision generators ---------------------------------------------------
# NAMED, NOT GLOBBED. A glob on "^rev_fig" empties the moment that prefix is
# retired, and this runner would still exit 0 having generated nothing --
# the assembler would then copy the previous run's figures and report success.
# The list below was generated from the glob and asserted set-identical to it.
REVISION <- FIGURE_SCRIPTS_REVISION
stopifnot(all(file.exists(REVISION)))
cat(sprintf("\n== REVISION GENERATORS (%d) ==\n", length(REVISION)))
for (p in REVISION) run_one(p, "revision")

# --- 3) assemble --------------------------------------------------------------
cat("\n== ASSEMBLE ==\n")
run_one("code/08_figures/zz_assemble_figures.R", "assemble")

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
