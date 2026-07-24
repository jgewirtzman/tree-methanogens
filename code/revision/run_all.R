#!/usr/bin/env Rscript
# ==============================================================================
# Revision pipeline — regenerate every final revision output and assemble the
# numbered manuscript figure set. Run from the repo root:
#     Rscript code/revision/run_all.R
# Reads data/ (the Zenodo drop-in archive) + code/ (git); writes to
# outputs/revision/ and outputs/revision/figures/{main,SI}/.
# Order: stats first (a few produce CSVs the figures consume), then figures, then
# the assembler. Exploratory scripts (code/revision/exploratory/) are NOT run.
# ==============================================================================
run <- function(f) {
  cat("\n>>>", basename(f), "\n")
  st <- tryCatch(system2("Rscript", f, stdout = "", stderr = ""), warning = function(w) 1L)
  if (!identical(st, 0L)) cat("   [non-zero exit — check above]\n")
}
# 0) unchanged original-pipeline figures (fig8, several SI) — the assembler copies
#    these through; the revision generators below override the changed ones. Two
#    scripts fail here by design (fig5, S12 methanome) — both replaced by revision
#    versions, so the assembler still finds them.
cat(">>> generate_all_figures.R (original pipeline — populates unchanged figures)\n")
run("code/generate_all_figures.R")

stats <- sort(list.files("code/revision", "^rev_stat_.*\\.R$", full.names = TRUE))
figs  <- sort(list.files("code/revision", "^rev_fig.*\\.R$",   full.names = TRUE))
cat(sprintf("== Revision pipeline: %d stat + %d figure generators ==\n", length(stats), length(figs)))
for (f in c(stats, figs)) run(f)
cat("\n>>> assembling numbered manuscript figures\n")
source("code/revision/rev_00_assemble_figures.R")
cat("\nDone. Numbered figures in outputs/revision/figures/{main,SI}/ (see MANIFEST.md).\n")
