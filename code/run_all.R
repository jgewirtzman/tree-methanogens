# ==============================================================================
# run_all.R -- the whole pipeline, in order, from code/pipeline.csv
# ------------------------------------------------------------------------------
# The manifest is the single list of what runs and in what order. Each row is one
# script with its stage, working directory, whether a failure stops the run, and
# what it is for. This file only reads the manifest and runs it; to add, remove or
# reorder a step, edit code/pipeline.csv.
#
# Stages
#   A  raw -> processed      import, flux fitting (goFlux), cleaning, harmonisation
#   B  model                 random-forest training tables, fit, export, held-out skill
#   C  drivers and scaling   climatologies, moisture surface, predictions, budget, grid
#   D  analyses              model audits, molecular summaries, statistics and tables
#   E  figures               one generator per paper figure, then the assembler
#   F  archive               data/compiled/ datasets and their README
#   G  gate                  check_consistency.R
#
# Usage (from the repository root)
#   Rscript code/run_all.R                    # stages B-G: everything from data/processed/
#   Rscript code/run_all.R --from raw         # stages A-G: everything from data/raw/
#   Rscript code/run_all.R --from C           # skip the model fit (uses outputs/models/)
#   Rscript code/run_all.R --only E           # one stage (make_figures.R does this)
#
# Each script's console output goes to outputs/logs/<script>.txt; a summary of the
# run is written to outputs/logs/pipeline_run.csv. Scripts whose manifest row says
# workdir = script are run from their own folder (the older processing scripts use
# relative paths); all others run from the repository root.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
STAGES <- c("A", "B", "C", "D", "E", "F", "G")
from <- toupper(opt("--from", "B")); if (from == "RAW") from <- "A"; if (from == "PROCESSED") from <- "B"
to   <- toupper(opt("--to", "G"))
only <- opt("--only", NA)
run_stages <- if (!is.na(only)) toupper(only) else STAGES[match(from, STAGES):match(to, STAGES)]
stopifnot(all(run_stages %in% STAGES))

if (!file.exists("code/pipeline.csv"))
  stop("Run from the repository root:  Rscript code/run_all.R", call. = FALSE)
P <- read.csv("code/pipeline.csv", stringsAsFactors = FALSE)
stopifnot(all(P$stage %in% STAGES), all(file.exists(P$script)))
P <- P[P$stage %in% run_stages, ]

# Starting after stage B needs a fitted model on disk.
if (!"B" %in% run_stages && any(c("C", "D", "E") %in% run_stages)) {
  need <- c("outputs/models/RF_MODELS.RData", "outputs/models/TRAINING_DATA.RData")
  if (!all(file.exists(need)))
    stop("no fitted model in outputs/models/; run stage B first (Rscript code/run_all.R --from B)",
         call. = FALSE)
}

ROOT <- normalizePath(".")
LOGDIR <- file.path(ROOT, "outputs/logs")
dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)
# Marker for the figure assembler's staleness check: a figure older than this was
# not produced by this run.
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), file.path(ROOT, "outputs/.pipeline_run_started"))

run_one <- function(script, workdir) {
  logf <- file.path(LOGDIR, sub("\\.R$", ".txt", basename(script)))
  dir  <- if (workdir == "script") file.path(ROOT, dirname(script)) else ROOT
  file <- if (workdir == "script") basename(script) else script
  t0 <- Sys.time(); owd <- setwd(dir); on.exit(setwd(owd))
  st <- tryCatch(system2("Rscript", file, stdout = logf, stderr = logf),
                 warning = function(w) 1L, error = function(e) 1L)
  list(status = if (identical(st, 0L)) "ok" else "FAIL",
       secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

cat(sprintf("Pipeline: stages %s, %d scripts\n", paste(run_stages, collapse = ""), nrow(P)))
res <- data.frame(stage = P$stage, script = P$script, status = NA_character_, secs = NA_real_)
for (i in seq_len(nrow(P))) {
  if (i == 1 || P$stage[i] != P$stage[i - 1]) cat(sprintf("\n== stage %s ==\n", P$stage[i]))
  r <- run_one(P$script[i], P$workdir[i])
  res$status[i] <- r$status; res$secs[i] <- r$secs
  cat(sprintf("  %-55s %s (%ss)\n", basename(P$script[i]), r$status, r$secs))
  if (r$status == "FAIL" && isTRUE(as.logical(P$fatal[i]))) {
    write.csv(res, file.path(LOGDIR, "pipeline_run.csv"), row.names = FALSE)
    stop(sprintf("required step failed: %s -- see outputs/logs/%s", P$script[i],
                 sub("\\.R$", ".txt", basename(P$script[i]))), call. = FALSE)
  }
}
write.csv(res, file.path(LOGDIR, "pipeline_run.csv"), row.names = FALSE)
fails <- res$script[res$status == "FAIL"]
cat(sprintf("\n%d scripts, %d failed, %.0f min\n", nrow(res), length(fails), sum(res$secs) / 60))
if (length(fails)) cat("failed (non-fatal):\n", paste("  ", fails, collapse = "\n"), "\n")
