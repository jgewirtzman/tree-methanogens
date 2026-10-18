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
#   B  tables and model      canonical tables (data/tables/), then random-forest training tables, fit, export, held-out skill
#   C  drivers and scaling   climatologies, moisture surface, predictions, budget, grid
#   D  analyses              model audits, molecular summaries, statistics and tables
#   E  figures               one generator per paper figure, then the assembler
#   F  archive               data/compiled/: result tables and the README of data/tables/
#   G  gate                  check_consistency.R
#
# Usage (from the repository root)
#   Rscript code/run_all.R                    # stages B-G: from data/processed/ and data/raw/
#   Rscript code/run_all.R --from raw         # stages A-G: everything from data/raw/
#   Rscript code/run_all.R --from C           # skip the model fit (uses outputs/models/)
#   Rscript code/run_all.R --only E           # one stage (make_figures.R does this)
#
# Full and lite runs
#   Rscript code/run_all.R --lite --from raw  # every script, in minutes, outside the repository
#   The run sizes (forest sizes, cross-validation repeats, resampling loops) are in
#   code/lib/run_mode.R. The default is the full run, the canonical one. --lite builds a
#   fresh copy of the repository at --lite-root (default: <temp>/tree-methanogens-lite,
#   replaced on every lite run): code/, renv settings and the folder layout are copied,
#   data/raw/ is linked (read only by the pipeline; checked unchanged afterwards) and, for
#   a run that starts after stage A, the existing data/processed/ (and, after B,
#   data/tables/, data/compiled/ and outputs/) are copied in; a run from raw gets only the
#   archived PICRUSt2 pathway tables (data/processed/molecular/picrust/). The pipeline then runs inside
#   the copy with PIPELINE_RUN_MODE=lite, so every write lands in the copy and the
#   canonical outputs are never touched. Lite numbers are not results.
#
# Each script's console output goes to outputs/logs/<script>.txt; a summary of the
# run is written to outputs/logs/pipeline_run.csv. Scripts whose manifest row says
# workdir = script are run from their own folder (the older processing scripts use
# relative paths); all others run from the repository root.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
LITE <- "--lite" %in% args
STAGES <- c("A", "B", "C", "D", "E", "F", "G")
from <- toupper(opt("--from", "B")); if (from == "RAW") from <- "A"; if (from == "PROCESSED") from <- "B"
to   <- toupper(opt("--to", "G"))
only <- opt("--only", NA)
run_stages <- if (!is.na(only)) toupper(only) else STAGES[match(from, STAGES):match(to, STAGES)]
stopifnot(all(run_stages %in% STAGES))

if (!file.exists("code/pipeline.csv"))
  stop("Run from the repository root:  Rscript code/run_all.R", call. = FALSE)
# Several scripts read column names and text with non-ASCII characters (µ, °), which
# fail under a C locale. The scripts inherit this process's locale.
if (!isTRUE(l10n_info()[["UTF-8"]]))
  stop("run under a UTF-8 locale, e.g.  LC_ALL=en_US.UTF-8 Rscript code/run_all.R", call. = FALSE)
source("code/lib/run_mode.R")   # RUN_MODE; stops on lite settings outside a lite copy

# ---- --lite: build the copy and run this script inside it ---------------------------
if (LITE) {
  if (RUN_MODE == "lite") stop("--lite inside a lite copy", call. = FALSE)
  ROOT <- normalizePath(".")
  LR <- normalizePath(opt("--lite-root", file.path(dirname(tempdir()), "tree-methanogens-lite")), mustWork = FALSE)
  inside <- function(a, b) startsWith(paste0(a, "/"), paste0(b, "/"))
  if (inside(LR, ROOT) || inside(ROOT, LR)) stop("--lite-root must lie outside the repository: ", LR, call. = FALSE)
  if (file.exists(LR)) {            # only a previous lite copy is ever replaced
    if (!file.exists(file.path(LR, ".lite_run"))) stop(LR, " exists and is not a lite copy; not touching it", call. = FALSE)
    if (nzchar(Sys.readlink(file.path(LR, "data/raw")))) file.remove(file.path(LR, "data/raw"))   # the link, never its target
    unlink(LR, recursive = TRUE)
  }
  dir.create(LR, recursive = TRUE)
  rel <- function(d) substring(list.dirs(file.path(ROOT, d)), nchar(ROOT) + 2)
  # folder layout of data/ and outputs/ (not data/raw/, which is linked, or the archives)
  lay <- c(rel("data"), rel("outputs"))
  lay <- lay[!grepl("^data/(raw|deprecated)(/|$)|^outputs/(snapshots|scratch)(/|$)", lay)]
  for (d in lay) dir.create(file.path(LR, d), recursive = TRUE, showWarnings = FALSE)
  # the code (not goFlux import intermediates, Finder files or the revision notes link)
  cf <- list.files(file.path(ROOT, "code"), recursive = TRUE, all.files = TRUE)
  cf <- file.path("code", cf[!grepl("(^|/)RData/|(^|/)\\.DS_Store$|^revision/", cf)])
  copy_in <- function(f) {
    for (d in unique(dirname(f))) dir.create(file.path(LR, d), recursive = TRUE, showWarnings = FALSE)
    stopifnot(all(file.copy(file.path(ROOT, f), file.path(LR, f), copy.date = TRUE)))
  }
  copy_in(c(cf, ".Rprofile", "renv.lock", "renv/activate.R", "renv/settings.json"))
  file.symlink(file.path(ROOT, "renv/library"), file.path(LR, "renv/library"))
  file.symlink(file.path(ROOT, "data/raw"), file.path(LR, "data/raw"))
  # a run from raw needs, of data/processed/, only the archived PICRUSt2 pathway tables, and
  # those only when data/raw/picrust/ is absent (01_import/06_picrust_pathways.R re-derives
  # them otherwise; either way it checks them against their recorded checksums)
  if (run_stages[1] == "A") {
    pp <- "data/processed/molecular/picrust"
    copy_in(file.path(pp, list.files(file.path(ROOT, pp), "[.]csv$")))
  }
  # a run that starts after a stage reads what that stage left in the repository
  if (run_stages[1] != "A") {
    keep <- "data/processed"
    if (!run_stages[1] %in% c("A", "B")) keep <- c(keep, "data/tables", "data/compiled", "outputs")
    f <- unlist(lapply(keep, function(d) file.path(d, list.files(file.path(ROOT, d), recursive = TRUE))))
    copy_in(f[!grepl("^outputs/(snapshots|scratch|logs)/", f)])
  }
  writeLines(c("lite copy made by code/run_all.R --lite", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ROOT),
             file.path(LR, ".lite_run"))
  pass <- args[args != "--lite"]
  i <- match("--lite-root", pass); if (!is.na(i)) pass <- pass[-c(i, i + 1)]
  cat(sprintf("Lite run in %s (data/raw linked from the repository)\n", LR))
  t0 <- Sys.time(); owd <- setwd(LR)
  st <- system2("Rscript", c("code/run_all.R", pass),
                env = c("PIPELINE_RUN_MODE=lite", paste0("RENV_PROJECT=", shQuote(LR))))
  setwd(owd)
  raw <- list.files(file.path(ROOT, "data/raw"), recursive = TRUE, full.names = TRUE)
  touched <- raw[file.mtime(raw) >= t0]
  if (length(touched)) { cat("data/raw/ files changed during the lite run:\n", paste(" ", touched, collapse = "\n"), "\n"); st <- 1L }
  cat(sprintf("\nLite run %s in %.1f min; its outputs are in %s\n", if (identical(st, 0L)) "finished" else "FAILED",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), LR))
  quit(status = if (identical(st, 0L)) 0L else 1L)
}
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

# renv checks the library against renv.lock at every R start-up (~4 s). It has run once,
# for this process; the scripts skip it (time only: the library they load is the same).
Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = "FALSE")

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

cat(sprintf("Pipeline (%s run): stages %s, %d scripts\n", RUN_MODE, paste(run_stages, collapse = ""), nrow(P)))
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
if (length(fails)) {
  cat("failed (non-fatal):\n", paste("  ", fails, collapse = "\n"), "\n")
  quit(status = 1)
}
