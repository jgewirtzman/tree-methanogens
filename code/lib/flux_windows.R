# ==============================================================================
# flux_windows.R -- saved closure windows and model choices for the goFlux scripts
# ------------------------------------------------------------------------------
# Closure windows are picked by hand with goFlux::click.peak2() and are part of the
# record: every goFlux fitting script loads its saved window file and re-fits from it.
# The picker opens only when no saved file exists, or when FLUX_REPICK=1 is set, and
# only in an interactive session -- a pipeline run never stops at a prompt.
#
# Model choices: goFlux 0.2.0 (renv.lock) converges the Hutchinson-Mosier fit on some
# deployments where the version used when these data were processed (2025) did not,
# and so selects differently for ~9% of deployments. The selection made at the time
# is recorded in data/processed/flux/flux_model_choices.csv and applied here, so a
# re-fit reproduces the stored fluxes. (Built by code/tools/rebuild_closure_windows.R.)
# ==============================================================================
use_saved_windows <- function(path) {
  if (Sys.getenv("FLUX_REPICK") == "1") {
    if (!interactive()) stop("FLUX_REPICK=1 needs an interactive session (the picker is graphical)", call. = FALSE)
    return(FALSE)
  }
  if (file.exists(path)) return(TRUE)
  if (!interactive())
    stop("no saved closure windows at ", path, " and a non-interactive run cannot open the picker; ",
         "see code/tools/README.md", call. = FALSE)
  FALSE
}

load_windows <- function(path) {
  m <- as.data.frame(readr::read_csv(path, show_col_types = FALSE, guess_max = 1e6))
  for (v in intersect(c("POSIX.time", "start.time", "start.time_corr", "end.time_corr"), names(m)))
    m[[v]] <- as.POSIXct(m[[v]], tz = "UTC")
  cat(sprintf("loaded saved closure windows: %s (%d windows, %d rows)\n",
              basename(path), length(unique(m$UniqueID)), nrow(m)))
  m
}

apply_model_choices <- function(best, gas, path) {
  if (!file.exists(path)) { cat("no recorded model choices; using best.flux() as run\n"); return(best) }
  mc <- read.csv(path, stringsAsFactors = FALSE); mc <- mc[mc$gas == gas, ]
  k <- match(best$UniqueID, mc$UniqueID); has <- !is.na(k)
  ch <- mc$model[k[has]]
  changed <- sum(best$model[has] != ch, na.rm = TRUE)
  best$model[has] <- ch
  best$best.flux[has] <- ifelse(ch == "HM", best$HM.flux[has], best$LM.flux[has])
  cat(sprintf("%s: applied %d recorded model choices (%d differ from this goFlux version's choice)\n",
              gas, sum(has), changed))
  best
}

# Diagnostic plots (one PDF page per deployment) are slow and are not a data product:
# made only with FLUX_PLOTS=1, and a plotting failure never stops a run.
maybe_flux_plot <- function(...) {
  if (Sys.getenv("FLUX_PLOTS") != "1") return(list())
  tryCatch(goFlux::flux.plot(...), error = function(e) { message("flux.plot skipped: ", conditionMessage(e)); list() })
}
maybe_flux2pdf <- function(plot.list, ...) {
  if (!length(plot.list)) { cat("diagnostic flux plots skipped (set FLUX_PLOTS=1 to make them)\n"); return(invisible()) }
  tryCatch(goFlux::flux2pdf(plot.list = plot.list, ...), error = function(e) message("flux2pdf skipped: ", conditionMessage(e)))
}
