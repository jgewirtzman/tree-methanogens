# ==============================================================================
# flux_windows.R -- saved closure windows and the fitting precision for the goFlux scripts
# ------------------------------------------------------------------------------
# Closure windows are picked by hand with goFlux::click.peak2() and are part of the
# record: every goFlux fitting script loads its saved window file and re-fits from it.
# The picker opens only when no saved file exists, or when FLUX_REPICK=1 is set, and
# only in an interactive session -- a pipeline run never stops at a prompt.
#
# Fitting precision: every fit uses the analyzer precision (sigma) of its field day,
# measured from the closures by 02_flux/02b_instrument_precision.R, with conf = 0.95, so
# that goFlux's MDF column (and with it the HM curvature bound kappa.max and best.flux's
# MDF criterion) is 1.96 sigma / t x flux term, while best.flux's RMSE, MAE and SE
# criteria compare the fit with sigma itself. best.flux() then selects LM or HM.
#
# Where they live: saved windows are inputs, in data/raw/flux_windows/, and no pipeline
# script writes there. A re-pick (FLUX_REPICK=1) writes its windows to
# data/processed/flux/repicked/; copy them into data/raw/flux_windows/ to adopt them.
# ==============================================================================
# Paths from the repository root, valid from the root or a script's own folder.
repo_path <- function(rel) {
  for (up in c("", "../", "../../", "../../../"))
    if (file.exists(paste0(up, "code/pipeline.csv"))) return(paste0(up, rel))
  stop("cannot find the repository root from ", getwd(), call. = FALSE)
}
window_path <- function(name) repo_path(file.path("data/raw/flux_windows", name))

# goFlux must be renv.lock's version: later versions compute the closure time, and so
# the MDF, differently and converge Hutchinson-Mosier fits on other deployments
# (code/lib/project_library.R).
local({
  want <- jsonlite::read_json(repo_path("renv.lock"))$Packages$goFlux$Version
  got  <- if ("goFlux" %in% loadedNamespaces()) getNamespaceVersion("goFlux") else as.character(utils::packageVersion("goFlux"))
  if (!identical(unname(got), want))
    stop("goFlux ", got, " is loaded but renv.lock has ", want, "; run from the repository's ",
         "package library (source code/lib/project_library.R first)", call. = FALSE)
})
repick_path <- function(name) {
  d <- repo_path("data/processed/flux/repicked"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(d, name)
  message("re-picked windows written to ", f, "; copy into data/raw/flux_windows/ to adopt them")
  f
}

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

# Set each closure's CH4 and CO2 precision to its field day's sigma (02b) before goFlux().
SIGMA_CONF <- 0.95
set_day_sigma <- function(manID) {
  S <- read.csv(repo_path("data/processed/flux/analyzer_sigma_by_day.csv"), stringsAsFactors = FALSE)
  k <- match(format(manID$POSIX.time, "%Y-%m-%d", tz = "UTC"), S$day)
  if (anyNA(k[manID$flag == 1]))
    stop("no analyzer sigma for ", length(unique(manID$UniqueID[is.na(k) & manID$flag == 1])),
         " closures; run code/02_flux/02b_instrument_precision.R", call. = FALSE)
  manID$CH4_prec <- S$sigma_CH4_ppb[k]
  manID$CO2_prec <- S$sigma_CO2_ppm[k]
  manID
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
