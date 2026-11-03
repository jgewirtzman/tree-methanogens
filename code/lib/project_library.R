# ==============================================================================
# project_library.R -- the project's package library (renv.lock), however a script is started
# ------------------------------------------------------------------------------
# The repository's .Rprofile activates renv only for an R session started at the
# repository root. A script that runs from its own folder (manifest workdir = script)
# inherits the project library when run_all.R starts it, but started by hand from that
# folder (cd code/02_flux/static; Rscript 03_goflux_trees_2021.R) it used to load the
# user library instead, whose goFlux can differ from renv.lock's (0.5.0.9003, the
# jgewirtzman fork): goFlux versions differ in how they take the closure time of a
# record and in which model best.flux selects, so a fit must run on the pinned version.
#
# Sourced first, before any library(), by every workdir = script script: activates the
# project library when it is not the active one (nothing happens under run_all.R).
# code/lib/flux_windows.R also checks that goFlux is renv.lock's version.
# ==============================================================================
local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "code/pipeline.csv")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "renv/activate.R"))) stop("cannot find the repository root from ", getwd(), call. = FALSE)
  active <- Sys.getenv("RENV_PROJECT")
  if (!nzchar(active) || normalizePath(active, mustWork = FALSE) != d) {
    owd <- setwd(d); on.exit(setwd(owd))
    Sys.setenv(RENV_PROJECT = d)
    source("renv/activate.R")
  }
})
