#!/usr/bin/env Rscript
# ==============================================================================
# rev_rescue_untagged_click.R  (revision "rescue" — INTERACTIVE, run in RStudio)
# Computes CH4 (+CO2) fluxes for the 49 monthly untagged/dead-snag measurements
# that were dropped before flux calc. Follows the standard semirigid workflow:
# import LGR (import2RData) -> obs.win -> click.peak2 (YOU click ~49 peaks) -> goFlux.
# Prereq: run code/revision/rev_rescue_untagged_prep.R first (writes untagged_auxfile.csv).
# Run from repo ROOT, interactively (needs a graphics device for clicking).
# Writes outputs/revision/untagged_fluxes.csv
# ==============================================================================
suppressMessages({library(goFlux);library(dplyr);library(readr)})

# ---- 1. import LGR3 traces (same settings as 04_goflux_trees.R) --------------
lgr_path <- "data/raw/lgr/semirigid_2020-2021"
cat("Importing LGR3 traces from", lgr_path, "...\n")
lgr_data <- import2RData(path = lgr_path, instrument = "UGGA", date.format = "mdy",
                         timezone = "UTC", keep_all = FALSE, prec = c(0.35, 0.9, 200), merge = TRUE)
cat("LGR loaded:", nrow(lgr_data), "rows;",
    "time range", as.character(min(lgr_data$POSIX.time)), "to", as.character(max(lgr_data$POSIX.time)), "\n")

# ---- 2. untagged auxfile (from rev_rescue_untagged_prep.R) -------------------
auxfile <- read.csv("outputs/revision/untagged_auxfile.csv")
auxfile$start.time <- as.POSIXct(auxfile$start.time_formatted, tz = "UTC")
cat("Auxfile:", nrow(auxfile), "untagged measurements\n")
stopifnot(all(c("UniqueID","start.time","Area","Vtot","Tcham","Pcham") %in% names(auxfile)))

# ---- 3. observation windows --------------------------------------------------
ow <- obs.win(inputfile = lgr_data, auxfile = auxfile, gastype = "CH4dry_ppb",
              obs.length = 600, shoulder = 300)
cat("Created", length(ow), "windows. First window has", nrow(ow[[1]]), "points.\n")
if (Sys.info()["sysname"] == "Darwin") options(device = "quartz")

# ---- 4. INTERACTIVE: click start/end of each CH4 peak ------------------------
# Click START then END of the chamber-closure change on each plot (~49 total).
# Do it in batches to avoid fatigue, e.g.:
#   manID   <- click.peak2(ow, gastype="CH4dry_ppb", seq = 1:25,
#                          plot.lim=c(1800,2500), save.plots="outputs/revision/untagged_ID_b1")
#   manID2  <- click.peak2(ow, gastype="CH4dry_ppb", seq = 26:length(ow),
#                          plot.lim=c(1800,2500), save.plots="outputs/revision/untagged_ID_b2")
#   manID   <- rbind(manID, manID2)
cat("\n>>> Now run click.peak2() interactively (see commented lines in this script).\n",
    ">>> When done, you'll have `manID`; then run the flux + save block below.\n")

# ---- 5. fluxes + save (run after clicking) ----------------------------------
compute_and_save <- function(manID) {
  ch4 <- goFlux(manID, "CH4dry_ppb"); co2 <- goFlux(manID, "CO2dry_ppm")
  out <- auxfile %>%
    select(UniqueID, site, species, Sample, Dstem, Area, Vtot) %>%
    left_join(ch4 %>% transmute(UniqueID, CH4_best.flux = best.flux, CH4_model = model,
                                CH4_r2 = if ("HM.r2" %in% names(.)) HM.r2 else NA), by = "UniqueID") %>%
    left_join(co2 %>% transmute(UniqueID, CO2_best.flux = best.flux), by = "UniqueID") %>%
    mutate(dead = grepl("dead|snag", Sample, ignore.case = TRUE))
  write.csv(out, "outputs/revision/untagged_fluxes.csv", row.names = FALSE)
  cat("Wrote outputs/revision/untagged_fluxes.csv (", nrow(out), "trees-measurements)\n")
  print(out %>% group_by(site, dead) %>%
          summarise(n = n(), CH4_median = round(median(CH4_best.flux, na.rm = TRUE), 3), .groups = "drop"))
  invisible(out)
}
# after clicking:  compute_and_save(manID)
