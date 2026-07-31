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

# ---- 1. import LGR3 traces (same as 04_goflux_trees.R: flatten subfolders + unzip) ----
# cached after first import so re-runs (after a clicking mistake) are instant
cache <- "outputs/revision/.lgr_cache.rds"
if (file.exists(cache)) {
  lgr_data <- readRDS(cache); cat("Loaded cached LGR (delete", cache, "to force re-import)\n")
} else {
  data_path <- "data/raw/lgr/semirigid_2020-2021"
  cat("Gathering LGR3 traces from", data_path, "...\n")
  zip_files <- list.files(data_path, recursive = TRUE, pattern = "\\.zip$", full.names = TRUE)
  temp_extract <- tempfile("lgr3_extract_"); dir.create(temp_extract, recursive = TRUE)
  for (z in zip_files) try(unzip(z, exdir = temp_extract, overwrite = TRUE), silent = TRUE)
  txt1 <- list.files(data_path, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE); txt1 <- txt1[file.size(txt1) > 0]
  txt2 <- list.files(temp_extract, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE)
  flat <- tempfile("lgr3_flat_"); dir.create(flat, recursive = TRUE); file.copy(c(txt1, txt2), flat)
  cat("Flattened", length(list.files(flat, pattern = "\\.txt$")), "trace files. Importing...\n")
  lgr_data <- import2RData(path = flat, instrument = "UGGA", date.format = "mdy",
                           timezone = "UTC", keep_all = FALSE, prec = c(0.35, 0.9, 200), merge = TRUE)
  unlink(temp_extract, recursive = TRUE); unlink(flat, recursive = TRUE); saveRDS(lgr_data, cache)
}
cat("LGR loaded:", nrow(lgr_data), "rows;",
    "time range", as.character(min(lgr_data$POSIX.time)), "to", as.character(max(lgr_data$POSIX.time)), "\n")

# ---- 2. untagged auxfile (from rev_rescue_untagged_prep.R) -------------------
auxfile <- read.csv("outputs/revision/untagged_auxfile.csv")
auxfile$start.time <- as.POSIXct(auxfile$start.time_formatted, tz = "UTC")
cat("Auxfile:", nrow(auxfile), "untagged measurements\n")
stopifnot(all(c("UniqueID","start.time","Area","Vtot","Tcham","Pcham") %in% names(auxfile)))

# ---- 3. observation windows (drop the 4 with no LGR trace: 07-18, 02-27, 03-03) ----
ow_all <- obs.win(inputfile = lgr_data, auxfile = auxfile, gastype = "CH4dry_ppb",
                  obs.length = 600, shoulder = 300)
sizes <- sapply(ow_all, nrow); ow <- ow_all[sizes > 30]
cat("Windows with data:", length(ow), "of", length(ow_all),
    "(", sum(sizes <= 30), "empty — LGR trace missing for those dates, unrecoverable)\n")
if (Sys.info()["sysname"] == "Darwin") options(device = "quartz")

# ---- 4. INTERACTIVE: click start/end of each CH4 peak ------------------------
# Click START then END of the chamber-closure change on each plot (~49 total).
# Do it in batches to avoid fatigue, e.g.:
#   manID   <- click.peak2(ow, gastype="CH4dry_ppb", seq = 1:25,
#                          plot.lim=c(1800,2500), save.plots="outputs/revision/untagged_ID_b1")
#   manID2  <- click.peak2(ow, gastype="CH4dry_ppb", seq = 26:length(ow),
#                          plot.lim=c(1800,2500), save.plots="outputs/revision/untagged_ID_b2")
#   manID   <- rbind(manID, manID2)
cat("\n>>> Now click the", length(ow), "windows interactively, in 2 batches, e.g.:\n",
    '    manID  <- click.peak2(ow, gastype="CH4dry_ppb", seq=1:23,            plot.lim=c(1800,6000), save.plots="outputs/revision/untagged_ID_b1")\n',
    '    manID2 <- click.peak2(ow, gastype="CH4dry_ppb", seq=24:length(ow),   plot.lim=c(1800,6000), save.plots="outputs/revision/untagged_ID_b2")\n',
    "    manID  <- rbind(manID, manID2); compute_and_save(manID)\n",
    ">>> Widen plot.lim if a peak is clipped. The 4 dropped measurements (missing LGR) appear as NA in the output.\n")

# ---- 5. fluxes + save (run after clicking) ----------------------------------
compute_and_save <- function(manID) {
  # goFlux 0.2.0: goFlux() gives LM/HM fluxes; best.flux() selects the best model
  ch4 <- best.flux(goFlux(manID, "CH4dry_ppb")); co2 <- best.flux(goFlux(manID, "CO2dry_ppm"))
  out <- auxfile %>%
    select(UniqueID, site, species, Sample, Dstem, Area, Vtot) %>%
    left_join(ch4 %>% transmute(UniqueID, CH4_best.flux = best.flux, CH4_model = model, CH4_quality = quality.check), by = "UniqueID") %>%
    left_join(co2 %>% transmute(UniqueID, CO2_best.flux = best.flux), by = "UniqueID") %>%
    mutate(dead = grepl("dead|snag", Sample, ignore.case = TRUE))
  write.csv(out, "outputs/revision/untagged_fluxes.csv", row.names = FALSE)
  cat("Wrote outputs/revision/untagged_fluxes.csv (", nrow(out), "trees-measurements)\n")
  print(out %>% group_by(site, dead) %>%
          summarise(n = n(), CH4_median = round(median(CH4_best.flux, na.rm = TRUE), 3), .groups = "drop"))
  invisible(out)
}
# after clicking:  compute_and_save(manID)
