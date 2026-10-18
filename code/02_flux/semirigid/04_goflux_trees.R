# ==============================================================================
# 04_goflux_trees.R -- monthly semi-rigid stem CO2 and CH4 fluxes (goFlux), 2020-21
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/02_flux/semirigid/).
#
# Fits goFlux (linear and Hutchinson-Mosier) to each deployment's hand-picked closure
# window and selects the model with best.flux(); the model selected when the data were
# processed is then applied from flux_model_choices.csv (code/lib/flux_windows.R).
# Chamber temperature and geometry (Tcham, Area, Vtot) come from the current auxfile
# (station air in true time; the chamber used that day, code/lib/stem_chamber_geometry.R),
# not from the values saved with the windows.
#
# Closure windows: data/raw/flux_windows/lgr_manual_identification_semirigid_tree.csv.
# Only when it is absent, or with FLUX_REPICK=1 in an interactive session, are the raw
# LGR3 records read and the picker (goFlux::click.peak2) opened; the new windows go to
# data/processed/flux/repicked/ (code/lib/flux_windows.R, code/tools/README.md).
# Diagnostic plots only with FLUX_PLOTS=1. The campaign table is assembled by
# code/02_flux/helper_assemble_campaign_flux.R.
#
# Inputs:  data/raw/flux_windows/lgr_manual_identification_semirigid_tree.csv
#          data/raw/flux_windows/flux_model_choices.csv
#          data/processed/flux/auxfile_goFlux_with_weather.csv (03_prep_tree_auxfile.R)
#          [re-pick only] data/raw/lgr/semirigid_2020-2021/ (raw LGR3 records)
# Outputs (data/processed/flux/):
#          {CO2,CH4}_flux_lgr_results_semirigid_tree.csv       goFlux results
#          {CO2,CH4}_best_flux_lgr_results_semirigid_tree.csv  best.flux results
#          semirigid_tree_final_complete_dataset.csv  (assemble_campaign("semirigid_tree"))
# ==============================================================================

source("../../lib/project_library.R")   # the project package library (renv.lock), however this script is started
library(goFlux)
library(dplyr)
library(readr)

source("../../lib/flux_windows.R")
# Date-times are written with microseconds, the setting goFlux's raw-data import
# (import.UGGA) leaves behind; set here so the tables do not depend on whether the raw
# records were read in this session.
options(digits.secs = 6)

data_path <- "../../../data/raw/lgr/semirigid_2020-2021"
MC_PATH <- window_path("flux_model_choices.csv")
WIN <- window_path("lgr_manual_identification_semirigid_tree.csv")

lgr3_auxfile <- read_csv("../../../data/processed/flux/auxfile_goFlux_with_weather.csv", show_col_types = FALSE)
cat("Stem auxfile:", nrow(lgr3_auxfile), "deployments\n")

# ---- closure windows ------------------------------------------------------------------
if (use_saved_windows(WIN)) {
  manID.lgr3 <- load_windows(WIN)
  ai <- match(manID.lgr3$UniqueID, lgr3_auxfile$UniqueID)
  manID.lgr3$Tcham <- lgr3_auxfile$Tcham[ai]
  manID.lgr3$Area <- lgr3_auxfile$Area[ai]
  manID.lgr3$Vtot <- lgr3_auxfile$Vtot[ai]
  stopifnot(!anyNA(manID.lgr3$Tcham), !anyNA(manID.lgr3$Area), !anyNA(manID.lgr3$Vtot))
} else {
  # interactive re-pick: raw LGR3 records (loose .txt files and .zip archives) ->
  # observation windows -> click.peak2 in batches of 20
  temp_extract_dir <- tempfile("lgr3_extract_"); dir.create(temp_extract_dir, recursive = TRUE)
  for (zip_file in list.files(data_path, recursive = TRUE, pattern = "\\.zip$", full.names = TRUE))
    unzip(zip_file, exdir = temp_extract_dir, overwrite = TRUE)
  txt_files <- list.files(data_path, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE)
  txt_files <- c(txt_files[file.size(txt_files) > 0],
                 list.files(temp_extract_dir, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE))
  complete_lgr3_dir <- tempfile("lgr3_complete_"); dir.create(complete_lgr3_dir, recursive = TRUE)
  file.copy(txt_files, complete_lgr3_dir)
  unlink("RData", recursive = TRUE)  # import2RData(merge = TRUE) merges every file already in ./RData
  lgr3_data_complete <- import2RData(path = complete_lgr3_dir, instrument = "UGGA", date.format = "mdy",
                                     timezone = "UTC", keep_all = FALSE, prec = c(0.35, 0.9, 200), merge = TRUE)
  unlink(c(temp_extract_dir, complete_lgr3_dir), recursive = TRUE)
  tz <- attr(lgr3_data_complete$POSIX.time, "tzone"); if (is.null(tz) || tz == "") tz <- "UTC"
  ow.lgr3.complete <- obs.win(inputfile = lgr3_data_complete,
                              auxfile = mutate(lgr3_auxfile, start.time = as.POSIXct(start.time, tz = tz)),
                              gastype = "CO2dry_ppm", obs.length = 600, shoulder = 300)
  names(ow.lgr3.complete) <- sapply(ow.lgr3.complete, function(x) unique(x$UniqueID)[1])
  good_windows <- which(sapply(ow.lgr3.complete, nrow) >= 30)
  n_batches <- ceiling(length(good_windows) / 20)
  manID_batches <- list()
  for (batch_num in seq_len(n_batches)) {
    batch_windows <- good_windows[((batch_num - 1) * 20 + 1):min(batch_num * 20, length(good_windows))]
    cat("\nBatch", batch_num, "of", n_batches, "\n")
    manID_batches[[batch_num]] <- click.peak2(ow.lgr3.complete, seq = batch_windows, gastype = "CO2dry_ppm",
                                              sleep = 3, plot.lim = c(200, 5000),
                                              save.plots = paste0("../../../data/processed/flux/lgr3_batch_", batch_num, "_plots"))
    while (dev.cur() > 1) dev.off()
    if (batch_num < n_batches) { cat("Press Enter for the next batch...\n"); readline() }
  }
  manID.lgr3 <- do.call(rbind, manID_batches)
  write_csv(manID.lgr3, repick_path(basename(WIN)))
}

# ---- fluxes ---------------------------------------------------------------------------
CO2_flux_lgr3 <- goFlux(dataframe = manID.lgr3, gastype = "CO2dry_ppm", H2O_col = "H2O_ppm", warn.length = 60)
CH4_flux_lgr3 <- goFlux(dataframe = manID.lgr3, gastype = "CH4dry_ppb", H2O_col = "H2O_ppm", warn.length = 60)
write_csv(CO2_flux_lgr3, "../../../data/processed/flux/CO2_flux_lgr_results_semirigid_tree.csv")
write_csv(CH4_flux_lgr3, "../../../data/processed/flux/CH4_flux_lgr_results_semirigid_tree.csv")

best_flux <- function(flux_result)
  best.flux(flux.result = flux_result,
            criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", "MDF", "nb.obs", "intercept", "p-value"),
            intercept.lim = NULL, g.limit = 2, p.val = 0.05, k.ratio = 1, warn.length = 60)
CO2_best_lgr3 <- apply_model_choices(best_flux(CO2_flux_lgr3), "CO2", MC_PATH)
CH4_best_lgr3 <- apply_model_choices(best_flux(CH4_flux_lgr3), "CH4", MC_PATH)
write_csv(CO2_best_lgr3, "../../../data/processed/flux/CO2_best_flux_lgr_results_semirigid_tree.csv")
write_csv(CH4_best_lgr3, "../../../data/processed/flux/CH4_best_flux_lgr_results_semirigid_tree.csv")
cat("\nCO2 quality: OK", sum(CO2_best_lgr3$quality.check == "OK", na.rm = TRUE),
    "| flagged", sum(CO2_best_lgr3$quality.check != "OK", na.rm = TRUE), "\n")

# ---- diagnostic plots (FLUX_PLOTS=1 only) ---------------------------------------------
plot_args <- list(shoulder = 30, plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
                  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
                  quality.check = TRUE, best.model = TRUE, p.val.disp = "round")
maybe_flux2pdf(
  plot.list = c(
    do.call(maybe_flux_plot, c(list(flux.results = CO2_best_lgr3, dataframe = manID.lgr3, gastype = "CO2dry_ppm"), plot_args)),
    do.call(maybe_flux_plot, c(list(flux.results = CH4_best_lgr3, dataframe = manID.lgr3, gastype = "CH4dry_ppb"), plot_args))),
  outfile = "../../../data/processed/flux/LGR3_flux_plots_complete.pdf", width = 11.6, height = 8.2)

# ---- campaign table -------------------------------------------------------------------
source("../helper_assemble_campaign_flux.R")
assemble_campaign("semirigid_tree")
