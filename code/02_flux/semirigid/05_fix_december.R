# ==============================================================================
# 05_fix_december.R -- December 2020 soil fluxes, from their own closure windows
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/02_flux/semirigid/). Runs after
# 04_goflux_soils.R.
#
# The December 2020 soil deployments have their own window file. This script fits
# them (goFlux + best.flux + the recorded model choices, as 04_goflux_soils.R) and, in
# every soil result table, replaces the rows of 8 December 2020 (UniqueID containing
# "20201208") and any earlier copy of these deployments with the new fits, so it can run
# more than once without duplicating December.
# Geometry: Area = SOIL_COLLAR_AREA_CM2, Vtot = SOIL_VTOT_L (code/lib/chamber_constants.R);
# the December windows were picked with the pre-correction collar volume, which is not
# used. Tcham is the soil auxfile's (station air in true time).
#
# Closure windows: data/raw/flux_windows/lgr_manual_identification_results_december_soil.csv.
# Only when it is absent, or with FLUX_REPICK=1 in an interactive session, are the raw
# LGR3 records read and the picker opened (code/lib/flux_windows.R); the re-picked
# December windows, and the soil window file with its December windows replaced by
# them, go to data/processed/flux/repicked/. Diagnostic plots only with FLUX_PLOTS=1.
#
# Inputs:  data/raw/flux_windows/lgr_manual_identification_results_december_soil.csv
#          data/raw/flux_windows/flux_model_choices.csv
#          data/processed/flux/auxfile_goFlux_soilflux_with_weather.csv
#          data/processed/flux/{CO2,CH4}_flux_lgr_results_soil.csv,
#          {CO2,CH4}_best_flux_lgr_results_soil.csv,
#          semirigid_tree_final_complete_dataset_soil.csv (04_goflux_soils.R)
#          [re-pick only] data/raw/lgr/semirigid_2020-2021/, data/raw/flux_windows/
#          lgr_manual_identification_results_soil.csv
# Outputs (data/processed/flux/):
#          the five soil tables above, rewritten with the December rows replaced
#          {CO2,CH4}_best_flux_december_lgr_results_soil.csv (December only)
# ==============================================================================

source("../../lib/project_library.R")   # the project package library (renv.lock), however this script is started
library(goFlux)
library(dplyr)
library(readr)
library(lubridate)

source("../../lib/flux_windows.R")
source("../../lib/chamber_constants.R")   # SOIL_VTOT_L, SOIL_COLLAR_AREA_CM2
# Date-times are written with microseconds, the setting goFlux's raw-data import
# (import.UGGA) leaves behind; set here so the tables do not depend on whether the raw
# records were read in this session.
options(digits.secs = 6)

FLUX <- "../../../data/processed/flux"
data_path <- "../../../data/raw/lgr/semirigid_2020-2021"
MC_PATH <- window_path("flux_model_choices.csv")
WIN <- window_path("lgr_manual_identification_results_december_soil.csv")
SOIL_WIN <- window_path("lgr_manual_identification_results_soil.csv")

lgr3_auxfile_full <- read_csv(file.path(FLUX, "auxfile_goFlux_soilflux_with_weather.csv"), show_col_types = FALSE) %>%
  mutate(start.time = as.POSIXct(start.time, tz = "UTC"),
         month = month(start.time))
lgr3_auxfile_december <- lgr3_auxfile_full %>% filter(month == 12)
cat("Soil deployments:", nrow(lgr3_auxfile_full), "| in December:", nrow(lgr3_auxfile_december), "\n")

# ---- closure windows ------------------------------------------------------------------
SAVED_WINDOWS <- use_saved_windows(WIN)
if (SAVED_WINDOWS) {
  manID.lgr3.december <- load_windows(WIN) %>%
    mutate(Tcham = lgr3_auxfile_full$Tcham[match(UniqueID, lgr3_auxfile_full$UniqueID)])
  stopifnot(!anyNA(manID.lgr3.december$Tcham))
} else {
  # interactive re-pick: raw LGR3 records -> December observation windows -> click.peak2
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
  ow.lgr3.december <- obs.win(inputfile = lgr3_data_complete, auxfile = lgr3_auxfile_december,
                              gastype = "CO2dry_ppm", obs.length = 300, shoulder = 300)
  names(ow.lgr3.december) <- sapply(ow.lgr3.december, function(x) unique(x$UniqueID)[1])
  good_windows <- which(sapply(ow.lgr3.december, nrow) >= 30)
  n_batches <- ceiling(length(good_windows) / 20)
  manID_batches <- list()
  for (batch_num in seq_len(n_batches)) {
    batch_windows <- good_windows[((batch_num - 1) * 20 + 1):min(batch_num * 20, length(good_windows))]
    cat("\nDecember batch", batch_num, "of", n_batches, "\n")
    manID_batches[[batch_num]] <- click.peak2(ow.lgr3.december, seq = batch_windows, gastype = "CO2dry_ppm",
                                              sleep = 3, plot.lim = c(200, 5000),
                                              save.plots = paste0(FLUX, "/lgr3_december_batch_", batch_num, "_plots"))
    while (dev.cur() > 1) dev.off()
    if (batch_num < n_batches) { cat("Press Enter for the next batch...\n"); readline() }
  }
  manID.lgr3.december <- do.call(rbind, manID_batches)
  write_csv(manID.lgr3.december, repick_path(basename(WIN)))
}
# one collar geometry for every fit (see header)
manID.lgr3.december <- manID.lgr3.december %>% mutate(Vtot = SOIL_VTOT_L, Area = SOIL_COLLAR_AREA_CM2)

# ---- fluxes ---------------------------------------------------------------------------
CO2_flux_lgr3_december <- goFlux(dataframe = manID.lgr3.december, gastype = "CO2dry_ppm", H2O_col = "H2O_ppm", warn.length = 60)
CH4_flux_lgr3_december <- goFlux(dataframe = manID.lgr3.december, gastype = "CH4dry_ppb", H2O_col = "H2O_ppm", warn.length = 60)

best_flux <- function(flux_result)
  best.flux(flux.result = flux_result,
            criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", "MDF", "nb.obs", "intercept", "p-value"),
            intercept.lim = NULL, g.limit = 2, p.val = 0.05, k.ratio = 1, warn.length = 60)
CO2_best_lgr3_december <- apply_model_choices(best_flux(CO2_flux_lgr3_december), "CO2", MC_PATH)
CH4_best_lgr3_december <- apply_model_choices(best_flux(CH4_flux_lgr3_december), "CH4", MC_PATH)
cat("December CO2 quality: OK", sum(CO2_best_lgr3_december$quality.check == "OK", na.rm = TRUE),
    "| flagged", sum(CO2_best_lgr3_december$quality.check != "OK", na.rm = TRUE), "\n")

plot_args <- list(shoulder = 30, plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
                  plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
                  quality.check = TRUE, best.model = TRUE, p.val.disp = "round")
maybe_flux2pdf(
  plot.list = c(
    do.call(maybe_flux_plot, c(list(flux.results = CO2_best_lgr3_december, dataframe = manID.lgr3.december, gastype = "CO2dry_ppm"), plot_args)),
    do.call(maybe_flux_plot, c(list(flux.results = CH4_best_lgr3_december, dataframe = manID.lgr3.december, gastype = "CH4dry_ppb"), plot_args))),
  outfile = file.path(FLUX, "LGR3_flux_plots_december_soil.pdf"), width = 11.6, height = 8.2)

# ---- replace the December rows in every soil table ----------------------------------------
december_8_pattern <- "20201208"
december_unique_ids <- unique(manID.lgr3.december$UniqueID)
drop_december <- function(d) d %>% filter(!grepl(december_8_pattern, UniqueID), !UniqueID %in% december_unique_ids)

if (!SAVED_WINDOWS) {
  # the soil window file with its December windows replaced by the re-picked ones (the
  # other months are the saved soil windows, unchanged)
  existing_manID <- read_csv(SOIL_WIN, show_col_types = FALSE) %>% mutate(DATE = as.character(DATE))
  updated_manID <- bind_rows(drop_december(existing_manID),
                             manID.lgr3.december %>% mutate(DATE = as.character(DATE)))
  soil_out <- file.path(repo_path("data/processed/flux/repicked"), basename(SOIL_WIN))
  write_csv(updated_manID, soil_out)
  message("soil windows with the re-picked December windows written to ", soil_out,
          " (other months: the saved soil windows); copy into data/raw/flux_windows/ to adopt them")
}

replace_december <- function(file, december) {
  existing <- read_csv(file.path(FLUX, file), show_col_types = FALSE)
  updated <- bind_rows(drop_december(existing), december)
  write_csv(updated, file.path(FLUX, file))
  cat(sprintf("  %-45s %d rows (%d removed, %d added)\n", file, nrow(updated),
              nrow(existing) - nrow(drop_december(existing)), nrow(december)))
}
cat("Replacing the December rows:\n")
replace_december("CO2_flux_lgr_results_soil.csv", CO2_flux_lgr3_december)
replace_december("CO2_best_flux_lgr_results_soil.csv", CO2_best_lgr3_december)
replace_december("CH4_flux_lgr_results_soil.csv", CH4_flux_lgr3_december)
replace_december("CH4_best_flux_lgr_results_soil.csv", CH4_best_lgr3_december)

december_final_data <- lgr3_auxfile_december %>%
  mutate(Vtot = SOIL_VTOT_L, Area = SOIL_COLLAR_AREA_CM2) %>%
  left_join(CO2_best_lgr3_december %>% rename_with(~ paste0("CO2_", .), -UniqueID), by = "UniqueID") %>%
  left_join(CH4_best_lgr3_december %>% rename_with(~ paste0("CH4_", .), -UniqueID), by = "UniqueID")
replace_december("semirigid_tree_final_complete_dataset_soil.csv", december_final_data)

write_csv(CO2_best_lgr3_december, file.path(FLUX, "CO2_best_flux_december_lgr_results_soil.csv"))
write_csv(CH4_best_lgr3_december, file.path(FLUX, "CH4_best_flux_december_lgr_results_soil.csv"))
cat("Wrote {CO2,CH4}_best_flux_december_lgr_results_soil.csv (", nrow(december_final_data), "December deployments )\n")
