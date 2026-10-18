# ==============================================================================
# 01_fix_soil_flux.R -- soil CH4 fluxes re-fitted from the saved soil windows
# ------------------------------------------------------------------------------
# Stage A; workdir: this script's folder (code/03_merge/). Runs after
# 02_flux/semirigid/05_fix_december.R and before the screens and detection limits that
# read the soil fluxes.
#
# Re-fits CH4 (goFlux, then best.flux) for every soil deployment from the saved soil
# closure windows, with the one collar geometry (Area = SOIL_COLLAR_AREA_CM2,
# Vtot = SOIL_VTOT_L; code/lib/chamber_constants.R) and the auxfile's chamber
# temperature (station air in true time), and rewrites the soil CH4 result files.
# The CO2 results of 04_goflux_soils.R / 05_fix_december.R are kept. The model (linear
# or non-linear) of each deployment is the recorded choice
# (data/raw/flux_windows/flux_model_choices.csv), as in every other fit.
# Diagnostic plots only with FLUX_PLOTS=1.
#
# Inputs:  data/raw/flux_windows/lgr_manual_identification_results_soil.csv
#          data/processed/flux/auxfile_goFlux_soilflux_with_weather.csv (Tcham)
#          data/processed/flux/auxfile_goFlux_soilflux_with_weather_formatted.csv (table base)
#          data/processed/flux/CO2_best_flux_lgr_results_soil.csv
#          data/raw/flux_windows/flux_model_choices.csv
# Outputs (data/processed/flux/):
#          CH4_flux_lgr_results_soil.csv, CH4_best_flux_lgr_results_soil.csv (read by the
#          screens, detection and the soil flux table)
#          semirigid_tree_final_complete_dataset_soil_CORRECTED.csv  auxfile + CO2_/CH4_ columns
# ==============================================================================

source("../lib/project_library.R")   # the project package library (renv.lock), however this script is started
library(goFlux)
library(dplyr)
library(readr)

source("../lib/flux_windows.R")       # window_path(), maybe_flux_plot(), maybe_flux2pdf()
source("../lib/chamber_constants.R")  # SOIL_VTOT_L, SOIL_COLLAR_AREA_CM2

FLUX <- "../../data/processed/flux"
cat("Soil collar geometry: Area", SOIL_COLLAR_AREA_CM2, "cm2 | Vtot", round(SOIL_VTOT_L, 3), "L\n")

soil_aux <- read_csv(file.path(FLUX, "auxfile_goFlux_soilflux_with_weather.csv"), show_col_types = FALSE)
manID.lgr3 <- read_csv(window_path("lgr_manual_identification_results_soil.csv"), show_col_types = FALSE) %>%
  mutate(
    Vtot = SOIL_VTOT_L,
    Area = SOIL_COLLAR_AREA_CM2,
    Tcham = soil_aux$Tcham[match(UniqueID, soil_aux$UniqueID)]
  )
stopifnot(!anyNA(manID.lgr3$Tcham))
cat("Saved soil windows:", length(unique(manID.lgr3$UniqueID)), "deployments | CH4 range",
    round(min(manID.lgr3$CH4dry_ppb, na.rm = TRUE), 1), "-", round(max(manID.lgr3$CH4dry_ppb, na.rm = TRUE), 1), "ppb\n")

CH4_flux_lgr3 <- goFlux(dataframe = manID.lgr3, gastype = "CH4dry_ppb", H2O_col = "H2O_ppm", warn.length = 60)
write_csv(CH4_flux_lgr3, file.path(FLUX, "CH4_flux_lgr_results_soil.csv"))

CH4_best_lgr3 <- apply_model_choices(best.flux(
  flux.result = CH4_flux_lgr3,
  criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa", "MDF", "nb.obs", "intercept", "p-value"),
  intercept.lim = NULL, g.limit = 2, p.val = 0.05, k.ratio = 1, warn.length = 60),
  "CH4", window_path("flux_model_choices.csv"))
write_csv(CH4_best_lgr3, file.path(FLUX, "CH4_best_flux_lgr_results_soil.csv"))
cat(sprintf("CH4 best flux: n %d | mean %.4f | median %.4f | range %.4f to %.4f nmol m-2 s-1\n",
            nrow(CH4_best_lgr3), mean(CH4_best_lgr3$best.flux, na.rm = TRUE),
            median(CH4_best_lgr3$best.flux, na.rm = TRUE),
            min(CH4_best_lgr3$best.flux, na.rm = TRUE), max(CH4_best_lgr3$best.flux, na.rm = TRUE)))

# ---- assembled table ------------------------------------------------------------------
CO2_best_lgr3 <- read_csv(file.path(FLUX, "CO2_best_flux_lgr_results_soil.csv"), show_col_types = FALSE)
final_dataset <- read_csv(file.path(FLUX, "auxfile_goFlux_soilflux_with_weather_formatted.csv"),
                          show_col_types = FALSE) %>%
  mutate(Vtot = SOIL_VTOT_L, Area = SOIL_COLLAR_AREA_CM2) %>%
  left_join(CO2_best_lgr3 %>% rename_with(~ paste0("CO2_", .), -UniqueID), by = "UniqueID") %>%
  left_join(CH4_best_lgr3 %>% rename_with(~ paste0("CH4_", .), -UniqueID), by = "UniqueID")
write_csv(final_dataset, file.path(FLUX, "semirigid_tree_final_complete_dataset_soil_CORRECTED.csv"))

maybe_flux2pdf(
  plot.list = maybe_flux_plot(
    flux.results = CH4_best_lgr3, dataframe = manID.lgr3, gastype = "CH4dry_ppb", shoulder = 30,
    plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"), plot.display = c("MDF", "prec", "nb.obs", "flux.term"),
    quality.check = TRUE, best.model = TRUE, p.val.disp = "round"),
  outfile = "../../outputs/figures/CH4_flux_plots_soil_CORRECTED.pdf", width = 11.6, height = 8.2)

cat("Wrote data/processed/flux/: CH4_flux_lgr_results_soil.csv, CH4_best_flux_lgr_results_soil.csv,",
    "semirigid_tree_final_complete_dataset_soil_CORRECTED.csv (", nrow(final_dataset), "rows )\n")
