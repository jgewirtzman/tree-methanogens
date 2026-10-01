#!/usr/bin/env Rscript
# ==============================================================================
# 01_apply_auxfile_vtot.R -- carry a changed system volume (Vtot) through the stored
# goFlux results without repeating the interactive window selection.
# ------------------------------------------------------------------------------
# WHY THIS EXISTS (2026-09-30). The analyzer volume in Vtot was corrected from
# 70 to 28 cm3 (code/lib/chamber_constants.R). The auxfile scripts pick that up
# on a rerun, but the goFlux scripts that turn an auxfile into fluxes begin with
# click.peak2(), and the saved window selections survive for only two of the
# five campaigns. goFlux does not need rerunning for a volume change, though:
#
#   flux = (fitted slope, ppb s-1) x flux.term,   flux.term = Vtot*P*(1-H2O)/(R*T*Area)
#
# The fits are made on concentration, so a new Vtot rescales flux.term and every
# quantity in flux units by Vtot_new/Vtot_old and leaves the rest alone (slopes,
# C0, r2, AICc, k, g.fact, the model chosen by best.flux, quality.check). Checked
# against a real goFlux + best.flux rerun from the saved windows of the 2021
# multi-height (441) and soil (288) campaigns: best.flux, MDF and flux.term agree
# to < 2e-9 and the standard errors to < 2e-6 (relative), with the same model and
# quality flags on every row. A rerun under the installed goFlux (0.2.0) also
# moves 195 of the 441 stored 2021 fluxes by up to 0.15 % with no volume change
# at all, so rescaling is the only way to change the volume and nothing else.
#
# Run from the repo root AFTER the auxfile scripts:
#   (cd code/02_flux/semirigid && Rscript 03_prep_tree_auxfile.R && Rscript 03_prep_soil_auxfile.R)
#   (cd code/02_flux/static    && Rscript 01_prep_auxfile_2021.R && Rscript 02_prep_auxfile_2023.R)
#   Rscript code/02_flux/01_apply_auxfile_vtot.R
#
# Old Vtot is whatever the campaign's assembled table holds; new Vtot is what its
# regenerated auxfile holds (for the untagged stems, also rerun the archived
# rev_rescue_untagged_prep.R with its output pointed at outputs/data/). A campaign whose two already agree is left alone, so
# the script is safe to rerun. Everything is computed before anything is written.
# ==============================================================================
suppressMessages({ library(readr); library(dplyr) })
source("code/lib/chamber_constants.R")
source("code/lib/monthly_cols.R")

FLUX <- "data/processed/flux"
RESC <- "data/processed/flux/untagged_rescue"   # the rescued untagged stems (02_flux/rescue/)
num  <- function(x) suppressWarnings(as.numeric(x))
rd   <- function(f) read_csv(f, col_types = cols(.default = "c"), na = character(), progress = FALSE)

# Soil collars: the soil fits use one fixed geometry (SOIL_VTOT_L, chamber_constants.R),
# so the target is that constant, not the auxfile.

# --- campaign -> where its old (state) and new (target) Vtot live ---------------
CAMPAIGNS <- list(
  semirigid_tree = list(state  = file.path(FLUX, "semirigid_tree_final_complete_dataset.csv"),
                        target = file.path(FLUX, "auxfile_goFlux_with_weather.csv")),
  soil           = list(state  = file.path(FLUX, "semirigid_tree_final_complete_dataset_soil.csv"),
                        target = SOIL_VTOT_L),
  `2021_multiheight`   = list(state  = file.path(FLUX, "tree_flux_2021_multiheight.csv"),
                              target = file.path(FLUX, "goflux_auxfile.csv")),
  `2023_cross_species` = list(state  = file.path(FLUX, "tree_flux_2023_cross_species.csv"),
                              target = file.path(FLUX, "ymf2023_goflux_auxfile.csv")))

# UniqueIDs are unique across campaigns (asserted), so one lookup serves every file.
vt <- bind_rows(lapply(names(CAMPAIGNS), function(cp) {
  cfg <- CAMPAIGNS[[cp]]
  st  <- rd(cfg$state) %>% transmute(UniqueID, V_old = num(Vtot))
  if (is.numeric(cfg$target)) st$V_new <- cfg$target
  else st <- left_join(st, rd(cfg$target) %>% transmute(UniqueID, V_new = num(Vtot)), by = "UniqueID")
  mutate(st, campaign = cp)
}))
stopifnot(!anyNA(vt$V_old), !anyNA(vt$V_new))

# The only thing allowed to differ between state and target is the analyzer term,
# i.e. one constant offset shared by every chamber of every campaign that moved.
moved <- abs(vt$V_new / vt$V_old - 1) > 1e-12
dV <- (vt$V_old - vt$V_new)[moved] * 1000
if (length(dV) && diff(range(dV)) > 1e-6)
  stop("state and target Vtot differ by more than one constant offset (", paste(signif(range(dV), 6), collapse = " to "),
       " cm3): something other than the analyzer volume changed. Rerun goFlux for that campaign instead.", call. = FALSE)

# The rescued untagged monthly stems. Their flux table keeps the UniqueIDs of the
# click step, which the auxfile has since re-keyed (species corrected from field
# notes), so the two cannot be joined on UniqueID; rev_rescue_untagged_finalize.R
# joined them on start time, which the flux table does not carry. Take the offset
# found above instead and require every result to be a volume the regenerated
# auxfile holds. (That auxfile comes from archive/revision/rev_rescue_untagged_prep.R;
# untagged_fluxes.csv and untagged_manID.rds predate the corrected geometry and
# are left alone.)
un  <- rd(file.path(RESC, "untagged_monthly_fluxes.csv")) %>% transmute(UniqueID, V_old = num(Vtot), campaign = "untagged")
aux <- round(num(rd(file.path(RESC, "untagged_auxfile.csv"))$Vtot), 9)
off <- if (length(dV)) mean(dV) / 1000 else 0
un$V_new <- if (all(round(un$V_old, 9) %in% aux)) un$V_old else un$V_old - off
if (!all(round(un$V_new, 9) %in% aux))
  stop("untagged_monthly_fluxes.csv and untagged_auxfile.csv disagree by something other than the analyzer offset", call. = FALSE)
vt <- bind_rows(vt, un)
stopifnot(!anyDuplicated(vt$UniqueID))
vt$r <- vt$V_new / vt$V_old
moved <- abs(vt$r - 1) > 1e-12
cat("Vtot change by campaign (old -> new, L; flux ratio):\n")
vt %>% group_by(campaign) %>%
  summarise(n = n(), old = sprintf("%.3f-%.3f", min(V_old), max(V_old)),
            new = sprintf("%.3f-%.3f", min(V_new), max(V_new)),
            ratio = sprintf("%.4f-%.4f", min(r), max(r)), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
if (!any(moved)) { cat("\nEvery stored table already holds its auxfile's Vtot. Nothing to do.\n"); quit(status = 0) }
cat(sprintf("\nOffset: %.1f cm3 on %d of %d measurements.\n\n", off * 1000, sum(moved), nrow(vt)))

# --- what scales ----------------------------------------------------------------
# Columns in flux units, bare (goFlux result files) or with a gas prefix and an
# optional .x/.y merge suffix (assembled tables). MDF.lim is best.flux's display
# copy of MDF, rounded to 2 significant figures, so it is re-derived, not scaled.
SCALED <- "^((CO2|CH4)_)?(LM\\.flux|LM\\.SE|HM\\.flux|HM\\.SE|MDF|flux\\.term|best\\.flux)(\\.[xy])?$"
rescale <- function(x, f) {
  i <- match(x$UniqueID, vt$UniqueID); hit <- !is.na(i)
  r <- vt$r[i]; r[!hit] <- 1
  sc <- grep(SCALED, names(x), value = TRUE)
  for (v in sc) { y <- num(x[[v]]); x[[v]] <- ifelse(hit & !is.na(y), y * r, x[[v]]) }
  for (v in grep("MDF\\.lim", names(x), value = TRUE)) {
    y <- num(x[[sub("MDF\\.lim", "MDF", v)]]); x[[v]] <- ifelse(hit & !is.na(y), signif(y, 2), x[[v]])
  }
  if ("Vtot" %in% names(x)) {
    # a file's own Vtot can sit in another state (the December soil windows carry the
    # pre-correction 17.75 L collar volume): move it by the same offset, never overwrite
    y <- num(x$Vtot); x$Vtot <- ifelse(hit & !is.na(y), y - (vt$V_old[i] - vt$V_new[i]), x$Vtot)
  }
  # the 2023 table carries the auxfile's volume breakdown alongside Vtot
  if (all(c("analyzer_cell_volume_cm3", "total_system_volume_cm3", "total_system_volume_L") %in% names(x))) {
    x$analyzer_cell_volume_cm3 <- ifelse(hit, ANALYZER_VOLUME_CM3, x$analyzer_cell_volume_cm3)
    x$total_system_volume_cm3  <- ifelse(hit, vt$V_new[i] * 1000,   x$total_system_volume_cm3)
    x$total_system_volume_L    <- ifelse(hit, vt$V_new[i],          x$total_system_volume_L)
  }
  cat(sprintf("  %-62s %5d of %5d rows, %2d flux columns\n", basename(f), sum(hit), nrow(x), length(sc)))
  x
}

csv_files <- c(
  file.path(FLUX, c(
    # goFlux / best.flux result files, per campaign
    "CH4_best_flux_lgr_results_soil.csv", "CO2_best_flux_lgr_results_soil.csv",
    "CH4_flux_lgr_results_soil.csv",      "CO2_flux_lgr_results_soil.csv",
    "CH4_best_flux_lgr_results_2021_multiheight.csv", "CO2_best_flux_lgr_results_2021_multiheight.csv",
    "CH4_flux_lgr_results_2021_multiheight.csv",
    "CH4_best_flux_lgr_results_2023_cross_species.csv", "CO2_flux_lgr_results_2023_cross_species.csv",
    # saved window selections: Vtot only, so a later goFlux run from them is right
    "lgr_manual_identification_results_soil.csv", "lgr_manual_identification_results_december_soil.csv",
    "lgr_manual_identification_2021_multiheight.csv", "lgr_manual_identification_2021_multiheight_final.csv",
    # assembled tables (the state files come last within each campaign)
    "semirigid_tree_final_complete_dataset_soil_CORRECTED.csv", "semirigid_tree_final_complete_dataset_soil.csv",
    "tree_flux_2021_multiheight.csv", "tree_flux_2023_cross_species.csv",
    "semirigid_tree_final_complete_dataset.csv")),
  file.path(RESC, "untagged_monthly_fluxes.csv"))
csv_files <- csv_files[file.exists(csv_files)]
new <- setNames(lapply(csv_files, function(f) rescale(rd(f), f)), csv_files)

# The monthly table with the rescued untagged stems. Its tagged rows carry a
# UniqueID; the 45 untagged rows were appended without one, in the row order of
# untagged_monthly_fluxes.csv (archive/revision/rev_merge_untagged_monthly.R), and
# hold only CH4_best.flux (CH4_best.flux.x/.y in files written before 2026-10-01;
# monthly_plain() reads either). Match them by position and prove it by value.
f_wu <- file.path(FLUX, "semirigid_tree_final_complete_dataset_with_untagged.csv")
wu <- monthly_plain(rescale(rd(f_wu), f_wu))
um_old <- rd(file.path(RESC, "untagged_monthly_fluxes.csv")); um_new <- new[[file.path(RESC, "untagged_monthly_fluxes.csv")]]
k <- which(wu$UniqueID == "NA" | wu$UniqueID == "")
stopifnot(length(k) == nrow(um_old),
          isTRUE(all.equal(num(wu$CH4_best.flux[k]), num(um_old$CH4_best.flux))))
wu$CH4_best.flux[k] <- um_new$CH4_best.flux
cat(sprintf("  %-62s %5d untagged rows matched by position\n", "", length(k)))
new[[f_wu]] <- wu

# --- write -----------------------------------------------------------------------
for (f in names(new)) write_csv(new[[f]], f, na = "NA", progress = FALSE)
cat(sprintf("\nRewrote %d tables at analyzer volume %g cm3.\n", length(new), ANALYZER_VOLUME_CM3))
cat("Downstream, in order: 03_merge/04_harmonize_all_data.R, 05_model/01_load_and_prep_data.R,\n",
    " 05_model/02_rf_models.R, code/run_all.R, zenodo/01_compile_datasets.R, check_consistency.R\n", sep = "")
