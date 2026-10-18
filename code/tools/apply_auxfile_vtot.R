#!/usr/bin/env Rscript
# ==============================================================================
# apply_auxfile_vtot.R -- carry a changed system volume (Vtot) through stored goFlux
#                         results without re-fitting
# ------------------------------------------------------------------------------
# A hand-run tool, not a pipeline step (retired from stage A 2026-10-17). In a run from
# raw every goFlux fit is made with its auxfile's Vtot, which carries ANALYZER_VOLUME_CM3
# (code/lib/chamber_constants.R), so a change of the analyzer volume reaches every flux
# by re-fitting, and this script found nothing to do (the 2026-10-17 run: ratio 1.0000
# on all five campaigns). It remains for stored results that are not re-fitted, such as
# an archive of data/processed/ written at the old 70 cm3: run from the repository root
# after the auxfile scripts. Its file list does not cover every result file (not the
# 2021 CO2 goFlux results, the 2023 CO2 best.flux and CH4 goFlux results, or the
# semirigid stem result files); re-fitting is the complete route.
#
# Compares the Vtot each campaign's assembled table was fitted with ("state") against
# the Vtot its auxfile now holds ("target"; for the soil collars, SOIL_VTOT_L in
# code/lib/chamber_constants.R). In a run from raw data they agree and nothing is
# written. If they differ by one constant offset (a change of the analyzer volume,
# ANALYZER_VOLUME_CM3), every flux-unit column of every stored result is rescaled by
# Vtot_new / Vtot_old instead of re-fitting:
#
#   flux = (fitted slope, ppb s-1) x flux.term,   flux.term = Vtot*P*(1-H2O)/(R*T*Area)
#
# The fits are made on concentration, so a new Vtot rescales flux.term and every
# quantity in flux units and leaves the rest alone (slopes, C0, r2, AICc, k, g.fact,
# the model chosen by best.flux, quality.check). Checked against a goFlux + best.flux
# re-fit from the saved windows of the 2021 multi-height (441) and soil (288)
# campaigns: best.flux, MDF and flux.term agree to < 2e-9 and the standard errors to
# < 2e-6 (relative), with the same model and quality flags on every row.
# Any other difference between state and target stops the script (re-fit instead).
# Everything is computed before anything is written.
#
# Inputs:  data/processed/flux/ assembled tables and auxfiles (CAMPAIGNS below), the
#          goFlux result files (csv_files below), untagged_rescue/untagged_monthly_fluxes.csv
#          and untagged_auxfile.csv
# Outputs: the same files, rewritten only when a Vtot changed
# ==============================================================================
suppressMessages({ library(readr); library(dplyr) })
source("code/lib/chamber_constants.R")
source("code/lib/monthly_cols.R")

FLUX <- "data/processed/flux"
RESC <- "data/processed/flux/untagged_rescue"   # the rescued untagged stems (02_flux/rescue/)
num  <- function(x) suppressWarnings(as.numeric(x))
rd   <- function(f) read_csv(f, col_types = cols(.default = "c"), na = character(), progress = FALSE)

# --- campaign -> where its old (state) and new (target) Vtot live ---------------
CAMPAIGNS <- list(
  semirigid_tree = list(state  = file.path(FLUX, "semirigid_tree_final_complete_dataset.csv"),
                        target = file.path(FLUX, "auxfile_goFlux_with_weather.csv")),
  soil           = list(state  = file.path(FLUX, "semirigid_tree_final_complete_dataset_soil.csv"),
                        target = SOIL_VTOT_L),     # the soil fits' one collar geometry, not the auxfile
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

# The untagged monthly stems (02_flux/rescue/): their flux table carries the auxfile's
# UniqueIDs. Take the offset found above and require every result to be a volume the
# auxfile holds. (The picked windows, untagged_manID.rds, are an input and are not
# rewritten.)
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
    # (saved closure windows are inputs in data/raw/flux_windows/ and are never rewritten:
    # this step rescales the regenerated tables after every fit, so a rerun is identical)
    # assembled tables (the state files come last within each campaign)
    "semirigid_tree_final_complete_dataset_soil.csv",
    "tree_flux_2021_multiheight.csv", "tree_flux_2023_cross_species.csv",
    "semirigid_tree_final_complete_dataset.csv")),
  file.path(RESC, "untagged_monthly_fluxes.csv"))
stopifnot(all(file.exists(csv_files)))
new <- setNames(lapply(csv_files, function(f) rescale(rd(f), f)), csv_files)

# The monthly table with the untagged stems (03_merge_untagged_monthly.R): every row,
# tagged or untagged, carries its UniqueID, so it rescales like the tables above.
f_wu <- file.path(FLUX, "semirigid_tree_final_complete_dataset_with_untagged.csv")
new[[f_wu]] <- monthly_plain(rescale(rd(f_wu), f_wu))

# --- write -----------------------------------------------------------------------
for (f in names(new)) write_csv(new[[f]], f, na = "NA", progress = FALSE)
cat(sprintf("\nRewrote %d tables at analyzer volume %g cm3.\n", length(new), ANALYZER_VOLUME_CM3))
