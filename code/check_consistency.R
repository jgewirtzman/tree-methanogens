#!/usr/bin/env Rscript
# ==============================================================================
# check_consistency.R   --   fast invariant check over the canonical outputs
# ------------------------------------------------------------------------------
# run_all.R takes ~40 minutes and is the only end-to-end check the repo has, which
# means in practice it does not get run and inconsistencies survive. This script
# re-derives every canonical quantity from the files on disk and asserts the
# invariants, in seconds. It RUNS NOTHING and WRITES NOTHING -- it only reads.
#
# It is not a substitute for regenerating outputs. It answers a narrower question:
# given whatever is on disk right now, do the pieces agree with each other?
#
# Each check prints PASS or FAIL with the numbers, and the exit status is non-zero
# if anything failed, so it can gate a commit.
#
#   Rscript code/check_consistency.R
# ==============================================================================
source("code/lib/geometry.R")

FAILS <- 0L
chk <- function(label, ok, detail = "") {
  cat(sprintf("  [%s] %-56s %s\n", if (ok) "PASS" else "FAIL", label, detail))
  if (!ok) FAILS <<- FAILS + 1L
  invisible(ok)
}
near <- function(a, b, tol = 1e-6) is.finite(a) && is.finite(b) && abs(a - b) <= tol * max(1, abs(b))
rd <- function(p) if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL

B   <- rd("outputs/data/canonical_budget.csv")
TRP <- rd("outputs/tables/tree_flux_predictions.csv")
SOA <- rd("outputs/tables/soil_surface_annual.csv")
G   <- rd("outputs/data/scaling_full_grid.csv")
H   <- rd("outputs/data/scaling_headline.csv")
W   <- rd("outputs/data/wai_bottomup.csv")
PRF <- rd("outputs/tables/tree_band_profile.csv")
GCV <- rd("outputs/data/rf_grouped_cv.csv")
if (is.null(B)) stop("no canonical_budget.csv -- nothing to check")
v <- function(q) { x <- B$value[B$quantity == q]; if (length(x)) x[1] else NA_real_ }

SEC_YR <- 86400 * 365.25; NMOL_TO_MG <- 16e-6
CONV <- SEC_YR * NMOL_TO_MG

cat("\n== geometry ==\n")
chk("stand = nominal - gap", near(v("stand_area_m2"), v("plot_nominal_area_m2") - v("notch_area_m2")),
    sprintf("%.0f = %.0f - %.0f", v("stand_area_m2"), v("plot_nominal_area_m2"), v("notch_area_m2")))
chk("budget stand area == geometry.R", near(v("stand_area_m2"), STAND_AREA_M2),
    sprintf("%.0f vs %.0f", v("stand_area_m2"), STAND_AREA_M2))
chk("soil area = stand - basal area", near(v("soil_area_m2"), v("stand_area_m2") - v("basal_area_m2")),
    sprintf("%.3f", v("soil_area_m2")))
chk("stem area index = stem area / stand", near(v("stem_area_index"), v("stem_area_0_2m_m2")/v("stand_area_m2")))

cat("\n== tree term, recomputed from the per-stem file ==\n")
if (!is.null(TRP)) {
  S <- TRP[TRP$in_stand %in% c(TRUE, "TRUE", 1), ]
  chk("in-stand stem count matches", nrow(S) == v("n_stems"), sprintf("%d vs %.0f", nrow(S), v("n_stems")))
  chk("stem band area matches", near(sum(S$A_stem_m2), v("stem_area_0_2m_m2")),
      sprintf("%.4f vs %.4f", sum(S$A_stem_m2), v("stem_area_0_2m_m2")))
  tree <- sum(S$flux_band_nmol_m2_s * S$A_stem_m2) / v("stand_area_m2") * CONV
  chk("tree_measured recomputes", near(tree, v("tree_measured_mg_m2_yr"), 1e-4),
      sprintf("%.6f vs %.6f", tree, v("tree_measured_mg_m2_yr")))
  t2m <- sum(S$flux_2m_nmol_m2_s * S$A_stem_m2) / v("stand_area_m2") * CONV
  chk("tree_at_2m_basis recomputes", near(t2m, v("tree_at_2m_basis_mg_m2_yr"), 1e-4),
      sprintf("%.6f vs %.6f", t2m, v("tree_at_2m_basis_mg_m2_yr")))
  chk("no GEN_* level lost (species ladder reached the inventory)",
      any(grepl("^GEN_", S$species_model)),
      sprintf("%d GEN_ stems", sum(grepl("^GEN_", S$species_model))))
}

cat("\n== soil term, and the basis both terms share ==\n")
if (!is.null(SOA)) {
  raw <- mean(SOA$flux_nmol_m2_s) * CONV                       # per m2 of SOIL
  grd <- raw * v("soil_area_m2") / v("stand_area_m2")          # per m2 of GROUND
  chk("soil term is on the GROUND basis (x soil fraction)", near(grd, v("soil_mg_m2_yr"), 1e-4),
      sprintf("%.4f vs %.4f (per-soil would be %.4f)", grd, v("soil_mg_m2_yr"), raw))
}
chk("net = tree + soil", near(v("tree_measured_mg_m2_yr") + v("soil_mg_m2_yr"), v("net_measured_mg_m2_yr"), 1e-9),
    sprintf("%.4f", v("net_measured_mg_m2_yr")))
# The MONTHLY series must be on the same basis as the annual scalars. It was not: the
# soil-fraction rescaling was applied to the annual term only, so canonical_monthly.csv
# stayed per-m2-of-soil while the tree column beside it was per-m2-of-ground, and
# fig09_budget.R plots that series labelled "per sq.m ground". 0.45%, and invisible
# to every other check here, which is why this one exists.
MON <- rd("outputs/data/canonical_monthly.csv")
if (!is.null(MON)) {
  DPM <- 365.25/12
  chk("monthly soil sums to the annual soil term",
      near(sum(MON$soil_mg_m2_d)*DPM, v("soil_mg_m2_yr"), 5e-3),
      sprintf("%.4f vs %.4f", sum(MON$soil_mg_m2_d)*DPM, v("soil_mg_m2_yr")))
  chk("monthly tree sums to the annual tree term",
      near(sum(MON$tree_mg_m2_d)*DPM, v("tree_measured_mg_m2_yr"), 5e-3),
      sprintf("%.4f vs %.4f", sum(MON$tree_mg_m2_d)*DPM, v("tree_measured_mg_m2_yr")))
  chk("monthly plot = monthly tree + monthly soil (same basis)",
      all(abs(MON$plot_mg_m2_d - (MON$tree_mg_m2_d + MON$soil_mg_m2_d)) < 1e-12))
}
chk("tree %% of soil consistent",
    near(v("tree_pct_of_soil"), 100*abs(v("tree_measured_mg_m2_yr"))/abs(v("soil_mg_m2_yr")), 1e-6),
    sprintf("%.4f%%", v("tree_pct_of_soil")))

cat("\n== scaling grid ==\n")
if (!is.null(G)) {
  chk("grid is 6 x 5 x 4 x 2 = 240 rows", nrow(G) == 240, sprintf("%d rows", nrow(G)))
  mm <- unique(round(G$measured_mg, 9))
  chk("measured band invariant across the grid", length(mm) == 1, sprintf("%d distinct", length(mm)))
  chk("grid band == canonical tree term", near(mm[1], v("tree_measured_mg_m2_yr"), 1e-6),
      sprintf("%.6f vs %.6f", mm[1], v("tree_measured_mg_m2_yr")))
  chk("total = measured + extrapolated (all rows)",
      all(abs(G$total_mg - (G$measured_mg + G$extrapolated_mg)) < 1e-9))
  chk("net = soil + total (all rows)",
      all(abs(G$net_mg - (v("soil_mg_m2_yr") + G$total_mg)) < 1e-6))
  chk("sink in every combination", all(G$net_mg < 0), sprintf("%d/%d", sum(G$net_mg < 0), nrow(G)))
  if (!is.null(H)) chk("headline row is present in the grid",
    any(G$flux == H$flux[1] & G$branch == H$branch[1] & G$bole == H$bole[1] & G$WAI == H$WAI[1]),
    sprintf("%s | %s", H$flux[1], H$WAI[1]))
}

cat("\n== WAI provenance ==\n")
if (!is.null(W)) {
  chk("wai_bottomup computed over the canonical stand", near(W$stand_area_m2[1], STAND_AREA_M2),
      sprintf("%.0f m2, %d stems", W$stand_area_m2[1], W$n_stems[1]))
  chk("wai_bottomup stem count == budget", W$n_stems[1] == v("n_stems"),
      sprintf("%d vs %.0f", W$n_stems[1], v("n_stems")))
  if (!is.null(G)) chk("grid WAI axis carries the bottom-up values",
    all(sapply(W$wai, function(x) any(abs(as.numeric(sub(" .*","",unique(G$WAI))) - x) < 0.005))),
    paste(sprintf("%.2f", W$wai), collapse=" / "))
}

cat("\n== band profile export (the grid's 2 m anchor) ==\n")
if (!is.null(PRF) && !is.null(TRP)) {
  Ps <- PRF[PRF$in_stand %in% c(TRUE,"TRUE",1), ]
  Ts <- TRP[TRP$in_stand %in% c(TRUE,"TRUE",1), ]
  fcol <- grep("^f_i[0-9]+$", names(Ps), value = TRUE)
  chk("profile rows match the per-stem file", nrow(Ps) == nrow(Ts), sprintf("%d vs %d", nrow(Ps), nrow(Ts)))
  if (nrow(Ps) == nrow(Ts) && length(fcol))
    chk("profile's last interval == calibrated flux_2m",
        max(abs(Ps[[tail(fcol,1)]] - Ts$flux_2m_nmol_m2_s)) < 1e-9,
        sprintf("max diff %.2e", max(abs(Ps[[tail(fcol,1)]] - Ts$flux_2m_nmol_m2_s))))
}

cat("\n== model skill ==\n")
if (!is.null(GCV)) {
  tg <- GCV$grouped_r2[GCV$model == "TreeRF"]; sg <- GCV$grouped_r2[GCV$model == "SoilRF"]
  chk("grouped CV in budget matches rf_grouped_cv.csv",
      near(v("tree_r2_grouped_cv"), tg[1], 1e-6) && near(v("soil_r2_grouped_cv"), sg[1], 1e-6),
      sprintf("tree %.4f, soil %.4f", tg[1], sg[1]))
  chk("OOB exceeds grouped CV (expected; flags if reversed)",
      v("tree_r2_oob") > v("tree_r2_grouped_cv"),
      sprintf("OOB %.4f vs grouped %.4f, ratio %.2fx",
              v("tree_r2_oob"), v("tree_r2_grouped_cv"), v("tree_r2_oob")/v("tree_r2_grouped_cv")))
}

# ------------------------------------------------------------------------------
# DATA PLUMBING (added 2026-09-30). The checks above test whether outputs agree with
# each other. These test the failure modes that slipped past them:
#   - two scripts writing one file, so whichever ran last silently won (the 2021 and
#     2023 tree-flux tables shared a name; the 2023 run overwrote the 2021 one)
#   - a table whose rows quietly stop accounting for the measurements (Figure 3 read
#     the model's training rows as if they were every deployment, and lost 17 ash trees)
#   - a derived copy left behind when its source changed
# ------------------------------------------------------------------------------
cat("\n== data plumbing ==\n")

# 1. one writer per file. Static scan of live code (archive excluded). A file may
#    appear here only with the reason it is allowed.
SEQUENTIAL_STAGES <- c(   # interactive 2020-21 soil processing: goFlux, then the December
  # fix merges into the same files, then the geometry fix rewrites CH4. Run once, in that
  # order, by hand; not part of run_all.R. The model reads the end state.
  "CH4_best_flux_lgr_results_soil.csv", "CH4_flux_lgr_results_soil.csv",
  "CO2_best_flux_lgr_results_soil.csv", "CO2_flux_lgr_results_soil.csv",
  "lgr_manual_identification_results_soil.csv", "semirigid_tree_final_complete_dataset_soil.csv")
live <- setdiff(list.files("code", "\\.R$", recursive = TRUE, full.names = TRUE),
                list.files("code/archive", "\\.R$", recursive = TRUE, full.names = TRUE))
wpat <- "(write\\.csv|write_csv|write\\.table|saveRDS|save|ggsave|fwrite|write_tsv)\\s*\\("
writers <- list()
for (f in live) for (ln in readLines(f, warn = FALSE)) {
  if (grepl("^\\s*#", ln) || !grepl(wpat, ln)) next
  for (p in regmatches(ln, gregexpr("[\"'][^\"']+\\.(csv|rds|RData|rda|tsv|png|pdf|txt)[\"']", ln))[[1]]) {
    k <- basename(gsub("[\"']", "", p)); writers[[k]] <- union(writers[[k]], f)
  }
}
multi <- Filter(function(w) length(w) > 1, writers)
multi <- multi[!names(multi) %in% SEQUENTIAL_STAGES]
chk("every output file has exactly one writing script", !length(multi),
    if (length(multi)) paste(sprintf("%s <- %s", names(multi),
                                     vapply(multi, paste, "", collapse = " + ")), collapse = "; ") else
    sprintf("%d files scanned", length(writers)))

# 2. the retired shared names must not come back, in code or on disk
RETIRED <- c("methanogen_tree_flux_complete_dataset.csv", "CH4_best_flux_lgr_results.csv",
             "CO2_best_flux_lgr_results.csv", "CH4_flux_lgr_results.csv", "CO2_flux_lgr_results.csv",
             "lgr_manual_identification_results.csv")
hits <- live[vapply(live, function(f) any(grepl(paste(gsub(".", "\\.", RETIRED, fixed = TRUE),
                                                      collapse = "|"), readLines(f, warn = FALSE))), TRUE)]
hits <- setdiff(hits, "code/check_consistency.R")
hits <- hits[!grepl("assemble_campaign_flux\\.R$", hits)]   # names them in its history note
chk("retired shared flux filenames absent from code", !length(hits), paste(hits, collapse = ", "))
chk("retired shared flux filenames absent from data/processed/flux",
    !any(file.exists(file.path("data/processed/flux", RETIRED))),
    paste(RETIRED[file.exists(file.path("data/processed/flux", RETIRED))], collapse = ", "))

# 3. the measurement table accounts for every deployment, and its training subset
#    is exactly what the model was fitted on
MT <- rd("outputs/data/flux_measurements_tree.csv")
if (!is.null(MT)) {
  fd <- "data/processed/flux"; nn <- function(x) sum(!is.na(x))
  src <- c(monthly_2020_2021 = nn(rd(file.path(fd, "semirigid_tree_final_complete_dataset_with_untagged.csv"))$CH4_best.flux.x),
           `2021_multiheight` = nn(rd(file.path(fd, "tree_flux_2021_multiheight.csv"))$CH4_best.flux),
           `2023_cross_species` = nn(rd(file.path(fd, "tree_flux_2023_cross_species.csv"))$CH4_best.flux))
  got <- table(factor(MT$campaign, levels = names(src)))
  chk("measurement table = every deployment in the campaign files", all(got == src),
      paste(sprintf("%s %d/%d", names(src), as.integer(got), src), collapse = ", "))
  if (file.exists("outputs/models/TRAINING_DATA.RData")) {
    e <- new.env(); load("outputs/models/TRAINING_DATA.RData", envir = e)
    chk("in_rf_training rows == the frozen model's training rows",
        sum(MT$in_rf_training %in% c(TRUE, "TRUE")) == nrow(e$tree_train_complete),
        sprintf("%d vs %d", sum(MT$in_rf_training %in% c(TRUE, "TRUE")), nrow(e$tree_train_complete)))
  }
  chk("every excluded deployment states its reason",
      all(!is.na(MT$exclusion_reason[!MT$in_rf_training %in% c(TRUE, "TRUE")])))
  if (file.exists("data/compiled/flux_measurements_tree.csv"))
    chk("data/compiled copy matches outputs/data", identical(unname(tools::md5sum("data/compiled/flux_measurements_tree.csv")),
                                                           unname(tools::md5sum("outputs/data/flux_measurements_tree.csv"))),
        "copy outputs/data/flux_measurements_tree.csv after rerunning the export")
}

# 4. internal gas: recalibrated (no zero clamp) and carried through to the merged table
GS <- rd("data/processed/internal_gas/sample_data_only.csv")
MG <- rd("data/processed/integrated/merged_tree_dataset_final.csv")
if (!is.null(GS)) chk("no internal-gas CH4 is exactly 0 (the old clamp)", !any(GS$CH4_concentration == 0, na.rm = TRUE),
                      sprintf("%d zeros", sum(GS$CH4_concentration == 0, na.rm = TRUE)))
if (!is.null(GS) && !is.null(MG))
  chk("merged table carries the current gas calibration",
      isTRUE(all.equal(sort(GS$CH4_concentration), sort(MG$CH4_concentration[!is.na(MG$CH4_concentration)]))),
      "rerun code/03_merge/02_harmonize_all_data.R after 03_process_internal_gas.R")

cat(sprintf("\n%s  %d check(s) failed\n\n", if (FAILS == 0L) "ALL CONSISTENT." else "INCONSISTENT.", FAILS))
if (FAILS > 0L) quit(status = 1L)
