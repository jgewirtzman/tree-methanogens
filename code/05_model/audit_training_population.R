# ==============================================================================
# audit_training_population.R -- which stems should the TreeRF be trained on?
# ------------------------------------------------------------------------------
# Records the evidence behind two 2026-09-30 decisions, so it can be rerun rather
# than remembered:
#   1. dead stems are excluded from training (the inventory is live stems only), and
#      a live/dead predictor does not rescue them;
#   2. the stand total moves with the restored deployments (2023 non-numeric tags,
#      untagged monthly stems) by more than random-forest seed noise.
#
# For each training population x seed it refits the TreeRF with the canonical
# hyper-parameters, scores it through predict_tree_flux_current.R in SANDBOX mode
# (temporary directory; outputs/models and outputs/tables are never touched), and
# reports out-of-bag R2 and the measured-band stand total.
#
# Run from the repo root after 02_rf_models.R and the CORE chain (it needs the
# inventory, moisture surface and climatologies):
#   Rscript code/05_model/audit_training_population.R
# ~1 min per fit. Writes outputs/audit/training_population_audit.{txt,csv}
# ==============================================================================
source("code/lib/outputs.R"); source("code/lib/geometry.R")
suppressMessages({ library(ranger); library(dplyr) })

load("outputs/models/TRAINING_CANDIDATES.RData")     # tree_train_candidates, dead flagged
e <- new.env(); load("outputs/models/RF_MODELS.RData", envir = e)
CONV  <- 86400 * 365.25 * 16e-6
SEEDS <- c(42, 1, 2)
cand  <- tree_train_candidates
restored <- grepl("^2023_|^UNTAG", cand$tree_id)

design <- function(d, with_dead = FALSE) {
  X <- data.frame(species = factor(d$species_clean), dbh_m = d$dbh_m,
                  soil_moisture_at_tree = d$soil_moisture_at_tree,
                  soil_temp_C_mean = d$soil_temp_C_mean, air_temp_C_mean = d$air_temp_C_mean,
                  height_cm = ifelse(is.na(d$measurement_height_cm), 125, d$measurement_height_cm))
  if (with_dead) X$dead_stem <- as.numeric(d$dead_stem)
  X
}

VARIANTS <- list(
  live_only           = list(keep = !cand$dead_stem,             dead_pred = FALSE),
  all_stems           = list(keep = rep(TRUE, nrow(cand)),       dead_pred = FALSE),
  live_dead_predictor = list(keep = rep(TRUE, nrow(cand)),       dead_pred = TRUE),
  live_no_restored    = list(keep = !cand$dead_stem & !restored, dead_pred = FALSE))

rows <- list()
for (v in names(VARIANTS)) for (s in SEEDS) {
  d <- cand[VARIANTS[[v]]$keep, ]; X <- design(d, VARIANTS[[v]]$dead_pred)
  TreeRF <- ranger(x = X, y = d$y_asinh, num.trees = e$TreeRF$num.trees, min.node.size = e$TreeRF$min.node.size,
                   mtry = floor(sqrt(ncol(X))), importance = "impurity", num.threads = 1, oob.error = TRUE, seed = s)
  SoilRF <- e$SoilRF
  sb <- tempfile("treepred_"); dir.create(sb)
  tree_train_complete <- d; X_tree <- X
  save(TreeRF, SoilRF, file = file.path(sb, "RF_MODELS.RData"))
  save(tree_train_complete, X_tree, file = file.path(sb, "TRAINING_DATA.RData"))
  st <- system2("Rscript", "code/06_upscale/predict_tree_flux_current.R", stdout = FALSE, stderr = FALSE,
                env = paste0("TREE_PRED_SANDBOX=", sb))
  P <- read.csv(file.path(sb, "tree_flux_predictions.csv"))
  P <- P[P$in_stand %in% c(TRUE, "TRUE"), ]
  band <- sum(P$flux_band_nmol_m2_s * P$A_stem_m2) / STAND_AREA_M2 * CONV
  rows[[length(rows) + 1]] <- data.frame(variant = v, seed = s, rows = nrow(d), dead_rows = sum(d$dead_stem),
                                         oob_r2 = round(TreeRF$r.squared, 4), band_mg = round(band, 3), status = st)
  unlink(sb, recursive = TRUE)
  cat(sprintf("%-20s seed %2d  rows %4d  OOB %.3f  band %.3f mg\n", v, s, nrow(d), TreeRF$r.squared, band))
}
R <- do.call(rbind, rows)
write.csv(R, out_path("training_population_audit.csv"), row.names = FALSE)

S <- R %>% group_by(variant) %>% summarise(rows = first(rows), dead = first(dead_rows),
       oob = sprintf("%.3f-%.3f", min(oob_r2), max(oob_r2)),
       band = sprintf("%.2f-%.2f", min(band_mg), max(band_mg)), .groups = "drop")
sink(out_path("training_population_audit.txt"))
cat("TRAINING POPULATION AUDIT (code/05_model/audit_training_population.R)\n")
cat(sprintf("built %s; seeds %s; measured-band stand total in mg CH4 m-2 yr-1\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M"), paste(SEEDS, collapse = ", ")))
print(as.data.frame(S), row.names = FALSE)
cat("\nlive_only is the canonical training set (seed 42 reproduces 02_rf_models.R).\n")
cat("A live/dead predictor is only worth keeping if it moves live_dead_predictor away from\n")
cat("all_stems; restored rows matter if live_no_restored lies outside live_only's seed range.\n")
sink()
cat(readLines(out_path("training_population_audit.txt")), sep = "\n")
