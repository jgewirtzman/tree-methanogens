# ==============================================================================
# 02_rf_models.R -- fit the stem and soil random forests
# ------------------------------------------------------------------------------
# Stage B. Run from code/05_model/ (workdir = script), after 01_load_and_prep_data.R.
#
# Input
#   data/processed/integrated/rf_workflow_input_data_with_2023.RData (01)
# Outputs
#   outputs/models/RF_MODELS.RData            TreeRF, SoilRF
#   outputs/models/TRAINING_DATA.RData        training rows and predictor frames
#                                             (tree_train_complete, X_tree,
#                                             soil_train_complete, X_soil)
#   outputs/models/TRAINING_CANDIDATES.RData  every trainable stem row, dead stems
#                                             flagged (05_audit_training_population.R)
#
# Response: CH4 flux in nmol m-2 s-1 on the measured scale (no transform). The
# stem response is chamber-harmonised (section 2).
# Stem forest: species (taxonomic ladder, section 3), DBH, soil moisture, soil
#   temperature, measurement height.
# Soil forest: soil moisture, soil temperature, month.
# Both: 10,000 trees (code/lib/rf_spec.R), min.node.size 5, mtry = floor(sqrt(p)),
#   seed 42. Air temperature is not a predictor; the with-air forests are the
#   sensitivity 05_model/15_rf_air_temperature_test.R.
# ==============================================================================

source("../lib/rf_spec.R")   # RF_NUM_TREES, RF_THREADS
source("../lib/dead_stems.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ranger)
})

cat("=== RANDOM FORESTS: STEM AND SOIL CH4 FLUX ===\n")
cat("Started at:", format(Sys.time()), "\n\n")

# =============================================================================
# 1. INPUT
# =============================================================================
load("../../data/processed/integrated/rf_workflow_input_data_with_2023.RData")
TREE_JULY <- rf_workflow_data$TREE_JULY     # 2021 multi-height + 2023 cross-species (rigid)
TREE_YEAR <- rf_workflow_data$TREE_YEAR     # monthly survey 2020-21 (semi-rigid)
SOIL_YEAR <- rf_workflow_data$SOIL_YEAR     # monthly soil collars 2020-21

cat("Flux ranges (nmol m-2 s-1):\n")
cat("  TREE_JULY: [", round(min(TREE_JULY$stem_flux_nmol_m2_s), 6), ",",
    round(max(TREE_JULY$stem_flux_nmol_m2_s), 6), "]\n")
cat("  TREE_YEAR: [", round(min(TREE_YEAR$stem_flux_nmol_m2_s), 6), ",",
    round(max(TREE_YEAR$stem_flux_nmol_m2_s), 6), "]\n")
cat("  SOIL_YEAR: [", round(min(SOIL_YEAR$soil_flux_nmol_m2_s), 6), ",",
    round(max(SOIL_YEAR$soil_flux_nmol_m2_s), 6), "]\n")
# Species, DBH and position of every stem row are its tree's (data/tables/trees.csv,
# by tree_id; 01_load_and_prep_data.R).
stopifnot(!anyNA(TREE_JULY$species), !anyNA(TREE_YEAR$species))

# =============================================================================
# 2. STEM ROWS AND CHAMBER HARMONISATION
# =============================================================================
tree_combined <- bind_rows(TREE_JULY, TREE_YEAR) %>% filter(!is.na(stem_flux_nmol_m2_s))
stopifnot(all(tree_combined$chamber_type %in% c("rigid", "semirigid")))

# DBH as measured (01 from trees.csv). Trees with no DBH measurement (six 2021 trees)
# keep dbh_m NA; ranger (0.17) handles a missing value natively, as it does for any
# missing covariate.
stopifnot(all(tree_combined$dbh_m > 0 & tree_combined$dbh_m < 2, na.rm = TRUE))
cat(sprintf("  DBH missing (not measured; kept, ranger handles NA): %d deployments on %d trees\n",
            sum(is.na(tree_combined$dbh_m)), n_distinct(tree_combined$tree_id[is.na(tree_combined$dbh_m)])))

# Chamber offset. On the asinh scale, z ~ I_semirigid + the stem forest's predictors
# (species, DBH, soil moisture, soil temperature, measurement height): beta1 is the
# systematic offset of the monthly (semi-rigid) chamber relative to the rigid one.
# Only semi-rigid rows are corrected, z_corr = z - beta1, and the response is
# sinh(z_corr). Species enters only where it has >= 5 rows; a missing height (75 cm
# and root crown) is 125 cm, as in the forest; rows missing DBH or soil values are
# left out of this regression only (lm complete cases). The log also reports two
# alternatives (with air temperature; conditions only).
cal_df <- tree_combined %>%
  mutate(z = asinh(stem_flux_nmol_m2_s),
         I_semirigid = as.numeric(chamber_type == "semirigid"),
         height_cm = ifelse(is.na(measurement_height_cm), 125, measurement_height_cm))
stopifnot(sum(cal_df$I_semirigid == 0) >= 10, sum(cal_df$I_semirigid == 1) >= 10)
sp_ok <- cal_df %>% count(species) %>% filter(n >= 5, !is.na(species)) %>% pull(species)
cal_df$species_ctl <- ifelse(cal_df$species %in% sp_ok, cal_df$species, "OTHER")
chamber_fit <- lm(z ~ I_semirigid + species_ctl + dbh_m + soil_moisture_abs + soil_temp_C + height_cm,
                  data = cal_df)
for (alt in list(c("with air temperature", "z ~ I_semirigid + air_temp_C + soil_temp_C + soil_moisture_abs + species_ctl"),
                 c("conditions only", "z ~ I_semirigid + soil_temp_C + soil_moisture_abs + species_ctl"))) {
  fa <- lm(as.formula(alt[2]), data = cal_df); ca <- summary(fa)$coefficients
  cat(sprintf("  (alternative, not used) %s: beta1 = %.4f (SE %.4f), n = %d\n", alt[1],
              ca["I_semirigid", 1], ca["I_semirigid", 2], nobs(fa)))
}
cf <- summary(chamber_fit)$coefficients
beta1 <- unname(cf["I_semirigid", "Estimate"]); beta1_se <- unname(cf["I_semirigid", "Std. Error"])
ci <- beta1 + c(-1.96, 1.96) * beta1_se
cat(sprintf("  Chamber bias beta1 = %.4f (SE %.4f, 95%% CI %.4f to %.4f), model R2 = %.3f, n = %d\n",
            beta1, beta1_se, ci[1], ci[2], summary(chamber_fit)$r.squared, nobs(chamber_fit)))
cat(sprintf("  n rigid = %d, n semirigid = %d\n", sum(cal_df$I_semirigid == 0), sum(cal_df$I_semirigid == 1)))

z <- asinh(tree_combined$stem_flux_nmol_m2_s)
tree_combined$stem_flux_corrected <- sinh(ifelse(tree_combined$chamber_type == "semirigid", z - beta1, z))

# =============================================================================
# 3. TRAINING FRAMES
# =============================================================================
# Temperatures at the measurement: air_temp_C is the air temperature measured at the
# deployment (carried, not a predictor); soil_temp_C the soil temperature read at the
# tree or collar (attach_soil(); NA where there is none within 3 days). At prediction
# the same names take the monthly 24-h climatologies. A value outside the physical
# range stops the build (corrections belong in the raw files, data/raw/CORRECTIONS.md).
check_temps <- function(df) {
  df$air_temp_C  <- suppressWarnings(as.numeric(df$air_temp_C))
  df$soil_temp_C <- suppressWarnings(as.numeric(df$soil_temp_C))
  stopifnot(all(df$air_temp_C >= -25 & df$air_temp_C <= 45, na.rm = TRUE),
            all(df$soil_temp_C >= -5 & df$soil_temp_C <= 35, na.rm = TRUE))
  df
}

# Species levels: a taxonomic ladder. A species with >= 10 records keeps its own
# level; below that it falls back to its genus (GEN_<genus>) if the genus was
# measured, and only otherwise to SPECIES_OTHER (unidentified: missing, "unknown*",
# or a single token). Chosen by repeated 5-fold CV with rare-species rows scored
# separately (10_rf_species_pooling.R): the ladder has the lowest rare-species RMSE
# and bias and ties for best overall. A species with no records cannot be predicted
# well by any scheme (09_rf_species_fallback_loso.R). The prediction side applies
# the same ladder through code/lib/species_levels.R.
species_ladder <- function(species) {
  sp_raw <- as.character(species)
  unident <- is.na(sp_raw) | grepl("^unknown", sp_raw, ignore.case = TRUE) | !grepl(" ", trimws(sp_raw))
  sp_raw[unident] <- "SPECIES_OTHER"
  gen <- ifelse(unident, "SPECIES_OTHER", sub(" .*", "", sp_raw))
  cnt <- table(sp_raw); big <- names(cnt)[cnt >= 10]
  gen_seen <- unique(gen[!unident])
  out <- ifelse(sp_raw %in% big, sp_raw,
         ifelse(gen %in% gen_seen & !unident, paste0("GEN_", gen), "SPECIES_OTHER"))
  cat(sprintf("  Species levels: %d own + %d genus + SPECIES_OTHER (%d unidentified)\n",
              sum(!grepl("^GEN_|SPECIES_OTHER", unique(out))), sum(grepl("^GEN_", unique(out))), sum(unident)))
  out
}

tree_train <- check_temps(tree_combined)
cat(sprintf("  tree temperatures at the measurement: air %d of %d rows (%d distinct), soil %d (%d distinct)\n",
            sum(!is.na(tree_train$air_temp_C)), nrow(tree_train), length(unique(na.omit(tree_train$air_temp_C))),
            sum(!is.na(tree_train$soil_temp_C)), length(unique(na.omit(tree_train$soil_temp_C)))))
tree_train$dbh_m <- as.numeric(tree_train$dbh_m)
tree_train$soil_moisture_at_tree <- as.numeric(tree_train$soil_moisture_abs)
tree_train$species_clean <- species_ladder(tree_train$species)
tree_train$y_flux <- as.numeric(tree_train$stem_flux_corrected)

soil_train <- check_temps(SOIL_YEAR)
soil_train$soil_moisture_at_site <- as.numeric(soil_train$soil_moisture_abs)
soil_train$y_flux <- as.numeric(soil_train$soil_flux_nmol_m2_s)

# One fixed row order: a seeded random forest still depends on the order of its rows.
# radix sorts text in C-locale order, so the order does not depend on the locale.
tree_train <- tree_train[order(as.character(tree_train$tree_id), tree_train$Date,
                               tree_train$measurement_height_cm, tree_train$chamber_type,
                               tree_train$stem_flux_nmol_m2_s, method = "radix"), ]
soil_train <- soil_train[order(as.character(soil_train$site_id), soil_train$Date,
                               soil_train$soil_flux_nmol_m2_s, method = "radix"), ]
rownames(tree_train) <- NULL; rownames(soil_train) <- NULL

# Training set = flux_stem.in_rf_training (code/tables/10_flux_stem.R, carried through
# 01): passed the closure screen, not a dead stem (the inventory is live stems only;
# code/lib/dead_stems.R), not the 2021 root-crown deployment, not the felled oak.
# Excluded rows remain in every descriptive analysis.
stopifnot(!anyNA(tree_train$measurement_id), !anyNA(tree_train$in_rf_training),
          !anyNA(tree_train$y_flux), !anyNA(tree_train$species_clean))
tree_train$dead_stem <- tree_train$measurement_id %in% dead_stem_ids()
stopifnot(!any(tree_train$dead_stem & tree_train$in_rf_training))
cat("  Deployments outside the training set (flux_stem.in_rf_training FALSE):", sum(!tree_train$in_rf_training),
    "of which dead stems", sum(tree_train$dead_stem), "\n")
# every trainable row, dead stems included and flagged (05_audit_training_population.R)
tree_train_candidates <- tree_train
save(tree_train_candidates, file = "../../outputs/models/TRAINING_CANDIDATES.RData")

tree_train_complete <- tree_train[tree_train$in_rf_training, ]
cat("  Training rows:", nrow(tree_train_complete), "of", nrow(tree_train), "\n")

# =============================================================================
# 4. STEM FOREST
# =============================================================================
# Predictors (each beat a random reference column in permutation importance over 20
# fits when chosen): DBH, species (a native factor), soil moisture, soil temperature,
# and measurement height. Height matters because the 2021 campaign measured 50, 125
# and 200 cm on the same trees, whose rows are otherwise identical in every
# predictor; monthly and 2023 rows (and the 75 cm and root-crown rows, which carry no
# height) are 125 cm. Rejected: within-species standardised DBH (standardised within
# different populations for fitting and for the inventory), elevation, local basal
# area and distance to river (static per tree), month, air temperature (the totals
# then depend on the diurnal basis of the scaling drivers).
# Importance is read as permutation importance elsewhere (06_rf_model_diagnostics.R,
# figS24); the stored forests carry impurity importance.
X_tree <- data.frame(
  species               = factor(tree_train_complete$species_clean),
  dbh_m                 = tree_train_complete$dbh_m,
  soil_moisture_at_tree = tree_train_complete$soil_moisture_at_tree,
  soil_temp_C           = tree_train_complete$soil_temp_C,
  height_cm             = ifelse(is.na(tree_train_complete$measurement_height_cm), 125,
                                 tree_train_complete$measurement_height_cm))
no_variance <- function(X) sapply(X, function(x) {
  if (is.factor(x)) return(length(unique(x[!is.na(x)])) < 2)
  all(is.na(x)) || isTRUE(var(x, na.rm = TRUE) == 0)
})
stopifnot(!any(no_variance(X_tree)))
cat("  Tree feature matrix:", nrow(X_tree), "obs x", ncol(X_tree), "features\n")

TreeRF <- ranger(
  x = X_tree, y = tree_train_complete$y_flux,
  num.trees = RF_NUM_TREES, min.node.size = 5, mtry = floor(sqrt(ncol(X_tree))),
  importance = "impurity", num.threads = RF_THREADS, oob.error = TRUE, seed = 42)
cat("\nTreeRF trained:\n")
cat("  OOB R²:", round(TreeRF$r.squared, 3), "\n")
cat("  OOB RMSE:", round(sqrt(TreeRF$prediction.error), 4), "\n")

# =============================================================================
# 5. SOIL FOREST
# =============================================================================
# Predictors: soil moisture (diffusion-limited uptake), soil temperature, month.
# Excluded: distance to river, elevation and local basal area. Each is static per
# collar (28 values over 266 rows), i.e. collar identity: with them OOB R2 rises to
# 0.659 while grouped-CV R2 is 0.451.
X_soil <- soil_train[, c("soil_moisture_at_site", "soil_temp_C", "month"), drop = FALSE]
stopifnot(!any(no_variance(X_soil)))
soil_train_complete <- soil_train
cat("  Soil feature matrix:", nrow(X_soil), "obs x", ncol(X_soil), "features\n")

SoilRF <- ranger(
  x = as.data.frame(X_soil), y = soil_train_complete$y_flux,
  num.trees = RF_NUM_TREES, min.node.size = 5, mtry = floor(sqrt(ncol(X_soil))),
  importance = "impurity", num.threads = RF_THREADS, oob.error = TRUE, seed = 42)
cat("\nSoilRF trained:\n")
cat("  OOB R²:", round(SoilRF$r.squared, 3), "\n")
cat("  OOB RMSE:", round(sqrt(SoilRF$prediction.error), 4), "\n\n")

# =============================================================================
# 6. SAVE
# =============================================================================
# pred_flux: the out-of-bag prediction of each training row (nmol m-2 s-1)
tree_train_complete$pred_flux <- TreeRF$predictions
soil_train_complete$pred_flux <- SoilRF$predictions
save(tree_train_complete, X_tree, X_soil, soil_train_complete,
     file = "../../outputs/models/TRAINING_DATA.RData")
save(TreeRF, SoilRF, file = "../../outputs/models/RF_MODELS.RData")
cat("Saved outputs/models/RF_MODELS.RData, TRAINING_DATA.RData, TRAINING_CANDIDATES.RData\n")
cat("\n=== 02_rf_models.R COMPLETE ===\n")
