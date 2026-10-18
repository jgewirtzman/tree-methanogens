#!/usr/bin/env Rscript
# ==============================================================================
# 03_export_canonical_tables.R -- write the model's merged tables as CSV
#
# Stage B. Writes out the rows the forests were fitted on, plus every other stem
# deployment with its exclusion reason, so there is one measurement-level flux table
# with one producer (read by fig03 and 11_variance-partition.R). It merges nothing:
# the training rows are 01_load_and_prep_data.R's (flux_stem deployments with trees.csv
# species, DBH and position and the conditions at the deployment).
#
# Run after 02_rf_models.R, from the repository root:
#   Rscript code/05_model/03_export_canonical_tables.R
#
# Outputs:
#   outputs/data/flux_measurements_tree.csv   every stem deployment with a CH4 flux;
#                                             in_rf_training marks the rows the model used
#   outputs/data/flux_measurements_soil.csv     266 soil measurements
# ==============================================================================
source("code/lib/outputs.R")

NEED <- "outputs/models/TRAINING_DATA.RData"
if (!file.exists(NEED))
  stop("missing ", NEED, "\nBuild it first (from code/05_model/):\n",
       "  Rscript 01_load_and_prep_data.R && Rscript 02_rf_models.R", call. = FALSE)

e <- new.env(); load(NEED, envir = e)
tree <- get("tree_train_complete", envir = e)
soil <- get("soil_train_complete", envir = e)

# Model columns are not measurements: the response (y_flux, the chamber-harmonised
# flux) and the out-of-bag prediction (pred_flux) are reproduced by refitting.
drop_pred <- function(d) d[, !names(d) %in% c("pred_flux", "y_flux"), drop = FALSE]
tree_out <- drop_pred(tree)
soil_out <- drop_pred(soil)

# ---- every stem measurement, flagged -- not only the modelled ones -------------
# Every stem deployment of the three surveys in data/tables/flux_stem.csv is listed,
# matched to the training rows by measurement_id: the training rows as the model
# saw them, the others (dead stems; the root-crown deployment of the Kalmia that
# also has a 75 cm one) from flux_stem and trees.csv (DBH: recorded with the
# deployment, else the tree's), with in_rf_training
# and exclusion_reason as flux_stem states them. The table refuses to write unless
# the two account exactly for each other.
#   Refit the model:            filter(in_rf_training)
#   Describe the measurements:  all rows (covariates are NA where none were measured)
suppressMessages(library(dplyr))
source("code/lib/soil_fill.R")      # soil values of the excluded deployments: attach_soil()
rdt <- function(t) read.csv(file.path("data/tables", paste0(t, ".csv")), colClasses = "character", na.strings = "")
FS <- rdt("flux_stem") %>% filter(campaign != "2022_felled_oak")
TR <- rdt("trees"); SP <- rdt("species")
num <- function(x) suppressWarnings(as.numeric(x))

i <- match(tree_out$measurement_id, FS$measurement_id)
stopifnot(!anyNA(i), !anyDuplicated(tree_out$measurement_id), all(FS$in_rf_training[i] == "TRUE"),
          identical(tree_out$tree_id, FS$tree_id[i]))
tree_out$campaign         <- FS$campaign[i]
stopifnot(isTRUE(all.equal(tree_out$air_temp_C, num(FS$air_temp_C[i]))))   # the measured air the model used
tree_out$in_rf_training   <- TRUE
tree_out$exclusion_reason <- NA_character_

ex <- FS[!FS$measurement_id %in% tree_out$measurement_id, ]
stopifnot(all(ex$in_rf_training == "FALSE"), !anyNA(ex$exclusion_reason),
          nrow(tree_out) + nrow(ex) == nrow(FS))
t <- match(ex$tree_id, TR$tree_id)
d <- as.Date(ex$date)
exs <- attach_soil(ex[, c("tree_id", "start_time")], "tree_id", "start_time", "tree")
extra <- data.frame(measurement_id = ex$measurement_id, tree_id = ex$tree_id,
                    species = SP$species[match(TR$species_code[t], SP$species_code)], species_code = TR$species_code[t],
                    Date = as.character(d), month = as.integer(format(d, "%m")), year = as.integer(format(d, "%Y")),
                    stem_flux_nmol_m2_s = num(ex$CH4_flux_nmol_m2_s), air_temp_C = num(ex$air_temp_C),
                    soil_temp_C = exs$soil_temp_C, soil_moisture_abs = exs$soil_vwc_pct / 100, soil_match = exs$soil_match,
                    dbh_m = coalesce(num(ex$dbh_recorded_cm), num(TR$dbh_cm[t])) / 100, chamber_type = ex$chamber_type,
                    measurement_height_cm = num(ex$height_cm), campaign = ex$campaign, in_rf_training = FALSE,
                    dead_stem = ex$dead_stem == "TRUE", exclusion_reason = ex$exclusion_reason,
                    stringsAsFactors = FALSE)
for (cmp in sort(unique(FS$campaign)))
  cat(sprintf("  %-20s %4d deployments = %4d in training + %3d excluded\n", cmp, sum(FS$campaign == cmp),
              sum(tree_out$campaign == cmp), sum(extra$campaign == cmp)))

tree_out$Date <- as.character(tree_out$Date)
tree_out <- bind_rows(tree_out, extra)
cat(sprintf("  all stem deployments: %d (%d in training, %d excluded)\n",
            nrow(tree_out), sum(tree_out$in_rf_training), sum(!tree_out$in_rf_training)))

write.csv(tree_out, out_path("flux_measurements_tree.csv"), row.names = FALSE)
write.csv(soil_out, out_path("flux_measurements_soil.csv"), row.names = FALSE)

cat(sprintf("flux_measurements_tree.csv  %d measurements, %d trees, %d cols\n",
            nrow(tree_out), length(unique(tree_out$tree_id)), ncol(tree_out)))
cat("  years:      ", paste(sprintf("%s=%d", names(table(tree_out$year)),
                                    table(tree_out$year)), collapse = "  "), "\n")
cat("  chambers:   ", paste(sprintf("%s=%d", names(table(tree_out$chamber_type)),
                                    table(tree_out$chamber_type)), collapse = "  "), "\n")
cat("  heights (cm):", paste(range(tree_out$measurement_height_cm, na.rm = TRUE),
                             collapse = " - "), "\n")
cat(sprintf("flux_measurements_soil.csv  %d measurements, %d cols\n",
            nrow(soil_out), ncol(soil_out)))
