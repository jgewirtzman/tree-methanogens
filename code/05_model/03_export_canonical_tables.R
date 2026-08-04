#!/usr/bin/env Rscript
# ==============================================================================
# 03_export_canonical_tables.R -- write the model's merged tables as CSV
#
# The flux data reaches the model through a merge built in 01_load_and_prep_data.R:
# the 2021 multi-height campaign is assembled from five per-height pieces
# (125 cm, 50 cm, 200 cm, Kalmia, other 75 cm), bound to the 2023 cross-species
# campaign and the 2020-2021 semi-rigid monthly survey, then joined to
# environmental drivers. The result is 1,130 stem measurements from 478 trees
# across 2020, 2021 and 2023.
#
# That merged table existed only inside outputs/models/TRAINING_DATA.RData, so
# every other consumer rebuilt its own version from the component files:
# rev_stat_campaign_counts.R reads five, 04_variance_partition.R reads two and
# covers 2023 only, and the Zenodo archive was compiled from a different subset
# again -- it carried the 2020-2021 semi-rigid and 2023 campaigns but not the
# 2021 rigid one, so 328 of the 1,130 measurements the model uses had no
# archived source and the model could not be refit from the archive.
#
# This script does not merge anything. It loads what the model was actually
# fitted on and writes it out, so there is one measurement-level flux table with
# one producer.
#
# Run after 02_rf_models.R, from the repository root:
#   Rscript code/05_model/03_export_canonical_tables.R
#
# Outputs:
#   outputs/data/flux_measurements_tree.csv   1,130 stem measurements
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

# Prediction columns are model output, not measurements; they are reproduced by
# refitting and would otherwise read as observed data in the archive.
drop_pred <- function(d) d[, !names(d) %in%
  c("pred_asinh", "pred_flux", "pred_flux_nmol", "y_asinh", "obs_flux_nmol"), drop = FALSE]

tree_out <- drop_pred(tree)
soil_out <- drop_pred(soil)

# UNIT MISNOMER, corrected at the archive boundary. stem_flux_umol_m2_s and
# soil_flux_umol_m2_s hold nmol m-2 s-1, not umol; manuscript_statistics.R:1447
# records this for the equivalent Phi_* columns. Verified against the budget:
# the mean stem flux read as nmol gives 6.0 mg CH4 m-2 ground yr-1 against the
# canonical 4.912, while reading it as umol gives 6,005 -- three orders out.
# Internally the misnomer is survivable because the comment sits beside the code.
# In a published dataset it is not: nothing travels with the column but its name.
names(tree_out)[names(tree_out) == "stem_flux_umol_m2_s"] <- "stem_flux_nmol_m2_s"
names(soil_out)[names(soil_out) == "soil_flux_umol_m2_s"] <- "soil_flux_nmol_m2_s"

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
