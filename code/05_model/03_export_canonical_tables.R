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
# 02_campaign_counts.R reads five, 04_variance_partition.R reads two and
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
#   outputs/data/flux_measurements_tree.csv   every stem deployment with a CH4 flux;
#                                             in_rf_training marks the 1,130 the model used
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
# soil_flux_umol_m2_s hold nmol m-2 s-1, not umol; 20_manuscript_statistics.R:1447
# records this for the equivalent Phi_* columns. Verified against the budget:
# the mean stem flux read as nmol gives 6.0 mg CH4 m-2 ground yr-1 against the
# canonical 4.912, while reading it as umol gives 6,005 -- three orders out.
# Internally the misnomer is survivable because the comment sits beside the code.
# In a published dataset it is not: nothing travels with the column but its name.
names(tree_out)[names(tree_out) == "stem_flux_umol_m2_s"] <- "stem_flux_nmol_m2_s"
names(soil_out)[names(soil_out) == "soil_flux_umol_m2_s"] <- "soil_flux_nmol_m2_s"

# ---- every stem measurement, flagged -- not only the modelled ones -------------
# 2026-09-30. This table used to hold ONLY the model's training rows, under a name
# that said "measurements". 61 of the 1,191 stem deployments with a CH4 flux never
# reach training, because the covariate assembly in 01_load_and_prep_data.R links
# rows through tree tags and the inventory:
#   2023 survey   22  tags that are not numbers ("Ash1".."Ash16", "untagged")
#   monthly       38  recovered untagged / dead-snag stems ("UNTAG_...")
#   2021 heights   1  a Kalmia root-crown deployment (no stem height)
# The rebuilt Figure 3 read this file as if it were complete; losing 17 low-emitting
# ash trees made white ash look like the top emitter. The model is frozen, so the
# training set is not changed here. Instead every deployment is listed, with
# in_rf_training and exclusion_reason, and the table refuses to write unless each
# campaign's rows account exactly for its source file.
#   Refit the model:            filter(in_rf_training)
#   Describe the measurements:  all rows (covariates are NA where none were measured)
suppressMessages(library(dplyr))
source("code/lib/dead_stems.R"); DEAD <- dead_stem_flux()
fd   <- "data/processed/flux"
key  <- function(x) round(as.numeric(x), 8)
LATIN <- setNames(as.character(tree_out$species), tree_out$species_code)
LATIN <- LATIN[!duplicated(names(LATIN))]

tree_out$campaign <- ifelse(tree_out$chamber_type == "semirigid", "monthly_2020_2021",
                     ifelse(tree_out$year == 2023, "2023_cross_species", "2021_multiheight"))
tree_out$in_rf_training   <- TRUE
tree_out$exclusion_reason <- NA_character_

src <- list(
  monthly_2020_2021 = read.csv(file.path(fd, "semirigid_tree_final_complete_dataset_with_untagged.csv"),
                               check.names = FALSE) %>%
    filter(!is.na(CH4_best.flux.x)) %>%
    transmute(tree_id = Plot.Tag, flux = CH4_best.flux.x, Date = as.character(Date), height = NA_real_,
              species_code = sub("^UNTAG_[A-Za-z]+_([A-Z]{4}).*$", "\\1", Plot.Tag),
              dbh_m = NA_real_, air_temp_C = suppressWarnings(as.numeric(Tcham)),
              soil_temp_C = NA_real_, soil_moisture_abs = NA_real_, chamber_type = "semirigid"),
  `2021_multiheight` = read.csv(file.path(fd, "tree_flux_2021_multiheight.csv"), check.names = FALSE) %>%
    filter(!is.na(CH4_best.flux)) %>%
    transmute(tree_id = as.character(tree_id), flux = CH4_best.flux,   # date is the UniqueID prefix
              Date = as.character(as.Date(substr(UniqueID, 1, 8), "%Y%m%d")),
              height = suppressWarnings(as.numeric(measurement_height)), species_code = species,
              dbh_m = NA_real_, air_temp_C = suppressWarnings(as.numeric(Tcham)),
              soil_temp_C = NA_real_, soil_moisture_abs = NA_real_, chamber_type = "rigid"),
  `2023_cross_species` = (function(r) {
    pick <- function(p) suppressWarnings(as.numeric(r[[grep(p, names(r), ignore.case = TRUE, value = TRUE)[1]]]))
    data.frame(tree_id = as.character(r$`Tree Tag`), flux = r$CH4_best.flux, Date = as.character(r$Date),
               height = NA_real_, species_code = r$`Species Code`, dbh_m = pick("^DBH") / 100,
               air_temp_C = pick("^air_temp_C$"), soil_temp_C = pick("^Soil.?Temp"),
               soil_moisture_abs = pick("^vwc_mean$") / 100, chamber_type = "rigid") %>%
      filter(!is.na(flux))
  })(read.csv(file.path(fd, "tree_flux_2023_cross_species.csv"), check.names = FALSE))
)

reason <- c(monthly_2020_2021  = "untagged or dead-snag stem: no tag link for the model's covariates",
            `2021_multiheight` = "root-crown deployment: no stem height",
            `2023_cross_species` = "non-numeric tree tag: no link for the model's covariates")

extra <- bind_rows(lapply(names(src), function(cmp) {
  s <- src[[cmp]]; trained <- key(tree_out$stem_flux_nmol_m2_s[tree_out$campaign == cmp])
  stopifnot(!anyDuplicated(key(s$flux)))            # flux values identify deployments within a campaign
  miss <- s[!key(s$flux) %in% trained, ]
  # ACCOUNTING: every source deployment is either in training or listed here, and
  # every training row traces back to a source deployment. Anything else stops.
  if (nrow(s) != sum(tree_out$campaign == cmp) + nrow(miss) || !all(trained %in% key(s$flux)))
    stop(sprintf("%s: %d source deployments, %d in training, %d unmatched -- rows unaccounted for",
                 cmp, nrow(s), sum(tree_out$campaign == cmp), nrow(miss)), call. = FALSE)
  cat(sprintf("  %-20s %4d deployments = %4d in training + %3d excluded\n", cmp, nrow(s),
              sum(tree_out$campaign == cmp), nrow(miss)))
  if (!nrow(miss)) return(NULL)
  d <- as.Date(substr(miss$Date, 1, 10), tryFormats = c("%Y-%m-%d", "%m/%d/%y", "%m/%d/%Y"))
  why <- ifelse(is_dead_stem_flux(miss$flux, DEAD), "dead stem: the inventory the model is applied to is live stems only",
                reason[[cmp]])
  data.frame(tree_id = paste0(sub("_.*", "", cmp), "_", gsub("[^A-Za-z0-9]", "", miss$tree_id)),
             species = unname(LATIN[miss$species_code]), species_code = miss$species_code,
             Date = as.character(d), month = as.integer(format(d, "%m")), year = as.integer(format(d, "%Y")),
             stem_flux_nmol_m2_s = miss$flux, air_temp_C = miss$air_temp_C, soil_temp_C = miss$soil_temp_C,
             soil_moisture_abs = miss$soil_moisture_abs, dbh_m = miss$dbh_m, chamber_type = miss$chamber_type,
             measurement_height_cm = miss$height, campaign = cmp, in_rf_training = FALSE,
             dead_stem = is_dead_stem_flux(miss$flux, DEAD),   # was left NA for excluded rows
             exclusion_reason = why, stringsAsFactors = FALSE)
}))
# Untagged stems share labels ("untagged" x5 in 2023); number them so each deployment
# keeps its own id. Monthly UNTAG_ labels are per-stem and stay as they are.
dup <- extra$tree_id %in% extra$tree_id[duplicated(extra$tree_id)] & extra$campaign == "2023_cross_species"
extra$tree_id[dup] <- paste0(extra$tree_id[dup], "_", ave(seq_along(extra$tree_id[dup]), extra$tree_id[dup], FUN = seq_along))

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
