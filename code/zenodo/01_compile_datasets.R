# ==============================================================================
# 01_compile_datasets.R -- the archived datasets in data/compiled/
# ------------------------------------------------------------------------------
# Runs at the end of stage B (after cleaning, the model fit and the measurement-table
# export), so every later stage reads the same archive a reader downloads. Results
# computed later are archived separately by 02_compile_results.R.
#
# One table per kind of observation, each at a stated grain, each row kept and
# flagged rather than deleted (code/DATA_HYGIENE.md lists every rule):
#   flux_stem.csv          one stem chamber deployment (all four campaigns)
#   flux_soil.csv          one soil chamber deployment
#   internal_gas.csv       one internal-gas sample (2021-22 survey, 2024 run)
#   isotopes.csv           one Picarro run, with the whole-tree sample rule applied as flags
#   ddpcr_gene_abundances.csv, otu_table_16S.csv, taxonomy_key_16S.csv,
#   sample_metadata_16S.csv, tree_properties.csv, forest_inventory_stems.csv,
#   environmental_timeseries.csv, black_oak_experiment.csv, black_oak_its_load.csv,
#   methanotroph_definitions.csv
#
# Replaces (2026-10-01) the five overlapping flux files: semirigid_chamber_flux,
# static_chamber_flux, height_chamber_flux, flux_measurements_tree, flux_measurements_soil.
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/dead_stems.R"); source("code/lib/soil_collars.R"); source("code/lib/isotope_samples.R")
out_dir <- "data/compiled"
dir.create(out_dir, showWarnings = FALSE)
old <- c("semirigid_chamber_flux.csv", "static_chamber_flux.csv", "height_chamber_flux.csv",
         "flux_measurements_tree.csv", "flux_measurements_soil.csv", "picrust_pathway_associations.csv")
unlink(file.path(out_dir, old))                 # superseded files must not linger in the archive
fd <- "data/processed/flux"
key <- function(x) sprintf("%.8f", as.numeric(x))  # deployments are matched by campaign + flux value
wr <- function(d, f) { write_csv(d, file.path(out_dir, f)); cat(sprintf("  %-32s %6d rows x %2d cols\n", f, nrow(d), ncol(d))) }

# ---------------------------------------------------------------- flux_stem --
mon_all <- read.csv(file.path(fd, "semirigid_tree_final_complete_dataset_with_untagged.csv"), check.names = FALSE)
mon_tag <- read.csv(file.path(fd, "semirigid_tree_final_complete_dataset.csv"), check.names = FALSE)
stopifnot(isTRUE(all.equal(mon_all$CH4_best.flux.x, mon_all$CH4_best.flux.y)))
mon <- mon_all %>% filter(!is.na(CH4_best.flux.x)) %>% transmute(
  unique_id = UniqueID, campaign = "monthly_2020_2021", date = as.Date(Date),
  tree_tag = as.character(Plot.Tag), landscape_position = Plot.Letter, chamber_type = "semirigid",
  height_label = "125", chamber_id = as.character(Chamber_ID), chamber_area_m2 = Area, chamber_vol_L = Vtot,
  chamber_temp_C = Tcham, CH4_flux_nmol_m2_s = CH4_best.flux.x, CH4_model = CH4_model.x,
  CH4_quality = CH4_quality.check.x, CH4_LM_r2 = CH4_LM.r2.x, CH4_LM_pval = CH4_LM.p.val.x,
  CO2_flux_umol_m2_s = CO2_best.flux, notes = Notes,
  untagged_rescue = !UniqueID %in% mon_tag$UniqueID)
h21 <- read.csv(file.path(fd, "tree_flux_2021_multiheight.csv"), check.names = FALSE) %>%
  filter(!is.na(CH4_best.flux)) %>% transmute(
  unique_id = UniqueID, campaign = "2021_multiheight", date = as.Date(substr(start.time, 1, 10)),
  tree_tag = as.character(tree_id), landscape_position = as.character(plot), chamber_type = "rigid",
  height_label = as.character(measurement_height), chamber_id = NA_character_, chamber_area_m2 = Area,
  chamber_vol_L = Vtot, chamber_temp_C = Tcham, CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model,
  CH4_quality = CH4_quality.check, CH4_LM_r2 = CH4_LM.r2, CH4_LM_pval = CH4_LM.p.val,
  CO2_flux_umol_m2_s = CO2_best.flux, notes = NA_character_, untagged_rescue = FALSE)
y23 <- read.csv(file.path(fd, "tree_flux_2023_cross_species.csv"), check.names = FALSE)
names(y23)[grepl("^Bark", names(y23))][1] -> .bark; names(y23)[grepl("^Wounding", names(y23))][1] -> .wnd
c23 <- y23 %>% filter(!is.na(CH4_best.flux)) %>% transmute(
  unique_id = UniqueID, campaign = "2023_cross_species", date = as.Date(date_clean),
  tree_tag = as.character(`Tree Tag`), landscape_position = NA_character_, chamber_type = "rigid",
  height_label = "125", chamber_id = as.character(`Chamber ID`), chamber_area_m2 = Area, chamber_vol_L = Vtot,
  chamber_temp_C = Tcham, CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_quality = CH4_quality.check,
  CH4_LM_r2 = CH4_LM.r2, CH4_LM_pval = CH4_LM.p.val, CO2_flux_umol_m2_s = CO2_best.flux,
  notes = as.character(Notes), untagged_rescue = FALSE,
  bark_missing = as.character(.data[[.bark]]), wounding = as.character(.data[[.wnd]]))
oak <- read.csv("data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv", check.names = FALSE) %>%
  filter(!is.na(CH4_best.flux)) %>% transmute(
  unique_id = UniqueID, campaign = "2022_felled_oak", date = as.Date("2022-10-04"),
  tree_tag = "felled black oak", landscape_position = NA_character_, chamber_type = "rigid",
  height_label = paste(as.character(Height_m), "m"), chamber_id = as.character(Chamber), chamber_area_m2 = NA_real_,
  chamber_vol_L = NA_real_, chamber_temp_C = Stem_Temp_C, CH4_flux_nmol_m2_s = CH4_best.flux,
  CH4_model, CH4_quality = CH4_quality.check, CH4_LM_r2 = CH4_LM.r2, CH4_LM_pval = CH4_LM.p.val,
  CO2_flux_umol_m2_s = CO2_best.flux, notes = as.character(Notes), untagged_rescue = FALSE,
  oak_mdf = CH4_MDF)
# NOTE: felled-oak date is the felling date; the standing-tree flux was measured before it.
stem <- bind_rows(mon, h21, c23, oak)
# The hand-rescued untagged deployments carry no goFlux UniqueID; name them by stem
# label and date so every archived deployment has an identifier.
stem$unique_id <- ifelse(is.na(stem$unique_id), paste0("untagged_", stem$tree_tag, "_", stem$date), stem$unique_id)
stopifnot(!anyDuplicated(stem$unique_id))

# detection (MDF, analyzer precision) by deployment ID
DET <- read.csv("outputs/data/flux_FINAL.csv") %>% filter(type == "stem") %>%
  transmute(unique_id = UniqueID, CH4_MDF_nmol_m2_s = MDF, analyzer_sigma_ppb = sigma,
            closure_s = t_sec, detected = detected, detection_class = class)
stem <- stem %>% left_join(DET, by = "unique_id") %>%
  mutate(CH4_MDF_nmol_m2_s = coalesce(CH4_MDF_nmol_m2_s, oak_mdf)) %>% select(-oak_mdf)
# no detection class for the 45 rescued untagged stems (hand-picked windows have no
# per-record precision estimate) or the felled oak (MDF from its own goFlux run only)

# the modelled measurement table: canonical tree, species, drivers, training flag
M <- read.csv("outputs/data/flux_measurements_tree.csv") %>% transmute(
  campaign, k = key(stem_flux_nmol_m2_s), tree_id, species, species_code, dbh_m,
  measurement_height_cm, air_temp_C, soil_temp_C, soil_moisture_vwc = soil_moisture_abs,
  dead_stem, in_rf_training, exclusion_reason)
stopifnot(!anyDuplicated(M[, c("campaign", "k")]))
stem <- stem %>% mutate(k = key(CH4_flux_nmol_m2_s)) %>% left_join(M, by = c("campaign", "k")) %>% select(-k)
surveyed <- stem$campaign != "2022_felled_oak"
stopifnot(sum(surveyed) == nrow(M), all(!is.na(stem$in_rf_training[surveyed])))
stem <- stem %>% mutate(
  in_rf_training = ifelse(surveyed, in_rf_training, FALSE),
  exclusion_reason = ifelse(surveyed, exclusion_reason, "felled black oak: analyzed separately"),
  dead_stem = ifelse(surveyed, dead_stem, FALSE),
  qc_pass = TRUE, qc_reason = NA_character_)          # no stem failed the closure screen
wr(stem, "flux_stem.csv")

# ---------------------------------------------------------------- flux_soil --
QX <- read.csv("outputs/data/qc_excluded_measurements.csv", stringsAsFactors = FALSE)
SOIL_TRAIN <- read.csv("outputs/data/flux_measurements_soil.csv") %>% transmute(site_id, k = key(soil_flux_nmol_m2_s), in_rf_training = TRUE)
soil <- read.csv(file.path(fd, "semirigid_tree_final_complete_dataset_soil.csv"), check.names = FALSE) %>%
  filter(!is.na(CH4_best.flux)) %>% transmute(
  unique_id = UniqueID, campaign = "monthly_2020_2021", date = as.Date(Date),
  site_id = paste0(`Plot letter`, "_", `Plot Tag`), landscape_position = `Plot letter`,
  chamber_area_m2 = Area, chamber_vol_L = Vtot, chamber_temp_C = Tcham,
  CH4_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_quality = CH4_quality.check,
  CH4_LM_r2 = CH4_LM.r2, CH4_LM_pval = CH4_LM.p.val, CO2_flux_umol_m2_s = CO2_best.flux, notes = Notes) %>%
  left_join(read.csv("outputs/data/flux_FINAL.csv") %>% filter(type == "soil") %>%
              transmute(unique_id = UniqueID, CH4_MDF_nmol_m2_s = MDF, analyzer_sigma_ppb = sigma,
                        closure_s = t_sec, detected, detection_class = class), by = "unique_id") %>%
  mutate(k = key(CH4_flux_nmol_m2_s), qc_pass = !unique_id %in% QX$UniqueID,
         qc_reason = QX$reason[match(unique_id, QX$UniqueID)]) %>%
  left_join(SOIL_TRAIN, by = c("site_id", "k")) %>% select(-k) %>%
  mutate(in_rf_training = coalesce(in_rf_training, FALSE),
         exclusion_reason = case_when(in_rf_training ~ NA_character_,
                                      !qc_pass ~ "chamber not at ambient at closure",
                                      site_id %in% OUT_OF_STAND_COLLARS ~ "collar outside the modelled stand",
                                      TRUE ~ "no soil temperature or moisture record for the measurement"))
stopifnot(sum(soil$in_rf_training) == nrow(SOIL_TRAIN))
wr(soil, "flux_soil.csv")
cat("    soil exclusions:", paste(names(table(soil$exclusion_reason)), table(soil$exclusion_reason), collapse = "; "), "\n")

# ------------------------------------------------------------- internal_gas --
g21 <- read.csv("data/processed/internal_gas/sample_data_only.csv") %>% transmute(
  campaign = "2021_survey", sample_id = Lab.ID, tree_id = Tree.ID, species_code = Species.ID,
  tissue = "whole stem (bark to pith)", analysis_date = as.character(FID.Date),
  CH4_ppm = CH4_concentration, CO2_ppm = CO2_concentration, N2O_ppm = N2O_concentration,
  O2_ppm = O2_concentration, CH4_below_lod, CO2_below_lod, N2O_below_lod,
  CH4_above_top_standard = CH4_above_SB6)
g24 <- read.csv("data/processed/internal_gas/stem_gas_2024_calibrated.csv") %>% filter(Sample.Type == "Sample") %>%
  transmute(campaign = "2024_heartwood_sapwood", sample_id = Sample.ID, tree_id = as.character(Tree.No),
            species_code = Species, tissue = recode(Tissue, H = "heartwood", S = "sapwood"),
            analysis_date = as.character(FID_Date), CH4_ppm = CH4_calibrated_ppm, CO2_ppm = CO2_calibrated_ppm)
wr(bind_rows(g21, g24), "internal_gas.csv")

# ----------------------------------------------------------------- isotopes --
raw <- picarro_runs()
keep <- isotope_whole_tree_samples(raw)
iso <- raw %>% transmute(sample_name = SampleName, d13CH4_permil = HR_Delta_iCH4_Raw_mean,
                         CH4_ppm = HR_12CH4_dry_mean, d13CO2_permil = Delta_Raw_iCO2_mean, CO2_ppm = `12CO2_mean`) %>%
  mutate(sample_type = case_when(
           grepl("^Amb|room[[:space:]]*air", sample_name, ignore.case = TRUE) ~ "atmosphere or room air",
           grepl("^V", sample_name) ~ "incubation",
           grepl("[HS]$", sample_name) ~ "heartwood/sapwood pair",
           sample_name %in% ISO_STANDARDS ~ "standard",
           TRUE ~ "whole tree"),
         redone = sample_name %in% sub("_Redo$", "", sample_name[grepl("_Redo$", sample_name)]),
         in_whole_tree_set = sample_name %in% keep$SampleName,
         exclusion_reason = case_when(
           in_whole_tree_set ~ NA_character_,
           sample_type != "whole tree" ~ "not a whole-tree sample",
           redone ~ "replaced by its re-run",
           is.na(d13CH4_permil) ~ "no delta13C value",
           CH4_ppm < ISO_CH4_FLOOR ~ "CH4 below 1.5 ppm",
           TRUE ~ "delta13CH4 outside -115 to +25 permil")) %>%
  left_join(keep %>% transmute(sample_name = SampleName, tree_id, species), by = "sample_name")
stopifnot(sum(iso$in_whole_tree_set) == nrow(keep))
wr(iso, "isotopes.csv")

# ------------------------------------------------------------ ddPCR, 16S -----
ddpcr <- read_csv("data/processed/molecular/processed_ddpcr_data.csv", show_col_types = FALSE)
names(ddpcr)[26] <- "Conc_copies_per_uL"; names(ddpcr)[35] <- "Copies_20uL_Well"
wr(ddpcr %>% transmute(sample_id = Inner.Core.Sample.ID, plate_id = plate_identifier, species, material,
                       core_type, target_gene = Target, concentration_copies_per_uL = Conc_copies_per_uL,
                       copies_per_20uL_well = Copies_20uL_Well, accepted_droplets = Accepted.Droplets,
                       positives = Positives, negatives = Negatives, status = Status, analysis_type,
                       sample_mass_mg = Sample.Mass.Added.to.Tube..mg., extraction_plate = Extraction.Plate.ID),
   "ddpcr_gene_abundances.csv")
tax <- read_tsv("data/raw/16s/taxonomy_table.txt", show_col_types = FALSE) %>%
  rename(feature_id = `Feature ID`, taxonomy = Taxon) %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = "; ", fill = "right", remove = FALSE)
wr(tax, "taxonomy_key_16S.csv")
wr(read_csv("data/raw/16s/16s_w_metadata.csv", show_col_types = FALSE), "sample_metadata_16S.csv")
wr(read_tsv("data/raw/16s/OTU_table.txt", show_col_types = FALSE), "otu_table_16S.csv")

# ------------------------------------------------------- trees and the stand --
wr(read_csv("data/processed/integrated/merged_tree_dataset_final.csv", show_col_types = FALSE) %>%
     rename(landscape_position = plot), "tree_properties.csv")
inv <- read_csv("data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_2019.csv", show_col_types = FALSE)
names(inv) <- gsub(" ", "_", gsub("﻿", "", tolower(trimws(names(inv)))))
census <- inv %>% transmute(tag = as.character(tag), stem = as.character(stem_tag),
            census_section_id = section_id, census_quadrat = quadrat, census_sub_quadrat = sub_quadrat,
            census_status = status, census_pom_m = pom, census_codes = codes, census_notes = notes,
            latitude, longitude) %>% distinct(tag, stem, .keep_all = TRUE)
wr(read_csv("outputs/tables/inventory_stems.csv", show_col_types = FALSE) %>%
     mutate(tag = as.character(tag), stem = as.character(stem)) %>% left_join(census, by = c("tag", "stem")),
   "forest_inventory_stems.csv")

# ------------------------------------------------- environment, felled oak ---
wr(read_csv("data/raw/weather/ymf_clean_sorted.csv", show_col_types = FALSE) %>% transmute(
     timestamp = TIMESTAMP, year, month, day, hour, air_temp_C = Tair, air_temp_avg_C = Tair_Avg,
     air_temp_max_C = Tair_Max, air_temp_min_C = Tair_Min, soil_temp_avg_C = Tsoil_Avg,
     relative_humidity_pct = RH, vwc = VWC, vwc_avg = VWC_Avg, precip_mm = Rain_mm_Tot,
     wind_speed_avg_ms = WindSpeed_ms_Avg, wind_speed_max_ms = WindSpeed_ms_Max, wind_dir_deg = WindDir,
     solar_rad_avg_kW = Solar_Rad_kW_Avg, solar_rad_total_MJ = Solar_Rad_Tot_MJ_Tot,
     dewpoint_C = TdewPointC, ET_ref = ETos), "environmental_timeseries.csv")
wr(read_csv("data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv", show_col_types = FALSE) %>% transmute(
     unique_id = UniqueID, height_m = Height_m, chamber = Chamber, stem_temp_C = Stem_Temp_C,
     stem_diam_mm = Stem_Diam_mm, air_temp_C = Air_Temp_C, obs_length_s = obs_length_sec,
     CH4_best_flux_nmol_m2_s = CH4_best.flux, CH4_model, CH4_quality = CH4_quality.check,
     CH4_LM_flux = CH4_LM.flux, CH4_LM_r2 = CH4_LM.r2, CO2_best_flux_umol_m2_s = CO2_best.flux,
     CO2_model, CO2_quality = CO2_quality.check, notes = Notes), "black_oak_experiment.csv")
wr(read_csv("data/processed/molecular/black_oak/bo_its_load.csv", show_col_types = FALSE) %>%
     transmute(sample_id = `Sample ID`, ITS_copies_uL = ITS_per_ul, material = Material), "black_oak_its_load.csv")
wr(read_csv("code/lib/methanotroph_definitions.csv", show_col_types = FALSE), "methanotroph_definitions.csv")
cat(sprintf("\n%d datasets in %s/\n", length(list.files(out_dir, "\\.csv$")), out_dir))
