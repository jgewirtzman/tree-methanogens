# ==============================================================================
# Compile Clean Datasets for Zenodo
# ==============================================================================
# Purpose: Creates analysis-ready compiled datasets from processed pipeline
#   outputs for public archival on Zenodo. Each dataset is simplified to
#   essential columns with clear naming and documented units.
#
# Outputs: data/compiled/*.csv
# ==============================================================================

library(tidyverse)

out_dir <- "data/compiled"
dir.create(out_dir, showWarnings = FALSE)

cat("============================================================\n")
cat("Compiling clean datasets for Zenodo\n")
cat("============================================================\n\n")

# ============================================================================
# 1. SEMIRIGID CHAMBER FLUX — Tree stems (2020–2021)
# ============================================================================
cat("--- 1a. Semirigid tree stem flux ---\n")
tree_flux <- read_csv("data/processed/flux/semirigid_tree_final_complete_dataset.csv",
                       show_col_types = FALSE)

tree_flux_clean <- tree_flux %>%
  transmute(
    unique_id        = UniqueID,
    date             = as.Date(Date),
    tree_tag         = as.character(`Plot Tag`),
    landscape_position = `Plot Letter`,
    chamber_id       = Chamber_ID,
    chamber_area_m2  = Area,
    chamber_vol_L    = Vtot,
    chamber_temp_C   = Tcham,
    CH4_best_flux_nmol_m2_s = CH4_best.flux.y,
    CH4_model        = CH4_model.y,
    CH4_quality      = CH4_quality.check.y,
    CH4_LM_r2        = CH4_LM.r2.y,
    CH4_LM_pval      = CH4_LM.p.val.y,
    CO2_best_flux_nmol_m2_s = CO2_best.flux,
    CO2_model        = CO2_model,
    CO2_quality      = CO2_quality.check,
    measurement_type = "tree_stem",
    notes            = Notes
  )

cat(sprintf("  %d tree stem flux measurements\n", nrow(tree_flux_clean)))

# ============================================================================
# 1b. SEMIRIGID CHAMBER FLUX — Soil (2020–2021)
# ============================================================================
cat("--- 1b. Semirigid soil flux ---\n")
soil_flux <- read_csv("data/processed/flux/semirigid_tree_final_complete_dataset_soil_CORRECTED.csv",
                       show_col_types = FALSE)

soil_flux_clean <- soil_flux %>%
  transmute(
    unique_id        = UniqueID,
    date             = as.Date(Date),
    tree_tag         = as.character(`Plot Tag`),
    landscape_position = `Plot letter`,
    chamber_id       = NA_character_,
    chamber_area_m2  = Area,
    chamber_vol_L    = Vtot,
    chamber_temp_C   = Tcham,
    CH4_best_flux_nmol_m2_s = CH4_best.flux,
    CH4_model        = CH4_model,
    CH4_quality      = CH4_quality.check,
    CH4_LM_r2        = CH4_LM.r2,
    CH4_LM_pval      = CH4_LM.p.val,
    CO2_best_flux_nmol_m2_s = CO2_best.flux,
    CO2_model        = CO2_model,
    CO2_quality      = CO2_quality.check,
    measurement_type = "soil",
    notes            = Notes
  )

cat(sprintf("  %d soil flux measurements\n", nrow(soil_flux_clean)))

# Combine tree + soil
semirigid_flux <- bind_rows(tree_flux_clean, soil_flux_clean)
write_csv(semirigid_flux, file.path(out_dir, "semirigid_chamber_flux.csv"))
cat(sprintf("  -> Wrote semirigid_chamber_flux.csv (%d rows)\n\n", nrow(semirigid_flux)))


# ============================================================================
# 2. STATIC CHAMBER FLUX (2021 + 2023)
# ============================================================================
cat("--- 2. Static chamber flux ---\n")
static_flux <- read_csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv",
                         show_col_types = FALSE)

# Fix special-char column names (degree symbols, etc.)
names(static_flux)[7]  <- "DBH_cm"
names(static_flux)[8]  <- "Air_Temp_F"
names(static_flux)[9]  <- "Soil_Temp_C_raw"
names(static_flux)[10] <- "Stem_Temp_C_raw"
names(static_flux)[11] <- "VWC_1_pct"
names(static_flux)[12] <- "VWC_2_pct"
names(static_flux)[13] <- "VWC_3_pct"
names(static_flux)[18] <- "Bark_Missing"
names(static_flux)[19] <- "Visible_fungus"
names(static_flux)[24] <- "Bole_Damage"

static_flux_clean <- static_flux %>%
  transmute(
    unique_id        = UniqueID,
    date             = as.Date(date_clean),
    tree_tag         = `Tree Tag`,
    species_code     = `Species Code`,
    species          = Species,
    dbh_cm           = DBH_cm,
    air_temp_C       = air_temp_C,
    soil_temp_C      = Soil_Temp_C_raw,
    stem_temp_C      = stem_temp_C,
    vwc_mean_pct     = vwc_mean,
    chamber_id       = `Chamber ID`,
    chamber_area_m2  = Area,
    chamber_vol_L    = Vtot,
    chamber_temp_C   = Tcham,
    obs_length_s     = obs_length_seconds,
    CH4_best_flux_nmol_m2_s = CH4_best.flux,
    CH4_model        = CH4_model,
    CH4_quality      = CH4_quality.check,
    CH4_LM_r2        = CH4_LM.r2,
    CH4_LM_pval      = CH4_LM.p.val,
    CO2_best_flux_nmol_m2_s = CO2_best.flux,
    CO2_model        = CO2_model,
    CO2_quality      = CO2_quality.check,
    bark_missing     = Bark_Missing,
    visible_fungus   = Visible_fungus,
    bole_damage      = Bole_Damage,
    notes            = Notes
  )

write_csv(static_flux_clean, file.path(out_dir, "static_chamber_flux.csv"))
cat(sprintf("  -> Wrote static_chamber_flux.csv (%d rows)\n\n", nrow(static_flux_clean)))


# ============================================================================
# 3. ddPCR GENE ABUNDANCES
# ============================================================================
cat("--- 3. ddPCR gene abundances ---\n")
ddpcr <- read_csv("data/processed/molecular/processed_ddpcr_data.csv",
                   show_col_types = FALSE)

# Fix special-char column names (µ symbol)
names(ddpcr)[26] <- "Conc_copies_per_uL"
names(ddpcr)[35] <- "Copies_20uL_Well"

ddpcr_clean <- ddpcr %>%
  transmute(
    sample_id        = Inner.Core.Sample.ID,
    plate_id         = plate_identifier,
    species          = species,
    material         = material,
    core_type        = core_type,
    target_gene      = Target,
    concentration_copies_per_uL = Conc_copies_per_uL,
    copies_per_20uL_well = Copies_20uL_Well,
    accepted_droplets = Accepted.Droplets,
    positives        = Positives,
    negatives        = Negatives,
    status           = Status,
    analysis_type    = analysis_type,
    sample_mass_mg   = Sample.Mass.Added.to.Tube..mg.,
    extraction_plate = Extraction.Plate.ID
  )

write_csv(ddpcr_clean, file.path(out_dir, "ddpcr_gene_abundances.csv"))
cat(sprintf("  -> Wrote ddpcr_gene_abundances.csv (%d rows)\n\n", nrow(ddpcr_clean)))


# ============================================================================
# 4. 16S COMMUNITY COMPOSITION
# ============================================================================
cat("--- 4. 16S community composition ---\n")
# OTU table is tab-delimited with samples as columns
otu <- read_tsv("data/raw/16s/OTU_table.txt", show_col_types = FALSE)
tax <- read_tsv("data/raw/16s/taxonomy_table.txt", show_col_types = FALSE)
meta_16s <- read_csv("data/raw/16s/16s_w_metadata.csv", show_col_types = FALSE)

# Write taxonomy key
tax_clean <- tax %>%
  rename(feature_id = `Feature ID`, taxonomy = Taxon) %>%
  separate(taxonomy,
           into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = "; ", fill = "right", remove = FALSE)

write_csv(tax_clean, file.path(out_dir, "taxonomy_key_16S.csv"))
cat(sprintf("  -> Wrote taxonomy_key_16S.csv (%d ASVs)\n", nrow(tax_clean)))

# Write sample metadata
write_csv(meta_16s, file.path(out_dir, "sample_metadata_16S.csv"))
cat(sprintf("  -> Wrote sample_metadata_16S.csv (%d samples)\n", nrow(meta_16s)))

# OTU table: too large to reshape, write as-is (samples x ASVs matrix)
write_csv(otu, file.path(out_dir, "otu_table_16S.csv"))
cat(sprintf("  -> Wrote otu_table_16S.csv (%d ASVs x %d samples)\n\n",
            nrow(otu), ncol(otu) - 1))


# ============================================================================
# 5. PICRUSt PATHWAY ASSOCIATIONS
# ============================================================================
cat("--- 5. PICRUSt pathway associations ---\n")
picrust <- read_csv("data/processed/molecular/picrust/pathway_associations_combined.csv",
                     show_col_types = FALSE)

# Also read pmoA-specific associations
picrust_pmoa <- read_csv("data/processed/molecular/picrust/pathway_associations_pmoa.csv",
                          show_col_types = FALSE) %>%
  mutate(target_gene = "pmoA")

picrust_combined <- picrust %>%
  rename(
    p_value_all    = p.all,
    t_stat_all     = t.all,
    FDR_all        = FDR.all,
    p_value_no_mcra = p.no_mcra,
    t_stat_no_mcra = t.no_mcra,
    FDR_no_mcra    = FDR.no_mcra
  )

write_csv(picrust_combined, file.path(out_dir, "picrust_pathway_associations.csv"))
cat(sprintf("  -> Wrote picrust_pathway_associations.csv (%d pathways)\n\n", nrow(picrust_combined)))


# ============================================================================
# 6. TREE PROPERTIES (Master integrated dataset)
# ============================================================================
cat("--- 6. Tree properties (master dataset) ---\n")
merged <- read_csv("data/processed/integrated/merged_tree_dataset_final.csv",
                    show_col_types = FALSE)

# This is already a clean integrated dataset — write with improved names
merged_clean <- merged %>%
  rename(
    landscape_position = plot
  )

write_csv(merged_clean, file.path(out_dir, "tree_properties.csv"))
cat(sprintf("  -> Wrote tree_properties.csv (%d trees x %d variables)\n\n",
            nrow(merged_clean), ncol(merged_clean)))


# ============================================================================
# 7. FOREST INVENTORY (ForestGEO)
# ============================================================================
cat("--- 7. Forest inventory ---\n")
inventory <- read_csv("data/raw/inventory/ForestGEO_data2021UPDATE_6_21_DW_2019.csv",
                       show_col_types = FALSE)

# Clean column names (handles BOM and spaces)
inv_names <- tolower(trimws(names(inventory)))
inv_names <- gsub("\uFEFF", "", inv_names)  # remove BOM
inv_names <- gsub(" ", "_", inv_names)
names(inventory) <- inv_names

inventory_clean <- inventory %>%
  rename(
    dbh_mm       = dbh,
    x_m          = px,
    y_m          = py
  )

write_csv(inventory_clean, file.path(out_dir, "forest_inventory.csv"))
cat(sprintf("  -> Wrote forest_inventory.csv (%d stems)\n\n", nrow(inventory_clean)))


# ============================================================================
# 8. ENVIRONMENTAL TIME SERIES
# ============================================================================
cat("--- 8. Environmental time series ---\n")
weather <- read_csv("data/raw/weather/ymf_clean_sorted.csv",
                     show_col_types = FALSE)

weather_clean <- weather %>%
  transmute(
    timestamp        = TIMESTAMP,
    year             = year,
    month            = month,
    day              = day,
    hour             = hour,
    air_temp_C       = Tair,
    air_temp_avg_C   = Tair_Avg,
    air_temp_max_C   = Tair_Max,
    air_temp_min_C   = Tair_Min,
    soil_temp_avg_C  = Tsoil_Avg,
    relative_humidity_pct = RH,
    vwc              = VWC,
    vwc_avg          = VWC_Avg,
    precip_mm        = Rain_mm_Tot,
    wind_speed_avg_ms = WindSpeed_ms_Avg,
    wind_speed_max_ms = WindSpeed_ms_Max,
    wind_dir_deg     = WindDir,
    solar_rad_avg_kW = Solar_Rad_kW_Avg,
    solar_rad_total_MJ = Solar_Rad_Tot_MJ_Tot,
    dewpoint_C       = TdewPointC,
    ET_ref           = ETos
  )

write_csv(weather_clean, file.path(out_dir, "environmental_timeseries.csv"))
cat(sprintf("  -> Wrote environmental_timeseries.csv (%d records)\n\n", nrow(weather_clean)))


# ============================================================================
# 9. BLACK OAK FELLED TREE EXPERIMENT
# ============================================================================
cat("--- 9. Black oak experiment ---\n")
black_oak <- read_csv("data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv",
                       show_col_types = FALSE)

black_oak_clean <- black_oak %>%
  transmute(
    unique_id        = UniqueID,
    height_m         = Height_m,
    chamber          = Chamber,
    stem_temp_C      = Stem_Temp_C,
    stem_diam_mm     = Stem_Diam_mm,
    air_temp_C       = Air_Temp_C,
    obs_length_s     = obs_length_sec,
    CH4_best_flux_nmol_m2_s = CH4_best.flux,
    CH4_model        = CH4_model,
    CH4_quality      = CH4_quality.check,
    CH4_LM_flux      = CH4_LM.flux,
    CH4_LM_r2        = CH4_LM.r2,
    CO2_best_flux_nmol_m2_s = CO2_best.flux,
    CO2_model        = CO2_model,
    CO2_quality      = CO2_quality.check,
    notes            = Notes
  )

write_csv(black_oak_clean, file.path(out_dir, "black_oak_experiment.csv"))
cat(sprintf("  -> Wrote black_oak_experiment.csv (%d rows)\n\n", nrow(black_oak_clean)))


# ============================================================================
# 10. METHANOTROPH DEFINITIONS (curated taxonomy lookup)
# ============================================================================
cat("--- 10. Methanotroph definitions ---\n")
methano_def <- read_csv("data/processed/molecular/methanotroph_definitions.csv",
                         show_col_types = FALSE)
write_csv(methano_def, file.path(out_dir, "methanotroph_definitions.csv"))
cat(sprintf("  -> Wrote methanotroph_definitions.csv (%d taxa)\n\n", nrow(methano_def)))


# ============================================================================
# SUMMARY
# ============================================================================
cat("============================================================\n")
compiled_files <- list.files(out_dir, pattern = "\\.csv$")
cat(sprintf("Done! %d compiled datasets written to %s/\n", length(compiled_files), out_dir))
for (f in compiled_files) {
  size <- file.size(file.path(out_dir, f))
  size_str <- if (size > 1e6) sprintf("%.1f MB", size / 1e6) else sprintf("%.0f KB", size / 1e3)
  cat(sprintf("  %-45s %s\n", f, size_str))
}
cat("============================================================\n")
