# Compiled Datasets for Zenodo

These datasets are analysis-ready versions of the processed data from the Yale-Myers Forest tree methane project. Each CSV is self-contained with clear column names and documented units, intended for reanalysis and meta-analysis.

**Citation**: Gewirtzman et al. (in review). DOI: [pending]

**License**: CC-BY 4.0

**Related**: Raw instrument data and full processing pipeline are available in the `data/raw/` and `data/processed/` directories. See the project `README.md` for the complete processing workflow.

---

## 1. `semirigid_chamber_flux.csv` (692 rows)

Semirigid dynamic chamber CH4 and CO2 flux measurements from tree stems and adjacent soil, Yale-Myers Forest, 2020-2021.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| unique_id | character | | Unique measurement identifier (date_tag_position_chamber) |
| date | date | YYYY-MM-DD | Measurement date |
| tree_tag | character | | ForestGEO tree tag number |
| landscape_position | character | | Landscape position: U (upland), I (intermediate), W (wetland) |
| chamber_id | character | | Chamber identifier |
| chamber_area_m2 | numeric | m2 | Chamber surface area |
| chamber_vol_L | numeric | L | Total system volume (chamber + tubing + analyzer) |
| chamber_temp_C | numeric | C | Chamber air temperature |
| CH4_best_flux_nmol_m2_s | numeric | nmol m-2 s-1 | Best-fit CH4 flux (goFlux model selection) |
| CH4_model | character | | Selected model: LM (linear) or HM (Hutchinson-Mosier) |
| CH4_quality | character | | Quality flag from goFlux |
| CH4_LM_r2 | numeric | | Linear model R-squared |
| CH4_LM_pval | numeric | | Linear model p-value |
| CO2_best_flux_nmol_m2_s | numeric | nmol m-2 s-1 | Best-fit CO2 flux |
| CO2_model | character | | Selected CO2 model |
| CO2_quality | character | | CO2 quality flag |
| measurement_type | character | | "tree_stem" or "soil" |
| notes | character | | Field notes |

---

## 2. `static_chamber_flux.csv` (406 rows)

Static chamber CH4 and CO2 flux measurements from tree stems with associated tree and environmental metadata, 2021-2023.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| unique_id | character | | Unique measurement identifier |
| date | date | YYYY-MM-DD | Measurement date |
| tree_tag | character | | ForestGEO tree tag number |
| species_code | character | | USDA species code (e.g., ACRU, QURU) |
| species | character | | Full species name |
| dbh_cm | numeric | cm | Diameter at breast height |
| air_temp_C | numeric | C | Air temperature |
| soil_temp_C | numeric | C | Soil temperature |
| stem_temp_C | numeric | C | Stem surface temperature |
| vwc_mean_pct | numeric | % | Mean volumetric water content (3 probes) |
| chamber_id | character | | Chamber identifier |
| chamber_area_m2 | numeric | m2 | Chamber surface area |
| chamber_vol_L | numeric | L | Total system volume |
| chamber_temp_C | numeric | C | Chamber temperature |
| obs_length_s | numeric | s | Observation length |
| CH4_best_flux_nmol_m2_s | numeric | nmol m-2 s-1 | Best-fit CH4 flux |
| CH4_model | character | | Selected model |
| CH4_quality | character | | Quality flag |
| CH4_LM_r2 | numeric | | Linear model R-squared |
| CH4_LM_pval | numeric | | Linear model p-value |
| CO2_best_flux_nmol_m2_s | numeric | nmol m-2 s-1 | Best-fit CO2 flux |
| CO2_model | character | | Selected CO2 model |
| CO2_quality | character | | CO2 quality flag |
| bark_missing | character | | Bark condition notes |
| visible_fungus | character | | Fungal presence observed |
| bole_damage | character | | Bole damage severity (1-3 scale) |
| notes | character | | Field notes |

---

## 3. `ddpcr_gene_abundances.csv` (4318 rows)

Digital droplet PCR (ddPCR) quantification of methane-cycling functional genes (mcrA, pmoA, mmoX) and 16S rRNA from tree wood cores and associated soil.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| sample_id | character | | Inner core sample identifier |
| plate_id | character | | ddPCR plate identifier |
| species | character | | Tree species |
| material | character | | Sample material (inner bark, outer bark, mineral soil, organic soil) |
| core_type | character | | Core type |
| target_gene | character | | Target gene: mcrA, pmoA, mmoxY, 16S_arc, 16S_bact, mcrA_probe |
| concentration_copies_per_uL | numeric | copies/uL | Gene concentration |
| copies_per_20uL_well | numeric | copies/well | Copies per 20 uL reaction well |
| accepted_droplets | integer | | Number of accepted droplets in partition |
| positives | integer | | Number of positive droplets |
| negatives | integer | | Number of negative droplets |
| status | character | | QC status |
| analysis_type | character | | Analysis classification (loose/strict threshold) |
| sample_mass_mg | numeric | mg | Sample mass added to extraction tube |
| extraction_plate | character | | DNA extraction plate identifier |

---

## 4. `otu_table_16S.csv` (79769 ASVs x 596 samples)

16S rRNA amplicon sequence variant (ASV) abundance table. Rows are ASVs, columns are samples. Values are read counts.

## 5. `taxonomy_key_16S.csv` (79769 ASVs)

Taxonomic assignments for each ASV in the OTU table.

| Column | Type | Description |
|--------|------|-------------|
| feature_id | character | ASV identifier (matches OTU table row names) |
| taxonomy | character | Full semicolon-delimited taxonomic string |
| kingdom | character | Kingdom assignment |
| phylum | character | Phylum assignment |
| class | character | Class assignment |
| order | character | Order assignment |
| family | character | Family assignment |
| genus | character | Genus assignment |
| species | character | Species assignment (often "none") |

## 6. `sample_metadata_16S.csv` (896 samples)

Sample metadata linking 16S sequencing IDs to experimental variables.

| Column | Type | Description |
|--------|------|-------------|
| Sample ID | character | Sample identifier |
| Well | character | Sequencing plate well |
| 16S_per_ul | numeric | 16S rRNA gene concentration (copies/uL) |
| core_type | character | Core type |
| seq_id_raw | character | Raw sequencing identifier |
| seq_id | character | Cleaned sequencing identifier |
| Material | character | Sample material type |

---

## 7. `picrust_pathway_associations.csv` (407 pathways)

PICRUSt2-predicted MetaCyc pathway associations with mcrA-harboring ASVs.

| Column | Type | Description |
|--------|------|-------------|
| pathway | character | MetaCyc pathway identifier |
| p_value_all | numeric | P-value (all ASVs included) |
| t_stat_all | numeric | T-statistic (all ASVs) |
| FDR_all | numeric | FDR-corrected p-value (all ASVs) |
| description | character | MetaCyc pathway description |
| mean_counts_per_million | numeric | Mean pathway abundance (CPM) |
| p_value_no_mcra | numeric | P-value (excluding mcrA-harboring ASVs) |
| t_stat_no_mcra | numeric | T-statistic (no mcrA ASVs) |
| FDR_no_mcra | numeric | FDR-corrected p-value (no mcrA ASVs) |
| mean_percent_from_mcra | numeric | Mean % of pathway contributed by mcrA ASVs |

---

## 8. `tree_properties.csv` (235 trees x 73 variables)

Master tree-level integrated dataset linking tree identity, DBH, wood properties, internal gas concentrations, gene abundances, and soil characteristics.

Key variable groups:
- **Identity**: tree_id, species_id, landscape_position (plot)
- **Size**: dbh (mm)
- **Internal gas**: CO2/CH4/N2O/O2 concentration (ppm)
- **Stem flux**: CH4/CO2 best flux at multiple heights (nmol m-2 s-1)
- **Wood moisture**: outer/inner/middle moisture (% fresh weight and dry weight basis)
- **Wood density**: outer/inner/middle density (g/cm3)
- **Gene abundances**: ddPCR copies for mcrA, pmoA, mmoX, 16S archaeal/bacterial by compartment and threshold (loose/strict)
- **Soil**: VWC (%), ORP (mV), soil temp (C), organic/mineral depth (cm)

---

## 9. `forest_inventory.csv` (7225 stems)

ForestGEO forest census data for the Yale-Myers Forest 25.6-ha plot.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| section_id | character | | Forest plot section |
| quadrat | character | | 20x20 m quadrat identifier |
| sub_quadrat | character | | 5x5 m sub-quadrat |
| tag | integer | | Tree tag number |
| stem_tag | integer | | Individual stem tag |
| species_code | character | | USDA species code |
| x_m | numeric | m | X coordinate within plot |
| y_m | numeric | m | Y coordinate within plot |
| dbh_mm | numeric | mm | Diameter at breast height |
| status | character | | Tree status (alive/dead) |
| pom | numeric | m | Point of measurement height |
| codes | character | | Condition codes |
| notes | character | | Census notes |
| secondary | character | | Secondary stem flag |
| latitude | numeric | degrees | Geographic latitude |
| longitude | numeric | degrees | Geographic longitude |

---

## 10. `environmental_timeseries.csv` (40458 records)

Meteorological tower and environmental sensor data from Yale-Myers Forest.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| timestamp | datetime | | Measurement timestamp |
| year, month, day, hour | integer | | Date/time components |
| air_temp_C | numeric | C | Instantaneous air temperature |
| air_temp_avg_C | numeric | C | Average air temperature |
| air_temp_max_C | numeric | C | Maximum air temperature |
| air_temp_min_C | numeric | C | Minimum air temperature |
| soil_temp_avg_C | numeric | C | Average soil temperature |
| relative_humidity_pct | numeric | % | Relative humidity |
| vwc | numeric | m3/m3 | Instantaneous volumetric water content |
| vwc_avg | numeric | m3/m3 | Average volumetric water content |
| precip_mm | numeric | mm | Total precipitation |
| wind_speed_avg_ms | numeric | m/s | Average wind speed |
| wind_speed_max_ms | numeric | m/s | Maximum wind speed |
| wind_dir_deg | numeric | degrees | Wind direction |
| solar_rad_avg_kW | numeric | kW/m2 | Average solar radiation |
| solar_rad_total_MJ | numeric | MJ/m2 | Total solar radiation |
| dewpoint_C | numeric | C | Dew point temperature |
| ET_ref | numeric | mm | Reference evapotranspiration |

---

## 11. `black_oak_experiment.csv` (9 rows)

CH4 and CO2 flux measurements from a felled black oak (Quercus velutina) at multiple stem heights.

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| unique_id | character | | Measurement identifier |
| height_m | numeric | m | Measurement height on stem |
| chamber | character | | Chamber identifier |
| stem_temp_C | numeric | C | Stem temperature |
| stem_diam_mm | numeric | mm | Stem diameter at measurement point |
| air_temp_C | numeric | C | Air temperature |
| obs_length_s | numeric | s | Observation length |
| CH4_best_flux_nmol_m2_s | numeric | nmol m-2 s-1 | Best-fit CH4 flux |
| CH4_model | character | | Selected model |
| CH4_quality | character | | Quality flag |
| CH4_LM_flux | numeric | nmol m-2 s-1 | Linear model CH4 flux |
| CH4_LM_r2 | numeric | | Linear model R-squared |
| CO2_best_flux_nmol_m2_s | numeric | nmol m-2 s-1 | Best-fit CO2 flux |
| CO2_model | character | | Selected CO2 model |
| CO2_quality | character | | Quality flag |
| notes | character | | Field notes |

---

## 12. `methanotroph_definitions.csv` (37 taxa)

Curated lookup table defining methanotrophic taxa used in 16S community analyses.

---

## Flux Units

All flux values are reported in **nmol m-2 s-1** (nanomoles per square meter per second). Positive values indicate emission to the atmosphere; negative values indicate uptake.

To convert to common alternative units:
- **umol m-2 hr-1**: multiply by 3.6
- **mg CH4 m-2 hr-1**: multiply by 3.6 x 16.04 / 1000 = 0.0577

## Flux Model Selection

Fluxes were calculated using the [goFlux](https://github.com/Qepanna/goFlux) R package, which fits both linear (LM) and Hutchinson-Mosier (HM) nonlinear models and selects the best fit based on AICc and diagnostic criteria. The `*_model` column indicates which model was selected; `*_quality` provides the quality flag.
