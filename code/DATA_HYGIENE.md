# Data hygiene ledger

Every rule that removes, flags, corrects or defines data between the raw files and the
results, in pipeline order. One row per rule: what it does, why, where it is applied,
how many rows it touches, and what kind of action it is.

**Principle.** Measurements are flagged, not deleted. A measurement is removed only for
procedural failure (the instrument or chamber was not measuring what it was meant to).
Statistical filtering of plausible values is not used anywhere: it biases means
(an r² > 0.7 screen would raise the mean stem flux by 144%; Methods S2).

Action types: **removed** (row dropped from analysis, still in the archive with a flag
where possible) · **flagged** (row kept, column says why) · **corrected** (value changed,
original recoverable from raw) · **defined** (a convention every downstream step follows).

Counts are as of 2026-10-01. *Status* marks rules whose effect is not yet logged by the
pipeline itself; those are being wired so the counts are produced, not typed.

## 1. Stem and soil flux

| # | rule | action | applied in | rows affected | why |
|---|---|---|---|---|---|
| F1 | July 2021 campaign dates: the import hard-coded one base date onto a 19-day campaign | corrected | 05_model/01_load_and_prep_data.R (per-deployment dates from the field sheets; base date is a logged fallback only) | 46% of July 2021 rows had August dates labelled July | wrong dates mis-assign temperature and moisture |
| F2 | December 2020 monthly-survey files handled separately (date format) | corrected | 02_flux/semirigid/05_fix_december.R | December 2020 campaign | otherwise not matched to its aux file |
| F3 | Failed goFlux fits reprocessed with manual window corrections | corrected | 02_flux/static/03_fix_failed.R | *status: count to be logged* | fits that failed on default windows |
| F4 | Untagged monthly stems: closure windows chosen by hand, stems given species | corrected | interactive rescue; products in data/processed/flux/untagged_rescue/ | 49 fluxes on 7 stems | otherwise silently dropped; cannot be regenerated without repeating the clicking |
| F5 | Analyzer cell volume 28 cm³ (was goFlux's 70 cm³ example value) | corrected | lib/chamber_constants.R → 02_flux/apply_auxfile_vtot.R | every flux; −0.5 to −6.8% | wrong system volume scales every flux |
| F6 | Closure duration in seconds for the minimum detectable flux (2021 logged every 5 s) | corrected | 02_flux/mdf_FINAL_precision_and_detection.R | MDF of every 2021 measurement (was 5× too high) | MDF used sample count in place of seconds |
| F7 | Chamber not at ambient at closure: starting CH₄ > 1.5× the campaign median | removed | 02_flux/qc_c0_screen.R → outputs/data/qc_excluded_measurements.csv | 2 soil deployments (5.7× and 24× median); 0 stems | procedural failure |
| F8 | Detection status from the MDF (90%) and slope p-value | flagged | 02_flux/mdf_FINAL_precision_and_detection.R | 772 of 1,146 stem and 260 of 286 soil measurements detected | reported, never used to exclude |
| F9 | Chamber-type offset: semirigid (monthly) relative to rigid, estimated on the asinh scale with environment and species | corrected | 05_model/02_rf_models.R | all semirigid stem rows entering the model | systematic chamber difference; chamber type has no counterpart on inventory stems. *Status: fitted β is printed but not saved — to be written to outputs* |
| F10 | Per-campaign flux filenames (three campaigns had overwritten one shared name) | defined | 02_flux/assemble_campaign_flux.R, migrate_campaign_filenames.R | all campaigns | a filename collision silently dropped campaigns |

## 2. Training population for the stem model

| # | rule | action | applied in | rows affected | why |
|---|---|---|---|---|---|
| T1 | Dead stems | flagged; excluded from model training | lib/dead_stems.R | 43 deployments on 10 stems | the inventory contains no dead stems; kept in all descriptive analyses |
| T2 | Root-crown measurement | excluded from model training | 05_model/03_export_canonical_tables.R | 1 deployment | not a stem surface |
| T3 | Every deployment carried in one table with `in_rf_training` and `exclusion_reason` | defined | 05_model/03_export_canonical_tables.R | 1,191 deployments; 1,147 train | row accounting is checked by check_consistency.R |
| T4 | Species calibration needs ≥ 5 training measurements | defined | lib/species_calibration.R | levels with < 5 get ratio 1 | one measurement cannot calibrate a species |

## 3. Internal gas and isotopes

| # | rule | action | applied in | rows affected | why |
|---|---|---|---|---|---|
| G1 | GC calibration: fits weighted to minimise relative error; linear CH₄ to 5,029 ppm with interpolation above; quadratic CO₂, N₂O, O₂; no truncation at zero | corrected | 01_import/03_process_internal_gas.R, 03b_process_internal_gas_2024.R | all internal-gas samples; 1 sample above the top standard reported at it | the earlier calibration clamped low CH₄ to 0 |
| I1 | Not a whole-tree sample: atmosphere and room air, incubations, H/S tissue pairs, standards | removed (from the whole-tree set) | lib/isotope_samples.R | 18 + 93 + 122 + 12 of 408 Picarro rows | different sample types |
| I2 | A re-run replaces its original | corrected | lib/isotope_samples.R | 2 | the redo is the valid measurement |
| I3 | CH₄ below 1.5 ppm | removed | lib/isotope_samples.R | 33 of 161 candidates | δ¹³C unreliable near ambient |
| I4 | δ¹³CH₄ outside −115 to +25‰ | removed | lib/isotope_samples.R | 13 | not CH₄ from any source |

## 4. Inventory

| # | rule | action | applied in | rows affected | why |
|---|---|---|---|---|---|
| V1 | 2018 census positions recomputed from quadrat and local coordinates | corrected | 01_import/inventory_build.R | 29 stems in quadrat 907 had been placed 20 m south | geometry |
| V2 | Subquadrat-only positions placed at the subquadrat centre (± 2.5 m) | corrected | inventory_build.R | 38 stems | otherwise unlocated |
| V3 | Stems far from their own quadrat label | flagged unlocated, kept | inventory_build.R | 14 | position contradicts the census |
| V4 | Repeated tags kept (tags are reused; 70 carry different species) | corrected | inventory_build.R | 169 stems restored | had been de-duplicated away |
| V5 | 2018 multi-stem trees from the census "multiple stems" sheet | corrected | inventory_build.R, data/raw/inventory/fg_2018_multiple_stems.csv | 40 stems on 17 trees | diameters live only on that sheet |
| V6 | Diameter decimal-shift repair | corrected | inventory_build.R | 8 of 8,208 stems | data-entry errors |
| V7 | Stems in the uncensused notch | removed | inventory_build.R | 3 | outside the censused stand |
| V8 | Stems without coordinates kept in the stand total | flagged | inventory_build.R | 34 | they are in the stand; only their position is unknown |
| V9 | No dead stems: neither census records a dead status | defined | — | — | why the model is trained on live stems only |

## 5. Molecular data

| # | rule | action | applied in | rows affected | why |
|---|---|---|---|---|---|
| M1 | 16S: chloroplast and mitochondrial reads removed, rarefied to 3,500 | removed | 16S processing (lib/build_phyloseq.R) | 35 wood samples fall below 3,500 after plastid removal | host DNA; dropout disclosed in Methods S3 |
| M2 | ddPCR copies g⁻¹ = copies µL⁻¹ × 75 µL elution ÷ mass; dry mass for wood, fresh mass for soil | defined | 03_merge/02_harmonize_all_data.R | all ddPCR values | *Status: possible missing ×10 template dilution, pending confirmation* |
| M3 | ddPCR "loose" vs "strict" positive calls; analyses use loose | defined | harmonize_all_data.R | all ddPCR values | stated convention |
| M4 | Methanotroph classification by methane-monooxygenase capacity (Known / Putative / listed, not counted) | defined | lib/load_methanotroph_definitions.R, methanotroph_definitions.csv | all 16S methanotroph ASVs | growth phenotype is not the criterion |
| M5 | 16S and ITS copy numbers are the sequencing facility's qPCR values | defined | sample metadata | — | not ddPCR; the ddPCR conversion does not apply |

## 6. Environment

| # | rule | action | applied in | rows affected | why |
|---|---|---|---|---|---|
| E1 | Soil environment archive is the union of the per-date sheets and the 2020 summary | corrected | 03_merge/compile_soil_env.R | 291 records, 22 dates | earlier builds read one source and lost dates |
| E2 | Stream points in the December moisture survey set to 100% VWC | defined | 04_drivers/moisture_surface.R | stream boundary points | an assumption, not a measurement; anchors the wet end |
