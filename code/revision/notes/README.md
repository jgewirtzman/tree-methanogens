# Revision workspace — NPH-MS-2026-56441

**➡️ See [STATUS.md](STATUS.md) for the full stocktake** (every analysis, result,
decision, and what's still needed). This README is the file index.

**Rule:** everything here is NEW. No existing pipeline script or manuscript file
is modified. Merge into the main pipeline / manuscript later only with Jon's
explicit approval.

Decision: Reject, resubmission encouraged (09-Jun-2026). Due 07-Dec-2027; target
fall 2026. All 3 referees positive. Load-bearing issue (editor + R3): plant framing.

## Contents

- `reviews/reviews_verbatim.md` — full verbatim editor + R1/R2/R3 comments (reference; do not edit)
- `response_to_referees_skeleton.md` — point-by-point plan, every comment, with status tags
- `analysis/_prep_species_data.R` — shared, side-effect-free data prep (reproduces script 02's aggregation; sourced by the analyses below)
- `analysis/R2_species_aggregation_RMA.R` — reproduces the species aggregation (ratio R²=0.513, r=0.717) and adds OLS-vs-SMA/RMA robustness + leave-one-species-out jackknife + variance partition
- `analysis/R3_aggregate_x_keep_y.R` — Mark's aggregate-x/keep-individual-y model vs individual vs species-medians, on raw AND pseudo-log flux; makes `aggregate_x_ratio_panels.png`
- `analysis/R_regression_figures.R` — per-regression 6-panel figures (raw/pseudo-log x individual/aggregate-x/species) for each predictor → `regfig_*.png`
- `analysis/R_variance_partition_review.R` — reconciles 82.9% (manuscript Method 2) vs 91% (species-only) + mixed-model ICC (raw/arcsinh) → `variance_partition_review.txt`
- `analysis/R_isotopes_co2.R` — d13C-CO2, CH4-vs-CO2 fractionation (alpha_C/eps_C), Keeling; validates d13CH4=-63.7 → `iso_co2_*`
- `analysis/R_tree_distribution.R` — R3 L286 DBH x landscape table + species moisture niches (yellow birch wettest) → `tree_*` + `tree_distribution_panels.png`
- `analysis/R_copies_per_gram.R` — real-mass copies/g (fixes 100mg proxy) + wood dry-weight factors → `copies_per_g_*`
- `analysis/R_isotopes_provenance_hwsw.R` — CORRECTED isotope provenance (whole-tree vs paired HW/SW vs incubation) + heartwood→sapwood enrichment test → `iso_provenance_summary.txt`, `iso_hwsw_*`. Supersedes the tissue split in R_isotopes_co2.R.
- `analysis/R_height_slope_vs_moisture.R` — R3's question: per-tree/species flux-height slope vs soil VWC + birch wettest-soil check → `height_slope_*`
- **`analysis/ISOTOPES_final.R` — CANONICAL settled isotope workflow.** Emits `ISOTOPES_methods.md`, `ISOTOPES_results.md`, `ISOTOPES_fig1_source.png`, `ISOTOPES_fig2_crossplot.png`, `ISOTOPES_sample_table.csv`. The `R_isotopes_*` scripts below are exploratory/superseded (kept for provenance).
- `analysis/R_isotopes_wholetree_sources.R` — whole-tree CH4/CO2 source terms via Keeling AND Miller-Tans (loads all 4 Picarro runs) → `iso_wholetree_sources.*`
- `analysis/R_fig7b_copies_per_gram.R` — Fig 7b fix: black-oak mcrA in copies/g dry using masses from Black Oak xlsx → `fig7b_copies_per_gram.*`
- `analysis/R_isotopes_by_species.R` — is the isotope source one pool or species-specific? Threshold-sensitivity shows species difference is low-CH4 noise → aggregate is appropriate → `iso_by_species.*`
- `analysis/R_isotopes_hwsw_explore.R` — all 61 HW/SW pairs, no conc filter (still no enrichment) → `iso_hwsw_explore.*`
- `analysis/R_scaling_wood_soil_pool.R` — standing methanogen gene pool wood vs soil per m² (wood ~20-80% of soil) → `scaling_wood_soil_pool.*`
- `data/black_oak/` — extracted from Black Oak Tree Project Data.xlsx (extraction masses, core fresh/dry masses, soil GWC)
- `outputs/` — all results (tables, figures, summaries)
- `text/aggregation_framing.md` — draft discussion/methods text answering R3 L664/L668-674 (the "circular"/aggregation asks)
- `text/aggregate_x_finding.md` — what the aggregate-x/keep-y test showed + decisions for Jon/Mark (raw vs pseudo-log; whether to retire the "exactly null individually" claim)
- `text/isotope_co2_finding.md` — CO2 isotope findings + response-to-R2 language
- `text/felled_tree_methods.md` — draft felled-tree Methods subsection (⚠️ needs field facts confirmed)
- `text/microbiology_corrections.md` — every R2 correction: current text → proposed fix

## Open items needing Jon (see response skeleton [NEEDS JON] tags)
1. **TRY trait data** — not in this repo; bring over from the Covey internal-gas manuscript for the plant-trait reframing.
2. **Felled-tree field facts** — same tree climbed then felled, or two trees? dates/method/tissue protocol (see felled_tree_methods.md CONFIRM list).
3. **ddPCR negative-control / LOD numbers** — to report explicitly for R2 #1.
4. Title final form (with/without species/functional-group language).
5. Which figure-panel fixes (Fig 2A, 4, 7b/c) and analysis re-tabs (pmoA/mmoX separate; putative-taxa table; species×DBH table) to build next — all as new scripts.

## Reproduce the analysis
```
Rscript revision/analysis/R2_species_aggregation_RMA.R
```
