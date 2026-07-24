# Revision — Stocktake (state of all analyses)

> ## 🚨 OPEN CRITICAL ISSUE — ddPCR ~10× underestimate
> All ddPCR gene abundances (mcrA/pmoA/mmoX) may be **~10× underreported**: the
> pipeline (`02_harmonize_all_data.R:317`, `Conc × 75 / mass`) omits the 10×
> sample→reaction dilution (2.5 µL in 25 µL, per methods-paper SI Table 1).
> **Jon to ask Wyatt Arnold to confirm.** If confirmed, correct ALL ddPCR copies/g ×10.
> Uniform → ratios & relative results UNCHANGED. 16S (qPCR) NOT affected. LoD unaffected.
> Detail + plain walkthrough: `text/ddpcr_controls_LoD.md`.


NPH-MS-2026-56441 resubmission. All scripts in `revision/analysis/` reproduce
clean from scratch (12/12 pass). Inputs are the project's Zenodo/local data
(compiled + raw, git-ignored per repo convention). Nothing here edits existing
pipeline code or the manuscript.

**Git status:** the entire `revision/` folder is currently UNTRACKED (not
committed). Files are saved on disk. Committing is Jon's call — see end.

---

## A. Analyses DONE (script → key result → outputs)

| # | Analysis | Key result | Script | Outputs |
|---|---|---|---|---|
| 1 | **RMA/SMA robustness** of species ratio–flux | reproduces R²=0.513 (r=0.717); SMA slope 0.073 vs OLS 0.052; jackknife stable, sig 9/10 folds | `R2_species_aggregation_RMA.R` | `rma_*` |
| 2 | **Aggregate-x / keep-y** (Mark) | raw fails (skew); on pseudo-log, aggregate-x slope 0.076≈SMA, p=0.028; individual marginal (p=0.08) | `R3_aggregate_x_keep_y.R` | `aggregate_x_*` |
| 3 | **Per-regression figures** | 6-panel per predictor (raw/pseudo-log × indiv/agg-x/species) | `R_regression_figures.R` | `regfig_*.png` |
| 4 | **Variance partition** reconcile | manuscript 82.9% = OLS Method-2 (reproduced); mixed-model ICC 7.3% raw / 15.6% arcsinh | `R_variance_partition_review.R` | `variance_partition_review.txt` |
| 5 | **Whole-tree isotope sources** | δ¹³CH₄ −63‰, δ¹³CO₂ −20‰; Miller-Tans source CH₄ −65‰/CO₂ −35‰ (R²0.93) vs Keeling (R²0.13) | `R_isotopes_wholetree_sources.R` | `iso_wholetree_sources.*` |
| 6 | **Isotope provenance + HW/SW** | whole-tree=manuscript basis; tissue set=61 paired trees; NO sig HW→SW enrichment (31/61, p=0.58) | `R_isotopes_provenance_hwsw.R`, `R_isotopes_hwsw_explore.R` | `iso_hwsw_*`, `iso_provenance_summary.txt` |
| 7 | **Height-slope vs soil moisture** (R3) | steeper decline in wetter soil: per-tree r=−0.33 p=0.0002; per-species r=−0.65 p=0.03; birch wettest+steepest | `R_height_slope_vs_moisture.R` | `height_slope_*` |
| 8 | **Tree distribution + niches** (R3 L286) | DBH×landscape table; yellow birch wettest (VWC 42%) | `R_tree_distribution.R` | `tree_*`, `forest_inventory_*` |
| 9 | **Copies/g (real mass) + dry/wet** | wood mass=dry(freeze-dried), soil=fresh; heartwood:soil 34× (mixed) / ~21–25× (consistent dry); soil dry est via GWC | `R_copies_per_gram.R` | `copies_per_g_*` |
| 10 | **Fig 7b copies/g fix** | black-oak mcrA in copies/g dry (21/21 matched); profile shape preserved | `R_fig7b_copies_per_gram.R` | `fig7b_copies_per_gram.*` |
| — | Shared species prep helper | reproduces script-02 aggregation | `_prep_species_data.R` | (sourced) |
| — | ⚠️ SUPERSEDED (spurious tissue labels) | whole-tree bulk ok; tissue split invalid | `R_isotopes_co2.R` | `iso_co2_*` |

## B. Confirmed facts (from data / Black Oak xlsx)

- H = heartwood, S = sapwood (Jon confirmed).
- `sample_mass_mg` = mass weighed into extraction tube (xlsx: "TARGET WOOD 100mg / SOIL 250mg").
- **Wood freeze-dried before weighing → wood copies/g is DRY; soil fresh → WET.** (= R2's point.)
- No DNA extraction-efficiency correction anywhere → copies/g are conservative LOWER BOUNDS.
- Isotope runs: original script loaded 3 of 4 Picarro runs (omitted `205833`, ~all tissue). Whole-tree unaffected.
- Black-oak abundances (~10²) are NORMAL — manuscript heartwood median is 10²·⁷.
- Black Oak Tree Project Data.xlsx has the felled-tree destructive data: height cores (fresh+dry mass), bark, foliage/roots ("Misc Components - DNA"), stem fluxes, soil moisture/GWC → resolves R3's felled-tree/"AI hallucination" methods gap.

## C. Text drafts ready

`reviews/reviews_verbatim.md`, `response_to_referees_skeleton.md`,
`text/aggregation_framing.md`, `text/aggregate_x_finding.md`,
`text/felled_tree_methods.md`, `text/microbiology_corrections.md`,
`text/isotope_co2_finding.md` (with correction note).

## D. Decisions Jon has made

- Drop the H/S tissue isotope dataset; focus on whole-tree.
- Report copies/g on FRESH basis for all (common unit); also dry for wood, estimated dry for soil (organic/mineral) as SI note.
- Keep stand-level upscaling as a bounding exercise; don't lean on preprints.
- Title: "Tree stem methane exchange in upland forests reflects a balance of microbial production and consumption" (± species/functional-group language).

## E0. MANUSCRIPT RESTRUCTURING (text-assembly phase, in progress)

- **Decisions (Jon):** (1) NEW synthesis main fig = convergent hydrogenotrophy evidence (key taxonomy assoc + key pathway/function + isotopes), full heatmaps/isotopes → SI; absorbs S6-taxonomy + F6-highlights + §7 isotopes. (2) F9 upscaling → reworked: keep soil map + tree map, species rates + seasonal → SI, ADD component budget; reframe as bounding exercise for OUR stand (clarify not soften). (3) Methods S1–S3 → main. (4) Plant traits → SI (done). (5) Fixed panels (F2a/F4/F7b/F7c) swap in place.
- **Component budget figure — DRAFTED** `R_fig_component_budget.R` → `fig_budget_maps.png` (soil map + tree map + budget) & `fig_budget_panel.png`. Budget = measured extent (1.25 mg, 0.14%) + extrapolated extent (to 114 mg, ~13% full-WAI SCENARIO, unconstrained) as distinct components w/ zoom inset; flux SE propagated (small); surface-area/constant-emission flagged unconstrained. Numbers MATCH manuscript (L229/269: 0.14%→~13%). Maps reproduce existing F9 (can restyle). Seasonal bars → SI.
- **⚠️ CODE BUG found (not manuscript):** `manuscript_statistics.R:1514` computes WAI extrapolation as per-GROUND tree flux (1.25) × WAI 3.07 = 3.85 mg / 0.42% — a UNIT CONFLATION (Phi_tree is already per-ground, `02_rf_models.R:919`). Correct full-WAI = per-stem (1.25/0.0347) × WAI = ~114 mg / 13%. The manuscript TEXT is correct (114/13%); only the stale stat is wrong. Flag for code cleanup; not manuscript-affecting.
- **Budget figure — DESIGN LOCKED (v3).** color=process (blue sink/red source, matches Figs 6-8); (a)soil map +(b)tree map aligned same extent w/ bottom legends, tree color QUANTILE-spread to show skewed range; (c) 4-season bars (soil blue/tree red, own scales); (d) net waterfall (0→soil→tree offset→net −790), no zoom. Units: maps nmol m⁻²s⁻¹, budget mg m⁻²yr⁻¹. Honesty caption (net sink, RF R²=0.15/0.28, offset 0.14→13% scenario unconstrained). Micro-polish left: caption CH₄ subscript, soil-raster interp artifacts, opt. scale bar.
- **STILL TO BUILD:** the synthesis (hydrogenotrophy) main figure — convergent taxonomy + pathway/function + isotopes.

## E. New calculations / outputs STILL NEEDED

1. **pmoA and mmoX reported separately** in main text (R2 #1c) — small re-tab from ddPCR. [not built]
2. **Known/Putative methanotroph & methanogen taxa table** (R2 #2) — from taxonomy_key_16S + methanotroph_definitions. [not built]
3. **Fig 2A axis fix, Fig 4 legends, Fig 7c flux unit** (R2 #6) — DONE. Fig 2A: shared signed pseudo-log axis + defined zero-flux dashed line (`R_fig2a_consistent_axis.R`). Fig 4: legends already present after figure restructure (phylo-distance colorbar + Sample Type) — caption clarification only. Fig 7c: unit added, `CH4 flux (nmol m-2 s-1)` (`R_fig7c_flux_unit.R`; one-line patch at 09_felled_oak_profiles.R:221). Note: `text/figure_panel_fixes.md`.
4. **Plant-trait reframing** (editor + R3) — **UNBLOCKED & FIRST-PASS DONE.** TRY/wood-trait data were NOT missing: the sibling repo `../tree-gas-traits` already assembled per-species traits (TRY+GWDD+BAAD+wood pH+USDA porosity) and a full trait→internal-CH4 atlas (porosity R²=0.41; phylo λ=0.97). NEW bridge `R_traits_flux_microbial_bridge.R` joins those traits (by species code) to THIS paper's canonical species vars (via `_prep_species_data.R`; validity check reproduces ratio→flux R²=0.47≈0.513). Result: **wood & bark density suppress methanotrophs (ρ=−0.65/−0.73, p≈0.05/0.03) → tip methanogen:methanotroph balance (ρ=+0.73, p=0.03) → raise flux (ρ=+0.63, p=0.08)** — a plant-structural (O₂-barrier) mechanism = the editor's ask. Porosity predicts internal *concentration* but NOT flux (nice nuance). EXPANDED (`R_traits_multitrait.R`): 20 traits incl. ROOTS + trait & OUR-realized moisture (VWC niche) + wood chemistry; LOO-CV + FDR + collinearity. Validity check now reproduces ratio→flux EXACTLY (r=0.717, R²=0.513; the "0.47" was a pseudo-log-transform artifact in the check, not a data issue). **Honest verdict: the 80-test atlas survives NO FDR correction (n=10 underpowered) — so we lead with (a) out-of-sample LOO-CV and (b) pre-specified mechanism, NOT the scan.** What holds under LOO-CV: realized soil moisture→methanogen pool (R²=0.23, independent of density); wood/bark density→balance via methanotroph suppression (R²=0.20–0.29); **net flux itself only weakly trait-predictable (≤0.13)** — components predictable, sum noisy (reframes R3's low flux-R²). Robust NULLS answer R3: rooting depth→flux ρ=−0.19 p=0.61; realized moisture→flux ≈0 p=0.97; porosity→flux ≈0. Powered anchor = sibling n=29 concentration (λ=0.97, porosity R²=0.41). DESCRIPTIVE HEATMAP (`R_traits_heatmap.R` → `traits_heatmap.png`): traits × {mcrA, methanotroph, balance, flux}, grouped by category. Building it CAUGHT an overclaim — the O₂-diffusion-barrier mechanism I'd asserted is NOT supported (wood density→O₂-deficit ρ=−0.43 wrong sign; O₂-deficit→methanotroph n.s.), so density→methanotroph is a robust ASSOCIATION with OPEN mechanism (withdrawn the O₂ story). Internal-gas cols dropped (61% near-atmospheric). DIMENSION REDUCTION (`R_traits_dimreduction.R`): PCA PC1=35% wood-economics axis; response arrows point mechanistically-sensibly BUT no PC sig. correlates w/ any response; **supervised PLS LOO-CV: NOTHING predicts out-of-sample (flux −0.29, balance −0.05, mcrA −0.02, meth +0.07).** GYMNOSPERM CHECK (`R_gymnosperm_sensitivity.R`) — IMPORTANT: the density→methanotroph/balance pattern is driven by the 2 conifers (PIST, TSCA: low-density + high-methanotroph + consumption-tilted); it VANISHES angiosperm-only (density→balance +0.26→−0.05). So it's a tentative conifer/hardwood CATEGORICAL contrast (n=2, hypothesis-only), NOT a continuous density effect — corrects the earlier "robust mixed-model" claim (those p's were conifer-driven too). The ONE continuous effect surviving angiosperm-only = moisture→methanogen (mcrA~VWC 0.42→0.37). Don't code wood-type as ordinal porosity. Verdict: plant-trait signal is DESCRIPTIVE/DIRECTIONAL + a-priori single-trait hypotheses, NOT a predictive multivariable model. **SCOPE (Jon): this is an SI-figure + discussion add to appease the editor's plant-framing note ONLY — not central, not a repositioning. Do NOT bring in the separate n=29 cross-site concentration paper (Jon's separate work).** **CLEAN DELIVERABLE (agreed framing): `text/plant_traits_SI_and_discussion.md`** — a walk-through: (1) heatmap = apparent structure → (2) ordination (PCA biplot + Procrustes) = moderate/non-sig co-structure, conifers as misfits → (3) breakdown dumbbell (`traits_breakdown.png`) = density effects collapse without the 2 conifers → CONCLUSION: measured traits explain only a coarse hydrological + conifer/hardwood signal; **most among-species variation (esp. net flux) is not yet explainable.** Includes SI-methods + discussion paragraph + editor/R3 response, all drafted. Working log with all numbers/corrections: `text/plant_traits_bridge_finding.md`. Figures: traits_heatmap, traits_pca_biplot, traits_procrustes, traits_breakdown. **Remaining: at text-assembly, pick composite-vs-separate SI figs + fold paragraph in.**
5. **ddPCR negative controls / LOD** (R2 #1) — DONE. Via Arnold et al. 2024 (MEE methods paper) + Dryad raw data: probe-mcrA NTC = true zero (6/6), mmoX/pmoA NTC~0, LoD ~500 cells/100mg dry (empirical: 6/6 detect at 3 copies/rxn, 2/6 at 0.33), recovery ~20% (uncorrected → conservative). Final R2 response drafted (LoD + NTC + one log-scale sentence; **LoQ deliberately omitted**; values-below-LoD explained as probabilistic detection, not error). `text/ddpcr_controls_LoD.md`. **Remaining for Jon: (i) confirm the ×10 conversion with Wyatt [see 🚨 banner]; (ii) pull pmoA/mmoX primer sequences from Luton 2002 / Fuse 1998 / Bourne 2001 for the Methods table.**
6. **Isotopes (R2 δ¹³C)** — DONE. `ISOTOPES_final.R`; report δ¹³CH₄ −63.7 + δ¹³CO₂ −20 + ε_C ~54-64 (corroborating). Tissue/incubation sets confirmed to belong to Arnold et al. 2025.
7. **Felled-tree Methods subsection** — DONE (`text/felled_tree_methods.md`). Same tree as Arnold & Gewirtzman 2025 (Nature); climbed standing for flux/gas, then felled for destructive sampling. Confirm flux/gas were pre-felling.
8. Manuscript-finalization sweeps: gene-name italics, primer sequences, unit clause, define ddPCR, move Fig S6 to main, Ranniku 2023 cite, Gauci discussion, mechanistic sentence fixes (drafts in `text/microbiology_corrections.md`). [text-only]

- **R3 L342 (defend exponential expression) — DONE.** Finding: NO exponential model is fit anywhere (height=linear mixed-effects, flux=goFlux, upscaling=RF). The "exponential" is only the CLAIM that soil-transported CH₄ declines exponentially with height in wetlands (contrast to our upland uniform-profile hypothesis; cited Pangala 2013/Barba 2019/Jeffrey/Anttila). Defense: it's empirically established + mechanistically expected (first-order radial diffusive loss → C(h)=C₀e^(−kh)); and we DON'T fit it — clarify we only use it as the predicted soil-transport signature, testing our data with assumption-free LME. Drop-in clarification + R3 response: `text/exponential_height_defense.md`. (PDF line #s don't map exactly; confirm at finalization — but it's the only 'exponential' in the ms.)

- **R2 #1e primer sequences — RESOLVED (Jon provided SI Table 2), pending Wyatt confirm.** Sequences now in hand: pmoA 189f `GGNGACTGGGACTTCTGG` / mb661r `CCGGMGCAACGTCYTTACC` (**Bourne 2001**); mmoX 536f `CGCTGTGGAAGGGCATGAAGCG` / 898r `GCTCGACCTTGAACTTGGAGCC` (**Fuse 1998**); mcrA = probe assay (L117); ML mcrA (Luton 2002) was superseded. **CITATION FIX: ms L117 currently cites Luesken 2011/McDonald&Murrell 1997/McDonald 1995 for pmoA/mmoX — WRONG, replace with Bourne 2001 + Fuse 1998.** Full drop-in Methods text + SI table + R2 response: `text/primer_sequences.md`. TO CONFIRM w/ Wyatt: exact primers + whether thermocycling follows Arnold 2024 (MEE) SI Table 2 or Arnold 2023 (Microbiol Spectr).

- **R2 #3 (16S functional-prediction caveats) — DONE (measured stance, per Jon).** `R_faprotax_caveat_metrics.R` quantifies limits: FAPROTAX annotates only **46.9%** of HW abundance; **32% of ASVs** genus-resolved; categories non-additive (dark-H₂ 3.6% + hydrogenotrophic methanogenesis 3.2% → total H₂ ≥~7%). STANCE: functional inference is standard, NOT retracting; stop over-interpreting ABSOLUTE %s; foreground the robust **between-tissue** differential (HW enriched in methanogenesis 32×/H₂ 3.1×/fermentation 2.5×; SW enriched in aerobic 1.5×/methanotrophy 2.3×) = expected core/outer gradient matching ddPCR + isotopes. CORRECTED my earlier error: 'aerobic>fermentation in HW is backwards' is WRONG — heartwood O₂ is not systematically low but VARIABLE (Arnold Nat. 2025), and aerobic is relatively LOWER in HW (0.7×). H₂ non-additivity CORROBORATES hydrogenotrophy (δ¹³CH₄ −63.7 + Methanobacteriaceae). Frame as corroborated-where-independent-evidence-agrees, hypothesis-generating otherwise. **PROVENANCE (Jon's catch): paired depth-resolved HW/SW O₂ profiles = companion Nature paper (Arnold & Gewirtzman 2025), NOT measured here — cite as prior work. THIS study has only single internal (bark-to-pith) gas per tree (n=157, ~1–22% O₂, variable) + internal CH₄–O₂ neg. corr. (Fig S8). In-study redox corroboration = OUR ddPCR compartment split (methanogens Inner / methanotrophs Outer, Figs 4,5), NOT our own O₂ profiles.** **DISCLOSURE (Jon, keep minimal): TISSUE-level FAPROTAX was in the Nature paper; Fig S3 here adds the SPECIES-level view; PICRUSt2 (Figs 6,S4,S5) is new. Minor point / one SI fig → ONE parenthetical at L185 or Fig S3 caption ("tissue-level predictions also reported in Arnold & Gewirtzman 2025, extended here to species level"); NO methods paragraph. 16S-data overlap already at L115.** Caveat + R2 response + provenance: `text/faprotax_picrust_caveats.md`.
9. **Upscaling → bounding-exercise demotion** (R3); **wood-vs-soil pool scenario** available (`R_scaling_wood_soil_pool.R`). [text + decided]
10. **Full response-to-referees assembly** + **title** finalization. [text]

## F. Save / commit

Recommend: create a branch and commit `revision/` so nothing is lost. Not done
yet (awaiting Jon's go-ahead). Input data stays Zenodo-only per repo convention;
the black-oak extracts in `revision/data/black_oak/` WILL commit (small, derived).
