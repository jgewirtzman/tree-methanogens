# Response to Referees — skeleton (NPH-MS-2026-56441)

Point-by-point, one entry per reviewer comment. Status tags:
**[DONE]** drafted in revision/ · **[EASY]** quick fix at finalization ·
**[NEEDS JON]** decision/field-fact required · **[NEW SCRIPT]** re-analysis needed.

> Note (NP process): resubmission = new submission, likely new referees, usually
> one shot. Every response must let the manuscript stand on its own. Keep tone
> gracious even where reviewers missed content already present (data
> availability). Editor + R3 agree on the load-bearing issue: **plant framing.**

---

## Editor (Martin De Kauwe) — plant framing

> "greater attention… to the framing of the manuscript and its linkage to plant
> sciences… how do differences among plant species drive variation in stem
> methane flux or associated microbial activity?"

**Response plan [NEEDS JON + NEW SCRIPT]:**
- Reframe intro/objectives around plant traits (angiosperm vs gymnosperm; ring-
  vs diffuse-porous; wood chemistry, moisture, heartwood formation, rooting
  depth) as controls on flux and on the microbial balance.
- Fold in TRY trait data (from the Covey internal-gas manuscript — **must be
  brought into this repo**) as an SI table + brief discussion.
- Reorganize the species results around traits; connect the "confounds" of the
  aggregation analysis (anatomy/chemistry) explicitly as the *plant-side
  mechanism* that sets up the microbial balance — this answers De Kauwe and R3's
  circularity worry at once.
- New title (see below) puts scope/mechanism forward; consider adding species/
  functional-group language per Mark.

---

## Referee 1 (near-accept)

| # | Comment | Status | Plan |
|---|---|---|---|
| 1.1 | Discuss contrast with other stem-flux results (Gauci et al. 2024 Nature) | [EASY] | Extend existing L93-96 / Discussion: net small source vs Gauci net-uptake; reframe as net-vs-gross, not contradiction. Gauci already cited. |
| 1.2 | Cryptogamic/*Melaleuca* microbiome work (Jeffrey) | [EASY] | Add 1–2 sentences; Jeffrey already cited elsewhere. |
| 1.3 | Higher birch flux consistent with birches on organic/swampy soils; cite Ranniku et al. 2023 | [EASY] | Add Ranniku et al. 2023 (doi:10.1088/2515-7620/acd7c7) — the one genuinely new citation — to the *B. alleghaniensis* discussion. |

---

## Referee 2 (rigorous microbiologist)

| # | Comment | Status | Plan (detail in revision/text/microbiology_corrections.md) |
|---|---|---|---|
| 2.1a | Low-biomass qPCR: confidence of detection; **negative controls reported explicitly** | [NEEDS JON] | Our claims rest on **ddPCR**, which had negative controls reading **0 copies**; report them + LOD/LOQ explicitly, and point to the methods paper for how thresholds are set (novel low-false-positive mcrA probe). This is the single highest-stakes add. Confirm the NTC/LOD numbers to state. |
| 2.1b | Units dry vs wet weight | [EASY] | Standardize + state basis. |
| 2.1c | Report pmoA and mmoX **separately** | [NEW SCRIPT?] | Report separately in Section 4 (already separated in SI). Possibly `revision/analysis/pmoa_mmox_separate.R`. |
| 2.1d | Italicize gene names | [EASY] | Global sweep. |
| 2.1e | L248-9: report pmoA/mmoX primer sequences in Methods | [EASY] | Add primer table to Methods/SI. |
| 2.2 | 16S methanotroph classification (Beijerinckiaceae, Methylacidiphilaceae); provide **supp. table of known/putative** taxa | [NEW SCRIPT] | Build Table S(new) of taxa → Known/Putative methanotroph & methanogen. Data exist (methanotroph_definitions.csv; taxonomy_key_16S.csv). `revision/analysis/putative_taxa_table.R`. |
| 2.2b | Methylacidiphilaceae not "verified" methanotrophs (mesophilic lack genes) | [DONE-draft] | Soften L177; cite peatland + tree-stem genomes. |
| 2.3 | Caveats of 16S-based functional prediction (fermentation, H₂ oxidation poorly predicted; "only minor fraction fermenting in hypoxic habitat" implausible; 3.6% dark-H₂ vs hydrogenotrophic methanogens) | [NEEDS JON] | Add explicit caveat paragraph on PICRUSt/FAPROTAX limits; reconcile the 3.6% dark-H₂ figure with hydrogenotroph abundance (reframe as lower bound / inference, not measurement). |
| 2.4a | L412-3 *Methylocapsa* sMMO/multi-carbon overstated | [DONE-draft] | Corrected. |
| 2.4b | L623 acetyl-CoA & RuMP ≠ fermentation | [DONE-draft] | Corrected at L253; L187 already fine. |
| 2.4c | L695-700 sMMO low-affinity (>500 ppm) | [DONE-draft] | Added; reframed as support for net-source interpretation. |
| 2.5 | **Open data** — "none of the raw data were made available" | [EASY-but-graceful] | **Data ARE released**: full Data Availability Statement (SRA BioProject for 16S; Zenodo DOI for flux/ddPCR/gas/characterization; GitHub for code). Note neutrally that the statement likely didn't render in the review PDF; confirm links resolve. Don't be smug. |
| 2.6 | Fig 2A axis; Fig 4a/c legends; Fig 7b copies g⁻¹; Fig 7c unit | [NEW SCRIPT] | Figure-panel fixes; new scripts only. |
| 2.7a | **Title** too review-like | [DONE-draft] | New title (below). |
| 2.7b | L607-10 δ¹³C caution + source CO₂ not measured + cite | [DONE-draft] | Strengthen L251; add Whiticar 1999 / Conrad 2005. |

---

## Referee 3 (most consequential; aligned with editor)

| # | Comment | Status | Plan |
|---|---|---|---|
| 3.1 | **Plant-light / microbe-heavy**; "species" used 7× with no account of what about a species drives flux | [NEEDS JON + NEW SCRIPT] | Same as editor — plant-trait reframing + TRY data. Load-bearing. |
| 3.2 | Multi-year campaigns hard to track; **felled tree not in methods**; methods in SI; length/figures | [NEEDS JON] | Add felled-tree Methods subsection (revision/text/felled_tree_methods.md — needs field facts); add campaign-structure overview (L157); move core Methods (S1-S3) into main text (L263); add species×DBH×location table (L286). |
| 3.3 | Weak RF upscaling "unhelpful"; overstating global magnitude | [NEEDS JON] | Keep stand-level scaling but **demote to a clearly-labeled bounding scenario** within the budget discussion, not a standalone numbered Result; foreground OOB R² honestly; play R3 (don't scale) vs R1 (reconcile w/ Gauci) against each other to justify a hedged retention. Don't lean on companion preprints. |
| 3.4 | Height-decline ≠ soil source (logical gap); *B. alleghaniensis* has 4th-highest heartwood mcrA; are highest birch fluxes in wettest soil? (Fig 3 suggests not) | [NEEDS JON] | Argue, don't assert. Foreground the cross-species correlation (r = −0.76 between height-effect and soil methanogens) as the actual evidence; concede internal-production alternative for birch; address the "other wet-soil species don't show the trend" point directly. This is also the R1↔R3 birch tension — resolve without contradiction. |
| 3.5 | RF explains only 28% soil / 15% stem flux variability | [EASY] | Acknowledge; frame RF as descriptive/bounding not predictive; ties to 3.3. |

### R3 line comments

| Line | Comment | Status | Plan |
|---|---|---|---|
| L32 | define ddPCR | [EASY] | Define at first use. |
| L38 | which species? | [EASY] | Name species in abstract result. |
| L66 | claims need citations | [EASY] | Add cites. |
| L103-110 | which species? | [EASY] | Specify. |
| L112 | reword to "We investigated net methane stem flux and associated microbial community…" | [EASY] | Adopt near-verbatim. |
| L125 | species distribution across transect? | [EASY] | Add (ties to L286 table). |
| L157 | describe campaigns + goals first, then protocol | [NEEDS JON] | Restructure Methods opening; helps 3.2. |
| L171 | felled tree / "AI hallucination?" | [NEEDS JON] | Felled-tree Methods subsection. |
| L189 | "fits" for flux calc, "models" for stats | [EASY] | Terminology sweep. |
| L215 | which data co-located? | [EASY] | Clarify. |
| L261 | moisture of what? | [EASY] | Specify soil vs stem. |
| L263 | why methods in SI? | [NEEDS JON] | Move core methods to main text. |
| L275,L563 | "transitional wetland" = wetland or intermediate? | [EASY] | Define term consistently. |
| L286 | table of #trees/species + DBH by location | [NEW SCRIPT] | Build it (also serves plant framing). `revision/analysis/species_dbh_table.R` from forest_inventory.csv / tree_properties.csv. |
| L342 | defend the exponential expression in Methods | [NEEDS JON] | Add justification/citation for the exponential height model. |
| L362 para | move to discussion | [EASY] | Relocate. |
| Results §8 | move discussion content to discussion | [EASY] | Relocate. |
| L599 | move Fig S6 to main text | [EASY] | Promote (load-bearing for syntrophy argument). |
| L664 | precedent for species aggregation? | [DONE-draft] | revision/text/aggregation_framing.md. |
| L668-674 | "seems circular" | [DONE-draft] | Rewritten as positive scale-dependence argument + SMA robustness. |
| L741 | towers exist in upland forests, few CH₄ sensors | [EASY] | Reword. |

---

## Proposed title

**"Tree stem methane exchange in upland forests reflects a balance of microbial
production and consumption"**

Mark suggested optionally lengthening to work in species / functional groups for
the editor's plant fix, e.g.:
- "*Species differences in* tree stem methane exchange in upland forests reflect
  a balance of microbial production and consumption" — **[NEEDS JON/coauthors]**

---

## Cover letter (to editor, not referees) — bullets to make
- Thank editor; state the central change = plant-trait reframing (traits →
  microbial balance → flux), directly answering the scope concern.
- Note the new title, the demoted upscaling, the added felled-tree methods, and
  the explicit species-level trait framing.
- Neutrally note the data were already openly archived (in case it reassures).
