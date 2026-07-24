# Panel (b) — FAPROTAX function selection

## Provenance of the values

All HW/SW values are reproduced from the FAPROTAX pipeline (microeco
`cal_spe_func_perc(abundance_weighted = TRUE)` on the rarefied wood 16S table, merged by
core type Inner=heartwood / Outer=sapwood). The full 49-function output is saved to
`revision/outputs/FAPROTAX_all_functions_HW_SW.csv`; the figure reads the shown subset from it.
The six functions in the earlier draft were value-correct but a **non-systematic hand-pick**;
this version applies an explicit inclusion rule.

## Inclusion rule (applied to the full FAPROTAX output)

Show a FAPROTAX function iff:
1. it is an **energy / redox metabolism** function (not a host-association, pathogen, gut, or
   phototrophy annotation — those describe lifestyle/habitat, not wood biogeochemistry); AND
2. mean relative abundance **≥ 0.1%** in at least one compartment (drops near-zero noise); AND
3. it is not a **redundant parent / duplicate** annotation. Excluded on this ground:
   - `chemoheterotrophy` (parent of aerobic/anaerobic_chemoheterotrophy)
   - `methanogenesis`, `methanogenesis_by_CO2_reduction_with_H2`,
     `methanogenesis_by_reduction_of_methyl_compounds_with_H2` — all similarly HW-enriched;
     `hydrogenotrophic_methanogenesis` shown as the single representative methanogenesis route.

## Shown functions (10), with verified values

Heartwood-enriched (warm):
| function | HW % | SW % | log2(HW/SW) |
|---|---|---|---|
| hydrogenotrophic methanogenesis | 3.21 | 0.10 | +5.00 |
| dark hydrogen oxidation | 3.61 | 1.16 | +1.64 |
| fermentation | 16.77 | 6.65 | +1.33 |
| anaerobic chemoheterotrophy | 14.65 | 7.05 | +1.06 |
| methylotrophy | 0.88 | 0.53 | +0.73 |

Sapwood-enriched (cool):
| function | HW % | SW % | log2(HW/SW) |
|---|---|---|---|
| nitrate reduction | 0.84 | 1.06 | −0.34 |
| aerobic chemoheterotrophy | 28.90 | 42.53 | −0.56 |
| methanol oxidation | 0.22 | 0.49 | −1.16 |
| methanotrophy | 0.23 | 0.52 | −1.18 |
| cellulolysis | 0.84 | 2.47 | −1.56 |

## Consistency with panel (c) and the redox story

- **Nitrate reduction is the only substantial redox function** (≥0.1%) and it is *depleted* in
  heartwood — independently echoing panel (c)'s finding that assimilatory nitrate/sulfate reduction
  anticorrelate with *mcrA*.
- **Dissimilatory redox respiration is essentially absent:** denitrification (all variants),
  sulfate respiration, sulfur respiration, and manganese oxidation are all <0.03% in both
  compartments. Report this in text so the "no alternative-TEA respiration" claim is explicit and
  not an omission.

## Changes vs the earlier 6-function draft
- **Added** `anaerobic_chemoheterotrophy` (14.65% HW, the 2nd most abundant anaerobic function — a
  clear omission before), `methanol_oxidation`, `cellulolysis` (aerobic-C contrasts), and
  `nitrate_reduction` (redox, for panel-c consistency).
- All six original values were correct and are retained.
