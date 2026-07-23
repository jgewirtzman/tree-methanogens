# Isotopes — settled analysis & what to report (answers R2 L607-610)

Supersedes earlier drafts. Scripts (all exclude calibration standards; the
high-concentration standard S3c = −391‰ at 16,105 ppm otherwise hijacks
concentration-weighted fits):
- `R_isotopes_wholetree_sources.R` — source terms (Keeling / Miller-Tans)
- `R_isotopes_within_tree_fractionation.R` — per-tree apparent fractionation (primary for ε_C)
- `R_isotopes_by_species.R` — tests source homogeneity across species
- `R_isotopes_provenance_hwsw.R`, `R_isotopes_hwsw_explore.R` — tissue set (DROPPED)
- `R_isotopes_co2.R` — SUPERSEDED (invalid tissue labels); do not use

## Provenance (which data)
- **Whole-tree single-borehole samples = the manuscript basis** (~1 sample/tree,
  ~130–141 with reliable data). Median δ¹³CH₄ reproduces the manuscript −63.7‰.
- **Paired heartwood/sapwood tissue set (61 trees): DROPPED** — shows no
  significant HW→SW enrichment (31/61 pairs, p=0.58) even unfiltered.
- **Incubations: unrelated, disregarded.**

## The settled numbers (internally consistent — same samples per row)

| subset | δ¹³C-CH₄ | δ¹³C-CO₂ | within-tree ε_C |
|---|---|---|---|
| **bulk** (all, ≥1.5 ppm) | **−63.7‰** | −20‰ | **~44‰** |
| reliable (≥50 ppm CH₄) | −70‰ | −18‰ | ~55‰ |

- **Source estimates (whole-tree, standards excluded) — use Keeling / atmosphere
  correction, NOT Miller-Tans:**
  - **Miller-Tans is fragile here** (do not use): it swings −65 → −89‰ as the top
    1–5 highest-CH₄ points are dropped (33% of leverage in 3 points). Its R²=0.93
    is an artifact of concentration-weighting, not fit quality.
  - **Keeling is robust**: −79‰ [95% CI −92, −66], stable to dropping high- or
    low-CH₄ points (low R²=0.13 reflects low-CH₄ scatter, but the intercept is
    anchored).
  - **Atmosphere mass-balance (per-sample) ≈ −70‰**, insensitive to the assumed
    atmospheric endmember (−70 to −71 for δ_atm −47 to −35).
  - Interpretation: the *tree-produced* CH₄ source (~−70 to −79‰) is MORE depleted
    than the bulk measured −63.7‰ because the measured gas is diluted by
    atmospheric CH₄ (−47‰). This strengthens the hydrogenotrophy read.
- Species do **not** differ in source above the noise floor: apparent among-species
  δ¹³CH₄ difference (p=0.02 at ≥1.5 ppm) vanishes at ≥5–10 ppm — a low-CH₄ noise
  artifact. → the aggregate is appropriate.

## KEY conceptual points (hard-won)
1. **δ¹³C-CH₄ (a composition) and ε_C (a difference) are different quantities.**
   −63.7‰ CH₄ pairs with ε_C ≈ 44‰; ε_C ~55‰ goes with the more-depleted −70‰
   high-CH₄ CH₄. Never mix subsets.
2. **ε_C must be computed WITHIN each tree**, not from separately-pooled source
   terms. CH₄ and CO₂ concentrations barely co-vary (r=0.26), so a pooled CO₂
   source (−35‰, from high-CO₂ trees) is NOT the CO₂ that co-occurs with high CH₄
   (−18‰). Pooled ε_C (≈32‰) is wrong; within-tree ε_C (~55‰) is valid.
3. The pooled source term is a **between-tree, ecosystem-average** mixing line
   (1 sample/tree), not a within-tree source. Frame it as such if reported.

## WHAT TO REPORT
- **Primary:** δ¹³C-CH₄ = **−63.7‰** (bulk median) — depleted, hydrogenotrophic-
  leaning — together with the **Methanobacteriaceae** dominance. This is the
  load-bearing isotopic evidence. (Optionally note the concentration-weighted
  source estimate, −65 to −79‰, agrees, confirming a depleted internal source
  distinct from atmosphere.)
- **Now add δ¹³C-CO₂ ≈ −20‰** (this is what R2 asked for — the source CO₂ WAS
  measured) and the **within-tree apparent ε_C**: ~44‰ (bulk), rising to ~55‰
  as internal CH₄ increases and the signal sheds atmospheric dilution/oxidation.
- **Do NOT report** the pooled-source ε_C, and do NOT report the tissue (H/S)
  enrichment.

## What NOT to over-claim
- ε_C ≈ 44‰ (bulk) is *below* the classic ~55‰ CO₂-reduction line; it only
  reaches the boundary in the cleanest high-CH₄ samples. So **ε_C is
  corroborating, not decisive** for hydrogenotrophy.
- ε_C uses **bulk internal CO₂ (respiration-dominated)**, not the isolated
  methanogenic substrate — so it is an *apparent* fractionation. State this
  plainly (it is exactly R2's point).

## HOW HARD TO LEAN ON IT (Jon's call): light touch.
One line in text/SI + a response-to-referee entry. Do NOT make ε_C a centerpiece.

### One-line for text (SI methods/results or a clause in the isotope paragraph)
> Paired internal δ¹³C-CO₂ (median −20‰) implies an apparent CO₂–CH₄ carbon
> isotope fractionation (ε_C ≈ 44‰, rising to ~55‰ in higher-CH₄ samples)
> consistent with a hydrogenotrophic contribution, though bulk internal CO₂ is
> respiration-dominated rather than the isolated methanogenic substrate, so this
> fractionation is corroborative rather than diagnostic.

### Response to Referee 2 (L607-610) — final
> We thank the reviewer. The Picarro G2201-i recorded δ¹³C-CO₂ simultaneously with
> δ¹³CH₄, which we now report (median −20‰). Pairing the two within each tree gives
> an apparent CO₂–CH₄ fractionation of ε_C ≈ 44‰, rising to ~55‰ in samples with
> internal CH₄ high enough to minimize atmospheric dilution and oxidation (≥10 ppm;
> Picarro precision is <1.15‰ at 10 ppm). These are consistent with a substantial
> hydrogenotrophic (CO₂-reduction) contribution and reinforce the strongly depleted
> δ¹³CH₄ (−63.7‰; atmosphere-corrected tree-produced source ~−70‰) and the
> dominance of hydrogenotrophic Methanobacteriaceae. As the reviewer notes, the
> δ¹³C of the source CO₂ matters: we now state explicitly that the bulk internal
> CO₂ pool is respiration-dominated rather than the isolated methanogenic substrate,
> so we interpret ε_C as corroborating rather than diagnostic and have tempered the
> isotopic language accordingly (Whiticar 1999; Conrad 2005).
