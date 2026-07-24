# Plant traits & the stem methane cycle — SI figure + discussion (Editor + R3 #3.1)

**Clean deliverable** (the working log with all numbers/corrections is
`plant_traits_bridge_finding.md`). Framing agreed with Jon: **one SI figure + one
discussion paragraph** that answer the editor's plant-framing note honestly. Uses only
this study's data + per-species trait values; not central; does not import the separate
cross-site concentration paper. n ≤ 10 species — descriptive / hypothesis-generating.

## SI figure — one heatmap
`traits_heatmap_robust.png` (script `R_traits_heatmap_robust.R`). A single
publication-style heatmap: per-species traits (rows, grouped Structure / Chemistry /
Roots / Moisture / Whole-plant) × the four methane-cycling quantities (columns, ordered
along the microbial chain mcrA → methanotroph → balance → net flux). Cell = Spearman ρ.
Two markers carry the statistics without prose: a **graded significance ladder**
(`*` p<0.10, `**` p<0.05, `***` FDR<0.10, `****` FDR<0.05, Benjamini–Hochberg across all
cells) and a **bold outline** for associations that survive control for the
gymnosperm/angiosperm split (rank-partial). This replaces the earlier 3-part sequence
(heatmap + ordination + breakdown); the ordination and breakdown are available
(`traits_pca_biplot.png`, `traits_procrustes.png`) if a reviewer wants them, but the one
heatmap now folds in the robustness the breakdown used to show.

### Figure caption (drop-in)
> **Figure SX. Per-species plant traits and the stem methane-cycling community.** Spearman
> rank correlations between per-species functional traits (rows, grouped by category) and
> four species-level methane-cycling quantities (columns, ordered along the microbial chain:
> methanogen *mcrA* abundance → methanotroph (*pmoA*+*mmoX*) abundance →
> methanogen:methanotroph balance → net stem CH₄ flux), for the ten tree species with ≥5
> sampled trees and ≥5 flux measurements. Cell colour and value give Spearman ρ. Statistical
> support is shown as a graded ladder — \*, p<0.10; \*\*, p<0.05; \*\*\*, FDR<0.10; \*\*\*\*,
> FDR<0.05 (Benjamini–Hochberg across all cells). No cell reaches \*\*\*: **no association
> survives multiple-comparison correction at this species-level sample size**, so all are
> exploratory. A bold outline marks associations that remain (raw p<0.05, |ρ|≥0.3) after
> controlling for the gymnosperm/angiosperm split (rank-partial correlation), distinguishing
> clade-robust signals (plant longevity, realized soil moisture, bark density → balance) from
> those driven by the coniferous-vs-angiosperm contrast (wood/bark density → methanotroph).
> "Realized soil moisture" is each species' median volumetric water content across its stems
> at our site; all other traits are from external databases (see Methods).

## SI methods (drop-in)
> We compiled per-species functional traits — wood and bark density (GWDD), wood
> porosity/type (USDA), sapwood and heartwood pH, and TRY-derived anatomy, chemistry, and
> whole-plant traits — and each species' realized soil-moisture niche (median volumetric
> water content across its stems at our site), and related them to species-level stem CH₄
> flux and to area-weighted *mcrA*, methanotroph (*pmoA*+*mmoX*), and
> methanogen:methanotroph gene abundances (Spearman correlations). Because a trait varies
> only among species, species level is the correct unit; we confirmed the key effects on
> individual trees with mixed models [response ~ trait + (1|species)] rather than
> pseudoreplicated per-tree regressions. Given n≤10 species we report Benjamini–Hochberg
> false-discovery rates alongside raw p-values, and — because the two coniferous species are
> trait-space extremes — we tested robustness to the gymnosperm/angiosperm split with
> rank-partial correlations (controlling for a gymnosperm indicator) and with
> angiosperm-only subsets. We treat all trait–response associations as descriptive.

## Discussion paragraph (drop-in, ~200 words)
> To ask whether plant species properties account for among-species variation in stem CH₄
> flux and its microbial basis, we related per-species wood and whole-plant functional
> traits, and each species' realized soil-moisture niche, to species-level flux and gene
> abundances (Fig SX). Directional structure is apparent — denser-stemmed, longer-lived
> species harbour fewer methanotrophs and a more methanogen-tilted balance, and
> wetter-growing species carry more methanogens. This structure is exploratory: with ten
> species, no association survives false-discovery-rate correction. Probing robustness to the
> deepest phylogenetic split, the wood- and bark-density associations with methanotroph
> abundance prove to be largely a coniferous-versus-angiosperm contrast — they do not persist
> once the two conifers are controlled for. What does persist are a whole-plant **longevity**
> signal (longer-lived species show a more methanogen-tilted balance and fewer methanotrophs,
> a relationship also present among angiosperms alone) and the realized soil-moisture–
> methanogen relationship. Net stem CH₄ flux itself is only weakly related to any single
> trait, and no trait combination predicts it in cross-validation — as expected for the noisy
> integration of production, oxidation, and gas transport, which also accounts for the modest
> variance our flux models explain. We therefore present these as descriptive,
> hypothesis-generating patterns that implicate plant life-history and hydrology as controls
> on the stem methane-cycling community, while most among-species variation — and net flux in
> particular — remains unexplained: a priority for future, more highly-replicated trait work.

## Response to Editor / Referee 3 (#3.1)
> We agree that species-level controls deserved fuller treatment. We now integrate
> per-species wood and whole-plant functional traits, and each species' realized soil-
> moisture niche, with our flux and gene data, presented as a new SI figure and a discussion
> paragraph (Fig SX). Rather than assert a trait–flux mechanism, we report the associations
> transparently with both raw and false-discovery-rate-corrected significance and a test for
> robustness to the gymnosperm/angiosperm split. Directionally sensible structure is present
> (longer-lived, denser-stemmed species have a more methanogen-tilted balance; wetter-growing
> species carry more methanogens), but none survives multiple-comparison correction at ten
> species, and the apparent wood-density effects on methanotrophs are largely the two conifers
> (a functional-type contrast, not a continuous gradient); the clade-robust signals are
> whole-plant longevity and realized soil moisture. Net flux is not trait-predictable in
> cross-validation. We therefore frame plant traits as capturing only a coarse life-history
> and hydrological part of among-species variation, with most of it — and net flux especially —
> remaining unexplained, and identify this as a priority for future, more highly-replicated
> work. This candid treatment also speaks directly to the reviewer's related point that a
> substantial fraction of flux variance is not captured by our models.

## Scripts behind the figure
`R_traits_heatmap_robust.R` (the SI figure) · supporting/exploratory:
`R_traits_flux_microbial_bridge.R` · `R_traits_multitrait.R` · `R_traits_dimreduction.R`
(PCA/PLS) · `R_traits_dimreduction2.R` (RDA/Mantel/Procrustes) · `R_traits_mixed_robustness.R`.
Superseded by the single heatmap: `R_traits_heatmap.R`, `R_traits_breakdown_figure.R`
(the breakdown's "nothing significant among angiosperms" claim was wrong — plant longevity is,
p=0.005; the robust heatmap corrects and replaces it).
