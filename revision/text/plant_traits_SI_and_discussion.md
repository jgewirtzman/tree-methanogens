# Plant traits & the stem methane cycle — SI walk-through + discussion (Editor + R3 #3.1)

**Clean deliverable** (the working log with all numbers/corrections is
`plant_traits_bridge_finding.md`). Framing agreed with Jon: an **SI figure sequence +
one discussion paragraph** that *walks the reader through* the trait analysis honestly —
apparent structure → how it breaks → most variation not yet explainable. Uses only this
study's data + per-species trait values; not central; does not import the separate
cross-site concentration paper.

## SI figure sequence (the walk-through)
Present as a 3-part SI figure (or three panels a–d of one figure):

1. **Heatmap** — `traits_heatmap.png`. Descriptive Spearman correlations of per-species
   traits (structure, chemistry, roots, realized soil-moisture niche) against the four
   methane-cycling quantities (mcrA, methanotroph, balance, flux). Shows the *apparent*
   structure: denser stems ↔ fewer methanotrophs / methanogen-tilted balance; wetter niche
   ↔ more methanogens.
2. **Ordination** — `traits_pca_biplot.png` (+ optionally `traits_procrustes.png`). Trait-
   space PCA with the methane responses projected on; the responses point in sensible
   directions but no axis significantly tracks a response. The Procrustes superimposition
   (trait vs response ordination; PROTEST r=0.56, p≈0.10) shows the **two conifers as the
   principal misfits**.
3. **Breakdown** — `traits_breakdown.png`. Each association recomputed with all species vs
   angiosperms-only: the density effects on the balance and flux **collapse or flip** once
   the two conifers are removed; only moisture→methanogen (and a weaker density→methanotroph
   direction) persists, and **none is significant among angiosperms (n≤13).**

## SI methods (drop-in)
> We compiled per-species functional traits — wood and bark density (GWDD), wood porosity/
> type (USDA), sapwood and heartwood pH, and TRY-derived anatomy, chemistry, and rooting
> traits — and each species' realized soil-moisture niche (median volumetric water content
> across its stems at our site). These were related to species-level stem CH₄ flux and to
> area-weighted mcrA, methanotroph (pmoA+mmoX), and methanogen:methanotroph gene abundances
> (Spearman correlations; trait-space PCA with responses projected; redundancy analysis and
> Procrustes/PROTEST with 999 permutations; PLS regression with leave-one-out cross-
> validation; and, because a trait is a species-level constant, mixed models
> [response ~ trait + (1|species)] on individual trees rather than pseudoreplicated
> per-tree regressions). We assessed sensitivity to the two coniferous species by repeating
> analyses on angiosperms only. Given n≤16 species, we treat these as descriptive and
> report multiple-comparison-corrected as well as cross-validated results.

## Discussion paragraph (drop-in, ~190 words)
> To ask whether plant species properties account for the among-species variation in stem
> CH₄ flux and its microbial basis, we related a panel of wood and whole-plant functional
> traits, and each species' realized soil-moisture niche, to species-level flux and gene
> abundances (Fig SX). Correlational structure is apparent — denser-stemmed species harbour
> fewer methanotrophs and a more methanogen-tilted balance, and wetter-growing species carry
> more methanogens — and trait space and community space share moderate multivariate
> co-structure (redundancy analysis and Procrustes, both p≈0.1). This structure is, however,
> fragile. It is not statistically significant, no trait combination predicts stem flux in
> cross-validation, and a Procrustes superimposition identifies the two coniferous species as
> the principal outliers. Removing them collapses nearly all of the wood-density associations,
> revealing them to be a coarse conifer–hardwood contrast rather than a continuous trait
> gradient; among angiosperms, only the soil-moisture–methanogen relationship persists in
> direction, and none is individually significant at this sample size. We therefore conclude
> that measured plant traits capture only a coarse hydrological and functional-type signal,
> and that **most among-species variation in stem methane exchange and its microbial basis
> remains unexplained** — an important target for future, more highly replicated trait work.

## Response to Editor / Referee 3 (#3.1)
> We agree that species-level controls deserved fuller treatment. We now integrate per-
> species wood and whole-plant functional traits, and each species' realized soil-moisture
> niche, with our flux and gene data, and present the analysis as a new SI figure and
> discussion (Fig SX). Rather than assert a trait–flux mechanism, we walk through the
> evidence transparently: there is apparent and directionally-sensible trait structure in the
> methane-cycling community (e.g. moisture with methanogens; wood density with the
> methanotroph balance), but it is not statistically robust at our species-level sample size,
> does not predict net flux in cross-validation, and the apparent wood-density effects are
> driven by the two conifers (a functional-type contrast, not a continuous gradient). We
> therefore frame plant traits as explaining only a coarse part of among-species variation,
> with most of it — and net flux in particular — remaining unexplained, and we identify this
> as a priority for future, more highly replicated work. We believe this candid treatment is
> more useful to the field than an over-interpreted trait–flux model, and it directly speaks
> to the reviewer's related point that a substantial fraction of flux variance is not
> captured by our models.

## Scripts behind the figures
`R_traits_heatmap.R` · `R_traits_dimreduction.R` (PCA/PLS) · `R_traits_dimreduction2.R`
(RDA/Mantel/Procrustes) · `R_traits_procrustes_plot.R` · `R_gymnosperm_sensitivity.R` ·
`R_traits_breakdown_figure.R` · `R_traits_mixed_robustness.R`.
