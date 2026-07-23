# Captions for the two new/reworked main figures (revision)

## Synthesis figure — convergent evidence for hydrogenotrophic, syntrophy-associated methanogenesis
`revision/outputs/fig_hydrogenotrophy.png` (script `R_fig_hydrogenotrophy.R`)

> **Figure X. Convergent evidence for hydrogenotrophic, syntrophy-associated methanogenesis in tree
> heartwood.** Convergent taxonomic, functional-inference, and isotopic evidence. **(a)** Microbial families whose relative
> abundance is significantly associated with heartwood *mcrA* gene abundance (linear mixed-effects
> models controlling for compartment and 16S load; FDR < 0.05), grouped by evidence for feeding
> methanogenesis: methanogens (red; *Methanobacteriaceae*, *Methanomassiliicoccaceae*);
> fermentative associates with evidence of syntrophic support (gold; *Christensenellaceae*,
> demonstrated H₂ transfer to a methanogen; *Dysgonomonadaceae*, community-level association); and
> a co-occurring amino-acid fermenter without an established H₂ link (grey; *Eggerthellaceae*).
> **(b)** FAPROTAX-predicted
> metabolic functions as log₂(heartwood/sapwood) relative abundance (energy/redox functions with mean
> abundance ≥0.1% in a compartment; host-association, pathogen, and phototrophy annotations excluded).
> Anaerobic and H₂-based functions (hydrogenotrophic methanogenesis, dark hydrogen oxidation,
> fermentation, anaerobic chemoheterotrophy, methylotrophy) are enriched in heartwood (red); aerobic
> and oxidant-linked functions (aerobic chemoheterotrophy, methanol oxidation, methanotrophy,
> cellulolysis, nitrate reduction) are enriched in sapwood (blue); colour scales with log₂(HW/SW).
> Dissimilatory denitrification and sulfate/sulfur respiration were near-zero (<0.03%) in both
> compartments. **(c)** MetaCyc pathways
> (PICRUSt2) associated with *mcrA* abundance (linear mixed-effects t-statistic; FDR < 0.05, after
> removing predicted contributions from methanogen-classified ASVs and excluding pathways with >10%
> methanogen-ASV contribution). Pathways were classified by MetaCyc function *a priori*
> (independent of the association) and generic housekeeping pathways excluded; up to four representative
> pathways per category are shown (top by |t|). Anaerobic carbon flow tracks *mcrA* — C1/carbon-
> fixation (red; RuMP cycle, reductive acetyl-CoA, Calvin-Benson) and fermentation/carbohydrate
> degradation (gold) — while oxidant-linked metabolism is depleted: aerobic respiration (dark blue;
> aerobic respiration I, TCA-cycle variants) and sulfur/nitrogen metabolism (light blue; assimilatory
> sulfate and nitrate reduction, sulfate/cysteine and methionine biosynthesis). The joint depletion of
> aerobic respiration *and* nitrate/sulfate reduction where *mcrA* is high is consistent with methanogens
> occupying the most reducing, oxidant-poor microsites (the thermodynamic redox ladder); note these N/S
> pathways are assimilatory, so this is a community redox-state signal rather than direct evidence of
> terminal-electron-acceptor competition (no significant *dissimilatory* denitrification/sulfate-reduction
> signal). See Methods and SI for the full classification. **(d)** Stable carbon isotopic composition of internal stem CH₄
> (n internal samples; kernel density above, individual samples below sized by CH₄ concentration).
> The distribution is ¹³C-depleted, and the atmosphere-corrected tree-produced source (Keeling-plot
> intercept −79‰, 95% CI shown; two-endmember mixing correction −70‰) falls in the hydrogenotrophic
> (CO₂-reduction) range (brackets, Whiticar 1999); the apparent CO₂–CH₄ fractionation ε_C = 54–64‰ is
> likewise consistent with CO₂ reduction. Dashed line, atmospheric δ¹³CH₄ (−47‰).
> Colour scheme: warm = *mcrA*-associated / anaerobic (red = methanogenesis, C1/carbon fixation;
> gold = fermentative associates / fermentation); cool = depleted / oxidant-linked (dark blue = aerobic
> respiration; light blue = sulfur/nitrogen metabolism); grey = co-occurring fermenter without an H₂
> link (panel a); purple (d) = internal CH₄.
> The isotopic signal (d) is an independent geochemical line of evidence; the sequence-based lines
> (a–c) derive from the same 16S/*mcrA* data and are mutually corroborative rather than statistically
> independent. Taxonomy (a) and isotopes (d) are the primary evidence, with the functional predictions
> (b, c) subject to the known limitations of taxonomy-based functional inference (Methods).

## Reworked upscaling figure — plot CH₄ component budget (bounding exercise)
`revision/outputs/fig_budget_maps.png` (script `R_fig_component_budget.R`)

> **Figure Y. Plot-scale CH₄ budget at the upland study site (a bounding exercise for this stand).**
> **(a)** Annual mean soil CH₄ flux (nmol m⁻² s⁻¹, per m² ground) interpolated across the inventory
> plot; the soil is a net sink (blue = uptake). **(b)** Annual mean tree-stem CH₄ flux (nmol m⁻² s⁻¹,
> per m² woody surface) for each inventoried tree, point size proportional to DBH; stems are a net
> source (red = emission). Both maps are random-forest predictions (out-of-bag R² = 0.28 soil, 0.15
> tree) and should be read as model estimates. Scale bar 50 m; arrow points north. **(c)** Four-season
> mean soil (blue) and tree (red) CH₄ flux (per m² ground); soil uptake deepens through summer/autumn,
> tree emission peaks in summer. **(d)** Annual net CH₄ budget (mg m⁻² yr⁻¹, per m² ground): soil
> uptake (−904) dominates; tree emission measured to 2 m stem height offsets 0.14% of it, and a
> scaling scenario applying the measured per-area flux to the full woody-area index (WAI 3.07; Gauci
> et al. 2024, assuming constant emission across all woody surfaces, unconstrained above 2 m) offsets
> up to ~13%. The stand remains a strong net CH₄ sink under all scenarios.
