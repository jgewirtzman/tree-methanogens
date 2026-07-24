# Plant traits → stem CH₄ flux and the microbial balance (Editor + R3 #3.1)

> **UPDATE / CORRECTION (2026-07) — read first; supersedes parts of this log.**
> The final SI figure is now a SINGLE heatmap, `traits_heatmap_robust.png`
> (`R_traits_heatmap_robust.R`), which folds the robustness in: `*` p<0.10, `**` p<0.05,
> `***` FDR<0.10, `****` FDR<0.05 (BH), + bold outline = survives control for the
> gymnosperm/angiosperm split (rank-partial). It replaces the heatmap+ordination+breakdown
> sequence. Clean deliverable + caption: `plant_traits_SI_and_discussion.md`.
> Two substantive corrections to the analysis below: (1) the claim "**nothing significant
> among angiosperms**" (breakdown figure) is WRONG — **plant longevity** is significant among
> angiosperms (→methanotroph & →balance, p=0.005). (2) The clade-robust signals (survive the
> gymnosperm control) are **plant longevity → balance +0.67 / methanotroph −0.64**, **realized
> soil moisture → mcrA +0.72**, and **bark density → balance +0.73** — while density →
> methanotroph is NOT clade-robust (the conifer contrast, as this log already noted). Reframe
> the discussion around **longevity + moisture**, with density as a correlated-but-clade-
> confounded axis. FDR still kills everything at n=10 (min q≈0.5 full grid; PC axes don't
> rescue it) — descriptive only. The prose in this log below predates these corrections.

**Scope: SI figure(s) + a discussion paragraph to address the editor's plant-framing
note — NOT a repositioning of the paper.** Editor De Kauwe and R3 asked what *about a
species* relates to stem CH₄ flux / the microbial community. This adds a bounded,
honest answer using only this study's YMF data + per-species trait values. It does
**not** bring in the separate cross-site concentration paper, and it is not pitched as
a central result. Deliberately modest: descriptive structure + a couple of a-priori,
directional hypotheses, offered for the discussion.

Scripts: `R_traits_flux_microbial_bridge.R` (core), `R_traits_multitrait.R` (expanded),
`R_traits_heatmap.R` (descriptive figure). Start here: **`traits_heatmap.png`** — the
plain descriptive picture of how traits structure the methane-cycling community.
Outputs: `traits_bridge_*`, `traits_multitrait_*`, `traits_heatmap*`.

## Data & method
Per-species wood/plant traits from the sibling project `../tree-gas-traits` (TRY +
GWDD wood density + BAAD heartwood fraction + curated wood pH + USDA porosity + roots +
whole-plant), plus **our own realized soil-moisture niche** per species (`vwc_mean` from
`tree_species_moisture_niche.csv`), joined by species code to *this* paper's canonical
species-level flux and gene variables (`_prep_species_data.R`, which reproduces the
reviewed ratio→flux **exactly: r=0.717, R²=0.513** — the earlier "0.47" was a transform
artifact in a check, not a data discrepancy). Species-level, n = 10 (≥5 trees & ≥5 flux).

## The honest statistical situation — read this first
- The full exploratory atlas is **20 traits × 4 responses = 80 tests**. **Nothing survives
  FDR correction** (all q ≈ 0.61). At n=10 a broad scan is far too underpowered to yield
  corrected significance, so we do **not** present it as "many trait associations."
- We therefore lead with (i) **out-of-sample cross-validation** (LOO-CV R², which is not a
  p-value scan and directly measures generalization), and (ii) **pre-specified mechanistic
  hypotheses** grounded in theory + the well-powered (n=29) sibling concentration analysis —
  not with cherry-picked correlations. The atlas is reported in full for transparency and to
  guide future work, explicitly flagged as uncorrected/exploratory.

## What actually holds up (LOO-CV, out-of-sample)
| Relationship | LOO-CV R² | reading |
|---|---|---|
| Realized soil moisture → **methanogen (mcrA)** abundance | **0.23** | wetter niche ⇒ larger production pool |
| Wood/bark density → **methanogen:methanotroph balance** | **0.20–0.29** | denser stems ⇒ methanogen-tilted balance |
| (bark density → **methanotroph** abundance — what tilts the balance) | — (ρ=−0.73, p=0.03) | denser bark ⇒ fewer consumers |
| Any single trait → **net stem CH₄ flux** | **≤0.13** | flux is the noisy downstream sum |

**Two largely independent plant controls, each cross-validated:**
1. **Moisture → production.** Species growing in wetter soil (our realized VWC, *independent*
   of wood density: collinearity 0.19) carry more methanogens. Moisture drives the *pool*.
2. **Structure → balance.** Denser wood and bark are associated with fewer aerobic methanotrophs
   (bark→methanotroph ρ=−0.73, p=0.03; wood→methanotroph ρ=−0.65, p=0.05) and a more
   methanogen-tilted balance (bark→balance ρ=+0.73, p=0.03). Structure drives the *balance*.
   **Mechanism is OPEN — and it is NOT a simple internal-O₂ story.** We checked: wood
   density→O₂-deficit ρ=−0.43 (wrong sign for an O₂-barrier), and O₂-deficit→methanotroph
   ρ=+0.12 (n.s.). So denser stems are not measurably more anoxic, and anoxia does not track
   methanotroph loss. The density–methanotroph link is a robust association whose driver
   (wood chemistry, water content, antifungal compounds, substrate, pH?) we do not claim to
   have identified. (An earlier draft asserted an O₂-diffusion-barrier mechanism; the internal-gas
   data do not support it, so it is withdrawn.)

**But net flux itself is only weakly trait-predictable (LOO-CV R²≤0.13; multi-trait flux models
overfit to negative LOO).** This is the key honest point — and it *reframes R3's own critique*:
the microbial *components* (production pool, oxidation balance) are trait-structured, but net flux
is their noisy sum plus transport, which is exactly why the manuscript's RF explained only ~15% of
flux. Components predictable, sum noisy — biology, not a modeling failure.

## The informative NULLS (these answer R3 directly and are robust to multiple comparisons)
- **Rooting depth → flux: ρ=−0.19, p=0.61** — answers R3's "could this be rooting depth?": no.
- **Realized soil moisture → flux: ρ≈−0.02, p=0.97** — wetter-growing species are **not** higher-flux,
  even though they carry more methanogens. Directly engages the R1↔R3 birch/wet-soil tension:
  moisture raises the production pool but not net emission (consumption/transport intervene).
- **Wood porosity → flux: ρ≈0.06, p=0.88** — porosity, the top predictor of internal *concentration*
  in the sibling atlas (R²=0.41, n=29), does **not** predict realized flux. Different traits govern
  *accumulation* vs *net emission*.

## Dimension reduction — the robustness stress-test (and it tempers the claim)
`R_traits_dimreduction.R` (PCA + supervised PLS with LOO-CV). Two honest results:
- **PCA:** the dominant trait axis (**PC1, 35% variance**) is a wood-economics gradient
  (gymnosperm/tall/soft ↔ dense/ring-porous hardwood). The methane responses projected onto
  trait space point in the *mechanistically sensible* directions (flux & balance toward high
  wood/bark density; methanotroph opposite; mcrA toward realized VWC) — so the descriptive
  *directions* corroborate the heatmap. **But no PC significantly correlates with any response**
  (only a minor PC3→balance, one of 12 tests), and species do not separate cleanly by flux.
- **PLS (traits → response, LOO-CV):** **no trait combination predicts any response out of
  sample** — CV R²: flux **−0.29**, balance **−0.05**, mcrA **−0.02**, methanotroph **+0.07**.
  At n=10 with ~15 collinear traits, multivariable trait models do not generalize.
- **Permutation multivariate co-structure** (`R_traits_dimreduction2.R`): RDA — traits (PC1–3)
  explain 54.5% of joint response variance (adj R²=0.32) but the omnibus permutation is n.s.
  (F=2.4, **p=0.10**; only trait-PC2 significant, p=0.03); Mantel trait-dist vs response-dist
  r=0.19, **p=0.17**; Procrustes ordination alignment r=0.56, **p=0.10**. All three agree:
  moderate but **not statistically significant** joint co-structure at n=10 — matching the failed
  PLS. (t-SNE/UMAP deliberately NOT used — indefensible at n≤16; they manufacture artifacts.)

**What this means (and how to pitch it):** the plant-trait signal here is **descriptive and
directional, not a predictive multivariable model.** Only *simple, a-priori single-trait*
relationships have any out-of-sample support (VWC→mcrA single-predictor LOO R²≈0.23;
density→methanotroph≈0.20–0.29). This is exactly the right size for an **SI + discussion** add:
we show the structure, name two directional hypotheses, and explicitly do not claim a predictive
trait model. Appropriately modest for a point that is not central to the paper.

## IMPORTANT correction — the density effect is mostly the 2 conifers, not a trait gradient
`R_gymnosperm_sensitivity.R`. The two gymnosperms (PIST, TSCA) are the Procrustes misfits and
trait-space extremes (wood-density z = −4.4, −2.3) with relatively **high methanotrophs**
(z = +1.5, +0.55) and a **methanotroph-tilted balance** (z = −1.2, −0.4). Dropping them collapses
almost every density relationship:

| Relationship | ALL species | ANGIOSPERM-only |
|---|---|---|
| wood density → methanotroph | ρ=−0.48 (p=0.07) | −0.38 (p=0.20) |
| bark density → methanotroph | −0.47 | −0.23 (p=0.50) |
| bark density → balance | +0.26 | −0.05 |
| wood density → flux | +0.12 | −0.13 |
| **moisture → methanogen (mcrA)** | +0.42 | **+0.37 (holds)** |

**So the "wood/bark density → methanotroph suppression → balance" pattern is largely the
categorical conifer-vs-hardwood difference, NOT a continuous density gradient** — it does not hold
within angiosperms. This corrects the earlier claim (and the mixed-model p-values of 0.046/0.017,
which also used all species, are likewise conifer-driven). The **one continuous effect that survives
angiosperm-only is moisture → methanogen abundance.** Net flux is unrelated to traits at every level.

**Modeling implication:** (1) present moisture→methanogen as the defensible continuous relationship;
(2) demote density→methanotroph/balance to a *tentative conifer/hardwood contrast*, explicitly n=2
conifers and hypothesis-only; (3) do not code wood type as ordinal `porosity=1/2/3` (tracheid is a
categorical architecture, not a step below diffuse-porous) — treat clade/wood-type categorically and
gymnosperm as an outlier flag, not a modelled group (impossible at n=2).

## Aggregation: species-level vs individuals (robustness) — `R_traits_mixed_robustness.R`
Unlike the microbe→flux relationship (both vary within a species), a **trait has no
within-species variation** — it varies only between species. So species-level is the natural,
correct level, and a naive per-tree OLS would be **pseudoreplication**. The valid "use the
individuals" version is a mixed model `response_ij ~ trait + (1|species)` ("aggregate-x/keep-y"
done right). We ran both. The mixed model's Satterthwaite df come out **~5–9 (species-level, not
tree-level)**, so its p-values are honest, and the key effects **hold on individuals**:
wood density → methanotroph slope=−0.23, **p=0.046**; bark density → balance slope=+0.35,
**p=0.017**. So the density→consumer/balance signal is **not an aggregation artifact** (the mixed
model is if anything more supportive than species-medians). Stem flux stays null at both levels;
the VWC→methanogen effect is sample-dependent (stronger in the strict n=10 set, marginal at n=14).

## Caveats to state in the paper
- n=10 species; the exploratory atlas does not survive FDR. We claim only the CV-supported,
  a-priori relationships and present the rest as hypothesis-generating.
- Plant height→flux (ρ=−0.72) is **collinear with wood density** (−0.50) — not an independent effect;
  we fold it into the density axis rather than claim it separately.
- Absolute ddPCR loads are subject to the pending ×10 issue; the **balance is a ratio** and every
  correlation/CV here is **unaffected**.

## How to use it in the paper (SI + discussion — bounded)
A short discussion paragraph + 1–2 SI figures (heatmap; optionally the PCA biplot). The point to
make: species differences relate to the methane-cycling community along **two directional axes** —
(1) the soil-moisture niche a species occupies tracks its methanogen pool; (2) wood/bark **structure**
tracks the methanotroph-driven balance (mechanism unresolved — not simply internal O₂). Net flux is
the noisy integration of both, which is *why* single traits and even a multivariable model do not
predict it — consistent with the modest variance our flux models explain. This is enough to answer
the editor's plant-framing note and speak to R3's questions (it also shows *rooting depth* and
realized moisture do not predict flux) **without repositioning the paper** around traits.

## Response to Editor / Referee 3 (#3.1) — draft (discussion-scoped)
> To speak to how plant species properties relate to stem CH₄ dynamics, we compiled per-species wood
> and whole-plant functional traits (wood/bark density, wood pH, porosity, rooting depth, and related
> anatomy) together with each species' realized soil-moisture niche at our site, and examined their
> association with our flux and gene-abundance data (SI Fig SX). Two directional patterns emerge: a
> species' realized soil-moisture niche is associated with its methanogen abundance, while wood and
> bark density track the methanogen:methanotroph balance via fewer methanotrophs in denser-stemmed
> species. Net stem CH₄ flux itself is only weakly related to any single trait — and no trait
> combination predicts it in cross-validation — as expected for the noisy integration of production,
> oxidation, and transport, which also accounts for the modest variance our flux models explain.
> Notably, rooting depth and realized soil moisture do not predict net flux. We present these as
> descriptive, hypothesis-generating relationships (they do not survive multiple-comparison correction
> at this species-level sample size) that point to plant-structural and hydrological controls on the
> stem methane-cycling community as a target for future, more highly-replicated work.
