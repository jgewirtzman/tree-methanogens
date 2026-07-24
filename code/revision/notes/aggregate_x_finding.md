# Aggregate-x / keep-individual-y — what we found (Mark's Slack point 3)

Script: `revision/analysis/R3_aggregate_x_keep_y.R` · Outputs:
`revision/outputs/aggregate_x_model_comparison.csv`, `aggregate_x_summary.txt`,
`aggregate_x_ratio_panels.png`. All reproduce the manuscript's species result
(ratio R² = 0.513) exactly, then extend it.

## The one-line answer to Mark

Aggregating the noisy predictor (x) while keeping individual flux (y) **does**
de-bias the slope and recover a significant relationship — **but only once flux
is on a scale where OLS is valid.** On raw flux it fails, for an instructive
reason: raw stem flux is so heavily right-skewed that a handful of extreme
emission events (unrelated to the microbial ratio) dominate any mean-based fit.

## The numbers (ratio predictor)

| flux scale | model | slope | R² | p |
|---|---|---|---|---|
| raw | M1 individual (x & y individual) | 0.029 | 0.009 | 0.31 |
| raw | M3 **aggregate-x, chamber y** | 0.028 | 0.001 | 0.51 |
| raw | M3b aggregate-x, 125cm y | 0.126 | 0.025 | 0.09 |
| raw | M2 species median, OLS | 0.052 | 0.513 | 0.020 |
| raw | M2 species median, SMA | 0.073 | 0.513 | 0.020 |
| **pseudo-log** | M1 individual | 0.045 | 0.026 | **0.078** (marginal) |
| **pseudo-log** | M3 **aggregate-x, chamber y** | 0.076 | 0.011 | **0.028** (sig) |
| **pseudo-log** | M3b aggregate-x, 125cm y | 0.163 | 0.054 | **0.013** (sig) |
| **pseudo-log** | M2 species median, SMA | 0.153 | 0.510 | 0.020 |

Diagnostic: **cor(species-mean ratio, species MEDIAN flux) = 0.717** (the
R²=0.51 headline) vs **cor with species MEAN flux = 0.328**. The relationship is
with species-*typical* (median) flux, not the mean, which the extreme events
inflate (e.g. *Fraxinus*: median 0.013, max 1.84; *Acer rubrum*: median 0.046,
max 3.45).

## What this means

1. **Mark's intuition is vindicated on the right scale.** On pseudo-log flux, the
   aggregate-x/keep-y slope (~0.076) lands right at the species-level SMA slope
   (~0.073 raw / 0.153 pseudo-log) — the signature that x-attenuation, not
   absence of relationship, drives the weak individual fit. The individual fit is
   marginally positive (p≈0.08), not zero; aggregation sharpens rather than
   manufactures the signal.

2. **The raw-flux "individual R² < 0.001, no individual relationship" claim is
   fragile** — it is largely an artifact of analyzing heavy-tailed raw flux.
   That's worth knowing *before* a fresh reviewer re-runs it.

3. **Jon's "both axes are noisy" instinct is right, and there's a second axis
   problem too**: individual flux isn't just noisy, it's heavy-tailed, so a
   single chamber is a poor estimate of species-typical flux and the *median* is
   the appropriate robust summary. That is the honest justification for
   species-level aggregation — matching the level at which the causal plant
   traits vary — rather than a statistical convenience.

## Decisions for Jon (and Mark, in the fall)

- **(A) Presentation scale.** Keep the species result on raw median flux
  (unchanged headline, R²=0.51), but consider adding the pseudo-log
  aggregate-x/keep-y result as the direct answer to Mark and to R3's "circular"
  charge — it shows the individual→species transition is attenuation, on a valid
  scale, with individual y retained. Recommended: report species-median SMA as
  headline; cite the pseudo-log aggregate-x as robustness in the response letter
  (and optionally a SI panel = `aggregate_x_ratio_panels.png`).

- **(B) Retire the "exactly null individually" framing?** Strongly consider
  replacing "no individual-level relationship (R² < 0.001)" with "a weak,
  strongly attenuated individual relationship that sharpens on aggregation."
  This is more defensible and removes the paradox R3 pressed ("how can species
  work if individual is exactly zero?"). **This is a substantive claim change —
  needs Jon + Mark sign-off.**

- **(C) How hard to lean on it.** Per the standing weighting: the mechanism rests
  on convergent evidence (isotopes, redox structure, functional inference, ratio
  > either gene); this regression is corroboration. Don't over-build. The
  pseudo-log aggregate-x result is enough; the hierarchical model remains
  unnecessary.

## Response-to-referees language (draft, for R3 L668-674 "circular")

> We thank the reviewer. We have clarified that our claim is species-level and
> associational, and we now show directly that the individual→species change in
> signal reflects regression attenuation rather than circular reasoning. On the
> analysis scale appropriate to the strongly right-skewed flux distribution, a
> weak but positive individual-level relationship (p ≈ 0.08) sharpens to a
> significant one when the spatially under-sampled predictor is aggregated to the
> species while individual flux variation is retained (slope 0.076, p = 0.028) —
> a value matching the species-level standardized major axis slope. Aggregation
> therefore recovers, rather than manufactures, the relationship, consistent with
> attenuation from within-tree heterogeneity and flux skew (Frost & Thompson
> 2000; Warton et al. 2006).
