# Aggregation / RMA — draft text for the revision

Addresses **R3 L664** ("Averaging across individuals in a species… is there any
precedent for this?"), **R3 L668–674** ("This seems circular…"), and **Mark's**
ecological-correlation-vs-regression-dilution point.

**Strategy (per the thread):** don't lean on the individual-level null as
evidence; make the positive argument that (i) both core and single-position
chamber are spatially under-sampled point observations of a *species-level*
trait, so error lives on **both axes**; (ii) the claim is explicitly
species-level and associational, not individual-causal; (iii) the mechanism
rests on convergent evidence (isotopes, redox spatial structure, functional
inference, ratio > either gene alone), with the species regression as
*corroboration, not the load-bearing proof*. Add **one** SMA/RMA robustness
sentence. Do **not** build the hierarchical Bayesian model.

All numbers below are from `revision/outputs/` (script
`revision/analysis/R2_species_aggregation_RMA.R`), which reproduces the
as-reviewed species aggregation exactly (ratio R² = 0.513, r = 0.717).

---

## 1. Reframe the variance-partition sentence (currently L149)

> **Current:** "Variance partitioning showed that species identity explained
> only 5.3% of variance, species-environment interactions 8.7%, environmental
> factors alone <0.01%, leaving 82.9% unexplained at the individual tree level."

**Proposed (turn the 82.9% from confession into justification):**

> Variance partitioning showed that species identity explained only 5.3% of
> variance and environmental factors <0.01%, with 82.9% unexplained at the
> individual-tree level. This large within-species, within-tree component is the
> expected consequence of characterizing a spatially heterogeneous process from
> single-point measurements: a single chamber samples one patch of bark, and a
> single core <0.1% of stem volume (Fig. 7), so each is one draw from a broad
> within-tree distribution rather than a stable estimate of that tree. It is
> precisely this fine-scale heterogeneity — not an absence of biological signal —
> that motivates aggregating replicate observations to the species level, where
> the characteristic traits that structure methane cycling (wood anatomy,
> chemistry, heartwood formation) actually vary (see Discussion).

---

## 2. Replace / expand the aggregation-defense paragraph (currently L261)

The current text (L261) makes the right gestures (Levin 1992; Polussa et al.
2021; Jasienski & Bazzaz 1999; numerator/denominator terms) but deploys them
*defensively* as caveats. R3 read past it. Invert it: lead with the principle,
then apply it.

**Proposed replacement paragraph:**

> Methanotrophs compound this complexity through their own spatial heterogeneity
> and functional diversity. Their abundance also failed individually to predict
> flux, yet the species-level methanogen:methanotroph ratio was the strongest
> predictor of net flux (R² = 0.51; Fig. 8), outperforming methanogen abundance
> alone (R² = 0.39). We interpret this as a **species-level, trait-based
> relationship rather than an individual-tree causal one**, and the aggregation
> is a deliberate analytical choice rather than a post hoc rescue. Both of our
> individual-level measurements are spatially under-sampled point estimates of an
> underlying species characteristic: a single increment core captures <0.1% of
> stem volume, and a single chamber position integrates one patch of a stem whose
> internal methanogen abundance varies by >3 orders of magnitude within one tree
> (Fig. 7). With error on **both** the predictor and the response at the
> individual scale, individual-level correlations are attenuated toward zero
> (regression dilution; Frost & Thompson 2000), whereas averaging replicate
> observations to the species mean reduces that error and recovers the
> underlying relationship — aggregation cannot manufacture a correlation from a
> genuinely absent one, so its emergence on aggregation is itself evidence that a
> real, scale-dependent relationship exists (Levin 1992). This is the standard
> unit of analysis in trait-based and comparative ecology, where species-mean
> values are regressed against one another precisely because the causal traits —
> here wood chemistry, anatomy, and heartwood formation — are species-level
> properties. We therefore make the explicit and limited claim that species
> characterized by a higher production:consumption balance emit more, not that
> gene abundance predicts flux in any individual tree. Consistent with a genuine
> species-level signal rather than a small-sample or ratio artifact, the ratio
> outperformed either gene alone by AIC (−44.3 vs. −42.1), held with numerator
> and denominator entered as separate terms (Fig. 8; Jasienski & Bazzaz 1999),
> and was robust in leave-one-species-out jackknifing (slope 0.066–0.077 across
> all folds; significant in 9 of 10, and *strengthened* when the highest-leverage
> species was removed). Because both variables are estimated with error, we
> verified the relationship under a Model II standardized major axis (SMA/RMA)
> regression, which accounts for error on both axes and, as expected under
> attenuation, yields a steeper slope (SMA 0.073 vs. OLS 0.052; 95% CI
> 0.042–0.125) without altering its significance. Critically, the mechanistic
> interpretation does not rest on this regression alone but on its convergence
> with independent lines of evidence — depleted δ¹³CH₄, the inverse
> heartwood/sapwood spatial structuring of methanogens and methanotrophs,
> functional co-occurrence of fermentative and hydrogenotrophic pathways, and the
> ratio outperforming either gene — which together indicate that net flux
> reflects a within-stem balance of production and oxidation.

---

## 3. One-sentence Methods addition (SMA robustness)

Add to the species-level statistics methods (SI Methods S2 or main text):

> Because species-mean gene abundance and species-mean flux are both estimated
> with sampling error, we additionally fit the species-level gene–flux
> relationships with a Model II standardized major axis (SMA, equivalent to
> reduced major axis) regression, which is appropriate when both variables carry
> error and no strict predictor→response direction is assumed (Warton et al.
> 2006). SMA slopes and 95% confidence intervals were computed analytically
> (SMA slope = OLS slope / |r|); the correlation and its significance are
> unaffected by estimator choice.

**Reference to add:** Warton DI, Wright IJ, Falster DS, Westoby M. 2006.
Bivariate line-fitting methods for allometry. *Biological Reviews* 81: 259–291.

---

## 4. Response-to-referees language (concise)

**To R3 (L664):** "Species-mean values are the standard unit of analysis in
trait-based and comparative ecology, because the traits that plausibly drive
stem methane cycling — wood chemistry, anatomy, heartwood formation — are
species-level properties (Levin 1992). We now state this explicitly as a
deliberate choice (Discussion, p. X) and clarify that our claim is
species-level and associational, not a prediction for individual trees."

**To R3 (L668–674, "circular"):** "We agree the individual-level null should not
be used as evidence, and we have rewritten the passage to make the *positive*
argument: both core and single-position chamber are under-sampled point
estimates of a species trait, so the individual relationship is attenuated by
regression dilution (Frost & Thompson 2000), and aggregation recovers rather
than manufactures signal. We now report that the result is robust to
leave-one-species-out jackknifing and to a Model II (SMA/RMA) estimator that
accounts for error on both axes, and we make explicit that the mechanism rests
on convergent independent evidence rather than this regression alone."

**Optional, for Mark:** we can additionally present the "aggregate-x, retain-y"
model he suggested as a robustness row, but note that because the chamber flux
is itself a one-position under-sample, retaining individual y keeps a large
block of unmatched within-tree heterogeneity — which is why symmetric
aggregation is the more principled treatment here. (Not run yet; flagged in the
script as a possible add if Mark wants it.)

---

## Numbers appendix (from `revision/outputs/rma_summary.txt`)

| predictor | n | r | R² | p | OLS slope | SMA slope | SMA 95% CI | sig? |
|---|---|---|---|---|---|---|---|---|
| **ratio** | 10 | 0.717 | 0.513 | 0.020 | 0.052 | 0.073 | 0.042–0.125 | yes |
| mcrA (area-wt) | 10 | 0.628 | 0.394 | 0.052 | 0.098 | 0.156 | 0.086–0.283 | borderline |
| methanotroph | 10 | −0.480 | 0.230 | 0.161 | −0.042 | −0.087 | — | no |

- SMA correction factor = 1/|r|: ratio ×1.40, mcrA ×1.59, methanotroph ×2.09 —
  the OLS→SMA gap tracks 1/|r|, the signature of attenuation; the strong
  relationships move least.
- Jackknife (ratio): R² 0.333–0.725, SMA slope 0.066–0.077, significant 9/10
  (only non-sig dropping *A. saccharum*, p = 0.104); dropping highest-leverage
  *Q. rubra* → R² 0.725.
- **Significance comes from p/r, not the SMA CI** (an SMA slope CI is a
  magnitude interval, essentially never spanning zero). Report p for "is there a
  relationship," SMA for "how steep."
