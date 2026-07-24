# R3 L342 — "defend the exponential expression in the Methods"

## What R3 is actually pointing at
There is **no exponential model fit anywhere** in the manuscript (submitted or revised): height
effects use **linear mixed-effects models**, flux uses goFlux (LM/Hutchinson-Mosier), upscaling uses
Random Forest. The only "exponential" is the **stated expectation that soil-transported CH₄ declines
exponentially with height in wetland/flooded trees**, used as the a-priori contrast to our upland
"uniform vertical profile" hypothesis (objectives; the black-oak result "rather than … exponentially
decreasing patterns"; discussion). R3 wants that functional-form claim **justified/cited**, and wants
it clear we are not imposing it on our own data.

(Note: PDF line numbers from the review don't map 1:1 to the working file; this is the only
"exponential" in the manuscript, so it is what L342 refers to. Confirm when finalizing.)

## The defense (two parts)
**1. The exponential form is established for soil-transport in stems — empirical + mechanistic.**
- Empirically, basal/soil-derived CH₄ transported up flooded-tree stems shows approximately
  exponential attenuation with height (Pangala et al. 2013; Barba et al. 2019; Jeffrey et al. 2024;
  Anttila et al. 2024).
- Mechanistically this is expected: methane enters at the stem base and is lost radially through the
  bark/wood along the transport path. A constant *fractional* loss per unit height (first-order,
  Fickian radial diffusion out of the transport conduit) integrates to an exponential decline,
  C(h) = C₀·e^(−kh) — the standard result for advective/diffusive transport with distributed first-
  order loss. So "exponential" is not an arbitrary curve; it is the analytic expectation for a basal
  source attenuated by radial diffusion, which is why the wetland literature reports it.

**2. We do NOT fit an exponential to our data.** We invoke it only as the *predicted signature of a
soil/basal source*. Our observed height effects are tested with assumption-free linear mixed-effects
models (no functional form imposed), and the intensively-sampled black oak is described empirically
(mid-stem peak at 4–6 m, coinciding with heart rot). Finding *no* systematic decline (let alone an
exponential one) in upland stems is therefore evidence *against* a dominant basal-transport source —
the logic only requires that a soil source would produce a monotonic (approximately exponential)
decline, which we do not observe.

## Drop-in Methods/Results clarification (one sentence)
> We use the exponential decline of stem CH₄ flux with height only as the expected signature of a
> basal, soil-transported source — a pattern reported for wetland and flooded trees (Pangala et al.
> 2013; Barba et al. 2019; Jeffrey et al. 2024) and consistent with first-order radial diffusive loss
> along the transport path — and we do not impose this functional form on our data; height effects
> are evaluated with linear mixed-effects models.

## Response to Referee 3 (L342)
> To clarify: we do not fit an exponential model. We refer to the exponential decline of stem CH₄
> flux with height as the established signature of basal, soil-transported methane in wetland and
> flooded trees, which is both reported empirically (Pangala et al. 2013; Barba et al. 2019; Jeffrey
> et al. 2024; Anttila et al. 2024) and expected mechanistically from first-order radial diffusive
> loss of a basal source along the stem (yielding C(h) ≈ C₀e^(−kh)). We now state this explicitly and
> note that our own height effects are tested with assumption-free linear mixed-effects models; the
> absence of a systematic decline in upland stems is our evidence against a dominant soil-transport
> source. We have added the citations and the one-sentence clarification above to the Methods.
