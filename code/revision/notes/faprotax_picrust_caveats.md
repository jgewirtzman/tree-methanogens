# 16S-based functional inference — limitations & corroboration (Referee 2 #3)

Metrics: `revision/analysis/R_faprotax_caveat_metrics.R` → `faprotax_caveat_metrics.txt`
(reproduces the manuscript FAPROTAX chain, same seed, adds the numbers below).

## Stance (per Jon): acknowledge limits, don't over-concede
Taxonomy-based functional inference (FAPROTAX, PICRUSt2) is a **standard, widely-used** approach for
16S data; using it here is unremarkable. We are not retracting it. We (a) state its well-known
limitations plainly and quantify them for our data, (b) stop leaning on **absolute** predicted
fractions, (c) foreground the more robust **between-tissue** contrasts, and (d) note that the key
inferences are **independently corroborated** — and are otherwise flagged as hypothesis-generating.

## Provenance — one brief parenthetical (don't over-explain)
The **tissue-level** (heartwood/sapwood) FAPROTAX predictions were reported in the companion Nature
paper; here Fig S3 **adds the species-level view**, and PICRUSt2 (Figs 6, S4, S5) is new. This is one
SI figure, so a single parenthetical is enough — at L185 or in the Fig S3 caption:
> FAPROTAX analysis of 16S rRNA data (Fig. S3; tissue-level predictions also reported in Arnold &
> Gewirtzman et al. 2025, extended here to the species level) …

No separate methods paragraph needed. (The underlying 16S data overlap is already disclosed at L115.)

## Legitimate limitations (quantified, stated once)
- FAPROTAX assigns a function to only **46.9% of heartwood / 50.9% of sapwood** community abundance;
  the rest is functionally uncharacterised.
- Only **32% of ASVs** resolve to a named genus (FAPROTAX/PICRUSt need genus-level matches).
- Reported percentages are therefore fractions of the **whole** community and are **conservative
  floors**, not complete estimates (e.g. fermentation 16.8% of the community ≈ 36% of the *annotated*
  community).
- Categories are **non-exhaustive and non-additive** (H₂ example below).
=> we do not interpret absolute values as quantitative, and in particular do not read a modest
   fermentation % as evidence that few taxa ferment.

## The robust, more informative signal: BETWEEN-TISSUE contrasts
Absolute fractions are noisy, but the heartwood-vs-sapwood *differential* is internally coherent and
matches independent data:

| Function | HW% | SW% | HW/SW |
|---|---|---|---|
| hydrogenotrophic methanogenesis | 3.21 | 0.10 | **32×** |
| dark H₂ oxidation | 3.61 | 1.16 | 3.1× |
| fermentation | 16.8 | 6.7 | 2.5× |
| methylotrophy | 0.88 | 0.53 | 1.7× |
| aerobic chemoheterotrophy | 28.9 | 42.5 | **0.7×** |
| methanotrophy | 0.23 | 0.52 | 0.4× |

Anaerobic/production functions (methanogenesis, H₂ oxidation, fermentation) are **enriched in
heartwood**; aerobic/consumption functions (aerobic chemoheterotrophy, methanotrophy) are **enriched
in sapwood** — the expected core→outer gradient, mirroring the ddPCR (methanogens in heartwood,
methanotrophs in sapwood).

## On "aerobic in hypoxic heartwood" — heartwood is NOT uniformly hypoxic
R2's implicit premise (hypoxic heartwood ⇒ fermentation should dominate) is too strong for this
system. **Be precise about what is ours vs cited:** the paired, depth-resolved heartwood-vs-sapwood
O₂ profiles — showing heartwood O₂ **not systematically lower** but far **more variable** (near-anoxia
to near-atmospheric) — are from the **companion Nature paper (Arnold & Gewirtzman et al. 2025)**, cited
as prior work, **not measured here**. This study's own O₂ data are **single internal (bark-to-pith)
gas measurements per tree** (n=157; not depth-resolved HW/SW), which are themselves highly variable
(~1–22% O₂) and in which internal CH₄ correlates **negatively with internal O₂** (R²=0.07; Fig. S8) —
consistent with, but not a direct measure of, a heartwood redox gradient. Given this patchiness and
small samples, a substantial predicted aerobic fraction in heartwood is **biologically plausible**,
not prima-facie "backwards"; and aerobic chemoheterotrophy is nonetheless **relatively lower** in
heartwood than sapwood (0.7×). (Over-assignment of aerobic chemoheterotrophy remains possible; both
readings are compatible and we do not over-interpret the absolute value.)

## H₂ economy — R2 is right on non-additivity, and it CORROBORATES hydrogenotrophy
Hydrogenotrophic methanogens oxidise H₂ but are filed under *methanogenesis* (3.2% HW), not *dark
hydrogen oxidation* (3.6% HW); the categories are non-additive, so total H₂-based metabolism is
≥ ~7%, not 3.6%. Rather than a bookkeeping problem, this **reinforces** an H₂-based economy in
heartwood — the same conclusion independently supported by our depleted δ¹³CH₄ (−63.7‰) and the
dominance of hydrogenotrophic Methanobacteriaceae.

## Where the inference is corroborated (so it is more than hypothesis-generating)
| Inference | Independent support (★ = this study's own data) |
|---|---|
| Hydrogenotrophic / H₂-based methanogenesis in heartwood | ★δ¹³CH₄ = −63.7‰; ★Methanobacteriaceae dominance; ★ddPCR mcrA |
| Anaerobic (heartwood) → aerobic (sapwood) redox gradient | ★ddPCR methanogen(Inner)/methanotroph(Outer) split (Figs 4,5); ★internal CH₄–O₂ negative correlation (Fig S8); + heartwood O₂ variability from companion Nature paper (cited, not measured here) |
| Methanotrophy concentrated in sapwood | ★ddPCR pmoA/mmoX distribution |
The in-study corroboration is the **ddPCR compartment split and our internal-gas correlations**; the
depth-resolved O₂ profiles are cited prior work. Where a predicted function has no such corroboration,
we present it as hypothesis-generating only.

## Caveat paragraph (drop-in — measured)
> Taxonomy-based functional predictions (FAPROTAX; PICRUSt2) infer metabolic potential from 16S
> affiliation and have well-known limitations, which are relevant here: they annotate only taxa with
> cultured references at genus level, and in our wood communities ~32% of ASVs resolved to a named
> genus and FAPROTAX assigned any function to under half of community abundance (46.9% heartwood).
> Predicted percentages are thus fractions of the whole community and conservative floors, and the
> categories are non-exhaustive and non-additive; we therefore do not interpret absolute values
> quantitatively. We instead emphasise between-tissue contrasts, which are internally consistent:
> anaerobic and H₂-based functions (methanogenesis, hydrogen oxidation, fermentation) are enriched in
> heartwood and aerobic functions (aerobic chemoheterotrophy, methanotrophy) in sapwood, mirroring
> our ddPCR distribution of methanogens (heartwood) and methanotrophs (sapwood) and the negative
> internal CH₄–O₂ correlation in our gas measurements; the heartwood O₂ environment is itself locally
> anoxic and highly variable (Arnold & Gewirtzman et al. 2025). Consistent with this, hydrogenotrophic
> methanogens (classified under methanogenesis) and dark hydrogen oxidation are both enriched in
> heartwood, indicating an H₂-based economy also supported by our δ¹³CH₄ values and the dominance of
> Methanobacteriaceae. We treat functional predictions as corroborative where they agree with these
> independent lines of evidence and as hypothesis-generating otherwise; quantitative conclusions rest
> on ddPCR of functional genes and on stable isotopes.

## Response to Referee 2 (#3)
> We agree these predictions should be interpreted cautiously and have added an explicit treatment of
> their limitations, while noting that taxonomy-based functional inference is a standard tool for 16S
> data. Quantifying the limitation for our data: ~32% of wood ASVs resolved to genus and FAPROTAX
> annotated under half of community abundance, so the predicted percentages are fractions of the whole
> community and conservative floors — we now avoid interpreting absolute values (including the
> fermentation fraction, which is not evidence that few taxa ferment) and instead emphasise the more
> robust between-tissue contrasts. The reviewer is also correct that "dark hydrogen oxidation" (3.6%)
> understates H₂ metabolism, because hydrogenotrophic methanogens (≈3.2%, likewise H₂ oxidisers) are
> classified under methanogenesis; taken together these non-additive categories imply ≥~7% H₂-based
> metabolism, reinforcing rather than complicating our conclusion of an H₂-based, hydrogenotrophic
> economy in heartwood — a conclusion independently supported by our δ¹³CH₄ values and by
> Methanobacteriaceae dominance. On the balance of aerobic versus fermentative predictions, we note
> that heartwood in this system is not uniformly hypoxic but highly variable in O₂ (shown with
> depth-resolved profiles in our companion study, Arnold & Gewirtzman et al. 2025; our single internal
> gas measurements here span ~1–22% O₂), so mixed aerobic and anaerobic signals are expected;
> nonetheless aerobic chemoheterotrophy is relatively depleted, and anaerobic functions relatively
> enriched, in heartwood versus sapwood, consistent with our ddPCR compartment distribution and the
> negative internal CH₄–O₂ correlation we observe. We have revised the text to present
> these predictions as corroborative where they align with independent evidence and as
> hypothesis-generating otherwise, and we state that our quantitative conclusions derive from ddPCR
> and isotopes rather than from 16S-based functional inference.
