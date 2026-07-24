# ddPCR negative controls, LoD, recovery, primers (answers R2 #1)

> # 🚨 CRITICAL — MUST REVISIT BEFORE RESUBMISSION 🚨
> **All ddPCR gene abundances (mcrA, pmoA, mmoX — every methanogen and
> methanotroph load) may be ~10× UNDERESTIMATED** by the current pipeline.
> `02_harmonize_all_data.R:317` uses `Conc × 75 / mass`, which multiplies by the
> 75 µL elution volume but **omits the 10× dilution of sample into the ddPCR
> reaction** (2.5 µL sample in a 25 µL reaction, per methods-paper SI Table 1).
> The physically-complete conversion is `Conc × 750 / mass` = **10× higher**.
> - **ACTION: Jon to ask Wyatt (Arnold) to confirm the intended conversion / QuantaSoft
>   setup.** If confirmed, **ALL ddPCR copies/g in the manuscript must be corrected ×10.**
> - It is a UNIFORM 10×, so the **methanogen:methanotroph ratio and every relative
>   result are UNCHANGED**; only absolute loads shift.
> - **16S is separate (qPCR, UMGC) — NOT affected, leave as-is.**
> - The **LoD (500 cells/100 mg) is unaffected** (derived correctly in the methods paper).


Source: **Arnold, Gewirtzman, et al. 2024, "A method for sampling the living wood
microbiome," *Methods in Ecology and Evolution* 15:1084–1096, DOI
10.1111/2041-210X.14311** (the methods paper this study cites for ddPCR).
Data package: Dryad doi:10.5061/dryad.vq83bk3zj.

R2's asks: "What is the confidence of detection? Did negative controls produce
any quantifiable numbers? Negative control data should be reported explicitly."
Plus units (dry vs wet), and pmoA/mmoX primer sequences (L248-249).

## Citable facts (all from Arnold et al. 2024)

**Negative controls (NTC) — the direct answer:**
- The mcrA EvaGreen assay had a **persistent low rate of false positives in wood**
  (NTC ~10–20 copies; their Fig 2). This is *why* a degenerate **probe-based mcrA
  assay** was developed — the assay used for the wood/low-biomass detections here.
- The **probe-based mcrA assay produced ZERO false positives** across all 6 NTCs
  (their SI Fig 2 / Fig 2). mmoX and pmoA (EvaGreen) NTCs also read **~0** (Fig 2).
- Unspiked field controls confirmed **no pre-existing target sequences** capable
  of confounding results (their Fig S3).
- So: negative controls did **not** produce quantifiable signal for the assays
  used here. This is the clean answer to R2.

**Limit of detection (LoD):**
- Assay-level: the mcrA probe assay detects consistently to **3 copies/reaction
  (6/6)**, likely **~1 copy** (SI Fig 2); strong linearity (slope 0.99, R²=0.997)
  over 5 orders of magnitude.
- Method-level (incl. ~20% recovery through the full pipeline): **~300–600 cells
  per 100 mg dry wood** (conservative 500) — ~200× below typical soil methanogen
  abundance.

**Recovery / conservative bias:**
- Recovery efficiency averaged **~20%** (22.2% ± 3.6% across the three targets in
  field samples; 81.7% relative to an idealized liquid control). ⇒ reported
  abundances are **conservative underestimates** (no efficiency correction applied).
- PCR inhibition was minimal (dilution-series slopes ~1, R²>0.99; a Zymo inhibitor-
  removal kit gave no significant difference).

**Units (R2's dry-vs-wet point):** the wood method is explicitly **per 100 mg DRY
wood** (freeze-dried, dry-mass basis); soil extractions used **250 mg wet** soil.
This matches the copies/g dry-vs-wet issue we already flagged.

**Primers/probes (R2 L248-249):**
- mcrA probe assay (FAM, selected) — Arnold et al. 2024 SI Table 1:
  Forward `ACGACYTRCAGGAYCAGTGY`; Probe `WGGWCCWAACTAYCCBAACTACG`; (reverse in
  SI Table 1 — pull the exact sequence from the docx/table).
- pmoA and mmoX used EvaGreen with **literature primer sets**: Bourne et al. 2001;
  Fuse et al. 1998; Luton et al. 2002. R2 wants these *sequences in our Methods*,
  so pull them from those refs (e.g. pmoA A189/mb661-type; mmoX206/886-type) and
  tabulate — don't just cite.

## "Values below the LoD" — RESOLVED (see full explanation lower down)
Short version: the LoD is a *consistent-detection* threshold, not a floor. (The
methods paper set it empirically: 6/6 replicates detected at 3 copies/reaction,
2/6 at 0.33 copies/reaction -> "LoD below 3, near 1 copy" -> ~500 cells/100 mg. It
did NOT use a formal 95% criterion.) Reported values below it are real detections
of genuinely low-abundance samples
that fired a few droplets (NTC = true zero, so not contamination), caught in the
probabilistic zone between "always detectable" (~500 cells) and the ~1-copy floor.
Detection rate is therefore a conservative floor, and claims are order-of-magnitude
/ aggregated. No LoQ is introduced. (The earlier apparent "tension" was partly a
units artifact — comparing a cells LoD to copies/g — and partly the x10 conversion
issue; neither is a real problem. Detail below.)

## Is the LoD defensible & replicable? YES — reconstructed from the Dryad data

We replicated the two ingredients from the raw ddPCR (Dryad doi:10.5061/dryad.vq83bk3zj):

1. **Assay sensitivity + true-zero NTC** (`mcrA_probe_Stds_curve_*.csv`):
   - **6/6 no-template controls = 0 positive droplets** (of ~17,000–20,000 accepted
     droplets each) → genuine zeros, not ambiguous signal.
   - Detection floor: a standard well expected to hold ~1 copy shows exactly **1
     positive droplet**; consistent detection by ~3 copies/reaction (dMIQE-style
     LoD ≈ 1–3 copies/reaction).
   - **Clean amplitude separation**: positive droplets ~1,900 RFU vs negatives
     ~720 RFU — a wide gap, so the classifier threshold is unambiguous.

2. **Method LoD = assay LoD propagated back through the pipeline** (defensible arithmetic):
   ```
   method LoD (cells / 100 mg dry) = assay_LoD(copies/rxn) x (elution_vol / template_vol) / recovery
                                   ≈ 3 x (75 µL / 2 µL) / 0.20  ≈ 560 cells
   (range ~300–600 for assay_LoD 1–3 copies and recovery 20–25%)
   ```
   This is exactly the paper's "300–600 cells / 100 mg dry wood, based on 2 µL of
   eluent." So the LoD is a transparent, reproducible propagation of the empirical
   assay floor through the extraction recovery and template-subsampling — proper
   and defensible.

## How ddPCR gives TRUE zeros (for R2 — the mechanistic explanation)
- The reaction is partitioned into ~20,000 nL-scale **droplets**; each either
  contains ≥1 target (amplifies to high fluorescence) or not (stays at baseline).
- After endpoint PCR a reader measures each droplet's **fluorescence amplitude**;
  the data are **bimodal** (here ~1,900 vs ~720 RFU). A **threshold** is placed in
  the gap: droplets above = positive, below = negative. Concentration comes from
  the positive fraction by **Poisson** statistics — no standard curve.
- **Therefore a negative control with no droplets above threshold is a literal
  0/18,000 — a true zero**, not a late-cycle ambiguous read. This is the key
  contrast with qPCR, where a no-template control still eventually crosses the Cq
  threshold at a late cycle (primer-dimer / nonspecific amplification), yielding a
  number that cannot cleanly be called zero. R2's "did the negative controls
  produce quantifiable numbers?" has a categorical answer in ddPCR: no positive
  droplets = zero.
- The **probe** assay is why separation is clean: a hydrolysis probe fluoresces
  only on specific target amplification, so positives/negatives separate sharply.
  The intercalating-dye (EvaGreen) mcrA assay bound nonspecific products → poor
  separation → the false positives that motivated the probe assay.

## "500 cells -> 1–3 positives," and how it squares with 20% recovery
Yes — that is exactly the chain, and the 20% recovery is *baked into* it:
- 500 cells in 100 mg dry wood -> ~20% extraction recovery -> ~100 gene copies in
  the 75 µL eluent -> 2 µL of that eluent into the reaction = 100 x (2/75) ≈ **2.7
  copies in the reaction -> ~1–3 positive droplets** (above the zero NTC).
- So the 20% loss is *why* you need ~500 cells (not ~100) to reach the ~3-copy
  detection floor. Recovery and LoD are consistent, not contradictory: the LoD is
  the cell count that survives the loss and the subsampling to hit the assay floor.

## Recovery correction — is not-correcting standard? (yes)
- Reporting molecular abundances **without** an extraction-efficiency correction is
  standard for environmental qPCR/ddPCR: efficiency is sample-specific and usually
  unmeasured, so "copies per g" is reported **as recovered** — a conservative lower
  bound. Studies that do correct use internal spike-ins/whole-cell standards.
- **qPCR is handled the same way**: copies are read off a standard curve of known
  copies and reported as measured in the extract, typically not efficiency-
  corrected. (qPCR vs ddPCR differ in the *quantification* method — curve vs Poisson
  droplet counting — not in the recovery-correction convention.)
- This method went **beyond** the norm by measuring recovery (~20%) via dowel/cell
  spiking, so we can state the conservative bias explicitly: reported abundances are
  ~5x underestimates of the true in-wood amount. We keep values uncorrected (the
  defensible, conservative choice) and say so.

## Normalization — the copies/g conversion (the ~10x issue), walked through plainly

**Physical reality** (methods-paper SI Table 1 volumes + verified from the data):
- Extract a sample of mass M into **75 µL** elution buffer (all copies live here).
- Pipette **2.5 µL** of that into a **25 µL** ddPCR reaction (2.5 sample + 22.5 mix)
  → the sample is **diluted 10x** (2.5 into 25) going into the reaction.
- ddPCR reports `Conc` = copies per µL **of the reaction** (a concentration; it is
  independent of how much volume is dropleted — the 20-vs-25 µL "well" is only
  QuantaSoft bookkeeping and does NOT change the concentration).

**Correct copies/g** (worked example: Conc = 1 copy/µL, mass = 100 mg):
```
Conc (reaction)          = 1  copy/µL
x 10  (undo dilution)    = 10 copies/µL of ELUENT
x 75  (elution volume)   = 750 copies extracted total
/ 0.1 g (100 mg)         = 7,500 copies/g          => formula: Conc * 750 / mass
```
**Pipeline** (`02_harmonize_all_data.R:317`): `Conc * 75 / mass` -> 1 * 75 / 0.1 =
**750 copies/g** for the same sample. It applied the 75 µL elution but **skipped the
x10 sample->reaction dilution** -> **exactly 10x low** (7,500 vs 750).

The factor is a clean **10x** = total reaction volume / sample volume = 25 / 2.5.
(My earlier "8x" hedge was WRONG: I let QuantaSoft's 20 µL reporting convention leak
into the dilution factor; concentration is concentration, so droplet volume is
irrelevant. It is 10x.)

**What each place does:**
- *Methods paper LoD*: works forward cells->copies and correctly accounts for the
  2.5 µL-in-75 µL aliquot + ~20% recovery -> 500 cells/100 mg. CORRECT, unaffected.
- *Manuscript pipeline*: `Conc x 75` -> misses the x10 -> absolute copies/g ~10x low.
- *Correct*: `Conc x 750 / mass`.

**Status of our scripts:** left at x75 to MATCH the manuscript pending confirmation
(so nothing diverges silently). If Wyatt confirms the x10 is missing, switch the
manuscript AND our scripts to x750 (x10 correction), uniformly.

**Internal check (ddPCR-only, reassuring):** with x10 correction, heartwood mcrA
median ~5,000 copies/g and the LoD in the same recovered units ~1,000 copies/g ->
typical heartwood ~5x above detection (also retires the "near-LoD?" worry).

## Values BELOW the LoD — what they mean (resolved; no LoQ needed)
The LoD (~500 cells/100 mg) is a **consistent-detection threshold, NOT a floor** on
reported values. The methods paper set it empirically from the dilution series:
**6/6 replicates detected at 3 copies/reaction, 2/6 at 0.33 copies/reaction**, so
"LoD below 3 copies, near 1 copy" -> propagated to ~500 cells/100 mg. (No formal
95% criterion was used; that is just the generic dMIQE definition, not this paper's
method.) The absolute minimum detectable is ~1 copy (any time a real target
molecule fires a droplet). Between the ~1-copy floor and the consistent-detection
level is a **probabilistic zone — directly observed** as the 2/6 (33%) detection at
0.33 copies. So a reported value below the LoD = a **real, low-abundance
sample that fired a few droplets this time**: detection is genuine (NTC = true
zero), but at that level it is not guaranteed to re-detect and the exact number is
imprecise. This does NOT undermine anything: (1) our detection rate is a
conservative floor (some equally-low samples gave 0 droplets and were scored
non-detect); (2) claims are order-of-magnitude and species-aggregated, never
resting on an individual near-floor value. No LoQ is introduced (R2 did not ask
for one, and a CV threshold is arbitrary); detection is the robust claim, and
precision is adequate for log-scale comparisons (see below).

## Draft response to Referee 2 (#1, negative controls / LoD) — FINAL (LoD only)
> We now report the detection limit and controls explicitly. mcrA in wood was
> quantified with a degenerate probe-based ddPCR assay developed specifically to
> eliminate the low, persistent false-positive rate of intercalating-dye mcrA
> assays in wood (Arnold et al. 2024). Because ddPCR classifies each of ~18,000
> partitioned droplets as positive or negative, a no-template control that yields
> no positive droplets is a literal zero; across all no-template controls this
> probe assay produced no positive droplets, and mmoX and pmoA no-template controls
> were likewise negative, so contamination cannot account for our detections.
> Unspiked field controls contained no pre-existing target. The method limit of
> detection, incorporating extraction recovery (~20%), is ~500 cells per 100 mg
> dry wood — roughly two orders of magnitude below typical soil methanogen
> abundances — and ~97% of heartwood samples were positive. Reported abundances
> span ~5 orders of magnitude and are analysed on a log scale, where per-sample
> counting precision is small relative to that range; we therefore emphasise
> detection frequency and species-aggregated abundances, for which precision is
> more than sufficient. Because recovery is uncorrected, reported abundances are
> conservative. Full assay development, controls, and LoD analysis are in Arnold
> et al. (2024).

## Draft Methods clause
> mcrA was quantified using a degenerate probe-based ddPCR assay (FAM; primers/
> probe in Arnold et al. 2024, Table S1) designed to eliminate false positives in
> low-biomass wood; pmoA and mmoX were quantified by EvaGreen ddPCR using primers
> of [Luton et al. 2002 / Fuse et al. 1998 / Bourne et al. 2001] (sequences in
> Table S[x]). All assays included no-template and unspiked field controls
> (Arnold et al. 2024). Abundances are reported per gram dry wood; the method
> limit of detection is ~500 cells per 100 mg dry wood.
