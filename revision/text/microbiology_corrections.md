# Microbiology corrections (Referee 2) — current text → proposed fix

Each item quotes the current manuscript text (REVISED_manuscript_v2.md) and gives
a corrected version. All are cheap, correct, and make the paper more defensible.
These are DRAFTS for the revision — do not edit the manuscript yet (new-files-only
rule).

---

## R2 #4 (L412-413) — *Methylocapsa* / soluble MMO / multi-carbon growth

**Referee:** "Many Beijerinckiaceae genera especially *Methylocapsa* do not
encode soluble MMO and cannot grow on multi-carbon substrates."

**Current (manuscript L177):**
> "…while wood was dominated by Beijerinckiaceae genera (*Methylocella*,
> *Methylocapsa*) that use soluble MMO and can grow on multi-carbon substrates
> (Fig. S2)."

**Proposed:**
> "…while wood was dominated by Beijerinckiaceae, including *Methylocella* —
> facultative methanotrophs that possess soluble MMO and can grow on
> multi-carbon substrates — and *Methylocapsa*, most members of which lack
> soluble MMO and are not known to grow on multi-carbon substrates (Fig. S2)."

(R2 is right: the sMMO / multi-carbon attribute applies to *Methylocella*, not
to *Methylocapsa* generally. Don't attribute both traits to both genera.)

---

## R2 #2 — Methylacidiphilaceae "verified aerobic methanotrophs"

**Referee:** "For Methylacidiphilaceae, only geothermal Methylacidiphilaceae has
been known to be methanotrophic. All reported mesophilic Methylacidiphilaceae
genomes to date do not encode genes for methanotrophy" (cites peatland +
tree-stem genomes).

**Current (manuscript L177):**
> "…and Methylacidiphilaceae (mean 0.09% in heartwood), a family of verified
> aerobic methanotrophs (Fig. 5c-d; Fig. S2)."

**Proposed:**
> "…and Methylacidiphilaceae (mean 0.09% in heartwood). Although
> Methylacidiphilaceae are often cited as methanotrophs, verified methanotrophy
> in this family is currently restricted to thermoacidophilic geothermal
> members; mesophilic Methylacidiphilaceae genomes reported to date, including
> those from peatland and tree stems, lack the genes for methanotrophy (Sharp et
> al.; Arnold et al. 2025). We therefore treat these ASVs as *putative* rather
> than confirmed methanotrophs (see Table S[new])."

**Refs to add:** the two R2 pointed to — mesophilic peatland
Methylacidiphilaceae (microorganisms 9:2566, doi:10.3390/microorganisms9122566)
and the tree-stem *Science* paper (doi:10.1126/science.adu2182; = Arnold et al.
2025, already in refs). Verify exact citation.

---

## R2 #4 (L623) — reductive acetyl-CoA & RuMP are NOT fermentation

**Referee:** "Both reductive acetyl-CoA pathway and RuMP cycle are not directly
associated with fermentation. The former is an anaerobic carbon-fixation pathway
and the latter is associated with one-carbon compound assimilation."

**Two places in the manuscript conflate them with fermentation:**

**Current (L253, Discussion):**
> "Pathways associated with fermentation (reductive acetyl-CoA pathway, RuMP
> cycle), and anaerobic metabolism were positively correlated with mcrA…"

**Proposed (L253):**
> "Anaerobic and one-carbon pathways — including the reductive acetyl-CoA
> (Wood–Ljungdahl) carbon-fixation pathway and the RuMP cycle for one-carbon
> assimilation — as well as predicted fermentation pathways, were positively
> correlated with mcrA abundance in heartwood, while aerobic processes (TCA
> cycle, heme biosynthesis) were negatively associated (Fig. 6)."

**Current (L187, Results):** lists "the reductive acetyl coenzyme A pathway, and
formaldehyde assimilation (RuMP cycle), co-occurring with predicted fermentation
pathways (homolactic fermentation, mannan degradation)…" — this one is actually
**already correct** (it separates acetyl-CoA/RuMP from the fermentation
pathways). Just ensure the Discussion (L253) matches this careful Results
wording. Main fix is L253.

---

## R2 #4 (L695-700) — sMMO is low-affinity

**Referee:** "pMMO and sMMO are distinguished by their affinity towards CH4.
sMMO is [a] low-affinity enzyme… threshold usually > 500 ppm. Methanotrophs
that use sMMO likely only able to oxidize CH4 present at a significantly enriched
level."

**Current (manuscript L265):** discusses the pMMO:(pMMO+sMMO) shift and
copper/sMMO induction but does not mention the **affinity** distinction.

**Proposed — add a sentence at L265:**
> "Notably, sMMO is a low-affinity enzyme (CH₄ oxidation threshold typically
> >500 ppm), so sMMO-dependent oxidation is likely confined to the elevated
> internal CH₄ concentrations we observe within stems (often >1,000 ppm) rather
> than to atmospheric-concentration uptake at the bark surface (Semrau et al.
> 2010). This is consistent with our observation of net emission rather than
> atmospheric uptake, and implies that any high-affinity, atmospheric-CH₄
> oxidation at the stem surface would depend on pMMO-bearing methanotrophs."

This turns R2's correction into support for your net-source interpretation.

---

## R2 #1 — pmoA and mmoX reported separately

**Referee:** "pmoA and mmoX abundance should be reported separately as many
methanotrophs encode both."

**Action (not a sentence rewrite — an analysis/reporting change):**
- In the main abundance results (Section 4), report pmoA and mmoX **separately**
  in addition to the combined "pmoA + mmoX" methanotroph total.
- You already handle mmoX carefully in SI (Methods S1 — sapwood mmoX is your
  strongest single-tree predictor) and report the pMMO:(pMMO+sMMO) ratio (L265).
  Pull the separate values up into the main text.
- **NEW SCRIPT** needed if this requires re-tabulating from ddPCR — flag: I can
  write `revision/analysis/pmoa_mmox_separate.R` reading
  `data/compiled/ddpcr_gene_abundances.csv`. Confirm you want this.

---

## R2 #1 — units (dry vs wet weight)

**Referee:** "per gram dry weight or wet weight… dry weight for wood but wet
weight for soil. Please use the same unit."

**Action:** state the basis explicitly in Methods and standardize (or, if you
must keep both, add a clearly labeled conversion). Need to confirm which basis
each used. → prose fix in Methods + one clarifying clause in each figure/table
caption reporting copies g⁻¹.

---

## R2 #1 — italicize gene names globally

`mcrA`, `pmoA`, `mmoX` must be italic throughout (main text, captions, SI).
Referee 3 explicitly praised your taxonomic italics, so this is the one
inconsistency to sweep. → global find/replace at manuscript-finalization time
(not now). I can produce a checklist of every occurrence if useful.

---

## R2 #7 (L607-610) — δ¹³C-CH₄ caution + source CO₂

**Referee:** "interpretation of 13C-CH4 isotopic data requires more caution.
Methanogenesis can produce a wide range of δ13C-CH4 and the 13C depletion is
controlled also by the δ13C of the carbon source (e.g. CO2) which was not
measured… cite references."

**Current (manuscript L251):** already hedges ("Isotopic interpretation is
complicated by the co-occurrence of multiple methanogenic pathways and by
fractionation during oxidation… but the highly depleted values are difficult to
explain without substantial hydrogenotrophic production").

**Proposed — strengthen L251 by adding the source-CO₂ caveat + refs:**
> "…Isotopic interpretation is further complicated because the δ¹³C of CH₄
> depends not only on the methanogenic pathway but also on the δ¹³C of the
> substrate CO₂, which we did not measure; a given δ¹³CH₄ value is therefore
> consistent with a range of source conditions (Whiticar 1999; Conrad 2005).
> With that caveat, the strongly depleted values we observe (median −63.7‰) are
> difficult to explain without substantial hydrogenotrophic production, and are
> consistent with reported ranges for CO₂-reduction methanogenesis (refs)."

**Refs to add/confirm:** Whiticar 1999 (Chem. Geol.); Conrad 2005 (Org.
Geochem.) — standard δ¹³C-CH₄ pathway references.

---

## R2 #6 — figure fixes (quick)

- **Fig 2A:** make x-axis scale consistent across species panels (fix the
  Quercus rubra vs Carya ovata inconsistency); define the dotted line in caption;
  address outlier distortion (Acer rubrum, Quercus rubra) — consider consistent
  axis limits or noting clipped outliers. → figure-code task (new script).
- **Fig 4a, 4c:** add color legends to the stacked bar charts. → figure-code.
- **Fig 7b:** replace "copies µL⁻¹" with **copies g⁻¹** (R2: per-µL is
  sensitive to elution volume/extraction mass). This also aligns with the L149
  dry-weight standardization. → figure-code. NOTE this also affects Fig 7 as
  discussed; confirm copies g⁻¹ is available per sample.
- **Fig 7c:** add CH₄ flux **unit** (nmol m⁻² s⁻¹). → figure-code.

All figure fixes require NEW scripts in `revision/analysis/` that regenerate the
specific panels; none should edit existing figure code. Flag which you want
first.
