# Panel (c) — PICRUSt2 pathway associations, a-priori functional categories

## Why this replaces "top-8 by t"

The original Fig. 6 selection ("top positive/negative pathways by t") surfaced mostly
generic housekeeping biosynthesis (nucleotide, cofactor, cell-wall pathways) that co-vary
with *mcrA* because they track total anaerobic biomass, not because they are mechanistically
tied to methanogenesis. Panel (c) instead classifies the **full significant set** by MetaCyc
pathway function — so category membership is fixed by biology, independent of the association
statistic — then shows representative pathways per mechanistically-relevant category.

## Filters (identical stringency logic to Fig. 6)

1. Predicted contributions from methanogen-classified ASVs removed (`no_mcra_otus` table),
   then pathways with mean methanogen-ASV contribution ≥ 10% excluded (`gc < 0.10`). This is
   the same non-circularity filter as Fig. 6 — it removes the archaeal/methanogenesis-cofactor
   pathways (coenzyme M/B, tetrahydromethanopterin, methanogenesis-from-acetate, archaeal
   lipids), which are the methanogens' own genomes and would be tautological.
2. FDR < 0.05 (linear mixed-effects association with *mcrA*).

Result: 143 significant low-contribution pathways, of which ~69% are generic housekeeping (and
lipid/fatty-acid biosynthesis) that are **not shown**. The remaining pathways fall into the four
displayed functional categories below; up to four representative pathways per category are shown
(top by |t|).

**Display note (lipid + threshold):** the *displayed* Sulfur/nitrogen category excludes lipid/
fatty-acid biosynthesis so that the assimilatory nitrate- and sulfate-reduction pathways (ranks
#3–4 within an S/N-only category) are visible; in a merged lipid+S/N category they would rank #8–9
behind five lipid pathways. Lipid biosynthesis is retained in the SI classification table under its
own label. The FDR threshold is 0.05 (looser than Fig. 6's 0.01) so that the two reduction pathways
(FDR 0.015 and 0.020) are included; this is noted in the caption.

## Full-set category summary (FDR < 0.05, all pathways; independent of display)

| Category | n | sig POS | sig NEG | direction |
|---|---|---|---|---|
| Archaeal / methanogenesis (filtered out by gc<0.10) | 17 | 15 | 0 | + (circular — excluded) |
| C1 / carbon fixation | 5 | 4 | 0 | + |
| Fermentation / carbohydrate degradation | 42 | 10 | 4 | + (net) |
| Aerobic respiration / oxidative | 33 | 1 | 17 | − |
| Lipid / fatty-acid biosynthesis | 13 | 0 | 11 | − |
| Redox: sulfur (assimilatory) | 5 | 0 | 3 | − |
| Redox: nitrogen (assimilatory) | 3 | 0 | 1 | − |
| Other housekeeping (excluded) | 289 | 59 | 37 | mixed (noise) |

**Interpretation:** the association is a clean anaerobic-C-flow vs. aerobic-oxidative axis.
There is essentially **no dissimilatory N/S/Fe/Mn respiration signal** — the only redox pathways
present are *assimilatory* sulfate/nitrate and they trend negative (aerobe-associated). So we do
not frame this as metal/N/S redox cycling.

## Displayed pathways (top 4 by |t| per category, FDR < 0.05; colours as in figure)

- **C1 / carbon fixation** (red, positive): formaldehyde assimilation II (RuMP cycle, t=+11.5);
  reductive acetyl-CoA pathway (+7.6); Calvin-Benson-Bassham cycle (+6.0). [only 3 in category]
- **Fermentation / carbohydrate** (gold, positive): L-glutamate degradation VIII → propanoate
  (+6.3); mannan degradation (+4.4); homolactic fermentation (+4.3); glycolysis III (+3.4).
- **Aerobic respiration** (dark blue, negative): TCA cycle IV (−6.3); TCA cycle V
  (2-oxoglutarate:ferredoxin oxidoreductase) (−4.2); TCA cycle I, prokaryotic (−3.9);
  aerobic respiration I, cytochrome c (−3.9).
- **Sulfur / nitrogen metabolism** (light blue, negative): L-methionine biosynthesis by
  sulfhydrylation (−4.95); sulfate assimilation + cysteine biosynthesis (−4.31);
  **nitrate reduction VI, assimilatory (−2.84)**; **sulfate reduction I, assimilatory (−2.71)**.

## Redox-ladder interpretation (and its limit)

The joint depletion of aerobic respiration *and* assimilatory nitrate/sulfate reduction where
*mcrA* is high is consistent with methanogens occupying the most reducing, oxidant-poor microsites —
the classic thermodynamic redox ladder, in which methanogenesis is favoured only once higher-yield
electron acceptors (O₂ > NO₃⁻ > … > SO₄²⁻) are drawn down. **Important limit:** the significant N/S
pathways here are *assimilatory* (incorporating N/S into biomass), not *dissimilatory* (using
NO₃⁻/SO₄²⁻ as terminal electron acceptors in respiration). The only dissimilatory pathway present —
nitrate reduction I (denitrification) — is non-significant (t = −0.42, FDR = 0.74). So this is best
read as a **community redox-state proxy**, not direct evidence of terminal-electron-acceptor
competition. Phrase text accordingly ("consistent with," not "demonstrates").

## Defensibility notes
- Category membership is assigned from MetaCyc pathway descriptions (keyword rules, documented in
  `R_fig_hydrogenotrophy.R::classify_pw`), **before** examining the association — no outcome-based
  selection within the mechanistic categories.
- "H₂-producing" and "methanogen-syntrophic" are *not* MetaCyc attributes; they are ecological
  interpretations (H₂ production is inferred from fermentation). They are therefore discussed in
  text, not used as figure categories.
- MetaCyc's Expected-Taxonomic-Range / ontology (which would give rigorous "archaeal" and
  "obligate anaerobe/aerobe" tags) is license-restricted and not in the repo; the description-based
  keyword scheme is the reproducible substitute and is auditable via the full classification table.
- **SI:** export the full 407-pathway classification + association table so every assignment is
  inspectable.
