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
2. FDR < 0.01 (linear mixed-effects association with *mcrA*).

Result: 143 significant low-contribution pathways at FDR < 0.05 (109 at FDR < 0.01), of which
98 (~69%) are generic housekeeping and are **not shown**. The remaining pathways fall into the
functional categories below.

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

## Displayed pathways (top 3 by |t| per category; colours as in figure)

- **C1 / carbon fixation** (red, positive): formaldehyde assimilation II (RuMP cycle, t=+11.5);
  reductive acetyl-CoA pathway (+7.6); Calvin-Benson-Bassham cycle (+6.0).
- **Fermentation / carbohydrate** (gold, positive): L-glutamate degradation VIII → propanoate
  (+6.3); mannan degradation (+4.4); homolactic fermentation (+4.3).
- **Aerobic respiration** (dark blue, negative): TCA cycle IV (−6.3); TCA cycle V
  (2-oxoglutarate:ferredoxin oxidoreductase) (−4.2); TCA cycle I, prokaryotic (−3.9).
- **Oxidative biosynthesis — lipid, assimilatory S/N** (light blue, negative): superpathway of
  L-methionine biosynthesis by sulfhydrylation (−4.9); superpathway of sulfate assimilation +
  cysteine biosynthesis (−4.3); stearate biosynthesis II (−3.2).

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
