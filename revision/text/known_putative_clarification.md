# Known/Putative methanotroph & methanogen classification — clarified (R2 #2)

Table: `revision/outputs/known_putative_taxa_table.csv` · counts: `kp_counts.txt`
Script: `R_known_putative_table.R` · Original scheme: `data/compiled/methanotroph_definitions.csv`
Revised scheme (one change): `revision/data/methanotroph_definitions_revised.csv`

## Approach (decided): keep the existing scheme, document it, one within-scheme fix
The existing scheme (curated from **Knief 2015** + SILVA 138) is already more careful
than typical 16S studies: per-taxon Known/Putative flags with written rationale, an
explicit rule, and a genus whitelist. We keep it, document it, and make **one
correction that follows the scheme's own rule** — no overhaul, no new tiers, no
"known non-methanotroph" exclusion list (that would require asserting function from
16S, which is what we're avoiding; "putative" already means "not confirmed").

## The scheme's own rule (made explicit)
- **Known** = ASV genus is a cultivated methanotroph genus (whitelist), OR ASV family
  is an *exclusive* methanotroph family (all members are methanotrophs), genus unresolved.
- **Putative** = ASV family *contains methanotrophs and non-methanotrophs* (mixed
  family; e.g. Beijerinckiaceae, Methylocystaceae), genus unresolved → cannot confirm.
- **Methanogen** = one of 9 established methanogen families; Bathyarchaeia (debated) excluded.
- 16S cannot confirm methanotrophy (a functional trait); assignments reflect
  phylogenetic affiliation, and "putative" carries that uncertainty by design.

## The one fix — and it FOLLOWS your rule, it doesn't override it
Your rule: **exclusive families → Known; families that contain non-methanotrophs →
Putative only.** Methylacidiphilaceae was placed in the *exclusive* bucket on Knief
2015. But post-Knief genomes (R2's refs: peatland, tree-stem) show mesophilic
Methylacidiphilaceae **lack methanotrophy genes** — i.e. it is a **mixed family**.
So by your own criterion it moves to **Putative** (Include_known YES→NO), exactly like
Beijerinckiaceae. The genus *Methylacidiphilum* (cultivated geothermal methanotroph)
stays Known — parallel to Beijerinckiaceae-family-putative / Methylocapsa-genus-known.
=> same rule, one updated fact. (Also soften manuscript L177 "verified aerobic
methanotrophs".) Done in `methanotroph_definitions_revised.csv`.

## Bounded due-diligence (so we needn't re-audit the whole literature)
A stale fact can only wrongly inflate **Known** in two small, finite lists; we checked
only those (not all 79,769 ASVs):
- **Exclusive families detected:** Methylococcaceae (keep — established Type I family),
  Methylacidiphilaceae (→ Putative, the fix), **Methylomonadaceae (keep; FLAG for Jon's
  expert eye — newer family; named genera are methanotrophs but confirm SILVA 138 does
  not lump uncultured non-methanotroph lineages)**.
- **Known genera detected:** Methylocapsa, Methylocella, Methyloferula, Methylobacter,
  Methyloglobulus, Methylomagnum — all bona fide cultivated methanotroph genera; hold.
Everything in "Putative" is conservative by design, so stale facts there don't
over-claim. Net: **one change (Methylacidiphilaceae) + one flag (Methylomonadaceae).**

## Count effect (nothing rests on it, but for completeness)
Reclassifying family-level Methylacidiphilaceae Known→Putative moves its ASVs out of
Known (in the full unfiltered set, Known 465→77; re-run on the manuscript's filtered
ASV set for the exact reported 434/627 figures). Genuinely-Known methanotrophs are
then dominated by **Methylocapsa** (your wood sMMO genus); methanogens by
**Methanobacterium** — consistent with the paper. The quantitative methanotroph result
uses ddPCR pmoA/mmoX, not this 16S split.

## Methods paragraph (drop-in)
> Methanotroph and methanogen taxa were classified from 16S ASV taxonomy (SILVA 138)
> following a curated scheme based on Knief (2015) (Supplementary Table SX). Because
> 16S cannot confirm a functional trait, we distinguish "known" methanotrophs — ASVs
> assigned to a cultivated methanotroph genus, or to a family in which all members are
> methanotrophs — from "putative" methanotrophs, i.e. family-level assignments to
> methanotroph-containing families (e.g. Beijerinckiaceae) that also include
> non-methanotrophs and lack genus resolution. We reverified the exclusive-family
> designations against current genomic evidence; on this basis Methylacidiphilaceae,
> whose mesophilic members lack methanotrophy genes, is treated as putative at the
> family level (the cultivated geothermal genus Methylacidiphilum remains known).
> Methanogens comprise nine established methanogen families; Bathyarchaeia, with
> debated methanogenic capacity, are excluded.

## Response to Referee 2 (#2) — final
> We have added a supplementary table listing every detected taxon classified as a
> known or putative methanotroph or methanogen, with rank, classification, ASV count,
> relative abundance, and the basis for inclusion (Knief 2015; SILVA 138). We clarify
> our criteria: "known" methanotrophs are ASVs assigned to a cultivated methanotroph
> genus or to a family in which all members are methanotrophs; "putative" methanotrophs
> are family-level assignments to methanotroph-containing families (e.g. Beijerinckiaceae)
> that also include non-methanotrophs and lack genus resolution — so putative status
> already reflects that these are affiliations, not confirmed methanotrophs. Applying
> this same criterion with current evidence that mesophilic Methylacidiphilaceae lack
> methanotrophy genes, we now treat family-level Methylacidiphilaceae as putative rather
> than known, and have corrected the statement describing it as verified methanotrophs.
> We emphasise that our quantitative methanotroph results derive from ddPCR of the pmoA
> and mmoX functional genes, not from 16S-based classification.
