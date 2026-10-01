# MMO capacity screen, 2026-09-30

`mmo_capacity_screen_2026-09-30.csv` is the frozen evidence behind the genus rows in
`methanotroph_definitions.csv` marked "Capacity rule 2026-09-30". NCBI changes over time,
so the snapshot is committed rather than regenerated.

Reproduce with `code/tools/mmo_capacity_screen.R` (live NCBI E-utilities;
needs network and the Bioconductor package `pwalign`).

- `assemblies` — NCBI Assembly records for the genus.
- `*_by_name` — NCBI Protein records whose title or gene name places them in the
  copper-MMO (pMMO/AMO) or soluble di-iron MMO family. Name-based, so over-inclusive:
  the copper query also catches sMMO component names, and RefSeq's "ammonia monooxygenase"
  is also used for an unrelated membrane protein.
- `verified_mmo_capacity` — after aligning each hit (BLOSUM62, local) to M. capsulatus
  PmoA (Q607G3) and M. trichosporium MmoX (P27353). True subunit: >=40% id over >=150 aa to
  PmoA, or >=55% id over >=300 aa to MmoX.
- Positive controls: Methylocapsa (pMMO), Methylocella (sMMO).
