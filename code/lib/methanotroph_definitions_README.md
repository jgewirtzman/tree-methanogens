# Methanotroph definitions

`methanotroph_definitions.csv` is the one list every 16S analysis uses to call a taxon a
methanotroph; `load_methanotroph_definitions.R` reads it (`load_methanotroph_defs()`,
`classify_methanotrophs()`, `load_placed_asvs()`).

A taxon counts as a methanotroph if it carries methane monooxygenase (MMO; particulate or
soluble). The list follows Knief (2015), cross-referenced to SILVA 138; for genera in mixed
families the evidence is the NCBI genome screen of Table S3
(`mmo_capacity_screen_table_S3_2026-09-30.csv`, see `mmo_capacity_screen_README.md`).

Classification, first match wins:

1. **Known**: a genus that carries MMO (the cultivated methanotroph genera of Knief 2015,
   plus *Candidatus* Methylomirabilis, *Crenothrix* and *Candidatus* Clonothrix), or an ASV
   in a family whose characterized members all carry MMO (Methylococcaceae,
   Methylomonadaceae, Methylothermaceae, Methylomirabilaceae).
2. **Putative** (an upper bound): an ASV in a mixed family (Beijerinckiaceae,
   Methylocystaceae, Crenotrichaceae, Methylacidiphilaceae) whose genus is unresolved or an
   uncultured placeholder, or a genus whose MMO capacity is untested or found in only some
   genomes (*Methylovirgula*, *Rhodoblastus*, *Methylorosula*, *Psychroglaciecola*).
3. **Not counted**: resolved genera in those families whose genomes lack MMO (for example
   *Methylobacterium-Methylorubrum*, *Bosea*, *Roseiarcus*).
