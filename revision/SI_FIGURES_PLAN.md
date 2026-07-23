# Final SI figure list (revision) — decided

## Decisions (Jon)
- S6 taxonomy mcrA heatmap → **keep in SI** (Fig 6a distills it; full heatmap stays SI).
- Plant traits → **one composite** SI figure (heatmap + ordination + breakdown).
- Extended isotopes → **yes**, one composite; **fold in whole-tree sources**.
- Seasonal SI → **skip** (Fig 9 seasonal panel + RF results cover it).
- S14 mcrA vs methanotroph → **keep** (important).
- **Add:** regfig_meth + regfig_mcra (per-regression), height_slope_vs_moisture,
  iso_wholetree_sources (wrapped into the isotope composite).
- Methods S1–S3 → main text (separate from the figure series; doesn't remove any Fig S).

## Proposed organized SI list (numbering provisional; follows main-text flow)

### Flux & environment
| SI | Content | Source |
|---|---|---|
| 1 | Moisture overlay | old S1 |
| 2 | Height-decline slope vs soil moisture | NEW `height_slope_vs_moisture` (R3) |
| 3 | Random-forest flux predictions (upscaling support) | old S15 |

### Gene–flux relationships
| 4 | Per-regression: mcrA – flux (individual / aggregate-x / species) | NEW `regfig_mcra` |
| 5 | Per-regression: methanogenesis – flux | NEW `regfig_meth` |
| 6 | Scale-dependent gene patterns | old S11 |
| 7 | mcrA vs methanotroph | old S14 |

### Microbial abundance & composition
| 8 | Methanotroph abundance patterns | old S10 |
| 9 | pmoA vs mmoX by compartment | NEW `SI_fig_pmoa_mmox_separate` (R2 #1c) |
| 10 | Taxonomy mcrA heatmap (full) | old S6 |
| 11 | Taxonomy pmoA heatmap | old S2 |

### Functional / pathway inference
| 12 | FAPROTAX heatmaps | old S3 |
| 13 | PICRUSt mcrA-all heatmap | old S4 |
| 14 | PICRUSt mcrA (no-mcrA) heatmap | NEW — demoted main Fig 6 |
| 15 | PICRUSt pmoA heatmap | old S5 |

### Isotopes
| 16 | δ13CH4 raincloud (full) | old S9 |
| 17 | Extended isotopes composite: Keeling + CH4-vs-CO2 cross-plot + εC + whole-tree sources | NEW (`iso_co2_keeling`, `iso_co2_ch4_vs_co2`, `iso_epsC_weighted`/`iso_co2_epsC_distribution`, `iso_wholetree_sources`) |

### Internal gas & anatomy
| 18 | Internal-gas beeswarm | old S7 |
| 19 | Internal-gas profiles | old S8 |
| 20 | Black-oak methanome | old S12 |
| 21 | Tree radial sections | old S13 |

### Plant traits
| 22 | Plant-trait composite: heatmap + ordination (PCA/Procrustes) + density breakdown | NEW (`traits_heatmap`, `traits_pca_biplot`/`traits_procrustes`, `traits_breakdown`) |

**Total: 22 SI figures.**

## To build (new composites)
- **SI-17 extended isotopes** — assemble Keeling + CH4-vs-CO2 + εC + whole-tree sources into one 4-panel figure (new script).
- **SI-22 plant-trait composite** — assemble heatmap + ordination + breakdown into one figure (new script).
- SI-14 (PICRUSt no-mcrA) is the existing main Fig 6 image, just relabeled/renumbered.

## SI TABLES (separate track)
Known/putative taxa; pmoA/mmoX by compartment + by species; pmoA/mmoX robustness;
panel-c pathway classification; primer sequences.

## Open / to confirm
- 22 SI figs is on the larger side — could consolidate regfig_mcra + regfig_meth into one,
  or drop regfig_ratio (not included). Flag if leaner is wanted.
- Confirm which εC panel to use (`iso_epsC_weighted` vs `iso_co2_epsC_distribution`).
