# Final figure manifest (revision)

> Regenerate with `Rscript code/revision/run_all.R`. Authoritative dest→source mapping: `code/revision/rev_00_assemble_figures.R`.

Clean, correctly-numbered copies of the final figures. `main/` = main text only,
`SI/` = supplement only. Regenerate any file from the script in the Source column.
Numbers reflect the agreed restructure (hydrogenotrophy → main Fig 6; old PICRUSt-mcrA
heatmap → SI; SI resequenced per `SI_FIGURES_PLAN.md`).

## Main figures
| # | File | Content | Source script | Source PNG |
|---|---|---|---|---|
| 1 | Figure_1_temporal-flux | Temporal flux time series (soil + tree) | `rev_fig01_temporal-flux.R` | fig1_final.png |
| 2 | Figure_2_height-flux | Height-dependent flux | `rev_fig02_height-flux.R` | fig2_final.png |
| 3 | Figure_3_variance-partition | Variance partitioning | `rev_fig03_variance-partition.R` | fig3_final.png |
| 4 | Figure_4_gene-abundance | mcrA / methanotroph abundance by compartment (copies/g; ×10 toggle) | `rev_fig04_gene-abundance.R` | fig4_final.png |
| 5 | Figure_5_methane-cycling-composition | Combined methane-cycling composition (revised putative) | `rev_fig05_methane-cycling.R` | fig5_final.png |
| 6 | Figure_6_hydrogenotrophy | Convergent-hydrogenotrophy synthesis **(NEW — replaces old PICRUSt-mcrA heatmap)** | `rev_fig06_hydrogenotrophy.R` | fig_hydrogenotrophy.png |
| 7 | Figure_7_felled-oak | Felled-oak profiles (copies/g + flux unit) | `rev_fig07_felled-oak.R` | fig7_felled_oak_final.png |
| 8 | Figure_8_radial-species | Radial species comparison **(kept as-is)** | `05_gene_flux_analysis/04_species_gene_flux.R` | outputs/figures/main/fig8_radial_species_comparison.png |
| 9 | Figure_9_ch4-budget | CH₄ budget: maps + seasonal + net waterfall **(NEW — replaces old `fig9_upscaled_flux_seasonal`)** | `rev_fig09_budget.R` | fig_budget_maps.png |

## SI figures (18)
| # | File | Content | Source | Was |
|---|---|---|---|---|
| S01 | moisture-overlay | Moisture overlay | orig | old S1 |
| S02 | height-slope-moisture | Height-decline slope vs soil moisture (R3) | NEW `rev_figS02_height-slope-moisture.R` | — |
| S03 | mcra-vs-methanotroph | mcrA vs methanotroph | orig | old S14 |
| S04 | pmoa-mmox-coupling | pmoA/mmoX coupling & composition (4-compartment) | REBUILT `rev_figS04_pmoa-mmox-coupling.R` | old S10 |
| S05 | taxonomy-mcra | Taxonomy mcrA heatmap | orig | old S6 |
| S06 | taxonomy-pmoa | Taxonomy pmoA heatmap | orig | old S2 |
| S07 | faprotax | FAPROTAX heatmaps | orig | old S3 |
| S08 | picrust-mcra-all | PICRUSt mcrA-all heatmap | orig | old S4 |
| S09 | picrust-mcra-no-mcra | PICRUSt mcrA (no-mcrA) heatmap **(demoted main Fig 6)** | orig | main Fig 6 |
| S10 | picrust-pmoa | PICRUSt pmoA heatmap | orig | old S5 |
| S11 | scale-dependent-genes | Scale-dependent gene patterns, **3-column** | REBUILT `rev_figS11_scale-dependent.R` | old S11 |
| S12 | isotope-sources | Extended-isotope source composite | BUILT `rev_figS12_isotope-sources.R` | — |
| S13 | internal-gas-beeswarm | Internal-gas beeswarm | orig | old S7 |
| S14 | internal-gas-profiles | Internal-gas profiles | orig | old S8 |
| S15 | black-oak-methanome | Black-oak methanome **(revised putative defs)** | REBUILT `rev_figS15_black-oak-methanome.R` | old S12 |
| S16 | radial-sections | Tree radial sections | orig | old S13 |
| S17 | plant-traits | Plant-trait × methane-cycling heatmap (clade-robust) | NEW `rev_figS17_plant-traits.R` | — |
| S18 | rf-flux-predictions | Random-forest flux predictions | orig | old S15 |

### Retired from the SI set (moved to `exploratory/`)
- **pmoa-mmox-by-compartment** (`rev_figS07_pmoa-mmox-compartment.R`) — redundant with S04.
- **d13ch4-raincloud** (orig `figS9`) — redundant with main Fig 6.

## Assumptions / to confirm
- **Fig 9** is the new CH₄-budget figure (maps+seasonal+waterfall), replacing the old
  `fig9_upscaled_flux_seasonal`. If both are wanted, the budget becomes Fig 10 and everything holds.
- **Fig 8** is unchanged from the as-reviewed version (Jon: "I like Fig 8 as-is").
- SI numbering is sequential S01–S18 following the main-text flow (no gaps).
- **S15 (black-oak)** uses the revised methanotroph definitions (Methylacidiphilaceae
  Known→Putative), matching Fig 5 and the classification table.
- The 12 "orig" figures are pulled from the frozen original pipeline's `outputs/figures/`.
  Reproduce them by running the original pipeline first, then `run_all.R` re-assembles the set.
- Absolute-abundance figures (Fig 4, Fig 7b) still pending the ×10 flip and the soil
  dry-basis harmonization — re-render only, no renumbering.
