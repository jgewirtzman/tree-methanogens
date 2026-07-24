# Final figure manifest (revision)

Clean, correctly-numbered copies of the final figures. `main/` = main text only,
`SI/` = supplement only. Regenerate any file from the script in the Source column.
Numbers reflect the agreed restructure (hydrogenotrophy → main Fig 6; old PICRUSt-mcrA
heatmap → SI; SI resequenced per `SI_FIGURES_PLAN.md`).

## Main figures
| # | File | Content | Source script | Source PNG |
|---|---|---|---|---|
| 1 | Figure_1_temporal-flux | Temporal flux time series (soil + tree) | `R_fig1_final.R` | fig1_final.png |
| 2 | Figure_2_height-flux | Height-dependent flux | `R_fig2_final.R` | fig2_final.png |
| 3 | Figure_3_variance-partition | Variance partitioning | `R_fig3_final.R` | fig3_final.png |
| 4 | Figure_4_gene-abundance | mcrA / methanotroph abundance by compartment (copies/g; ×10 toggle) | `R_fig4_final.R` | fig4_final.png |
| 5 | Figure_5_methane-cycling-composition | Combined methane-cycling composition (revised putative) | `R_fig5_final.R` | fig5_final.png |
| 6 | Figure_6_hydrogenotrophy | Convergent-hydrogenotrophy synthesis **(NEW — replaces old PICRUSt-mcrA heatmap)** | `R_fig_hydrogenotrophy.R` | fig_hydrogenotrophy.png |
| 7 | Figure_7_felled-oak | Felled-oak profiles (copies/g + flux unit) | `R_fig7_final.R` | fig7_felled_oak_final.png |
| 8 | Figure_8_radial-species | Radial species comparison **(kept as-is)** | `05_gene_flux_analysis/04_species_gene_flux.R` | outputs/figures/main/fig8_radial_species_comparison.png |
| 9 | Figure_9_ch4-budget | CH₄ budget: maps + seasonal + net waterfall **(NEW — replaces old `fig9_upscaled_flux_seasonal`)** | `R_fig_component_budget.R` | fig_budget_maps.png |

## SI figures
| # | File | Content | Source | Was |
|---|---|---|---|---|
| S01 | moisture-overlay | Moisture overlay | orig | old S1 |
| S02 | height-slope-moisture | Height-decline slope vs soil moisture (R3) | NEW `R_height_slope_vs_moisture` | — |
| S03 | rf-flux-predictions | Random-forest flux predictions | orig | old S15 |
| S04 | scale-dependent-genes | Scale-dependent gene patterns, **3-column** | REBUILT `R_figS11_final.R` | old S11 |
| S05 | mcra-vs-methanotroph | mcrA vs methanotroph | orig | old S14 |
| S06 | pmoa-mmox-coupling | pmoA/mmoX coupling & composition (4-compartment) | REBUILT `R_figS10_final.R` | old S10 |
| S07 | pmoa-mmox-by-compartment | pmoA vs mmoX by compartment (abundances) | NEW `R_pmoa_mmox_separate.R` | — |
| S08 | taxonomy-mcra | Taxonomy mcrA heatmap | orig | old S6 |
| S09 | taxonomy-pmoa | Taxonomy pmoA heatmap | orig | old S2 |
| S10 | faprotax | FAPROTAX heatmaps | orig | old S3 |
| S11 | picrust-mcra-all | PICRUSt mcrA-all heatmap | orig | old S4 |
| S12 | picrust-mcra-no-mcra | PICRUSt mcrA (no-mcrA) heatmap **(demoted main Fig 6)** | orig | main Fig 6 |
| S13 | picrust-pmoa | PICRUSt pmoA heatmap | orig | old S5 |
| S14 | d13ch4-raincloud | δ¹³CH₄ raincloud | orig | old S9 |
| S15 | isotope-sources | Extended-isotope source composite | BUILT `R_SI_isotopes_composite.R` | — |
| S16 | internal-gas-beeswarm | Internal-gas beeswarm | orig | old S7 |
| S17 | internal-gas-profiles | Internal-gas profiles | orig | old S8 |
| S18 | black-oak-methanome | Black-oak methanome | orig | old S12 |
| S19 | radial-sections | Tree radial sections | orig | old S13 |
| S20 | plant-traits | Plant-trait × methane-cycling heatmap (clade-robust) | NEW `R_traits_heatmap_robust.R` | — |

## Assumptions / to confirm
- **Fig 9** is the new CH₄-budget figure (maps+seasonal+waterfall), replacing the old
  `fig9_upscaled_flux_seasonal`. If both are wanted, the budget becomes Fig 10 and everything holds.
- **Fig 8** is unchanged from the as-reviewed version (Jon: "I like Fig 8 as-is").
- SI numbering is sequential S01–S20 following the main-text flow in `SI_FIGURES_PLAN.md`
  (the retired regfig-combined row is dropped; no gaps).
- Absolute-abundance figures (Fig 4, Fig 7b, S07) still pending the ×10 flip and the soil
  dry-basis harmonization — re-render only, no renumbering.
