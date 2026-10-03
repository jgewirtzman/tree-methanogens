source("code/lib/outputs.R")
# ==============================================================================
# REVISION — Fig 4 final: methanogen (mcrA) / methanotroph (pmoA,mmoX) gene
# abundance by compartment (heartwood/sapwood/soil), (a) species barplot + (b) scatter.
# Copy of the original generator (code/archive/superseded_2026-10-01/08_figures/util_combined_plot.R).
#
# UNITS: the ddpcr_*_loose columns are ALREADY copies g^-1, converted once in
# code/03_merge/04_harmonize_all_data.R with code/lib/ddpcr_constants.R
# (Conc x (25/2.5) x 75 uL / mass). Basis: DRY for wood (freeze-dried cores); soil
# uses fresh sample mass.
# Output: outputs/figures/generated/fig4_final.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(cowplot); library(patchwork) })
out <- "outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

merged_final <- read_csv("data/processed/integrated/merged_tree_dataset_final.csv", show_col_types = FALSE)

source("code/08_figures/helper_species_barplots.R")   # defines create_mcra_barplot_by_species + species_mapping
source("code/07_molecular/helper_ridge_plots.R")       # defines create_gene_scatter_ggside_transformed_probe_mcra

result      <- create_mcra_barplot_by_species(merged_final, species_mapping)
scatterplot <- create_gene_scatter_ggside_transformed_probe_mcra(merged_final)
barplot_improved <- result$plot +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), panel.spacing = unit(1, "lines"))
combined_plot <- plot_grid(barplot_improved, scatterplot, ncol = 2,
                           labels = c("(a)", "(b)"), label_size = 11, label_fontface = "bold",
                           rel_widths = c(1.3, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(out_path("fig4_final.png"), combined_plot, width = 12, height = 7, dpi = 300, bg = "white")
cat("Wrote fig4_final.png\n")
