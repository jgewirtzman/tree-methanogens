# ==============================================================================
# REVISION — Fig 4 final: methanogen (mcrA) / methanotroph (pmoA,mmoX) gene
# abundance by compartment (heartwood/sapwood/soil), (a) species barplot + (b) scatter.
# Copy of the original generator (code/02_ddpcr/util_combined_plot.R) with a
# DILUTION_10X toggle (default 1) on the absolute ddPCR copies — flip to 10 if Wyatt
# confirms the dropped template->reaction dilution (A1). Original scripts untouched.
#
# UNITS: the ddpcr_*_loose columns are ALREADY copies g^-1 — the harmonization step
# (code/00_harmonization/02_harmonize_all_data.R:317) applies concentration_per_g =
# copies/uL * 75 / sample_mass * 1000 (elution / sample mass). Basis: DRY for wood
# (freeze-dried cores); soil uses fresh sample mass (the pending A2 dry-basis
# harmonization would move soil to a dry basis and needs soil moisture we don't have).
# The DILUTION_10X toggle multiplies these copies/g values, i.e. the A1 x10 correction
# (uniform log10 shift; relative/among-compartment comparisons unchanged).
# Output: outputs/revision/fig4_final.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(cowplot); library(patchwork) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
DILUTION_10X <- 1   # set to 10 when Wyatt confirms the dropped dilution

merged_final <- read_csv("data/processed/integrated/merged_tree_dataset_final.csv", show_col_types = FALSE)
ddpcr_cols <- grep("^ddpcr_.*_loose$", names(merged_final), value = TRUE)
merged_final[ddpcr_cols] <- merged_final[ddpcr_cols] * DILUTION_10X   # x10 toggle (absolute copies only)

source("code/02_ddpcr/04_species_barplots.R")   # defines create_mcra_barplot_by_species + species_mapping
source("code/02_ddpcr/util_ridge_plots.R")       # defines create_gene_scatter_ggside_transformed_probe_mcra

result      <- create_mcra_barplot_by_species(merged_final, species_mapping)
scatterplot <- create_gene_scatter_ggside_transformed_probe_mcra(merged_final)
barplot_improved <- result$plot +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), panel.spacing = unit(1, "lines"))
combined_plot <- plot_grid(barplot_improved, scatterplot, ncol = 2,
                           labels = c("(a)", "(b)"), label_size = 11, label_fontface = "bold",
                           rel_widths = c(1.3, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(out, "fig4_final.png"), combined_plot, width = 12, height = 7, dpi = 300, bg = "white")
cat("Wrote fig4_final.png (DILUTION_10X =", DILUTION_10X, ")\n")
