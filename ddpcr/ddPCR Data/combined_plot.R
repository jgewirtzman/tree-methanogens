# Generate your plots
result <- create_mcra_barplot_by_species(merged_final, species_mapping)
barplot <- result$plot

scatterplot <- create_gene_scatter_ggside_transformed_probe_mcra(merged_final)


# Option 1: Smaller text and better spacing
barplot_improved <- result$plot +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),  # Smaller text
    panel.spacing = unit(1, "lines")  # More space between facets
  )



# Then use the improved version:
combined_plot <- plot_grid(barplot_improved, scatterplot, ncol = 2)

combined_plot