# install.packages(c("tidyverse","ggside","scales")) # if needed
library(tidyverse)
library(ggside)
library(scales)

# Scatter + OUTWARD marginal densities, with EVERYTHING in pseudolog space
create_gene_scatter_ggside_transformed <- function(
    merged_final,
    stringency_type = "loose",
    sigma = 10.0,   # GREATLY INCREASED for much more expansion near zero
    base  = 10
) {
  # ---------- reshape to long ----------
  long_data <- merged_final %>%
    dplyr::select(tree_id, starts_with("ddpcr_")) %>%
    dplyr::select_if(~ !all(is.na(.))) %>%
    pivot_longer(
      cols = starts_with("ddpcr_"),
      names_to = "measurement_type",
      values_to = "gene_copies"
    ) %>%
    filter(!is.na(gene_copies)) %>%
    separate(
      measurement_type,
      into = c("method", "gene", "part1", "part2", "part3"),
      sep = "_", extra = "merge", fill = "right"
    ) %>%
    mutate(
      location   = if_else(part1 == "probe", part2, part1),
      stringency = if_else(part1 == "probe", part3, part2),
      location = str_remove(location, "probe"),
      location = case_when(
        location %in% c("Inner","inner")     ~ "Inner",
        location %in% c("Outer","outer")     ~ "Outer",
        location %in% c("Mineral","mineral") ~ "Mineral",
        location %in% c("Organic","organic") ~ "Organic",
        TRUE ~ location
      ),
      gene = case_when(
        gene == "mcra" ~ "mcrA",
        gene == "mmox" ~ "mmoX",
        gene == "pmoa" ~ "pmoA",
        TRUE ~ gene
      ),
      sample_type = case_when(
        location == "Inner"   ~ "Heartwood",
        location == "Outer"   ~ "Sapwood",
        location == "Mineral" ~ "Mineral",
        location == "Organic" ~ "Organic",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(stringency == stringency_type, !is.na(sample_type)) %>%
    filter(gene %in% c("mcrA", "pmoA", "mmoX")) %>%
    mutate(sample_id = paste(tree_id, sample_type, sep = "_"))
  
  # ---------- wide for scatter ----------
  mcra_data <- long_data %>%
    filter(gene == "mcrA") %>%
    dplyr::select(sample_id, tree_id, sample_type, mcra_copies = gene_copies)
  
  pmoa_mmox_data <- long_data %>%
    filter(gene %in% c("pmoA", "mmoX")) %>%
    group_by(sample_id, tree_id, sample_type) %>%
    summarise(pmoa_mmox_total = sum(gene_copies, na.rm = TRUE), .groups = "drop")
  
  scatter_data <- mcra_data %>%
    inner_join(pmoa_mmox_data, by = c("sample_id", "tree_id", "sample_type")) %>%
    filter(mcra_copies >= 0, pmoa_mmox_total >= 0)
  
  # ---------- palette & legend order ----------
  category_colors <- c(
    "Sapwood"   = "#dfc27d",
    "Heartwood" = "#a6611a",
    "Mineral"   = "#80cdc1",
    "Organic"   = "#018571"
  )
  legend_order <- c("Sapwood", "Heartwood", "Mineral", "Organic")
  scatter_data <- scatter_data %>% mutate(sample_type = factor(sample_type, levels = legend_order))
  
  # ---------- transform EVERYTHING to pseudolog space ----------
  trans <- scales::pseudo_log_trans(base = base, sigma = sigma)
  
  scatter_data <- scatter_data %>%
    mutate(
      mcra_log     = trans$transform(mcra_copies),
      pmoa_mmox_log = trans$transform(pmoa_mmox_total)
    )
  
  # centroids in transformed space
  centroids <- scatter_data %>%
    group_by(sample_type) %>%
    summarise(
      centroid_x = mean(mcra_log, na.rm = TRUE),
      centroid_y = mean(pmoa_mmox_log, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Updated breaks to show only powers of 10
  raw_x_breaks <- c(0, 0.1, 1, 10, 100, 1000, 100000)
  raw_y_breaks <- c(0, 0.1, 1, 10, 100, 1000)
  x_breaks <- trans$transform(raw_x_breaks)
  y_breaks <- trans$transform(raw_y_breaks)
  x_labels <- raw_x_breaks
  y_labels <- raw_y_breaks
  
  # limits on transformed axes
  x_limits <- range(scatter_data$mcra_log, na.rm = TRUE); x_limits[2] <- x_limits[2] * 1.02
  y_limits <- range(scatter_data$pmoa_mmox_log, na.rm = TRUE); y_limits[2] <- y_limits[2] * 1.02
  
  # ---------- plot (linear axes; data already transformed) ----------
  ggplot(scatter_data, aes(mcra_log, pmoa_mmox_log, color = sample_type)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_point(
      data = centroids,
      aes(x = centroid_x, y = centroid_y, color = sample_type),
      size = 6, shape = 3, stroke = 2, inherit.aes = FALSE
    ) +
    # side densities now run on transformed values -> true pseudolog KDE with smoothing
    ggside::geom_xsidedensity(
      aes(y = after_stat(density), fill = sample_type, color = NULL),
      alpha = 0.55, position = "identity", adjust = 2.0
    ) +
    ggside::geom_ysidedensity(
      aes(x = after_stat(density), fill = sample_type, color = NULL),
      alpha = 0.55, position = "identity", adjust = 2.0
    ) +
    # linear scales with transformed breaks/limits but raw labels
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits, expand = c(0, 0)) +
    scale_color_manual(values = category_colors, breaks = legend_order, name = "Sample Type") +
    scale_fill_manual(values = category_colors, breaks = legend_order, guide = "none") +
    labs(
      x = "mcrA copies (pseudolog ticks)",
      y = "pmoA + mmoX copies (pseudolog ticks)"
    ) +
    ggside::ggside(x.pos = "top", y.pos = "right") +
    theme_minimal(base_size = 12) +
    theme(
      ggside.panel.scale.x = 0.3,
      ggside.panel.scale.y = 0.3,
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),  # Remove all grid lines from main plot
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add border to main plot
      ggside.axis.text.x = element_blank(),  # Remove x-axis labels on side plots
      ggside.axis.text.y = element_blank(),  # Remove y-axis labels on side plots
      ggside.axis.ticks.x = element_blank(), # Remove x-axis ticks on side plots
      ggside.axis.ticks.y = element_blank(), # Remove y-axis ticks on side plots
      ggside.panel.border = element_blank()  # Remove borders from side plots
    )
}

# Example usage focused on loose stringency only:
p_loose_expanded <- create_gene_scatter_ggside_transformed(merged_final, "loose", sigma = 0.01)

# Alternative with even more expansion if needed:
p_loose_very_expanded <- create_gene_scatter_ggside_transformed(merged_final, "loose", sigma = 20.0)

print(p_loose_expanded)