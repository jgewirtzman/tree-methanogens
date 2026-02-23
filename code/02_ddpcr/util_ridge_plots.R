# ==============================================================================
# Gene Abundance Ridge Plots
# ==============================================================================
# Purpose: Ridge plots showing gene abundance distributions across trees.
# ==============================================================================

library(tidyverse)
library(ggside)
library(scales)

create_gene_scatter_ggside_transformed_probe_mcra <- function(
    merged_final,
    sigma = 10.0,
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
      is_probe   = (part1 == "probe"),
      location   = if_else(is_probe, part2, part1),
      stringency = if_else(is_probe, part3, part2),
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
      ),
      sample_id = paste(tree_id, sample_type, sep = "_")
    ) %>%
    filter(stringency == "loose", !is.na(sample_type)) %>%
    filter(gene %in% c("mcrA", "pmoA", "mmoX"))
  
  # ---------- mcrA (PROBE ONLY, 1:1 — NO AVERAGING) ----------
  mcra_probe_data <- long_data %>%
    filter(gene == "mcrA", is_probe) %>%
    select(sample_id, tree_id, sample_type, mcra_copies = gene_copies)
  
  # ---------- pmoA + mmoX (sum within sample) ----------
  pmoa_mmox_data <- long_data %>%
    filter(gene %in% c("pmoA", "mmoX")) %>%
    group_by(sample_id, tree_id, sample_type) %>%
    summarise(pmoa_mmox_total = sum(gene_copies, na.rm = TRUE), .groups = "drop")
  
  # ---------- join & sanity-check uniqueness ----------
  scatter_data <- mcra_probe_data %>%
    inner_join(pmoa_mmox_data, by = c("sample_id", "tree_id", "sample_type")) %>%
    filter(mcra_copies >= 0, pmoa_mmox_total >= 0)
  
  stopifnot(!any(duplicated(scatter_data["sample_id"])))  # 1 row per sample_id
  
  # ---------- aesthetics ----------
  category_colors <- c(
    "Sapwood"   = "#dfc27d",
    "Heartwood" = "#a6611a",
    "Mineral"   = "#80cdc1",
    "Organic"   = "#018571"
  )
  legend_order <- c("Sapwood", "Heartwood", "Mineral", "Organic")
  scatter_data <- scatter_data %>% mutate(sample_type = factor(sample_type, levels = legend_order))
  
  # ---------- pseudo-log transform ----------
  trans <- scales::pseudo_log_trans(base = base, sigma = sigma)
  scatter_data <- scatter_data %>%
    mutate(
      mcra_log      = trans$transform(mcra_copies),
      pmoa_mmox_log = trans$transform(pmoa_mmox_total)
    )
  
  centroids <- scatter_data %>%
    group_by(sample_type) %>%
    summarise(
      centroid_x = mean(mcra_log, na.rm = TRUE),
      centroid_y = mean(pmoa_mmox_log, na.rm = TRUE),
      .groups = "drop"
    )
  
  # axis breaks
  raw_x_breaks <- c(0, 100, 10000, 1e6, 1e8)
  raw_y_breaks <- c(0, 100, 10000, 1e6)
  x_breaks <- trans$transform(raw_x_breaks)
  y_breaks <- trans$transform(raw_y_breaks)
  x_labels <- c("0", expression(10^2), expression(10^4), expression(10^6), expression(10^8))
  y_labels <- c("0", expression(10^2), expression(10^4), expression(10^6))
  
  x_limits <- range(scatter_data$mcra_log, na.rm = TRUE); x_limits[2] <- x_limits[2] * 1.02
  y_limits <- range(scatter_data$pmoa_mmox_log, na.rm = TRUE); y_limits[2] <- y_limits[2] * 1.02
  
  # ---------- plot ----------
  ggplot(scatter_data, aes(mcra_log, pmoa_mmox_log, color = sample_type)) +
    geom_point(size = 3, alpha = 0.5) +
    # centroids with double ring
    geom_point(
      data = centroids,
      aes(x = centroid_x, y = centroid_y, color = sample_type),
      size = 4, shape = 21, fill = alpha("white", 0.6), stroke = 0, inherit.aes = FALSE
    ) +
    geom_point(
      data = centroids,
      aes(x = centroid_x, y = centroid_y),
      size = 4, shape = 21, color = "black", fill = "transparent",
      stroke = 3.0, inherit.aes = FALSE
    ) +
    geom_point(
      data = centroids,
      aes(x = centroid_x, y = centroid_y, color = sample_type),
      size = 4, shape = 21, fill = "transparent", stroke = 2.5, inherit.aes = FALSE
    ) +
    ggside::geom_xsidedensity(
      aes(y = after_stat(density), fill = sample_type, color = NULL),
      alpha = 0.55, position = "identity", adjust = 2.0
    ) +
    ggside::geom_ysidedensity(
      aes(x = after_stat(density), fill = sample_type, color = NULL),
      alpha = 0.55, position = "identity", adjust = 2.0
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits, expand = c(0, 0)) +
    scale_color_manual(values = category_colors, breaks = legend_order, name = "Sample Type") +
    scale_fill_manual(values = category_colors, breaks = legend_order, guide = "none") +
    labs(
      x = expression("mcrA (copies g"^-1*")"),
      y = expression("pmoA + mmoX (copies g"^-1*")")
    ) +
    ggside::ggside(x.pos = "top", y.pos = "right") +
    theme_minimal(base_size = 12) +
    theme(
      ggside.panel.scale.x = 0.3,
      ggside.panel.scale.y = 0.3,
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      ggside.axis.text.x = element_blank(),
      ggside.axis.text.y = element_blank(),
      ggside.axis.ticks.x = element_blank(),
      ggside.axis.ticks.y = element_blank(),
      ggside.panel.border = element_blank(),
      # Reduce legend spacing
      legend.margin = margin(t = 5, r = 0, b = 1, l = -50),  # Remove legend margins
      legend.box.margin = margin(t = 5, r = 0, b = 0, l = -50),  # Remove legend box margins
      legend.spacing = unit(0.5, "lines"),  # Reduce spacing between legend elements
      legend.box.spacing = unit(0, "lines")  # Reduce spacing around legend box
    )
}

# Example:
p_loose_probe_only <- create_gene_scatter_ggside_transformed_probe_mcra(merged_final)
print(p_loose_probe_only)







# Calculate % of samples with nonzero mcrA by sample type
mcra_nonzero_summary <- merged_final %>%
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
    is_probe   = (part1 == "probe"),
    location   = if_else(is_probe, part2, part1),
    stringency = if_else(is_probe, part3, part2),
    location = str_remove(location, "probe"),
    location = case_when(
      location %in% c("Inner","inner")     ~ "Inner",
      location %in% c("Outer","outer")     ~ "Outer",
      location %in% c("Mineral","mineral") ~ "Mineral",
      location %in% c("Organic","organic") ~ "Organic",
      TRUE ~ location
    ),
    sample_type = case_when(
      location == "Inner"   ~ "Heartwood",
      location == "Outer"   ~ "Sapwood",
      location == "Mineral" ~ "Mineral",
      location == "Organic" ~ "Organic",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(stringency == "loose", !is.na(sample_type), gene == "mcra", is_probe) %>%
  mutate(sample_id = paste(tree_id, sample_type, sep = "_")) %>%
  select(sample_id, sample_type, mcra_copies = gene_copies) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  group_by(sample_type) %>%
  summarise(
    n_total = n(),
    n_nonzero = sum(mcra_copies > 0),
    pct_nonzero = 100 * n_nonzero / n_total,
    .groups = "drop"
  )

cat("\n% of samples with nonzero mcrA:\n")
mcra_nonzero_summary %>%
  mutate(output = sprintf("  %s: %.1f%% (%d/%d)", 
                          sample_type, pct_nonzero, n_nonzero, n_total)) %>%
  pull(output) %>%
  cat(sep = "\n")

















# Calculate number of samples above various mcrA thresholds (copies per 2µL)
# Using probe columns only

# LOOSE stringency
mcra_threshold_loose <- merged_final %>%
  select(tree_id, 
         ddpcr_mcra_probe_Inner_loose,
         ddpcr_mcra_probe_Outer_loose,
         ddpcr_mcra_probe_Mineral_loose,
         ddpcr_mcra_probe_Organic_loose) %>%
  pivot_longer(
    cols = starts_with("ddpcr_"),
    names_to = "column",
    values_to = "mcra_copies"
  ) %>%
  filter(!is.na(mcra_copies)) %>%
  mutate(
    sample_type = case_when(
      str_detect(column, "Inner") ~ "Heartwood",
      str_detect(column, "Outer") ~ "Sapwood",
      str_detect(column, "Mineral") ~ "Mineral",
      str_detect(column, "Organic") ~ "Organic"
    )
  ) %>%
  group_by(sample_type) %>%
  summarise(
    n_total = n(),
    n_gt1 = sum(mcra_copies > 1),
    pct_gt1 = 100 * n_gt1 / n_total,
    n_gt3 = sum(mcra_copies > 3),
    pct_gt3 = 100 * n_gt3 / n_total,
    n_gt5 = sum(mcra_copies > 5),
    pct_gt5 = 100 * n_gt5 / n_total,
    n_gt10 = sum(mcra_copies > 10),
    pct_gt10 = 100 * n_gt10 / n_total,
    .groups = "drop"
  )

# STRICT stringency
mcra_threshold_strict <- merged_final %>%
  select(tree_id,
         ddpcr_mcra_probe_Inner_strict,
         ddpcr_mcra_probe_Outer_strict,
         ddpcr_mcra_probe_Mineral_strict,
         ddpcr_mcra_probe_Organic_strict) %>%
  pivot_longer(
    cols = starts_with("ddpcr_"),
    names_to = "column",
    values_to = "mcra_copies"
  ) %>%
  filter(!is.na(mcra_copies)) %>%
  mutate(
    sample_type = case_when(
      str_detect(column, "Inner") ~ "Heartwood",
      str_detect(column, "Outer") ~ "Sapwood",
      str_detect(column, "Mineral") ~ "Mineral",
      str_detect(column, "Organic") ~ "Organic"
    )
  ) %>%
  group_by(sample_type) %>%
  summarise(
    n_total = n(),
    n_gt1 = sum(mcra_copies > 1),
    pct_gt1 = 100 * n_gt1 / n_total,
    n_gt3 = sum(mcra_copies > 3),
    pct_gt3 = 100 * n_gt3 / n_total,
    n_gt5 = sum(mcra_copies > 5),
    pct_gt5 = 100 * n_gt5 / n_total,
    n_gt10 = sum(mcra_copies > 10),
    pct_gt10 = 100 * n_gt10 / n_total,
    .groups = "drop"
  )

cat("\n=== LOOSE stringency: mcrA above thresholds (copies per 2µL) ===\n")
mcra_threshold_loose %>%
  mutate(output = sprintf("  %s: >1: %d (%.1f%%), >3: %d (%.1f%%), >5: %d (%.1f%%), >10: %d (%.1f%%) [n=%d]", 
                          sample_type, 
                          n_gt1, pct_gt1,
                          n_gt3, pct_gt3,
                          n_gt5, pct_gt5,
                          n_gt10, pct_gt10,
                          n_total)) %>%
  pull(output) %>%
  cat(sep = "\n")

cat("\n\n=== STRICT stringency: mcrA above thresholds (copies per 2µL) ===\n")
mcra_threshold_strict %>%
  mutate(output = sprintf("  %s: >1: %d (%.1f%%), >3: %d (%.1f%%), >5: %d (%.1f%%), >10: %d (%.1f%%) [n=%d]", 
                          sample_type, 
                          n_gt1, pct_gt1,
                          n_gt3, pct_gt3,
                          n_gt5, pct_gt5,
                          n_gt10, pct_gt10,
                          n_total)) %>%
  pull(output) %>%
  cat(sep = "\n")
