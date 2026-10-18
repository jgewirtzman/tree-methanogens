# ==============================================================================
# species_gene_data.R -- shared species/gene/flux data prep (sourced)
# ==============================================================================
# The ONE copy of the tree- and species-level gene/flux prep, reading the canonical
# tables through code/lib/tables_read.R. Sourced by
# code/07_molecular/helper_scale_dependent_gene_patterns.R (and through it 14-16, 18,
# 19, Figs 4, 8, S21), code/07_molecular/helper_radial_gene_plots.R (the parser
# only), 09_tables_stats/17, 20 and 25, and Fig. S28. Builds, in the calling
# environment:
#   area_weighted_gene(), prepare_long_all_genes()
#   ymf2021, ymf2023     - survey_2021_trees() and the 2023 cross-species deployments
#   long_all             - wood genes, long (loose threshold; mcrA probe assay only)
#   tree_genes_weighted  - per-tree area-weighted genes (trees with a measured DBH)
#   tree_genes_complete  - ... with mcrA, pmoA and mmoX all present
#   tree_level_complete  - ... + 125 cm flux, log10(x + 1) genes and log_ratio
#   flux_all             - individual chamber CH4 fluxes (pooled 2021 125 cm + 2023) + species
#   flux_by_species      - species n_flux, median, quartiles of flux
#   analysis_ratio       - species-level: median_log_ratio, median_flux (n_trees>=5, n_flux>=5)
#   analysis_mcra        - species-level area-weighted mcrA (col 'value')
#   analysis_meth        - species-level methanotroph_total (col 'value')
# Side-effect free (no writes, no plots). Area weighting needs the tree's DBH: trees
# without a measured DBH (trees.csv dbh_cm NA) have no area-weighted value.
#
# Inputs (data/tables/): trees, species, flux_stem (2021_multiheight at 125 cm;
#   2023_cross_species), ddpcr + samples
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))
source("code/lib/tables_read.R")

ymf2021 <- survey_2021_trees()
ymf2023 <- flux_stem_campaign("2023_cross_species")

area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  # Cores were cut into an outer 5-cm segment (sapwood) and an inner segment (heartwood):
  # sapwood value in the outer 5 cm, heartwood value within min(5 cm, R - 5 cm) of the pith,
  # linear in radius between. A tree of radius <= 5 cm is all sapwood.
  r1 <- min(5, max(R - 5, 0)); r2 <- max(R - 5, r1)
  S <- r1^2 + r1 * r2 + r2^2
  gene_outer + (gene_inner - gene_outer) * (S / (3 * R^2))
}

prepare_long_all_genes <- function(df) {
  df %>%
    dplyr::select(tree_id, starts_with("ddpcr_")) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    pivot_longer(cols = starts_with("ddpcr_"),
                 names_to = "measurement_type", values_to = "gene_copies") %>%
    filter(!is.na(gene_copies)) %>%
    separate(measurement_type,
             into = c("method", "gene", "part1", "part2", "part3"),
             sep = "_", extra = "merge", fill = "right") %>%
    mutate(
      is_probe = (part1 == "probe"),
      location = if_else(is_probe, part2, part1),
      stringency = if_else(is_probe, part3, part2),
      location = str_remove(location, "probe"),
      location = case_when(location %in% c("Inner", "inner") ~ "Inner",
                           location %in% c("Outer", "outer") ~ "Outer",
                           TRUE ~ location),
      gene = case_when(gene == "mcra" ~ "mcrA", gene == "mmox" ~ "mmoX",
                       gene == "pmoa" ~ "pmoA", TRUE ~ gene),
      sample_type = case_when(location == "Inner" ~ "Heartwood",
                              location == "Outer" ~ "Sapwood",
                              TRUE ~ NA_character_)
    ) %>%
    filter(stringency == "loose", !is.na(sample_type),
           sample_type %in% c("Heartwood", "Sapwood"),
           gene != "mcrA" | is_probe)  # mcrA: probe assay only (the EvaGreen mcrA columns give false positives in wood; Fig. S4)
}

long_all <- prepare_long_all_genes(ymf2021)

tree_genes_weighted <- long_all %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
  filter(is.finite(dbh)) %>%
  group_by(tree_id, species_id, dbh, gene) %>%
  summarise(gene_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
            gene_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
            .groups = "drop") %>%
  filter(is.finite(gene_inner), is.finite(gene_outer)) %>%
  mutate(gene_area_weighted = mapply(area_weighted_gene, gene_inner, gene_outer, dbh),
         species = species_name(species_id)) %>%
  dplyr::select(tree_id, species_id, species, gene, gene_area_weighted) %>%
  pivot_wider(names_from = gene, values_from = gene_area_weighted) %>%
  mutate(methanotroph_total = pmoA + mmoX)   # NA where either gene is missing

tree_genes_complete <- tree_genes_weighted %>%
  filter(!is.na(mcrA), !is.na(pmoA), !is.na(mmoX), !is.na(methanotroph_total))

flux_all <- bind_rows(
  ymf2023 %>% dplyr::select(species_id = species_code, CH4_flux = CH4_flux_nmol_m2_s) %>%
    mutate(species = species_name(species_id)),
  ymf2021 %>% dplyr::select(species_id, CH4_flux = CH4_best.flux_125cm) %>%
    mutate(species = species_name(species_id))
) %>% filter(!is.na(CH4_flux), !is.nan(CH4_flux), !is.na(species))

flux_by_species <- flux_all %>%
  group_by(species, species_id) %>%
  summarise(n_flux = n(), median_flux = median(CH4_flux, na.rm = TRUE),
            q25_flux = quantile(CH4_flux, 0.25, na.rm = TRUE),
            q75_flux = quantile(CH4_flux, 0.75, na.rm = TRUE), .groups = "drop")

tree_level_complete <- tree_genes_complete %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
            by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux)) %>%
  mutate(log_tree_mcra = log10(mcrA + 1),
         log_tree_mmox = log10(mmoX + 1),
         log_pmoa = log10(pmoA + 1),
         log_methanotroph = log10(methanotroph_total + 1),
         ratio_mcra_methanotroph = (mcrA + 1) / (methanotroph_total + 1),
         log_ratio = log10(ratio_mcra_methanotroph))

.species_agg <- function(col) {
  tree_level_complete %>%
    group_by(species, species_id) %>%
    summarise(n_trees = n(), value = median(.data[[col]], na.rm = TRUE), .groups = "drop") %>%
    inner_join(flux_by_species, by = c("species", "species_id")) %>%
    filter(n_trees >= 5, n_flux >= 5)
}

analysis_ratio <- tree_level_complete %>%
  group_by(species, species_id) %>%
  summarise(n_trees = n(), median_log_ratio = median(log_ratio, na.rm = TRUE), .groups = "drop") %>%
  inner_join(flux_by_species, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

analysis_mcra <- .species_agg("mcrA")
analysis_meth <- .species_agg("methanotroph_total")
