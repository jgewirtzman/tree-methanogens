# ==============================================================================
# REVISION HELPER — shared species/gene/flux data prep (NEW file; sourced)
# ==============================================================================
# Reproduces the prep chain from
# code/05_gene_flux_analysis/02_scale_dependent_gene_patterns.R (lines 50-302),
# VERBATIM logic, so numbers match the as-reviewed manuscript
# (ratio R2 = 0.513, r = 0.717). Builds, in the calling environment:
#   tree_level_complete  - per-tree area-weighted genes + 125cm flux + log_ratio
#   flux_all             - individual chamber CH4 fluxes (pooled 2021+2023) + species
#   flux_by_species      - species median flux (+ n_flux)
#   analysis_ratio       - species-level: median_log_ratio, median_flux (n_trees>=5, n_flux>=5)
#   analysis_mcra        - species-level area-weighted mcrA (col 'value')
#   analysis_meth        - species-level methanotroph_total (col 'value')
# Side-effect free (no writes, no plots). Sourced by revision analysis scripts.
#
# Inputs:
#   data/processed/flux/methanogen_tree_flux_complete_dataset.csv
#   data/processed/integrated/merged_tree_dataset_final.csv
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

ymf2023 <- read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv")
ymf2021 <- read.csv("data/processed/integrated/merged_tree_dataset_final.csv")

species_mapping <- c(
  "ACRU" = "Acer rubrum", "ACSA" = "Acer saccharum",
  "BEAL" = "Betula alleghaniensis", "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera", "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana", "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra", "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata", "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina", "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina", "SAAL" = "Sassafras albidum"
)

area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  r1 <- min(5, R); r2 <- max(R - 5, r1)
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
           sample_type %in% c("Heartwood", "Sapwood"))
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
         species = species_mapping[species_id]) %>%
  dplyr::select(tree_id, species_id, species, gene, gene_area_weighted) %>%
  pivot_wider(names_from = gene, values_from = gene_area_weighted) %>%
  mutate(methanotroph_total = case_when(
    is.na(pmoA) & !is.na(mmoX) ~ mmoX,
    !is.na(pmoA) & is.na(mmoX) ~ pmoA,
    !is.na(pmoA) & !is.na(mmoX) ~ pmoA + mmoX,
    TRUE ~ NA_real_))

tree_genes_complete <- tree_genes_weighted %>%
  filter(!is.na(mcrA), !is.na(pmoA), !is.na(mmoX), !is.na(methanotroph_total))

flux_all <- bind_rows(
  ymf2023 %>% dplyr::select(species_id = Species.Code, CH4_flux = CH4_best.flux) %>%
    mutate(species = species_mapping[species_id]),
  ymf2021 %>% dplyr::select(species_id, CH4_flux = CH4_best.flux_125cm) %>%
    mutate(species = species_mapping[species_id])
) %>% filter(!is.na(CH4_flux), !is.nan(CH4_flux), !is.na(species))

flux_by_species <- flux_all %>%
  group_by(species, species_id) %>%
  summarise(n_flux = n(), median_flux = median(CH4_flux, na.rm = TRUE), .groups = "drop")

tree_level_complete <- tree_genes_complete %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm),
            by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux)) %>%
  mutate(ratio_mcra_methanotroph = (mcrA + 1) / (methanotroph_total + 1),
         log_ratio = log10(ratio_mcra_methanotroph),
         log_mcra = log10(mcrA + 1),
         log_meth = log10(methanotroph_total + 1))

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
