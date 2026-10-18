# ==============================================================================
# helper_radial_gene_plots.R -- species-averaged stem cross-sections (Figure 8a-c)
# ------------------------------------------------------------------------------
# SOURCED by 08_figures/fig08_radial-species.R (from the repository root); never run
# alone. For each species with >= 5 trees, a grid over a disc of the species' median
# DBH: heartwood value within 5 cm of the centre, sapwood value within 5 cm of the
# bark, log-linear in between (log10(copies g-1 + 1); species means of tree means).
#   mcra_cloud      mcrA (probe assay)
#   sum_cloud       pmoA + mmoX (where both were measured in a compartment)
#   overlay_cloud   mcrA and pmoA + mmoX on one grid (Clog_mcra, Clog_methan)
# Input: data/tables/ via code/lib/species_gene_data.R (wood genes, loose threshold).
# ==============================================================================

# Load libraries
library(dplyr)
library(tidyr)

# ============================================================
# LOAD DATA
# ============================================================

# 2021 survey per tree (species and DBH from trees.csv) and its wood genes, long
# (loose threshold; mcrA probe assay only): code/lib/species_gene_data.R, sourced into
# its own environment so the species-level objects of
# helper_scale_dependent_gene_patterns.R (sourced alongside in fig08) are not replaced.
.sg <- new.env()
sys.source("code/lib/species_gene_data.R", envir = .sg)
ymf2021 <- .sg$ymf2021

# ============================================================
# GENE DATA: mcrA, pmoA, mmoX
# ============================================================

tree_genes <- .sg$long_all %>% dplyr::filter(gene %in% c("mcrA", "pmoA", "mmoX"))

cat(sprintf("Total measurements: %d\n", nrow(tree_genes)))
print(table(tree_genes$gene))

# ============================================================
# GRID GENERATION FUNCTION
# ============================================================

generate_species_grid <- function(species_row, gene_inner, gene_outer, grid_resolution = 0.1) {
  # Use median DBH for this species
  R <- species_row$median_dbh / 2
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  
  x_seq <- seq(-R, R, by = grid_resolution)
  y_seq <- seq(-R, R, by = grid_resolution)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$r <- sqrt(grid$x^2 + grid$y^2)
  grid <- grid[grid$r <= R + grid_resolution/2, ]  # Include edge points
  
  # Log-scale interpolation
  log_inner <- log10(gene_inner + 1)
  log_outer <- log10(gene_outer + 1)
  
  grid$Clog <- ifelse(grid$r <= r1, log_inner,
                      ifelse(grid$r >= r2, log_outer,
                             log_inner + (log_outer - log_inner) *
                               (grid$r - r1) / max(r2 - r1, 1e-9)))
  
  grid$species_label <- species_row$species_label
  grid$R <- R
  
  grid
}

# ============================================================
# SPECIES-AVERAGED CROSS-SECTION (one grid per species)
# ============================================================

species_gene_cloud <- function(gene_data, gene_name, min_trees = 5) {
  
  # Calculate species averages
  # Check if species_id already exists in the data
  if(!"species_id" %in% names(gene_data)) {
    gene_data <- gene_data %>%
      dplyr::left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id")
  }
  
  species_data <- gene_data %>%
    dplyr::filter(gene == gene_name) %>%
    dplyr::mutate(species = species_name(species_id)) %>%
    dplyr::filter(is.finite(dbh)) %>%
    dplyr::group_by(tree_id, species_id, species, dbh, sample_type) %>%
    dplyr::summarise(gene_copies = mean(gene_copies, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample_type, values_from = gene_copies) %>%
    dplyr::group_by(species_id, species) %>%
    dplyr::summarise(
      n_trees = dplyr::n(),
      gene_inner = mean(Heartwood, na.rm = TRUE),
      gene_outer = mean(Sapwood, na.rm = TRUE),
      median_dbh = median(dbh, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_trees >= min_trees, is.finite(gene_inner), is.finite(gene_outer))
  
  # Add species label
  species_data$species_label <- ifelse(is.na(species_data$species) | species_data$species == "", 
                                       species_data$species_id, species_data$species)
  
  stopifnot(nrow(species_data) > 0)
  
  # Order alphabetically by species name
  species_data <- species_data %>%
    dplyr::mutate(mean_gene = (gene_inner + gene_outer) / 2) %>%
    dplyr::arrange(species_label)
  
  species_data$species_label <- factor(species_data$species_label,
                                       levels = species_data$species_label)
  
  cat(sprintf("\n%s: %d species (n≥%d trees)\n", 
              gene_name, nrow(species_data), min_trees))
  
  # Generate grids for all species
  all_cloud <- dplyr::bind_rows(lapply(1:nrow(species_data), function(i) {
    generate_species_grid(
      species_data[i,],
      species_data$gene_inner[i],
      species_data$gene_outer[i]
    )
  }))
  
  all_cloud
}

# ============================================================
# CLOUDS
# ============================================================
mcra_cloud <- species_gene_cloud(tree_genes, "mcrA")

# pmoA + mmoX per tree and compartment
tree_species <- ymf2021 %>% 
  dplyr::select(tree_id, species_id, dbh) %>%
  dplyr::distinct()

sum_data <- tree_genes %>%
  dplyr::filter(gene %in% c("pmoA", "mmoX")) %>%
  dplyr::left_join(tree_species, by = "tree_id") %>%
  dplyr::filter(is.finite(dbh)) %>%
  dplyr::group_by(tree_id, species_id, dbh, sample_type) %>%
  # pmoA + mmoX only where both genes were measured in that compartment (one missing -> NA)
  dplyr::summarise(gene_copies = if (dplyr::n_distinct(gene) == 2) sum(gene_copies) else NA_real_,
                   .groups = "drop") %>%
  dplyr::filter(!is.na(gene_copies)) %>%
  dplyr::mutate(gene = "pmoA+mmoX")

sum_cloud <- species_gene_cloud(sum_data, "pmoA+mmoX")

# Overlay: species with >= 5 trees for both mcrA and pmoA + mmoX, one grid each
mcra_species <- tree_genes %>%
  dplyr::filter(gene == "mcrA") %>%
  dplyr::left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
  dplyr::group_by(tree_id, species_id, dbh, sample_type) %>%
  dplyr::summarise(gene_copies = mean(gene_copies, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = sample_type, values_from = gene_copies) %>%
  dplyr::group_by(species_id) %>%
  dplyr::summarise(
    n_trees = dplyr::n(),
    mcra_inner = mean(Heartwood, na.rm = TRUE),
    mcra_outer = mean(Sapwood, na.rm = TRUE),
    median_dbh = median(dbh, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_trees >= 5)

methan_species <- tree_genes %>%
  dplyr::filter(gene %in% c("pmoA", "mmoX")) %>%
  dplyr::left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
  dplyr::group_by(tree_id, species_id, dbh, sample_type) %>%
  # pmoA + mmoX only where both genes were measured in that compartment (one missing -> NA)
  dplyr::summarise(gene_copies = if (dplyr::n_distinct(gene) == 2) sum(gene_copies) else NA_real_,
                   .groups = "drop") %>%
  dplyr::filter(!is.na(gene_copies)) %>%
  tidyr::pivot_wider(names_from = sample_type, values_from = gene_copies) %>%
  dplyr::group_by(species_id) %>%
  dplyr::summarise(
    n_trees = dplyr::n(),
    methan_inner = mean(Heartwood, na.rm = TRUE),
    methan_outer = mean(Sapwood, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_trees >= 5)

overlay_species <- mcra_species %>%
  dplyr::inner_join(methan_species, by = "species_id", suffix = c("_mcra", "_methan")) %>%
  dplyr::mutate(species = species_name(species_id))
overlay_species$species_label <- ifelse(is.na(overlay_species$species) | overlay_species$species == "",
                                        overlay_species$species_id, overlay_species$species)
cat(sprintf("Found %d species with both mcrA and methanotrophs\n", nrow(overlay_species)))
stopifnot(nrow(overlay_species) > 0)

ov_grid <- function(inner, outer, col) dplyr::bind_rows(lapply(1:nrow(overlay_species), function(i) {
  grid <- generate_species_grid(overlay_species[i,], overlay_species[[inner]][i], overlay_species[[outer]][i])
  grid[[col]] <- grid$Clog
  grid
}))
overlay_cloud <- ov_grid("mcra_inner", "mcra_outer", "Clog_mcra") %>%
  dplyr::select(x, y, species_label, R, Clog_mcra) %>%
  dplyr::inner_join(
    ov_grid("methan_inner", "methan_outer", "Clog_methan") %>% dplyr::select(x, y, species_label, Clog_methan),
    by = c("x", "y", "species_label"))
