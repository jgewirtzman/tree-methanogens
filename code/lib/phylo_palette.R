# ==============================================================================
# phylo_palette.R -- species colours in phylogenetic order (Figs S16-S17)
# ------------------------------------------------------------------------------
# Species are placed on the PhytoPhylo megatree (data/raw/external/phylogeny/,
# Qian & Jin 2016), ordered by mean cophenetic distance to the other species
# present, cut into up to 8 quantile groups, and coloured viridis by group.
# Names containing "Unknown" are grey and listed last.
#
#   phylo_species_colours(species) -> list(
#     colours   named colour vector for every species given
#     order     legend order (phylogenetic, then unknowns)
#     distances data.frame(species, distance, distance_group))
# ==============================================================================
PHYLO_TREE_FILE <- "data/raw/external/phylogeny/PhytoPhylo"

phylo_species_colours <- function(species) {
  if (!file.exists(PHYLO_TREE_FILE)) stop("phylogeny not found: ", PHYLO_TREE_FILE, call. = FALSE)
  all_species <- unique(species)
  known   <- all_species[!grepl("Unknown", all_species)]
  unknown <- all_species[grepl("Unknown", all_species)]

  tr <- ape::read.tree(PHYLO_TREE_FILE)
  tips <- intersect(gsub(" ", "_", known), tr$tip.label)
  stopifnot(length(tips) > 1)
  d <- rowMeans(stats::cophenetic(ape::keep.tip(tr, tips)))
  dd <- data.frame(species_clean = names(d), distance = as.numeric(d))
  dd$species <- gsub("_", " ", dd$species_clean)
  dd <- dd[order(dd$distance), ]
  dd$distance_group <- cut(dd$distance, include.lowest = TRUE, labels = FALSE,
                           breaks = stats::quantile(dd$distance, probs = seq(0, 1, length.out = min(8, nrow(dd)))))
  dd$color <- viridisLite::viridis(max(dd$distance_group, na.rm = TRUE))[dd$distance_group]

  colours <- c(stats::setNames(dd$color, dd$species),
               stats::setNames(rep("grey80", length(unknown)), unknown))
  list(colours = colours, order = c(dd$species, unknown),
       distances = dd[, c("species", "distance", "distance_group")])
}
