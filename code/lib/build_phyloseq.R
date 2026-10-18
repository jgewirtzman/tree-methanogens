# ==============================================================================
# build_phyloseq.R -- the single definition of the 16S phyloseq object, the
# organelle rule and the rarefaction (depth 3,500, seed 46814)
# ------------------------------------------------------------------------------
# Built from the canonical tables (see SOURCES below).
#
# Usage:
#   source("code/lib/build_phyloseq.R")
#   ps16     <- build_16s_phyloseq()
#   no_mito  <- ps16$ps          # unrarefied, mito/chloroplast removed
#   ps.ra    <- rarefy_16s(no_mito)
# Scripts that also need the intermediates take them from the same list:
#   phylo_tree <- ps16$tree      # the ASV phylogeny
#   ddpcr      <- ps16$ddpcr     # manuscript_statistics
# ==============================================================================

# Canonical rarefaction parameters. These were identical in all nine copies; they
# live here so a change is a one-line change rather than a nine-file sweep.
RAREFY_SEED  <- 46814
RAREFY_DEPTH <- 3500

# Organelle ASVs: SILVA family Mitochondria or order Chloroplast. The one rule for
# every 16S and PICRUSt2 analysis: these reads are removed before any relative
# abundance is formed (build_16s_phyloseq(), asv_reads() in tables_read.R, and the
# PICRUSt2 pathway table, code/01_import/06_picrust_pathways.R). tax: a data.frame
# (or matrix) with columns Family and Order.
is_organelle <- function(tax) {
  tax <- as.data.frame(tax, stringsAsFactors = FALSE)
  tax$Family %in% "Mitochondria" | tax$Order %in% "Chloroplast"
}

pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- phyloseq::taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  phyloseq::prune_taxa(allTaxa, physeq)
}

# SOURCES: the canonical tables only, plus the 16S tree.
#   asv_counts_16S, asv_taxonomy_16S   counts and SILVA ranks (main run)
#   samples_16S                        library -> sample_id (main run)
#   samples, trees                     sample -> tree, position, material; tree -> species
#   ddpcr                              *_loose / *_strict: mean over the sample's qc_pass
#                                      runs of copies per uL reaction x DDPCR_WELL_UL
#                                      (copies per 20 uL well, the unit the old raw
#                                      metadata columns held)
#   qpcr                               X16S_per_ul: facility 16S copies per uL extract
#   data/raw/16s/unrooted_tree.nwk     the ASV phylogeny (no table holds it)
# sample_data keeps the column names the analysis scripts use: seq_id (the tree_id;
# "MISC" for libraries with no survey sample), core_type (Inner, Outer, Mineral,
# Organic), species.x (trees.csv species_code), material (Wood, Soil), the ddPCR
# columns and X16S_per_ul; plus library_id, sample_id, tree_id, library_type.
build_16s_phyloseq <- function() {
  stopifnot(requireNamespace("phyloseq", quietly = TRUE))
  source("code/lib/ddpcr_constants.R", local = TRUE)      # DDPCR_WELL_UL
  rt <- function(t, ...) utils::read.csv(file.path("data/tables", paste0(t, ".csv")), stringsAsFactors = FALSE,
                                         check.names = FALSE, na.strings = "", ...)
  cnt <- rt("asv_counts_16S"); tax <- rt("asv_taxonomy_16S")
  stopifnot(identical(cnt$feature_id, tax$feature_id))
  counts <- as.matrix(cnt[, -1]); rownames(counts) <- cnt$feature_id
  taxm <- as.matrix(tax[, -1]); rownames(taxm) <- tax$feature_id
  colnames(taxm) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  L <- rt("samples_16S"); L <- L[L$run == "main", ]
  stopifnot(setequal(L$library_id, colnames(counts)), !anyDuplicated(L$library_id))
  L <- L[match(colnames(counts), L$library_id), ]
  S <- rt("samples"); TR <- rt("trees"); DD <- rt("ddpcr"); Q <- rt("qpcr")
  sm <- S[match(L$sample_id, S$sample_id), ]
  POS <- c(inner = "Inner", outer = "Outer", mineral = "Mineral", organic = "Organic")
  MAT <- c(wood = "Wood", soil = "Soil")
  sd <- data.frame(library_id = L$library_id, sample_id = L$sample_id, tree_id = sm$tree_id,
                   library_type = L$library_type,
                   seq_id = ifelse(is.na(sm$tree_id), "MISC", sm$tree_id),
                   core_type = unname(POS[sm$position]), species.x = TR$species_code[match(sm$tree_id, TR$tree_id)],
                   material = unname(MAT[sm$material]), stringsAsFactors = FALSE)
  DD <- DD[DD$qc_pass %in% TRUE, ]
  key <- c(mcra_probe = "mcrA|probe", mcra_eva = "mcrA|EvaGreen", pmoa = "pmoA|EvaGreen", mmox = "mmoX|EvaGreen",
           bact = "16S_bacteria|EvaGreen", arc = "16S_archaea|EvaGreen")
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  for (k in names(key)) {
    x <- DD[paste(DD$gene, DD$assay, sep = "|") == key[[k]], ]
    sd[[paste0(k, "_loose")]]  <- DDPCR_WELL_UL * vapply(sd$sample_id, function(i) mean_na(x$copies_per_uL_rxn[x$sample_id %in% i]), 0)
    sd[[paste0(k, "_strict")]] <- DDPCR_WELL_UL * vapply(sd$sample_id, function(i) mean_na(x$copies_per_uL_rxn_strict[x$sample_id %in% i]), 0)
  }
  q16 <- Q[Q$target == "16S", ]
  sd$X16S_per_ul <- q16$copies_per_uL_extract[match(sd$sample_id, q16$sample_id)]
  rownames(sd) <- sd$library_id

  phylo_tree <- phyloseq::read_tree("data/raw/16s/unrooted_tree.nwk")
  raw_ps <- phyloseq::phyloseq(phyloseq::tax_table(taxm), phyloseq::otu_table(counts, taxa_are_rows = TRUE), phylo_tree,
                               phyloseq::sample_data(sd))

  tt <- as.data.frame(as(phyloseq::tax_table(raw_ps), "matrix"), stringsAsFactors = FALSE)
  mitochondria <- rownames(tt)[tt$Family %in% "Mitochondria"]
  chloroplast  <- rownames(tt)[tt$Order %in% "Chloroplast"]
  no_mito <- pop_taxa(raw_ps, rownames(tt)[is_organelle(tt)])

  list(ps = no_mito, raw = raw_ps, tree = phylo_tree, ddpcr = sd[!is.na(sd$sample_id), ],
       samp_data = sd, n_mitochondria = length(mitochondria), n_chloroplast = length(chloroplast))
}

# Rarefaction is stochastic: the seed is what makes every figure reproducible and
# mutually consistent. Do not call rarefy_even_depth() directly.
rarefy_16s <- function(ps, depth = RAREFY_DEPTH, seed = RAREFY_SEED) {
  set.seed(seed)
  phyloseq::rarefy_even_depth(ps, sample.size = depth)
}
