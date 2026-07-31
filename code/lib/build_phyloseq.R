# ==============================================================================
# Single definition of the 16S phyloseq object.
#
# WHY THIS EXISTS. Nine scripts each rebuilt this object from the same five raw
# files with their own ~45-line copy of the assembly: the same hardcoded OTU
# column split (1:589 counts, 590:596 taxonomy), the same seq_id patches, the
# same single-fastq override, the same duplicate-RowName drop, the same
# mitochondria/chloroplast removal. Before consolidating, all nine preambles were
# executed in isolation and their objects fingerprinted: every one produced
# 71,511 taxa / 587 samples / 6,427,747 counts with identical taxa and sample
# names -- 1 distinct object out of 9. They had not drifted. The point of moving
# them here is that they cannot start.
#
# The column split is positional, so a differently-shaped OTU_table.txt would
# silently mis-split it. Asserted below rather than assumed.
#
# Usage:
#   source("code/lib/build_phyloseq.R")
#   ps16     <- build_16s_phyloseq()
#   no_mito  <- ps16$ps          # unrarefied, mito/chloroplast removed
#   ps.ra    <- rarefy_16s(no_mito)
# Scripts that also need the intermediates take them from the same list:
#   phylo_tree <- ps16$tree      # 08d, rev_faprotax_dump_HW_SW, manuscript_statistics
#   ddpcr      <- ps16$ddpcr     # manuscript_statistics
# ==============================================================================

# Canonical rarefaction parameters. These were identical in all nine copies; they
# live here so a change is a one-line change rather than a nine-file sweep.
RAREFY_SEED  <- 46814
RAREFY_DEPTH <- 3500

# OTU_table.txt carries counts and taxonomy in one file, split by position.
OTU_COUNT_COLS <- 1:589
OTU_TAX_COLS   <- 590:596

pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- phyloseq::taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  phyloseq::prune_taxa(allTaxa, physeq)
}

build_16s_phyloseq <- function() {
  stopifnot(requireNamespace("phyloseq", quietly = TRUE))

  ddpcr    <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv")
  water    <- read.csv("data/raw/tree_cores/Tree_Core_Sectioning_Data.csv")
  qpcr_16s <- read.csv("data/raw/16s/16s_w_metadata.csv")
  qpcr_16s <- subset(qpcr_16s, qpcr_16s$Sample.ID != "None")
  qpcr_16s <- qpcr_16s[, c(3, 4, 6)]

  water$seq_id <- toupper(water$seq_id)
  ddpcr <- merge(ddpcr, water,    by = c("core_type", "seq_id"), all.x = TRUE)
  ddpcr <- merge(ddpcr, qpcr_16s, by = c("core_type", "seq_id"), all.x = TRUE)

  otu_tab <- read.delim("data/raw/16s/OTU_table.txt", header = TRUE, row.names = 1)
  # Positional split -- fail loudly if the table is no longer 596 columns wide.
  if (ncol(otu_tab) != max(OTU_TAX_COLS))
    stop(sprintf("OTU_table.txt has %d columns, expected %d; the positional split at %s would mis-assign taxonomy.",
                 ncol(otu_tab), max(OTU_TAX_COLS), deparse(OTU_TAX_COLS)), call. = FALSE)

  bastard_tax <- otu_tab[, OTU_TAX_COLS]
  bastard_tax[bastard_tax == ""] <- NA
  tax_tab_pre <- phyloseq::tax_table(bastard_tax)
  phyloseq::taxa_names(tax_tab_pre) <- sub("sp", "seq", phyloseq::taxa_names(tax_tab_pre))

  otu_tab_corr  <- otu_tab[, OTU_COUNT_COLS]
  otu_table_pre <- phyloseq::otu_table(otu_tab_corr, taxa_are_rows = TRUE)

  phylo_tree <- phyloseq::read_tree("data/raw/16s/unrooted_tree.nwk")
  samp_data  <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", row.names = 1)
  samp_data$RowName <- row.names(samp_data)

  # seq_id patches: filesystem-safe substitutions reversed, plus one sample whose
  # mapping file entry is wrong at source.
  samp_data$seq_id <- sub("prime", "'", samp_data$seq_id)
  samp_data$seq_id <- sub("star",  "*", samp_data$seq_id)
  samp_data$seq_id <- sub("HM",    "H", samp_data$seq_id)
  samp_data$seq_id[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"]    <- "RO104"
  samp_data$core_type[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "Inner"

  samp_data_merged <- merge(ddpcr, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)

  dups <- which(duplicated(samp_data_merged$RowName))
  if (length(dups)) samp_data_merged <- samp_data_merged[-dups, ]
  row.names(samp_data_merged) <- samp_data_merged$RowName

  raw_ps <- phyloseq::phyloseq(tax_tab_pre, otu_table_pre, phylo_tree,
                               phyloseq::sample_data(samp_data_merged))

  tt <- phyloseq::tax_table(raw_ps)
  mitochondria <- rownames(tt)[which(tt[, 5] == "Mitochondria")]
  chloroplast  <- rownames(tt)[which(tt[, 4] == "Chloroplast")]
  no_mito <- pop_taxa(raw_ps, c(mitochondria, chloroplast))

  list(ps = no_mito, raw = raw_ps, tree = phylo_tree, ddpcr = ddpcr,
       samp_data = samp_data_merged,
       n_mitochondria = length(mitochondria), n_chloroplast = length(chloroplast))
}

# Rarefaction is stochastic: the seed is what makes every figure reproducible and
# mutually consistent. Do not call rarefy_even_depth() directly.
rarefy_16s <- function(ps, depth = RAREFY_DEPTH, seed = RAREFY_SEED) {
  set.seed(seed)
  phyloseq::rarefy_even_depth(ps, sample.size = depth)
}
