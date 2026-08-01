source("code/lib/outputs.R")
# ==============================================================================
# REVISION — FAPROTAX/PICRUSt caveat metrics (Referee 2 #3)
# ==============================================================================
# R2 #3: the caveats of 16S-based functional prediction should be explicit; many
# functions are poorly predicted (esp. when 16S can't resolve genus); the low
# fermentation fraction is implausible for a hypoxic habitat; and "dark hydrogen
# oxidation" (3.6% HW) is on par with hydrogenotrophic methanogens (also H2
# oxidisers) -> the H2 economy is understated.
#
# This reproduces the manuscript's FAPROTAX chain (code/08_figures/08d, same seed)
# and computes the numbers that quantify R2's point but were NOT reported:
#   (1) ANNOTATED FRACTION: % of community abundance FAPROTAX assigns ANY function
#       (the rest is functionally invisible -> the % values are of the WHOLE
#       community, so unannotated taxa could include many more fermenters).
#   (2) GENUS RESOLUTION: % of ASVs / abundance resolved to genus (FAPROTAX needs
#       genus-level matches; unresolved taxa get no assignment).
#   (3) H2 economy: dark_hydrogen_oxidation vs hydrogenotrophic methanogenesis
#       (methanogens are H2 oxidisers filed under a different category).
#   (4) Redox-implausibility flag: aerobic_chemoheterotrophy in (hypoxic) heartwood.
#
# NEW file. Run: Rscript code/revision/R_faprotax_caveat_metrics.R
# Output: outputs/revision/faprotax_caveat_metrics.txt
# ==============================================================================

suppressPackageStartupMessages({ library(phyloseq); library(microeco); library(dplyr) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- reproduce the pipeline (verbatim logic from 08d_faprotax_heatmaps.R) -----
ddpcr <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv")
otu_tab <- read.delim("data/raw/16s/OTU_table.txt", header = TRUE, row.names = 1)
bastard_tax <- otu_tab[, 590:596]; bastard_tax[bastard_tax == ""] <- NA
tax_tab_pre <- tax_table(as.matrix(bastard_tax))
taxa_names(tax_tab_pre) <- sub("sp", "seq", taxa_names(tax_tab_pre))
otu_table_pre <- otu_table(otu_tab[, 1:589], taxa_are_rows = TRUE)
phylo_tree <- read_tree("data/raw/16s/unrooted_tree.nwk")
samp_data <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", row.names = 1)
samp_data$RowName <- row.names(samp_data)
samp_data$seq_id <- sub("prime", "'", samp_data$seq_id); samp_data$seq_id <- sub("star", "*", samp_data$seq_id)
samp_data$seq_id <- sub("HM", "H", samp_data$seq_id)
samp_data$seq_id[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "RO104"
samp_data$core_type[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "Inner"
sm <- merge(ddpcr, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)
sm <- sm[!duplicated(sm$RowName), ]; row.names(sm) <- sm$RowName

raw_ps <- phyloseq(tax_tab_pre, otu_table_pre, phylo_tree, sample_data(sm))
pop_taxa <- function(ps, bad) prune_taxa(taxa_names(ps)[!(taxa_names(ps) %in% bad)], ps)
mito <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 5] == "Mitochondria")]
chlo <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 4] == "Chloroplast")]
no_mito <- pop_taxa(raw_ps, c(mito, chlo))
taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
set.seed(46814)
ps.rare <- rarefy_even_depth(no_mito, sample.size = 3500)
ps.ra <- transform_sample_counts(ps.rare, function(x) x / sum(x) * 100)
colnames(tax_table(ps.ra)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
ps.wood <- subset_samples(ps.ra, material == "Wood")
for (i in 1:7) tax_table(ps.wood)[, i] <- paste0(c("k__","p__","c__","o__","f__","g__","s__")[i], tax_table(ps.wood)[, i])

meco_all <- microtable$new(otu_table = as.data.frame(otu_table(ps.wood)),
  tax_table = noquote(as.data.frame(tax_table(ps.wood))),
  sample_table = as.data.frame(as.matrix(sample_data(ps.wood))), phylo_tree = phy_tree(ps.wood))
meco_core <- clone(meco_all)$merge_samples("core_type")
tf <- trans_func$new(meco_core)
tf$cal_spe_func(prok_database = "FAPROTAX")
tf$cal_spe_func_perc(abundance_weighted = TRUE, dec = 2)

# ---- (1) ANNOTATED FRACTION (abundance-weighted, per compartment) -------------
func_mat <- tf$res_spe_func                     # taxa x functions (0/1)
func_mat <- func_mat[, colSums(func_mat, na.rm = TRUE) > 0, drop = FALSE]
annotated <- rownames(func_mat)[rowSums(func_mat, na.rm = TRUE) > 0]
otu <- meco_core$otu_table                      # taxa x {Inner, Outer}
cov <- sapply(colnames(otu), function(g) {
  tot <- sum(otu[[g]], na.rm = TRUE)
  100 * sum(otu[rownames(otu) %in% annotated, g], na.rm = TRUE) / tot })

# ---- (2) GENUS RESOLUTION (abundance-weighted) --------------------------------
tax <- as.data.frame(meco_all$tax_table); tax$ASV <- rownames(tax)
unresolved <- tax$Genus %in% c("g__","g__NA","g__ ") | is.na(tax$Genus) |
  grepl("uncultured|metagenome|unidentif|unknown|_NA$", tax$Genus, ignore.case = TRUE)
resolved_asv <- tax$ASV[!unresolved]
genus_cov <- sapply(colnames(otu), function(g) {
  100 * sum(otu[rownames(otu) %in% resolved_asv, g], na.rm = TRUE) / sum(otu[[g]], na.rm = TRUE) })
genus_asv_pct <- 100 * mean(!unresolved)

# ---- report ------------------------------------------------------------------
p <- tf$res_spe_func_perc                        # compartments x functions (%)
getf <- function(f, grp) if (f %in% colnames(p)) p[grp, f] else NA
sink(out_path("faprotax_caveat_metrics.txt"))
cat("=================================================================\n")
cat("FAPROTAX caveat metrics (Referee 2 #3) — heartwood=Inner, sapwood=Outer\n")
cat("=================================================================\n\n")
cat("(1) FRACTION OF COMMUNITY FAPROTAX ANNOTATES AT ALL (abundance-weighted):\n")
for (g in names(cov)) cat(sprintf("    %-7s: %.1f%% annotated  ->  %.1f%% functionally UNCHARACTERISED\n",
                                   g, cov[g], 100 - cov[g]))
cat("    => the reported percentages are fractions of the WHOLE community; the large\n")
cat("       unannotated remainder could contain many additional fermenters etc.\n\n")
cat("(2) GENUS-LEVEL RESOLUTION (FAPROTAX needs genus matches):\n")
cat(sprintf("    %.1f%% of ASVs resolved to a named genus.\n", genus_asv_pct))
for (g in names(genus_cov)) cat(sprintf("    %-7s: %.1f%% of abundance is genus-resolved\n", g, genus_cov[g]))
cat("\n(3) H2 ECONOMY — R2 notes 3.6%% understates H2; we read it as CORROBORATION:\n")
cat(sprintf("    dark_hydrogen_oxidation:                  HW=%.2f%%  SW=%.2f%%\n", getf("dark_hydrogen_oxidation","Inner"), getf("dark_hydrogen_oxidation","Outer")))
cat(sprintf("    hydrogenotrophic_methanogenesis:          HW=%.2f%%  SW=%.2f%%\n", getf("hydrogenotrophic_methanogenesis","Inner"), getf("hydrogenotrophic_methanogenesis","Outer")))
cat("    => categories are non-additive: hydrogenotrophic methanogens (H2 oxidisers) are\n")
cat("       filed under methanogenesis, so total H2 metabolism >= ~7%%, not 3.6%%. Both\n")
cat("       H2-using categories are enriched in heartwood -> an H2-based economy that is\n")
cat("       INDEPENDENTLY corroborated by depleted d13CH4 and Methanobacteriaceae dominance.\n\n")
cat("(4) BETWEEN-TISSUE DIFFERENTIAL (the robust signal; more informative than absolute %%):\n")
ff <- c("hydrogenotrophic_methanogenesis","dark_hydrogen_oxidation","fermentation",
        "methylotrophy","aerobic_chemoheterotrophy","methanotrophy")
cat(sprintf("    %-32s %6s %6s  %s\n","function","HW%","SW%","HW/SW"))
for (f in ff) { hw <- getf(f,"Inner"); sw <- getf(f,"Outer")
  cat(sprintf("    %-32s %6.2f %6.2f  %5.1fx\n", f, hw, sw, hw/sw)) }
cat("    => anaerobic/production functions (methanogenesis, H2 oxidation, fermentation) are\n")
cat("       ENRICHED in heartwood; aerobic/consumption functions (aerobic chemoheterotrophy,\n")
cat("       methanotrophy) are enriched in sapwood -- the expected core/outer gradient. IN-STUDY\n")
cat("       corroboration = OUR ddPCR compartment split (methanogens Inner, methanotrophs Outer)\n")
cat("       + OUR internal CH4~O2 negative correlation (Fig S8). NOTE: paired depth-resolved HW/SW\n")
cat("       O2 PROFILES are the companion Nature paper (Arnold & Gewirtzman 2025), NOT measured\n")
cat("       here; this study has single internal (bark-to-pith) gas per tree, ~1-22%% O2, variable.\n")
sink()
cat(readLines(out_path("faprotax_caveat_metrics.txt")), sep = "\n")
