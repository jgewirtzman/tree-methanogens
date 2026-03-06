################################################################################
# 05_soil_microbiome_exploration.R
#
# Quick exploratory analysis: Is there enough signal in the soil microbiome
# (bacterial 16S + fungal ITS) for a standalone paper on tree species and
# soil moisture/redox effects?
#
# Key questions:
#   1. Do tree species, VWC, ORP, and/or soil horizon significantly structure
#      soil bacterial and fungal communities? (PERMANOVA)
#   2. How much variance does each factor explain?
#   3. Is there a species x moisture interaction?
#   4. After conditioning on sampling location (plot), do species effects hold?
#   5. What does the ordination look like?
################################################################################

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# ── Paths ────────────────────────────────────────────────────────────────────
base_dir <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane"
methanogens_dir <- file.path(base_dir, "tree-methanogens")
microbiome_dir  <- file.path(base_dir, "tree-microbiome")

out_dir <- file.path(methanogens_dir, "outputs/figures/new_paper")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# 0. BUILD PLOT (LOCATION) LOOKUP from merged tree dataset
################################################################################

cat("\n========== BUILDING PLOT LOOKUP ==========\n")

merged_trees <- read.csv(file.path(methanogens_dir,
                                   "data/processed/integrated/merged_tree_dataset_final.csv"),
                         stringsAsFactors = FALSE)

# tree_id in merged is lowercase; seq_id in 16S/ITS is uppercase → uppercase both
plot_lookup <- merged_trees %>%
  mutate(seq_id_upper = toupper(trimws(tree_id))) %>%
  filter(plot != "" & !is.na(plot) & plot != "NA") %>%
  # Rename to 'site' to avoid collision with base::plot()
  select(seq_id_upper, site = plot) %>%
  distinct(seq_id_upper, .keep_all = TRUE)

cat("Plot lookup:", nrow(plot_lookup), "trees with site assignment\n")
cat("Site breakdown:", paste(names(table(plot_lookup$site)),
                             table(plot_lookup$site), sep = "=", collapse = ", "), "\n")

################################################################################
# 1. BUILD 16S PHYLOSEQ (soil only)
################################################################################

cat("\n========== BUILDING 16S PHYLOSEQ ==========\n")

# Load OTU table (rows = ASVs, cols = samples + taxonomy)
otu_16s <- read.delim(file.path(methanogens_dir, "data/raw/16s/OTU_table.txt"),
                      header = TRUE, row.names = 1)

# Split taxonomy from OTU counts
tax_cols <- which(colnames(otu_16s) %in%
                    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
sample_cols <- setdiff(seq_len(ncol(otu_16s)), tax_cols)

tax_16s <- otu_16s[, tax_cols]
tax_16s[tax_16s == ""] <- NA
otu_counts_16s <- otu_16s[, sample_cols]

# Load mapping file + metadata
samp_16s <- read.delim(
  file.path(methanogens_dir, "data/raw/16s/tree_16s_mapping_dada2_corrected.txt"),
  row.names = 1
)
samp_16s$RowName <- rownames(samp_16s)

# Fix known ID issues (from existing pipeline)
samp_16s$seq_id <- sub("prime", "'", samp_16s$seq_id)
samp_16s$seq_id <- sub("star", "*", samp_16s$seq_id)
samp_16s$seq_id <- sub("HM", "H", samp_16s$seq_id)
samp_16s$seq_id[samp_16s$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "RO104"
samp_16s$core_type[samp_16s$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "Inner"

# Merge ddPCR + environmental metadata
ddpcr <- read.csv(file.path(methanogens_dir, "data/raw/ddpcr/ddPCR_meta_all_data.csv"))
water <- read.csv(file.path(methanogens_dir, "data/raw/tree_cores/Tree_Core_Sectioning_Data.csv"))
water$seq_id <- toupper(water$seq_id)
ddpcr <- merge(ddpcr, water, by = c("core_type", "seq_id"), all.x = TRUE)

samp_16s_merged <- merge(ddpcr, samp_16s, by = c("seq_id", "core_type"), all.y = TRUE)
dups <- duplicated(samp_16s_merged$RowName)
samp_16s_merged <- samp_16s_merged[!dups, ]
rownames(samp_16s_merged) <- samp_16s_merged$RowName

# Join plot assignment via seq_id
samp_16s_merged$seq_id_upper <- toupper(trimws(samp_16s_merged$seq_id))
samp_16s_merged <- merge(samp_16s_merged, plot_lookup,
                         by = "seq_id_upper", all.x = TRUE)
rownames(samp_16s_merged) <- samp_16s_merged$RowName

# Build phyloseq
ps_16s_raw <- phyloseq(
  tax_table(as.matrix(tax_16s)),
  otu_table(as.matrix(otu_counts_16s), taxa_are_rows = TRUE),
  sample_data(samp_16s_merged)
)

# Remove contaminants
mito <- taxa_names(ps_16s_raw)[which(tax_table(ps_16s_raw)[, "Family"] == "Mitochondria")]
chloro <- taxa_names(ps_16s_raw)[which(tax_table(ps_16s_raw)[, "Order"] == "Chloroplast")]
ps_16s_clean <- prune_taxa(!(taxa_names(ps_16s_raw) %in% c(mito, chloro)), ps_16s_raw)

# Rarefy
taxa_names(ps_16s_clean) <- paste0("ASV", seq(ntaxa(ps_16s_clean)))
set.seed(46814)
ps_16s_rare <- rarefy_even_depth(ps_16s_clean, sample.size = 3500, verbose = FALSE)
colnames(tax_table(ps_16s_rare)) <- c("Kingdom", "Phylum", "Class", "Order",
                                       "Family", "Genus", "Species")

# Filter to soil, well-sampled species
well_sampled <- c("ACRU", "ACSA", "BEAL", "BELE", "BEPA",
                  "FAGR", "FRAM", "PIST", "QURU", "TSCA")

sd_16s <- data.frame(sample_data(ps_16s_rare), stringsAsFactors = FALSE)
sd_16s$horizon <- case_when(
  sd_16s$core_type == "Mineral" ~ "Mineral",
  sd_16s$core_type == "Organic" ~ "Organic",
  TRUE ~ NA_character_
)
# Write horizon back to phyloseq
sample_data(ps_16s_rare) <- sample_data(sd_16s)

soil_16s_idx <- which(sd_16s$core_type %in% c("Mineral", "Organic") &
                        sd_16s$species.x %in% well_sampled)
ps_16s_soil <- prune_samples(rownames(sd_16s)[soil_16s_idx], ps_16s_rare)
ps_16s_soil <- prune_taxa(taxa_sums(ps_16s_soil) > 0, ps_16s_soil)

cat("16S soil samples:", nsamples(ps_16s_soil), "\n")
cat("16S soil ASVs:", ntaxa(ps_16s_soil), "\n")
cat("Species represented:\n")
print(table(sample_data(ps_16s_soil)$species.x))

################################################################################
# 2. BUILD ITS PHYLOSEQ (soil only)
################################################################################

cat("\n========== BUILDING ITS PHYLOSEQ ==========\n")

otu_its <- read.delim(file.path(microbiome_dir, "ITS/OTU_table.txt"),
                      header = TRUE, row.names = 1)

tax_cols_its <- which(colnames(otu_its) %in%
                        c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
sample_cols_its <- setdiff(seq_len(ncol(otu_its)), tax_cols_its)

tax_its <- otu_its[, tax_cols_its]
tax_its[tax_its == ""] <- NA
otu_counts_its <- otu_its[, sample_cols_its]

# Use annotated metadata (already has core_type, species, VWC, ORP merged)
samp_its <- read.csv(
  file.path(microbiome_dir, "Other-Metadata/annotated_metadata/ITS_tree_sample_table_with_meta.csv"),
  row.names = 1, check.names = FALSE
)
samp_its$RowName <- rownames(samp_its)

# Join plot assignment via seq_id
samp_its$seq_id_upper <- toupper(trimws(samp_its$seq_id))
samp_its <- merge(samp_its, plot_lookup, by = "seq_id_upper", all.x = TRUE)
rownames(samp_its) <- samp_its$RowName

# Keep only samples present in OTU table
common_its <- intersect(rownames(samp_its), colnames(otu_counts_its))
cat("ITS samples in both OTU table and metadata:", length(common_its), "\n")
samp_its <- samp_its[common_its, ]
otu_counts_its <- otu_counts_its[, common_its]

# Build phyloseq
ps_its_raw <- phyloseq(
  tax_table(as.matrix(tax_its)),
  otu_table(as.matrix(otu_counts_its), taxa_are_rows = TRUE),
  sample_data(samp_its)
)

# Rarefy ITS
taxa_names(ps_its_raw) <- paste0("OTU", seq(ntaxa(ps_its_raw)))
set.seed(46814)
ps_its_rare <- rarefy_even_depth(ps_its_raw, sample.size = 3500, verbose = FALSE)
colnames(tax_table(ps_its_rare)) <- c("Kingdom", "Phylum", "Class", "Order",
                                       "Family", "Genus", "Species")

# Filter to soil, well-sampled species
sd_its <- data.frame(sample_data(ps_its_rare), stringsAsFactors = FALSE)
sd_its$horizon <- case_when(
  sd_its$core_type == "Mineral" ~ "Mineral",
  sd_its$core_type == "Organic" ~ "Organic",
  TRUE ~ NA_character_
)
# Write horizon back to phyloseq
sample_data(ps_its_rare) <- sample_data(sd_its)

soil_its_idx <- which(sd_its$core_type %in% c("Mineral", "Organic") &
                        sd_its$species.x %in% well_sampled)
ps_its_soil <- prune_samples(rownames(sd_its)[soil_its_idx], ps_its_rare)
ps_its_soil <- prune_taxa(taxa_sums(ps_its_soil) > 0, ps_its_soil)

cat("ITS soil samples:", nsamples(ps_its_soil), "\n")
cat("ITS soil ASVs:", ntaxa(ps_its_soil), "\n")
cat("Species represented:\n")
print(table(sample_data(ps_its_soil)$species.x))

################################################################################
# 2b. CONFOUNDING DIAGNOSTIC: Species × Plot
################################################################################

cat("\n========== SPECIES × PLOT CONFOUNDING ==========\n")

sd16_diag <- data.frame(sample_data(ps_16s_soil), stringsAsFactors = FALSE)
cat("\n--- 16S Soil: Species × Plot ---\n")
print(table(sd16_diag$species.x, sd16_diag$site, useNA = "ifany"))

cat("\nPlot assignment coverage:", sum(!is.na(sd16_diag$site)), "of",
    nrow(sd16_diag), "samples\n")

# Identify species present in ≥2 plots (for sensitivity analysis)
sp_plot_tab <- table(sd16_diag$species.x, sd16_diag$site)
sp_plot_tab <- sp_plot_tab[, colnames(sp_plot_tab) != "NA", drop = FALSE]
plots_per_sp <- apply(sp_plot_tab > 0, 1, sum)
multi_plot_species <- names(plots_per_sp[plots_per_sp >= 2])

cat("\nSpecies in ≥2 plots (usable for plot-controlled analysis):\n")
for (sp in multi_plot_species) {
  cat(sprintf("  %s: %s\n", sp,
              paste(colnames(sp_plot_tab)[sp_plot_tab[sp, ] > 0],
                    sp_plot_tab[sp, sp_plot_tab[sp, ] > 0],
                    sep = "=", collapse = ", ")))
}

cat("\nSpecies confined to 1 plot (confounded):\n")
single_plot_species <- names(plots_per_sp[plots_per_sp <= 1])
for (sp in single_plot_species) {
  nonzero <- sp_plot_tab[sp, sp_plot_tab[sp, ] > 0, drop = FALSE]
  cat(sprintf("  %s: %s\n", sp,
              paste(colnames(nonzero), nonzero, sep = "=", collapse = ", ")))
}

################################################################################
# 3. PERMANOVA — 16S Bacteria (with plot)
################################################################################

cat("\n========== PERMANOVA: 16S BACTERIA ==========\n")

sd16 <- data.frame(sample_data(ps_16s_soil), stringsAsFactors = FALSE)
sd16$VWC <- as.numeric(sd16$mean_vwc)
sd16$ORP <- as.numeric(sd16$mean_orp)

# Subset to samples with complete environmental data AND plot
complete_16s <- !is.na(sd16$VWC) & !is.na(sd16$ORP) & !is.na(sd16$site)
cat("16S samples with complete data (env + plot):", sum(complete_16s), "of", nrow(sd16), "\n")

ps_16s_complete <- prune_samples(rownames(sd16)[complete_16s], ps_16s_soil)
ps_16s_complete <- prune_taxa(taxa_sums(ps_16s_complete) > 0, ps_16s_complete)
sd16c <- data.frame(sample_data(ps_16s_complete), stringsAsFactors = FALSE)
sd16c$VWC <- as.numeric(sd16c$mean_vwc)
sd16c$ORP <- as.numeric(sd16c$mean_orp)

bc_16s <- phyloseq::distance(ps_16s_complete, method = "bray")

# A) Full model WITHOUT plot conditioning (baseline)
set.seed(999)
perm_16s_noplot <- adonis2(
  bc_16s ~ species.x + VWC + ORP + horizon,
  data = sd16c, permutations = 999, by = "margin"
)
cat("\n--- 16S PERMANOVA: NO plot conditioning (marginal) ---\n")
print(perm_16s_noplot)

# B) Model WITH plot conditioning (partial out spatial effect)
set.seed(999)
perm_16s_cond_plot <- adonis2(
  bc_16s ~ species.x + VWC + ORP + horizon,
  data = sd16c, permutations = 999, by = "margin",
  strata = sd16c$site
)
cat("\n--- 16S PERMANOVA: Stratified by SITE (marginal) ---\n")
print(perm_16s_cond_plot)

# C) Sequential with plot first, then species
set.seed(999)
perm_16s_seq <- adonis2(
  bc_16s ~ site + species.x + VWC + ORP + horizon,
  data = sd16c, permutations = 999, by = "terms"
)
cat("\n--- 16S PERMANOVA: Sequential (plot first) ---\n")
print(perm_16s_seq)

# D) Interaction model
set.seed(999)
perm_16s_interact <- adonis2(
  bc_16s ~ species.x + VWC + ORP + horizon + species.x:VWC,
  data = sd16c, permutations = 999, by = "margin",
  strata = sd16c$site
)
cat("\n--- 16S PERMANOVA: Conditioned on plot + interaction ---\n")
print(perm_16s_interact)

################################################################################
# 4. PERMANOVA — ITS Fungi (with plot)
################################################################################

cat("\n========== PERMANOVA: ITS FUNGI ==========\n")

sd_it <- data.frame(sample_data(ps_its_soil), stringsAsFactors = FALSE)
sd_it$VWC <- as.numeric(sd_it$mean_vwc)
sd_it$ORP <- as.numeric(sd_it$mean_orp)

complete_its <- !is.na(sd_it$VWC) & !is.na(sd_it$ORP) & !is.na(sd_it$site)
cat("ITS samples with complete data (env + plot):", sum(complete_its), "of", nrow(sd_it), "\n")

ps_its_complete <- prune_samples(rownames(sd_it)[complete_its], ps_its_soil)
ps_its_complete <- prune_taxa(taxa_sums(ps_its_complete) > 0, ps_its_complete)
sditc <- data.frame(sample_data(ps_its_complete), stringsAsFactors = FALSE)
sditc$VWC <- as.numeric(sditc$mean_vwc)
sditc$ORP <- as.numeric(sditc$mean_orp)

bc_its <- phyloseq::distance(ps_its_complete, method = "bray")

# A) No plot conditioning
set.seed(999)
perm_its_noplot <- adonis2(
  bc_its ~ species.x + VWC + ORP + horizon,
  data = sditc, permutations = 999, by = "margin"
)
cat("\n--- ITS PERMANOVA: NO plot conditioning (marginal) ---\n")
print(perm_its_noplot)

# B) Conditioned on plot
set.seed(999)
perm_its_cond_plot <- adonis2(
  bc_its ~ species.x + VWC + ORP + horizon,
  data = sditc, permutations = 999, by = "margin",
  strata = sditc$site
)
cat("\n--- ITS PERMANOVA: Stratified by SITE (marginal) ---\n")
print(perm_its_cond_plot)

# C) Sequential with site first
set.seed(999)
perm_its_seq <- adonis2(
  bc_its ~ site + species.x + VWC + ORP + horizon,
  data = sditc, permutations = 999, by = "terms"
)
cat("\n--- ITS PERMANOVA: Sequential (site first) ---\n")
print(perm_its_seq)

# D) Interaction model
set.seed(999)
perm_its_interact <- adonis2(
  bc_its ~ species.x + VWC + ORP + horizon + species.x:VWC,
  data = sditc, permutations = 999, by = "margin",
  strata = sditc$site
)
cat("\n--- ITS PERMANOVA: Stratified by site + interaction ---\n")
print(perm_its_interact)

################################################################################
# 5. SENSITIVITY: Multi-plot species subset
################################################################################

cat("\n========== SENSITIVITY: MULTI-PLOT SPECIES ==========\n")

# Only species found in ≥2 plots
cat("Multi-plot species:", paste(multi_plot_species, collapse = ", "), "\n")

# 16S subset
ps_16s_multi <- prune_samples(sd16c$species.x %in% multi_plot_species, ps_16s_complete)
ps_16s_multi <- prune_taxa(taxa_sums(ps_16s_multi) > 0, ps_16s_multi)
sd16m <- data.frame(sample_data(ps_16s_multi), stringsAsFactors = FALSE)
sd16m$VWC <- as.numeric(sd16m$mean_vwc)
sd16m$ORP <- as.numeric(sd16m$mean_orp)
bc_16s_multi <- phyloseq::distance(ps_16s_multi, method = "bray")

cat("\n16S multi-plot subset:", nsamples(ps_16s_multi), "samples,",
    length(unique(sd16m$species.x)), "species\n")
print(table(sd16m$species.x, sd16m$site))

set.seed(999)
perm_16s_multi <- adonis2(
  bc_16s_multi ~ species.x + VWC + ORP + horizon,
  data = sd16m, permutations = 999, by = "margin",
  strata = sd16m$site
)
cat("\n--- 16S MULTI-PLOT: Conditioned on plot (marginal) ---\n")
print(perm_16s_multi)

# ITS subset
ps_its_multi <- prune_samples(sditc$species.x %in% multi_plot_species, ps_its_complete)
ps_its_multi <- prune_taxa(taxa_sums(ps_its_multi) > 0, ps_its_multi)
sditm <- data.frame(sample_data(ps_its_multi), stringsAsFactors = FALSE)
sditm$VWC <- as.numeric(sditm$mean_vwc)
sditm$ORP <- as.numeric(sditm$mean_orp)
bc_its_multi <- phyloseq::distance(ps_its_multi, method = "bray")

cat("\nITS multi-plot subset:", nsamples(ps_its_multi), "samples,",
    length(unique(sditm$species.x)), "species\n")
print(table(sditm$species.x, sditm$site))

set.seed(999)
perm_its_multi <- adonis2(
  bc_its_multi ~ species.x + VWC + ORP + horizon,
  data = sditm, permutations = 999, by = "margin",
  strata = sditm$site
)
cat("\n--- ITS MULTI-PLOT: Conditioned on plot (marginal) ---\n")
print(perm_its_multi)

################################################################################
# 6. ORDINATION — PCoA with Bray-Curtis
################################################################################

cat("\n========== ORDINATION ==========\n")

# 16S PCoA
ord_16s <- ordinate(ps_16s_complete, method = "PCoA", distance = bc_16s)
eig_16s <- ord_16s$values$Relative_eig[1:2] * 100

# ITS PCoA
ord_its <- ordinate(ps_its_complete, method = "PCoA", distance = bc_its)
eig_its <- ord_its$values$Relative_eig[1:2] * 100

# Mycorrhizal type assignment
myco_map <- c(
  ACRU = "AM/ECM", ACSA = "AM", BEAL = "ECM", BELE = "ECM", BEPA = "ECM",
  FAGR = "ECM", FRAM = "AM", PIST = "ECM", QURU = "ECM", TSCA = "ECM"
)

# Build plot dataframes
make_ord_df <- function(ord, ps_complete) {
  sc <- data.frame(ord$vectors[, 1:2])
  colnames(sc) <- c("Axis1", "Axis2")
  sd <- data.frame(sample_data(ps_complete), stringsAsFactors = FALSE)
  sc$species <- sd$species.x
  sc$VWC <- as.numeric(sd$mean_vwc)
  sc$ORP <- as.numeric(sd$mean_orp)
  sc$horizon <- ifelse(sd$core_type == "Mineral", "Mineral", "Organic")
  sc$mycorrhiza <- myco_map[sc$species]
  sc$site <- sd$site
  sc
}

df_16s <- make_ord_df(ord_16s, ps_16s_complete)
df_its <- make_ord_df(ord_its, ps_its_complete)

# ── Plot: colored by VWC ──
theme_ord <- theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        legend.position = "right")

p1 <- ggplot(df_16s, aes(Axis1, Axis2, color = VWC)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "viridis") +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_16s[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_16s[2]),
       title = "16S Bacteria — Soil VWC",
       color = "VWC (%)", shape = "Horizon") +
  theme_ord

p2 <- ggplot(df_its, aes(Axis1, Axis2, color = VWC)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "viridis") +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_its[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_its[2]),
       title = "ITS Fungi — Soil VWC",
       color = "VWC (%)", shape = "Horizon") +
  theme_ord

# ── Plot: colored by species ──
sp_colors <- c(
  ACRU = "#E41A1C", ACSA = "#FF7F00", BEAL = "#4DAF4A", BELE = "#377EB8",
  BEPA = "#984EA3", FAGR = "#A65628", FRAM = "#F781BF", PIST = "#66C2A5",
  QURU = "#E7298A", TSCA = "#1B9E77"
)

p3 <- ggplot(df_16s, aes(Axis1, Axis2, color = species)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.7) +
  scale_color_manual(values = sp_colors) +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_16s[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_16s[2]),
       title = "16S Bacteria — Tree Species",
       color = "Species", shape = "Horizon") +
  theme_ord

p4 <- ggplot(df_its, aes(Axis1, Axis2, color = species)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.7) +
  scale_color_manual(values = sp_colors) +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_its[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_its[2]),
       title = "ITS Fungi — Tree Species",
       color = "Species", shape = "Horizon") +
  theme_ord

# ── Plot: colored by sampling location (plot) ──
site_colors <- c("Lab" = "#E41A1C", "Lowland" = "#377EB8", "Ridge" = "#4DAF4A")

p5 <- ggplot(df_16s, aes(Axis1, Axis2, color = site)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.7) +
  scale_color_manual(values = site_colors, na.value = "grey50") +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_16s[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_16s[2]),
       title = "16S Bacteria — Sampling Location",
       color = "Site", shape = "Horizon") +
  theme_ord

p6 <- ggplot(df_its, aes(Axis1, Axis2, color = site)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.7) +
  scale_color_manual(values = site_colors, na.value = "grey50") +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_its[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_its[2]),
       title = "ITS Fungi — Sampling Location",
       color = "Site", shape = "Horizon") +
  theme_ord

# ── Plot: colored by mycorrhizal type ──
p7 <- ggplot(df_16s, aes(Axis1, Axis2, color = mycorrhiza)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("AM" = "#D95F02", "ECM" = "#1B9E77", "AM/ECM" = "#7570B3")) +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_16s[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_16s[2]),
       title = "16S Bacteria — Mycorrhizal Type",
       color = "Mycorrhiza", shape = "Horizon") +
  theme_ord

p8 <- ggplot(df_its, aes(Axis1, Axis2, color = mycorrhiza)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("AM" = "#D95F02", "ECM" = "#1B9E77", "AM/ECM" = "#7570B3")) +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_its[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_its[2]),
       title = "ITS Fungi — Mycorrhizal Type",
       color = "Mycorrhiza", shape = "Horizon") +
  theme_ord

# ── Plot: colored by ORP ──
p9 <- ggplot(df_16s, aes(Axis1, Axis2, color = ORP)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "magma") +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_16s[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_16s[2]),
       title = "16S Bacteria — Soil ORP",
       color = "ORP (mV)", shape = "Horizon") +
  theme_ord

p10 <- ggplot(df_its, aes(Axis1, Axis2, color = ORP)) +
  geom_point(aes(shape = horizon), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "magma") +
  labs(x = sprintf("PCoA1 (%.1f%%)", eig_its[1]),
       y = sprintf("PCoA2 (%.1f%%)", eig_its[2]),
       title = "ITS Fungi — Soil ORP",
       color = "ORP (mV)", shape = "Horizon") +
  theme_ord

# Save combined figure
fig_combined <- plot_grid(
  p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
  ncol = 2, align = "hv"
)
ggsave(file.path(out_dir, "fig_soil_exploration_ordination.png"),
       fig_combined, width = 16, height = 30, dpi = 200)
cat("Saved ordination figure.\n")

################################################################################
# 7. VARIANCE PARTITIONING (4-way: species, VWC+ORP, horizon, plot)
################################################################################

cat("\n========== VARIANCE PARTITIONING ==========\n")

bc_16s_mat <- as.matrix(bc_16s)
bc_its_mat <- as.matrix(bc_its)

# 4-way: Species, VWC+ORP, Horizon, Plot
varpart_16s <- varpart(
  bc_16s_mat,
  ~ species.x,          # [a] species
  ~ VWC + ORP,          # [b] moisture/redox
  ~ horizon,            # [c] horizon
  ~ site,               # [d] sampling location
  data = sd16c
)
cat("\n--- 16S Variance Partitioning (4-way) ---\n")
print(varpart_16s)

varpart_its <- varpart(
  bc_its_mat,
  ~ species.x,
  ~ VWC + ORP,
  ~ horizon,
  ~ site,
  data = sditc
)
cat("\n--- ITS Variance Partitioning (4-way) ---\n")
print(varpart_its)

# Save varpart plots
png(file.path(out_dir, "fig_soil_exploration_varpart.png"),
    width = 14, height = 7, units = "in", res = 200)
par(mfrow = c(1, 2))
plot(varpart_16s, Xnames = c("Species", "VWC+ORP", "Horizon", "Plot"),
     main = "16S Bacteria",
     bg = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"))
plot(varpart_its, Xnames = c("Species", "VWC+ORP", "Horizon", "Plot"),
     main = "ITS Fungi",
     bg = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00"))
dev.off()
cat("Saved variance partitioning figure.\n")

################################################################################
# 8. ALPHA DIVERSITY vs VWC
################################################################################

cat("\n========== ALPHA DIVERSITY ==========\n")

# 16S
alpha_16s <- estimate_richness(ps_16s_complete, measures = c("Shannon", "Observed"))
alpha_16s$VWC <- as.numeric(sd16c$mean_vwc)
alpha_16s$species <- sd16c$species.x
alpha_16s$horizon <- ifelse(sd16c$core_type == "Mineral", "Mineral", "Organic")
alpha_16s$mycorrhiza <- myco_map[alpha_16s$species]
alpha_16s$site <- sd16c$site

# ITS
alpha_its <- estimate_richness(ps_its_complete, measures = c("Shannon", "Observed"))
alpha_its$VWC <- as.numeric(sditc$mean_vwc)
alpha_its$species <- sditc$species.x
alpha_its$horizon <- ifelse(sditc$core_type == "Mineral", "Mineral", "Organic")
alpha_its$mycorrhiza <- myco_map[alpha_its$species]
alpha_its$site <- sditc$site

# Shannon vs VWC (colored by plot)
pa1 <- ggplot(alpha_16s, aes(VWC, Shannon, color = site)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = site_colors, na.value = "grey50") +
  labs(title = "16S Shannon vs. VWC", x = "VWC (%)", y = "Shannon Index") +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank())

pa2 <- ggplot(alpha_its, aes(VWC, Shannon, color = site)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = site_colors, na.value = "grey50") +
  labs(title = "ITS Shannon vs. VWC", x = "VWC (%)", y = "Shannon Index") +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank())

# Richness vs VWC
pa3 <- ggplot(alpha_16s, aes(VWC, Observed, color = site)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = site_colors, na.value = "grey50") +
  labs(title = "16S Richness vs. VWC", x = "VWC (%)", y = "Observed ASVs") +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank())

pa4 <- ggplot(alpha_its, aes(VWC, Observed, color = site)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_manual(values = site_colors, na.value = "grey50") +
  labs(title = "ITS Richness vs. VWC", x = "VWC (%)", y = "Observed ASVs") +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank())

# Shannon by species (colored by mycorrhiza)
pa5 <- ggplot(alpha_16s, aes(reorder(species, Shannon, median), Shannon, fill = mycorrhiza)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c("AM" = "#D95F02", "ECM" = "#1B9E77", "AM/ECM" = "#7570B3")) +
  coord_flip() +
  labs(title = "16S Shannon by Tree Species", x = "", y = "Shannon Index", fill = "Mycorrhiza") +
  theme_bw(base_size = 11) + theme(panel.grid.major.y = element_blank())

pa6 <- ggplot(alpha_its, aes(reorder(species, Shannon, median), Shannon, fill = mycorrhiza)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  scale_fill_manual(values = c("AM" = "#D95F02", "ECM" = "#1B9E77", "AM/ECM" = "#7570B3")) +
  coord_flip() +
  labs(title = "ITS Shannon by Tree Species", x = "", y = "Shannon Index", fill = "Mycorrhiza") +
  theme_bw(base_size = 11) + theme(panel.grid.major.y = element_blank())

fig_alpha <- plot_grid(pa1, pa2, pa3, pa4, pa5, pa6, ncol = 2, align = "hv")
ggsave(file.path(out_dir, "fig_soil_exploration_alpha.png"),
       fig_alpha, width = 14, height = 16, dpi = 200)
cat("Saved alpha diversity figure.\n")

################################################################################
# 9. BETA-DISPERSION
################################################################################

cat("\n========== BETA-DISPERSION ==========\n")

bd_16s_sp <- betadisper(bc_16s, sd16c$species.x)
bd_its_sp <- betadisper(bc_its, sditc$species.x)

cat("\n--- 16S betadisper by species ---\n")
print(permutest(bd_16s_sp, pairwise = FALSE))
cat("\n--- ITS betadisper by species ---\n")
print(permutest(bd_its_sp, pairwise = FALSE))

# By plot
bd_16s_pl <- betadisper(bc_16s, sd16c$site)
bd_its_pl <- betadisper(bc_its, sditc$site)

cat("\n--- 16S betadisper by plot ---\n")
print(permutest(bd_16s_pl, pairwise = FALSE))
cat("\n--- ITS betadisper by plot ---\n")
print(permutest(bd_its_pl, pairwise = FALSE))

# By horizon
bd_16s_hz <- betadisper(bc_16s, sd16c$core_type)
bd_its_hz <- betadisper(bc_its, sditc$core_type)

cat("\n--- 16S betadisper by horizon ---\n")
print(permutest(bd_16s_hz, pairwise = FALSE))
cat("\n--- ITS betadisper by horizon ---\n")
print(permutest(bd_its_hz, pairwise = FALSE))

################################################################################
# 10. MYCORRHIZAL TYPE (with plot conditioning)
################################################################################

cat("\n========== MYCORRHIZAL TYPE PERMANOVA ==========\n")

sd16c$mycorrhiza <- myco_map[sd16c$species.x]
sditc$mycorrhiza <- myco_map[sditc$species.x]

set.seed(999)
perm_16s_myco <- adonis2(
  bc_16s ~ mycorrhiza + VWC + ORP + horizon,
  data = sd16c, permutations = 999, by = "margin",
  strata = sd16c$site
)
cat("\n--- 16S: Mycorrhizal type (conditioned on plot) ---\n")
print(perm_16s_myco)

set.seed(999)
perm_its_myco <- adonis2(
  bc_its ~ mycorrhiza + VWC + ORP + horizon,
  data = sditc, permutations = 999, by = "margin",
  strata = sditc$site
)
cat("\n--- ITS: Mycorrhizal type (conditioned on plot) ---\n")
print(perm_its_myco)

################################################################################
# 11. SUMMARY
################################################################################

cat("\n\n")
cat("================================================================\n")
cat("                    EXPLORATION SUMMARY                         \n")
cat("================================================================\n\n")

cat("DATASET SIZE:\n")
cat(sprintf("  16S soil samples (complete): %d\n", nsamples(ps_16s_complete)))
cat(sprintf("  ITS soil samples (complete): %d\n", nsamples(ps_its_complete)))
cat(sprintf("  Tree species (all): %d\n", length(well_sampled)))
cat(sprintf("  Tree species (multi-plot): %d (%s)\n",
            length(multi_plot_species), paste(multi_plot_species, collapse = ", ")))

cat("\n--- KEY COMPARISON: Species effect WITHOUT vs WITH plot conditioning ---\n\n")
cat("  16S Bacteria:\n")
cat(sprintf("    Without plot: R²=%.3f, p=%.3f\n",
            perm_16s_noplot["species.x", "R2"], perm_16s_noplot["species.x", "Pr(>F)"]))
cat(sprintf("    With stratified by site: R²=%.3f, p=%.3f\n",
            perm_16s_cond_plot["species.x", "R2"], perm_16s_cond_plot["species.x", "Pr(>F)"]))
cat(sprintf("    Multi-plot subset + stratified by site: R²=%.3f, p=%.3f\n",
            perm_16s_multi["species.x", "R2"], perm_16s_multi["species.x", "Pr(>F)"]))

cat("\n  ITS Fungi:\n")
cat(sprintf("    Without plot: R²=%.3f, p=%.3f\n",
            perm_its_noplot["species.x", "R2"], perm_its_noplot["species.x", "Pr(>F)"]))
cat(sprintf("    With stratified by site: R²=%.3f, p=%.3f\n",
            perm_its_cond_plot["species.x", "R2"], perm_its_cond_plot["species.x", "Pr(>F)"]))
cat(sprintf("    Multi-plot subset + stratified by site: R²=%.3f, p=%.3f\n",
            perm_its_multi["species.x", "R2"], perm_its_multi["species.x", "Pr(>F)"]))

cat("\n--- VWC/ORP (conditioned on plot) ---\n")
cat(sprintf("  16S VWC: R²=%.3f (p=%.3f)  |  ITS VWC: R²=%.3f (p=%.3f)\n",
            perm_16s_cond_plot["VWC", "R2"], perm_16s_cond_plot["VWC", "Pr(>F)"],
            perm_its_cond_plot["VWC", "R2"], perm_its_cond_plot["VWC", "Pr(>F)"]))
cat(sprintf("  16S ORP: R²=%.3f (p=%.3f)  |  ITS ORP: R²=%.3f (p=%.3f)\n",
            perm_16s_cond_plot["ORP", "R2"], perm_16s_cond_plot["ORP", "Pr(>F)"],
            perm_its_cond_plot["ORP", "R2"], perm_its_cond_plot["ORP", "Pr(>F)"]))

cat("\n--- Mycorrhizal type (conditioned on plot) ---\n")
cat(sprintf("  16S: R²=%.3f (p=%.3f)  |  ITS: R²=%.3f (p=%.3f)\n",
            perm_16s_myco["mycorrhiza", "R2"], perm_16s_myco["mycorrhiza", "Pr(>F)"],
            perm_its_myco["mycorrhiza", "R2"], perm_its_myco["mycorrhiza", "Pr(>F)"]))

cat("\n--- Species×VWC interaction (conditioned on plot) ---\n")
if ("species.x:VWC" %in% rownames(perm_16s_interact)) {
  cat(sprintf("  16S: R²=%.3f (p=%.3f)  |  ITS: R²=%.3f (p=%.3f)\n",
              perm_16s_interact["species.x:VWC", "R2"],
              perm_16s_interact["species.x:VWC", "Pr(>F)"],
              perm_its_interact["species.x:VWC", "R2"],
              perm_its_interact["species.x:VWC", "Pr(>F)"]))
}

cat("\nFigures saved to:", out_dir, "\n")
cat("  - fig_soil_exploration_ordination.png (10-panel PCoA incl. plot)\n")
cat("  - fig_soil_exploration_varpart.png (Venn diagrams with plot)\n")
cat("  - fig_soil_exploration_alpha.png (alpha diversity)\n")
cat("\nDone!\n")
