# ==============================================================================
# 04_known-putative-table.R -- Table S2: methane-cycling taxa and their classification
# ------------------------------------------------------------------------------
# Stage D; run from the repository root.
# Supplementary table of the taxa classified as known and putative methanotrophs and
# methanogens. Methylacidiphilaceae: only geothermal members are methanotrophic;
# mesophilic ones are not.
#
# CLASSIFICATION LOGIC: the project's shared, genus-aware classifier
# (code/lib/load_methanotroph_definitions.R) -- the rule Figure 5 and the Results use:
#   METHANOTROPH, "Known"    = ASV genus is a cultivated methanotroph genus, OR ASV
#                              family is an EXCLUSIVE methanotroph family (all members
#                              are methanotrophs) when genus is unresolved.
#   METHANOTROPH, "Putative" = ASV family CONTAINS methanotrophs but also non-
#                              methanotrophs (e.g. Beijerinckiaceae, Methylocystaceae)
#                              and genus is unresolved -> cannot be confirmed.
#   Methylotroph-only families (Methylophilaceae, Methyloligellaceae, Methylopilaceae)
#                              are NOT methanotrophs and are excluded.
#   METHANOGEN               = ASV family is in code/lib/methanogen_families.R.
#
# METHYLACIDIPHILACEAE: reclassify from Known -> Putative,
# since verified methanotrophy is restricted to thermoacidophilic geothermal members;
# temperate-forest (mesophilic) Methylacidiphilaceae are not known to oxidize CH4.
#
# Reads:   asv_reads("main", samples_only = TRUE) (data/tables/asv_counts_16S.csv +
#          asv_taxonomy_16S.csv, via code/lib/tables_read.R): the sample libraries of the main
#          16S run (no mock, no "other" libraries); organelle ASVs removed
# Outputs: outputs/data/known_putative_taxa_table.csv, kp_counts.txt
# ==============================================================================
source("code/lib/outputs.R")

suppressPackageStartupMessages({ library(tidyverse) })

source("code/lib/load_methanotroph_definitions.R")
defs_rev  <- load_methanotroph_defs()                                     # revised (current)
# Previous rule, for the Methylacidiphilaceae comparison: Methylacidiphilaceae counted as an exclusive
# (Known) family. Built from the revised file rather than read from
# data/compiled/methanotroph_definitions.csv, which is now byte-identical to the revised file
# and so made the comparison report "moves 0 ASVs".
defs_orig <- defs_rev
defs_orig$Include_known[defs_orig$Taxon_rank == "Family" & defs_orig$Taxon == "Methylacidiphilaceae"] <- "YES"
source("code/lib/tables_read.R")
# real samples only (samples_16S library_type "sample": no Zymo mock, no "other"
# libraries) and organelle ASVs out of the denominator (asv_reads()); ASVs with no read
# in those libraries are not detected and are left out
A <- asv_reads("main", samples_only = TRUE)
otu_mat <- A$cnt[rowSums(A$cnt, na.rm = TRUE) > 0, , drop = FALSE]
# per-ASV total relative abundance (summed over samples, as % of all their non-organelle reads)
asv_relabund <- (rowSums(otu_mat, na.rm = TRUE) / sum(A$tot)) * 100
tx <- A$tax[rownames(otu_mat), ]
tax <- tibble(feature_id = rownames(otu_mat), phylum = tx$Phylum, family = tx$Family, genus = tx$Genus,
              relabund = asv_relabund[feature_id])

source("code/lib/methanogen_families.R")
methanogen_families <- METHANOGEN_FAMILIES

# shared classifier expects Family / Genus / Phylum
tdf <- data.frame(Family = tax$family, Genus = tax$genus, Phylum = tax$phylum, stringsAsFactors = FALSE)
label <- function(mt_status) {
  out <- ifelse(is.na(mt_status), NA_character_, paste0("Methanotroph_", mt_status))
  ifelse(is.na(out) & tax$family %in% methanogen_families, "Methanogen", out)
}
tax$class_orig <- label(classify_methanotrophs(tdf, defs_orig))
tax$class_rev  <- label(classify_methanotrophs(tdf, defs_rev))

# ---- per-taxon supplementary table (detected taxa) ---------------------------
tab <- tax %>% filter(!is.na(class_rev)) %>%
  mutate(genus_res = genus,                                       # NA = unresolved
         display_taxon = coalesce(genus_res, family),
         level = if_else(!is.na(genus_res), "Genus", "Family")) %>%
  group_by(classification = class_rev, family, display_taxon, level) %>%
  summarise(n_ASVs = n(), mean_relabund_pct = round(sum(relabund), 4), .groups = "drop") %>%
  left_join(defs_rev %>% transmute(display_taxon = Taxon, source = Primary_source, note = Notes),
            by = "display_taxon") %>%
  arrange(classification, desc(mean_relabund_pct))
write.csv(tab, out_path("known_putative_taxa_table.csv"), row.names = FALSE)

# ---- counts: original vs Methylacidiphilaceae-revised ------------------------
cnt <- function(col) tax %>% filter(!is.na(.data[[col]])) %>% count(.data[[col]]) %>% deframe()
c0 <- cnt("class_orig"); c1 <- cnt("class_rev")

sink(out_path("kp_counts.txt"))
cat("=================================================================\n")
cat("Known/Putative methanotroph & methanogen classification\n")
cat("=================================================================\n\n")
cat("Methanotroph ASV counts (n ASVs):\n")
cat(sprintf("  Original scheme : Known = %d, Putative = %d  (total %d)\n",
            c0["Methanotroph_Known"], c0["Methanotroph_Putative"],
            c0["Methanotroph_Known"]+c0["Methanotroph_Putative"]))
cat(sprintf("  Methylacidi-REVISED: Known = %d, Putative = %d  (total %d)\n",
            c1["Methanotroph_Known"], c1["Methanotroph_Putative"],
            c1["Methanotroph_Known"]+c1["Methanotroph_Putative"]))
cat(sprintf("  -> reclassifying Methylacidiphilaceae Known->Putative moves %d ASVs.\n\n",
            c0["Methanotroph_Known"]-c1["Methanotroph_Known"]))
cat(sprintf("Methanogen ASVs (n): %d\n\n", c1["Methanogen"]))
cat("Top methanotroph taxa (revised classification), by mean relative abundance:\n")
print(as.data.frame(tab %>% filter(grepl("Methanotroph", classification)) %>%
  transmute(classification = sub("Methanotroph_","",classification), display_taxon, level, n_ASVs, mean_relabund_pct) %>%
  head(12)), row.names = FALSE)
cat("\nMethanogen taxa:\n")
print(as.data.frame(tab %>% filter(classification=="Methanogen") %>%
  transmute(display_taxon, level, n_ASVs, mean_relabund_pct)), row.names = FALSE)
cat("\nNOTE: Methylacidiphilaceae now PUTATIVE (was Known) -- verified methanotrophy\n")
cat("is restricted to thermoacidophilic geothermal members; mesophilic forest members are not\n")
cat("known to oxidize CH4 (peatland & tree-stem genomes lack methanotrophy genes).\n")
sink()
cat(readLines(out_path("kp_counts.txt")), sep="\n")
