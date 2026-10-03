source("code/lib/outputs.R")
# ==============================================================================
# REVISION — Known/Putative methanotroph & methanogen taxa table (Referee 2 #2)
# ==============================================================================
# R2: "provide a supplementary table for the taxa classified as known and putative
# methanotrophs and methanogens", and (R2 #2b) Methylacidiphilaceae: only geothermal
# members are methanotrophic; mesophilic ones are not.
#
# CLASSIFICATION LOGIC: the project's shared, genus-aware classifier
# (code/lib/load_methanotroph_definitions.R, revised definitions file) -- the SAME
# rule Figure 5 and the Results use. Until 2026-09-30 this script carried its own
# classify(), which read the pre-revision definitions and marked every mixed-family
# ASV Putative WITHOUT checking the genus, so resolved non-methanotroph genera
# (Methylobacterium, Roseiarcus, Bosea, ...) were listed as putative methanotrophs:
#   METHANOTROPH, "Known"    = ASV genus is a cultivated methanotroph genus, OR ASV
#                              family is an EXCLUSIVE methanotroph family (all members
#                              are methanotrophs) when genus is unresolved.
#   METHANOTROPH, "Putative" = ASV family CONTAINS methanotrophs but also non-
#                              methanotrophs (e.g. Beijerinckiaceae, Methylocystaceae)
#                              and genus is unresolved -> cannot be confirmed.
#   Methylotroph-only families (Methylophilaceae, Methyloligellaceae, Methylopilaceae)
#                              are NOT methanotrophs and are excluded.
#   METHANOGEN               = ASV family is one of 9 established methanogen families.
#
# METHYLACIDIPHILACEAE REVISION (answers R2 #2b): reclassify from Known -> Putative,
# since verified methanotrophy is restricted to thermoacidophilic geothermal members;
# temperate-forest (mesophilic) Methylacidiphilaceae are not known to oxidize CH4.
#
# NEW file. Run: Rscript code/revision/R_known_putative_table.R
# Outputs: outputs/data/known_putative_taxa_table.csv, kp_counts.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out <- "outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

source("code/lib/load_methanotroph_definitions.R")
defs_rev  <- load_methanotroph_defs()                                     # revised (current)
# Pre-revision rule, for the R2 #2b comparison: Methylacidiphilaceae counted as an exclusive
# (Known) family. Built from the revised file rather than read from
# data/compiled/methanotroph_definitions.csv, which is now byte-identical to the revised file
# and so made the comparison report "moves 0 ASVs".
defs_orig <- defs_rev
defs_orig$Include_known[defs_orig$Taxon_rank == "Family" & defs_orig$Taxon == "Methylacidiphilaceae"] <- "YES"
tax  <- read_csv("data/compiled/taxonomy_key_16S.csv", show_col_types = FALSE)
otu  <- read_csv("data/compiled/otu_table_16S.csv", show_col_types = FALSE)
otu_mat <- otu %>% column_to_rownames(names(otu)[1])
otu_mat <- otu_mat[, sapply(otu_mat, is.numeric), drop = FALSE]   # keep numeric sample cols only
# per-ASV total relative abundance (summed over samples, as % of grand total)
asv_relabund <- (rowSums(otu_mat, na.rm = TRUE) / sum(otu_mat, na.rm = TRUE)) * 100
tax <- tax %>% mutate(relabund = asv_relabund[feature_id])

methanogen_families <- c("Methanobacteriaceae","Methanomassiliicoccaceae","Methanoregulaceae",
  "Methanocellaceae","Methanosaetaceae","Methanomicrobiaceae","Methanosarcinaceae",
  "Methanomethyliaceae","Methanocorpusculaceae")

# shared classifier expects Family / Genus / Phylum
tdf <- data.frame(Family = tax$family, Genus = tax$genus, Phylum = tax$phylum, stringsAsFactors = FALSE)
rownames(tdf) <- tax$feature_id
PLACED <- load_placed_asvs()
label <- function(mt_status) {
  out <- ifelse(is.na(mt_status), NA_character_, paste0("Methanotroph_", mt_status))
  ifelse(is.na(out) & tax$family %in% methanogen_families, "Methanogen", out)
}
tax$class_orig <- label(classify_methanotrophs(tdf, defs_orig, placed = PLACED))
tax$class_rev  <- label(classify_methanotrophs(tdf, defs_rev, placed = PLACED))

# ---- per-taxon supplementary table (detected taxa) ---------------------------
tab <- tax %>% filter(!is.na(class_rev)) %>%
  mutate(genus_res = na_if(na_if(genus, ""), "none"),          # "none" = unresolved
         display_taxon = if_else(class_rev == "Methanotroph_Placed", paste0(family, " (placed with methanotroph genera)"),
                                 coalesce(genus_res, family)),
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
cat("Known/Putative methanotroph & methanogen classification (R2 #2)\n")
cat("=================================================================\n\n")
cat("Methanotroph ASV counts (n ASVs):\n")
cat(sprintf("  Original scheme : Known = %d, Putative = %d  (total %d)\n",
            c0["Methanotroph_Known"], c0["Methanotroph_Putative"],
            c0["Methanotroph_Known"]+c0["Methanotroph_Putative"]))
cat(sprintf("  Current rule: Known = %d, Placed = %d, Putative = %d  (total %d)\n",
            c1["Methanotroph_Known"], c1["Methanotroph_Placed"], c1["Methanotroph_Putative"],
            sum(c1[c("Methanotroph_Known","Methanotroph_Placed","Methanotroph_Putative")], na.rm = TRUE)))
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
cat("\nNOTE: genus-unresolved Methylacidiphilaceae are listed, not counted (amended 2026-10-02):\n")
cat("verified methanotrophy is restricted to thermoacidophilic geothermal members; mesophilic\n")
cat("members from peat and tree bark lack methanotrophy genes.\n")
sink()
cat(readLines(out_path("kp_counts.txt")), sep="\n")
