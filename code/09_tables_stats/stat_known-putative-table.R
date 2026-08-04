source("code/lib/outputs.R")
# ==============================================================================
# REVISION — Known/Putative methanotroph & methanogen taxa table (Referee 2 #2)
# ==============================================================================
# R2: "provide a supplementary table for the taxa classified as known and putative
# methanotrophs and methanogens", and (R2 #2b) Methylacidiphilaceae: only geothermal
# members are methanotrophic; mesophilic ones are not.
#
# CLASSIFICATION LOGIC (made explicit; from methanotroph_definitions.csv, Knief 2015
# + SILVA 138 audit):
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

defs <- read_csv("data/compiled/methanotroph_definitions.csv", show_col_types = FALSE) %>%
  mutate(across(c(Include_known, Include_putative), ~toupper(trimws(.))))
tax  <- read_csv("data/compiled/taxonomy_key_16S.csv", show_col_types = FALSE)
otu  <- read_csv("data/compiled/otu_table_16S.csv", show_col_types = FALSE)
otu_mat <- otu %>% column_to_rownames(names(otu)[1])
otu_mat <- otu_mat[, sapply(otu_mat, is.numeric), drop = FALSE]   # keep numeric sample cols only
# per-ASV total relative abundance (summed over samples, as % of grand total)
asv_relabund <- (rowSums(otu_mat, na.rm = TRUE) / sum(otu_mat, na.rm = TRUE)) * 100
tax <- tax %>% mutate(relabund = asv_relabund[feature_id])

known_genera <- defs %>% filter(Taxon_rank=="Genus", Include_known=="YES") %>% pull(Taxon)
excl_families <- defs %>% filter(Taxon_rank=="Family", Include_known=="YES") %>% pull(Taxon)
put_families  <- defs %>% filter(Taxon_rank=="Family", Include_putative %in% c("YES","CONDITIONAL")) %>% pull(Taxon)
methanogen_families <- c("Methanobacteriaceae","Methanomassiliicoccaceae","Methanoregulaceae",
  "Methanocellaceae","Methanosaetaceae","Methanomicrobiaceae","Methanosarcinaceae",
  "Methanomethyliaceae","Methanocorpusculaceae")

classify <- function(fam, gen, methylacidi_known = TRUE) {
  excl <- excl_families; if (!methylacidi_known) excl <- setdiff(excl, "Methylacidiphilaceae")
  put  <- union(put_families, if (!methylacidi_known) "Methylacidiphilaceae" else character(0))
  if (!is.na(gen) & gen %in% known_genera) return("Methanotroph_Known")
  if (!is.na(fam) & fam %in% excl)         return("Methanotroph_Known")
  if (!is.na(fam) & fam %in% put)          return("Methanotroph_Putative")
  if (!is.na(fam) & fam %in% methanogen_families) return("Methanogen")
  NA_character_
}

tax <- tax %>% rowwise() %>%
  mutate(class_orig = classify(family, genus, TRUE),
         class_rev  = classify(family, genus, FALSE)) %>% ungroup()

# ---- per-taxon supplementary table (detected taxa) ---------------------------
tab <- tax %>% filter(!is.na(class_rev)) %>%
  mutate(display_taxon = coalesce(na_if(genus,""), family),
         level = if_else(!is.na(genus) & genus != "", "Genus", "Family")) %>%
  group_by(classification = class_rev, family, display_taxon, level) %>%
  summarise(n_ASVs = n(), mean_relabund_pct = round(sum(relabund), 4), .groups = "drop") %>%
  left_join(defs %>% transmute(display_taxon = Taxon, source = Primary_source, note = Notes),
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
cat("\nNOTE: Methylacidiphilaceae now PUTATIVE (was Known) per R2 -- verified methanotrophy\n")
cat("is restricted to thermoacidophilic geothermal members; mesophilic forest members are not\n")
cat("known to oxidize CH4 (peatland & tree-stem genomes lack methanotrophy genes).\n")
sink()
cat(readLines(out_path("kp_counts.txt")), sep="\n")
