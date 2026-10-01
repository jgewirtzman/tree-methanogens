source("code/lib/outputs.R")
# ==============================================================================
# stat_clade-census.R
# ------------------------------------------------------------------------------
# For every category the paper defines -- each methanotroph rule in the definitions
# file, the groups listed but not counted, and each methanogen family -- how many
# ASVs actually fall in it, and where. Also sweeps ALL archaeal families present,
# so a group we define nothing for (e.g. ANME, which carry mcrA) cannot hide.
#
# Basis: Figure 5's own pipeline (build_16s_phyloseq, seed 46814, rarefied to
# 3,500; heartwood / sapwood / mineral / organic samples with a species). A second
# pass counts ASVs anywhere in the full dataset, including samples outside the
# analysis (black-oak tissues, muck, aquatic), unrarefied.
#
# NEW file; edits nothing. Output: outputs/audit/clade_census.txt
# ==============================================================================
suppressMessages({library(phyloseq); library(dplyr)})
source("code/lib/build_phyloseq.R"); source("code/lib/load_methanotroph_definitions.R")
defs <- load_methanotroph_defs()

ps <- suppressWarnings(build_16s_phyloseq()$ps); taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
TX_ALL <- data.frame(tax_table(ps), stringsAsFactors = FALSE)
names(TX_ALL)[1:7] <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
reads_all <- taxa_sums(ps)

set.seed(46814); pr <- suppressMessages(rarefy_even_depth(ps, sample.size = 3500, verbose = FALSE))
sd <- data.frame(sample_data(pr))
sd$comp <- c(Inner="Heartwood", Outer="Sapwood", Mineral="Mineral", Organic="Organic")[sd$core_type]
keep <- !is.na(sd$comp) & !is.na(sd$species.x) & sd$species.x != ""
pr <- prune_samples(rownames(sd)[keep], pr); sd <- sd[keep, ]
otu <- as(otu_table(pr), "matrix"); if (!taxa_are_rows(pr)) otu <- t(otu)
otu <- otu / 3500 * 100                                     # % of each sample
TX <- TX_ALL[rownames(otu), ]
unres <- function(g) is.na(g) | g == "" | tolower(g) %in% c("none","unclassified")
COMPS <- c("Heartwood","Sapwood","Mineral","Organic")

census <- function(label, idx_all, tier) {
  idx_all[is.na(idx_all)] <- FALSE; names(idx_all) <- rownames(TX_ALL)   # align by ASV name
  idx <- idx_all[rownames(otu)]
  out <- data.frame(tier = tier, category = label,
                    ASVs_anywhere = sum(idx_all & reads_all > 0),
                    ASVs_analysed = sum(idx & rowSums(otu) > 0))
  for (k in COMPS) {
    s <- sd$comp == k
    m <- if (any(idx)) colSums(otu[idx, s, drop = FALSE]) else rep(0, sum(s))
    out[[paste0(k, "_pct")]]  <- round(mean(m), 3)
    out[[paste0(k, "_prev")]] <- round(100 * mean(m > 0))
  }
  out
}
G <- function(x) TX_ALL$Genus == x & !is.na(TX_ALL$Genus)
F <- function(x) TX_ALL$Family == x & !is.na(TX_ALL$Family)
st <- classify_methanotrophs(TX_ALL, defs, include_conditional = FALSE)

rows <- list()
# ---- every rule in the definitions file ----
for (i in seq_len(nrow(defs))) {
  r <- defs[i, ]; tier <- if (r$Include_known == "YES") "Known rule" else if (r$Include_putative == "YES") "Putative rule" else "Conditional / excluded rule"
  idx <- switch(r$Taxon_rank,
    Genus  = G(r$Taxon),
    Family = if (r$Include_known == "YES") F(r$Taxon) else F(r$Taxon) & unres(TX_ALL$Genus),
    Phylum = TX_ALL$Phylum == r$Taxon & !is.na(TX_ALL$Phylum))
  lab <- paste0(r$Taxon_rank, ": ", r$Taxon, if (r$Taxon_rank == "Family" && r$Include_known != "YES") " (genus unresolved)" else "")
  rows[[length(rows) + 1]] <- census(lab, idx, tier)
}
# ---- resolved non-methanotroph genera inside mixed families (listed, not counted) ----
# (conditional families -- NC10 -- are handled separately below: their resolved genera are
#  habitat-conditional methanotrophs, not known non-methanotrophs)
mixed <- defs$Taxon[defs$Taxon_rank == "Family" & defs$Include_known != "YES" & defs$Include_putative == "YES"]
res <- TX_ALL$Family %in% mixed & !unres(TX_ALL$Genus) & is.na(st)
for (g in names(sort(table(TX_ALL$Genus[res]), decreasing = TRUE)))
  rows[[length(rows) + 1]] <- census(paste0("Genus: ", g, " (", TX_ALL$Family[res & TX_ALL$Genus == g][1], ")"),
                                     res & TX_ALL$Genus == g, "Listed, not counted: resolved non-methanotroph")
# ---- resolved genera in conditional families (NC10), not already a rule ----
cond <- defs$Taxon[defs$Taxon_rank == "Family" & defs$Include_putative == "CONDITIONAL"]
cres <- TX_ALL$Family %in% cond & !unres(TX_ALL$Genus) & !(TX_ALL$Genus %in% defs$Taxon[defs$Taxon_rank == "Genus"])
for (g in names(sort(table(TX_ALL$Genus[cres]), decreasing = TRUE)))
  rows[[length(rows) + 1]] <- census(paste0("Genus: ", g, " (", TX_ALL$Family[cres & TX_ALL$Genus == g][1], ")"),
                                     cres & TX_ALL$Genus == g, "Conditional / excluded rule")
# ---- methanogen families ----
MG <- c("Methanobacteriaceae","Methanomassiliicoccaceae","Methanoregulaceae","Methanocellaceae",
        "Methanosaetaceae","Methanomicrobiaceae","Methanosarcinaceae","Methanomethyliaceae","Methanocorpusculaceae")
for (f in MG) rows[[length(rows) + 1]] <- census(paste0("Family: ", f), F(f), "Methanogen family")
# ---- all other archaeal families present (not on the methanogen list) ----
arch <- TX_ALL$Kingdom == "Archaea" & !is.na(TX_ALL$Kingdom) & !(TX_ALL$Family %in% MG)
for (f in names(sort(tapply(reads_all[arch], TX_ALL$Family[arch], sum), decreasing = TRUE)))
  rows[[length(rows) + 1]] <- census(paste0("Family: ", f, "  [", TX_ALL$Class[arch & TX_ALL$Family == f][1], "]"),
                                     arch & TX_ALL$Family == f & !is.na(TX_ALL$Family), "Other archaea (not defined)")
C <- bind_rows(rows)

sink(out_path("clade_census.txt"))
cat("CLADE CENSUS -- every defined category: ASVs, mean % of community, prevalence (% of samples)\n")
cat(sprintf("Analysed samples (rarefied 3,500): %s\n\n", paste(names(table(sd$comp)), table(sd$comp), sep = "=", collapse = ", ")))
for (t in unique(C$tier)) { cat("====", t, "====\n"); x <- C[C$tier == t, -1]
  print(x, row.names = FALSE); cat("\n") }
sink()
cat(readLines(out_path("clade_census.txt")), sep = "\n")
