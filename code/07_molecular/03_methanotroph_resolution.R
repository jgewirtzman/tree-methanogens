#!/usr/bin/env Rscript
# ==============================================================================
# 03_methanotroph_resolution.R -- methanotrophs in the 16S data at several
#   taxonomic resolutions, and which resolution tracks the ddPCR genes
# ------------------------------------------------------------------------------
# The main analysis counts methanotrophs at genus level (Known) with a family-level
# upper bound (Putative). Genus calls can miss methanotrophs that SILVA leaves
# unresolved, so this script asks the same question several ways:
#
#   Order   Methylococcales; Methylacidiphilales; Methylomirabilales (NC10)
#   Family  Methylomonadaceae + Methylococcaceae; Beijerinckiaceae (all);
#           Methylacidiphilaceae; Methylomirabilaceae
#   Genus   Known (the capacity rule); Known + Putative
#   Placed  Known, plus genus-unresolved Beijerinckiaceae ASVs whose nearest
#           genus-resolved Beijerinckiaceae neighbour on the ASV tree is a
#           methanotroph genus (Methylocapsa, Methylocella, Methyloferula), within
#           PLACE_MAX_DIST substitutions per site
#
# For each grouping and compartment it reports the mean share of the community,
# prevalence, and Spearman correlations with ddPCR pmoA and mmoX in the same
# extracts (all samples, zeros included). The felled black oak is summarised by
# tissue at the order, family and genus levels (its ASVs are not on the tree).
#
# Reads:  data/raw/16s/OTU_table.txt, unrooted_tree.nwk,
#         data/raw/external/tree-microbiome/tree_data_methanogen_group.csv,
#         data/raw/16s/black_oak/OTU_table.txt, code/lib/methanotroph_definitions.csv
# Writes: outputs/data/methanotroph_resolution.csv (survey, by compartment)
#         outputs/data/methanotroph_resolution_placement.csv (unresolved ASVs)
#         outputs/data/methanotroph_placed_asvs.csv (ASV ids for the Placed tier, read by
#           load_placed_asvs() in every downstream classifier call)
#         outputs/data/methanotroph_placement_sensitivity.csv (thresholds 0.02-0.04)
#         outputs/data/methanotroph_resolution_black_oak.csv (oak, by tissue)
# ==============================================================================
suppressMessages({ library(ape); library(dplyr) })
source("code/lib/outputs.R"); source("code/lib/load_methanotroph_definitions.R")
num <- function(x) suppressWarnings(as.numeric(x))
PLACE_MAX_DIST <- 0.03          # ~97% identity over V4
MT_GENERA <- c("Methylocapsa", "Methylocella", "Methyloferula")
NONMT_GENERA <- c("Methylobacterium-Methylorubrum", "1174-901-12", "Roseiarcus", "Bosea", "Microvirga")
unres <- function(g) is.na(g) | g == "" | tolower(g) %in% c("none", "unclassified", "uncultured")

# ---- survey 16S -------------------------------------------------------------
otu <- read.delim("data/raw/16s/OTU_table.txt", header = TRUE, row.names = 1, check.names = FALSE)
taxc <- intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), names(otu))
tax <- otu[, taxc]; for (k in taxc) tax[[k]] <- trimws(tax[[k]])
cnt <- as.matrix(sapply(otu[, setdiff(names(otu), taxc), drop = FALSE], num)); rownames(cnt) <- rownames(otu)
tot <- colSums(cnt, na.rm = TRUE)
keep_s <- tot >= 1000; cnt <- cnt[, keep_s]; tot <- tot[keep_s]
rel <- sweep(cnt, 2, tot, "/") * 100
defs <- load_methanotroph_defs()
tax$mt <- classify_methanotrophs(tax, defs, include_conditional = FALSE)

# ---- phylogenetic placement of genus-unresolved Beijerinckiaceae -----------
tr <- read.tree("data/raw/16s/unrooted_tree.nwk"); tr$tip.label <- gsub("'", "", tr$tip.label)
beij <- rownames(tax)[tax$Family %in% "Beijerinckiaceae" & rownames(tax) %in% tr$tip.label]
tb <- keep.tip(tr, beij); D <- cophenetic(tb)
resolved <- beij[!unres(tax[beij, "Genus"])]; unresolved <- beij[unres(tax[beij, "Genus"])]
P <- do.call(rbind, lapply(unresolved, function(a) {
  d <- D[a, resolved]; j <- which.min(d)
  data.frame(asv = a, nearest_genus = tax[resolved[j], "Genus"], dist = unname(d[j]))
}))
place <- function(thr) ifelse(P$dist > thr, "no close relative",
            ifelse(P$nearest_genus %in% MT_GENERA, "methanotroph genus",
            ifelse(P$nearest_genus %in% NONMT_GENERA, "non-methanotroph genus", "other genus")))
P$placed <- place(PLACE_MAX_DIST)
placed_mt <- P$asv[P$placed == "methanotroph genus"]
write.csv(data.frame(asv = placed_mt), out_path("methanotroph_placed_asvs.csv"), row.names = FALSE)

# ---- groupings ----------------------------------------------------------------
fam <- function(f) rownames(tax)[tax$Family %in% f]
ord <- function(o) rownames(tax)[tax$Order %in% o]
GROUPS <- list(
  "Order: Methylococcales"                     = ord("Methylococcales"),
  "Order: Methylacidiphilales"                 = ord("Methylacidiphilales"),
  "Order: Methylomirabilales (NC10)"           = ord("Methylomirabilales"),
  "Family: Methylomonadaceae + Methylococcaceae" = fam(c("Methylomonadaceae", "Methylococcaceae")),
  "Family: Beijerinckiaceae (all)"             = fam("Beijerinckiaceae"),
  "Family: Methylacidiphilaceae"               = fam("Methylacidiphilaceae"),
  "Family: Methylomirabilaceae"                = fam("Methylomirabilaceae"),
  "Genus: Known"                               = rownames(tax)[tax$mt == "Known"],
  "Genus: Known + Putative (all tiers)"        = rownames(tax)[tax$mt %in% c("Known", "Putative")],
  "Placed: Known + Placed"                     = union(rownames(tax)[tax$mt == "Known"], placed_mt))

# ---- pair with ddPCR ----------------------------------------------------------
key <- sub("[.]16S[.]S[0-9]*$", "", colnames(rel))
o <- read.csv("data/raw/external/tree-microbiome/tree_data_methanogen_group.csv", check.names = FALSE); o <- o[, names(o) != ""]
o$key <- paste0(o$seq_id, o$core_type)
COMP <- c(Inner = "Heartwood", Outer = "Sapwood", Mineral = "Mineral soil", Organic = "Organic soil")
o$comp <- COMP[o$core_type]
idx <- match(key, o$key)
smp <- data.frame(col = colnames(rel), comp = o$comp[idx], pmoa = num(o$pmoa_loose[idx]), mmox = num(o$mmox_loose[idx]))
smp <- smp[!is.na(smp$comp), ]

R <- do.call(rbind, lapply(names(GROUPS), function(g) {
  a <- intersect(GROUPS[[g]], rownames(rel))
  share <- if (length(a)) colSums(rel[a, smp$col, drop = FALSE]) else rep(0, nrow(smp))
  do.call(rbind, lapply(unname(COMP), function(k) {
    s <- smp$comp == k; x <- share[s]
    rho <- function(y) { ok <- is.finite(y[s]) & is.finite(x); if (sum(ok) < 8 || sd(x[ok]) == 0) return(c(NA, NA))
      ct <- suppressWarnings(cor.test(x[ok], y[s][ok], method = "spearman", exact = FALSE)); c(unname(ct$estimate), ct$p.value) }
    rp <- rho(smp$pmoa); rm <- rho(smp$mmox)
    data.frame(grouping = g, compartment = k, n = sum(s), ASVs = length(a),
               mean_pct = mean(x), prevalence_pct = 100 * mean(x > 0),
               rho_pmoA = rp[1], p_pmoA = rp[2], rho_mmoX = rm[1], p_mmoX = rm[2])
  }))
}))
write.csv(R, out_path("methanotroph_resolution.csv"), row.names = FALSE)

# placement summary, weighted by abundance in each compartment
P$heartwood <- rowMeans(rel[P$asv, smp$col[smp$comp == "Heartwood"], drop = FALSE])
P$sapwood   <- rowMeans(rel[P$asv, smp$col[smp$comp == "Sapwood"], drop = FALSE])
P$mineral   <- rowMeans(rel[P$asv, smp$col[smp$comp == "Mineral soil"], drop = FALSE])
P$organic   <- rowMeans(rel[P$asv, smp$col[smp$comp == "Organic soil"], drop = FALSE])
write.csv(P, out_path("methanotroph_resolution_placement.csv"), row.names = FALSE)

# sensitivity of the Placed tier to the distance threshold
known_ids <- rownames(tax)[tax$mt == "Known"]
SENS <- do.call(rbind, lapply(c(0.02, 0.03, 0.04), function(thr) {
  ids <- union(known_ids, P$asv[place(thr) == "methanotroph genus"])
  a <- intersect(ids, rownames(rel)); sh <- colSums(rel[a, smp$col, drop = FALSE])
  data.frame(threshold = thr, placed_ASVs = sum(place(thr) == "methanotroph genus"),
             t(sapply(unname(COMP), function(k) mean(sh[smp$comp == k]))), check.names = FALSE)
}))
write.csv(SENS, out_path("methanotroph_placement_sensitivity.csv"), row.names = FALSE)

# ---- black oak, by tissue -------------------------------------------------------
bo <- read.table("data/raw/16s/black_oak/OTU_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bnum <- names(bo)[sapply(bo, is.numeric)]; btot <- colSums(bo[bnum], na.rm = TRUE); bnum <- bnum[btot >= 1000]
brel <- sweep(as.matrix(bo[bnum]), 2, btot[bnum], "/") * 100
for (k in c("Order", "Family", "Genus")) bo[[k]] <- trimws(bo[[k]])
bo$mt <- classify_methanotrophs(data.frame(Family = bo$Family, Genus = bo$Genus, stringsAsFactors = FALSE), defs, include_conditional = FALSE)
tissue <- function(s) dplyr::case_when(grepl("HEART", s) ~ "Heartwood", grepl("SAP", s) ~ "Sapwood", grepl("BARK", s) ~ "Bark",
  grepl("FOLIAGE", s) ~ "Foliage", grepl("BRANCH", s) ~ "Branch", grepl("LITTER", s) ~ "Litter", grepl("MINERAL", s) ~ "Mineral soil",
  grepl("ORGANIC", s) ~ "Organic soil", grepl("ROT", s) ~ "Rot", grepl("COARSE", s) ~ "Coarse root", grepl("FINE", s) ~ "Fine root", TRUE ~ NA_character_)
BG <- list("Order: Methylococcales" = bo$Order %in% "Methylococcales",
           "Order: Methylacidiphilales" = bo$Order %in% "Methylacidiphilales",
           "Family: Beijerinckiaceae (all)" = bo$Family %in% "Beijerinckiaceae",
           "Genus: Known" = bo$mt %in% "Known", "Genus: Known + Putative" = bo$mt %in% c("Known", "Putative"),
           "Genus: Methylocella" = bo$Genus %in% "Methylocella", "Genus: Methylocapsa" = bo$Genus %in% "Methylocapsa")
tis <- tissue(bnum)
B <- do.call(rbind, lapply(names(BG), function(g) {
  sh <- colSums(brel[BG[[g]], , drop = FALSE])
  do.call(rbind, lapply(sort(unique(na.omit(tis))), function(t) {
    x <- sh[tis %in% t]
    data.frame(grouping = g, tissue = t, n = length(x), mean_pct = mean(x), samples_detected = sum(x > 0))
  }))
}))
write.csv(B, out_path("methanotroph_resolution_black_oak.csv"), row.names = FALSE)

# ---- print ----------------------------------------------------------------------
fmt <- R; fmt$mean_pct <- round(fmt$mean_pct, 3); fmt$prevalence_pct <- round(fmt$prevalence_pct)
for (v in c("rho_pmoA", "rho_mmoX")) fmt[[v]] <- round(fmt[[v]], 2)
for (v in c("p_pmoA", "p_mmoX")) fmt[[v]] <- signif(fmt[[v]], 2)
cat("\nSURVEY: methanotroph groupings by compartment (share of community; Spearman with ddPCR in the same extracts)\n")
print(fmt, row.names = FALSE)
cat("\nPLACEMENT of genus-unresolved Beijerinckiaceae (", length(unresolved), " ASVs; max distance ", PLACE_MAX_DIST, ")\n", sep = "")
print(P %>% group_by(placed) %>% summarise(ASVs = n(), heartwood = round(sum(heartwood), 3), sapwood = round(sum(sapwood), 3),
                                          mineral = round(sum(mineral), 3), organic = round(sum(organic), 3)) %>% as.data.frame(), row.names = FALSE)
cat("\nnearest resolved genus (all unresolved ASVs):\n"); print(sort(table(P$nearest_genus), decreasing = TRUE))
cat("\nSENSITIVITY: Known + Placed (% of community) by placement threshold\n"); print(SENS, row.names = FALSE, digits = 3)
cat("\nBLACK OAK by tissue:\n"); B$mean_pct <- round(B$mean_pct, 3); print(B, row.names = FALSE)
