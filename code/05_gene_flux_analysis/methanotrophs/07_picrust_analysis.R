# ==============================================================================
# PICRUSt2 Functional Prediction Analysis (Figures 6, S3)
# ==============================================================================
# Purpose: PICRUSt2 functional prediction analysis linking predicted metagenome
#   to mcrA abundance.
#
# Pipeline stage: 3 — Analysis
#
# Inputs:
#   - taxonomy_table.txt (from data/raw/16s/)
#   - OTU_table.txt (from data/raw/16s/)
#   - 16S_tree_sample_table_with_meta.csv (from data/processed/)
# ==============================================================================

library(tidyverse)
library(lme4)
library(pheatmap)

## ---- Read in taxonomy and metadata ----
taxonomy = read.table('../../../data/raw/16s/taxonomy_table.txt',
                      sep = "\t", header = T, row.names = 1, quote = "",
                      stringsAsFactors = FALSE)

## Read abundance table - MAKE SURE IT'S NUMERIC
otu_table = read.table('../../../data/raw/16s/OTU_table.txt',
                       sep = "\t", header = T, row.names = 1, quote = "",
                       check.names = FALSE)

## Check what we got
cat("OTU table dimensions:", dim(otu_table), "\n")
cat("First few entries:\n")
print(otu_table[1:5, 1:5])
cat("\nClass of first column:", class(otu_table[,1]), "\n")

## Convert to numeric matrix if needed
otu_matrix = as.matrix(otu_table)
## Force to numeric
otu_matrix = apply(otu_matrix, 2, as.numeric)
rownames(otu_matrix) = rownames(otu_table)

## Transpose so samples are rows
abundance = t(otu_matrix)

cat("Abundance matrix dimensions:", dim(abundance), "\n")
cat("Class check:", class(abundance[1,1]), "\n")

## Read metadata
meta = read.csv('../../../data/raw/picrust/16S_tree_sample_table_with_meta.csv',
                row.names = 1)
meta$log16S = log10(1 + meta$X16S_per_ul)
wood = meta[meta$material == "Wood" & !(is.na(meta$material)),]
wood = wood[!(is.na(wood$mcra_probe_loose)),]
wood = wood[wood$X16S_per_ul >= 100,]

## Match up samples
samples = intersect(rownames(abundance), rownames(wood))
wood = wood[samples,]
abundance = abundance[samples,]

## ---- Parse taxonomy to extract Family level ----
tax_split = str_split(taxonomy$Taxon, "; ", simplify = TRUE)
families = tax_split[, 5]
families = gsub("^\\s+|\\s+$", "", families)
families[families == "none" | families == ""] = NA
names(families) = rownames(taxonomy)

cat("\nNumber of families extracted:", length(families), "\n")
cat("Number of non-NA families:", sum(!is.na(families)), "\n")

## ---- Aggregate by Family level ----
common_otus = intersect(colnames(abundance), names(families))
cat("Common OTUs:", length(common_otus), "\n")

abundance = abundance[, common_otus]
families_filtered = families[common_otus]

## Create family-level abundance matrix
unique_families = unique(families_filtered)
unique_families = unique_families[!is.na(unique_families) & unique_families != "" & unique_families != "none"]

cat("Number of unique families:", length(unique_families), "\n")

## Initialize as numeric matrix
family_abundance = matrix(0, nrow = nrow(abundance), 
                          ncol = length(unique_families))
rownames(family_abundance) = rownames(abundance)
colnames(family_abundance) = unique_families

## Sum OTU abundances by family
for(i in 1:ncol(abundance)) {
  otu_id = colnames(abundance)[i]
  fam = families_filtered[otu_id]
  if(!is.na(fam) && fam != "" && fam != "none" && fam %in% colnames(family_abundance)) {
    ## Ensure we're adding numeric values
    family_abundance[, fam] = family_abundance[, fam] + as.numeric(abundance[, i])
  }
}

## Convert to relative abundance
percent = family_abundance / rowSums(family_abundance)
percent = percent[, colSums(percent > 0) > 20]

cat("Testing", ncol(percent), "families\n\n")

## ---- Test association with mcrA ----
pvals_list = lapply(colnames(percent), function(family){
  meta_copy = wood
  meta_copy$family = scale(percent[, family])
  meta_copy$mcra_probe_loose = scale(meta_copy$mcra_probe_loose)
  
  fit0 = lmer(mcra_probe_loose ~ core_type + log16S + 
                (1 | seq_id), data = meta_copy)
  fit1 = lmer(mcra_probe_loose ~ family + core_type + log16S + 
                (1 | seq_id), data = meta_copy)
  p = anova(fit0, fit1)[2, "Pr(>Chisq)"]
  t = coef(summary(fit1))["family", "t value"]
  return(data.frame(family = family, p = p, t = t, stringsAsFactors = FALSE))
})

pvals = do.call("rbind", pvals_list)
pvals$FDR = p.adjust(pvals$p, "BH")

cat("=== Top 20 families associated with mcrA ===\n")
print(head(pvals[order(pvals$p),], 20))

## ---- Create heatmap ----
sig_families = pvals[!is.na(pvals$p) & pvals$p < 0.01, ]
cat("\nFound", nrow(sig_families), "significant families at p < 0.01\n")

if(nrow(sig_families) > 0) {
  sig_families = sig_families[order(sig_families$t),]
  rownames(sig_families) = sig_families$family
  sig_families$neg_log_p = -log10(sig_families$p)
  
  sig_families_percent = percent[rownames(wood[order(wood$mcra_probe_loose),]), 
                                 sig_families$family]
  
  wood$log10_mcra_qpcr = log10(1 + wood$mcra_probe_loose)
  annotation_col = wood[rownames(sig_families_percent), 
                        c("core_type", "log16S", "log10_mcra_qpcr")]
  colnames(annotation_col) = c("Type", "Log10_16S", "Log10_MCRA")
  
  annotation_row = data.frame(
    Associate_w_MCRA = ifelse(sig_families$p < 0.01, "p<0.01", "N.S.")
  )
  rownames(annotation_row) = sig_families$family
  
  ann_colors = list(
    Type = c(Inner = "#E41A89", Outer = "#00BFC4"),
    Associate_w_MCRA = c("p<0.01" = "#FF6B6B", "N.S." = "#90EE90")
  )
  
  dir.create("../../../outputs/figures", showWarnings = FALSE, recursive = TRUE)
  
  pdf("../../../outputs/figures/taxonomy_mcra_heatmap.pdf", height = 10, width = 16)
  pheatmap(t(sig_families_percent), 
           color = colorRampPalette(c("dodgerblue", "white", "red"))(100),
           scale = "row", 
           annotation_col = annotation_col, 
           annotation_row = annotation_row,
           annotation_colors = ann_colors,
           show_colnames = FALSE, 
           treeheight_row = 0, 
           cluster_cols = FALSE, 
           cluster_rows = FALSE,
           labels_row = sig_families$family,
           annotation_names_row = FALSE)
  dev.off()
  
  cat("\nHeatmap saved!\n")
}

write.csv(pvals, "../../../outputs/tables/family_mcra_associations.csv", row.names = FALSE)



pheatmap(t(sig_families_percent), 
         color = colorRampPalette(c("dodgerblue", "white", "red"))(100),
         scale = "row", 
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_colnames = FALSE, 
         treeheight_row = 0, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         labels_row = sig_families$family,
         annotation_names_row = FALSE)
























library(tidyverse)
library(lme4)
library(pheatmap)

## ---- Read in taxonomy and metadata ----
taxonomy = read.table('../../../data/raw/16s/taxonomy_table.txt',
                      sep = "\t", header = T, row.names = 1, quote = "",
                      stringsAsFactors = FALSE)

## Read abundance table - MAKE SURE IT'S NUMERIC
otu_table = read.table('../../../data/raw/16s/OTU_table.txt',
                       sep = "\t", header = T, row.names = 1, quote = "",
                       check.names = FALSE)

## Convert to numeric matrix
otu_matrix = as.matrix(otu_table)
otu_matrix = apply(otu_matrix, 2, as.numeric)
rownames(otu_matrix) = rownames(otu_table)

## Transpose so samples are rows
abundance = t(otu_matrix)

## Read metadata
meta = read.csv('../../../data/raw/picrust/16S_tree_sample_table_with_meta.csv',
                row.names = 1)
meta$log16S = log10(1 + meta$X16S_per_ul)
wood = meta[meta$material == "Wood" & !(is.na(meta$material)),]
wood = wood[!(is.na(wood$mcra_probe_loose)),]
wood = wood[wood$X16S_per_ul >= 100,]

## Match up samples
samples = intersect(rownames(abundance), rownames(wood))
wood = wood[samples,]
abundance = abundance[samples,]

cat("Total samples:", nrow(wood), "\n")

## ---- Parse taxonomy to extract Family level ----
tax_split = str_split(taxonomy$Taxon, "; ", simplify = TRUE)
families = tax_split[, 5]
families = gsub("^\\s+|\\s+$", "", families)
families[families == "none" | families == ""] = NA
names(families) = rownames(taxonomy)

## ---- Aggregate by Family level ----
common_otus = intersect(colnames(abundance), names(families))
abundance = abundance[, common_otus]
families_filtered = families[common_otus]

## Create family-level abundance matrix
unique_families = unique(families_filtered)
unique_families = unique_families[!is.na(unique_families) & unique_families != "" & unique_families != "none"]

family_abundance = matrix(0, nrow = nrow(abundance), 
                          ncol = length(unique_families))
rownames(family_abundance) = rownames(abundance)
colnames(family_abundance) = unique_families

## Sum OTU abundances by family
for(i in 1:ncol(abundance)) {
  otu_id = colnames(abundance)[i]
  fam = families_filtered[otu_id]
  if(!is.na(fam) && fam != "" && fam != "none" && fam %in% colnames(family_abundance)) {
    family_abundance[, fam] = family_abundance[, fam] + as.numeric(abundance[, i])
  }
}

## Convert to relative abundance
percent = family_abundance / rowSums(family_abundance)
percent = percent[, colSums(percent > 0) > 2]  ## present in more than 2 samples

cat("Testing", ncol(percent), "families\n\n")

## ---- Test association with mcrA ----
pvals_list = lapply(colnames(percent), function(family){
  meta_copy = wood
  meta_copy$family = scale(percent[, family])
  meta_copy$mcra_probe_loose = scale(meta_copy$mcra_probe_loose)
  
  fit0 = lmer(mcra_probe_loose ~ core_type + log16S + 
                (1 | seq_id), data = meta_copy)
  fit1 = lmer(mcra_probe_loose ~ family + core_type + log16S + 
                (1 | seq_id), data = meta_copy)
  p = anova(fit0, fit1)[2, "Pr(>Chisq)"]
  t = coef(summary(fit1))["family", "t value"]
  return(data.frame(family = family, p = p, t = t, stringsAsFactors = FALSE))
})

pvals = do.call("rbind", pvals_list)
pvals$FDR = p.adjust(pvals$p, "BH")

cat("=== Top 20 families associated with mcrA ===\n")
print(head(pvals[order(pvals$p),], 20))

## ---- Create heatmap ----
sig_families = pvals[!is.na(pvals$p) & pvals$p < 0.05, ]  ## Changed to p < 0.05
cat("\nFound", nrow(sig_families), "significant families at p < 0.05\n")

if(nrow(sig_families) > 0) {
  sig_families = sig_families[order(sig_families$t),]
  rownames(sig_families) = sig_families$family
  sig_families$neg_log_p = -log10(sig_families$p)
  
  sig_families_percent = percent[rownames(wood[order(wood$mcra_probe_loose),]), 
                                 sig_families$family]
  
  wood$log10_mcra_qpcr = log10(1 + wood$mcra_probe_loose)
  annotation_col = wood[rownames(sig_families_percent), 
                        c("core_type", "log16S", "log10_mcra_qpcr")]
  colnames(annotation_col) = c("Type", "Log10_16S", "Log10_MCRA")
  
  annotation_row = data.frame(
    Associate_w_MCRA = ifelse(sig_families$p < 0.01, "p<0.01", 
                              ifelse(sig_families$p < 0.05, "p<0.05", "N.S."))
  )
  rownames(annotation_row) = sig_families$family
  
  ann_colors = list(
    Type = c(Inner = "#E41A89", Outer = "#00BFC4"),
    Associate_w_MCRA = c("p<0.01" = "#FF0000", "p<0.05" = "#FFA500", "N.S." = "#90EE90")
  )
  
  dir.create("../../../outputs/figures", showWarnings = FALSE, recursive = TRUE)
  
  pdf("../../../outputs/figures/taxonomy_mcra_heatmap.pdf", height = 12, width = 16)
  pheatmap(t(sig_families_percent), 
           color = colorRampPalette(c("dodgerblue", "white", "red"))(100),
           scale = "row", 
           annotation_col = annotation_col, 
           annotation_row = annotation_row,
           annotation_colors = ann_colors,
           show_colnames = FALSE, 
           treeheight_row = 0, 
           cluster_cols = FALSE, 
           cluster_rows = FALSE,
           labels_row = sig_families$family,
           annotation_names_row = FALSE)
  dev.off()
  
  cat("\nHeatmap saved!\n")
}

write.csv(pvals, "../../../outputs/tables/family_mcra_associations.csv", row.names = FALSE)