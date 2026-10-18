# ==============================================================================
# 01_faprotax_dump_HW_SW.R -- FAPROTAX functions in heartwood and sapwood (Fig 6b)
# ------------------------------------------------------------------------------
# Stage D; run from the repository root. The FAPROTAX calculation of
# 08_figures/figS10_faprotax.R (same phyloseq object and rarefaction), written out for
# every function, and the heartwood vs sapwood tests behind the asterisks of Figure 6b.
# Outputs: outputs/data/FAPROTAX_all_functions_HW_SW.csv  (pooled shares, every function)
#          outputs/data/FAPROTAX_HW_SW_tests.csv          (Fig 6b functions, per-library tests)
# ==============================================================================
suppressPackageStartupMessages({library(phyloseq);library(microeco);library(tidyverse)})

source("code/lib/build_phyloseq.R")
.ps16   <- build_16s_phyloseq()
no_mito <- .ps16$ps
taxa_names(no_mito)<-paste0("ASV",seq(ntaxa(no_mito))); ps.rare<-rarefy_16s(no_mito)
ps.ra<-transform_sample_counts(ps.rare,function(x) x/sum(x)*100); colnames(tax_table(ps.ra))<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
ps.wood<-subset_samples(ps.ra,material=="Wood")
for(i in 1:7) tax_table(ps.wood)[,i]<-paste0(c("k__","p__","c__","o__","f__","g__","s__")[i],tax_table(ps.wood)[,i])
meco_all<-microtable$new(otu_table=as.data.frame(otu_table(ps.wood)),tax_table=noquote(as.data.frame(tax_table(ps.wood))),sample_table=as.data.frame(as.matrix(sample_data(ps.wood))),phylo_tree=phy_tree(ps.wood))
meco_core<-clone(meco_all)$merge_samples("core_type")
t<-trans_func$new(meco_core); t$cal_spe_func(prok_database="FAPROTAX"); t$cal_spe_func_perc(abundance_weighted=TRUE,dec=2)
pc<-t$res_spe_func_perc
df<-data.frame(func=colnames(pc), HW=as.numeric(pc["Inner",]), SW=as.numeric(pc["Outer",]))
df$lr<-ifelse(df$SW>0,log2(df$HW/df$SW),NA)

# The same shares with methanogen ASVs removed. FAPROTAX assigns functions by genus
# name, and several categories recount the methanogens: "dark_hydrogen_oxidation"
# includes Methanobacteriaceae and Methanosarcinaceae, and "methylotrophy" includes
# the methyl-reducing Methanomassiliicoccaceae. Without them neither is enriched in
# heartwood (2026-10-02). Figure 6b uses these columns, as panel (c) removes methanogen
# contributions from the PICRUSt2 pathways.
source("code/lib/methanogen_families.R")
fm <- as.matrix(t$res_spe_func); ab <- as.matrix(t$otu_table)[rownames(fm), c("Inner", "Outer")]
fam <- sub("^f__", "", as.data.frame(t$tax_table)[rownames(fm), "Family"])
is_mg <- fam %in% METHANOGEN_FAMILIES
share <- function(keep) 100 * colSums(ab[keep, , drop = FALSE]) / colSums(as.matrix(t$otu_table)[, c("Inner", "Outer")])
chk <- t(sapply(df$func, function(f) share(fm[, f] == 1)))
stopifnot(max(abs(chk[, "Inner"] - df$HW), abs(chk[, "Outer"] - df$SW)) < 0.01)   # reproduces microeco's shares
nm <- t(sapply(df$func, function(f) share(fm[, f] == 1 & !is_mg)))
df$HW_nonmethanogen <- round(nm[, "Inner"], 2); df$SW_nonmethanogen <- round(nm[, "Outer"], 2)
df$lr_nonmethanogen <- ifelse(df$SW_nonmethanogen > 0, log2(df$HW_nonmethanogen / df$SW_nonmethanogen), NA)
df$methanogen_asvs <- sapply(df$func, function(f) sum(fm[, f] == 1 & is_mg))
df<-df[order(-df$HW),]
write.csv(df,"outputs/data/FAPROTAX_all_functions_HW_SW.csv",row.names=FALSE)
cat("=== ALL FAPROTAX functions (HW=Inner, SW=Outer, % relative abundance) ===\n")
for(i in seq_len(nrow(df))) cat(sprintf("%-52s HW=%6.2f SW=%6.2f lr=%s | without methanogens HW=%6.2f SW=%6.2f (%d methanogen ASVs)\n",df$func[i],df$HW[i],df$SW[i],ifelse(is.na(df$lr[i]),"  Inf",sprintf("%+.2f",df$lr[i])),df$HW_nonmethanogen[i],df$SW_nonmethanogen[i],df$methanogen_asvs[i]))

# ---- Per-library shares and tests for the functions of Figure 6b --------------------
# Each library's share of a function (methanogen ASVs removed except for the
# methanogenesis row, as in the panel). Pooled shares are the mean over libraries (each
# library sums to 100%), so the log2 ratio here is the unrounded one the panel plots.
# Paired test: heartwood vs sapwood of the same tree, replicate libraries of a tree and
# tissue averaged; Wilcoxon signed-rank, BH across the plotted functions. The unpaired
# test uses every library.
source("code/lib/fig6b_functions.R")
sel <- FIG6B_FUNCTIONS$func; stopifnot(all(sel %in% colnames(fm)))
ab_lib <- as.matrix(meco_all$otu_table)
stopifnot(all(rownames(fm) %in% rownames(ab_lib)),
          sum(ab_lib[setdiff(rownames(ab_lib), rownames(fm)), ]) == 0,   # ASVs outside fm are absent
          all(abs(colSums(ab_lib) - 100) < 1e-6))
ab_lib <- ab_lib[rownames(fm), , drop = FALSE]
lib_share <- sapply(sel, function(f) {
  keep <- fm[, f] == 1 & (f == "hydrogenotrophic_methanogenesis" | !is_mg)
  round(colSums(ab_lib[keep, , drop = FALSE]), 8) })   # shares are multiples of 1/35 %; rounding makes ties exact
smp <- as.data.frame(as.matrix(meco_all$sample_table), stringsAsFactors = FALSE)
L <- data.frame(tree_id = smp[rownames(lib_share), "tree_id"], core_type = smp[rownames(lib_share), "core_type"],
                lib_share, check.names = FALSE)
stopifnot(all(L$core_type %in% c("Inner", "Outer")))
Lt <- L %>% filter(!is.na(tree_id)) %>% group_by(tree_id, core_type) %>%
  summarise(across(all_of(sel), mean), .groups = "drop")
tests <- map_dfr(sel, function(f) {
  hw <- L[[f]][L$core_type == "Inner"]; sw <- L[[f]][L$core_type == "Outer"]
  w <- Lt %>% select(tree_id, core_type, v = all_of(f)) %>%
    pivot_wider(names_from = core_type, values_from = v) %>% filter(!is.na(Inner), !is.na(Outer))
  d <- w$Inner - w$Outer
  tibble(func = f, pooled_HW = mean(hw), pooled_SW = mean(sw), log2_ratio = log2(mean(hw) / mean(sw)),
         median_HW = median(hw), median_SW = median(sw), n_HW = length(hw), n_SW = length(sw),
         n_pairs = nrow(w), pairs_HW_gt = sum(d > 0), pairs_SW_gt = sum(d < 0), pairs_tied = sum(d == 0),
         paired_p = wilcox.test(w$Inner, w$Outer, paired = TRUE, exact = FALSE)$p.value,
         unpaired_p = wilcox.test(hw, sw, exact = FALSE)$p.value) }) %>%
  mutate(paired_FDR = p.adjust(paired_p, "BH"), unpaired_FDR = p.adjust(unpaired_p, "BH"))
# pooled shares agree with the table above (2-decimal columns)
pooled_chk <- df %>% filter(func %in% sel) %>%
  transmute(func, HWp = ifelse(func == "hydrogenotrophic_methanogenesis", HW, HW_nonmethanogen),
            SWp = ifelse(func == "hydrogenotrophic_methanogenesis", SW, SW_nonmethanogen)) %>%
  inner_join(tests, by = "func")
stopifnot(nrow(pooled_chk) == length(sel),
          max(abs(pooled_chk$HWp - pooled_chk$pooled_HW), abs(pooled_chk$SWp - pooled_chk$pooled_SW)) < 0.01)
write.csv(tests, "outputs/data/FAPROTAX_HW_SW_tests.csv", row.names = FALSE)
cat("\n=== Figure 6b functions: heartwood vs sapwood (", nrow(L), " libraries; ",
    sum(is.na(L$tree_id)), " without a tree) ===\n", sep = "")
print(as.data.frame(tests %>% mutate(across(where(is.double), ~ signif(.x, 3)))), row.names = FALSE)
