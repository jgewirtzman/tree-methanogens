# ==============================================================================
# REVISION helper: dump FAPROTAX HW/SW (Inner/Outer) for ALL functions.
# Replicates the phyloseq build + FAPROTAX calc of code/08_figures/08d verbatim,
# then writes the full per-function heartwood/sapwood table used by panel (b) of
# the hydrogenotrophy synthesis figure. NEW file (does not modify 08d).
# Output: outputs/data/FAPROTAX_all_functions_HW_SW.csv
# ==============================================================================
suppressPackageStartupMessages({library(phyloseq);library(microeco);library(tidyverse)})

# 16S phyloseq assembly: single definition in code/lib/build_phyloseq.R.
# All nine copies of this ~45-line block produced an identical object
# (71,511 taxa / 587 samples / 6,427,747 counts) before consolidation.
source("code/lib/build_phyloseq.R")
.ps16   <- build_16s_phyloseq()
no_mito <- .ps16$ps
phylo_tree <- .ps16$tree
taxa_names(no_mito)<-paste0("ASV",seq(ntaxa(no_mito))); set.seed(46814); ps.rare<-rarefy_even_depth(no_mito,sample.size=3500)
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
