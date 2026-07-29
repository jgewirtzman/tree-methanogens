# ==============================================================================
# REVISION helper: dump FAPROTAX HW/SW (Inner/Outer) for ALL functions.
# Replicates the phyloseq build + FAPROTAX calc of code/06_figures/08d verbatim,
# then writes the full per-function heartwood/sapwood table used by panel (b) of
# the hydrogenotrophy synthesis figure. NEW file (does not modify 08d).
# Output: outputs/revision/FAPROTAX_all_functions_HW_SW.csv
# ==============================================================================
suppressPackageStartupMessages({library(phyloseq);library(microeco);library(tidyverse)})
ddpcr <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv")
water <- read.csv("data/raw/tree_cores/Tree_Core_Sectioning_Data.csv")
qpcr_16s <- read.csv("data/raw/16s/16s_w_metadata.csv"); qpcr_16s<-subset(qpcr_16s,qpcr_16s$Sample.ID!="None"); qpcr_16s<-qpcr_16s[,c(3,4,6)]
water$seq_id<-toupper(water$seq_id)
ddpcr<-merge(ddpcr,water,by=c("core_type","seq_id"),all.x=TRUE); ddpcr<-merge(ddpcr,qpcr_16s,by=c("core_type","seq_id"),all.x=TRUE)
otu_tab<-read.delim("data/raw/16s/OTU_table.txt",header=TRUE,row.names=1)
bastard_tax<-otu_tab[,590:596]; bastard_tax[bastard_tax==""]<-NA
tax_tab_pre<-tax_table(bastard_tax); taxa_names(tax_tab_pre)<-sub("sp","seq",taxa_names(tax_tab_pre))
otu_tab_corr<-otu_tab[,1:589]; otu_table_pre<-otu_table(otu_tab_corr,taxa_are_rows=TRUE)
phylo_tree<-read_tree("data/raw/16s/unrooted_tree.nwk")
samp_data<-read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt",row.names=1); samp_data$RowName<-row.names(samp_data)
samp_data$seq_id<-sub("prime","'",samp_data$seq_id); samp_data$seq_id<-sub("star","*",samp_data$seq_id); samp_data$seq_id<-sub("HM","H",samp_data$seq_id)
samp_data$seq_id[samp_data$ForwardFastqFile=="206_B01_16S_S3_R1_001.fastq"]<-"RO104"; samp_data$core_type[samp_data$ForwardFastqFile=="206_B01_16S_S3_R1_001.fastq"]<-"Inner"
samp_data_merged<-merge(ddpcr,samp_data,by=c("seq_id","core_type"),all.y=TRUE)
dups<-which(duplicated(samp_data_merged$RowName)==TRUE); samp_data_merged<-samp_data_merged[-c(dups),]; row.names(samp_data_merged)<-samp_data_merged$RowName
raw_ps<-phyloseq(tax_tab_pre,otu_table_pre,phylo_tree,sample_data(samp_data_merged))
pop_taxa<-function(physeq,badTaxa){allTaxa<-taxa_names(physeq);allTaxa<-allTaxa[!(allTaxa%in%badTaxa)];return(prune_taxa(allTaxa,physeq))}
mitochondria<-rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[,5]=="Mitochondria")]; chloroplast<-rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[,4]=="Chloroplast")]
no_mito<-pop_taxa(raw_ps,c(mitochondria,chloroplast))
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
df<-df[order(-df$HW),]
write.csv(df,"outputs/revision/FAPROTAX_all_functions_HW_SW.csv",row.names=FALSE)
cat("=== ALL FAPROTAX functions (HW=Inner, SW=Outer, % relative abundance) ===\n")
for(i in seq_len(nrow(df))) cat(sprintf("%-52s HW=%6.2f SW=%6.2f lr=%s\n",df$func[i],df$HW[i],df$SW[i],ifelse(is.na(df$lr[i]),"  Inf",sprintf("%+.2f",df$lr[i]))))
