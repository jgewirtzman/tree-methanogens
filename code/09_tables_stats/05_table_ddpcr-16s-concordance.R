#!/usr/bin/env Rscript
# ==============================================================================
# 05_table_ddpcr-16s-concordance.R  (candidate SI TABLE replacing the S05 figure)
# Same data/definitions as rev_figS21 (taxonomy-based, Fig 5 defs), summarized as
# a concordance table: per gene x compartment, the 4 detection categories, ddPCR
# vs 16S detection rates, and the quantitative fit (r, p) where both detect.
# Reads data/tables/ via code/lib/tables_read.R: asv_table("main") (ASV counts + taxonomy),
# dna_samples("2021_survey") (16S library -> sample, material; ddPCR copies per uL extract);
# joined by library_id. Samples: 2021 DNA samples whose 16S library has at least RAREFY_DEPTH
# (3500) non-organelle reads, the libraries the main 16S analysis keeps (build_phyloseq.R).
# Writes outputs/data/tbl_ddpcr_16s_concordance.csv + .png
# ==============================================================================
suppressMessages({library(gridExtra);library(grid)}); options(warn=-1,stringsAsFactors=FALSE)
source("code/lib/tables_read.R"); source("code/lib/build_phyloseq.R")   # RAREFY_DEPTH
# denominator: the library's non-organelle reads (asv_reads(): mitochondria and
# chloroplast ASVs removed, as in the main 16S analysis)
A<-asv_reads("main"); cnt<-A$cnt; tax<-A$tax; tot<-A$tot
deep<-names(tot)[tot>=RAREFY_DEPTH]
source("code/lib/methanogen_families.R"); mgfam<-METHANOGEN_FAMILIES
mg_asv<-rownames(cnt)[tax$Family %in% mgfam]
source("code/lib/load_methanotroph_definitions.R")
dp<-"code/lib/methanotroph_definitions.csv"
tax$mt<-classify_methanotrophs(tax,load_methanotroph_defs(dp),include_conditional=TRUE)
mt_asv<-rownames(cnt)[tax$mt %in% c("Known","Putative")]
s16<-data.frame(library_id=colnames(cnt),
                mg16=100*colSums(cnt[mg_asv,,drop=FALSE],na.rm=TRUE)/tot,
                mt16=100*colSums(cnt[mt_asv,,drop=FALSE],na.rm=TRUE)/tot)
ds<-dna_samples("2021_survey"); ds<-ds[ds$library_id %in% deep,]
o<-data.frame(library_id=ds$library_id,material=ds$material,mcra=ds$mcra_probe_ul,pmoa=ds$pmoa_ul,mmox=ds$mmox_ul)
m<-merge(o,s16,by="library_id")
row1<-function(gene,func,glab,mat){
  s<-m[m$material==mat & is.finite(m[[gene]])&is.finite(m[[func]]),]; g<-s[[gene]]; f<-s[[func]]
  pos<-g>0&f>0; rp<-c(NA,NA); if(sum(pos)>=4){ct<-suppressWarnings(cor.test(log10(g[pos]),log10(f[pos]))); rp<-c(ct$estimate,ct$p.value)}
  data.frame(Gene=glab,Compartment=mat,n=nrow(s),
    `Both detect`=sum(g>0&f>0),`ddPCR+ / 16S=0`=sum(g>0&f==0),`ddPCR=0 / 16S+`=sum(g==0&f>0),`Both=0`=sum(g==0&f==0),
    `ddPCR %`=round(100*mean(g>0)),`16S %`=round(100*mean(f>0)),
    r=ifelse(is.na(rp[1]),"-",sprintf("%.2f",rp[1])),
    p=ifelse(is.na(rp[2]),"-",ifelse(rp[2]<0.001,"<0.001",sprintf("%.3f",rp[2]))),
    check.names=FALSE)
}
combos<-list(c("mcra","mg16","mcrA"),c("pmoa","mt16","pmoA"),c("mmox","mt16","mmoX"))
tab<-do.call(rbind,lapply(c("Wood","Soil"),function(mat) do.call(rbind,lapply(combos,function(x) row1(x[1],x[2],x[3],mat)))))
tab<-tab[order(tab$Gene,tab$Compartment),]
write.csv(tab,"outputs/data/tbl_ddpcr_16s_concordance.csv",row.names=FALSE)
th<-ttheme_minimal(base_size=10,core=list(fg_params=list(hjust=0.5,x=0.5)),
  colhead=list(fg_params=list(fontface="bold")))
g<-tableGrob(tab,rows=NULL,theme=th)
ggplot2::ggsave("outputs/figures/generated/tbl_ddpcr_16s_concordance.png",g,width=11,height=2.6,dpi=200,bg="white")
cat("wrote outputs/tbl_ddpcr_16s_concordance.{csv,png}\n"); print(tab,row.names=FALSE)
