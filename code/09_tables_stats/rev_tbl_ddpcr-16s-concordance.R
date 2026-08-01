#!/usr/bin/env Rscript
# ==============================================================================
# rev_tbl_ddpcr-16s-concordance.R  (candidate SI TABLE replacing the S05 figure)
# Same data/definitions as rev_figS21 (taxonomy-based, Fig 5 defs), summarized as
# a concordance table: per gene x compartment, the 4 detection categories, ddPCR
# vs 16S detection rates, and the quantitative fit (r, p) where both detect.
# Writes outputs/revision/tbl_ddpcr_16s_concordance.csv + .png
# ==============================================================================
suppressMessages({library(gridExtra);library(grid)}); options(warn=-1,stringsAsFactors=FALSE)
num<-function(x) suppressWarnings(as.numeric(x))
otu<-read.delim("data/raw/16s/OTU_table.txt",header=TRUE,row.names=1,check.names=FALSE)
taxc<-intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),names(otu))
tax<-otu[,taxc]; cnt<-as.matrix(sapply(otu[,setdiff(names(otu),taxc),drop=FALSE],num)); rownames(cnt)<-rownames(otu); tot<-colSums(cnt,na.rm=TRUE)
mgfam<-c("Methanobacteriaceae","Methanomassiliicoccaceae","Methanoregulaceae","Methanocellaceae",
         "Methanosaetaceae","Methanomicrobiaceae","Methanosarcinaceae","Methanomethyliaceae","Methanocorpusculaceae")
mg_asv<-rownames(otu)[tax$Family %in% mgfam]
source("code/lib/load_methanotroph_definitions.R")
dp<-if(file.exists("data/processed/molecular/methanotroph_definitions_revised.csv"))"data/processed/molecular/methanotroph_definitions_revised.csv" else "data/compiled/methanotroph_definitions.csv"
tax$mt<-classify_methanotrophs(tax,load_methanotroph_defs(dp),include_conditional=TRUE)
mt_asv<-rownames(otu)[tax$mt %in% c("Known","Putative")]
s16<-data.frame(key=sub("[.]16S[.]S[0-9]*$","",colnames(cnt)),
                mg16=100*colSums(cnt[mg_asv,,drop=FALSE],na.rm=TRUE)/tot,
                mt16=100*colSums(cnt[mt_asv,,drop=FALSE],na.rm=TRUE)/tot)
o<-read.csv("data/processed/molecular/tree_data_methanogen_group.csv",check.names=FALSE); o<-o[,names(o)!=""]
o$mcra<-num(o$mcra_probe_loose); o$pmoa<-num(o$pmoa_loose); o$mmox<-num(o$mmox_loose); o$key<-paste0(o$seq_id,o$core_type)
m<-merge(o[,c("key","material","mcra","pmoa","mmox")],s16,by="key")
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
ggplot2::ggsave("outputs/revision/tbl_ddpcr_16s_concordance.png",g,width=11,height=2.6,dpi=200,bg="white")
cat("wrote outputs/revision/tbl_ddpcr_16s_concordance.{csv,png}\n"); print(tab,row.names=FALSE)
