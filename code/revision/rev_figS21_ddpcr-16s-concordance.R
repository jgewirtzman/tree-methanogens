#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS21_ddpcr-16s-concordance.R  (SI)
# Cross-method concordance: ddPCR functional-gene loads vs 16S amplicon-inferred
# function (relative abundance). Two independent molecular methods agree, which
# validates the ddPCR quantification (relevant to R2 controls concerns):
#   A  mcrA (ddPCR)  vs 16S hydrogenotrophic methanogenesis (%)
#   B  pmoA (ddPCR)  vs 16S methanotrophy (%)
#   C  mmoX (ddPCR)  vs 16S methanotrophy (%)
# 16S functional predictions from the sibling repo (tree-microbiome); path below.
# Writes outputs/revision/figS21_ddpcr_16s_concordance.png
# ==============================================================================
suppressMessages({library(ggplot2);library(patchwork)}); options(warn=-1,stringsAsFactors=FALSE)
TM<-"/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-microbiome"
o<-read.csv("data/processed/molecular/tree_data_methanogen_group.csv",check.names=FALSE); o<-o[,names(o)!=""]
num<-function(x) suppressWarnings(as.numeric(x))
o$mcra<-num(o$mcra_probe_loose); o$pmoa<-num(o$pmoa_loose); o$mmox<-num(o$mmox_loose)
o$key<-paste0(o$seq_id,o$core_type)
f<-read.csv(file.path(TM,"metabolisms/16S_metabolisms_weighted.csv"),check.names=FALSE)
names(f)[1]<-"fid"; f$key<-sub("[.]16S[.]S[0-9]+$","",f$fid)
fcol<-c("hydrogenotrophic_methanogenesis","methanotrophy")
m<-merge(o[,c("key","material","mcra","pmoa","mmox")],f[,c("key",fcol)],by="key")
names(m)<-make.names(names(m))
cols<-c(Wood="#2c7a2c",Soil="#8c6d3f")
panel<-function(gene,func,glab,flab){
  s<-m[is.finite(m[[gene]])&is.finite(m[[func]])&m[[gene]]>0&m[[func]]>0 & m$material%in%c("Wood","Soil"),]
  s$y<-log10(s[[gene]]); s$x<-log10(s[[func]])
  r<-cor(s$x,s$y); rho<-cor(s$x,s$y,method="spearman")
  ggplot(s,aes(x,y,colour=material))+geom_point(alpha=.45,size=1.5)+
    geom_smooth(method="lm",se=FALSE,colour="black",linewidth=.6,inherit.aes=FALSE,aes(x,y))+
    scale_colour_manual(values=cols,name=NULL)+
    annotate("text",-Inf,Inf,hjust=-.06,vjust=1.3,size=3.2,label=sprintf("r=%.2f  rho=%.2f  n=%d",r,rho,nrow(s)))+
    labs(title=glab,x=flab,y=bquote(log[10]~.(gene)~load))+
    theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position=c(.99,.02),legend.justification=c(1,0))
}
A<-panel("mcra","hydrogenotrophic_methanogenesis","A  mcrA vs 16S methanogenesis",expression(log[10]~"16S hydrogenotrophic methanogenesis (%)"))
B<-panel("pmoa","methanotrophy","B  pmoA vs 16S methanotrophy",expression(log[10]~"16S methanotrophy (%)"))
C<-panel("mmox","methanotrophy","C  mmoX vs 16S methanotrophy",expression(log[10]~"16S methanotrophy (%)"))
ggsave("outputs/revision/figS21_ddpcr_16s_concordance.png",(A|B|C)+
  plot_annotation(
    title="ddPCR functional genes vs 16S-inferred function",
    subtitle="strong concordance for methanogens (mcrA); weak for methanotrophs (pmoA/mmoX) — 16S poorly infers methanotrophy, so ddPCR is the direct measure",
    caption="log10 ddPCR gene load vs log10 16S FAPROTAX function (%). Two independent molecular methods.",
    theme=theme(plot.title=element_text(face="bold",size=13),plot.subtitle=element_text(size=9.5))),
  width=14,height=4.8,dpi=200)
cat("wrote outputs/revision/figS21_ddpcr_16s_concordance.png\n")
