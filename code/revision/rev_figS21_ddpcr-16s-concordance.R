#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS21_ddpcr-16s-concordance.R  (SI)
# ddPCR functional-gene detection vs 16S amplicon-inferred function. ddPCR is far
# more sensitive, so the DIRECTIONAL prediction is: many samples ddPCR+ / 16S=0
# (ddPCR catches what 16S misses), almost none ddPCR=0 / 16S+ (16S should not see
# a function the more-sensitive assay cannot). Non-detects are shown (log axes
# with a labelled "0" floor), coloured by concordance; the black line is the
# quantitative fit where both detect. Compared to TOTAL 16S methanogenesis /
# methanotrophy. 16S functional data from sibling tree-microbiome repo (path below).
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
m<-merge(o[,c("key","material","mcra","pmoa","mmox")],f[,c("key","methanogenesis","methanotrophy")],by="key")
cols<-c("both detect"="#2c7fb8","ddPCR+ / 16S = 0"="#d95f02","ddPCR = 0 / 16S+"="#7570b3","both = 0"="grey65")
set.seed(1)
panel<-function(gene,func,glab,flab){
  s<-m[is.finite(m[[gene]])&is.finite(m[[func]]),]; s$g<-s[[gene]]; s$f<-s[[func]]
  s$cat<-with(s,ifelse(g>0&f>0,"both detect",ifelse(g>0&f==0,"ddPCR+ / 16S = 0",
             ifelse(g==0&f>0,"ddPCR = 0 / 16S+","both = 0"))))
  s$cat<-factor(s$cat,levels=names(cols))   # identical levels in every panel -> one merged legend
  fx<-function(v){p<-v[v>0]; fl<-10^(floor(log10(min(p)))-1); v2<-ifelse(v==0, fl*10^runif(length(v),-.2,.2), v); list(v=v2,fl=fl,top=max(p))}
  gg<-fx(s$g); ff<-fx(s$f); s$gp<-gg$v; s$fp<-ff$v
  br<-function(fl,top){b<-10^(seq(floor(log10(fl)),ceiling(log10(top)))); l<-ifelse(b==fl,"0",parse(text=paste0("10^",round(log10(b))))); list(b=b,l=l)}
  bg<-br(gg$fl,gg$top); bf<-br(ff$fl,ff$top)
  pos<-s[s$g>0&s$f>0,]; fit<-if(nrow(pos)>8){cf<-coef(lm(log10(g)~log10(f),data=pos)); data.frame(f=10^seq(log10(min(pos$f)),log10(max(pos$f)),length=50))}else NULL
  if(!is.null(fit)) fit$g<-10^(cf[1]+cf[2]*log10(fit$f))
  tab<-table(factor(s$cat,levels=names(cols)))
  lab<-sprintf("both detect: %d (r=%.2f)\nddPCR+ / 16S=0: %d\nddPCR=0 / 16S+: %d\nddPCR %.0f%% vs 16S %.0f%%",
     tab["both detect"], if(nrow(pos)>2)cor(log10(pos$g),log10(pos$f)) else NA,
     tab["ddPCR+ / 16S = 0"], tab["ddPCR = 0 / 16S+"], 100*mean(s$g>0),100*mean(s$f>0))
  p<-ggplot(s,aes(fp,gp,colour=cat))+
    annotate("rect",xmin=ff$fl/3,xmax=ff$fl*3,ymin=-Inf,ymax=Inf,fill="grey92",alpha=.6)+
    annotate("rect",ymin=gg$fl/3,ymax=gg$fl*3,xmin=-Inf,xmax=Inf,fill="grey92",alpha=.6)+
    geom_point(alpha=.6,size=1.5)+scale_colour_manual(values=cols,name=NULL,drop=FALSE)+
    scale_x_log10(breaks=bf$b,labels=bf$l)+scale_y_log10(breaks=bg$b,labels=bg$l)+
    annotate("text",ff$fl,gg$top,label=lab,hjust=0,vjust=1,size=2.8)+
    labs(title=glab,x=flab,y=bquote(.(gene)~"ddPCR load (copies, 0=nd)"))+
    theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="bottom")+
    guides(colour=guide_legend(override.aes=list(size=3,alpha=1),nrow=1))
  if(!is.null(fit)) p<-p+geom_line(data=fit,aes(f,g),colour="black",linewidth=.6,inherit.aes=FALSE)
  p
}
A<-panel("mcra","methanogenesis","A  mcrA vs 16S methanogenesis","16S methanogenesis (%, 0=nd)")
B<-panel("pmoa","methanotrophy","B  pmoA vs 16S methanotrophy","16S methanotrophy (%, 0=nd)")
C<-panel("mmox","methanotrophy","C  mmoX vs 16S methanotrophy","16S methanotrophy (%, 0=nd)")
ggsave("outputs/revision/figS21_ddpcr_16s_concordance.png",
  (A|B|C)+plot_layout(guides="collect")+
  plot_annotation(title="ddPCR detects methane-cycling genes where 16S is blind, and almost never the reverse",
    subtitle="non-detects shown at '0'; ddPCR far more sensitive (directional: many ddPCR+/16S=0, ~0 ddPCR=0/16S+)",
    theme=theme(plot.title=element_text(face="bold",size=13),plot.subtitle=element_text(size=9.5),legend.position="bottom")),
  width=14,height=5.4,dpi=200)
cat("wrote outputs/revision/figS21_ddpcr_16s_concordance.png\n")
