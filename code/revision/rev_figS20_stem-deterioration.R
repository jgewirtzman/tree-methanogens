#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS20_stem-deterioration.R  (SI)
# Stem deterioration vs CH4 flux (2023 cross-species survey, species-controlled).
# A PCA of four stem-condition scores (bark loss, wounding, moss/lichen, fungus)
# separates two orthogonal, mechanistically distinct axes:
#   PC1 = overall decay  -> INVERTED-U (early decay peaks; moisture/labile C/anoxia
#         before structural loss + O2 ingress) — quadratic p=0.002
#   PC2 = wounding        -> MONOTONIC increase (physical gas-escape) — linear p=0.021
# (C) the same hump on the intuitive bark-condition scale. Dead trees are too few
# and campaign-inconsistent to resolve the extreme; the signal is early/moderate.
# Reads data/; writes outputs/revision/figS20_stem_deterioration.png
# ==============================================================================
suppressMessages({library(ggplot2);library(patchwork);library(lme4);library(lmerTest)})
options(warn=-1, stringsAsFactors=FALSE)
asinh10<-function(x) asinh(x/0.1)/log(10)
ord<-function(v){v<-tolower(trimws(as.character(v)));x<-suppressWarnings(as.numeric(v));x[v=="dead"]<-4;x}
y<-read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv",check.names=FALSE)
names(y)<-make.names(names(y)); y<-y[!is.na(y$CH4_best.flux),]
y$fx<-asinh10(y$CH4_best.flux); y$sp<-as.factor(y$Species)
y$bark<-ord(y$Bark.Missing);y$wound<-ord(y$Wounding.Holes);y$moss<-ord(y$Moss.Lichen.Cover);y$fungus<-ord(y$Visible.fungus)
d<-y[complete.cases(y[,c("bark","wound","moss","fungus")]),]
pc<-prcomp(d[,c("bark","wound","moss","fungus")],scale.=TRUE)
d$PC1<-pc$x[,1]; d$PC2<-pc$x[,2]
ve<-round(100*summary(pc)$importance[2,1:2])

qfit<-function(v,quad){f<-if(quad)paste0("fx~",v,"+I(",v,"^2)+(1|sp)") else paste0("fx~",v,"+(1|sp)")
  m<-lmer(as.formula(f),data=d);co<-summary(m)$coefficients
  g<-data.frame(x=seq(min(d[[v]]),max(d[[v]]),length=60))
  g$y<-if(quad) co[1,1]+co[2,1]*g$x+co[3,1]*g$x^2 else co[1,1]+co[2,1]*g$x
  list(g=g, p=co[nrow(co),"Pr(>|t|)"])}
th<-theme_bw(base_size=11)
A<-qfit("PC1",TRUE); B<-qfit("PC2",FALSE)

pA<-ggplot(d,aes(PC1,fx))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
  geom_jitter(width=.05,height=0,alpha=.35,size=1.4,colour="#8c510a")+
  geom_line(data=A$g,aes(x,y),linewidth=1)+
  annotate("text",-Inf,Inf,hjust=-.08,vjust=1.4,size=3.2,label=sprintf("quadratic p=%.3f",A$p))+
  labs(title=sprintf("A  Overall decay axis (PC1, %d%%) — inverted-U",ve[1]),
       x="PC1  (bark loss + fungus + moss + wounding)",y=expression(CH[4]~flux~(arcsinh)))+th
pB<-ggplot(d,aes(PC2,fx))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
  geom_jitter(width=.05,height=0,alpha=.35,size=1.4,colour="#8c510a")+
  geom_line(data=B$g,aes(x,y),linewidth=1)+
  annotate("text",-Inf,Inf,hjust=-.08,vjust=1.4,size=3.2,label=sprintf("linear p=%.3f",B$p))+
  labs(title=sprintf("B  Wounding axis (PC2, %d%%) — monotonic",ve[2]),
       x="PC2  (wounding vs. surface colonization)",y=NULL)+th
d$cond<-cut(d$bark,c(0,1,2,4),labels=c("healthy","moderate","severe/dead"))
nD<-as.data.frame(table(d$cond));names(nD)<-c("cond","n")
pC<-ggplot(d,aes(cond,fx,fill=cond))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
  geom_boxplot(width=.55,alpha=.85,outlier.size=.6)+
  geom_text(data=nD,aes(cond,-2.1,label=paste0("n=",n)),inherit.aes=FALSE,size=3)+
  scale_fill_manual(values=c(healthy="#7fbf7b",moderate="#f0a202","severe/dead"="#8c510a"),guide="none")+
  labs(title="C  Bark condition (intuitive hump)",x=NULL,y=NULL)+th
ggsave("outputs/revision/figS20_stem_deterioration.png", pA+pB+pC+plot_layout(widths=c(1,1,.9)),
       width=14.5,height=4.6,dpi=300)
cat(sprintf("PC1 quad p=%.3f  PC2 linear p=%.3f\n",A$p,B$p))
cat("wrote outputs/revision/figS20_stem_deterioration.png\n")
