#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS20_stem-deterioration.R  (SI)
# Stem condition vs CH4 flux (2023 cross-species survey, species-controlled).
# Two stem-condition metrics, two distinct patterns:
#   (A) BARK LOSS (decay proxy) -> emission peaks at moderate loss and falls at
#       severe/dead: early decay (moisture, labile C, anoxia) before structural
#       loss + O2 ingress. Moderate vs healthy: +0.13, p<0.001 (mixed model).
#   (B) WOUNDING (physical damage) -> monotonic increase (gas-escape paths):
#       linear slope +0.05, p=0.04 (mixed model).
# Both from arcsinh(flux) ~ metric + (1|species). Dead trees are too few and
# campaign-inconsistent to resolve the severe extreme; the robust signal is the
# early/moderate range. Reads data/; writes outputs/revision/.
# ==============================================================================
suppressMessages({library(ggplot2);library(patchwork);library(lme4);library(lmerTest)})
options(warn=-1, stringsAsFactors=FALSE)
asinh10<-function(x) asinh(x/0.1)/log(10)
ord<-function(v){v<-tolower(trimws(as.character(v)));x<-suppressWarnings(as.numeric(v));x[v=="dead"]<-4;x}
y<-read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv",check.names=FALSE)
names(y)<-make.names(names(y)); y<-y[!is.na(y$CH4_best.flux),]
y$fx<-asinh10(y$CH4_best.flux); y$sp<-as.factor(y$Species)
y$bark<-ord(y$Bark.Missing); y$wound<-ord(y$Wounding.Holes)
th<-theme_bw(base_size=11)

# (A) bark loss -> hump, as 3-level condition
db<-y[!is.na(y$bark),]; db$cond<-cut(db$bark,c(0,1,2,4),labels=c("healthy","moderate","severe/dead"))
pbark<-summary(lmer(fx~cond+(1|sp),data=db))$coefficients["condmoderate","Pr(>|t|)"]
nB<-as.data.frame(table(db$cond)); names(nB)<-c("cond","n")
pA<-ggplot(db,aes(cond,fx,fill=cond))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
  geom_boxplot(width=.55,alpha=.85,outlier.size=.6)+
  geom_text(data=nB,aes(cond,-2.1,label=paste0("n=",n)),inherit.aes=FALSE,size=3.1)+
  scale_fill_manual(values=c(healthy="#7fbf7b",moderate="#f0a202","severe/dead"="#8c510a"),guide="none")+
  annotate("text",-Inf,Inf,hjust=-.06,vjust=1.4,size=3.2,
    label=paste0("moderate vs healthy: p", ifelse(pbark<0.001,"<0.001",sprintf("=%.3f",pbark))))+
  labs(title="A  Bark loss (decay): emission peaks at moderate",
       x="stem condition (bark loss)",y=expression(CH[4]~flux~(arcsinh)))+th

# (B) wounding -> monotonic
dw<-y[!is.na(y$wound)&y$wound<4,]; dw$wf<-factor(dw$wound)   # drop the single 'dead'-coded wound
pw<-summary(lmer(fx~wound+(1|sp),data=y[!is.na(y$wound),]))$coefficients["wound","Pr(>|t|)"]
nW<-as.data.frame(table(dw$wf)); names(nW)<-c("wf","n")
pB<-ggplot(dw,aes(wf,fx,fill=wf))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
  geom_boxplot(width=.55,alpha=.85,outlier.size=.6)+
  geom_text(data=nW,aes(wf,-2.1,label=paste0("n=",n)),inherit.aes=FALSE,size=3.1)+
  scale_fill_brewer(palette="Oranges",guide="none")+
  annotate("text",-Inf,Inf,hjust=-.06,vjust=1.4,size=3.2,label=sprintf("linear trend: p=%.3f",pw))+
  labs(title="B  Wounding (damage): monotonic increase",
       x="wounding / holes (severity)",y=NULL)+th

ggsave("outputs/revision/figS20_stem_deterioration.png", pA+pB, width=10, height=4.6, dpi=300)
cat(sprintf("bark moderate-vs-healthy p=%.3f  wounding linear p=%.3f\n", pbark, pw))
cat("wrote outputs/revision/figS20_stem_deterioration.png\n")
