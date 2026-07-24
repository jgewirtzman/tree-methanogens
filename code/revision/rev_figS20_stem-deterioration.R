#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS20_stem-deterioration.R  (SI, exploratory)
# CH4 flux vs stem condition (2023 cross-species survey), species-controlled.
# Two ordinal condition metrics, shown as species-controlled estimated marginal
# means +/- 95% CI PER category (mixed model arcsinh(flux) ~ category + (1|sp)),
# in ORIGINAL flux units (back-transformed axis). Both metrics agree that
# emission RISES into moderate deterioration (robust to mean/median/rank). They
# differ in COVERAGE, not necessarily mechanism:
#   (A) Bark loss spans the full gradient (healthy -> dead), so the decline at the
#       degraded end is observable: emission peaks at moderate (hump; quad p<0.001).
#   (B) Wounding is only sampled over its rising range (no well-populated degraded
#       category), so it shows the increase but cannot test whether it too declines.
# Categories use 3 levels (severe+dead merged) to avoid singleton categories.
# Reads data/; writes outputs/revision/figS20_stem_deterioration.png
# ==============================================================================
suppressMessages({library(ggplot2);library(patchwork);library(lme4);library(lmerTest);library(emmeans)})
options(warn=-1,stringsAsFactors=FALSE)
asinh10<-function(x) asinh(x/0.1)/log(10)
ord<-function(v){v<-tolower(trimws(as.character(v)));x<-suppressWarnings(as.numeric(v));x[v=="dead"]<-4;x}
y<-read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv",check.names=FALSE)
names(y)<-make.names(names(y)); y<-y[!is.na(y$CH4_best.flux),]
y$fx<-asinh10(y$CH4_best.flux); y$sp<-as.factor(y$Species)
y$bark<-ord(y$Bark.Missing); y$wound<-ord(y$Wounding.Holes)
obrk<-c(-0.2,0,0.05,0.2,1,5)                       # original-unit tick positions
pfmt<-function(p) ifelse(p<0.001,"p<0.001",sprintf("p=%.3f",p))
th<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.subtitle=element_text(size=9))

panel<-function(metric,labs,title,sub,note,trend){
  d<-y[!is.na(y[[metric]]),]; d$m<-pmin(d[[metric]],3); d$cf<-factor(d$m,labels=labs)
  m<-lmer(fx~cf+(1|sp),data=d); em<-as.data.frame(emmeans(m,~cf))
  d$mc<-d$m-mean(d$m)
  p_lin<-summary(lmer(fx~m+(1|sp),data=d))$coefficients["m","Pr(>|t|)"]
  p_quad<-summary(lmer(fx~mc+I(mc^2)+(1|sp),data=d))$coefficients["I(mc^2)","Pr(>|t|)"]
  lab<-if(trend=="hump") sprintf("hump: quad %s",pfmt(p_quad)) else sprintf("rise: linear %s",pfmt(p_lin))
  nl<-as.data.frame(table(d$cf)); names(nl)<-c("cf","n")
  ggplot(d,aes(cf,fx))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
    geom_jitter(width=.12,height=0,alpha=.16,size=1,colour="#b8732b")+
    geom_line(data=em,aes(as.integer(cf),emmean),colour="#5e3a06",linewidth=.8,inherit.aes=FALSE)+
    geom_errorbar(data=em,aes(cf,ymin=lower.CL,ymax=upper.CL),width=.13,linewidth=.7,inherit.aes=FALSE)+
    geom_point(data=em,aes(cf,emmean),size=2.8,colour="#5e3a06",inherit.aes=FALSE)+
    geom_text(data=nl,aes(cf,asinh10(-0.28),label=paste0("n=",n)),inherit.aes=FALSE,size=2.8,colour="grey35")+
    annotate("text",Inf,asinh10(5),hjust=1.05,vjust=1,size=3,label=lab)+
    annotate("text",-Inf,asinh10(-0.2),hjust=-.05,vjust=1,size=2.7,colour="grey40",label=note)+
    scale_y_continuous(breaks=asinh10(obrk),labels=obrk)+
    coord_cartesian(ylim=asinh10(c(-0.32,6)))+
    labs(title=title,subtitle=sub,x=NULL,y=expression(CH[4]~flux~(nmol~m^-2~s^-1)))+th
}
pA<-panel("bark",c("healthy","moderate","severe/\ndead"),
  "A  Bark loss: full gradient (healthy to dead)",
  "emission peaks at moderate, then declines","spans healthy to dead","hump")
pB<-panel("wound",c("1","2","3+"),
  "B  Wounding: rising range only",
  "emission rises; too few degraded stems to test a decline","no well-sampled degraded end","rise")
ggsave("outputs/revision/figS20_stem_deterioration.png",
  (pA | pB+labs(y=NULL)) + plot_annotation(
    caption="Species-controlled marginal means +/- 95% CI; pattern robust to mean/median/rank. 2023 survey (n=282)."),
  width=10.5,height=5,dpi=300)
cat("wrote outputs/revision/figS20_stem_deterioration.png\n")
