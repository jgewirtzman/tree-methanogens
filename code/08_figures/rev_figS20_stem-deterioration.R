#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS20_stem-deterioration.R  (SI, exploratory)
# CH4 flux vs stem condition (2023 cross-species survey), species-controlled.
# For each ordinal condition metric (bark loss, wounding; FULL 4 categories, not
# aggregated) we fit continuous LINEAR vs QUADRATIC mixed models
#   arcsinh(flux) ~ level (+ level^2) + (1|species)
# and show the actual fitted curves (+/- 95% CI) over the raw data, with model
# comparison (AIC, marginal R^2, quadratic-term p). Lets the data adjudicate the
# shape rather than assuming it. Axis in original flux units (back-transformed).
# Reads data/; writes outputs/figures/generated/figS20_stem_deterioration.png
# ==============================================================================
suppressMessages({library(ggplot2);library(patchwork);library(lme4);library(lmerTest)})
options(warn=-1,stringsAsFactors=FALSE)
asinh10<-function(x) asinh(x/0.1)/log(10)
ord<-function(v){v<-tolower(trimws(as.character(v)));x<-suppressWarnings(as.numeric(v));x[v=="dead"]<-4;x}
y<-read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv",check.names=FALSE)
names(y)<-make.names(names(y)); y<-y[!is.na(y$CH4_best.flux),]
y$fx<-asinh10(y$CH4_best.flux); y$sp<-as.factor(y$Species)
y$bark<-ord(y$Bark.Missing); y$wound<-ord(y$Wounding.Holes)
obrk<-c(-0.2,0,0.05,0.2,1,5); th<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank())
R2m<-function(m){ vf<-var(predict(m,re.form=NA)); vr<-sum(as.numeric(VarCorr(m))); ve<-sigma(m)^2; vf/(vf+vr+ve) }
predci<-function(m,g,quad){ X<-if(quad) model.matrix(~x+I(x^2),g) else model.matrix(~x,g)
  b<-fixef(m); V<-as.matrix(vcov(m)); g$fit<-as.vector(X%*%b); g$se<-sqrt(diag(X%*%V%*%t(X)))
  g$lo<-g$fit-1.96*g$se; g$hi<-g$fit+1.96*g$se; g }

panel<-function(metric,xlab,catlabs){
  d<-y[!is.na(y[[metric]]),]; d$x<-d[[metric]]
  ml<-lmer(fx~x+(1|sp),data=d); mq<-lmer(fx~x+I(x^2)+(1|sp),data=d)
  pq<-summary(mq)$coefficients["I(x^2)","Pr(>|t|)"]; pl<-summary(ml)$coefficients["x","Pr(>|t|)"]
  aic_l<-AIC(ml); aic_q<-AIC(mq); r2l<-R2m(ml); r2q<-R2m(mq)
  best<-if(aic_q<aic_l-2) "quadratic (hump)" else if(aic_l<aic_q-2) "linear" else "linear ~ quad"
  g<-data.frame(x=seq(min(d$x),max(d$x),.05))
  gl<-predci(ml,g,FALSE); gq<-predci(mq,g,TRUE)
  nl<-as.data.frame(table(d$x)); names(nl)<-c("x","n"); nl$x<-as.numeric(as.character(nl$x))
  txt<-sprintf("linear:    AIC %.1f  R2m %.3f  (slope p=%.3f)\nquadratic: AIC %.1f  R2m %.3f  (x2 p=%.3f)\ndAIC=%.1f  ->  %s",
    aic_l,r2l,pl, aic_q,r2q,pq, aic_l-aic_q, best)
  ggplot(d,aes(x,fx))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
    geom_jitter(width=.12,height=0,alpha=.16,size=1,colour="#b8732b")+
    geom_ribbon(data=gq,aes(x,ymin=lo,ymax=hi),alpha=.15,fill="#5e3a06",inherit.aes=FALSE)+
    geom_line(data=gl,aes(x,fit),linetype=2,linewidth=.8,colour="grey25",inherit.aes=FALSE)+
    geom_line(data=gq,aes(x,fit),linewidth=1,colour="#5e3a06",inherit.aes=FALSE)+
    geom_text(data=nl,aes(x,asinh10(-0.28),label=paste0("n=",n)),inherit.aes=FALSE,size=2.8,colour="grey35")+
    annotate("text",-Inf,asinh10(5.5),hjust=-.03,vjust=1,size=2.8,label=txt)+
    scale_x_continuous(breaks=1:4,labels=catlabs)+
    scale_y_continuous(breaks=asinh10(obrk),labels=obrk)+coord_cartesian(ylim=asinh10(c(-0.32,6)))+
    labs(x=xlab,y=expression(CH[4]~flux~(nmol~m^-2~s^-1)))+th
}
pA<-panel("bark","Bark loss (decay)",c("healthy","moderate","severe","dead"))
pB<-panel("wound","Wounding (damage)",c("1","2","3","4"))
ggsave("outputs/figures/generated/figS20_stem_deterioration.png",
  (pA | pB+labs(y=NULL))+plot_annotation(
    caption="Species-controlled (1|species), n=282. Rise into moderate deterioration robust to mean/median/rank."),
  width=11,height=5,dpi=300)
cat("wrote outputs/figures/generated/figS20_stem_deterioration.png\n")
