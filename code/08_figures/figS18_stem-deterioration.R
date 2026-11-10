#!/usr/bin/env Rscript
# ==============================================================================
# figS18_stem-deterioration.R -- Figure S19: stem flux against visible deterioration
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root. Exploratory.
# CH4 flux vs stem condition (2023 cross-species survey), species-controlled.
# For each ordinal condition metric (bark loss, wounding; FULL 4 categories, not
# aggregated) we fit continuous LINEAR vs QUADRATIC mixed models
#   arcsinh(flux) ~ level (+ level^2) + (1|species)
# and show the actual fitted curves (+/- 95% CI) over the raw data, with model
# comparison (AIC, marginal R^2, quadratic-term p). Lets the data adjudicate the
# shape rather than assuming it. Axis in original flux units (back-transformed).
# Inputs: data/tables/flux_stem.csv (2023_cross_species: bark_missing, wounding; scores
#   1-3, "dead" = 4) and trees.csv (species), via code/lib/tables_read.R
# Output: outputs/figures/generated/figS20_stem_deterioration.png
# ==============================================================================
suppressMessages({library(ggplot2);library(patchwork);library(lme4);library(lmerTest)})
source("code/lib/flux_transform.R")   # flux_asinh10()
ord<-function(v){v<-tolower(trimws(as.character(v)));x<-suppressWarnings(as.numeric(v));x[v=="dead"]<-4;x}
source("code/lib/fig_style.R")     # SI house style: theme_si(), si_save()
source("code/lib/tables_read.R")   # 2023 deployments (flux_stem.csv); species from trees.csv
y<-flux_stem_campaign("2023_cross_species"); y<-y[!is.na(y$CH4_flux_nmol_m2_s),]
y$fx<-flux_asinh10(y$CH4_flux_nmol_m2_s); y$sp<-as.factor(y$species_code)
y$bark<-ord(y$bark_missing); y$wound<-ord(y$wounding)
obrk<-c(-0.2,0,0.05,0.2,1,5); th<-theme_si()
R2m<-function(m){ vf<-var(predict(m,re.form=NA)); vr<-sum(as.numeric(VarCorr(m))); ve<-sigma(m)^2; vf/(vf+vr+ve) }
predci<-function(m,g,quad){ X<-if(quad) model.matrix(~x+I(x^2),g) else model.matrix(~x,g)
  b<-fixef(m); V<-as.matrix(vcov(m)); g$fit<-as.vector(X%*%b); g$se<-sqrt(diag(X%*%V%*%t(X)))
  g$lo<-g$fit-1.96*g$se; g$hi<-g$fit+1.96*g$se; g }

panel<-function(metric,xlab,catlabs){
  d<-y[!is.na(y[[metric]]),]; d$x<-d[[metric]]
  ml<-lmer(fx~x+(1|sp),data=d); mq<-lmer(fx~x+I(x^2)+(1|sp),data=d)
  pq<-summary(mq)$coefficients["I(x^2)","Pr(>|t|)"]; pl<-summary(ml)$coefficients["x","Pr(>|t|)"]
  aic_l<-AIC(update(ml,REML=FALSE)); aic_q<-AIC(update(mq,REML=FALSE)); r2l<-R2m(ml); r2q<-R2m(mq)
  best<-if(aic_q<aic_l-2) "quadratic preferred" else if(aic_l<aic_q-2) "linear preferred" else "linear and quadratic similar"
  g<-data.frame(x=seq(min(d$x),max(d$x),.05))
  gl<-predci(ml,g,FALSE); gq<-predci(mq,g,TRUE)
  nl<-as.data.frame(table(d$x)); names(nl)<-c("x","n"); nl$x<-as.numeric(as.character(nl$x))
  pf<-function(p) if(p<0.001) "< 0.001" else sprintf("%.3f",p)
  mn<-function(v,d=1) sub("^-","\u2212",sprintf(paste0("%.",d,"f"),v))   # true minus sign
  # model comparison as a small table at the top of the panel (npc positions)
  xl<-0.03; xa<-0.70; xb<-0.95; yr<-c(0.965,0.905,0.845,0.785); ya<-0.715
  tab<-list(
    annotate("text",x=I(c(xa,xb)),y=I(yr[1]),label=c("Linear","Quadratic"),hjust=1,vjust=1,size=2.8,fontface="bold"),
    annotate("text",x=I(xl),y=I(yr[2:4]),label=c("AIC","Marginal~R^2","italic(p)~'(slope; squared term)'"),
             hjust=0,vjust=1,size=2.8,parse=TRUE),
    annotate("text",x=I(rep(xa,3)),y=I(yr[2:4]),label=c(mn(aic_l),sprintf("%.3f",r2l),pf(pl)),hjust=1,vjust=1,size=2.8),
    annotate("text",x=I(rep(xb,3)),y=I(yr[2:4]),label=c(mn(aic_q),sprintf("%.3f",r2q),pf(pq)),hjust=1,vjust=1,size=2.8),
    annotate("text",x=I(xl),y=I(ya),label=sprintf("\u0394AIC = %s: %s",mn(aic_l-aic_q),best),hjust=0,vjust=1,size=2.8))
  ggplot(d,aes(x,fx))+geom_hline(yintercept=0,linetype=3,colour="grey60")+
    geom_jitter(width=.12,height=0,alpha=.16,size=1,colour="#b8732b")+
    geom_ribbon(data=gq,aes(x,ymin=lo,ymax=hi),alpha=.15,fill="#5e3a06",inherit.aes=FALSE)+
    geom_line(data=gl,aes(x,fit),linetype=2,linewidth=.8,colour="grey25",inherit.aes=FALSE)+
    geom_line(data=gq,aes(x,fit),linewidth=1,colour="#5e3a06",inherit.aes=FALSE)+
    geom_text(data=nl,aes(x,flux_asinh10(-0.28),label=paste0("n = ",n)),inherit.aes=FALSE,size=2.8,colour="grey35")+
    tab+
    scale_x_continuous(breaks=1:4,labels=catlabs)+
    scale_y_continuous(breaks=flux_asinh10(obrk),labels=obrk)+coord_cartesian(ylim=flux_asinh10(c(-0.32,30)))+   # headroom for the model table
    labs(x=xlab,y=expression(CH[4]~flux~(nmol~m^-2~s^-1)))+th
}
# both scores: healthy, moderate, severe (1-3) or dead stem (4) (SI Methods, Bark loss and wounding)
cl<-c("Healthy","Moderate","Severe","Dead")
pA<-panel("bark","Bark loss (decay)",cl)
pB<-panel("wound","Wounding (damage)",cl)
# model note (species random intercept, n, robustness) is in the SI caption
si_save("outputs/figures/generated/figS20_stem_deterioration.png", (pA | pB+labs(y=NULL)) + si_tags(), height=4.4)
cat("wrote outputs/figures/generated/figS20_stem_deterioration.png\n")
