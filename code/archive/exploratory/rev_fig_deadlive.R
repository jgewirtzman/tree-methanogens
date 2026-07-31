#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_deadlive.R   (REVISION — live vs dead stem CH4, done properly)
# Two complementary views of the dead-vs-live stem-flux contrast:
#   (A) emmeans-adjusted means (+/-95% CI) from a mixed model controlling for
#       species + plot + season, with a tree random effect (repeated measures).
#       Honest but underpowered: 8 dead trees across 6 species / 3 plots.
#   (B) matched within-species x plot dumbbells: each dead tree vs live trees of
#       the SAME species in the SAME plot (controls both) -> dead > matched live 8/8.
# Unit is the TREE (medians across each tree's QC'd measurements) for display;
# the model uses measurement level so season is a within-tree control.
# PARKED / EXPLORATORY: dead-vs-live is underpowered (n=8 dead; effect is a median/
# rank signal that vanishes on the budget-relevant mean, p=0.97). Not a manuscript
# result -> one descriptive sentence in the decay discussion. Kept as the record.
# Writes outputs/revision/exploratory/fig_deadlive.png (+ a stats sidecar .txt).
# ==============================================================================
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(lme4);library(lmerTest);library(emmeans)})
has_patch<-requireNamespace("patchwork",quietly=TRUE); if(has_patch) library(patchwork)
options(warn=-1); num<-function(x)suppressWarnings(as.numeric(x)); F<-function(v){v[is.na(v)]<-FALSE;v}
asinh10<-function(x)asinh(x/0.1)/log(10); inv10<-function(y)0.1*sinh(y*log(10))
sp_map<-c(ACRU="A. rubrum",BEAL="B. alleghaniensis",BELE="B. lenta",FRAM="F. americana",
          QURU="Q. rubra",TSCA="T. canadensis",PIST="P. strobus",QUAL="Q. alba")

## ---- data: QC'd monthly (merged), measurement level ----
ym<-read.csv("data/raw/field_data/static_chamber_field/YM_trees_measured.csv.csv",check.names=FALSE)
sp_lut<-setNames(toupper(as.character(ym$Species)),as.character(ym$Label))
y23<-read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv",check.names=FALSE)
sp23<-setNames(toupper(as.character(y23[["Species Code"]])),as.character(y23[["Tree Tag"]]))
m<-read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv")
m<-m[!is.na(num(m$CH4_best.flux.x)),]
m<-m[F(num(m$CO2_LM.r2)>=.7|num(m$CO2_HM.r2)>=.7|num(m$CH4_LM.r2.x)>=.7|num(m$CH4_HM.r2.x)>=.7)&
     F(num(m$CH4_best.flux.x)>=-100&num(m$CH4_best.flux.x)<=200),]
m$tag<-trimws(as.character(m$Plot.Tag)); m$plot<-recode(m$Plot.Letter,WD="W",WS="W",.default=m$Plot.Letter)
gsp<-function(t) ifelse(grepl("^UNTAG_",t), toupper(sub("^UNTAG_[A-Za-z]+_([A-Za-z]+).*$","\\1",t)),
                        ifelse(t%in%names(sp_lut),sp_lut[t],sp23[t]))
m$sp<-gsp(m$tag); m$dead<-F(m$dead); m$flux<-num(m$CH4_best.flux.x); m$fx<-asinh10(m$flux)
mo<-as.integer(format(as.Date(m$Date),"%m"))
m$season<-factor(c("Winter","Winter","Spring","Spring","Spring","Summer","Summer","Summer","Fall","Fall","Fall","Winter")[mo],
                 levels=c("Winter","Spring","Summer","Fall"))
m<-m[!is.na(m$sp)&m$sp!="",]
m$Status<-factor(ifelse(m$dead,"Dead","Live"),levels=c("Live","Dead"))
m$plotf<-recode(m$plot,U="Upland",I="Intermediate",W="Wetland")
tr<-m%>%group_by(tag,plotf,sp,Status)%>%summarise(med=median(flux),y=asinh10(median(flux)),.groups="drop")

## ---- model + emmeans (Panel A) ----
fit<-lmer(fx~dead+sp+plot+season+(1|tag),data=m)
em<-as.data.frame(summary(emmeans(fit,~dead)))
em$Status<-factor(ifelse(em$dead,"Dead","Live"),levels=c("Live","Dead"))
ctr<-summary(contrast(emmeans(fit,~dead),"revpairwise"))
pA<-ctr$p.value[1]

## ---- matched within-species x plot contrast (Panel B) ----
## each dead tree vs LIVE trees of the SAME species in the SAME plot (controls both).
library(tidyr)
cellL<-tr%>%filter(Status=="Live")%>%group_by(sp,plotf)%>%summarise(live=median(med),n_live=n(),.groups="drop")
dd<-tr%>%filter(Status=="Dead")%>%left_join(cellL,by=c("sp","plotf"))%>%filter(!is.na(live))
dd$Species<-ifelse(dd$sp%in%names(sp_map),sp_map[dd$sp],dd$sp)
dd<-dd%>%arrange(med-live)%>%mutate(rk=factor(row_number()))
rlab<-setNames(sprintf("%s · %s  (n=%d)",dd$Species,dd$plotf,dd$n_live),dd$rk)
signp<-binom.test(sum(dd$med>dd$live),nrow(dd))$p.value
wilp <-wilcox.test(dd$med-dd$live)$p.value
swl<-bind_rows(transmute(dd,rk,Status="Live",med=live),transmute(dd,rk,Status="Dead",med=med))
swl$Status<-factor(swl$Status,levels=c("Live","Dead"))

## ---- shared aesthetics ----
cols<-c(Live="#4575b4",Dead="#d73027")
brk<-c(0,0.01,0.05,0.1,0.2); ax<-scale_y_continuous(breaks=asinh10(brk),labels=brk)
th<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="top")

## Panel A: adjusted emmeans + raw tree points
pa<-ggplot(tr,aes(Status,y,color=Status))+
  geom_jitter(width=.14,height=0,alpha=.5,size=1.8,aes(shape=plotf))+
  geom_errorbar(data=em,aes(x=Status,ymin=asinh10(inv10(lower.CL)),ymax=asinh10(inv10(upper.CL))),
                width=.12,linewidth=.9,inherit.aes=FALSE,color="grey20")+
  geom_point(data=em,aes(Status,emmean),size=4,shape=18,inherit.aes=FALSE,color="grey15")+
  scale_color_manual(values=cols,guide="none")+scale_shape_manual(name="Plot",values=c(Upland=16,Intermediate=17,Wetland=15))+
  ax+labs(x=NULL,y=expression(CH[4]~flux~(nmol~m^-2~s^-1)),
    title="(A) Adjusted for species + plot + season",
    subtitle=sprintf("emmeans +/- 95%% CI (diamond); points = tree medians;  contrast p = %.2f (n.s.)",pA))+th

## Panel B: matched within-species x plot dumbbells (one row per dead tree)
pb<-ggplot(swl,aes(asinh10(med),rk))+
  geom_line(aes(group=rk),color="grey60",linewidth=.8)+
  geom_point(aes(color=Status),size=3.3)+
  scale_color_manual(name=NULL,values=cols)+
  scale_x_continuous(breaks=asinh10(brk),labels=brk)+
  scale_y_discrete(labels=rlab)+
  labs(x=expression(CH[4]~flux~(nmol~m^-2~s^-1)),y=NULL,
    title="(B) Matched within species × plot",
    subtitle=sprintf("each dead tree vs same-species same-plot live cohort\ndead > matched live in %d/%d;  Wilcoxon p = %.3f",sum(dd$med>dd$live),nrow(dd),wilp))+
  th+theme(axis.text.y=element_text(size=8))

if(has_patch){ g<-pa+pb+plot_layout(widths=c(1,1.15)) } else { g<-gridExtra::arrangeGrob(pa,pb,ncol=2) }
ggsave("outputs/revision/exploratory/fig_deadlive.png",g,width=11,height=5,dpi=200,bg="white")

## stats sidecar
sink("outputs/revision/exploratory/fig_deadlive_stats.txt")
cat("DEAD vs LIVE stem CH4 -- two views\n\n(A) lmer(asinh(flux) ~ dead + species + plot + season + (1|tree)); measurement level\n")
print(round(summary(fit)$coefficients[grep("dead",rownames(summary(fit)$coefficients)),,drop=FALSE],4))
cat("\nemmeans (back-transformed nmol m-2 s-1):  Live",round(inv10(em$emmean[em$Status=="Live"]),4),
    "| Dead",round(inv10(em$emmean[em$Status=="Dead"]),4)," contrast p =",round(pA,4),"\n\n")
cat("(B) matched within-species x plot (each dead tree vs same-sp same-plot live median), nmol m-2 s-1:\n")
print(as.data.frame(dd%>%transmute(Species,plot=plotf,n_live,dead=round(med,4),live=round(live,4),diff=round(med-live,4))),row.names=FALSE)
cat(sprintf("\ndead > matched live in %d/%d dead trees; sign test p = %.4f; Wilcoxon signed-rank p = %.4f\n",
    sum(dd$med>dd$live),nrow(dd),signp,wilp))
sink()
cat("wrote outputs/revision/exploratory/fig_deadlive.png + fig_deadlive_stats.txt\n")
cat(sprintf("A: adjusted contrast p=%.3f | B: matched %d/%d dead>live, Wilcoxon p=%.3f\n",pA,sum(dd$med>dd$live),nrow(dd),wilp))
