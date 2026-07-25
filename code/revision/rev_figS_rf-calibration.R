#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS_rf-calibration.R  (SI — RF budget calibration, R3.3/3.5)
# Out-of-bag aggregate accuracy of the upscaling RF (predicted on x throughout;
# binned means only — the raw point scatter is already in the RF-predictions SI fig):
#   TOP  — calibration: bin by PREDICTED flux (deciles); mean observed vs mean predicted
#          (a) stem  (b) soil
#   BOT  — budget-unit recovery: mean observed vs mean predicted within the units the
#          budget sums over — (c) per SPECIES (stem)  (d) per SOIL-MOISTURE bin
# Modest point-level OOB R2, but means/magnitudes recovered near 1:1 in every view.
# Reads outputs/models/TRAINING_DATA.RData. Writes outputs/revision/figS_rf_calibration.png.
# ==============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
has_patch <- requireNamespace("patchwork", quietly=TRUE); if (has_patch) library(patchwork)
options(warn=-1)
asinh10 <- function(x) asinh(x/0.1)/log(10)
ccc <- function(o,p) 2*cov(o,p)/(var(o)+var(p)+(mean(o)-mean(p))^2)
e <- new.env(); load("outputs/models/TRAINING_DATA.RData", envir=e)
th <- theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=8), plot.title=element_text(size=10.5))

flux <- function(df,col="y_asinh") list(o=sinh(df[[col]])*1000, p=sinh(df$pred_asinh)*1000)

# obs-vs-pred panel from a table of bin means (mp, mo, se); optional color aesthetic
obsvpred <- function(d, title, sub, col, brks, colour_by=NULL, col_name=NULL, col_grad=NULL){
  lim <- range(asinh10(c(d$mp,d$mo,d$mo+d$se,d$mo-d$se)),na.rm=TRUE)
  g <- ggplot(d, aes(asinh10(mp), asinh10(mo))) +
    geom_abline(slope=1,intercept=0,linetype="dashed",color="grey55") +
    geom_errorbar(aes(ymin=asinh10(mo-se),ymax=asinh10(mo+se)),width=0,color="grey65",linewidth=0.45)
  if(is.null(colour_by)) g <- g + geom_point(size=2.6,color=col) else
    g <- g + geom_point(aes(colour=.data[[colour_by]]),size=2.8) +
             scale_colour_gradientn(colours=col_grad,name=col_name)
  g + scale_x_continuous(breaks=asinh10(brks),labels=brks,limits=lim) +
      scale_y_continuous(breaks=asinh10(brks),labels=brks,limits=lim) + coord_equal() +
      labs(x=expression("Predicted CH"[4]*" (nmol m"^-2*" s"^-1*")"),
           y=expression("Observed CH"[4]), title=title, subtitle=sub) + th
}

## TOP — calibration by predicted decile
cal <- function(df, title, col, brks){
  f<-flux(df); ok<-is.finite(f$o)&is.finite(f$p); o<-f$o[ok]; p<-f$p[ok]
  d <- tibble(o,p) %>% mutate(b=ntile(p,10)) %>% group_by(b) %>%
       summarise(mp=mean(p),mo=mean(o),se=sd(o)/sqrt(n()),.groups="drop")
  obsvpred(d,title,sprintf("OOB R2 = %.2f | CCC = %.2f | mean pred/obs = %.2f",
           1-mean((o-p)^2)/var(o),ccc(o,p),mean(p)/mean(o)),col,brks)
}
pa <- cal(e$tree_train_complete,"(a) Stem: calibration (predicted deciles)","#b2182b",c(0,0.05,0.1,0.2,0.4))
pb <- cal(e$soil_train_complete,"(b) Soil: calibration (predicted deciles)","#2166ac",c(-3,-1,-0.3,0))

## BOT-c — per species (stem)
f<-flux(e$tree_train_complete)
spc <- tibble(o=f$o,p=f$p,sp=e$tree_train_complete$species) %>% filter(is.finite(o),is.finite(p)) %>%
  group_by(sp) %>% summarise(mp=mean(p),mo=mean(o),se=sd(o)/sqrt(n()),n=n(),.groups="drop") %>% filter(n>=3)
pc <- obsvpred(spc,"(c) Stem: per-species means",sprintf("%d species | CCC = %.2f",nrow(spc),ccc(asinh10(spc$mo),asinh10(spc$mp))),
               "#b2182b",c(0,0.05,0.1,0.2))

## BOT-d — per soil-moisture bin
f<-flux(e$soil_train_complete)
mob <- tibble(o=f$o,p=f$p,m=e$soil_train_complete$soil_moisture_at_site) %>% filter(is.finite(o),is.finite(p),is.finite(m)) %>%
  mutate(b=ntile(m,10)) %>% group_by(b) %>% summarise(mp=mean(p),mo=mean(o),se=sd(o)/sqrt(n()),moist=mean(m),.groups="drop")
pd <- obsvpred(mob,"(d) Soil: per-moisture-bin means",sprintf("10 moisture bins | CCC = %.2f",ccc(asinh10(mob$mo),asinh10(mob$mp))),
               NULL,c(-3,-1,-0.3,0),colour_by="moist",col_name="VWC",col_grad=c("#d1e5f0","#4393c3","#2166ac"))

g <- if (has_patch) (pa|pb)/(pc|pd) else gridExtra::arrangeGrob(pa,pb,pc,pd,ncol=2)
ggsave("outputs/revision/figS_rf_calibration.png", g, width=8.6, height=8.4, dpi=200, bg="white")
cat("wrote outputs/revision/figS_rf_calibration.png\n")
