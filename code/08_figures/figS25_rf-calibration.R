#!/usr/bin/env Rscript
# ==============================================================================
# figS25_rf-calibration.R -- Figure S27: calibration of the random-forest predictions
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
# Out-of-bag calibration of the upscaling forests: bins by PREDICTED flux (deciles),
# mean observed (+/- SE) against mean predicted, with the 1:1 line; (a) stem (b) soil.
# Binned means only (the point scatter is Figure S26). Per-species and per-moisture
# aggregations are not shown: species and moisture are model predictors, so they would
# be near-circular.
# Input: outputs/models/TRAINING_DATA.RData (y_flux, pred_flux; 05_model/02_rf_models.R)
# Output: outputs/figures/generated/figS_rf_calibration.png
# ==============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
source("code/lib/flux_transform.R")   # flux_asinh10()
ccc <- function(o,p) 2*cov(o,p)/(var(o)+var(p)+(mean(o)-mean(p))^2)
e <- new.env(); load("outputs/models/TRAINING_DATA.RData", envir=e)
th <- theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=8), plot.title=element_text(size=10.5))

# observed (y_flux) and out-of-bag predicted (pred_flux) flux, nmol m-2 s-1 (02_rf_models.R)
flux <- function(df) list(o=df$y_flux, p=df$pred_flux)

# obs-vs-pred panel from a table of bin means (mp, mo, se)
obsvpred <- function(d, title, sub, col, brks){
  lim <- range(flux_asinh10(c(d$mp,d$mo,d$mo+d$se,d$mo-d$se)),na.rm=TRUE)
  g <- ggplot(d, aes(flux_asinh10(mp), flux_asinh10(mo))) +
    geom_abline(slope=1,intercept=0,linetype="dashed",color="grey55") +
    geom_errorbar(aes(ymin=flux_asinh10(mo-se),ymax=flux_asinh10(mo+se)),width=0,color="grey65",linewidth=0.45)
  g + geom_point(size=2.6,color=col) +
      scale_x_continuous(breaks=flux_asinh10(brks),labels=brks,limits=lim) +
      scale_y_continuous(breaks=flux_asinh10(brks),labels=brks,limits=lim) + coord_equal() +
      labs(x=expression("Predicted CH"[4]*" (nmol m"^-2*" s"^-1*")"),
           y=expression("Observed CH"[4]), title=title, subtitle=sub) + th
}

## calibration by predicted decile
cal <- function(df, title, col, brks){
  f<-flux(df); ok<-is.finite(f$o)&is.finite(f$p); o<-f$o[ok]; p<-f$p[ok]
  d <- tibble(o,p) %>% mutate(b=ntile(p,10)) %>% group_by(b) %>%
       summarise(mp=mean(p),mo=mean(o),se=sd(o)/sqrt(n()),.groups="drop")
  obsvpred(d,title,sprintf("OOB R2 = %.2f | CCC = %.2f | mean pred/obs = %.2f",
           1-mean((o-p)^2)/var(o),ccc(o,p),mean(p)/mean(o)),col,brks)
}
pa <- cal(e$tree_train_complete,"(a) Stem flux","#b2182b",c(0,0.05,0.1,0.2,0.4))
pb <- cal(e$soil_train_complete,"(b) Soil flux","#2166ac",c(-3,-1,-0.3,0))

ggsave("outputs/figures/generated/figS_rf_calibration.png", pa|pb, width=8.8, height=4.6, dpi=200, bg="white")
cat("wrote outputs/figures/generated/figS_rf_calibration.png\n")
