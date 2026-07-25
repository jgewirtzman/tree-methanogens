#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS_rf-calibration.R  (SI — RF budget calibration, R3.3/3.5)
# Out-of-bag calibration of the upscaling RF: bin trees/soil plots by PREDICTED
# flux and plot mean OBSERVED vs mean PREDICTED per bin against the 1:1 line.
# Point-level OOB R2 is modest, but the budget-relevant quantities (means,
# magnitude-by-magnitude ordering) are recovered — this figure shows that.
# Reads outputs/models/TRAINING_DATA.RData (observed + ranger OOB predictions).
# Writes outputs/revision/figS_rf_calibration.png.
# ==============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
has_patch <- requireNamespace("patchwork", quietly=TRUE); if (has_patch) library(patchwork)
options(warn=-1)
asinh10 <- function(x) asinh(x/0.1)/log(10); inv10 <- function(y) 0.1*sinh(y*log(10))
ccc <- function(o,p) 2*cov(o,p)/(var(o)+var(p)+(mean(o)-mean(p))^2)

e <- new.env(); load("outputs/models/TRAINING_DATA.RData", envir=e)

panel <- function(df, ycol, title, col, brks){
  o <- sinh(df[[ycol]])*1000; p <- sinh(df$pred_asinh)*1000     # nmol m-2 s-1
  ok <- is.finite(o)&is.finite(p); o<-o[ok]; p<-p[ok]
  r2  <- 1 - mean((o-p)^2)/var(o); cc <- ccc(o,p); rat <- mean(p)/mean(o)
  d <- tibble(o,p) %>% mutate(bin = ntile(p, 10)) %>%
    group_by(bin) %>% summarise(mp=mean(p), mo=mean(o), se=sd(o)/sqrt(n()), .groups="drop")
  lim <- range(asinh10(c(d$mp, d$mo, d$mo+d$se, d$mo-d$se)), na.rm=TRUE)
  ggplot(d, aes(asinh10(mp), asinh10(mo))) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey55") +
    geom_errorbar(aes(ymin=asinh10(mo-se), ymax=asinh10(mo+se)), width=0, color="grey60", linewidth=0.5) +
    geom_point(size=2.6, color=col) +
    scale_x_continuous(breaks=asinh10(brks), labels=brks, limits=lim) +
    scale_y_continuous(breaks=asinh10(brks), labels=brks, limits=lim) +
    coord_equal() +
    labs(x=expression("Predicted CH"[4]*" (bin mean, nmol m"^-2*" s"^-1*")"),
         y=expression("Observed CH"[4]*" (bin mean)"),
         title=title,
         subtitle=sprintf("OOB R2 = %.2f  |  CCC = %.2f  |  mean pred/obs = %.2f", r2, cc, rat)) +
    theme_bw(base_size=10) +
    theme(panel.grid.minor=element_blank(), plot.subtitle=element_text(size=8.5))
}

pt <- panel(e$tree_train_complete, "y_asinh", "(a) Stem flux", "#b2182b", c(0,0.05,0.1,0.2,0.4))
ps <- panel(e$soil_train_complete, "y_asinh", "(b) Soil flux", "#2166ac", c(-3,-1,-0.3,0))

g <- if (has_patch) pt + ps else gridExtra::arrangeGrob(pt, ps, ncol=2)
ggsave("outputs/revision/figS_rf_calibration.png", g, width=9, height=4.7, dpi=200, bg="white")
cat("wrote outputs/revision/figS_rf_calibration.png\n")
