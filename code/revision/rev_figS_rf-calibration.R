#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS_rf-calibration.R  (SI — RF budget calibration, R3.3/3.5)
# Out-of-bag calibration of the upscaling RF, three views (predicted on x throughout):
#   (a) stem — continuous calibration: OOB points + LOESS smoother vs the 1:1 line
#   (b) soil — continuous calibration
#   (c) stem — aggregate recovery per SPECIES (the unit the tree budget sums over):
#       mean observed vs mean predicted, one point per species
# Point-level OOB R2 is modest, but the budget-relevant quantities (means, the
# magnitude ordering, and the per-species aggregates) are recovered near 1:1.
# Reads outputs/models/TRAINING_DATA.RData (observed + ranger OOB predictions).
# Writes outputs/revision/figS_rf_calibration.png.
# ==============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
has_patch <- requireNamespace("patchwork", quietly=TRUE); if (has_patch) library(patchwork)
options(warn=-1)
asinh10 <- function(x) asinh(x/0.1)/log(10)
ccc <- function(o,p) 2*cov(o,p)/(var(o)+var(p)+(mean(o)-mean(p))^2)
e <- new.env(); load("outputs/models/TRAINING_DATA.RData", envir=e)

of <- function(df, col="y_asinh") sinh(df[[col]])*1000          # nmol m-2 s-1
th <- theme_bw(base_size=10) + theme(panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=8), plot.title=element_text(size=10.5))

# continuous calibration: faint OOB points + LOESS + 1:1 (predicted = x)
# smoother fit only over the central 1-99% of predictions (avoids sparse-tail extrapolation)
cal <- function(df, title, col, brks){
  o <- of(df); p <- sinh(df$pred_asinh)*1000; ok<-is.finite(o)&is.finite(p); o<-o[ok]; p<-p[ok]
  r2 <- 1-mean((o-p)^2)/var(o); cc<-ccc(o,p); rat<-mean(p)/mean(o)
  d <- tibble(px=asinh10(p), oy=asinh10(o))
  qx <- quantile(d$px, c(0.01,0.99)); lim <- range(c(d$px,d$oy))
  ggplot(d, aes(px, oy)) +
    geom_abline(slope=1,intercept=0,linetype="dashed",color="grey55") +
    geom_point(alpha=0.16, size=1.1, color=col) +
    geom_smooth(data=subset(d, px>=qx[1] & px<=qx[2]), method="loess", se=TRUE,
                color=col, fill=col, alpha=0.18, linewidth=0.9, span=0.9) +
    scale_x_continuous(breaks=asinh10(brks), labels=brks, limits=lim) +
    scale_y_continuous(breaks=asinh10(brks), labels=brks, limits=lim) +
    coord_equal() +
    labs(x=expression("Predicted CH"[4]*" (nmol m"^-2*" s"^-1*")"), y=expression("Observed CH"[4]),
         title=title, subtitle=sprintf("OOB R2 = %.2f | CCC = %.2f | mean pred/obs = %.2f", r2,cc,rat)) + th
}

pa <- cal(e$tree_train_complete, "(a) Stem calibration", "#b2182b", c(-0.02,0,0.05,0.1,0.2,0.5))
pb <- cal(e$soil_train_complete, "(b) Soil calibration", "#2166ac", c(-3,-1,-0.3,0,0.3))

# (c) per-species aggregate recovery (stem)
sp <- e$tree_train_complete %>%
  mutate(o=of(.), p=sinh(pred_asinh)*1000) %>% filter(is.finite(o),is.finite(p)) %>%
  group_by(species) %>% summarise(mo=mean(o), mp=mean(p), n=n(), .groups="drop") %>% filter(n>=3)
lim3 <- range(asinh10(c(sp$mo,sp$mp)))
pc <- ggplot(sp, aes(asinh10(mp), asinh10(mo))) +
  geom_abline(slope=1,intercept=0,linetype="dashed",color="grey55") +
  geom_point(aes(size=n), alpha=0.8, color="#b2182b") +
  scale_size_continuous(range=c(1.6,5), name="n meas") +
  scale_x_continuous(breaks=asinh10(c(0,0.05,0.1,0.2)), labels=c(0,0.05,0.1,0.2), limits=lim3) +
  scale_y_continuous(breaks=asinh10(c(0,0.05,0.1,0.2)), labels=c(0,0.05,0.1,0.2), limits=lim3) +
  coord_equal() +
  labs(x=expression("Predicted CH"[4]*" (species mean)"), y=expression("Observed CH"[4]*" (species mean)"),
       title="(c) Per-species means (stem)",
       subtitle=sprintf("%d species | CCC = %.2f", nrow(sp), ccc(asinh10(sp$mo),asinh10(sp$mp)))) +
  th + theme(legend.position=c(0.82,0.22), legend.background=element_rect(fill=NA),
             legend.key.size=unit(0.3,"cm"), legend.title=element_text(size=7), legend.text=element_text(size=6.5))

g <- if (has_patch) pa + pb + pc + plot_layout(nrow=1) else gridExtra::arrangeGrob(pa,pb,pc,ncol=3)
ggsave("outputs/revision/figS_rf_calibration.png", g, width=12, height=4.4, dpi=200, bg="white")
cat("wrote outputs/revision/figS_rf_calibration.png\n")
