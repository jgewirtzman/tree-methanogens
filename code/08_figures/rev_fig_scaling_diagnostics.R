#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_scaling_diagnostics.R
# ------------------------------------------------------------------------------
# Diagnostic figure set for the height/area scaling. Three panels groups:
#
#   FIG 1  HEIGHT-FLUX PROFILES
#          a  observed, every three-height tree (paired within-tree)
#          b  observed by species, mean profiles
#          c  RF-learned profiles by species
#          d  functional forms, fitted to the band and projected to the canopy
#
#   FIG 2  SURFACE AREA WITH HEIGHT
#          a  stem height distribution from the allometry
#          b  per-tree area density scenarios (uniform / conic / bole+crown)
#          c  STAND-LEVEL woody area vs height -- the emergent distribution, which is
#             bottom-weighted even when each tree is uniform, because stems differ in
#             height
#          d  cumulative share of woody area below a given height
#
#   FIG 3  WHERE THE BUDGET COMES FROM
#          a  flux x area density vs height
#          b  cumulative share of the whole-surface budget below a given height,
#             under each flux form -- shows how much of the answer is extrapolated
#
# Also reports Kalmia latifolia heights under the unmodified allometry, for the
# decision on whether to apply a shrub cap at all.
# ==============================================================================

suppressMessages({library(ranger);library(dplyr);library(ggplot2);library(tidyr);library(patchwork)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
source("code/lib/rev_geometry.R")
source("code/lib/rev_species_levels.R")
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
d <- tree_train_complete
trained <- sort(unique(as.character(d$species_clean)))
CONV <- 86400*365.25*16e-6
A_PLOT <- STAND_AREA_M2   # censused ground; was hardcoded 200*200, the nominal square

# CANONICAL inventory, not rf_workflow_data$PLACEHOLDER_INVENTORY. That object is the
# pre-rebuild stem list (before the bytag coordinate rescue and the DBH unit fix), and
# the fx() stack that stood here was the over-correcting species-specific DBH repair
# that rev_inventory_build.R replaced with one documented rule. The species line also
# carried the ladder bug -- see rev_species_levels.R.
INVFILE <- "outputs/tables/inventory_stems.csv"
if (!file.exists(INVFILE))
  stop("missing ", INVFILE, " -- run rev_inventory_build.R first")
INV <- read.csv(INVFILE, stringsAsFactors=FALSE)
INV <- INV[INV$in_stand & is.finite(INV$dbh_m) & INV$dbh_m > 0, ]
INV$dbh <- INV$dbh_m          # already unit-checked and typo-repaired upstream
INV$sp  <- species_to_model_level(INV$species, trained)
GY <- c("Pinus strobus","Tsuga canadensis")
# One allometry, from rev_geometry.R. This was a copy-pasted duplicate that retyped
# the canopy anchor and both exponents and recomputed the 95th-percentile DBH
# locally, so changing CANOPY_H_M or the exponents would have moved none of them.
INV$H <- stem_height_m(INV$dbh, INV$species)      # NO shrub cap

# ---- Kalmia check -----------------------------------------------------------
# %in%, not ==. 41 inventory rows have species = NA (unidentified stems, retained by
# design), and == returns NA for those, so INV[NA, ] injects NA rows and every
# statistic below came out NA. The old PLACEHOLDER_INVENTORY had no NA species, so
# this was latent until the canonical inventory was wired in.
K <- INV[INV$species %in% "Kalmia latifolia", ]
cat(sprintf("
KALMIA LATIFOLIA under the unmodified allometry (n = %d)
  DBH  : median %.1f cm, range %.1f - %.1f cm
  H    : median %.2f m, mean %.2f m, range %.2f - %.2f m
  Field expectation for mountain laurel: 1-3 m, occasionally to 6 m.
  Share of stand: %.2f%% of basal area, %.2f%% of conic woody area.
  Effect of a 3 m cap on total woody area: %+.2f%%
",
 nrow(K), 100*median(K$dbh), 100*min(K$dbh), 100*max(K$dbh),
 median(K$H), mean(K$H), min(K$H), max(K$H),
 100*sum(pi*(K$dbh/2)^2)/sum(pi*(INV$dbh/2)^2),
 100*sum(pi*K$dbh*K$H/2)/sum(pi*INV$dbh*INV$H/2),
 100*(sum(pi*INV$dbh*ifelse(INV$species %in% "Kalmia latifolia",pmin(INV$H,3),INV$H)/2) /
      sum(pi*INV$dbh*INV$H/2) - 1)))

# ==============================================================================
# FIG 1  height-flux profiles
# ==============================================================================
pair <- d %>% filter(!is.na(measurement_height_cm), is.finite(stem_flux_corrected)) %>%
  group_by(tree_id, species_clean, h=measurement_height_cm/100) %>%
  summarise(f=mean(stem_flux_corrected), .groups="drop")
W <- pair %>% pivot_wider(names_from=h, values_from=f, names_prefix="h") %>%
  filter(is.finite(h0.5),is.finite(h1.25),is.finite(h2))
sig <- function(x) sign(x)*log10(1+abs(x)/0.01)

p1a <- ggplot(pair %>% semi_join(W,"tree_id"), aes(h, f, group=tree_id)) +
  geom_line(alpha=.16, colour="grey30") +
  stat_summary(aes(group=1), fun=mean, geom="line", colour="#B2182B", linewidth=1.2) +
  stat_summary(aes(group=1), fun=median, geom="line", colour="#2166AC", linewidth=1.2, linetype="22") +
  scale_y_continuous(trans=scales::pseudo_log_trans(sigma=.01, base=10)) +
  annotate("text",x=2,y=8,label="mean",colour="#B2182B",hjust=1,size=2.8) +
  annotate("text",x=2,y=3,label="median",colour="#2166AC",hjust=1,size=2.8) +
  labs(title="a  observed, all three-height trees", x="height (m)",
       y=expression(CH[4]~(nmol~m^-2~s^-1))) + theme_bw(base_size=8)

spmean <- pair %>% semi_join(W,"tree_id") %>% group_by(species_clean,h) %>%
  summarise(f=mean(f), n=n(), .groups="drop") %>% group_by(species_clean) %>%
  filter(min(n)>=3) %>% ungroup()
p1b <- ggplot(spmean, aes(h,f,colour=species_clean)) + geom_line(linewidth=.7) + geom_point(size=1.3) +
  scale_y_continuous(trans=scales::pseudo_log_trans(sigma=.01, base=10)) +
  labs(title="b  observed, species means (n>=3)", x="height (m)", y=NULL, colour=NULL) +
  theme_bw(base_size=8) + theme(legend.text=element_text(size=5.5), legend.key.height=unit(7,"pt"))

gr <- expand.grid(species=trained, height_cm=seq(50,200,12.5), stringsAsFactors=FALSE)
gr$dbh_m <- median(d$dbh_m,na.rm=TRUE)
gr$soil_moisture_at_tree <- 0.22; gr$soil_temp_C_mean <- 15; gr$air_temp_C_mean <- 17
gr$species <- factor(gr$species, levels=trained)
gr$pred <- predict(TreeRF, gr, num.threads=1)$predictions
p1c <- ggplot(gr, aes(height_cm/100, pred, colour=species)) + geom_line(linewidth=.7) +
  labs(title="c  RF-learned profiles by species", x="height (m)",
       y=expression(CH[4]~(nmol~m^-2~s^-1)), colour=NULL) +
  theme_bw(base_size=8) + theme(legend.text=element_text(size=5.5), legend.key.height=unit(7,"pt"))

hq <- c(.5,1.25,2); fb <- sapply(hq,function(z) mean(pair$f[abs(pair$h-z)<1e-6]))
bl <- coef(lm(log(pmax(fb,1e-9))~hq))[2]; bp <- coef(lm(log(pmax(fb,1e-9))~log(hq)))[2]
blin <- coef(lm(fb~hq)); f2b <- fb[3]
zz <- seq(.05,25,length.out=400)
FORMS <- list(constant=rep(1,length(zz)), `power`=exp(bp*(log(zz)-log(2))),
  `log-linear`=exp(bl*(zz-2)), `linear (floored)`=pmax(blin[1]+blin[2]*zz,0)/max(blin[1]+blin[2]*2,1e-9),
  `asymptote (50%)`=0.5+0.5*exp(bl*(zz-2)))
fdf <- bind_rows(lapply(names(FORMS),function(n) data.frame(z=zz,f=f2b*FORMS[[n]],form=n)))
p1d <- ggplot(fdf,aes(z,f,colour=form)) +
  annotate("rect",xmin=0,xmax=2,ymin=-Inf,ymax=Inf,fill="grey88") +
  geom_line(linewidth=.7) +
  geom_point(data=data.frame(z=hq,f=fb),aes(z,f),inherit.aes=FALSE,size=2,shape=24,fill="black") +
  scale_y_continuous(trans=scales::pseudo_log_trans(sigma=1e-4,base=10)) +
  labs(title="d  forms fitted to the band (grey), projected to canopy",
       x="height (m)", y=NULL, colour=NULL) + theme_bw(base_size=8) +
  theme(legend.text=element_text(size=6), legend.key.height=unit(8,"pt"))

ggsave(file.path(outdir,"fig_diag1_height_profiles.png"),
  (p1a|p1b)/(p1c|p1d) + plot_annotation(title="Height-flux profiles: observed, learned, and modelled",
   theme=theme(plot.title=element_text(size=10,face="bold"))),
  width=10,height=7,dpi=200,bg="white")

# ==============================================================================
# FIG 2  surface area with height
# ==============================================================================
p2a <- ggplot(INV,aes(H)) + geom_histogram(bins=60,fill="#4393C3",colour=NA) +
  geom_vline(xintercept=median(INV$H),linetype="22") +
  labs(title=sprintf("a  stem heights from allometry (median %.1f m)",median(INV$H)),
       x="height (m)", y="stems") + theme_bw(base_size=8)

Hc <- 25; zc <- seq(.05,Hc,length.out=400); zcb <- .5*Hc
A_bole <- sum(pi*INV$dbh*INV$H/2); A_tot <- 3.07*A_PLOT
dens <- list(uniform=rep(1,length(zc)), conic=pmax(Hc-zc,0),
  `bole+crown`={b<-pmax(Hc-zc,0);b<-b/sum(b)*A_bole
                c_<-as.numeric(zc>=zcb);c_<-c_/sum(c_)*(A_tot-A_bole); b+c_})
dens <- lapply(dens,function(v) v/max(v))
p2b <- ggplot(bind_rows(lapply(names(dens),function(n) data.frame(z=zc,a=dens[[n]],s=n))),
  aes(a,z,colour=s)) + annotate("rect",xmin=-Inf,xmax=Inf,ymin=0,ymax=2,fill="grey88") +
  geom_path(linewidth=.8) + labs(title="b  per-tree area density assumptions",
   x="relative density", y="height (m)", colour=NULL) + theme_bw(base_size=8) +
  theme(legend.text=element_text(size=6), legend.key.height=unit(8,"pt"))

# stand-level emergent distribution: each tree spreads its area over 0..H_i
zg <- seq(0,ceiling(max(INV$H)),by=.25)
Atree <- pi*INV$dbh*INV$H/2
stand <- sapply(zg, function(zz0) sum(Atree[INV$H>zz0]/INV$H[INV$H>zz0]))
sdf <- data.frame(z=zg, a=stand/max(stand))
p2c <- ggplot(sdf,aes(a,z)) + annotate("rect",xmin=-Inf,xmax=Inf,ymin=0,ymax=2,fill="grey88") +
  geom_path(colour="#B2182B",linewidth=.9) +
  labs(title="c  STAND-LEVEL woody area vs height", x="relative area density", y="height (m)",
       subtitle="bottom-weighted even though each tree is uniform,\nbecause stems differ in height") +
  theme_bw(base_size=8)

cum <- sapply(zg,function(zz0) sum(Atree*pmin(zz0,INV$H)/INV$H))/sum(Atree)
p2d <- ggplot(data.frame(z=zg,c=cum),aes(z,c)) + geom_line(colour="#2166AC",linewidth=.9) +
  geom_vline(xintercept=2,linetype="22",colour="grey40") +
  annotate("text",x=2.4,y=.06,label=sprintf("%.1f%% below 2 m",100*approx(zg,cum,2)$y),hjust=0,size=2.6) +
  scale_y_continuous(labels=scales::percent) +
  labs(title="d  cumulative woody area below a height", x="height (m)", y="share of woody area") +
  theme_bw(base_size=8)

ggsave(file.path(outdir,"fig_diag2_area_height.png"), (p2a|p2b)/(p2c|p2d) +
  plot_annotation(title="Vertical distribution of woody surface area",
   theme=theme(plot.title=element_text(size=10,face="bold"))),
  width=10,height=7,dpi=200,bg="white")

# ==============================================================================
# FIG 3  where the budget comes from
# ==============================================================================
adens <- approx(zg,c(stand[1],diff(cum)*sum(Atree)/diff(zg)[1]),xout=zz)$y
adens <- sapply(zz,function(z0) sum(Atree[INV$H>z0]/INV$H[INV$H>z0]))
fa <- bind_rows(lapply(names(FORMS),function(n)
  data.frame(z=zz, v=f2b*FORMS[[n]]*adens, form=n)))
p3a <- ggplot(fa,aes(z,v,colour=form)) +
  annotate("rect",xmin=0,xmax=2,ymin=-Inf,ymax=Inf,fill="grey88") +
  geom_line(linewidth=.7) +
  labs(title="a  flux x area density vs height", x="height (m)",
       y="contribution density", colour=NULL) + theme_bw(base_size=8) +
  theme(legend.text=element_text(size=6), legend.key.height=unit(8,"pt"))
cumf <- fa %>% group_by(form) %>% arrange(z) %>%
  mutate(c=cumsum(v)/sum(v)) %>% ungroup()
p3b <- ggplot(cumf,aes(z,c,colour=form)) + geom_line(linewidth=.7) +
  geom_vline(xintercept=2,linetype="22",colour="grey40") +
  scale_y_continuous(labels=scales::percent) +
  labs(title="b  cumulative share of the whole-surface budget",
       subtitle="everything right of the dashed line is extrapolated",
       x="height (m)", y="share of budget", colour=NULL) + theme_bw(base_size=8) +
  theme(legend.text=element_text(size=6), legend.key.height=unit(8,"pt"))
ggsave(file.path(outdir,"fig_diag3_budget_by_height.png"), p3a/p3b +
  plot_annotation(title="Where the whole-surface budget comes from",
   theme=theme(plot.title=element_text(size=10,face="bold"))),
  width=8,height=7,dpi=200,bg="white")

cat("\nShare of the whole-surface budget coming from ABOVE 2 m (extrapolated):\n")
for (n in names(FORMS)) {
  cc <- cumf %>% filter(form==n)
  cat(sprintf("  %-18s %.1f%%\n", n, 100*(1-approx(cc$z,cc$c,2)$y)))
}
cat("\nWritten: fig_diag1_height_profiles.png, fig_diag2_area_height.png, fig_diag3_budget_by_height.png\n")
