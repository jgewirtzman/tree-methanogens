#!/usr/bin/env Rscript
# ==============================================================================
# rev_area_distribution_scenarios.R
# ------------------------------------------------------------------------------
# HOW IS WOODY SURFACE AREA DISTRIBUTED WITH HEIGHT?
#
# This is an unconstrained assumption -- no LiDAR, no branch harvest at this site --
# and it was previously left implicit. The scaling had been assuming PER-TREE UNIFORM
# area density from 0 to each stem's height. That is documented here as one scenario
# among several rather than carried silently.
#
# The area of each stem is split into a BOLE and a BRANCH component, each with its own
# vertical shape, then summed across all stems. Because stem heights differ, the
# STAND-LEVEL distribution is bottom-weighted even when every individual tree is
# uniform -- only tall stems contribute area high up. That emergent effect is real and
# falls out of the allometry; it is not itself an assumption.
#
# BOLE SHAPES  (density of surface area with height, from the solid of revolution)
#   cylinder    d(z) = const          a(z) ∝ 1            no taper
#   cone        d(z) ∝ (H-z)          a(z) ∝ (H-z)        Whittaker & Woodwell conic
#   paraboloid  d(z) ∝ sqrt(H-z)      a(z) ∝ sqrt(H-z)    W&W parabolic form; they note
#                                                          true stems lie BETWEEN the
#                                                          conic and parabolic estimates
#
# BRANCH SHAPES
#   uniform_all     branch area spread evenly over 0 to H
#   uniform_top50   spread evenly over the top half, [0.5H, H]
#   uniform_crown   spread evenly over [0.6H, H]
#   gaussian_75     Gaussian centred at 0.75H, sd 0.15H
#
#   NONE of these is sourced. Branch area in real crowns is neither uniform nor
#   Gaussian; these bracket plausible shapes. The Gaussian is included because it is a
#   commonly assumed crown form, not because it is supported here.
#
# AREA SPLIT   bole  = pi*DBH*H/2 * taper_corr   (W&W conic surface, corrected 1.20-1.35)
#              branch = bole * branch:stem ratio (W&W: 2.7 Smokies deciduous, 4.0
#                       Brookhaven oak-pine)
#
# Output: outputs/revision/area_distribution_scenarios.csv / .txt
#         outputs/revision/fig_area_distribution_scenarios.png
# ==============================================================================

suppressMessages({library(ranger);library(dplyr);library(ggplot2);library(tidyr);library(patchwork)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
source("code/revision/rev_geometry.R")
source("code/revision/rev_species_levels.R")
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
d <- tree_train_complete; trained <- sort(unique(as.character(d$species_clean)))
CONV <- 86400*365.25*16e-6
A_PLOT <- STAND_AREA_M2   # censused ground. Was hardcoded 200*200 = 40,000, which
                          # is the NOMINAL square and 7.5% too large, so every WAI
                          # and every mg m-2 this script reported was low by that.

# CANONICAL INPUTS, not private re-derivations. This script used to carry its own
# copy of the whole upscaling front end, and every part of it had drifted:
#   * rf_workflow_data$PLACEHOLDER_INVENTORY -- the pre-rebuild inventory, before
#     the bytag coordinate rescue and the DBH unit fix (rev_inventory_build.R);
#   * an fx() stack of species-specific DBH "repairs" that rev_inventory_build.R
#     replaced with one documented rule, and which over-corrected;
#   * the campaign monthly drivers with approx(rule = 2) gapfilling, rather than
#     the climatologies;
#   * an UNCALIBRATED 2 m anchor predicted here from scratch;
#   * the species-ladder bug (see rev_species_levels.R).
# It now reads the same products rev_scaling_full_grid.R reads, so the two cannot
# disagree. One producer per quantity.
NEED <- c(inv="outputs/tables/inventory_stems.csv",
          tree="outputs/tables/tree_flux_predictions.csv")
if (any(!file.exists(NEED)))
  stop("missing ", paste(NEED[!file.exists(NEED)], collapse=", "),
       " -- run rev_inventory_build.R and rev_predict_tree_flux_current.R first")
INV <- read.csv(NEED[["inv"]],  stringsAsFactors=FALSE)
TRP <- read.csv(NEED[["tree"]], stringsAsFactors=FALSE)
INV <- INV[INV$in_stand, ]; TRP <- TRP[TRP$in_stand, ]
stopifnot(nrow(INV) == nrow(TRP))
INV$dbh <- INV$dbh_m          # already unit-checked and typo-repaired upstream
INV$sp  <- species_to_model_level(INV$species, trained)
GY <- c("Pinus strobus","Tsuga canadensis")
DA <- quantile(INV$dbh[INV$dbh>0.10],.95,na.rm=TRUE)
# NO shrub cap: the published allometry is applied uniformly (see scaling_parameters.md)
INV$H <- 1.37 + ifelse(INV$species %in% GY,(25-1.37)/DA^0.60,(25-1.37)/DA^0.53) *
         INV$dbh^ifelse(INV$species %in% GY,0.60,0.53)

TAPER <- 1.275     # midpoint of W&W 1.20-1.35
BRSTEM <- 3.35     # midpoint of 2.7-4.0
INV$A_bole   <- pi*INV$dbh*INV$H/2 * TAPER
INV$A_branch <- INV$A_bole * BRSTEM

# ---- 2 m anchor per tree: READ the canonical calibrated value ----------------
# This was re-predicted here from the campaign drivers with no per-species
# calibration, so it disagreed with the anchor the budget and the grid use.
f2 <- TRP$flux_2m_nmol_m2_s
stopifnot(all(is.finite(f2)))

# ---- vertical shapes --------------------------------------------------------
BOLE <- list(
  cylinder   = function(u) rep(1,length(u)),
  cone       = function(u) pmax(1-u,0),
  paraboloid = function(u) sqrt(pmax(1-u,0)))
BRANCH <- list(
  uniform_all   = function(u) rep(1,length(u)),
  uniform_top50 = function(u) as.numeric(u >= 0.50),
  uniform_crown = function(u) as.numeric(u >= 0.60),
  gaussian_75   = function(u) dnorm(u, 0.75, 0.15))

# flux forms, as ratios to the 2 m value
hq <- c(.5,1.25,2)
fbar <- sapply(hq,function(z) mean(d$stem_flux_corrected[
  !is.na(d$measurement_height_cm) & abs(d$measurement_height_cm/100-z)<1e-6], na.rm=TRUE))
b_log <- coef(lm(log(pmax(fbar,1e-9))~hq))[2]
b_pow <- coef(lm(log(pmax(fbar,1e-9))~log(hq)))[2]
FLUX <- list(constant=function(z) rep(1,length(z)),
             `asymptote 50%`=function(z) 0.5+0.5*exp(b_log*(z-2)),
             power=function(z) exp(b_pow*(log(pmax(z,.05))-log(2))),
             `log-linear`=function(z) exp(b_log*(z-2)))

# ---- evaluate: integrate each tree's own profile, sum over the stand ---------
NU <- 90; u <- (seq_len(NU)-0.5)/NU          # relative height within each tree
res <- list(); prof_store <- list()
for (bn in names(BOLE)) for (rn in names(BRANCH)) {
  wb <- BOLE[[bn]](u); wb <- wb/sum(wb)
  wr <- BRANCH[[rn]](u); wr <- wr/sum(wr)
  # per-tree absolute height of each slab, and the area in it
  Z <- outer(INV$H, u)                        # n x NU
  Aslab <- outer(INV$A_bole, wb) + outer(INV$A_branch, wr)
  frac_above2 <- sum(Aslab[Z>2])/sum(Aslab)
  row <- data.frame(bole=bn, branch=rn, WAI=sum(Aslab)/A_PLOT,
                    frac_area_above_2m=frac_above2,
                    mean_area_height=sum(Aslab*Z)/sum(Aslab))
  for (fn in names(FLUX)) {
    R <- matrix(FLUX[[fn]](as.vector(Z)), nrow=nrow(Z))
    row[[fn]] <- sum(f2 * rowSums(Aslab*R))/A_PLOT*CONV
  }
  res[[paste(bn,rn)]] <- row
  # stand-level density for plotting
  zg <- seq(0,ceiling(max(INV$H)),by=.25)
  dens <- sapply(seq_along(zg)[-1], function(k)
    sum(Aslab[Z>=zg[k-1] & Z<zg[k]]))
  prof_store[[paste(bn,rn)]] <- data.frame(z=zg[-1], a=dens/max(dens),
                                           bole=bn, branch=rn)
}
R <- bind_rows(res)
cat(sprintf("
Stand: %d stems, median height %.1f m, WAI %.2f (taper %.3f, branch:stem %.2f)
Per-tree area = bole (pi*DBH*H/2 x taper) + branch (bole x branch:stem).
'frac_area_above_2m' is the share of woody area OUTSIDE the measured band.\n",
  nrow(INV), median(INV$H), R$WAI[1], TAPER, BRSTEM))
cat("\n=== AREA DISTRIBUTION SCENARIOS ===\n\n")
print(R %>% mutate(across(where(is.numeric), ~round(.x,3))), row.names=FALSE)
write.csv(R, file.path(outdir,"area_distribution_scenarios.csv"), row.names=FALSE)

cat(sprintf("
  Share of woody area above 2 m ranges %.3f - %.3f across shapes.
  Whole-surface budget (mg CH4 m-2 yr-1):
    under CONSTANT flux      %.1f - %.1f   (spread %.2fx)  <- distribution barely matters
    under LOG-LINEAR decline %.1f - %.1f   (spread %.1fx)  <- distribution dominates
", min(R$frac_area_above_2m), max(R$frac_area_above_2m),
   min(R$constant), max(R$constant), max(R$constant)/min(R$constant),
   min(R$`log-linear`), max(R$`log-linear`), max(R$`log-linear`)/min(R$`log-linear`)))
cat("
  The constant-flux column is IDENTICAL across all twelve shapes, exactly. Under a
  constant flux the budget is f(2m) x total area, and total area is fixed by the
  allometry, taper and branch:stem ratio -- the vertical distribution cancels
  algebraically, it is not merely insensitive. Under any declining form it does not
  cancel, and the branch shape then dominates: branch area is ~77% of the total
  (branch:stem 3.35), so the bole shape matters far less than where branch area sits.
")

# ---- figure ------------------------------------------------------------------
P <- bind_rows(prof_store)
p1 <- ggplot(P, aes(a, z, colour=branch)) +
  annotate("rect",xmin=-Inf,xmax=Inf,ymin=0,ymax=2,fill="grey88") +
  geom_path(linewidth=.6) + facet_wrap(~bole, nrow=1) +
  labs(title="a  stand-level woody area density, by bole and branch shape",
       x="relative area density", y="height (m)", colour="branch shape") +
  theme_bw(base_size=8) + theme(legend.position="bottom",
    legend.text=element_text(size=6), legend.key.height=unit(7,"pt"))
RL <- R %>% pivot_longer(c(constant,`asymptote 50%`,power,`log-linear`),
                         names_to="flux", values_to="mg") %>%
  mutate(flux=factor(flux, levels=c("constant","asymptote 50%","power","log-linear")),
         combo=paste(bole, branch, sep=" / "))
p2 <- ggplot(RL, aes(reorder(combo,mg), mg, fill=branch)) +
  geom_col() + coord_flip() + facet_wrap(~flux, nrow=1, scales="free_x") +
  labs(title="b  whole-surface budget under each area shape x flux form",
       x=NULL, y=expression(mg~CH[4]~m^-2~yr^-1), fill="branch shape") +
  theme_bw(base_size=7) + theme(legend.position="bottom",
    legend.text=element_text(size=6), legend.key.height=unit(7,"pt"))
ggsave(file.path(outdir,"fig_area_distribution_scenarios.png"), p1/p2 +
  plot_layout(heights=c(1,1.15)) +
  plot_annotation(title="Vertical distribution of woody surface area: scenarios",
    subtitle="no shape is supported by data from this site; these bracket plausible geometries",
    theme=theme(plot.title=element_text(size=10,face="bold"),
                plot.subtitle=element_text(size=8))),
  width=11,height=8.5,dpi=200,bg="white")
cat("\nWritten: outputs/revision/fig_area_distribution_scenarios.png\n")
