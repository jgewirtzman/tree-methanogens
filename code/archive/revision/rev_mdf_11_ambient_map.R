#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_11_ambient_map.R      STEP 11
# ------------------------------------------------------------------------------
# ALTERNATIVE HYPOTHESIS: is elevated chamber C0 real local atmosphere?
#
# Steps 9-10 showed that apparent stem uptake scales with how far the chamber's
# starting concentration sits above LOCAL ambient. That was interpreted as chamber
# carryover. A competing explanation is that the air itself is genuinely enriched over
# wet, CH4-emitting ground -- in which case elevated C0 is a real property of the site,
# not a measurement fault, and the negatives might reflect something biological.
#
# THE TWO HYPOTHESES MAKE DIFFERENT PREDICTIONS
#   carryover        excess is a property of the MEASUREMENT (what happened just
#                    before it) and should be spatially unstructured
#   real atmosphere  ambient is a property of the PLACE, so it should be spatially
#                    coherent and should track soil moisture, which drives the surface
#                    source
#
# Note the excess metric already differences against local ambient measured minutes
# either side, so a spatially elevated background should largely cancel. The test is
# therefore run on BOTH quantities: absolute local ambient (does it map?) and excess
# (does it map?).
#
# CONTROLLING FOR CHANGING BACKGROUND
#   Ambient CH4 varies with synoptic conditions between days and with boundary-layer
#   development through the day. Both are removed before any spatial or moisture
#   analysis by taking residuals from  ambient ~ factor(date) + hour + hour^2.
#   Without this, a day when sampling happened to run through a wet patch in the early
#   morning would masquerade as a moisture effect.
#
# Output: outputs/revision/ambient_map.txt
#         outputs/revision/fig_ambient_map.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

L <- read.csv("outputs/revision/local_ambient_check.csv", stringsAsFactors=FALSE)
a <- read.csv("data/processed/flux/ymf2023_goflux_auxfile.csv", stringsAsFactors=FALSE,
              check.names=TRUE)
cat(sprintf("local-ambient records: %d\n", nrow(L)))

# --- attach VWC, temperature, tag, time --------------------------------------
vw <- grep("^VWC", names(a), value=TRUE)
cat(sprintf("VWC columns found: %s\n", paste(vw, collapse=", ")))
a$vwc <- if (length(vw)) rowMeans(sapply(vw, function(c) suppressWarnings(as.numeric(a[[c]]))),
                                  na.rm=TRUE) else NA_real_
tagc <- grep("Tree.Tag|Tree_Tag|tree_tag", names(a), value=TRUE)[1]
stc  <- grep("Soil.Temp|Soil_Temp", names(a), value=TRUE)[1]
a$st <- as.POSIXct(a$start.time, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
A <- a %>% transmute(UniqueID, tag=as.character(.data[[tagc]]), vwc,
                     soil_temp=if(!is.na(stc)) suppressWarnings(as.numeric(.data[[stc]])) else NA_real_,
                     st, date=as.Date(st), hour=as.numeric(format(st,"%H")) +
                       as.numeric(format(st,"%M"))/60)
D <- L %>% inner_join(A, by="UniqueID")
cat(sprintf("joined: %d | with VWC: %d | with tag: %d\n",
            nrow(D), sum(is.finite(D$vwc)), sum(!is.na(D$tag))))

# --- attach inventory coordinates ---------------------------------------------
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
INV <- rf_workflow_data$PLACEHOLDER_INVENTORY
INV$tag <- as.character(INV$tree_id)
D <- D %>% left_join(INV %>% dplyr::select(tag, x, y, PX, PY, species) %>% distinct(tag,.keep_all=TRUE),
                     by="tag")
cat(sprintf("with coordinates: %d\n", sum(is.finite(D$x))))

# ==============================================================================
rule("1. REMOVING THE CHANGING BACKGROUND")
# ==============================================================================
cat(sprintf("
  raw local ambient: median %.1f ppb, range %.1f - %.1f, SD %.1f
  between-day SD of the daily mean: %.1f ppb
", median(D$local_ambient), min(D$local_ambient), max(D$local_ambient), sd(D$local_ambient),
   sd(tapply(D$local_ambient, D$date, mean))))
m <- lm(local_ambient ~ factor(date) + hour + I(hour^2), data=D)
D$amb_resid <- resid(m)
cat(sprintf("  ambient ~ date + hour + hour^2 :  R2 = %.3f  (that much is background, not place)\n",
            summary(m)$r.squared))
cat(sprintf("  residual SD after detrending: %.1f ppb\n", sd(D$amb_resid)))

# ==============================================================================
rule("2. DOES AMBIENT TRACK SOIL MOISTURE?")
# ==============================================================================
k <- is.finite(D$vwc)
if (sum(k) > 10) {
  for (v in c("local_ambient","amb_resid","excess")) {
    ct <- suppressWarnings(cor.test(D$vwc[k], D[[v]][k], method="spearman"))
    cat(sprintf("  %-16s vs VWC : Spearman rho %+.3f  p = %.4g   (n=%d)\n",
                v, ct$estimate, ct$p.value, sum(k)))
  }
  cat("\n  by VWC tercile:\n\n")
  D$vt <- cut(D$vwc, quantile(D$vwc, c(0,1/3,2/3,1), na.rm=TRUE),
              labels=c("dry","mid","wet"), include.lowest=TRUE)
  print(D %>% filter(!is.na(vt)) %>% group_by(vt) %>%
    summarise(n=n(), vwc=round(median(vwc),1), ambient=round(median(local_ambient),1),
              ambient_detrended=round(median(amb_resid),1), excess=round(median(excess),1),
              pct_neg=round(100*mean(flux<0),1), .groups="drop") %>% as.data.frame(), row.names=FALSE)
} else cat("  insufficient VWC data\n")

# ==============================================================================
rule("3. IS AMBIENT SPATIALLY COHERENT?")
# ==============================================================================
kk <- is.finite(D$x) & is.finite(D$y)
if (sum(kk) > 20) {
  MLON <- 111320*cos(41.99*pi/180); MLAT <- 111320
  P <- D[kk,]; P$xm <- P$x*MLON; P$ym <- P$y*MLAT
  # Moran-style check: correlation of pairwise value difference with distance
  n <- nrow(P); idx <- expand.grid(i=1:n, j=1:n); idx <- idx[idx$i<idx$j,]
  dd <- sqrt((P$xm[idx$i]-P$xm[idx$j])^2 + (P$ym[idx$i]-P$ym[idx$j])^2)
  for (v in c("amb_resid","excess")) {
    vv <- abs(P[[v]][idx$i] - P[[v]][idx$j])
    ct <- suppressWarnings(cor.test(dd, vv, method="spearman"))
    cat(sprintf("  |difference in %-10s| vs pairwise distance : rho %+.3f  p = %.3g\n",
                v, ct$estimate, ct$p.value))
  }
  cat("
  A POSITIVE rho means nearby measurements are more alike than distant ones, i.e. the
  quantity is spatially structured -- the signature of a real property of the place.
  A rho near zero means no spatial structure, which is what carryover would look like.
")
  cat("\n  repeatability within trees measured more than once:\n")
  rp <- P %>% group_by(tag) %>% filter(n()>=3) %>%
    summarise(n=n(), sd_within=sd(amb_resid), .groups="drop")
  if (nrow(rp)) cat(sprintf("    %d trees with >=3 visits | median within-tree SD %.1f ppb vs overall SD %.1f ppb\n",
      nrow(rp), median(rp$sd_within), sd(P$amb_resid)))
}

# ==============================================================================
rule("4. FIGURE")
# ==============================================================================
if (sum(kk) > 20) {
  P <- D[kk,]
  pa <- ggplot(P, aes(x, y, colour=amb_resid, size=abs(excess))) + geom_point(alpha=.85) +
    scale_colour_gradient2(low="#2166AC", mid="grey90", high="#B2182B", midpoint=0,
                           name="ambient\n(detrended, ppb)") +
    scale_size_continuous(range=c(.7,4), name="|excess|") +
    coord_quickmap() +
    labs(title="a  detrended local ambient across the plot", x=NULL, y=NULL) +
    theme_bw(base_size=8) + theme(axis.text=element_text(size=5))
  pb <- ggplot(P[is.finite(P$vwc),], aes(vwc, amb_resid)) +
    geom_point(aes(colour=flux<0), size=1.3, alpha=.8) +
    geom_smooth(method="lm", se=TRUE, colour="grey25", linewidth=.5, formula=y~x) +
    scale_colour_manual(values=c(`FALSE`="grey60",`TRUE`="#B2182B"),
                        labels=c("positive","negative"), name="flux") +
    labs(title="b  detrended ambient vs soil moisture", x="VWC (%)", y="ambient residual (ppb)") +
    theme_bw(base_size=8) + theme(legend.position="bottom", legend.text=element_text(size=6))
  pc <- ggplot(P[is.finite(P$vwc),], aes(vwc, excess)) +
    geom_point(aes(colour=flux<0), size=1.3, alpha=.8) +
    geom_smooth(method="lm", se=TRUE, colour="grey25", linewidth=.5, formula=y~x) +
    scale_y_continuous(trans=scales::pseudo_log_trans(sigma=5, base=10)) +
    scale_colour_manual(values=c(`FALSE`="grey60",`TRUE`="#B2182B"), guide="none") +
    labs(title="c  excess vs soil moisture", x="VWC (%)", y="C0 - local ambient (ppb)") +
    theme_bw(base_size=8)
  pd <- ggplot(D, aes(hour, local_ambient, colour=factor(date))) +
    geom_point(size=.8, alpha=.7) +
    labs(title="d  raw ambient by hour and day (the background being removed)",
         x="hour of day", y="local ambient (ppb)") +
    theme_bw(base_size=8) + theme(legend.position="none")
  ggsave(file.path(outdir,"fig_ambient_map.png"), (pa|pb)/(pc|pd) +
    plot_annotation(title="Is elevated chamber C0 real local atmosphere?",
      subtitle="detrended for date and time of day before spatial / moisture analysis",
      theme=theme(plot.title=element_text(size=11,face="bold"))),
    width=11, height=8, dpi=200, bg="white")
  cat("  Written: outputs/revision/fig_ambient_map.png\n")
}
write.csv(D, file.path(outdir,"ambient_map_data.csv"), row.names=FALSE)
