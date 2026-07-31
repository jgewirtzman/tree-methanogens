#!/usr/bin/env Rscript
# ==============================================================================
# rev_negative_flux_quality.R
# ------------------------------------------------------------------------------
# ARE ANY OF OUR NEGATIVE STEM FLUXES REAL?
#
# 8.6% of stem flux measurements are negative, and they are NOT randomly distributed:
# essentially all occur in the June and July campaigns, with none at all in eight
# other months. Before a negative value is allowed to set a bound on the vertical
# extrapolation (the "decline through zero into uptake" scenario), it needs a quality
# assessment.
#
# Three criteria are available, though not for every campaign:
#   MDF        minimal detectable flux, from goFlux (precision / t x flux.term).
#              Available only for the 2021 height campaign, which is the only campaign
#              reprocessed through goFlux with the diagnostic columns retained.
#   LM r2      linear-model fit quality, available for all campaigns.
#   quality.check / CH4_quality   goFlux diagnostic flags.
#
# A negative flux is treated as CREDIBLE only if it is (a) more negative than -MDF
# where MDF exists, and (b) supported by an acceptable model fit.
#
# Output: outputs/revision/negative_flux_quality.txt
#         outputs/revision/fig_negative_flux_quality.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

# ==============================================================================
rule("1. HEIGHT CAMPAIGN (2021) -- the only data with MDF")
# ==============================================================================
m <- read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE)
m$height_cm <- suppressWarnings(as.numeric(sub(".*_(\\d+)\\s*(cm)?_\\d+_\\d+$","\\1",m$UniqueID)))
m <- m[is.finite(m$best.flux),]
cat(sprintf("
  n = %d measurements
  precision used by goFlux : %.2f ppb  (MANUFACTURER spec, constant across rows)
  MDF        median %.4f   IQR %.4f - %.4f   nmol m-2 s-1
  |flux|     median %.4f
  ratio |flux| / MDF, median %.1f  -- typical measurement is well above detection
", nrow(m), unique(m$prec)[1], median(m$MDF), quantile(m$MDF,.25), quantile(m$MDF,.75),
   median(abs(m$best.flux)), median(abs(m$best.flux)/m$MDF)))

cat(sprintf("
  BELOW DETECTION      |flux| < MDF : %3d (%.1f%%)
  NEGATIVE                 flux < 0 : %3d (%.1f%%)
  SIGNIFICANTLY NEGATIVE  flux < -MDF: %3d (%.1f%%)
", sum(abs(m$best.flux)<m$MDF), 100*mean(abs(m$best.flux)<m$MDF),
   sum(m$best.flux<0), 100*mean(m$best.flux<0),
   sum(m$best.flux < -m$MDF), 100*mean(m$best.flux < -m$MDF)))

sig <- m[m$best.flux < -m$MDF,]
cat(sprintf("\n  Of the %d significantly-negative measurements, model fit quality:\n", nrow(sig)))
cat(sprintf("    LM r2 > 0.7 : %d\n    LM r2 < 0.3 : %d\n    with a goFlux quality flag : %d\n",
    sum(sig$LM.r2>0.7,na.rm=TRUE), sum(sig$LM.r2<0.3,na.rm=TRUE),
    sum(nchar(trimws(sig$quality.check))>0,na.rm=TRUE)))
cat(sprintf("    most negative CREDIBLE (exceeds MDF and r2 > 0.7): %s\n",
    { z <- sig$best.flux[sig$LM.r2>0.7]; if(length(z)) sprintf("%.4f", min(z)) else "NONE" }))
cat(sprintf("    most negative of any kind in this campaign: %.4f\n", min(m$best.flux)))

cat("\n  By measurement height:\n")
hh <- m %>% filter(height_cm %in% c(50,125,200)) %>% group_by(height_cm) %>%
  summarise(n=n(), n_neg=sum(best.flux<0), pct_neg=round(100*mean(best.flux<0),1),
            n_sig_neg=sum(best.flux < -MDF), below_MDF=sum(abs(best.flux)<MDF),
            min_flux=round(min(best.flux),4), .groups="drop")
print(as.data.frame(hh), row.names=FALSE)

# ==============================================================================
rule("2. ALL CAMPAIGNS -- where do the negatives come from?")
# ==============================================================================
grab <- function(f, fluxcol, r2col, qcol, lab) {
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors=FALSE)
  if (!fluxcol %in% names(d)) return(NULL)
  data.frame(campaign=lab, flux=d[[fluxcol]],
             r2=if(r2col %in% names(d)) d[[r2col]] else NA_real_,
             qual=if(qcol %in% names(d)) as.character(d[[qcol]]) else NA_character_,
             date=if("date" %in% names(d)) as.character(d$date) else NA_character_,
             stringsAsFactors=FALSE)
}
A <- bind_rows(
  grab("data/compiled/semirigid_chamber_flux.csv","CH4_best_flux_nmol_m2_s","CH4_LM_r2","CH4_quality","semirigid"),
  grab("data/compiled/static_chamber_flux.csv","CH4_best_flux_nmol_m2_s","CH4_LM_r2","CH4_quality","static"))
A <- A[is.finite(A$flux),]
cat("\n  Negatives by campaign source:\n\n")
cat(sprintf("  %-12s %6s %8s %7s %10s %10s\n","campaign","n","n_neg","%neg","min flux","median neg"))
for (cp in unique(A$campaign)) {
  q <- A[A$campaign==cp,]; ng <- q$flux[q$flux<0]
  cat(sprintf("  %-12s %6d %8d %6.1f%% %10.4f %10.4f\n", cp, nrow(q), length(ng),
      100*length(ng)/nrow(q), min(q$flux), if(length(ng)) median(ng) else NA))
}
cat("\n  Fit quality of NEGATIVE measurements vs POSITIVE, by campaign:\n\n")
cat(sprintf("  %-12s %-9s %6s %10s %10s\n","campaign","sign","n","median r2","%% r2>0.7"))
for (cp in unique(A$campaign)) for (s in c("negative","positive")) {
  q <- A[A$campaign==cp & (if(s=="negative") A$flux<0 else A$flux>=0),]
  if (!nrow(q)) next
  cat(sprintf("  %-12s %-9s %6d %10.3f %9.0f%%\n", cp, s, nrow(q),
      median(q$r2,na.rm=TRUE), 100*mean(q$r2>0.7,na.rm=TRUE)))
}

# the extreme value
ext <- A[which.min(A$flux),]
cat(sprintf("\n  MOST NEGATIVE MEASUREMENT IN THE WHOLE DATASET
    campaign %s | flux %.4f | LM r2 %.3f | quality '%s' | date %s
", ext$campaign, ext$flux, ext$r2, ext$qual, ext$date))

# ==============================================================================
rule("3. IS THERE A DEFENSIBLE NEGATIVE FLOOR?")
# ==============================================================================
credible <- m$best.flux[m$best.flux < -m$MDF & m$LM.r2 > 0.7]
A_cred <- A$flux[A$flux < 0 & is.finite(A$r2) & A$r2 > 0.7]
cat(sprintf("
  Candidate floors, most to least conservative:

    most negative measurement of any quality          %.4f
    most negative with LM r2 > 0.7 (all campaigns)    %.4f
    most negative exceeding MDF with r2 > 0.7         %s
    5th percentile of all negative measurements       %.4f
    median negative measurement                       %.4f
",
 min(A$flux,na.rm=TRUE),
 if(length(A_cred)) min(A_cred) else NA_real_,
 if(length(credible)) sprintf("%.4f", min(credible)) else "NONE EXIST",
 quantile(A$flux[A$flux<0], .05, na.rm=TRUE),
 median(A$flux[A$flux<0], na.rm=TRUE)))

# --- figure ------------------------------------------------------------------
p1 <- ggplot(m, aes(abs(best.flux), MDF)) +
  geom_point(aes(colour=best.flux<0), size=1.1, alpha=.7) +
  geom_abline(slope=1, intercept=0, linetype="22") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=c("FALSE"="grey55","TRUE"="#B2182B"),
                      labels=c("positive","negative"), name=NULL) +
  labs(title="a  |flux| vs minimal detectable flux (height campaign)",
       subtitle="points below the 1:1 line are detectable; above it, indistinguishable from zero",
       x=expression("|"*CH[4]*" flux| (nmol m"^-2*" s"^-1*")"), y="MDF") +
  theme_bw(base_size=8)
p2 <- ggplot(A[A$flux<0,], aes(flux, fill=campaign)) +
  geom_histogram(bins=40) +
  labs(title="b  distribution of negative fluxes by campaign",
       x=expression(CH[4]~flux~(nmol~m^-2~s^-1)), y="count", fill=NULL) +
  theme_bw(base_size=8)
p3 <- ggplot(A, aes(r2, fill=flux<0)) + geom_histogram(bins=30, position="identity", alpha=.6) +
  geom_vline(xintercept=0.7, linetype="22") +
  scale_fill_manual(values=c("FALSE"="grey55","TRUE"="#B2182B"),
                    labels=c("positive","negative"), name=NULL) +
  labs(title="c  model fit quality, negative vs positive measurements",
       x=expression(LM~r^2), y="count") + theme_bw(base_size=8)
ggsave(file.path(outdir,"fig_negative_flux_quality.png"), (p1|p2)/p3 +
  plot_annotation(title="Quality assessment of negative stem CH4 fluxes",
    theme=theme(plot.title=element_text(size=10,face="bold"))),
  width=10,height=7,dpi=200,bg="white")
cat("\nWritten: outputs/revision/fig_negative_flux_quality.png\n")
