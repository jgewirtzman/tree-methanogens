#!/usr/bin/env Rscript
# ==============================================================================
# rev_extreme_outliers.R
# ------------------------------------------------------------------------------
# BIOLOGICALLY IMPOSSIBLE VALUES: identify and quantify.
#
# This is a Tier 1 (hard QC) question, distinct from detection limits. Tier 1 removes
# measurements that cannot be real, on physical grounds. It does NOT apply statistical
# thresholds to plausible flux values -- doing that biases the mean (see
# ch4-data-filtering, sec 5.3.4).
#
# Distributions are examined SEPARATELY BY MEASUREMENT TYPE, because soil and tree
# stems have completely different expected ranges and signs. Pooling them makes soil
# uptake look like a stem artefact, and vice versa.
#
# PLAUSIBILITY BOUNDS FROM THE LITERATURE (nmol CH4 m-2 s-1)
#   upland forest soil uptake   typically -0.1 to -3; extreme well-drained soils to -10
#   upland tree stem emission   typically 0 to a few; hotspots to tens
#   upland tree stem uptake     small, order -0.1
# Values far outside these are flagged as candidates for hard removal, with the
# decision left explicit rather than automatic.
#
# Output: outputs/revision/extreme_outliers.txt
#         outputs/revision/fig_extreme_outliers.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

sr <- read.csv("data/compiled/semirigid_chamber_flux.csv", stringsAsFactors=FALSE)
st <- read.csv("data/compiled/static_chamber_flux.csv", stringsAsFactors=FALSE)
D <- bind_rows(
  sr %>% transmute(src=paste0("semirigid/", measurement_type), flux=CH4_best_flux_nmol_m2_s,
                   r2=CH4_LM_r2, qual=CH4_quality, id=unique_id, date=as.character(date)),
  st %>% transmute(src="static/tree_stem", flux=CH4_best_flux_nmol_m2_s,
                   r2=CH4_LM_r2, qual=CH4_quality, id=unique_id, date=as.character(date)))
D <- D[is.finite(D$flux),]

cat("=== DISTRIBUTION BY MEASUREMENT TYPE ===\n\n")
qs <- c(0,.001,.01,.05,.25,.5,.75,.95,.99,.999,1)
for (s in unique(D$src)) {
  v <- D$flux[D$src==s]
  cat(sprintf("%-24s n = %d\n", s, length(v)))
  print(round(quantile(v, qs), 4)); cat("\n")
}

# --- Tukey fences on each type, plus literature plausibility ------------------
cat("=== EXTREME OUTLIERS: Tukey far-out fence (Q1/Q3 -/+ 3*IQR) ===\n\n")
cat(sprintf("  %-24s %8s %10s %10s %8s %8s\n","type","n","lower fence","upper fence","n below","n above"))
FEN <- D %>% group_by(src) %>% summarise(n=n(),
  lo=quantile(flux,.25)-3*IQR(flux), hi=quantile(flux,.75)+3*IQR(flux),
  n_lo=sum(flux < quantile(flux,.25)-3*IQR(flux)),
  n_hi=sum(flux > quantile(flux,.75)+3*IQR(flux)), .groups="drop")
for (i in seq_len(nrow(FEN))) with(FEN[i,],
  cat(sprintf("  %-24s %8d %10.4f %10.4f %8d %8d\n", src, n, lo, hi, n_lo, n_hi)))

LIT <- list(`semirigid/soil`=c(-10, 2), `semirigid/tree_stem`=c(-1, 100),
            `static/tree_stem`=c(-1, 100))
cat("\n=== BEYOND LITERATURE-PLAUSIBLE RANGE ===\n\n")
cat(sprintf("  %-24s %14s %8s %s\n","type","plausible range","n outside","values outside"))
flagged <- list()
for (s in names(LIT)) {
  v <- D[D$src==s,]; b <- LIT[[s]]
  out <- v[v$flux < b[1] | v$flux > b[2],]
  flagged[[s]] <- out
  cat(sprintf("  %-24s %14s %8d  %s\n", s, sprintf("[%.0f, %.0f]", b[1], b[2]), nrow(out),
      if (nrow(out)) paste(sprintf("%.2f", sort(out$flux)), collapse=", ") else "--"))
}

cat("\n=== THE FLAGGED MEASUREMENTS IN DETAIL ===\n\n")
FL <- bind_rows(flagged)
if (nrow(FL)) {
  FL <- FL %>% arrange(flux)
  cat(sprintf("  %-26s %10s %7s %-22s %s\n","id","flux","LM r2","quality","date"))
  for (i in seq_len(nrow(FL))) with(FL[i,],
    cat(sprintf("  %-26s %10.3f %7.3f %-22s %s\n", substr(id,1,26), flux,
        ifelse(is.na(r2),NA,r2), substr(ifelse(is.na(qual),"",qual),1,22), date)))
} else cat("  none\n")

cat("\n=== EFFECT OF REMOVING THEM ===\n\n")
for (s in names(LIT)) {
  v <- D$flux[D$src==s]; b <- LIT[[s]]; k <- v >= b[1] & v <= b[2]
  cat(sprintf("  %-24s mean %10.4f -> %8.4f (%+.1f%%)   median %8.4f -> %8.4f   n %d -> %d\n",
      s, mean(v), mean(v[k]), 100*(mean(v[k])/mean(v)-1), median(v), median(v[k]), length(v), sum(k)))
}
cat("
  NOTE the asymmetry: removing a handful of impossible values changes the MEAN
  substantially where they exist, and the MEDIAN not at all. That is the signature of
  a genuine outlier problem rather than a heavy-tailed but real distribution, and it is
  why Tier 1 removal is justified here while statistical filtering of plausible values
  is not.
")

p <- ggplot(D, aes(flux, src, fill=src)) +
  geom_violin(scale="width", linewidth=.2, alpha=.7) +
  geom_point(data=FL, aes(flux, src), inherit.aes=FALSE, colour="#B2182B", size=2, shape=4, stroke=1) +
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=.01, base=10)) +
  labs(title="CH4 flux distributions by measurement type",
       subtitle="pseudo-log axis; red crosses = values outside the literature-plausible range",
       x=expression(CH[4]~flux~(nmol~m^-2~s^-1)), y=NULL) +
  theme_bw(base_size=9) + theme(legend.position="none")
ggsave(file.path(outdir,"fig_extreme_outliers.png"), p, width=9, height=4.5, dpi=200, bg="white")
cat("\nWritten: outputs/revision/fig_extreme_outliers.png\n")
