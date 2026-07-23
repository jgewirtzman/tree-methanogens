# ==============================================================================
# REVISION — Per-regression figures so each specification is legible
# ==============================================================================
# For each predictor (ratio, mcrA, methanotroph) produces ONE 6-panel figure:
#   rows   = flux scale (raw ; arcsinh/pseudo-log)
#   cols   = M1 individual (x & y individual)
#            M3 aggregate-x, individual chamber y  (Mark's model)
#            M2 species medians (both axes; OLS red, SMA blue dashed)
# Each panel titled with slope / R2 / p. Lets Jon see WHY raw fails (skew) and
# how aggregation changes the fit. NEW file; edits nothing.
#
# Run: Rscript revision/analysis/R_regression_figures.R
# Outputs: revision/outputs/regfig_<predictor>.png  (+ combined contact sheet)
# ==============================================================================

source("revision/analysis/_prep_species_data.R")
suppressPackageStartupMessages({ library(ggplot2); library(gridExtra) })
out_dir <- "revision/outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

psl <- function(x) asinh(x / 0.2) / log(10)
sma_fit <- function(x, y) { ok<-is.finite(x)&is.finite(y); x<-x[ok]; y<-y[ok]
  r<-cor(x,y); list(r=r,r2=r^2,p=cor.test(x,y)$p.value,
    b_ols=unname(coef(lm(y~x))[2]), b_sma=sign(r)*sd(y)/sd(x)) }

fit_lab <- function(m) { s<-summary(m); sprintf("slope=%.3f  R2=%.3f  p=%.3f",
  coef(m)[2], s$r.squared, coef(s)[2,4]) }

flux_pool <- flux_all %>% dplyr::select(species, CH4_flux)

# predictor configs: individual x column + species-level x
cfgs <- list(
  list(key="ratio", lab="mcrA:methanotroph ratio",
       xind=tree_level_complete %>% transmute(species, x=log_ratio, CH4_flux),
       xsp = analysis_ratio %>% transmute(species, x_sp=median_log_ratio),
       xlab="log10 mcrA:methanotroph"),
  list(key="mcra", lab="mcrA (area-weighted)",
       xind=tree_level_complete %>% transmute(species, x=log_mcra, CH4_flux),
       xsp = analysis_mcra %>% transmute(species, x_sp=log10(value+1)),
       xlab="log10 mcrA"),
  list(key="meth", lab="methanotroph (pmoA+mmoX)",
       xind=tree_level_complete %>% transmute(species, x=log_meth, CH4_flux),
       xsp = analysis_meth %>% transmute(species, x_sp=log10(value+1)),
       xlab="log10 methanotroph")
)

make_row <- function(cfg, ytrans, scale_name, ylab) {
  ty <- ytrans
  ind <- cfg$xind %>% mutate(y = ty(CH4_flux))
  cham <- flux_pool %>% inner_join(cfg$xsp, by="species") %>% mutate(y = ty(CH4_flux))
  spm <- flux_pool %>% inner_join(cfg$xsp, by="species") %>%
    group_by(species, x_sp) %>% summarise(y=median(ty(CH4_flux)), .groups="drop")

  pA <- ggplot(ind, aes(x, y)) + geom_point(alpha=.35, size=1) +
    geom_smooth(method="lm", se=FALSE, color="firebrick", formula=y~x) +
    labs(title=paste0("M1 individual  [", scale_name, "]"),
         subtitle=fit_lab(lm(y~x, ind)), x=cfg$xlab, y=ylab) + theme_bw(base_size=9)
  pB <- ggplot(cham, aes(x_sp, y)) +
    geom_point(alpha=.30, size=1, position=position_jitter(width=.02)) +
    geom_smooth(method="lm", se=FALSE, color="firebrick", formula=y~x) +
    labs(title=paste0("M3 aggregate-x, individual y  [", scale_name, "]"),
         subtitle=fit_lab(lm(y~x_sp, cham)), x=paste(cfg$xlab,"(species median)"), y=NULL) +
    theme_bw(base_size=9)
  f <- sma_fit(spm$x_sp, spm$y)
  pC <- ggplot(spm, aes(x_sp, y)) + geom_point(size=2.3) +
    geom_abline(intercept=mean(spm$y)-f$b_ols*mean(spm$x_sp), slope=f$b_ols, color="firebrick") +
    geom_abline(intercept=mean(spm$y)-f$b_sma*mean(spm$x_sp), slope=f$b_sma, color="steelblue", linetype="dashed") +
    labs(title=paste0("M2 species medians  [", scale_name, "]"),
         subtitle=sprintf("OLS %.3f / SMA %.3f  R2=%.3f  p=%.3f", f$b_ols,f$b_sma,f$r2,f$p),
         x=paste(cfg$xlab,"(species median)"), y=ylab) + theme_bw(base_size=9)
  list(pA,pB,pC)
}

for (cfg in cfgs) {
  raw <- make_row(cfg, identity, "raw", "CH4 flux (nmol/m2/s)")
  pl  <- make_row(cfg, psl,      "pseudo-log", "pseudo-log10 CH4 flux")
  g <- arrangeGrob(grobs=c(raw, pl), nrow=2,
                   top=grid::textGrob(cfg$lab, gp=grid::gpar(fontsize=13, fontface="bold")))
  ggsave(file.path(out_dir, paste0("regfig_", cfg$key, ".png")), g, width=13, height=8, dpi=140)
  cat("wrote regfig_", cfg$key, ".png\n", sep="")
}
cat("Done.\n")
