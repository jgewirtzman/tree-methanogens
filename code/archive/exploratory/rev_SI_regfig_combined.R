# ==============================================================================
# REVISION — SI: combined gene-flux regression figure (mcrA + methanotroph).
# Merges regfig_mcra + regfig_meth into ONE figure: two labelled blocks, each
# 6 panels = aggregation level (individual / aggregate-x / species-median)
# x flux scale (raw / arcsinh pseudo-log). Shows the relationship emerges at the
# species-median scale for BOTH genes (mcrA +, methanotroph -).
# Reuses helpers/data from R_regression_figures.R. NEW file. Edits nothing.
# Output: outputs/revision/SI_regfig_combined.png
# ==============================================================================
source("code/revision/rev_prep_species_data.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(gridExtra) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

psl <- function(x) asinh(x / 0.2) / log(10)
sma_fit <- function(x, y) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  r <- cor(x, y); list(r = r, r2 = r^2, p = cor.test(x, y)$p.value,
    b_ols = unname(coef(lm(y ~ x))[2]), b_sma = sign(r)*sd(y)/sd(x)) }
fit_lab <- function(m) { s <- summary(m); sprintf("slope=%.3f  R2=%.3f  p=%.3f",
  coef(m)[2], s$r.squared, coef(s)[2,4]) }
flux_pool <- flux_all %>% dplyr::select(species, CH4_flux)

cfgs <- list(
  mcra = list(key="mcra", lab="mcrA (area-weighted)",
       xind = tree_level_complete %>% transmute(species, x = log_mcra, CH4_flux),
       xsp  = analysis_mcra %>% transmute(species, x_sp = log10(value+1)),
       xlab = "log10 mcrA"),
  meth = list(key="meth", lab="methanotroph (pmoA+mmoX)",
       xind = tree_level_complete %>% transmute(species, x = log_meth, CH4_flux),
       xsp  = analysis_meth %>% transmute(species, x_sp = log10(value+1)),
       xlab = "log10 methanotroph"))

make_row <- function(cfg, ytrans, scale_name, ylab) {
  ty <- ytrans
  ind  <- cfg$xind %>% mutate(y = ty(CH4_flux))
  cham <- flux_pool %>% inner_join(cfg$xsp, by = "species") %>% mutate(y = ty(CH4_flux))
  spm  <- flux_pool %>% inner_join(cfg$xsp, by = "species") %>%
    group_by(species, x_sp) %>% summarise(y = median(ty(CH4_flux)), .groups = "drop")
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
  list(pA, pB, pC)
}

block <- function(cfg) arrangeGrob(
  grobs = c(make_row(cfg, identity, "raw", "CH4 flux (nmol/m2/s)"),
            make_row(cfg, psl, "pseudo-log", "pseudo-log10 CH4 flux")),
  nrow = 2, top = grid::textGrob(cfg$lab, gp = grid::gpar(fontsize = 13, fontface = "bold")))

g <- arrangeGrob(block(cfgs$mcra), block(cfgs$meth), ncol = 1)
ggsave(file.path(out, "SI_regfig_combined.png"), g, width = 12, height = 15, dpi = 150)
cat("Wrote SI_regfig_combined.png\n")
