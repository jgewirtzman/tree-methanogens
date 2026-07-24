# ==============================================================================
# REVISION — S11 rebuilt: 4-column scale-dependent gene-flux figure.
# Columns = aggregation ladder + a median-vs-mean robustness check:
#   M1 Individual         : individual x (gene), individual y (flux)
#   M3 Aggregate predictor: species-MEDIAN x, INDIVIDUAL y                 [bridge]
#   M2 Species-median     : species-MEDIAN x & y
#   M2m Species-mean      : species-MEAN x & y                             [robustness]
# Rows = 5 predictors (mcrA, pmoA, mmoX, pmoA+mmoX, ratio). Bottom row = R2 bars.
#
# Every panel reports slope + R2 + p from a GENE-ONLY simple model fit to the plotted
# response (so stats match the drawn line, and R2 bars = the gene's R2, not the old
# full species+gene model). M1/M3 responses on the pseudolog individual-flux axis;
# M2/M2m on species median/mean flux. The two species columns test whether the
# aggregation result is robust to central-tendency choice (median vs mean).
#
# Sources the original generator for data/theme/palette/pseudolog only; builds ALL
# panels fresh. ORIGINAL UNTOUCHED; its two re-rendered outputs are git-restored by
# the wrapper. Tags a..t (row-major) + u,v,w,x bars.
# Output: revision/outputs/figS11_final.png
# ==============================================================================
suppressWarnings(suppressMessages(
  source("code/05_gene_flux_analysis/02_scale_dependent_gene_patterns.R")
))
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(grid) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

PRED <- tibble(
  key   = c("mcrA", "pmoA", "mmoX", "methanotroph", "ratio"),
  label = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  ix    = c("log_tree_mcra", "log_pmoa", "log_tree_mmox", "log_methanotroph", "log_ratio"),
  col   = c("#E74C3C", "#3498DB", "#5DADE2", "#1F618D", "#9B59B6"))
xlab_i <- list(mcrA=expression(log[10]~mcrA), pmoA=expression(log[10]~pmoA), mmoX=expression(log[10]~mmoX), methanotroph=expression(log[10]~(pmoA+mmoX)), ratio=expression(log[10]~ratio))
xlab_md<- list(mcrA=expression(log[10]~median~mcrA), pmoA=expression(log[10]~median~pmoA), mmoX=expression(log[10]~median~mmoX), methanotroph=expression(log[10]~median~(pmoA+mmoX)), ratio=expression(log[10]~median~ratio))
xlab_mn<- list(mcrA=expression(log[10]~mean~mcrA), pmoA=expression(log[10]~mean~pmoA), mmoX=expression(log[10]~mean~mmoX), methanotroph=expression(log[10]~mean~(pmoA+mmoX)), ratio=expression(log[10]~mean~ratio))

# species-MEDIAN predictor + flux (from the original analysis_* objects)
A <- list(
  mcrA=analysis_mcra %>% mutate(sx=log10(median_mcra+1)), pmoA=analysis_pmoa %>% mutate(sx=log10(median_pmoa+1)),
  mmoX=analysis_mmox %>% mutate(sx=log10(median_mmox+1)), methanotroph=analysis_methanotroph %>% mutate(sx=log10(median_methanotroph+1)),
  ratio=analysis_ratio %>% mutate(sx=median_log_ratio))
spx <- lapply(A, function(x) x %>% transmute(species, sx))

# species-MEAN predictor + flux (built the same way, mean instead of median)
flux_mean <- flux_all %>% group_by(species, species_id) %>% summarise(mean_flux=mean(CH4_flux, na.rm=TRUE), n_flux=n(), .groups="drop")
gm <- tree_level_complete %>% group_by(species, species_id) %>% summarise(
  n_trees=n(), mcrA=mean(mcrA, na.rm=TRUE), pmoA=mean(pmoA, na.rm=TRUE), mmoX=mean(mmoX, na.rm=TRUE),
  methanotroph=mean(methanotroph_total, na.rm=TRUE), ratio=mean(log_ratio, na.rm=TRUE), .groups="drop") %>%
  inner_join(flux_mean, by=c("species","species_id")) %>% filter(n_trees>=5, n_flux>=5)
Am <- list(
  mcrA=gm %>% transmute(species, mx=log10(mcrA+1), mean_flux), pmoA=gm %>% transmute(species, mx=log10(pmoA+1), mean_flux),
  mmoX=gm %>% transmute(species, mx=log10(mmoX+1), mean_flux), methanotroph=gm %>% transmute(species, mx=log10(methanotroph+1), mean_flux),
  ratio=gm %>% transmute(species, mx=ratio, mean_flux))

mk <- function(k, level, tag, ylab="") {
  col <- PRED$col[match(k, PRED$key)]
  if (level=="ind") { d <- tree_level_complete %>% mutate(xv=.data[[PRED$ix[match(k,PRED$key)]]], yv=pseudolog10_individual(CH4_flux)) %>% filter(is.finite(xv), is.finite(yv)); xlab<-xlab_i[[k]] }
  else if (level=="agg") { d <- tree_level_complete %>% inner_join(spx[[k]], by="species") %>% mutate(xv=sx, yv=pseudolog10_individual(CH4_flux)) %>% filter(is.finite(xv), is.finite(yv)); xlab<-xlab_md[[k]] }
  else if (level=="sp")  { d <- A[[k]]  %>% mutate(xv=sx, yv=median_flux) %>% filter(is.finite(xv), is.finite(yv)); xlab<-xlab_md[[k]] }
  else                   { d <- Am[[k]] %>% mutate(xv=mx, yv=mean_flux)   %>% filter(is.finite(xv), is.finite(yv)); xlab<-xlab_mn[[k]] }
  f <- lm(yv ~ xv, d); s<-coef(f)[2]; r2<-summary(f)$r.squared; p<-summary(f)$coef[2,4]
  g <- ggplot(d, aes(xv, yv, color=species)) +
    geom_point(size=if (level %in% c("sp","spm")) 3.1 else 1.9, alpha=if (level %in% c("sp","spm")) 0.85 else 0.6) +
    geom_smooth(aes(group=1), method="lm", formula=y~x, se=TRUE, color=col, fill=col, alpha=0.15, linewidth=1) +
    annotate("label", x=Inf, y=Inf, label=sprintf("slope=%.3f\nR2=%.3f  p=%.3f", s, r2, p), hjust=1.04, vjust=1.04, size=2.5, fill="white", alpha=0.9, label.size=0.15, lineheight=0.92) +
    scale_color_manual(values=species_palette, name="Species") +
    labs(tag=tag, x=xlab, y=ylab) + theme_pub_gene + theme(legend.position="none")
  if (level %in% c("sp","spm")) g <- g + geom_hline(yintercept=0, linetype="dashed", color="gray50")
  else g <- g + geom_hline(yintercept=pseudolog10_individual(0), linetype="dashed", color="gray50") +
    scale_y_continuous(breaks=pseudolog10_individual(c(0,0.1,1)), labels=c("0","0.1","1"))
  list(panel=g, r2=r2, p=p)
}

tags <- matrix(c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)",
                 "(m)","(n)","(o)","(p)","(q)","(r)","(s)","(t)"), ncol=4, byrow=TRUE)
levs <- c("ind","agg","sp","spm")
ylabs <- list(ind=y_lab_individual, agg="", sp=expression(Median~CH[4]~flux~(nmol~m^{-2}~s^{-1})), spm=expression(Mean~CH[4]~flux~(nmol~m^{-2}~s^{-1})))
res <- list()
for (r in 1:5) for (cc in 1:4) { k<-PRED$key[r]; lv<-levs[cc]; res[[paste(k,lv)]] <- mk(k, lv, tags[r,cc], ylab=ylabs[[lv]]) }

# ---- R2 bars (gene-only) per column, shared y-limit --------------------------
bar_df <- function(lv) data.frame(Model=PRED$label, R2=sapply(PRED$key, function(k) res[[paste(k,lv)]]$r2),
  P=sapply(PRED$key, function(k) res[[paste(k,lv)]]$p)) %>% mutate(Significant=P<0.05, Model=factor(Model, levels=Model[order(R2)]))
D <- lapply(levs, bar_df); names(D) <- levs
y_lim <- max(sapply(D, function(x) max(x$R2))) * 1.2
bar_panel <- function(dat, tag, ylab="") ggplot(dat, aes(Model, R2, fill=Significant)) +
  geom_col(alpha=0.85) + geom_text(aes(label=sprintf("%.3f", R2)), vjust=-0.3, size=2.7) +
  scale_fill_manual(values=c("FALSE"="gray70","TRUE"="#27AE60"), labels=c("NS","p < 0.05"), name="") +
  ylim(0, y_lim) + labs(tag=tag, x="", y=ylab) + theme_pub_gene +
  theme(legend.position="none", axis.text.x=element_text(angle=35, hjust=1, size=7.5))
bars <- list(bar_panel(D$ind,"(u)",expression(italic(R)^2~"(gene only)")), bar_panel(D$agg,"(v)"), bar_panel(D$sp,"(w)"), bar_panel(D$spm,"(x)"))

# ---- assemble 4 columns ------------------------------------------------------
hh <- function(txt) wrap_elements(full=textGrob(txt, gp=gpar(fontface="bold", fontsize=11)))
H <- list(hh("Individual\n(individual x, individual y)"), hh("Aggregate predictor\n(species-median x, individual y)"),
          hh("Species-median\n(median x & y)"), hh("Species-mean\n(mean x & y)"))
hts <- c(0.13, 1, 1, 1, 1, 1, 0.92)
P <- function(k, lv) res[[paste(k, lv)]]$panel
mkcol <- function(cc) H[[cc]] / P("mcrA",levs[cc]) / P("pmoA",levs[cc]) / P("mmoX",levs[cc]) / P("methanotroph",levs[cc]) / P("ratio",levs[cc]) / bars[[cc]] + plot_layout(ncol=1, heights=hts)
combined <- (mkcol(1) | mkcol(2) | mkcol(3) | mkcol(4)) / legend_panel + plot_layout(heights=c(1, 0.05))
ggsave(file.path(out, "figS11_final.png"), combined, width=21, height=16, dpi=300, bg="white")
for (lv in levs) cat(sprintf("%-4s gene R2: %s\n", lv, paste(sprintf("%s=%.3f", D[[lv]]$Model, D[[lv]]$R2), collapse="  ")))
cat("Wrote figS11_final.png\n")
