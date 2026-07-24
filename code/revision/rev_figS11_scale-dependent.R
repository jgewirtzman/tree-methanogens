# ==============================================================================
# REVISION — S11 rebuilt: 4-column scale-dependent gene-flux figure.
# ALL columns on a common arcsinh flux scale (asinh(x/0.1)/log10; cofactor matches Fig 1/2/3),
# and aggregation is TRANSFORM-THEN-AGGREGATE on both axes, so the median-vs-mean
# comparison reflects real aggregation, not raw right-skew sensitivity.
#   M1 Individual         : individual log-gene, individual arcsinh-flux
#   M3 Aggregate predictor: species-MEDIAN log-gene, INDIVIDUAL arcsinh-flux    [bridge]
#   M2 Species-median     : per species, MEDIAN of (log-gene, arcsinh-flux)
#   M2m Species-mean      : per species, MEAN of (log-gene, arcsinh-flux)       [robustness]
# Rows = 5 predictors (mcrA, pmoA, mmoX, pmoA+mmoX, ratio). Bottom row = R2 bars.
# Every panel reports slope + R2 + p from a GENE-ONLY simple model on the plotted
# (arcsinh) response, so stats match the drawn line and R2 bars are the gene's R2.
#
# Sources the original generator for data/theme/palette/pseudolog only; builds ALL
# panels fresh. ORIGINAL UNTOUCHED; its two re-rendered outputs git-restored by the
# wrapper. Tags a..t (row-major) + u,v,w,x bars.
# Output: outputs/revision/figS11_final.png
# ==============================================================================
suppressWarnings(suppressMessages(
  source("code/05_gene_flux_analysis/02_scale_dependent_gene_patterns.R")
))
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(grid) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
COFAC <- 0.1                             # arcsinh cofactor — matches Fig 1/2/3; conclusions robust across 0.05-2 (sweep)
asf <- function(x) asinh(x / COFAC) / log(10)
ybrk <- asf(c(0, 0.1, 1))

PRED <- tibble(
  key   = c("mcrA", "pmoA", "mmoX", "methanotroph", "ratio"),
  label = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  lg    = c("lg_mcra", "lg_pmoa", "lg_mmox", "lg_meth", "lg_ratio"),
  col   = c("#E74C3C", "#3498DB", "#5DADE2", "#1F618D", "#9B59B6"))
xlab_i  <- list(mcrA=expression(log[10]~mcrA), pmoA=expression(log[10]~pmoA), mmoX=expression(log[10]~mmoX), methanotroph=expression(log[10]~(pmoA+mmoX)), ratio=expression(log[10]~ratio))
xlab_md <- list(mcrA=expression(median~log[10]~mcrA), pmoA=expression(median~log[10]~pmoA), mmoX=expression(median~log[10]~mmoX), methanotroph=expression(median~log[10]~(pmoA+mmoX)), ratio=expression(median~log[10]~ratio))
xlab_mn <- list(mcrA=expression(mean~log[10]~mcrA), pmoA=expression(mean~log[10]~pmoA), mmoX=expression(mean~log[10]~mmoX), methanotroph=expression(mean~log[10]~(pmoA+mmoX)), ratio=expression(mean~log[10]~ratio))

# per-individual log-gene + arcsinh flux
TL <- tree_level_complete %>% mutate(
  lg_mcra=log10(mcrA+1), lg_pmoa=log10(pmoA+1), lg_mmox=log10(mmoX+1),
  lg_meth=log10(methanotroph_total+1), lg_ratio=log_ratio, af=asf(CH4_flux))
# species aggregates: transform THEN aggregate (median & mean of each transformed value)
gene_agg <- TL %>% group_by(species, species_id) %>%
  summarise(across(c(lg_mcra, lg_pmoa, lg_mmox, lg_meth, lg_ratio), list(med=~median(.x, na.rm=TRUE), mean=~mean(.x, na.rm=TRUE)), .names="{.col}_{.fn}"),
            n_trees=n(), .groups="drop")
flux_agg <- flux_all %>% mutate(af=asf(CH4_flux)) %>% group_by(species, species_id) %>%
  summarise(med_af=median(af, na.rm=TRUE), mean_af=mean(af, na.rm=TRUE), n_flux=n(), .groups="drop")
spp <- gene_agg %>% inner_join(flux_agg, by=c("species","species_id")) %>% filter(n_trees>=5, n_flux>=5)

mk <- function(k, level, tag, ylab="") {
  col <- PRED$col[match(k, PRED$key)]; lg <- PRED$lg[match(k, PRED$key)]
  if (level=="ind")      { d <- TL %>% transmute(species, xv=.data[[lg]], yv=af); xlab<-xlab_i[[k]] }
  else if (level=="agg") { lk <- spp %>% transmute(species, xv=.data[[paste0(lg,"_med")]]); d <- TL %>% transmute(species, yv=af) %>% inner_join(lk, by="species"); xlab<-xlab_md[[k]] }
  else if (level=="sp")  { d <- spp %>% transmute(species, xv=.data[[paste0(lg,"_med")]],  yv=med_af);  xlab<-xlab_md[[k]] }
  else                   { d <- spp %>% transmute(species, xv=.data[[paste0(lg,"_mean")]], yv=mean_af); xlab<-xlab_mn[[k]] }
  d <- d %>% filter(is.finite(xv), is.finite(yv))
  f <- lm(yv ~ xv, d); s<-coef(f)[2]; r2<-summary(f)$r.squared; p<-summary(f)$coef[2,4]
  ggplot(d, aes(xv, yv, color=species)) +
    geom_point(size=if (level %in% c("sp","spm")) 3.1 else 1.9, alpha=if (level %in% c("sp","spm")) 0.85 else 0.6) +
    geom_smooth(aes(group=1), method="lm", formula=y~x, se=TRUE, color=col, fill=col, alpha=0.15, linewidth=1) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50") +
    annotate("label", x=Inf, y=Inf, label=sprintf("slope=%.3f\nR2=%.3f  p=%.3f", s, r2, p), hjust=1.04, vjust=1.04, size=2.5, fill="white", alpha=0.9, label.size=0.15, lineheight=0.92) +
    scale_y_continuous(breaks=ybrk, labels=c("0","0.1","1")) +
    scale_color_manual(values=species_palette, name="Species") +
    labs(tag=tag, x=xlab, y=ylab) + theme_pub_gene + theme(legend.position="none") -> g
  list(panel=g, r2=r2, p=p)
}

tags <- matrix(c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)",
                 "(m)","(n)","(o)","(p)","(q)","(r)","(s)","(t)"), ncol=4, byrow=TRUE)
levs <- c("ind","agg","sp","spm")
ylabs <- list(ind=expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})), agg="",
              sp=expression(Median~CH[4]~flux), spm=expression(Mean~CH[4]~flux))
res <- list()
for (r in 1:5) for (cc in 1:4) { k<-PRED$key[r]; lv<-levs[cc]; res[[paste(k,lv)]] <- mk(k, lv, tags[r,cc], ylab=ylabs[[lv]]) }

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

hh <- function(txt) wrap_elements(full=textGrob(txt, gp=gpar(fontface="bold", fontsize=11)))
H <- list(hh("Individual\n(individual x, individual y)"), hh("Aggregate predictor\n(species-median x, individual y)"),
          hh("Species-median\n(median of arcsinh flux & log-gene)"), hh("Species-mean\n(mean of arcsinh flux & log-gene)"))
hts <- c(0.13, 1, 1, 1, 1, 1, 0.92)
P <- function(k, lv) res[[paste(k, lv)]]$panel
mkcol <- function(cc) H[[cc]] / P("mcrA",levs[cc]) / P("pmoA",levs[cc]) / P("mmoX",levs[cc]) / P("methanotroph",levs[cc]) / P("ratio",levs[cc]) / bars[[cc]] + plot_layout(ncol=1, heights=hts)
combined <- (mkcol(1) | mkcol(2) | mkcol(3) | mkcol(4)) / legend_panel + plot_layout(heights=c(1, 0.05))
ggsave(file.path(out, "figS11_final.png"), combined, width=21, height=16, dpi=300, bg="white")
for (lv in levs) cat(sprintf("%-4s gene R2: %s\n", lv, paste(sprintf("%s=%.3f(p%.3f)", D[[lv]]$Model, D[[lv]]$R2, D[[lv]]$P), collapse="  ")))
cat("Wrote figS11_final.png\n")
