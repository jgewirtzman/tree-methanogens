# ==============================================================================
# figS21_scale-dependent-genes.R -- Figure S21: gene-flux relationships across scales
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
# Four columns, one flux scale (arcsinh, cofactor 0.1, log10 units, as Figs 1-3);
# every aggregate is transform-then-aggregate on both axes, so median vs mean
# reflects aggregation rather than raw right skew.
#   M1  Individual          individual log-gene, individual arcsinh flux
#   M3  Aggregate predictor species-median log-gene, individual arcsinh flux
#   M2  Species median      per species, median of (log-gene, arcsinh flux)
#   M2m Species mean        per species, mean of (log-gene, arcsinh flux)
# Rows: mcrA, pmoA, mmoX, pmoA+mmoX, ratio; bottom row: R2 bars. Each panel reports
# slope, R2 and p of the gene-only model on the plotted response.
# Data: code/07_molecular/helper_scale_dependent_gene_patterns.R (tables only).
# Output: outputs/figures/generated/figS11_final.png (assembled as Figure S21)
# ==============================================================================
source("code/lib/outputs.R")
suppressWarnings(suppressMessages(
  source("code/07_molecular/helper_scale_dependent_gene_patterns.R")
))
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(grid) })
source("code/lib/fig_style.R")      # SI house style
source("code/lib/flux_transform.R")   # flux_asinh10(), cofactor 0.1 (conclusions robust across 0.05-2)
ybrk <- flux_asinh10(c(0, 0.1, 1))
si_panel <- theme(legend.position="none", plot.tag=element_text(size=SI_BASE+1, face="bold"),
                  axis.title.x=element_text(size=rel(0.9)), plot.margin=margin(2, 4, 2, 2))

# about four x ticks per small panel, at least three inside the range
xbrk <- function(lims) {
  b <- scales::breaks_extended(n=4)(lims)
  if (sum(b >= lims[1] & b <= lims[2]) < 3) b <- scales::breaks_extended(n=5)(lims)
  b
}

PRED <- tibble(
  key   = c("mcrA", "pmoA", "mmoX", "methanotroph", "ratio"),
  label = c("mcrA", "pmoA", "mmoX", "pmoA+mmoX", "Ratio"),
  lg    = c("lg_mcra", "lg_pmoa", "lg_mmox", "lg_meth", "lg_ratio"),
  col   = c("#C03221", "#4A6FA5", "#7FA7D0", "#2C4A73", "#6B5B95"))   # as in Figure 8
xlab_i  <- list(mcrA=expression(log[10]~italic(mcrA)), pmoA=expression(log[10]~italic(pmoA)), mmoX=expression(log[10]~italic(mmoX)), methanotroph=expression(log[10]*(italic(pmoA)+italic(mmoX))), ratio=expression(log[10]~ratio))
# medians / means: the column header names the statistic, so every column uses the same x titles
xlab_md <- xlab_i; xlab_mn <- xlab_i

# per-individual log-gene + arcsinh flux
TL <- tree_level_complete %>% mutate(
  lg_mcra=log10(mcrA+1), lg_pmoa=log10(pmoA+1), lg_mmox=log10(mmoX+1),
  lg_meth=log10(methanotroph_total+1), lg_ratio=log_ratio, af=flux_asinh10(CH4_flux))
# species aggregates: transform THEN aggregate (median & mean of each transformed value)
gene_agg <- TL %>% group_by(species, species_id) %>%
  summarise(across(c(lg_mcra, lg_pmoa, lg_mmox, lg_meth, lg_ratio), list(med=~median(.x, na.rm=TRUE), mean=~mean(.x, na.rm=TRUE)), .names="{.col}_{.fn}"),
            n_trees=n(), .groups="drop")
flux_agg <- flux_all %>% mutate(af=flux_asinh10(CH4_flux)) %>% group_by(species, species_id) %>%
  summarise(med_af=median(af, na.rm=TRUE), mean_af=mean(af, na.rm=TRUE), n_flux=n(), .groups="drop")
spp <- gene_agg %>% inner_join(flux_agg, by=c("species","species_id")) %>% filter(n_trees>=5, n_flux>=5)

mk <- function(k, level, tag, ylab=NULL) {
  col <- PRED$col[match(k, PRED$key)]; lg <- PRED$lg[match(k, PRED$key)]
  if (level=="ind")      { d <- TL %>% transmute(species, xv=.data[[lg]], yv=af); xlab<-xlab_i[[k]] }
  else if (level=="agg") { lk <- spp %>% transmute(species, xv=.data[[paste0(lg,"_med")]]); d <- TL %>% transmute(species, yv=af) %>% inner_join(lk, by="species"); xlab<-xlab_md[[k]] }
  else if (level=="sp")  { d <- spp %>% transmute(species, xv=.data[[paste0(lg,"_med")]],  yv=med_af);  xlab<-xlab_md[[k]] }
  else                   { d <- spp %>% transmute(species, xv=.data[[paste0(lg,"_mean")]], yv=mean_af); xlab<-xlab_mn[[k]] }
  d <- d %>% filter(is.finite(xv), is.finite(yv))
  f <- lm(yv ~ xv, d); s<-coef(f)[2]; r2<-summary(f)$r.squared; p<-summary(f)$coef[2,4]
  ggplot(d, aes(xv, yv)) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50", linewidth=0.4) +
    geom_point(size=if (level %in% c("sp","spm")) 1.8 else 1, alpha=if (level %in% c("sp","spm")) 0.9 else 0.45,
               colour="grey35", stroke=0) +
    geom_smooth(aes(group=1), method="lm", formula=y~x, se=TRUE, color=col, fill=col, alpha=0.15, linewidth=0.8) +
    scale_y_continuous(breaks=ybrk, labels=c("0","0.1","1")) +
    scale_x_continuous(breaks=xbrk, labels=function(x) sub("^-", "\u2212", as.character(signif(x, 6)))) +
    labs(tag=tag, x=xlab, y=ylab) + theme_si() + si_panel -> g
  list(panel=g, r2=r2, p=p)
}

tags <- matrix(c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)",
                 "(m)","(n)","(o)","(p)","(q)","(r)","(s)","(t)"), ncol=4, byrow=TRUE)
levs <- c("ind","agg","sp","spm")
# one y title for the five gene rows (collected by patchwork); the column headers say whether
# flux is per tree, a species median or a species mean
ylab_flux <- expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))
res <- list()
for (r in 1:5) for (cc in 1:4) { k<-PRED$key[r]; lv<-levs[cc]; res[[paste(k,lv)]] <- mk(k, lv, tags[r,cc], ylab=if (cc==1) ylab_flux else NULL) }

bar_df <- function(lv) data.frame(Model=PRED$label, R2=sapply(PRED$key, function(k) res[[paste(k,lv)]]$r2),
  P=sapply(PRED$key, function(k) res[[paste(k,lv)]]$p)) %>% mutate(Significant=P<0.05, Model=factor(Model, levels=Model[order(R2)]))
D <- lapply(levs, bar_df); names(D) <- levs
y_lim <- max(sapply(D, function(x) max(x$R2))) * 1.2
bar_panel <- function(dat, tag, ylab=NULL) ggplot(dat, aes(Model, R2, fill=Significant)) +
  geom_col(width=0.75) + geom_text(aes(label=sprintf("%.2f", R2)), vjust=-0.3, size=2.8) +
  scale_fill_manual(values=c("FALSE"="grey75","TRUE"="#285238"), labels=c("Not significant","p < 0.05"), name=NULL) +
  scale_x_discrete(labels=c(mcrA=expression(italic(mcrA)), pmoA=expression(italic(pmoA)), mmoX=expression(italic(mmoX)),
                            `pmoA+mmoX`=expression(italic(pmoA)+italic(mmoX)), Ratio="Ratio")) +
  scale_y_continuous(limits=c(0, y_lim), expand=expansion(mult=c(0, 0.02))) +
  labs(tag=tag, x=NULL, y=ylab) + theme_si() + si_panel +
  theme(axis.text.x=element_text(angle=40, hjust=1, vjust=1), panel.grid.major.x=element_blank())
bars <- list(bar_panel(D$ind,"(u)",expression(R^2~"(gene only)")), bar_panel(D$agg,"(v)"), bar_panel(D$sp,"(w)"), bar_panel(D$spm,"(x)"))

# column headers: plain sentence case above the top row
hh <- function(txt) wrap_elements(full=textGrob(txt, gp=gpar(fontsize=SI_BASE, lineheight=0.9), vjust=0.5))
H <- list(hh("Individual trees"), hh("Species median gene,\nindividual flux"), hh("Species medians"), hh("Species means"))
P <- function(k, lv) res[[paste(k, lv)]]$panel
cells <- c(H, unlist(lapply(PRED$key, function(k) lapply(levs, function(lv) P(k, lv))), recursive=FALSE), bars)
combined <- wrap_plots(cells, ncol=4, byrow=TRUE, heights=c(0.16, 1, 1, 1, 1, 1, 1.15)) +
  plot_layout(axis_titles="collect_y")   # no species legend: points are grey
si_save(out_path("figS11_final.png"), combined, height=9.8)
for (lv in levs) cat(sprintf("%-4s gene R2: %s\n", lv, paste(sprintf("%s=%.3f(p%.3f)", D[[lv]]$Model, D[[lv]]$R2, D[[lv]]$P), collapse="  ")))
cat("Wrote figS11_final.png\n")
