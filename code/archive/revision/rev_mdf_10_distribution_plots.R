#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_10_distribution_plots.R      STEP 10
# ------------------------------------------------------------------------------
# Flux distributions and their relationship to chamber starting concentration.
#
# FIGURE 1  distributions
#   a  all stem fluxes, by campaign          (pseudo-log, so negatives are visible)
#   b  positive fluxes only                  (log)
#   c  negative fluxes only, |flux|          (log)
#   d  detection class composition by campaign
#
# FIGURE 2  relationship to chamber starting concentration
#   a  flux vs C0, all campaigns             (C0 is absolute; campaigns differ in ambient)
#   b  flux vs excess above LOCAL ambient    (2023 only, where traces were matched)
#   c  |negative flux| vs excess
#   d  prevalence of negatives by excess bin  -- the clearest view of the effect
#
# 'excess' = C0 minus the 10th percentile of the 3 min of trace either side of the
# enclosure, i.e. how far above local ambient the chamber opened. See step 9.
#
# Output: outputs/revision/fig_flux_distributions.png
#         outputs/revision/fig_flux_vs_c0.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

ALL <- read.csv("outputs/revision/mdf_allan_all_campaigns.csv", stringsAsFactors=FALSE)
G <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="height_2021"),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="semirigid_soil"),
  read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="static_2023")) %>%
  filter(is.finite(best.flux), is.finite(C0), is.finite(nb.obs), nb.obs>2)
SIG <- ALL %>% mutate(campaign=ifelse(campaign=="untagged","semirigid_tree",campaign)) %>%
  group_by(campaign) %>% summarise(sigma=median(allan_sd_CH4,na.rm=TRUE), .groups="drop")
SIG <- bind_rows(SIG, data.frame(campaign="static_2023", sigma=1.298/0.86))
G <- G %>% left_join(SIG,"campaign") %>% filter(is.finite(sigma))
G$MDF95 <- G$sigma*qnorm(.975)/G$nb.obs*G$flux.term
G$class <- ifelse(G$best.flux < -G$MDF95, "sig. negative",
           ifelse(G$best.flux >  G$MDF95, "sig. positive", "below MDF"))
G$class <- factor(G$class, levels=c("sig. negative","below MDF","sig. positive"))
G$type <- ifelse(G$campaign=="semirigid_soil", "soil", "stem")
CL <- c(`sig. negative`="#B2182B", `below MDF`="grey65", `sig. positive`="#2166AC")
psl <- scales::pseudo_log_trans(sigma=.005, base=10)

L <- read.csv("outputs/revision/local_ambient_check.csv", stringsAsFactors=FALSE)

# ============================== FIGURE 1 =====================================
S <- G[G$type=="stem",]
p1a <- ggplot(S, aes(best.flux, fill=class)) +
  geom_histogram(bins=70) + facet_wrap(~campaign, ncol=1, scales="free_y") +
  geom_vline(xintercept=0, linetype="22", linewidth=.3) +
  scale_x_continuous(trans=psl) + scale_fill_manual(values=CL, name=NULL) +
  labs(title="a  all stem fluxes", x=expression(CH[4]~(nmol~m^-2~s^-1)), y="count") +
  theme_bw(base_size=8) + theme(legend.position="bottom", legend.text=element_text(size=6))

p1b <- ggplot(S[S$best.flux>0,], aes(best.flux, colour=campaign)) +
  stat_ecdf(linewidth=.6) + scale_x_log10() +
  labs(title="b  positive fluxes (ECDF)", x=expression(CH[4]~(nmol~m^-2~s^-1)),
       y="cumulative fraction", colour=NULL) +
  theme_bw(base_size=8) + theme(legend.position="bottom", legend.text=element_text(size=6))

p1c <- ggplot(S[S$best.flux<0,], aes(abs(best.flux), colour=campaign)) +
  stat_ecdf(linewidth=.6) + scale_x_log10() +
  labs(title="c  negative fluxes, |flux| (ECDF)", x=expression("|"*CH[4]*"| (nmol m"^-2*" s"^-1*")"),
       y="cumulative fraction", colour=NULL) +
  theme_bw(base_size=8) + theme(legend.position="bottom", legend.text=element_text(size=6))

p1d <- G %>% count(campaign, class) %>% group_by(campaign) %>% mutate(f=n/sum(n)) %>%
  ggplot(aes(f, campaign, fill=class)) + geom_col() +
  scale_fill_manual(values=CL, name=NULL) + scale_x_continuous(labels=scales::percent) +
  labs(title="d  detection class composition", x=NULL, y=NULL) +
  theme_bw(base_size=8) + theme(legend.position="bottom", legend.text=element_text(size=6))

ggsave(file.path(outdir,"fig_flux_distributions.png"),
  (p1a | (p1b/p1c)) / p1d + plot_layout(heights=c(2.4,1)) +
  plot_annotation(title="Stem CH4 flux distributions",
    subtitle="detection class at Wassmann 95% with campaign-level empirical precision",
    theme=theme(plot.title=element_text(size=11,face="bold"))),
  width=11, height=9, dpi=200, bg="white")

# ============================== FIGURE 2 =====================================
p2a <- ggplot(G, aes(C0, best.flux, colour=class)) +
  geom_hline(yintercept=0, linetype="22", colour="grey50", linewidth=.3) +
  geom_point(size=.9, alpha=.7) + facet_wrap(~campaign, scales="free_x", nrow=1) +
  scale_y_continuous(trans=psl) + scale_colour_manual(values=CL, name=NULL) +
  labs(title="a  flux vs chamber starting concentration",
       x=expression(C[0]~(ppb)), y=expression(CH[4]~(nmol~m^-2~s^-1))) +
  theme_bw(base_size=8) + theme(legend.position="none")

LL <- L %>% mutate(class=factor(ifelse(flux < -MDF95, "sig. negative",
                    ifelse(flux > MDF95, "sig. positive","below MDF")),
                    levels=names(CL)))
p2b <- ggplot(LL, aes(excess, flux, colour=class)) +
  geom_hline(yintercept=0, linetype="22", colour="grey50", linewidth=.3) +
  geom_vline(xintercept=50, linetype="33", colour="grey40", linewidth=.3) +
  geom_point(size=1.2, alpha=.8) +
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=5, base=10)) +
  scale_y_continuous(trans=psl) + scale_colour_manual(values=CL, name=NULL) +
  labs(title="b  flux vs excess above LOCAL ambient (2023)",
       subtitle="dotted line = 50 ppb",
       x="C0 - local ambient (ppb)", y=expression(CH[4]~(nmol~m^-2~s^-1))) +
  theme_bw(base_size=8) + theme(legend.position="bottom", legend.text=element_text(size=6))

p2c <- ggplot(LL[LL$flux<0,], aes(excess, abs(flux))) +
  geom_point(aes(colour=class), size=1.4, alpha=.85) +
  geom_smooth(method="lm", se=TRUE, colour="grey25", linewidth=.5, formula=y~x) +
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=5, base=10)) + scale_y_log10() +
  scale_colour_manual(values=CL, guide="none") +
  labs(title="c  |negative flux| vs excess", x="C0 - local ambient (ppb)",
       y=expression("|"*CH[4]*"| (nmol m"^-2*" s"^-1*")")) +
  theme_bw(base_size=8)

LL$bin <- cut(LL$excess, c(-Inf,0,10,25,50,100,Inf),
              labels=c("<0","0-10","10-25","25-50","50-100",">100"))
bb <- LL %>% filter(!is.na(bin)) %>% group_by(bin) %>%
  summarise(n=n(), pct_neg=100*mean(flux<0), pct_signeg=100*mean(flux < -MDF95), .groups="drop")
p2d <- bb %>% pivot_longer(c(pct_neg,pct_signeg)) %>%
  mutate(name=ifelse(name=="pct_neg","any negative","significantly negative")) %>%
  ggplot(aes(bin, value, fill=name)) +
  geom_col(position="dodge") +
  geom_text(data=bb, aes(bin, -2, label=paste0("n=",n)), inherit.aes=FALSE, size=2.2, colour="grey35") +
  scale_fill_manual(values=c(`any negative`="#F4A582", `significantly negative`="#B2182B"), name=NULL) +
  labs(title="d  negative prevalence by excess bin", x="C0 - local ambient (ppb)", y="% of measurements") +
  theme_bw(base_size=8) + theme(legend.position="bottom", legend.text=element_text(size=6))

ggsave(file.path(outdir,"fig_flux_vs_c0.png"),
  p2a / (p2b | p2c | p2d) + plot_layout(heights=c(1,1.15)) +
  plot_annotation(title="Flux vs chamber starting concentration",
    theme=theme(plot.title=element_text(size=11,face="bold"))),
  width=12, height=8, dpi=200, bg="white")

cat("Written: fig_flux_distributions.png, fig_flux_vs_c0.png\n\n")
cat("Negative prevalence by excess bin:\n"); print(as.data.frame(bb), row.names=FALSE, digits=3)
cat("\nCorrelations (2023, matched to trace):\n")
n <- LL[LL$flux<0,]
cat(sprintf("  |negative flux| vs excess : Spearman rho %.3f (p=%.3g, n=%d)\n",
  suppressWarnings(cor.test(n$excess, abs(n$flux), method="spearman")$estimate),
  suppressWarnings(cor.test(n$excess, abs(n$flux), method="spearman")$p.value), nrow(n)))
cat(sprintf("  all flux vs excess        : Spearman rho %.3f (p=%.3g, n=%d)\n",
  suppressWarnings(cor.test(LL$excess, LL$flux, method="spearman")$estimate),
  suppressWarnings(cor.test(LL$excess, LL$flux, method="spearman")$p.value), nrow(LL)))
for (cp in unique(G$campaign)) { d <- G[G$campaign==cp,]
  cat(sprintf("  %-16s flux vs C0 : rho %.3f (p=%.3g)\n", cp,
    suppressWarnings(cor.test(d$C0,d$best.flux,method="spearman")$estimate),
    suppressWarnings(cor.test(d$C0,d$best.flux,method="spearman")$p.value))) }
