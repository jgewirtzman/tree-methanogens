#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_scaling_heatmap.R
# ------------------------------------------------------------------------------
# Summary view of the 280-combination scaling grid.
#
#   a  HEATMAP  flux form x WAI, cell = median across the 8 area shapes.
#      These are the two axes with real leverage; the area shapes are collapsed into
#      the cell and their spread is shown in panel b.
#   b  spread contributed by the 8 area shapes within each cell
#   c  measured vs extrapolated decomposition per flux form
#
# Diverging colour scale centred on zero: with the empirically bounded negative forms
# the tree term can be a net sink, so sign matters and a sequential scale would hide it.
# ==============================================================================
suppressMessages({library(dplyr);library(ggplot2);library(tidyr);library(patchwork)})
outdir <- "outputs/revision"
R <- read.csv(file.path(outdir,"scaling_full_grid.csv"), stringsAsFactors=FALSE)
FO <- c("constant","rf_learned_capped","power","exponential","linear_floored",
        "linear_bounded_median","linear_bounded_p05")
R$flux <- factor(R$flux, levels=rev(FO))
# linear_bounded_median is CONTRADICTED by the only direct evidence above 2 m: all four
# climbed-tree measurements above 2 m are positive (0.031-0.354 nmol m-2 s-1), each
# larger in magnitude than the assumed sink of -0.0247, and the profile's only negative
# is at 0.5 m where goFlux flags it below detection. Retained because n = 1 tree cannot
# exclude uptake in other species or conditions, but marked so it is not read as an
# equal member of the set.
CONTRADICTED <- "linear_bounded_median"
R$WAI  <- factor(R$WAI, levels=sort(unique(R$WAI)))

H <- R %>% group_by(flux, WAI) %>%
  summarise(med=median(total_mg), lo=min(total_mg), hi=max(total_mg),
            spread=max(total_mg)-min(total_mg), .groups="drop")
lim <- max(abs(H$med))
pa <- ggplot(H, aes(WAI, flux, fill=med)) +
  geom_tile(colour="white", linewidth=.6) +
  geom_text(aes(label=sprintf("%.0f", med)), size=2.9,
            colour=ifelse(abs(H$med)>0.55*lim, "white", "grey15")) +
  scale_fill_gradient2(low="#2166AC", mid="grey95", high="#B2182B", midpoint=0,
                       limits=c(-lim,lim), name=expression(mg~CH[4]~m^-2~yr^-1)) +
  labs(title="a  whole-woody-surface tree budget",
       subtitle="median across the 8 area-distribution shapes; blue = net sink\n* contradicted by the climbed tree: all 4 measurements above 2 m were positive",
       x="woody area index", y=NULL) +
  scale_y_discrete(labels=function(l) ifelse(l==CONTRADICTED, paste0(l, " *"), l)) +
  theme_minimal(base_size=8) +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="bottom",
        legend.key.height=unit(.25,"cm"), panel.grid=element_blank())

pb <- ggplot(H, aes(WAI, flux, fill=spread)) +
  geom_tile(colour="white", linewidth=.6) +
  geom_text(aes(label=sprintf("%.0f", spread)), size=2.7, colour="grey15") +
  scale_fill_gradient(low="grey96", high="#F4A582", name="range (mg)") +
  labs(title="b  spread contributed by the area shapes",
       subtitle="zero for constant flux: the vertical distribution cancels algebraically",
       x="woody area index", y=NULL) +
  theme_minimal(base_size=8) +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="bottom",
        legend.key.height=unit(.25,"cm"), panel.grid=element_blank())

D <- R %>% group_by(flux) %>%
  summarise(measured=median(measured_mg), extrapolated=median(extrapolated_mg), .groups="drop") %>%
  pivot_longer(c(measured,extrapolated))
pc <- ggplot(D, aes(value, flux, fill=name)) +
  geom_col() + geom_vline(xintercept=0, colour="grey40", linewidth=.3) +
  scale_fill_manual(values=c(measured="#2166AC", extrapolated="#F4A582"), name=NULL) +
  labs(title="c  measured (0-2 m) vs extrapolated (>2 m)",
       subtitle="the measured band is 5.8 mg and identical in every combination",
       x=expression(mg~CH[4]~m^-2~yr^-1), y=NULL) +
  theme_bw(base_size=8) + theme(legend.position="bottom")

ggsave(file.path(outdir,"fig_scaling_heatmap.png"), (pa|pb)/pc + plot_layout(heights=c(1.5,1)) +
  plot_annotation(title="Stem CH4 scaling: 280 combinations",
    subtitle="5 WAI x 2 bole shapes x 4 branch shapes x 7 vertical flux forms",
    theme=theme(plot.title=element_text(size=11,face="bold"))),
  width=11, height=8.5, dpi=200, bg="white")
cat("Written: outputs/revision/fig_scaling_heatmap.png\n\n")
cat("Cell medians (mg CH4 m-2 yr-1):\n\n")
print(H %>% dplyr::select(flux,WAI,med) %>% pivot_wider(names_from=WAI, values_from=med) %>%
  as.data.frame(), row.names=FALSE, digits=3)
cat("\nRange by flux form (across all WAI and shapes):\n\n")
print(R %>% group_by(flux) %>% summarise(lo=min(total_mg), med=median(total_mg),
  hi=max(total_mg), pct_soil_lo=min(pct_of_soil), pct_soil_hi=max(pct_of_soil),
  .groups="drop") %>% arrange(desc(med)) %>% as.data.frame(), row.names=FALSE, digits=3)
