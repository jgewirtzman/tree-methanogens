#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_scaling_heatmap.R
# ------------------------------------------------------------------------------
# Summary view of the 240-combination scaling grid (5 WAI x 2 bole x 4 branch x 6 flux).
#
#   a  HEATMAP  flux form x WAI, cell = median across the 8 area shapes
#      (2 bole x 4 branch). These are the two axes with real leverage: across the
#      grid the flux form moves the median total by 91.5 mg CH4 m-2 yr-1 and WAI by
#      13.8, against 7.8 for branch placement and 0.6 for bole shape.
#   b  spread contributed by the 8 area shapes within each cell
#   c  measured vs extrapolated decomposition per flux form
#
# Diverging colour scale centred on zero: with the empirically bounded negative forms
# the tree term can be a net sink, so sign matters and a sequential scale would hide it.
# ==============================================================================
suppressMessages({library(dplyr);library(ggplot2);library(tidyr);library(patchwork)})
outdir <- "outputs/revision"
R <- read.csv(file.path(outdir,"scaling_full_grid.csv"), stringsAsFactors=FALSE)
# Forms carried in the current grid; linear_bounded_p05 was retired.
FO <- c("constant","rf_learned_capped","power","exponential","linear_floored",
        "linear_bounded_median")
R$flux <- factor(R$flux, levels=rev(FO))
# linear_bounded_median is CONTRADICTED by the only direct evidence above 2 m: all four
# climbed-tree measurements above 2 m are positive (0.031-0.354 nmol m-2 s-1), each
# larger in magnitude than the assumed sink of -0.0262, and the profile's only negative
# is at 0.5 m where goFlux flags it below detection. Retained because n = 1 tree cannot
# exclude uptake in other species or conditions, but marked so it is not read as an
# equal member of the set.
CONTRADICTED <- "linear_bounded_median"
R$WAI  <- factor(R$WAI, levels=sort(unique(R$WAI)))
# Every count and value in the titles is READ FROM THE GRID. These were hardcoded
# ("280 combinations", "7 vertical flux forms", "the measured band is 5.8 mg") and
# had drifted from a grid that is now 240 combinations, 6 forms and a 4.91 mg band
# -- the same failure as the "216" that was hardcoded in the Figure 9 caption.
N_COMB   <- nrow(R)
N_WAI    <- length(unique(R$WAI))
N_BOLE   <- length(unique(R$bole))
N_BRANCH <- length(unique(R$branch))
N_FLUX   <- length(unique(R$flux))
N_SHAPE  <- N_BOLE * N_BRANCH
MEAS     <- unique(round(R$measured_mg, 2))[1]

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
       subtitle=sprintf("median across the %d area-distribution shapes; blue = net sink\n* contradicted by the climbed tree: all 4 measurements above 2 m were positive", N_SHAPE),
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
       subtitle=sprintf("the measured band is %.2f mg and identical in every combination", MEAS),
       x=expression(mg~CH[4]~m^-2~yr^-1), y=NULL) +
  theme_bw(base_size=8) + theme(legend.position="bottom")

ggsave(file.path(outdir,"fig_scaling_heatmap.png"), (pa|pb)/pc + plot_layout(heights=c(1.5,1)) +
  plot_annotation(title=sprintf("Stem CH4 scaling: %d combinations", N_COMB),
    subtitle=sprintf("%d WAI x %d bole shapes x %d branch shapes x %d vertical flux forms",
                     N_WAI, N_BOLE, N_BRANCH, N_FLUX),
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
