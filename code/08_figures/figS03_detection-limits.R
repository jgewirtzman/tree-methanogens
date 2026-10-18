#!/usr/bin/env Rscript
# ==============================================================================
# figS03_detection-limits.R -- Figure S3: detection limits and negative stem fluxes
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) |flux| against its own minimum detectable flux (MDF), by campaign and type
#   (b) stem flux distributions by detection class
#   (c) stem measurements detected as uptake against instrument precision, by campaign
#   (d) effect of each screening criterion on retained n, uptake share and mean flux
# Precision sigma = MAD(dx)/sqrt(2) per field period; MDF = z*sigma/t*flux.term at 90%
# (code/02_flux/03_precision_and_detection.R, which sets flux, MDF, sigma and class).
# Inputs: data/tables/flux_stem.csv and flux_soil.csv, every surveyed deployment with a
#   detection limit (the felled oak has none; nor have the two soil deployments that
#   failed the closure screen)
# Output: outputs/figures/generated/fig_SI_detection.png; the criterion table to the console
# ==============================================================================
suppressMessages({library(dplyr);library(ggplot2);library(tidyr);library(patchwork)})
source("code/lib/outputs.R")
source("code/lib/tables_read.R")
CAMP <- c(`2021_multiheight`="Height+molecular", `2023_cross_species`="Cross-species", monthly_2020_2021="Monthly survey")
.g <- function(d, type) data.frame(camp=unname(CAMP[d$campaign]), type=type, best.flux=d$CH4_flux_nmol_m2_s,
  MDF=d$CH4_MDF_nmol_m2_s, sigma=d$analyzer_sigma_ppb, class=d$detection_class, LM.p.val=d$CH4_LM_pval,
  LM.r2=d$CH4_LM_r2, stringsAsFactors=FALSE)
G <- rbind(.g(subset(read_table("flux_stem"), campaign %in% names(CAMP)), "stem"),
           .g(subset(read_table("flux_soil"), qc_pass), "soil"))
G$camp <- factor(G$camp, levels=c("Height+molecular","Cross-species","Monthly survey"))
G$class <- factor(G$class, levels=c("uptake","below detection","emission"))
# sink blue / source red, the convention of Figure 9 and the rest of the paper
CL <- c(uptake="#2166AC", `below detection`="grey70", emission="#B2182B")
S <- G[G$type=="stem",]

pa <- ggplot(G, aes(abs(best.flux), MDF, colour=class)) +
  geom_abline(slope=1, intercept=0, linetype="22", colour="grey40", linewidth=.4) +
  geom_point(size=.8, alpha=.7) +
  facet_wrap(~paste0(camp, " (", type, ")"), nrow=1) +
  scale_x_log10(labels=scales::label_log()) + scale_y_log10(labels=scales::label_log()) + scale_colour_manual(values=CL, name=NULL) +
  labs(x=expression("|"*CH[4]*" flux| (nmol m"^-2*" s"^-1*")"),
       y=expression("MDF (nmol m"^-2*" s"^-1*")")) +
  theme_bw(base_size=11) + theme(legend.position="bottom", strip.text=element_text(size=9.5))

pb <- ggplot(S, aes(best.flux, fill=class)) +
  geom_histogram(bins=70) + facet_wrap(~camp, ncol=1, scales="free_y") +
  geom_vline(xintercept=0, linetype="22", linewidth=.3) +
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=.005, base=10)) +
  scale_fill_manual(values=CL, name=NULL) +
  labs(x=expression(CH[4]~flux~(nmol~m^-2~s^-1)), y="Count") +
  theme_bw(base_size=11) + theme(legend.position="none", strip.text=element_text(size=9.5))

cs <- S %>% group_by(camp) %>% summarise(sigma=sigma[1], n=n(),
  pct=100*mean(class=="uptake"), .groups="drop")
pc <- ggplot(cs, aes(sigma, pct)) +
  geom_line(colour="grey60", linewidth=.4) +
  geom_point(aes(size=n), colour="#2166AC") +
  geom_text(aes(label=camp), vjust=-1.4, size=3.4) +
  scale_size_continuous(range=c(2.5,5), guide="none") +
  expand_limits(y=c(-1,13), x=c(1.0,2.5)) +
  labs(x=expression("Instrument precision, "*sigma~(ppb~CH[4])), y="Stem measurements detected\nas uptake (%)") +
  theme_bw(base_size=11)

CR <- list(`Retain all`=rep(TRUE,nrow(S)), `Slope p < 0.05`=S$LM.p.val<0.05,
           `Above MDF (90%)`=abs(S$best.flux)>S$MDF,
           `r² > 0.7`=S$LM.r2>0.7, `r² > 0.9`=S$LM.r2>0.9)
D <- bind_rows(lapply(names(CR), function(nm){ m<-CR[[nm]]; m[is.na(m)]<-FALSE; f<-S$best.flux[m]
  data.frame(criterion=nm, kept=sum(m), pct_kept=100*mean(m), pct_uptake=100*mean(f<0),
             bias=100*(mean(f)/mean(S$best.flux)-1))}))
D$criterion <- factor(D$criterion, levels=rev(names(CR)))
pd <- D %>% pivot_longer(c(pct_kept,pct_uptake,bias)) %>%
  mutate(name=factor(recode(name, pct_kept="Retained (%)", pct_uptake="Retained that are uptake (%)",
                     bias="Change in mean flux (%)"), levels=c("Retained (%)","Retained that are uptake (%)","Change in mean flux (%)"))) %>%
  ggplot(aes(value, criterion, fill=name)) + geom_col(show.legend=FALSE) +
  facet_wrap(~name, scales="free_x", strip.position="bottom") +
  scale_fill_manual(values=c("grey45","grey65","grey80")) +
  labs(x=NULL, y=NULL) + theme_bw(base_size=11) +
  theme(strip.placement="outside", strip.background=element_blank(), strip.text=element_text(size=10))

ggsave(out_path("fig_SI_detection.png"), pa/(pb|pc)/pd +
  plot_layout(heights=c(1,1.25,0.9)) +
  plot_annotation(tag_levels="a", tag_prefix="(", tag_suffix=")") &
  theme(plot.tag=element_text(size=14, face="bold")),
  width=12, height=13, dpi=300, bg="white")
cat("Written: outputs/figures/generated/fig_SI_detection.png\n\n")
print(as.data.frame(D), row.names=FALSE, digits=3)
