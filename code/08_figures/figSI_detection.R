source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# figSI_detection.R
# ------------------------------------------------------------------------------
# SI figure set for precision, detection, and the stem-uptake question.
#
#   a  measured flux against its own detection limit, by campaign
#   b  flux distributions with detection class
#   c  detected uptake prevalence against instrument precision, by campaign
#   d  effect of each screening criterion on retained n and apparent uptake
#
# Uses the final method: sigma = MAD(dx)/sqrt(2) per field period,
# MDF = z*sigma/t*flux.term at 90%.
# ==============================================================================
suppressMessages({library(dplyr);library(ggplot2);library(tidyr);library(patchwork)})
outdir <- "outputs"
G <- read.csv(out_path("flux_FINAL.csv"), stringsAsFactors=FALSE)
G$camp <- factor(G$camp, levels=c("Height+molecular","Cross-species","Monthly survey"))
G$class <- factor(G$class, levels=c("uptake","below detection","emission"))
CL <- c(uptake="#B2182B", `below detection`="grey70", emission="#2166AC")
S <- G[G$type=="stem",]

pa <- ggplot(G, aes(abs(best.flux), MDF, colour=class)) +
  geom_abline(slope=1, intercept=0, linetype="22", colour="grey40", linewidth=.4) +
  geom_point(size=.8, alpha=.7) +
  facet_wrap(~paste0(camp, " (", type, ")"), nrow=1) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values=CL, name=NULL) +
  labs(title="a  measured flux against its detection limit",
       subtitle="points below the 1:1 line are resolved; above it, indistinguishable from zero",
       x=expression("|"*CH[4]*" flux| (nmol m"^-2*" s"^-1*")"), y="MDF") +
  theme_bw(base_size=8) + theme(legend.position="bottom", strip.text=element_text(size=6.5))

pb <- ggplot(S, aes(best.flux, fill=class)) +
  geom_histogram(bins=70) + facet_wrap(~camp, ncol=1, scales="free_y") +
  geom_vline(xintercept=0, linetype="22", linewidth=.3) +
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=.005, base=10)) +
  scale_fill_manual(values=CL, name=NULL) +
  labs(title="b  stem flux distribution by campaign",
       x=expression(CH[4]~(nmol~m^-2~s^-1)), y="count") +
  theme_bw(base_size=8) + theme(legend.position="none")

cs <- S %>% group_by(camp) %>% summarise(sigma=sigma[1], n=n(),
  pct=100*mean(class=="uptake"), .groups="drop")
pc <- ggplot(cs, aes(sigma, pct)) +
  geom_line(colour="grey60", linewidth=.4) +
  geom_point(aes(size=n), colour="#B2182B") +
  geom_text(aes(label=camp), vjust=-1.4, size=2.4) +
  scale_size_continuous(range=c(2.5,5), guide="none") +
  expand_limits(y=c(-1,13), x=c(1.0,2.5)) +
  labs(title="c  apparent stem uptake tracks instrument precision",
       subtitle="the quietest period detects none",
       x=expression(sigma~(ppb~CH[4])), y="% of stem measurements detected as uptake") +
  theme_bw(base_size=8)

CR <- list(`retain all`=rep(TRUE,nrow(S)), `p < 0.05`=S$LM.p.val<0.05,
           `above MDF (90%)`=abs(S$best.flux)>S$MDF,
           `r2 > 0.7`=S$LM.r2>0.7, `r2 > 0.9`=S$LM.r2>0.9)
D <- bind_rows(lapply(names(CR), function(nm){ m<-CR[[nm]]; m[is.na(m)]<-FALSE; f<-S$best.flux[m]
  data.frame(criterion=nm, kept=sum(m), pct_kept=100*mean(m), pct_uptake=100*mean(f<0),
             bias=100*(mean(f)/mean(S$best.flux)-1))}))
D$criterion <- factor(D$criterion, levels=rev(names(CR)))
pd <- D %>% pivot_longer(c(pct_kept,pct_uptake,bias)) %>%
  mutate(name=recode(name, pct_kept="% retained", pct_uptake="% of retained that are uptake",
                     bias="change in mean flux (%)")) %>%
  ggplot(aes(value, criterion, fill=name)) + geom_col(show.legend=FALSE) +
  facet_wrap(~name, scales="free_x") +
  scale_fill_manual(values=c("#4393C3","#B2182B","#F4A582")) +
  labs(title="d  what each screening criterion does",
       subtitle="an r2 threshold removes nearly all apparent uptake and inflates the mean",
       x=NULL, y=NULL) + theme_bw(base_size=7.5)

ggsave(out_path("fig_SI_detection.png"), pa/(pb|pc)/pd +
  plot_layout(heights=c(1,1.25,0.9)) +
  plot_annotation(title="Detection limits and the interpretation of negative stem fluxes",
    theme=theme(plot.title=element_text(size=11,face="bold"))),
  width=10.5, height=11, dpi=200, bg="white")
cat("Written: outputs/figures/generated/fig_SI_detection.png\n\n")
print(as.data.frame(D), row.names=FALSE, digits=3)
