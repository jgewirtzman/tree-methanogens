#!/usr/bin/env Rscript
# ==============================================================================
# figS04_mcra-probe-validation.R -- Figure S4: probe mcrA assay against EvaGreen
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
# The probe-based mcrA ddPCR assay (this study) against the EvaGreen assay of the
# literature on the same 2021 wood extracts, copies per uL of reaction (loose
# threshold, runs that passed QC; a sample run more than once is its mean). Pseudo-log
# axes, so non-detects sit at 0 (jittered within the grey band) and low detections
# rise continuously. EvaGreen on x: samples EvaGreen calls positive and the probe scores
# zero (its false-positive signature; Arnold et al. 2024) form a band at probe = 0.
# Line: log-log fit over samples both assays detect; dashed: 1:1.
# Inputs: data/tables/ddpcr.csv, samples.csv (code/lib/tables_read.R)
# Output: outputs/figures/generated/figS19_mcra_probe_validation.png
# ==============================================================================
suppressMessages({library(ggplot2); library(scales)})
source("code/lib/tables_read.R")
dd <- read_table("ddpcr"); smp <- read_table("samples")
wood <- smp$sample_id[smp$campaign == "2021_survey" & smp$material == "wood"]
dd <- dd[dd$qc_pass & dd$gene == "mcrA" & dd$sample_id %in% wood & !is.na(dd$copies_per_uL_rxn), ]
per_sample <- function(assay) { x <- dd[dd$assay == assay, ]
  data.frame(sample_id = sort(unique(x$sample_id)), conc = as.vector(tapply(x$copies_per_uL_rxn, x$sample_id, mean)[sort(unique(x$sample_id))])) }
a <- per_sample("EvaGreen"); b <- per_sample("probe")
names(a)<-c("sample_id","eg"); names(b)<-c("sample_id","pr")
m <- merge(a,b,by="sample_id"); m <- m[is.finite(m$eg)&is.finite(m$pr),]
m$cat <- with(m, ifelse(eg>0 & pr>0, "both detect",
                 ifelse(eg>0 & pr==0, "EvaGreen only (probe = 0)",
                 ifelse(eg==0 & pr>0, "probe only", "both = 0"))))
sig <- 0.05; pslog <- pseudo_log_trans(sigma=sig, base=10)
set.seed(1); zj <- function(n) runif(n, -sig*0.5, sig*0.5)
m$egp <- ifelse(m$eg==0, zj(nrow(m)), m$eg); m$prp <- ifelse(m$pr==0, zj(nrow(m)), m$pr)
topmax <- max(c(m$eg, m$pr))
pos <- m[m$eg>0 & m$pr>0, ]
lab <- sprintf("both detect: n=%d, r=%.2f, slope=%.2f\nEvaGreen only (false +): %d\nprobe only: %d",
               nrow(pos), cor(log10(pos$eg),log10(pos$pr)), coef(lm(log10(pr)~log10(eg),data=pos))[2],
               sum(m$cat=="EvaGreen only (probe = 0)"), sum(m$cat=="probe only"))
cols <- c("both detect"="#2c7fb8","EvaGreen only (probe = 0)"="#d95f02",
          "probe only"="#7570b3","both = 0"="grey55")
brk <- c(0,0.1,1,10,100,1e3,1e4,1e5,1e6)
lbl <- c("0",expression(10^-1),expression(10^0),expression(10^1),expression(10^2),
         expression(10^3),expression(10^4),expression(10^5),expression(10^6))
eglab <- expression(EvaGreen~italic(mcrA)~assay~"(literature)"~~copies~mu*L^{-1})
prlab <- expression(Probe~italic(mcrA)~assay~"(this study)"~~copies~mu*L^{-1})

# x/y: plotted columns (non-detects jittered); xr/yr: the measured columns for the fit
build <- function(xp, yp, xr, yr, xlab, ylab, out){
  px<-pos[[xr]]; py<-pos[[yr]]; cf<-coef(lm(log10(py)~log10(px)))
  fx<-10^seq(log10(min(px)),log10(max(px)),length=100); fit<-data.frame(x=fx, y=10^(cf[1]+cf[2]*log10(fx)))
  p <- ggplot(m, aes(.data[[xp]], .data[[yp]], colour=cat)) +
    annotate("rect", xmin=-sig, xmax=sig, ymin=-Inf, ymax=Inf, fill="grey92", alpha=.6) +
    annotate("rect", ymin=-sig, ymax=sig, xmin=-Inf, xmax=Inf, fill="grey92", alpha=.6) +
    geom_abline(slope=1, intercept=0, linetype=2, colour="grey45") +
    geom_line(data=fit, aes(x,y), colour="black", linewidth=.6, inherit.aes=FALSE) +
    geom_point(alpha=.7, size=1.9) +
    scale_x_continuous(trans=pslog, breaks=brk, labels=lbl, limits=c(-sig,topmax)) +
    scale_y_continuous(trans=pslog, breaks=brk, labels=lbl, limits=c(-sig,topmax)) +
    scale_colour_manual(values=cols, name=NULL) +
    annotate("text", x=0, y=topmax, label=lab, hjust=0, vjust=1, size=3.2) +
    labs(x=xlab, y=ylab) +
    theme_bw(base_size=11) + theme(legend.position=c(.99,.02), legend.justification=c(1,0),
                                   legend.background=element_rect(fill=alpha("white",.7), colour=NA))
  ggsave(out, p, width=7.4, height=6.3, dpi=300)
}
build("egp","prp","eg","pr", eglab, prlab, "outputs/figures/generated/figS19_mcra_probe_validation.png")
cat("wrote outputs/figures/generated/figS19_mcra_probe_validation.png\n")
