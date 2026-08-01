#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS19_mcra-probe-validation.R  (SI)
# Validate the probe-based mcrA ddPCR assay against the original EvaGreen/ML assay
# on the same wood samples. Pseudo-log axes so non-detects sit at true 0 and low
# detections rise continuously (no artificial gap). EvaGreen calls many positives
# the specific probe scores zero (its false-positive signature; Arnold et al. 2024),
# while agreeing quantitatively where both detect. Compared on copies/uL (direct
# measurement). Writes BOTH axis orientations for choosing; the canonical
# figS19_mcra_probe_validation.png = probe-on-x (trusted-assay-as-reference).
# ==============================================================================
suppressMessages({library(ggplot2); library(scales)}); options(warn=-1, stringsAsFactors=FALSE)
dd <- read.csv("data/compiled/ddpcr_gene_abundances.csv")
dd$conc <- suppressWarnings(as.numeric(dd$concentration_copies_per_uL))   # copies/uL (direct)
AT <- "loose"
a <- dd[dd$target_gene=="mcra"       & dd$analysis_type==AT, c("sample_id","conc")]
b <- dd[dd$target_gene=="mcra_probe" & dd$analysis_type==AT & dd$material=="Wood", c("sample_id","conc")]
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
  ggsave(out, p, width=7.4, height=6.3, dpi=300); out
}
build("prp","egp","pr","eg", prlab, eglab, "outputs/revision/figS19_mcra_probe-x.png")      # probe on X
build("egp","prp","eg","pr", eglab, prlab, "outputs/revision/figS19_mcra_evagreen-x.png")   # EvaGreen on X (chosen)
# canonical = EvaGreen-on-x: false positives form a band pinned to 0 on OUR probe measure
file.copy("outputs/revision/figS19_mcra_evagreen-x.png", "outputs/revision/figS19_mcra_probe_validation.png", overwrite=TRUE)
cat("wrote figS19_mcra_{probe-x,evagreen-x}.png (canonical = EvaGreen-x)\n")
