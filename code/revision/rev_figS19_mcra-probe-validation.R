#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS19_mcra-probe-validation.R
# SI figure: validation of the probe-based mcrA ddPCR assay against the original
# intercalating-dye (EvaGreen/ML) assay on the same wood samples. Shows both
# non-detects and detections on a 1:1 plot so the false-positive pattern is
# explicit: EvaGreen calls many samples positive that the specific probe scores
# zero (its false-positive signature; Arnold et al. 2024), while agreeing
# quantitatively where both detect. Supports the R2 controls/LoD response.
# Units: copies g-1 dry wood. Reads data/; writes outputs/revision/.
# ==============================================================================
suppressMessages({library(ggplot2); library(scales)}); options(warn=-1, stringsAsFactors=FALSE)
dd <- read.csv("data/compiled/ddpcr_gene_abundances.csv")
# Compare on copies/uL (the direct ddPCR measurement). copies/g would multiply BOTH
# assays of a sample by the same 75/mass factor, artificially inflating agreement;
# copies/uL removes that shared mass term (fairer) and tightens the detection gap.
dd$copies_g <- suppressWarnings(as.numeric(dd$concentration_copies_per_uL))

AT <- "loose"
a <- dd[dd$target_gene=="mcra"       & dd$analysis_type==AT, c("sample_id","copies_g")]
b <- dd[dd$target_gene=="mcra_probe" & dd$analysis_type==AT & dd$material=="Wood", c("sample_id","copies_g")]
names(a)<-c("sample_id","eg"); names(b)<-c("sample_id","pr")
m <- merge(a,b,by="sample_id"); m <- m[is.finite(m$eg)&is.finite(m$pr),]
m$cat <- with(m, ifelse(eg>0 & pr>0, "both detect",
                 ifelse(eg>0 & pr==0, "EvaGreen only (probe = 0)",
                 ifelse(eg==0 & pr>0, "probe only", "both = 0"))))
# axes: probe (trusted/reference) on X, EvaGreen (assay under test) on Y, so its
# false positives read as a vertical stripe rising from probe = 0 (non-detect).

# ---- place non-detects (0) at a floor a decade below the lowest detected value ----
posmin <- min(c(m$eg[m$eg>0], m$pr[m$pr>0]))
floor  <- 10^(floor(log10(posmin)) - 1)                       # one clean decade below
topmax <- max(c(m$eg, m$pr))
set.seed(1)
jit <- function(x) x * 10^runif(length(x), -0.18, 0.18)       # spread the piled-up zeros
m$egp <- ifelse(m$eg==0, jit(rep(floor, nrow(m))), m$eg)
m$prp <- ifelse(m$pr==0, jit(rep(floor, nrow(m))), m$pr)

pos <- m[m$eg>0 & m$pr>0, ]; fit <- lm(log10(pr)~log10(eg), data=pos)
lab <- sprintf("both detect: n=%d, r=%.2f, slope=%.2f\nEvaGreen only (false +): %d\nprobe only: %d",
               nrow(pos), cor(log10(pos$eg),log10(pos$pr)), coef(fit)[2],
               sum(m$cat=="EvaGreen only (probe = 0)"), sum(m$cat=="probe only"))

cols <- c("both detect"="#2c7fb8","EvaGreen only (probe = 0)"="#d95f02",
          "probe only"="#7570b3","both = 0"="grey55")
brk <- 10^(seq(floor(log10(floor)), ceiling(log10(topmax))))
lbl <- ifelse(brk==floor, "0 (nd)", parse(text=paste0("10^", round(log10(brk)))))

p <- ggplot(m, aes(prp, egp, colour=cat)) +
  annotate("rect", xmin=floor/3.5, xmax=floor*3.5, ymin=-Inf, ymax=Inf, fill="grey92", alpha=.5) +
  annotate("rect", ymin=floor/3.5, ymax=floor*3.5, xmin=-Inf, xmax=Inf, fill="grey92", alpha=.5) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey45") +
  geom_smooth(data=pos, aes(pr,eg), method="lm", se=FALSE, colour="black", linewidth=.6, inherit.aes=FALSE) +
  geom_point(alpha=.7, size=1.9) +
  scale_x_log10(breaks=brk, labels=lbl) + scale_y_log10(breaks=brk, labels=lbl) +
  annotation_logticks(sides="bl", colour="grey60") +
  scale_colour_manual(values=cols, name=NULL) +
  annotate("text", x=floor, y=topmax, label=lab, hjust=0, vjust=1, size=3.2) +
  labs(x=expression(Probe~italic(mcrA)~(copies~mu*L^{-1}~reaction)~-~trusted~assay),
       y=expression(EvaGreen~italic(mcrA)~(copies~mu*L^{-1}~reaction)),
       title=expression("EvaGreen vs. probe "*italic(mcrA)*" ddPCR (paired wood samples; 0 = non-detect)")) +
  theme_bw(base_size=11) + theme(legend.position=c(.99,.02), legend.justification=c(1,0),
                                 legend.background=element_rect(fill=alpha("white",.7), colour=NA))
ggsave("outputs/revision/figS19_mcra_probe_validation.png", p, width=7.6, height=6.4, dpi=300)
cat("n paired:",nrow(m)," EvaGreen-only(false+):",sum(m$cat=="EvaGreen only (probe = 0)"),
    " probe-only:",sum(m$cat=="probe only")," both0:",sum(m$cat=="both = 0"),
    " | floor=",floor,"\n")
cat("wrote outputs/revision/figS19_mcra_probe_validation.png\n")
