#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS19_mcra-probe-validation.R
# SI figure: validation of the probe-based mcrA ddPCR assay against the original
# intercalating-dye (EvaGreen/ML) assay on the same wood samples. Shows (A) the
# two assays agree quantitatively where both detect, and (B) EvaGreen calls extra
# positives the specific probe rejects (its false-positive signature — the reason
# the probe assay was developed; Arnold et al. 2024). Supports the R2 controls/LoD
# response. Reads data/; writes outputs/revision/figS19_mcra_probe_validation.png
# ==============================================================================
suppressMessages({library(ggplot2); library(patchwork)}); options(warn=-1, stringsAsFactors=FALSE)
dd <- read.csv("data/compiled/ddpcr_gene_abundances.csv")
dd$conc <- suppressWarnings(as.numeric(dd$concentration_copies_per_uL))
dd$copies_g <- dd$conc * 75 / dd$sample_mass_mg * 1000        # per g dry (harmonization)

AT <- "loose"                                                 # standard detection basis
a <- dd[dd$target_gene=="mcra"       & dd$analysis_type==AT, c("sample_id","copies_g")]
b <- dd[dd$target_gene=="mcra_probe" & dd$analysis_type==AT & dd$material=="Wood", c("sample_id","copies_g")]
names(a)<-c("sample_id","eg"); names(b)<-c("sample_id","pr")
m <- merge(a,b,by="sample_id"); m <- m[is.finite(m$eg)&is.finite(m$pr),]
m$cat <- with(m, ifelse(eg>0 & pr>0, "both detect",
                 ifelse(eg>0 & pr==0, "EvaGreen only (probe rejects)",
                 ifelse(eg==0 & pr>0, "probe only", "both zero"))))
pos <- m[m$eg>0 & m$pr>0, ]
fit <- lm(log10(pr)~log10(eg), data=pos)
r  <- cor(log10(pos$eg), log10(pos$pr))
lab <- sprintf("both-positive n=%d\nPearson r=%.2f\nslope=%.2f", nrow(pos), r, coef(fit)[2])

cols <- c("both detect"="#2c7fb8","EvaGreen only (probe rejects)"="#d95f02",
          "probe only"="#7570b3","both zero"="grey70")
rng <- range(c(pos$eg,pos$pr)); rng[1]<-max(rng[1], min(pos$eg[pos$eg>0]))

pA <- ggplot(pos, aes(eg, pr)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey50") +
  geom_point(colour=cols[["both detect"]], alpha=.6, size=2) +
  geom_smooth(method="lm", se=FALSE, colour="black", linewidth=.6) +
  scale_x_log10() + scale_y_log10() +
  annotate("text", x=rng[1], y=max(pos$pr), label=lab, hjust=0, vjust=1, size=3.4) +
  labs(x=expression(EvaGreen~italic(mcrA)~(copies~g^{-1})),
       y=expression(Probe~italic(mcrA)~(copies~g^{-1})),
       title="A  Quantitative agreement (both assays detect)") +
  theme_bw(base_size=11)

tab <- as.data.frame(table(m$cat)); names(tab)<-c("cat","n"); tab$cat<-factor(tab$cat, levels=names(cols))
pB <- ggplot(tab, aes(reorder(cat,-n), n, fill=cat)) +
  geom_col(width=.7) + geom_text(aes(label=n), vjust=-0.3, size=3.4) +
  scale_fill_manual(values=cols, guide="none") +
  labs(x=NULL, y="wood samples",
       title=sprintf("B  Detection concordance (n=%d paired wood samples)", nrow(m))) +
  theme_bw(base_size=11) + theme(axis.text.x=element_text(angle=20, hjust=1))

p <- pA + pB + plot_layout(widths=c(1,1.1))
ggsave("outputs/revision/figS19_mcra_probe_validation.png", p, width=11, height=4.6, dpi=300)
cat("detection: EvaGreen", sum(m$eg>0),"/",nrow(m)," probe", sum(m$pr>0),"/",nrow(m),
    " | EvaGreen-only(false+):", sum(m$eg>0&m$pr==0), " probe-only:", sum(m$eg==0&m$pr>0),"\n")
cat("wrote outputs/revision/figS19_mcra_probe_validation.png\n")
