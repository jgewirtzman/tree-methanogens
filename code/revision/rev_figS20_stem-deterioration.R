#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS20_stem-deterioration.R  (exploratory / SI candidate)
# Stem deterioration vs CH4 flux (2023 cross-species survey, species-controlled).
# Tests the early-decay hypothesis: emission should peak at *moderate* stem
# deterioration (moisture, labile C, anoxia) then fall as structure is lost and
# O2 enters. (A) arcsinh flux by 3-level bark condition (the hump); (B) composite
# deterioration score (bark+wounding) with a species-controlled quadratic fit.
# Reads data/; writes outputs/revision/figS20_stem_deterioration.png
# ==============================================================================
suppressMessages({library(ggplot2); library(patchwork); library(lme4); library(lmerTest)})
options(warn=-1, stringsAsFactors=FALSE)
asinh10<-function(x) asinh(x/0.1)/log(10)
y<-read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv", check.names=FALSE)
names(y)<-make.names(names(y)); y<-y[!is.na(y$CH4_best.flux),]
y$fx<-asinh10(y$CH4_best.flux); y$sp<-as.factor(y$Species)
ord<-function(v){v<-tolower(trimws(as.character(v))); x<-suppressWarnings(as.numeric(v)); x[v=="dead"]<-4; x}
y$bark<-ord(y$Bark.Missing); y$wound<-ord(y$Wounding.Holes)
d<-y[complete.cases(y[,c("bark","wound")]),]
d$cond<-cut(d$bark, c(0,1,2,4), labels=c("healthy","moderate","severe/dead"))
d$score<-d$bark+d$wound

# ---- Panel A: hump by bark condition ----
nlab<-as.data.frame(table(d$cond)); names(nlab)<-c("cond","n")
pA<-ggplot(d, aes(cond, fx, fill=cond)) +
  geom_hline(yintercept=0, linetype=3, colour="grey60") +
  geom_boxplot(width=.55, outlier.size=.7, alpha=.85) +
  geom_text(data=nlab, aes(cond, y=-2.2, label=paste0("n=",n)), inherit.aes=FALSE, size=3.2) +
  scale_fill_manual(values=c(healthy="#7fbf7b",moderate="#f0a202","severe/dead"="#8c510a"), guide="none") +
  labs(x="stem condition (bark loss)", y=expression(CH[4]~flux~(arcsinh)),
       title="A  Emission peaks at moderate deterioration") +
  theme_bw(base_size=11)

# ---- Panel B: composite score + species-controlled quadratic fit ----
d$sc<-scale(d$score, scale=FALSE)[,1]
m<-lmer(fx~sc+I(sc^2)+(1|sp), data=d); co<-summary(m)$coefficients
pq<-co["I(sc^2)","Pr(>|t|)"]
gr<-data.frame(score=seq(min(d$score),max(d$score),.1)); gr$sc<-gr$score-mean(d$score)
gr$fx<-co["(Intercept)","Estimate"]+co["sc","Estimate"]*gr$sc+co["I(sc^2)","Estimate"]*gr$sc^2
pB<-ggplot(d, aes(score, fx)) +
  geom_hline(yintercept=0, linetype=3, colour="grey60") +
  geom_jitter(width=.15, height=0, alpha=.35, size=1.4, colour="#8c510a") +
  geom_line(data=gr, linewidth=1, colour="black") +
  annotate("text", x=min(d$score), y=max(d$fx), hjust=0, vjust=1, size=3.3,
    label=sprintf("quadratic (hump): p=%.3f\nspecies-controlled", pq)) +
  scale_x_continuous(breaks=2:7) +
  labs(x="deterioration score (bark + wounding)", y=expression(CH[4]~flux~(arcsinh)),
       title="B  Inverted-U on a continuous deterioration axis") +
  theme_bw(base_size=11)

ggsave("outputs/revision/figS20_stem_deterioration.png", pA+pB+plot_layout(widths=c(1,1.15)),
       width=11, height=4.6, dpi=300)
cat(sprintf("quadratic p=%.3f | means by cond: %s\n", pq,
  paste(sprintf("%s=%.3f",levels(d$cond),tapply(d$CH4_best.flux,d$cond,mean)),collapse="  ")))
cat("wrote outputs/revision/figS20_stem_deterioration.png\n")
