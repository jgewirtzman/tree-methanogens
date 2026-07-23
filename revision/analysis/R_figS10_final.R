# ==============================================================================
# REVISION — S10 rebuilt: pmoA/mmoX coupling & composition, WOOD vs SOIL.
# Core finding: the two methanotroph marker genes behave differently by habitat.
#   (a) Wood decoupling: pmoA vs mmoX — no co-variation (r2 ~ 0).
#   (b) Soil decoupling: pmoA vs mmoX — co-vary (r2 ~ 0.26); genes scale in tandem.
#   (c) Composition: balance log(pmoA/mmoX) vs size 1/2[log pmoA + log mmoX].
#       Wood slope ~1.4 (larger pools more pmoA-dominated); soil slope ~0 (stable).
#       Each slope tested against a 5000x independence permutation null: wood
#       matches the null (perm p~0.4) -> the shift reflects pmoA's larger dynamic
#       range (mmoX = fixed low background), NOT co-regulation; soil perm p~1.0
#       (isometric, composition invariant with abundance).
#   (d) Scaling vs total 16S (wood only; soil has no 16S): pmoA tracks community
#       abundance, mmoX is flat -> mechanism behind the wood decoupling.
#
# UNITS: copies uL^-1. This is a correlation/composition analysis, so copies uL^-1
#   is used deliberately: converting to copies g^-1 divides pmoA and mmoX by the
#   SAME sample mass and manufactures spurious co-variation (wood: r2 0.000 -> 0.019
#   from a 3% mass CV alone). Balances/ratios are unit-invariant. Absolute copies g^-1
#   dry abundances are reported separately (SI_fig_pmoa_mmox_separate + table).
# NONDETECTS: kept via log10(x+1) in the decoupling/scaling panels (wood ND 3-5%,
#   soil 0%). The composition panel needs both genes detected (a ratio with a zero
#   is undefined), so it is complete-case; this is stated.
# NOTE size axis = geometric-mean center (the 1/2 is the Bland-Altman mean); this is
#   NOT "total abundance" (that is the arithmetic sum pmoA+mmoX used elsewhere).
# Absolute axes subject to pending x10 (DILUTION_10X); ratios/slopes are not.
# NEW file; original S10 generator untouched. Output: revision/outputs/figS10_final.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
MAT  <- c(Wood = "#8c510a", Soil = "#01665e")
GENE <- c("pmoA (pMMO)" = "#2166ac", "mmoX (sMMO)" = "#b2182b")

d0 <- read.csv("data/compiled/ddpcr_gene_abundances.csv") %>%
  filter(analysis_type == "loose", target_gene %in% c("pmoa", "mmox", "16s_bact"))
meta <- d0 %>% group_by(sample_id) %>%
  summarise(material = first(material), mass = first(sample_mass_mg), .groups = "drop")
W <- d0 %>% select(sample_id, target_gene, conc = concentration_copies_per_uL, pos = positives) %>%
  pivot_wider(names_from = target_gene, values_from = c(conc, pos)) %>%
  left_join(meta, by = "sample_id") %>%
  filter(material %in% c("Wood", "Soil")) %>%
  mutate(material = factor(material, levels = c("Wood", "Soil")))

# ---- decoupling data: keep nondetects, +1, copies/uL -------------------------
dec <- W %>% mutate(pmoa = ifelse(is.na(conc_pmoa), 0, conc_pmoa),
                    mmox = ifelse(is.na(conc_mmox), 0, conc_mmox),
                    lp1 = log10(pmoa + 1), lm1 = log10(mmox + 1))
lim <- range(c(dec$lp1, dec$lm1))
decstat <- dec %>% group_by(material) %>%
  summarise(r2 = cor(lp1, lm1)^2, slope = coef(lm(lm1 ~ lp1))[2],
            p = summary(lm(lm1 ~ lp1))$coef[2, 4], n = n(), .groups = "drop")

# ---- composition data: complete-case (both detected), copies/uL --------------
cc <- W %>% filter(!is.na(pos_pmoa), pos_pmoa > 0, !is.na(pos_mmox), pos_mmox > 0) %>%
  mutate(lp = log10(conc_pmoa), lmm = log10(conc_mmox), size = 0.5*(lp + lmm), bal = lp - lmm)

perm_slope <- function(lp, lm, n = 5000) {
  m0  <- stats::lm((lp - lm) ~ I(0.5*(lp + lm)))
  o   <- unname(coef(m0)[2])
  nul <- replicate(n, { p <- sample(lp); m <- sample(lm); coef(stats::lm((p - m) ~ I(0.5*(p + m))))[2] })
  list(obs = o, ci = confint(m0)[2, ], r2 = summary(m0)$r.squared,
       p = summary(m0)$coef[2, 4], permp = mean(nul >= o))
}
set.seed(1)
pw <- perm_slope(cc$lp[cc$material == "Wood"], cc$lmm[cc$material == "Wood"])
ps <- perm_slope(cc$lp[cc$material == "Soil"], cc$lmm[cc$material == "Soil"])

th <- theme_bw(base_size = 11) +
  theme(aspect.ratio = 1,
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 8.3), panel.grid.minor = element_blank())
fmtp <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.2f", p))
fmtr <- function(r) ifelse(r < 0.001, "< 0.001", sprintf("= %.3f", r))

# ---- (a),(b) decoupling by habitat -------------------------------------------
decouple_panel <- function(mat) {
  dd <- dec %>% filter(material == mat); st <- decstat %>% filter(material == mat)
  ggplot(dd, aes(lp1, lm1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(color = MAT[[mat]], alpha = 0.55, size = 1.5) +
    geom_smooth(method = "lm", formula = y~x, se = FALSE, color = "grey20", linewidth = 0.9) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.35, size = 2.9, lineheight = 0.95,
             label = sprintf("%s (n = %d)\nslope = %.2f   R2 %s   p %s", mat, st$n, st$slope, fmtr(st$r2), fmtp(st$p))) +
    coord_fixed(ratio = 1, xlim = lim, ylim = lim) +
    labs(x = expression(log[10]~"pmoA (copies "*mu*"L"^-1*" + 1)"),
         y = expression(log[10]~"mmoX (copies "*mu*"L"^-1*" + 1)")) + th
}
pa <- decouple_panel("Wood"); pb <- decouple_panel("Soil")

# ---- (c) composition vs size, wood + soil overlaid ---------------------------
pc <- ggplot(cc, aes(size, bal, color = material)) +
  geom_smooth(method = "lm", formula = y~x, se = TRUE, linewidth = 0.9) +
  geom_point(alpha = 0.5, size = 1.4) +
  scale_color_manual(values = MAT, name = NULL) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 2.8, color = MAT[["Wood"]],
           label = sprintf("Wood: slope = %.2f  R2 = %.2f  perm p = %.2f", pw$obs, pw$r2, pw$permp)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 3.0, size = 2.8, color = MAT[["Soil"]],
           label = sprintf("Soil: slope = %.2f  R2 = %.2f  perm p = %.2f", ps$obs, ps$r2, ps$permp)) +
  labs(x = expression("size  "*frac(1,2)*"[log"[10]*" pmoA + log"[10]*" mmoX]"),
       y = expression("balance  log"[10]*"(pmoA / mmoX)")) +
  th + theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0))

# ---- (d) scaling vs independent 16S (wood only) ------------------------------
wc <- W %>% filter(material == "Wood", pos_pmoa > 0, pos_mmox > 0, !is.na(pos_16s_bact), pos_16s_bact > 0) %>%
  mutate(g16 = log10(conc_16s_bact), lp = log10(conc_pmoa), lmm = log10(conc_mmox))
fcp <- lm(lp ~ g16, wc); fcm <- lm(lmm ~ g16, wc)
long <- wc %>% transmute(g16, `pmoA (pMMO)` = lp, `mmoX (sMMO)` = lmm) %>%
  pivot_longer(-g16, names_to = "gene", values_to = "val")
pd <- ggplot(long, aes(g16, val, color = gene)) +
  geom_point(alpha = 0.5, size = 1.4) +
  geom_smooth(method = "lm", formula = y~x, se = TRUE, linewidth = 0.9) +
  scale_color_manual(values = GENE, name = NULL) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.04, vjust = -0.5, size = 2.8, lineheight = 0.95,
           label = sprintf("Wood (n = %d)\npmoA: slope=%.2f  R2=%.2f  p=%.0e\nmmoX: slope=%.2f  R2=%.2f  p=%.2f",
                           nrow(wc), coef(fcp)[2], summary(fcp)$r.squared, summary(fcp)$coef[2,4],
                           coef(fcm)[2], summary(fcm)$r.squared, summary(fcm)$coef[2,4])) +
  labs(x = expression(log[10]~"16S bacterial (copies "*mu*"L"^-1*")"),
       y = expression(log[10]~"gene (copies "*mu*"L"^-1*")")) +
  th + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

fig <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(out, "figS10_final.png"), fig, width = 10, height = 9.4, dpi = 300, bg = "white")
cat(sprintf("Wood: decouple r2=%.3f | comp slope=%.2f perm_p=%.2f || Soil: decouple r2=%.3f | comp slope=%.2f perm_p=%.2f\n",
            decstat$r2[decstat$material=="Wood"], pw$obs, pw$permp,
            decstat$r2[decstat$material=="Soil"], ps$obs, ps$permp))
cat("Wrote figS10_final.png\n")
