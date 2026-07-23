# ==============================================================================
# REVISION — S10 rebuilt: pmoA/mmoX coupling & composition by COMPARTMENT.
# Four compartments (BrBG palette): Heartwood, Sapwood (wood) / Organic, Mineral
# (soil), each with its own fit. Core finding: the two methanotroph marker genes
# are DECOUPLED in wood but COUPLED in soil.
#   (a) Wood decoupling: pmoA vs mmoX by heartwood/sapwood — no co-variation.
#   (b) Soil decoupling: pmoA vs mmoX by organic/mineral — co-vary.
#   (c) Composition: balance log(pmoA/mmoX) vs size 1/2[log pmoA + log mmoX], all
#       four compartments with individual fits + 3000x independence permutation per
#       compartment. Wood compartments steep (perm p ~ 0.4 -> matches independence
#       null: mmoX = fixed low background, pmoA carries dynamics); soil compartments
#       flat (perm p ~ 1 -> isometric, composition invariant with abundance).
#   (d) Scaling vs total 16S (wood only; soil has no 16S): pmoA tracks community
#       abundance, mmoX is flat -> mechanism behind the wood decoupling.
#
# UNITS: copies g^-1 dry (environmental concentration = copies/uL x elution/mass;
#   what is actually in the wood/soil, not the eluate concentration). Robust to
#   normalization: within-habitat sample-mass CV is only 3-4%, so copies/uL gives
#   near-identical structure -> shared-mass normalization is NOT the driver. (Only
#   POOLING wood+soil is unsafe: cross-habitat mass CV = 44%; habitats kept separate.)
#   Balances/ratios are unit-invariant regardless.
# NONDETECTS: kept via log10(x+1) in the decoupling panels (wood ND 3-5%, soil 0%).
#   The composition/scaling panels need genes detected (a ratio with a zero is
#   undefined), so they are complete-case; this is stated.
# NOTE size axis = geometric-mean center (the 1/2 is the Bland-Altman mean); NOT
#   "total abundance" (that is the arithmetic sum pmoA+mmoX used elsewhere).
# Absolute axes subject to pending x10 (DILUTION_10X); ratios/slopes are not.
# NEW file; original S10 generator untouched. Output: revision/outputs/figS10_final.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
E <- 75; DILUTION_10X <- 1
COMP4 <- c(Heartwood = "#a6611a", Sapwood = "#dfc27d", Organic = "#018571", Mineral = "#80cdc1")
GENE  <- c("pmoA (pMMO)" = "#2166ac", "mmoX (sMMO)" = "#b2182b")
WOODC <- c("Heartwood", "Sapwood"); SOILC <- c("Organic", "Mineral")

d0 <- read.csv("data/compiled/ddpcr_gene_abundances.csv") %>%
  filter(analysis_type == "loose", target_gene %in% c("pmoa", "mmox", "16s_bact"))
meta <- d0 %>% group_by(sample_id) %>%
  summarise(material = first(material), core_type = first(core_type), mass = first(sample_mass_mg), .groups = "drop")
W <- d0 %>% select(sample_id, target_gene, conc = concentration_copies_per_uL, pos = positives) %>%
  pivot_wider(names_from = target_gene, values_from = c(conc, pos)) %>%
  left_join(meta, by = "sample_id") %>%
  filter(material %in% c("Wood", "Soil"), !is.na(mass), mass > 0) %>%
  mutate(comp = factor(recode(core_type, Inner = "Heartwood", Outer = "Sapwood"), levels = names(COMP4)),
         cg = E / (mass/1000) * DILUTION_10X) %>%
  filter(!is.na(comp))

# ---- decoupling data: keep nondetects, +1, copies/g --------------------------
dec <- W %>% mutate(pmoa_g = ifelse(is.na(conc_pmoa), 0, conc_pmoa) * cg,
                    mmox_g = ifelse(is.na(conc_mmox), 0, conc_mmox) * cg,
                    lp1 = log10(pmoa_g + 1), lm1 = log10(mmox_g + 1))
lim <- range(c(dec$lp1, dec$lm1))
decstat <- dec %>% group_by(comp) %>%
  summarise(slope = coef(lm(lm1 ~ lp1))[2], r2 = cor(lp1, lm1)^2,
            p = summary(lm(lm1 ~ lp1))$coef[2, 4], n = n(), .groups = "drop")

# ---- composition data: complete-case (both detected), copies/g ---------------
cc <- W %>% filter(pos_pmoa > 0, pos_mmox > 0, !is.na(conc_pmoa), !is.na(conc_mmox)) %>%
  mutate(lp = log10(conc_pmoa * cg), lmm = log10(conc_mmox * cg), size = 0.5*(lp + lmm), bal = lp - lmm)
set.seed(1)
compstat <- lapply(names(COMP4), function(k) {
  s <- cc %>% filter(comp == k); m0 <- stats::lm(bal ~ size, s); o <- unname(coef(m0)[2])
  nul <- replicate(3000, { p <- sample(s$lp); m <- sample(s$lmm); coef(stats::lm((p - m) ~ I(0.5*(p + m))))[2] })
  permp <- min(1, 2 * min(mean(nul >= o), mean(nul <= o)))  # two-sided: deviation from independence
  data.frame(comp = k, slope = o, r2 = summary(m0)$r.squared, p = summary(m0)$coef[2, 4], permp = permp, n = nrow(s))
}) %>% bind_rows()

th <- theme_bw(base_size = 11) +
  theme(aspect.ratio = 1,
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 8.3), panel.grid.minor = element_blank())
fmtp <- function(p) ifelse(p < 0.001, "<0.001", sprintf("=%.2f", p))
fmtr <- function(r) ifelse(r < 0.001, "<0.001", sprintf("=%.3f", r))

# ---- (a),(b) decoupling by compartment ---------------------------------------
decouple_panel <- function(comps) {
  dd <- dec %>% filter(comp %in% comps)
  p <- ggplot(dd, aes(lp1, lm1, color = comp)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.55, size = 1.5) +
    geom_smooth(method = "lm", formula = y~x, se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = COMP4[comps], name = NULL)
  for (i in seq_along(comps)) {
    k <- comps[i]; r <- decstat[decstat$comp == k, ]
    p <- p + annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.3 + (i-1)*1.7, size = 2.7, color = COMP4[[k]],
                      label = sprintf("%s (n=%d): slope=%.2f  R2%s  p%s", k, r$n, r$slope, fmtr(r$r2), fmtp(r$p)))
  }
  p + coord_fixed(ratio = 1, xlim = lim, ylim = lim) +
    labs(x = expression(log[10]~"pmoA (copies g"^-1*" dry + 1)"),
         y = expression(log[10]~"mmoX (copies g"^-1*" dry + 1)")) + th + theme(legend.position = "none")
}
pa <- decouple_panel(WOODC); pb <- decouple_panel(SOILC)

# ---- (c) composition vs size, all four compartments --------------------------
pc <- ggplot(cc, aes(size, bal, color = comp)) +
  geom_smooth(method = "lm", formula = y~x, se = TRUE, linewidth = 0.9) +
  geom_point(alpha = 0.5, size = 1.3) +
  scale_color_manual(values = COMP4, name = NULL)
for (i in seq_len(nrow(compstat))) {
  r <- compstat[i, ]
  pc <- pc + annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.3 + (i-1)*1.5, size = 2.45, color = COMP4[[r$comp]],
                      label = sprintf("%s: slope=%.2f  perm p %s", r$comp, r$slope, ifelse(r$permp < 0.001, "<0.001", sprintf("=%.2f", r$permp))))
}
pc <- pc + labs(x = expression("size  "*frac(1,2)*"[log"[10]*" pmoA + log"[10]*" mmoX]"),
                y = expression("balance  log"[10]*"(pmoA / mmoX)")) + th + theme(legend.position = "none")

# ---- (d) scaling vs independent 16S (wood only) ------------------------------
wc <- W %>% filter(material == "Wood", pos_pmoa > 0, pos_mmox > 0, !is.na(pos_16s_bact), pos_16s_bact > 0) %>%
  mutate(g16 = log10(conc_16s_bact * cg), lp = log10(conc_pmoa * cg), lmm = log10(conc_mmox * cg))
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
  labs(x = expression(log[10]~"16S bacterial (copies g"^-1*" dry)"),
       y = expression(log[10]~"gene (copies g"^-1*" dry)")) +
  th + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

fig <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(out, "figS10_final.png"), fig, width = 10, height = 9.4, dpi = 300, bg = "white")
print(decstat); print(compstat)
cat("Wrote figS10_final.png\n")
