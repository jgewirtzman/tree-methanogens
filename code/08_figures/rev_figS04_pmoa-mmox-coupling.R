# ==============================================================================
# REVISION — S10: pmoA/mmoX coupling & composition, WOOD (focus) vs SOIL (contrast).
# Tree-focused; four compartments (BrBG palette): Heartwood, Sapwood (wood) /
# Organic, Mineral (soil), each with its own fit.
#   (a) Wood decoupling: pmoA vs mmoX (heartwood/sapwood) — no co-variation (r2~0).
#   (b) Soil decoupling: pmoA vs mmoX (organic/mineral) — co-vary.
#   (c) Wood composition: balance log(pmoA/mmoX) vs size 1/2[log pmoA + log mmoX],
#       heartwood/sapwood with individual fits + independence permutation.
#   (d) Soil composition: same axes for organic/mineral.
#       PERMUTATION NULL: shuffle pmoA and mmoX INDEPENDENTLY across samples (preserves
#       each gene's marginal, breaks within-sample pairing), recompute the balance-size
#       slope; two-sided p = deviation from that null. This controls the built-in
#       balance-vs-size geometry (a steep slope is mechanically expected when one gene
#       has restricted variance). WOOD matches the null (perm p 0.47/0.32) -> the ~1.4
#       slope IS that mechanical expectation (mmoX = fixed low background, pmoA carries
#       the variance), NOT co-regulation. SOIL deviates (perm p < 0.001) -> real coupling.
#
# UNITS: copies g^-1 (= copies/uL x elution/mass). BASIS DIFFERS BY MATERIAL: WOOD tube
#   mass is a subsample of freeze-dried core -> DRY; SOIL has no per-sample dry mass or
#   moisture -> FIELD-MOIST -> FRESH (labeled as such; dry-basis harmonization pending).
#   Conclusions basis-invariant (ratio/slope-based; per-sample mass scalar cancels) and
#   normalization-robust (within-habitat mass CV 3-4%). Do NOT pool wood+soil (cross-
#   habitat CV 44%); habitats kept separate. Balances/ratios are unit-invariant.
# NONDETECTS: kept via log10(x+1) in the decoupling panels (wood ND 3-5%, soil 0%).
#   Composition panels need both genes detected (ratio undefined at a zero) -> complete-case.
# NOTE size axis = geometric-mean center (the 1/2 is the Bland-Altman mean); NOT the
#   arithmetic-sum "total abundance" used elsewhere.
# Absolute axes subject to pending x10 (DILUTION_10X); ratios/slopes are not.
# NEW file; original S10 generator untouched. Output: outputs/revision/figS10_final.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
E <- 75; DILUTION_10X <- 1
COMP4 <- c(Heartwood = "#a6611a", Sapwood = "#dfc27d", Organic = "#018571", Mineral = "#80cdc1")
WOODC <- c("Heartwood", "Sapwood"); SOILC <- c("Organic", "Mineral")

d0 <- read.csv("data/compiled/ddpcr_gene_abundances.csv") %>%
  filter(analysis_type == "loose", target_gene %in% c("pmoa", "mmox"))
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
  theme(aspect.ratio = 1, panel.grid.minor = element_blank())
fmtp <- function(p) ifelse(p < 0.001, "<0.001", sprintf("=%.2f", p))
fmtr <- function(r) ifelse(r < 0.001, "<0.001", sprintf("=%.3f", r))
permlab <- function(p) ifelse(p < 0.001, "<0.001", sprintf("=%.2f", p))

# ---- (a),(b) decoupling by compartment ---------------------------------------
decouple_panel <- function(comps, basis) {
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
    labs(x = bquote(log[10]~"pmoA (copies g"^-1*" "*.(basis)*" + 1)"),
         y = bquote(log[10]~"mmoX (copies g"^-1*" "*.(basis)*" + 1)")) + th + theme(legend.position = "none")
}
pa <- decouple_panel(WOODC, "dry"); pb <- decouple_panel(SOILC, "fresh")

# ---- (c),(d) composition by compartment --------------------------------------
comp_panel <- function(comps) {
  d <- cc %>% filter(comp %in% comps); st <- compstat %>% filter(comp %in% comps)
  p <- ggplot(d, aes(size, bal, color = comp)) +
    geom_smooth(method = "lm", formula = y~x, se = TRUE, linewidth = 0.9) +
    geom_point(alpha = 0.5, size = 1.4) +
    scale_color_manual(values = COMP4[comps], name = NULL)
  for (i in seq_len(nrow(st))) {
    r <- st[i, ]
    p <- p + annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3 + (i-1)*1.6, size = 2.6, color = COMP4[[r$comp]],
                      label = sprintf("%s: slope=%.2f  R2=%.2f  perm p %s", r$comp, r$slope, r$r2, permlab(r$permp)))
  }
  p + labs(x = expression("size  "*frac(1,2)*"[log"[10]*" pmoA + log"[10]*" mmoX]"),
           y = expression("balance  log"[10]*"(pmoA / mmoX)")) + th + theme(legend.position = "none")
}
pc <- comp_panel(WOODC); pd <- comp_panel(SOILC)

fig <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(out, "figS10_final.png"), fig, width = 9.4, height = 9.4, dpi = 300, bg = "white")
print(decstat); print(compstat[, c("comp","slope","r2","permp","n")])
cat("Wrote figS10_final.png\n")
