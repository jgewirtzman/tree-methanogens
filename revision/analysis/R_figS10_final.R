# ==============================================================================
# REVISION — S10 rebuilt: pmoA/mmoX coupling & composition by COMPARTMENT.
# Four compartments (BrBG palette): Heartwood, Sapwood (wood) / Organic, Mineral
# (soil), each with its own fit. Core finding: the two methanotroph marker genes
# are DECOUPLED in wood but COUPLED in soil.
#   (a) Wood decoupling: pmoA vs mmoX by heartwood/sapwood — no co-variation.
#   (b) Soil decoupling: pmoA vs mmoX, points colored by SOIL MOISTURE (mean_vwc).
#       Genes co-vary (slope ~0.48), and the wet samples are the mmoX-dominated
#       cluster near the 1:1 line: wetter -> more mmoX (balance~VWC p~1e-4; mmoX~VWC
#       R2 0.29 vs pmoA~VWC 0.07). Toeslope/wetland-margin signal (sMMO Type II favored
#       in wet, organic, low-Cu soils). moisture joined by treecode+horizon (239/265).
#   (c) Composition: balance log(pmoA/mmoX) vs size 1/2[log pmoA + log mmoX], all
#       four compartments with individual fits + independence permutation (3000x).
#       NULL: shuffle pmoA and mmoX INDEPENDENTLY across samples (preserves each
#       gene's marginal distribution, breaks the within-sample pairing), recompute
#       the balance-size slope; two-sided p = deviation from that null. This is the
#       correct control for the built-in balance-vs-size geometry: when one gene has
#       restricted variance, a steep slope is mechanically expected. WOOD matches the
#       null (perm p 0.47/0.32) -> the ~1.4 slope IS that mechanical expectation
#       (mmoX = fixed low background, pmoA carries the variance), NOT co-regulation.
#       SOIL deviates from the null (perm p < 0.001) -> genuine coupling.
#   (d) Scaling vs total 16S (wood AND soil, faceted): pmoA & mmoX vs total 16S
#       rRNA. WOOD: pmoA weakly tracks community, mmoX is FLAT/negative (fixed
#       background) -> mechanism behind the wood decoupling. SOIL: BOTH genes track
#       total community strongly (mmoX most of all) -> mechanism behind soil coupling.
#       Uses X16S_per_ul (16S_tree_sample_table_with_meta.csv), the project-standard
#       total-16S normalizer and the ONLY 16S available for soil (ddPCR 16s_bact is
#       wood-only). 16S axis is a normalizer (copies/uL), so the SLOPE is the result
#       and is unit-invariant; complete-case (both genes detected), 16S >= 100.
#
# UNITS: copies g^-1 (= copies/uL x elution/mass; what is in the wood/soil). BASIS
#   DIFFERS BY MATERIAL: WOOD tube mass is a subsample of freeze-dried core -> DRY
#   basis; SOIL has no per-sample dry mass or moisture (only bulk mean_vwc), so the
#   tube mass is FIELD-MOIST -> FRESH basis (labeled as such; dry-basis harmonization
#   is pending and needs per-sample soil moisture we do not have). Conclusions are
#   basis-invariant: decoupling/composition are ratio- and slope-based, so a per-sample
#   mass scalar cancels; also robust to normalization (within-habitat mass CV 3-4%,
#   copies/uL near-identical). Only POOLING wood+soil is unsafe (cross-habitat mass
#   CV = 44%); habitats kept separate. Balances/ratios are unit-invariant regardless.
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

# join soil moisture (mean_vwc) for panel (b) recolor; key = treecode + horizon
vwc_lu <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv", row.names = 1, check.names = FALSE)
vwc_df <- tibble(joinkey = sub("[.]16S[.].*$", "", rownames(vwc_lu)), vwc = vwc_lu$mean_vwc) %>%
  filter(!is.na(vwc)) %>% distinct(joinkey, .keep_all = TRUE)
W <- W %>% mutate(joinkey = paste0(sub("\n.*", "", sub("^[A-Z]+_", "", sample_id)), core_type)) %>%
  left_join(vwc_df, by = "joinkey")

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
pa <- decouple_panel(WOODC, "dry")

# ---- (b) soil decoupling, recolored by soil moisture (mean_vwc) ---------------
soil <- dec %>% filter(comp %in% SOILC)
fs  <- lm(lm1 ~ lp1, soil)
fbv <- lm((lp1 - lm1) ~ vwc, soil %>% filter(!is.na(vwc)))
cat(sprintf("Soil moisture join: %d/%d matched | balance~VWC slope=%.3f p=%.1e | mmoX~VWC R2=%.2f pmoA~VWC R2=%.2f\n",
            sum(!is.na(soil$vwc)), nrow(soil), coef(fbv)[2], summary(fbv)$coef[2,4],
            summary(lm(lm1~vwc, soil))$r.squared, summary(lm(lp1~vwc, soil))$r.squared))
pb <- ggplot(soil, aes(lp1, lm1)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = vwc), alpha = 0.8, size = 1.7) +
  geom_smooth(method = "lm", formula = y~x, se = FALSE, color = "grey25", linewidth = 0.9) +
  scale_color_distiller(palette = "YlGnBu", direction = 1, na.value = "grey80", name = "soil moisture\n(% VWC)") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.3, size = 2.7, color = "grey15",
           label = sprintf("Soil (n=%d): slope=%.2f  R2=%.2f", nrow(soil), coef(fs)[2], summary(fs)$r.squared)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 2.75, size = 2.6, color = "#08519c",
           label = sprintf("wetter -> mmoX-dominated\n(balance~VWC p=%.0e)", summary(fbv)$coef[2,4])) +
  coord_fixed(ratio = 1, xlim = lim, ylim = lim) +
  labs(x = bquote(log[10]~"pmoA (copies g"^-1*" fresh + 1)"),
       y = bquote(log[10]~"mmoX (copies g"^-1*" fresh + 1)")) + th +
  theme(legend.position = c(0.985, 0.03), legend.justification = c(1, 0),
        legend.key.height = unit(0.30, "cm"), legend.key.width = unit(0.28, "cm"),
        legend.title = element_text(size = 7.3), legend.text = element_text(size = 6.8),
        legend.background = element_rect(fill = "white", color = "grey80"))

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

# ---- (d) scaling vs total 16S (wood AND soil, faceted) -----------------------
tm <- read.csv("data/raw/picrust/16S_tree_sample_table_with_meta.csv", row.names = 1, check.names = FALSE) %>%
  filter(material %in% c("Wood", "Soil"), pmoa_loose > 0, mmox_loose > 0, X16S_per_ul >= 100) %>%
  mutate(material = factor(material, levels = c("Wood", "Soil")),
         g16 = log10(X16S_per_ul), lp = log10(pmoa_loose), lmm = log10(mmox_loose))
sfit <- tm %>% group_by(material) %>% summarise(
  n = n(),
  ps = coef(lm(lp ~ g16))[2], pr = summary(lm(lp ~ g16))$r.squared,
  ms = coef(lm(lmm ~ g16))[2], mr = summary(lm(lmm ~ g16))$r.squared, .groups = "drop") %>%
  mutate(lab = sprintf("n=%d\npmoA: slope=%.2f  R2=%.2f\nmmoX: slope=%.2f  R2=%.2f", n, ps, pr, ms, mr))
longd <- tm %>% transmute(material, g16, `pmoA (pMMO)` = lp, `mmoX (sMMO)` = lmm) %>%
  pivot_longer(c(-material, -g16), names_to = "gene", values_to = "val")
pd <- ggplot(longd, aes(g16, val, color = gene)) +
  geom_point(alpha = 0.4, size = 1.1) +
  geom_smooth(method = "lm", formula = y~x, se = TRUE, linewidth = 0.9) +
  geom_text(data = sfit, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.07, vjust = 1.12, size = 2.35, lineheight = 0.95) +
  facet_wrap(~ material, scales = "free_x") +
  scale_color_manual(values = GENE, name = NULL) +
  labs(x = expression(log[10]~"total 16S rRNA (copies "*mu*"L"^-1*")"),
       y = expression(log[10]~"gene abundance")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"), legend.text = element_text(size = 8.3))

fig <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(out, "figS10_final.png"), fig, width = 10, height = 9.4, dpi = 300, bg = "white")
print(decstat); print(compstat)
cat("Wrote figS10_final.png\n")
