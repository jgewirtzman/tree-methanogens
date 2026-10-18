# ==============================================================================
# figS08_pmoa-mmox.R -- Figure S9: pmoA and mmoX by compartment, and their coupling
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) pmoA and mmoX by compartment: median and IQR, copies g-1 (non-detects as 0)
#   (b, c) coupling: log10(mmoX + 1) against log10(pmoA + 1), wood (heartwood, sapwood)
#       and soil (organic, mineral), one linear fit per compartment
#   (d, e) composition: balance log10(pmoA/mmoX) against size
#       1/2 (log10 pmoA + log10 mmoX), wood and soil, samples with both genes detected
# PERMUTATION NULL for (d, e): pmoA and mmoX are shuffled INDEPENDENTLY across samples
# (each gene's marginal kept, the within-sample pairing broken) and the balance-size
# slope recomputed, 3000 times; two-sided p = deviation from that null. A steep slope
# is expected mechanically when one gene has restricted variance, so only a deviation
# from the null indicates coupling.
# Units: copies g-1, dry mass for wood, fresh mass for soil (data/tables/ddpcr.csv,
# code/lib/ddpcr_constants.R); slopes and ratios are unit-free. Wood and soil are never
# pooled. The size axis is the geometric-mean centre, not the arithmetic total
# abundance used elsewhere.
# Samples: 2021 DNA samples with both genes measured (dna_samples(), mean over a
#   sample's runs that passed QC); a sample with one gene missing is left out.
# Outputs: outputs/figures/generated/figS10_final.png, outputs/data/pmoa_mmox_by_compartment.csv
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
source("code/lib/outputs.R")
source("code/lib/ddpcr_constants.R")
COMP4 <- c(Heartwood = "#a6611a", Sapwood = "#dfc27d", Organic = "#018571", Mineral = "#80cdc1")
WOODC <- c("Heartwood", "Sapwood"); SOILC <- c("Organic", "Mineral")

source("code/lib/tables_read.R")
W <- dna_samples("2021_survey") %>%
  filter(!is.na(pmoa_g), !is.na(mmox_g)) %>%
  mutate(comp = factor(recode(core_type, Inner = "Heartwood", Outer = "Sapwood"), levels = names(COMP4))) %>%
  filter(!is.na(comp))

# ---- decoupling data: keep nondetects (0 copies), +1, copies/g ----------------
dec <- W %>% mutate(lp1 = log10(pmoa_g + 1), lm1 = log10(mmox_g + 1))
lim <- range(c(dec$lp1, dec$lm1))
decstat <- dec %>% group_by(comp) %>%
  summarise(slope = coef(lm(lm1 ~ lp1))[2], r2 = cor(lp1, lm1)^2,
            p = summary(lm(lm1 ~ lp1))$coef[2, 4], n = n(), .groups = "drop")

# ---- composition data: complete-case (both detected), copies/g ---------------
cc <- W %>% filter(pmoa_g > 0, mmox_g > 0) %>%
  mutate(lp = log10(pmoa_g), lmm = log10(mmox_g), size = 0.5*(lp + lmm), bal = lp - lmm)
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

# ---- (b),(c) coupling by compartment ------------------------------------------
decouple_panel <- function(comps, basis) {
  dd <- dec %>% filter(comp %in% comps)
  p <- ggplot(dd, aes(lp1, lm1, color = comp)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.55, size = 1.5) +
    geom_smooth(method = "lm", formula = y~x, se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = COMP4[comps], name = NULL)
  for (i in seq_along(comps)) {
    k <- comps[i]; r <- decstat[decstat$comp == k, ]
    p <- p + annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.3 + (i-1)*1.7, size = 3.3, color = COMP4[[k]],
                      label = gsub("([=<])(?=[0-9])", "\\1 ", sprintf("%s (n = %d): slope %.2f, R\u00b2 %s, p %s", k, r$n, r$slope, fmtr(r$r2), fmtp(r$p)), perl = TRUE))
  }
  p + coord_fixed(ratio = 1, xlim = lim, ylim = lim) +
    labs(x = bquote(log[10]~italic(pmoA)~"(copies g"^-1*" "*.(basis)*" + 1)"),
         y = bquote(log[10]~italic(mmoX)~"(copies g"^-1*" "*.(basis)*" + 1)")) + th + theme(legend.position = "none")
}

# ---- (d),(e) composition by compartment --------------------------------------
comp_panel <- function(comps) {
  d <- cc %>% filter(comp %in% comps); st <- compstat %>% filter(comp %in% comps)
  p <- ggplot(d, aes(size, bal, color = comp)) +
    geom_smooth(method = "lm", formula = y~x, se = TRUE, linewidth = 0.9) +
    geom_point(alpha = 0.5, size = 1.4) +
    scale_color_manual(values = COMP4[comps], name = NULL)
  for (i in seq_len(nrow(st))) {
    r <- st[i, ]
    p <- p + annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3 + (i-1)*1.6, size = 3.3, color = COMP4[[r$comp]],
                      label = gsub("([=<])(?=[0-9])", "\\1 ", sprintf("%s: slope %.2f, R\u00b2 = %.2f, permutation p %s", r$comp, r$slope, r$r2, permlab(r$permp)), perl = TRUE))
  }
  p + labs(x = expression("Size, "*frac(1,2)*"(log"[10]*" "*italic(pmoA)*" + log"[10]*" "*italic(mmoX)*")"),
           y = expression("log"[10]*"("*italic(pmoA)*" / "*italic(mmoX)*")")) + th + theme(legend.position = "none")
}

# ---- (a) pmoA and mmoX separately, by compartment (median + IQR) -------------
# Same samples and copies/g as panels (b)-(e); non-detects enter as zero.
sep <- dec %>% select(comp, pmoA = pmoa_g, mmoX = mmox_g) %>%
  pivot_longer(c(pmoA, mmoX), names_to = "gene", values_to = "cg") %>%
  group_by(comp, gene) %>%
  summarise(med = median(cg), lo = quantile(cg, .25), hi = quantile(cg, .75), n = n(), .groups = "drop") %>%
  mutate(across(c(med, lo, hi), ~ pmax(.x, 1)),
         gene = factor(gene, levels = c("pmoA", "mmoX"), labels = c("pmoA (pMMO)", "mmoX (sMMO)")),
         comp = factor(comp, levels = names(COMP4)))
p_sep <- ggplot(sep, aes(comp, med, color = comp, shape = gene)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5),
                  size = 0.55, linewidth = 0.8) +
  scale_color_manual(values = COMP4, guide = "none") +
  scale_shape_manual(values = c("pmoA (pMMO)" = 16, "mmoX (sMMO)" = 17), name = NULL,
                     labels = c(expression(italic(pmoA)*" (pMMO)"), expression(italic(mmoX)*" (sMMO)"))) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  labs(x = NULL, y = expression("Copies g"^-1*" (median, IQR)")) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "top")
write.csv(sep, out_path("pmoa_mmox_by_compartment.csv"), row.names = FALSE)

pa <- p_sep; pb0 <- decouple_panel(WOODC, "dry"); pc0 <- decouple_panel(SOILC, "fresh")
pd0 <- comp_panel(WOODC); pe0 <- comp_panel(SOILC)
fig <- pa / ((pb0 | pc0) / (pd0 | pe0)) + plot_layout(heights = c(0.45, 2)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(out_path("figS10_final.png"), fig, width = 9.4, height = 12.2, dpi = 300, bg = "white")
print(decstat); print(compstat[, c("comp","slope","r2","permp","n")])
cat("Wrote figS10_final.png\n")
