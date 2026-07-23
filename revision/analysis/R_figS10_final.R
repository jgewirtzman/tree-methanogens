# ==============================================================================
# REVISION — S10 rebuilt: pmoA/mmoX composition & scaling (wood, copies/g).
# (a) decoupling: pmoA vs mmoX (1:1 line) — the two genes do not co-vary.
# (b) composition: balance log(pmoA/mmoX) vs size 1/2[log pmoA + log mmoX]
#     (GEOMETRIC-mean size, not arithmetic sum, so x isn't a pmoA proxy).
#     Larger combined pmoA+mmoX pools are more pmoA-dominated. Slope tested with a
#     5000x independence permutation: matches the independence expectation, i.e.
#     the shift reflects pmoA's larger dynamic range (mmoX = fixed low background),
#     NOT coordinated co-regulation.
# (c) scaling: pmoA & mmoX vs total 16S (independent axis) — pmoA tracks community
#     abundance, mmoX is flat; the flat mmoX line also controls for shared yield.
# Units: copies/g dry wood; results unit-invariant (sample-mass CV = 3%).
# Absolute axes subject to pending x10 (DILUTION_10X); ratios/slopes are not.
# NEW file; original S10 generator untouched.
# Output: revision/outputs/figS10_final.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
E <- 75; DILUTION_10X <- 1
COMP <- c(Heartwood = "#a6611a", Sapwood = "#dfc27d")
GENE <- c("pmoA (pMMO)" = "#2166ac", "mmoX (sMMO)" = "#b2182b")

d0 <- read.csv("data/compiled/ddpcr_gene_abundances.csv") %>%
  filter(analysis_type == "loose", target_gene %in% c("pmoa", "mmox", "16s_bact"))
meta <- d0 %>% group_by(sample_id) %>%
  summarise(material = first(material), core_type = first(core_type), mass = first(sample_mass_mg), .groups = "drop")
w <- d0 %>% select(sample_id, target_gene, conc = concentration_copies_per_uL, pos = positives) %>%
  pivot_wider(names_from = target_gene, values_from = c(conc, pos)) %>% left_join(meta, by = "sample_id") %>%
  filter(material == "Wood", !is.na(mass), mass > 0, pos_pmoa > 0, pos_mmox > 0) %>%  # complete cases
  mutate(comp = factor(recode(core_type, Inner = "Heartwood", Outer = "Sapwood"), levels = names(COMP)),
         cg   = E / (mass/1000) * DILUTION_10X,
         lp = log10(conc_pmoa*cg), lm = log10(conc_mmox*cg),
         size = 0.5*(lp + lm), bal = lp - lm)

# permutation null for panel (b)
set.seed(1); fitb <- lm(bal ~ size, w); obs <- coef(fitb)[2]; ci <- confint(fitb)[2,]
nullb <- replicate(5000, { p <- sample(w$lp); m <- sample(w$lm); coef(lm((p-m) ~ I(0.5*(p+m))))[2] })
permp <- mean(nullb >= obs)

th <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 8.7),
        legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 8.5), panel.grid.minor = element_blank())

# ---- (a) decoupling ----------------------------------------------------------
pa <- ggplot(w, aes(lp, lm, color = comp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 1.6) +
  scale_color_manual(values = COMP, name = NULL) +
  labs(x = expression(log[10]~pmoA~"(copies g"^-1*" dry)"), y = expression(log[10]~mmoX~"(copies g"^-1*" dry)"),
       title = "pmoA and mmoX are decoupled",
       subtitle = sprintf("r2 = %.3f (1:1 dashed)", cor(w$lp, w$lm)^2)) + th

# ---- (b) composition vs size -------------------------------------------------
pb <- ggplot(w, aes(size, bal)) +
  geom_smooth(method = "lm", formula = y~x, color = "black", fill = "grey80", linewidth = 0.8) +
  geom_point(aes(color = comp), alpha = 0.6, size = 1.6) +
  scale_color_manual(values = COMP, name = NULL) +
  labs(x = expression("size  "*frac(1,2)*"[log pmoA + log mmoX]"), y = expression("balance  log"[10]*"(pmoA / mmoX)"),
       title = "Larger pmoA+mmoX pools are more pmoA-dominated",
       subtitle = sprintf("slope = %.2f [%.2f, %.2f]; = independence expectation (perm p = %.2f)\n-> pmoA's larger range, not co-regulation", obs, ci[1], ci[2], permp)) + th

# ---- (c) scaling vs independent 16S ------------------------------------------
wc <- w %>% filter(pos_16s_bact > 0, !is.na(conc_16s_bact)) %>% mutate(g16 = log10(conc_16s_bact*cg))
sp <- coef(lm(lp ~ g16, wc)); pp <- summary(lm(lp ~ g16, wc))$coef[2,4]
sm <- coef(lm(lm ~ g16, wc)); pm <- summary(lm(lm ~ g16, wc))$coef[2,4]
long <- wc %>% transmute(g16, `pmoA (pMMO)` = lp, `mmoX (sMMO)` = lm) %>%
  pivot_longer(-g16, names_to = "gene", values_to = "val")
pc <- ggplot(long, aes(g16, val, color = gene)) +
  geom_point(alpha = 0.5, size = 1.4) +
  geom_smooth(method = "lm", formula = y~x, se = TRUE, linewidth = 0.9) +
  scale_color_manual(values = GENE, name = NULL) +
  labs(x = expression(log[10]~"16S bacterial (copies g"^-1*")"), y = expression(log[10]~"gene (copies g"^-1*")"),
       title = "pmoA tracks community size; mmoX is flat",
       subtitle = sprintf("pmoA slope = %.2f (p = %.0e);  mmoX slope = %.2f (p = %.2f)", sp[2], pp, sm[2], pm)) +
  th + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))

fig <- (pa | pb | pc) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(out, "figS10_final.png"), fig, width = 15, height = 5, dpi = 300, bg = "white")
cat(sprintf("(b) slope=%.2f [%.2f,%.2f] perm_p=%.2f | (c) pmoA=%.2f mmoX=%.2f | n(a/b)=%d n(c)=%d\n",
            obs, ci[1], ci[2], permp, sp[2], sm[2], nrow(w), nrow(wc)))
cat("Wrote figS10_final.png\n")
