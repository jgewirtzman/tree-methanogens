source("code/lib/outputs.R")
# ==============================================================================
# REVISION — SI extended-isotope source composite (4 panels, 2x2).
# Design (Jon): CH4 Keeling      | CO2 Keeling
#               d13CH4 convergence | eps_C convergence
# (Miller-Tans dropped: leverage-biased, redundant with Keeling.)
# Data load replicated verbatim from ISOTOPES_final.R (settled QC). NEW file.
# Output: outputs/figures/generated/SI_isotopes_source_composite.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
out <- "outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
RED <- "#b2182b"; BLU <- "#2166ac"

# ---- data (verbatim from ISOTOPES_final.R) -----------------------------------
source("code/lib/isotope_samples.R")   # the one sample-selection rule (see that file)
# reference values READ from 12_isotopes-canonical.R; panels (c) and (d) typed -70 and
# 54-64 by hand, and the band outlived the reconciliation that moved eps_C to ~51-53
.IS <- with(read.csv("outputs/data/ISOTOPES_summary.csv"), setNames(value, quantity))
ATM_SRC <- .IS[["atm_corrected_source"]]
EPS_LO <- min(.IS[["eps_C_source_atm"]], .IS[["eps_C_source_keeling"]])
EPS_HI <- max(.IS[["eps_C_source_atm"]], .IS[["eps_C_source_keeling"]])
d <- isotope_whole_tree_samples()
wt <- d %>% filter(!is.na(d13CO2), co2_ppm > 0) %>%
  mutate(eps_C = ((d13CO2 + 1000)/(d13CH4 + 1000) - 1)*1000)

th <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
# each panel's result as a small in-panel label (titles are in the SI caption)
res_lab <- function(txt, right = FALSE)
  annotate("label", x = if (right) Inf else -Inf, y = Inf, label = txt, hjust = if (right) 1.03 else -0.03,
           vjust = 1.25, size = 3.4, fill = "white", alpha = 0.85, label.size = 0)
pt <- function(dat, x, y) ggplot(dat, aes({{x}}, {{y}})) +
  geom_point(alpha = 0.45, size = 1.2, color = "grey30")

# ---- source estimates --------------------------------------------------------
kc <- lm(d13CH4 ~ I(1/ch4_ppm), data = d);  s_kc <- coef(kc)[1];  ci_kc <- confint(kc)[1,]
ko <- lm(d13CO2 ~ I(1/co2_ppm), data = wt); s_ko <- coef(ko)[1];  ci_ko <- confint(ko)[1,]

# ---- (a) CH4 Keeling ---------------------------------------------------------
pa <- pt(d, 1/ch4_ppm, d13CH4) +
  geom_smooth(method = "lm", formula = y~x, color = "black", fill = "grey80", linewidth = 0.7) +
  annotate("point", x = 0, y = s_kc, color = RED, size = 2.8) +
  res_lab(sprintf("Intercept %.0f\u2030 (95%% CI %.0f to %.0f)", s_kc, ci_kc[1], ci_kc[2]), right = TRUE) +
  labs(x = expression(1/CH[4]~"(ppm"^-1*")"), y = expression(delta^13*"C-CH"[4]~"(\u2030)")) + th
# ---- d13CH4 convergence ------------------------------------------------------
pc <- pt(d, ch4_ppm, d13CH4) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", linewidth = 0.7) +
  geom_hline(yintercept = ATM_SRC, linetype = "dashed", color = RED, linewidth = 0.6) +
  scale_x_log10() +
  res_lab(sprintf("Plateau %.0f\u2030 (two-endmember mixing)", ATM_SRC), right = TRUE) +
  labs(x = expression(CH[4]~"(ppm)"), y = expression(delta^13*"C-CH"[4]~"(\u2030)")) + th
# ---- (d) CO2 Keeling ---------------------------------------------------------
pd <- pt(wt, 1/co2_ppm, d13CO2) +
  geom_smooth(method = "lm", formula = y~x, color = "black", fill = "grey80", linewidth = 0.7) +
  annotate("point", x = 0, y = s_ko, color = BLU, size = 2.8) +
  res_lab(sprintf("Intercept %.0f\u2030 (95%% CI %.0f to %.0f)", s_ko, ci_ko[1], ci_ko[2])) +
  labs(x = expression(1/CO[2]~"(ppm"^-1*")"), y = expression(delta^13*"C-CO"[2]~"(\u2030)")) + th
# ---- eps_C convergence -------------------------------------------------------
pf <- pt(wt, ch4_ppm, eps_C) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = EPS_LO, ymax = EPS_HI, fill = RED, alpha = 0.10) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", linewidth = 0.7) +
  scale_x_log10() +
  res_lab(sprintf("Source-based \u03b5C %.0f\u2013%.0f\u2030", EPS_LO, EPS_HI), right = TRUE) +
  labs(x = expression(CH[4]~"(ppm)"), y = expression(epsilon[C]~"(\u2030)")) + th

fig <- (pa | pd) / (pc | pf) +          # top row = Keeling source plots; bottom = convergence
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(out_path("SI_isotopes_source_composite.png"), fig, width = 9.5, height = 8, dpi = 300, bg = "white")
cat(sprintf("CH4 Keeling %.0f | CO2 Keeling %.0f | source-based eps_C = %.0f | n=%d\n",
            s_kc, s_ko, s_ko - s_kc, nrow(wt)))
cat("Wrote SI_isotopes_source_composite.png\n")
