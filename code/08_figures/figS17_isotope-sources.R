# ==============================================================================
# figS17_isotope-sources.R -- Figure S18: source signatures of internal CH4 and CO2
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) CH4 Keeling plot, MM robust regression (as 09_tables_stats/12_isotopes-canonical.R);
#       open points: robustness weight < 0.5
#   (b) CO2 Keeling plot, OLS
#   (c) d13C-CH4 against CH4 (loess), with the atmosphere-corrected source
#   (d) eps_C against CH4 (loess), with the band of the source-based eps_C estimates
# Inputs: data/tables/isotopes.csv (isotope_whole_tree(): the whole-tree set, rule in
#   code/lib/isotope_samples.R); outputs/data/ISOTOPES_summary.csv (12_isotopes-canonical.R)
# Output: outputs/figures/generated/SI_isotopes_source_composite.png
# SI house style: code/lib/fig_style.R
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
source("code/lib/outputs.R")
source("code/lib/fig_style.R")
RED <- "#b2182b"; BLU <- "#2166ac"

# ---- data ---------------------------------------------------------------------
source("code/lib/tables_read.R")
# reference values read from 12_isotopes-canonical.R, never typed
.IS <- with(read.csv("outputs/data/ISOTOPES_summary.csv"), setNames(value, quantity))
ATM_SRC <- .IS[["atm_corrected_source"]]
EPS_LO <- min(.IS[["eps_C_source_atm"]], .IS[["eps_C_source_keeling"]])
EPS_HI <- max(.IS[["eps_C_source_atm"]], .IS[["eps_C_source_keeling"]])
d <- isotope_whole_tree()   # species from trees.csv
wt <- d %>% filter(!is.na(d13CO2), co2_ppm > 0) %>%
  mutate(eps_C = ((d13CO2 + 1000)/(d13CH4 + 1000) - 1)*1000)

th <- theme_si()
y_num <- scale_y_continuous(labels = scales::label_number(style_negative = "minus"))
neg <- function(x) sub("^-", "\u2212", sprintf("%.0f", x))   # true minus sign in annotations
# each panel's result as a small in-panel label (titles are in the SI caption): plain
# text on a borderless white backing so it stays legible over points
res_lab <- function(txt, right = FALSE, parse = FALSE)
  annotate("label", x = if (right) Inf else -Inf, y = Inf, label = txt, hjust = if (right) 1.02 else -0.02,
           vjust = 1.2, size = 3, fill = "white", alpha = 0.85, border.colour = NA, parse = parse)
pt <- function(dat, x, y) ggplot(dat, aes({{x}}, {{y}})) +
  geom_point(alpha = 0.45, size = 1.2, color = "grey30")

# ---- source estimates --------------------------------------------------------
kc <- MASS::rlm(d13CH4 ~ I(1/ch4_ppm), data = d, method = "MM", maxit = 100)   # as 12_isotopes-canonical.R
s_kc <- coef(kc)[1]; ci_kc <- s_kc + c(-1, 1) * qnorm(0.975) * summary(kc)$coefficients[1, 2]
d$mm_w <- kc$w
ko <- lm(d13CO2 ~ I(1/co2_ppm), data = wt); s_ko <- coef(ko)[1];  ci_ko <- confint(ko)[1,]

# ---- (a) CH4 Keeling ---------------------------------------------------------
pa <- ggplot(d, aes(1/ch4_ppm, d13CH4)) +
  geom_point(aes(shape = mm_w < 0.5), alpha = 0.55, size = 1.3, color = "grey30") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  geom_abline(intercept = coef(kc)[1], slope = coef(kc)[2], color = "black", linewidth = 0.7) +
  annotate("point", x = 0, y = s_kc, color = RED, size = 2.8) +
  res_lab(sprintf("Intercept %s\u2030 (95%% CI %s to %s)", neg(s_kc), neg(ci_kc[1]), neg(ci_kc[2])), right = TRUE) +
  y_num + labs(x = expression(1/CH[4]~(ppm^{-1})), y = expression(delta^{13}*CH[4]~"(\u2030)")) + th
# ---- (c) d13CH4 convergence --------------------------------------------------
pc <- pt(d, ch4_ppm, d13CH4) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", linewidth = 0.7) +
  geom_hline(yintercept = ATM_SRC, linetype = "dashed", color = RED, linewidth = 0.6) +
  scale_x_log10(labels = scales::label_comma()) + y_num +
  res_lab(sprintf("Plateau %s\u2030 (two-endmember mixing)", neg(ATM_SRC)), right = TRUE) +
  labs(x = expression(CH[4]~(ppm)), y = expression(delta^{13}*CH[4]~"(\u2030)")) + th
# ---- (b) CO2 Keeling ---------------------------------------------------------
pd <- pt(wt, 1/co2_ppm, d13CO2) +
  geom_smooth(method = "lm", formula = y~x, color = "black", fill = "grey80", linewidth = 0.7) +
  annotate("point", x = 0, y = s_ko, color = BLU, size = 2.8) +
  res_lab(sprintf("Intercept %s\u2030 (95%% CI %s to %s)", neg(s_ko), neg(ci_ko[1]), neg(ci_ko[2]))) +
  y_num + labs(x = expression(1/CO[2]~(ppm^{-1})), y = expression(delta^{13}*CO[2]~"(\u2030)")) + th
# ---- (d) eps_C convergence ---------------------------------------------------
pf <- pt(wt, ch4_ppm, eps_C) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = EPS_LO, ymax = EPS_HI, fill = RED, alpha = 0.10) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", linewidth = 0.7) +
  scale_x_log10(labels = scales::label_comma()) + y_num +
  res_lab(if (round(EPS_LO) == round(EPS_HI)) sprintf("'Source-based' ~ epsilon[C] ~ '%s\u2030'", neg(EPS_LO))
          else sprintf("'Source-based' ~ epsilon[C] ~ '%s\u2013%s\u2030'", neg(EPS_LO), neg(EPS_HI)),
          right = TRUE, parse = TRUE) +
  labs(x = expression(CH[4]~(ppm)), y = expression(epsilon[C]~"(\u2030)")) + th

# tags as si_tags() intends (si_tags() itself returns NULL: plot_annotation & theme)
fig <- (pa | pd) / (pc | pf) +          # top row = Keeling source plots; bottom = convergence
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = SI_BASE + 1, face = "bold"))
si_save(out_path("SI_isotopes_source_composite.png"), fig, height = 6.6)
cat(sprintf("CH4 Keeling %.0f | CO2 Keeling %.0f | source-based eps_C = %.0f | n=%d\n",
            s_kc, s_ko, s_ko - s_kc, nrow(wt)))
cat("Wrote SI_isotopes_source_composite.png\n")
