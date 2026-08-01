source("code/lib/outputs.R")
# ==============================================================================
# REVISION — SI extended-isotope source composite (4 panels, 2x2).
# Design (Jon): CH4 Keeling      | CO2 Keeling
#               d13CH4 convergence | eps_C convergence
# (Miller-Tans dropped: leverage-biased, redundant with Keeling.)
# Data load replicated verbatim from ISOTOPES_final.R (settled QC). NEW file.
# Output: outputs/revision/SI_isotopes_source_composite.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
RED <- "#b2182b"; BLU <- "#2166ac"

# ---- data (verbatim from ISOTOPES_final.R) -----------------------------------
runs <- list.files("data/raw/internal_gas/picarro", pattern = "_results.csv$", full.names = TRUE)
runs <- runs[!grepl("merged", runs)]
raw  <- map_dfr(runs, ~read_csv(.x, show_col_types = FALSE))
STANDARDS <- c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c","SA1")
d <- raw %>%
  transmute(SampleName,
            d13CH4 = HR_Delta_iCH4_Raw_mean, ch4_ppm = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean,     co2_ppm = `12CO2_mean`) %>%
  filter(!str_detect(SampleName, "^Amb"), !str_detect(SampleName, "^V"),
         !str_detect(SampleName, "[HS]$"), !SampleName %in% STANDARDS,
         !is.na(d13CH4), ch4_ppm >= 1.5)
wt <- d %>% filter(!is.na(d13CO2), co2_ppm > 0) %>%
  mutate(eps_C = ((d13CO2 + 1000)/(d13CH4 + 1000) - 1)*1000)

th <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 8), panel.grid.minor = element_blank())
pt <- function(dat, x, y) ggplot(dat, aes({{x}}, {{y}})) +
  geom_point(alpha = 0.45, size = 1.2, color = "grey30")

# ---- source estimates --------------------------------------------------------
kc <- lm(d13CH4 ~ I(1/ch4_ppm), data = d);  s_kc <- coef(kc)[1];  ci_kc <- confint(kc)[1,]
ko <- lm(d13CO2 ~ I(1/co2_ppm), data = wt); s_ko <- coef(ko)[1];  ci_ko <- confint(ko)[1,]

# ---- (a) CH4 Keeling ---------------------------------------------------------
pa <- pt(d, 1/ch4_ppm, d13CH4) +
  geom_smooth(method = "lm", formula = y~x, color = "black", fill = "grey80", linewidth = 0.7) +
  annotate("point", x = 0, y = s_kc, color = RED, size = 2.8) +
  labs(x = expression(1/CH[4]~"(ppm"^-1*")"), y = expression(delta^13*"C-CH"[4]~"(per mil)"),
       title = "CH4 Keeling",
       subtitle = sprintf("source intercept = %.0f [%.0f, %.0f] per mil", s_kc, ci_kc[1], ci_kc[2])) + th
# ---- d13CH4 convergence ------------------------------------------------------
pc <- pt(d, ch4_ppm, d13CH4) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", linewidth = 0.7) +
  geom_hline(yintercept = -70, linetype = "dashed", color = RED, linewidth = 0.6) +
  scale_x_log10() +
  labs(x = expression(CH[4]~"(ppm, log)"), y = expression(delta^13*"C-CH"[4]~"(per mil)"),
       title = "d13CH4 source convergence", subtitle = "plateau ~ -70 per mil (atmosphere-corrected)") + th
# ---- (d) CO2 Keeling ---------------------------------------------------------
pd <- pt(wt, 1/co2_ppm, d13CO2) +
  geom_smooth(method = "lm", formula = y~x, color = "black", fill = "grey80", linewidth = 0.7) +
  annotate("point", x = 0, y = s_ko, color = BLU, size = 2.8) +
  labs(x = expression(1/CO[2]~"(ppm"^-1*")"), y = expression(delta^13*"C-CO"[2]~"(per mil)"),
       title = "CO2 Keeling",
       subtitle = sprintf("source intercept = %.0f [%.0f, %.0f] per mil", s_ko, ci_ko[1], ci_ko[2])) + th
# ---- eps_C convergence -------------------------------------------------------
pf <- pt(wt, ch4_ppm, eps_C) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = 54, ymax = 64, fill = RED, alpha = 0.10) +
  geom_smooth(method = "loess", se = TRUE, color = "black", fill = "grey80", linewidth = 0.7) +
  scale_x_log10() +
  labs(x = expression(CH[4]~"(ppm, log)"), y = expression(epsilon[C]~"(per mil)"),
       title = "Apparent fractionation convergence",
       subtitle = expression(paste("plateau ", epsilon[C], " ~ 54-64 per mil (hydrogenotrophic)"))) + th

fig <- (pa | pd) / (pc | pf) +          # top row = Keeling source plots; bottom = convergence
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(out_path("SI_isotopes_source_composite.png"), fig, width = 9.5, height = 8, dpi = 300, bg = "white")
cat(sprintf("CH4 Keeling %.0f | CO2 Keeling %.0f | source-based eps_C = %.0f | n=%d\n",
            s_kc, s_ko, s_ko - s_kc, nrow(wt)))
cat("Wrote SI_isotopes_source_composite.png\n")
