source("code/lib/outputs.R")
# ==============================================================================
# REVISION — Robustness of the central gene-flux results to splitting pmoA/mmoX
# ==============================================================================
# R2 #1c asked pmoA and mmoX be treated separately (many methanotrophs encode both).
# The central Fig-8 inferences use the COMBINED methanotroph total (pmoA+mmoX):
#   (i)  median_flux ~ log10(gene)         [abundance -> flux]
#   (ii) median_flux ~ log10(mcrA:methanotroph ratio)   [the strong R2~0.51 result]
# Here we re-run BOTH with pmoA-only and mmoX-only in place of the combined total,
# to check whether any major conclusion changes.
#
# Reuses the exact per-species objects from 02_scale_dependent_gene_patterns.R
# (analysis_mcra/pmoa/mmox/methanotroph/ratio; 10 species). NEW file; edits nothing.
# Output: outputs/revision/pmoa_mmox_robustness.txt
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

invisible(capture.output(suppressMessages(
  source("code/07_molecular/02_scale_dependent_gene_patterns.R"))))

# assemble one per-species table (median gene abundances + median flux)
d <- flux_by_species %>% dplyr::select(species, median_flux) %>%
  left_join(analysis_mcra         %>% dplyr::select(species, median_mcra),         by = "species") %>%
  left_join(analysis_pmoa         %>% dplyr::select(species, median_pmoa),         by = "species") %>%
  left_join(analysis_mmox         %>% dplyr::select(species, median_mmox),         by = "species") %>%
  left_join(analysis_methanotroph %>% dplyr::select(species, median_methanotroph), by = "species") %>%
  filter(!is.na(median_mcra), !is.na(median_pmoa), !is.na(median_mmox), !is.na(median_flux))

fit_g <- function(g) { m <- lm(d$median_flux ~ log10(g + 1)); s <- summary(m)
  c(R2 = s$r.squared, p = coef(s)[2, 4]) }
fit_r <- function(num, den) { r <- log10((num + 1) / (den + 1)); m <- lm(d$median_flux ~ r); s <- summary(m)
  c(R2 = s$r.squared, p = coef(s)[2, 4]) }

rows <- tribble(~model, ~term, ~stats,
  "Abundance -> flux", "mcrA",                     fit_g(d$median_mcra),
  "Abundance -> flux", "pmoA+mmoX (combined)",     fit_g(d$median_methanotroph),
  "Abundance -> flux", "pmoA only",                fit_g(d$median_pmoa),
  "Abundance -> flux", "mmoX only",                fit_g(d$median_mmox),
  "Ratio -> flux",     "mcrA : (pmoA+mmoX)",       fit_r(d$median_mcra, d$median_methanotroph),
  "Ratio -> flux",     "mcrA : pmoA",              fit_r(d$median_mcra, d$median_pmoa),
  "Ratio -> flux",     "mcrA : mmoX",              fit_r(d$median_mcra, d$median_mmox)) %>%
  mutate(R2 = round(map_dbl(stats, "R2"), 3), p = round(map_dbl(stats, "p"), 3)) %>%
  dplyr::select(model, term, R2, p)

sink(out_path("pmoa_mmox_robustness.txt"))
cat("Robustness of central gene-flux results to pmoA/mmoX split (n =", nrow(d), "species)\n")
cat("=====================================================================\n")
print(as.data.frame(rows), row.names = FALSE)
cat("\nReading:\n")
cat(" * Combined pmoA+mmoX abundance -> flux R2 =", rows$R2[rows$term=="pmoA+mmoX (combined)"], "\n")
cat("   pmoA-only R2 =", rows$R2[rows$term=="pmoA only"],
    "| mmoX-only R2 =", rows$R2[rows$term=="mmoX only"], "\n")
cat(" * KEY ratio mcrA:(pmoA+mmoX) -> flux R2 =", rows$R2[rows$term=="mcrA : (pmoA+mmoX)"],
    "(p =", rows$p[rows$term=="mcrA : (pmoA+mmoX)"], ")\n")
cat("   using mcrA:pmoA  R2 =", rows$R2[rows$term=="mcrA : pmoA"],
    "(p =", rows$p[rows$term=="mcrA : pmoA"], ")\n")
cat("   using mcrA:mmoX  R2 =", rows$R2[rows$term=="mcrA : mmoX"],
    "(p =", rows$p[rows$term=="mcrA : mmoX"], ")\n")
sink()
cat(readLines(out_path("pmoa_mmox_robustness.txt")), sep = "\n")
