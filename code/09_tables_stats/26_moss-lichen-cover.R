#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 26_moss-lichen-cover.R -- moss and lichen cover at the chamber vs stem CH4 flux
# (Methods S2, "Moss and lichen cover"; Discussion, bark paragraph)
# ------------------------------------------------------------------------------
# The 2023 cross-species datasheet scored moss/lichen cover within each chamber
# footprint visually on a 1-3 scale of relative cover (1 = least); taxa and
# functional groups were not identified. Scores join to fluxes by tag and date.
# Tests: Spearman and Kruskal-Wallis across species; within species, a linear
# model of asinh(flux x 100) on species plus score; Spearman within each species
# with at least 10 scored trees.
# Output: outputs/audit/moss_lichen_flux.txt
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))

f <- read_csv("outputs/data/flux_measurements_tree.csv", show_col_types = FALSE) %>%
  filter(campaign == "2023_cross_species", !dead_stem) %>%
  mutate(tag = as.character(tree_id), day = as.Date(Date))
r <- read_csv("data/raw/flux_windows/lgr_manual_identification_2023_cross_species.csv",
              show_col_types = FALSE) %>%
  distinct(UniqueID, .keep_all = TRUE) %>%
  transmute(tag = as.character(`Tree Tag`), day = as.Date(datetime_start), moss = `Moss/Lichen Cover`) %>%
  distinct(tag, day, .keep_all = TRUE)

j <- f %>% left_join(r, by = c("tag", "day")) %>%
  filter(moss %in% c("1", "2", "3")) %>%
  mutate(moss = as.integer(moss), y = asinh(stem_flux_nmol_m2_s * 100))

out <- c(sprintf("Live-stem 2023 fluxes: %d; with a moss/lichen score: %d", nrow(f), nrow(j)), "")
tab <- j %>% group_by(moss) %>%
  summarise(n = n(), median = median(stem_flux_nmol_m2_s), mean = mean(stem_flux_nmol_m2_s),
            pct_positive = 100 * mean(stem_flux_nmol_m2_s > 0))
out <- c(out, capture.output(print(tab)), "")

sp <- cor.test(j$moss, j$stem_flux_nmol_m2_s, method = "spearman", exact = FALSE)
kw <- kruskal.test(stem_flux_nmol_m2_s ~ factor(moss), j)
out <- c(out, sprintf("Across species: Spearman rho = %.2f, p = %.3f; Kruskal-Wallis p = %.3f",
                      sp$estimate, sp$p.value, kw$p.value))

m0 <- lm(y ~ species_clean, j); m1 <- lm(y ~ species_clean + moss, j)
cf <- summary(m1)$coefficients["moss", ]
out <- c(out, sprintf("Within species (asinh(flux x 100) ~ species + score): slope = %.2f (SE %.2f), p = %.2f",
                      cf[1], cf[2], anova(m0, m1)$`Pr(>F)`[2]), "")

bysp <- j %>% group_by(species_clean) %>% filter(n() >= 10, n_distinct(moss) > 1) %>%
  summarise(n = n(), rho = cor(moss, stem_flux_nmol_m2_s, method = "spearman"),
            p = cor.test(moss, stem_flux_nmol_m2_s, method = "spearman", exact = FALSE)$p.value,
            share_score1 = mean(moss == 1), share_score3 = mean(moss == 3)) %>% arrange(p)
out <- c(out, "By species (at least 10 scored trees):", capture.output(print(bysp)))

writeLines(out, out_path("moss_lichen_flux.txt"))
cat(out, sep = "\n")
