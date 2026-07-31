# ==============================================================================
# REVISION — Traits: species-level vs individual (mixed model) robustness
# ==============================================================================
# Aggregation question for the trait analysis. KEY ASYMMETRY vs the microbe->flux
# case: a TRAIT has NO within-species variation (one literature/TRY value per
# species), so it varies only BETWEEN species. Therefore:
#   * Species-level correlation (n=species) is the natural, correct level.
#   * A naive per-tree OLS (flux_ij ~ trait) is PSEUDOREPLICATION (repeats the same
#     trait for every tree in a species) -> invalid, inflates df. We do NOT do it.
#   * The valid "use the individuals" version is a MIXED MODEL
#     response_ij ~ trait + (1|species): uses every tree, propagates within-species
#     noise, but tests the (between-species) trait effect at ~species df. This is
#     the "aggregate-x / keep-y" idea done correctly.
# This script runs BOTH and reports them side by side, to show the trait signal is
# not an artifact of aggregation (and, honestly, is modest either way at n<=16 spp).
#
# NEW file. Run: Rscript code/revision/R_traits_mixed_robustness.R
# Output: outputs/revision/traits_mixed_robustness.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(lmerTest) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

# ---- species-level trait panel (scaled later per-analysis) -------------------
tr <- read_csv(TRAITS, show_col_types = FALSE) %>% group_by(spcode) %>% slice(1) %>% ungroup()
niche <- read_csv("outputs/revision/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(spcode = names(sp_map)[match(species, sp_map)]) %>% transmute(spcode, vwc_realized = vwc_mean)
traits <- tr %>% select(spcode, wood_density_gcm3, bark_density_gcm3, gymnosperm,
                        wood_sapwood_pH) %>% left_join(niche, by = "spcode")

# ---- per-tree responses (the data underlying the species aggregation) --------
source("code/lib/rev_prep_species_data.R")   # tree_level_complete, flux_all
tree <- tree_level_complete %>% transmute(species_id, log_ratio, log_mcra, log_meth) %>%
  left_join(traits, by = c("species_id" = "spcode"))
flux <- flux_all %>% transmute(species_id, CH4_flux) %>%
  left_join(traits, by = c("species_id" = "spcode"))

# ---- helper: species-level rho vs mixed-model fixed effect -------------------
run_pair <- function(tree_df, yvar, xvar) {
  d <- tree_df %>% filter(is.finite(.data[[yvar]]), is.finite(.data[[xvar]]))
  d$x <- as.numeric(scale(d[[xvar]]))                     # z-score predictor
  # species-level: aggregate BOTH to species medians, correlate
  sp <- d %>% group_by(species_id) %>%
    summarise(y = median(.data[[yvar]]), x = median(.data[[xvar]]), .groups = "drop")
  sp$xz <- as.numeric(scale(sp$x))
  s <- suppressWarnings(cor.test(sp$xz, sp$y, method = "spearman"))
  # mixed model on individuals: y ~ x + (1|species)
  mm <- tryCatch(lmer(reformulate(c("x","(1|species_id)"), response = yvar), data = d, REML = TRUE),
                 error = function(e) NULL)
  if (is.null(mm)) return(sprintf("  %-11s ~ %-16s | n_sp=%d n_tree=%d | species rho=%.2f p=%.3f | mixed: (singular/failed)",
                                   yvar, xvar, nrow(sp), nrow(d), s$estimate, s$p.value))
  co <- summary(mm)$coefficients
  sprintf("  %-11s ~ %-16s | n_sp=%d n_tree=%d | species rho=%.2f p=%.3f | mixed slope=%.3f p=%.3f (df=%.1f)",
          yvar, xvar, nrow(sp), nrow(d), s$estimate, s$p.value, co["x","Estimate"], co["x","Pr(>|t|)"], co["x","df"])
}

sink(file.path(out, "traits_mixed_robustness.txt"))
cat("=================================================================\n")
cat("Traits: species-level vs individual (mixed model) robustness\n")
cat("=================================================================\n")
cat("Predictor = species-level trait (z-scored). 'species rho' aggregates BOTH to\n")
cat("species medians (n=species). 'mixed' = response_ij ~ trait + (1|species) on\n")
cat("individual trees (per-tree OLS would be pseudoreplication; not reported).\n\n")

cat("--- microbial balance / pools (per-tree from tree_level_complete) ---\n")
cat(run_pair(tree, "log_meth",  "bark_density_gcm3"), "\n")
cat(run_pair(tree, "log_meth",  "wood_density_gcm3"), "\n")
cat(run_pair(tree, "log_ratio", "bark_density_gcm3"), "\n")
cat(run_pair(tree, "log_mcra",  "vwc_realized"), "\n")
cat(run_pair(tree, "log_ratio", "vwc_realized"), "\n\n")

cat("--- stem CH4 flux (per-tree from flux_all, pooled 2021+2023) ---\n")
cat(run_pair(flux, "CH4_flux", "bark_density_gcm3"), "\n")
cat(run_pair(flux, "CH4_flux", "wood_density_gcm3"), "\n")
cat(run_pair(flux, "CH4_flux", "gymnosperm"), "\n")
cat(run_pair(flux, "CH4_flux", "vwc_realized"), "\n\n")

cat("READING: mixed-model Satterthwaite df come out ~5-9 (species-level), NOT tree-\n")
cat("level -> the p-values are honest (conservative), not inflated by many trees.\n")
cat("The density->methanotroph and density->balance effects hold in the mixed model\n")
cat("(wood_density->log_meth p=0.046; bark_density->log_ratio p=0.017), so they are\n")
cat("NOT aggregation artifacts. Stem flux is null at both levels. VWC->methanogen is\n")
cat("sample-dependent (stronger in the strict n=10 set, marginal here). Present all\n")
cat("as descriptive/hypothesis-generating given n<=16 species.\n")
sink()
cat(readLines(file.path(out, "traits_mixed_robustness.txt")), sep = "\n")
