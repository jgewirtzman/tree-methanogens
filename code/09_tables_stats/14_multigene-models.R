# ==============================================================================
# 14_multigene-models.R -- species-level multi-gene flux models
# ------------------------------------------------------------------------------
# Does adding a methanotroph term to mcrA improve the species-level fit? Species
# with >= 5 trees and >= 5 fluxes (n = 10). Response: arcsinh flux (cofactor 0.1,
# log10 units), transformed before aggregating; genes log10(x + 1). Every model is
# reported with R2, adjusted R2, AIC and each term's sign and p, for median and mean
# aggregation.
# Data: code/07_molecular/helper_scale_dependent_gene_patterns.R (tables only).
# Output: outputs/audit/species_multigene_models.txt
# ==============================================================================
source("code/lib/outputs.R")
suppressWarnings(suppressMessages(source("code/07_molecular/helper_scale_dependent_gene_patterns.R")))
suppressPackageStartupMessages(library(tidyverse))
sink(out_path("species_multigene_models.txt"), split = TRUE)
source("code/lib/flux_transform.R")   # flux_asinh10()

# species aggregates: transform THEN aggregate (median & mean), arcsinh flux + log genes
TL <- tree_level_complete %>% mutate(lg_mcra = log10(mcrA+1), lg_pmoa = log10(pmoA+1),
        lg_mmox = log10(mmoX+1), lg_meth = log10(methanotroph_total+1), lg_ratio = log_ratio)
G <- TL %>% group_by(species, species_id) %>%
  summarise(across(c(lg_mcra, lg_pmoa, lg_mmox, lg_meth, lg_ratio),
                   list(med = ~median(.x, na.rm=TRUE), mn = ~mean(.x, na.rm=TRUE)), .names = "{.col}_{.fn}"),
            n_trees = n(), .groups = "drop")
F <- flux_all %>% mutate(af = flux_asinh10(CH4_flux)) %>% group_by(species, species_id) %>%
  summarise(fmed = median(af, na.rm=TRUE), fmn = mean(af, na.rm=TRUE), n_flux = n(), .groups = "drop")
S <- G %>% inner_join(F, by = c("species","species_id")) %>% filter(n_trees >= 5, n_flux >= 5)
cat("n species =", nrow(S), "  (response = arcsinh flux, cofactor 0.1)\n\n")

# models: name -> predictor columns (suffix _med / _mn filled per aggregation)
MODELS <- list(
  "mcrA"                = c("lg_mcra"),
  "pmoA"                = c("lg_pmoa"),
  "mmoX"                = c("lg_mmox"),
  "pmoA+mmoX (methanotroph)" = c("lg_meth"),
  "ratio (mcrA:MOB)"    = c("lg_ratio"),
  "mcrA + pmoA"         = c("lg_mcra","lg_pmoa"),
  "mcrA + mmoX"         = c("lg_mcra","lg_mmox"),
  "mcrA + methanotroph" = c("lg_mcra","lg_meth"),
  "mcrA + pmoA + mmoX"  = c("lg_mcra","lg_pmoa","lg_mmox"))

fit_report <- function(sfx, yv) {
  cat(sprintf("=== %s aggregation ===\n", ifelse(sfx=="med","SPECIES-MEDIAN","SPECIES-MEAN")))
  cat(sprintf("%-26s %6s %7s %7s   %s\n", "model", "R2", "adjR2", "AIC", "term coefficients (sign, p)"))
  cat(strrep("-", 92), "\n")
  for (nm in names(MODELS)) {
    xs <- paste0(MODELS[[nm]], "_", sfx)
    d <- S %>% transmute(y = .data[[yv]], across(all_of(xs))) %>% filter(complete.cases(.))
    m <- lm(as.formula(paste("y ~", paste(xs, collapse=" + "))), d); sm <- summary(m)
    terms <- rownames(sm$coefficients)[-1]
    tl <- paste(sapply(terms, function(t) sprintf("%s=%+.2f(p=%.3f)", sub(paste0("_",sfx),"",t),
                       sm$coefficients[t,1], sm$coefficients[t,4])), collapse="  ")
    cat(sprintf("%-26s %6.3f %7.3f %7.1f   %s\n", nm, sm$r.squared, sm$adj.r.squared, AIC(m), tl))
  }
  cat("\n")
}
fit_report("med", "fmed")
fit_report("mn",  "fmn")

cat("READING:\n")
cat("- Compare on ADJUSTED R2 / AIC (plain R2 always rises with more predictors at n=10).\n")
cat("- mcrA enters POSITIVE (production); pmoA/mmoX/methanotroph enter NEGATIVE (oxidation).\n")
cat("- Adding a methanotroph term to mcrA improves adj-R2 under MEDIAN but not MEAN; no\n")
cat("  methanotroph term is individually significant (n=10) -> suggestive, not decisive.\n")
sink()
cat("Wrote species_multigene_models.txt\n")
