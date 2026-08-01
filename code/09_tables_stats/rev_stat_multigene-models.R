source("code/lib/outputs.R")
# ==============================================================================
# REVISION — species-level multi-gene flux models (does adding a methanotroph term
# to mcrA improve the fit?). Companion to R_S1S2_arcsinh.R / R_figS11_final.R.
# Species level (n=10, >=5 trees & >=5 flux). Response = arcsinh flux (cofactor 0.1),
# transform-then-aggregate on both axes. Reports every combination with R2, ADJUSTED
# R2, AIC (the fair small-n metrics), and the per-term coefficient sign + p, for BOTH
# median and mean aggregation.
# Key question (Jon): flux ~ mcrA + pmoA vs flux ~ mcrA. Honest reading: adding the
# oxidation term improves adj-R2 under MEDIAN (mcrA enters +, methanotroph -), but not
# under MEAN, and no methanotroph term is individually significant at n=10.
# Sources 02_scale for data only (outputs git-restored by wrapper). NEW file.
# Output: outputs/revision/species_multigene_models.txt
# ==============================================================================
suppressWarnings(suppressMessages(source("code/07_molecular/02_scale_dependent_gene_patterns.R")))
suppressPackageStartupMessages(library(tidyverse))
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(out_path("species_multigene_models.txt"), split = TRUE)
asf <- function(x) asinh(x / 0.1) / log(10)

# species aggregates: transform THEN aggregate (median & mean), arcsinh flux + log genes
TL <- tree_level_complete %>% mutate(lg_mcra = log10(mcrA+1), lg_pmoa = log10(pmoA+1),
        lg_mmox = log10(mmoX+1), lg_meth = log10(methanotroph_total+1), lg_ratio = log_ratio)
G <- TL %>% group_by(species, species_id) %>%
  summarise(across(c(lg_mcra, lg_pmoa, lg_mmox, lg_meth, lg_ratio),
                   list(med = ~median(.x, na.rm=TRUE), mn = ~mean(.x, na.rm=TRUE)), .names = "{.col}_{.fn}"),
            n_trees = n(), .groups = "drop")
F <- flux_all %>% mutate(af = asf(CH4_flux)) %>% group_by(species, species_id) %>%
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
