# ==============================================================================
# 16_s1-rf-soil-arcsinh.R -- tree-level gene random forest and soil-gene additions
# ------------------------------------------------------------------------------
# (1) An area-weighted random forest (species + DBH + area-weighted genes), with the
#     response raw and arcsinh (asinh(flux / 0.1) / log(10)): variance explained and
#     importance.
# (2) Soil genes (mcrA, pmoA, mmoX per g) added to species + tree mmoX, raw and arcsinh.
# Data: code/07_molecular/helper_scale_dependent_gene_patterns.R (tables only).
# Output: outputs/audit/S1_rf_soil_arcsinh_report.txt
# ==============================================================================
source("code/lib/outputs.R")
suppressWarnings(suppressMessages(source("code/07_molecular/helper_scale_dependent_gene_patterns.R")))
suppressPackageStartupMessages({ library(tidyverse); library(randomForest) })
sink(out_path("S1_rf_soil_arcsinh_report.txt"), split = TRUE)
source("code/lib/flux_transform.R")   # flux_asinh10()
line <- function() cat(strrep("-", 78), "\n")

# ---- assemble area-weighted RF data: tree_level_complete genes + dbh ----------
y <- ymf2021   # survey_2021_trees() (data/tables/), sourced with the helper above
rf <- tree_level_complete %>%
  left_join(y %>% dplyr::select(tree_id, dbh), by = "tree_id") %>%
  dplyr::select(CH4_flux, species, dbh, log_tree_mcra, log_pmoa, log_tree_mmox) %>%
  filter(complete.cases(.)) %>% mutate(species = factor(species))

cat("### S1 (RF). Area-weighted random forest: species + dbh + area-weighted genes\n"); line()
cat("n trees =", nrow(rf), "\n")
run_rf <- function(resp_raw) {
  d <- rf; d$y <- if (resp_raw) d$CH4_flux else flux_asinh10(d$CH4_flux); d$CH4_flux <- NULL
  set.seed(42)
  m <- randomForest(y ~ ., data = d, importance = TRUE, ntree = 1000)
  ve <- 100 * m$rsq[length(m$rsq)]
  imp <- importance(m)[, "%IncMSE"]; imp <- sort(imp, decreasing = TRUE)
  list(ve = ve, imp = imp)
}
for (resp in c("RAW", "ARCSINH")) {
  r <- run_rf(resp == "RAW")
  cat(sprintf("\n[%s]  variance explained = %.1f%%\n   importance (%%IncMSE, ranked):\n", resp, r$ve))
  for (nm in names(r$imp)) cat(sprintf("      %-16s %.2f\n", nm, r$imp[nm]))
}

# ---- soil-gene additions (linear): species + tree mmoX + soil gene ------------
cat("\n\n### S1 (soil). Soil-gene additions to species + tree mmoX (linear)\n"); line()
# soil genes at the tree: copies per g (fresh soil; data/tables/ddpcr.csv via survey_2021_trees()),
# mean of the organic and mineral horizons (NA where neither was measured)
soil_mean <- function(g) { o <- y[[paste0("ddpcr_", g, "_Organic_loose")]]; m <- y[[paste0("ddpcr_", g, "_Mineral_loose")]]
  ifelse(is.na(o) & is.na(m), NA_real_, rowMeans(cbind(o, m), na.rm = TRUE)) }
sd <- tibble(tree_id = y$tree_id, soil_mcra_probe = soil_mean("mcra_probe"), soil_pmoa = soil_mean("pmoa"),
             soil_mmox = soil_mean("mmox"))
base <- tree_level_complete %>% left_join(sd, by = "tree_id")
soil_genes <- c(mcrA = "soil_mcra_probe", pmoA = "soil_pmoa", mmoX = "soil_mmox")
for (nm in names(soil_genes)) {
  g <- soil_genes[[nm]]; if (!g %in% names(base)) { cat(sprintf("  soil %s: column missing\n", nm)); next }
  d <- base %>% filter(!is.na(.data[[g]]), !is.na(CH4_flux)) %>% mutate(lsg = log10(.data[[g]] + 1))
  for (resp in c("RAW", "ARCSINH")) {
    d$yv <- if (resp == "RAW") d$CH4_flux else flux_asinh10(d$CH4_flux)
    m0 <- lm(yv ~ species + log_tree_mmox, d); m1 <- lm(yv ~ species + log_tree_mmox + lsg, d)
    dR2 <- summary(m1)$r.squared - summary(m0)$r.squared; p <- summary(m1)$coef["lsg", 4]
    cat(sprintf("  soil %-5s [%-7s] n=%d  dR2=%+.4f  p(soil)=%.3f\n", nm, resp, nrow(d), dR2, p))
  }
}
# ---- concentration RF: species + dbh + tissue-specific + soil genes ----------
cat("\n\n### S1 (RF conc). Concentration random forest: species + dbh + hw/sw genes + soil genes\n"); line()
cat("(tests the 'heartwood mcrA is top predictor' ranking; absolute variance not comparable to paper\n")
cat(" because environmental predictors from the original pipeline are not included)\n")
crf <- y %>% filter(!is.na(CH4_best.flux_125cm), !is.na(species_id)) %>%
  transmute(tree_id, species = factor(species_name(species_id)), dbh, CH4_flux = CH4_best.flux_125cm,
    hw_mcra=log10(ddpcr_mcra_probe_Inner_loose+1), sw_mcra=log10(ddpcr_mcra_probe_Outer_loose+1),
    hw_pmoa=log10(ddpcr_pmoa_Inner_loose+1), sw_pmoa=log10(ddpcr_pmoa_Outer_loose+1),
    hw_mmox=log10(ddpcr_mmox_Inner_loose+1), sw_mmox=log10(ddpcr_mmox_Outer_loose+1)) %>%
  left_join(sd, by = "tree_id") %>%
  mutate(soil_mcra=log10(soil_mcra_probe+1), soil_pmoa=log10(soil_pmoa+1), soil_mmox=log10(soil_mmox+1)) %>%
  dplyr::select(CH4_flux, species, dbh, hw_mcra, sw_mcra, hw_pmoa, sw_pmoa, hw_mmox, sw_mmox, soil_mcra, soil_pmoa, soil_mmox) %>%
  filter(complete.cases(.))
cat("n trees =", nrow(crf), "\n")
for (resp in c("RAW", "ARCSINH")) {
  d <- crf; d$yv <- if (resp == "RAW") d$CH4_flux else flux_asinh10(d$CH4_flux); d$CH4_flux <- NULL
  set.seed(42); m <- randomForest(yv ~ ., data = d, importance = TRUE, ntree = 1000)
  imp <- sort(importance(m)[, "%IncMSE"], decreasing = TRUE)
  cat(sprintf("\n[%s]  variance explained = %.1f%%   top predictors (%%IncMSE):\n", resp, 100*m$rsq[length(m$rsq)]))
  for (nm in names(imp)[1:6]) cat(sprintf("      %-12s %.2f\n", nm, imp[nm]))
}
sink(); cat("Wrote S1_rf_soil_arcsinh_report.txt\n")
