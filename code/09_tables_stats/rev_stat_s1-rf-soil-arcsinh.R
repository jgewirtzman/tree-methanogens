source("code/lib/outputs.R")
# ==============================================================================
# REVISION — S1 area-weighted RANDOM FOREST + SOIL-gene additions, redone with an
# arcsinh RESPONSE (asinh(flux/0.1)/log10). Companion to R_S1S2_arcsinh.R.
# Area-weighted RF (species + dbh + area-weighted genes) is reliably assemblable from
# tree_level_complete + dbh; reports RAW vs ARCSINH variance-explained + importance.
# Soil-gene additions: soil mcrA/pmoA/mmoX added to species + tree mmoX (RAW vs ARCSINH).
# NOTE: the CONCENTRATION RF and the full partial-dependence / H-statistic diagnostics
# live in a 5000-line path-dependent pipeline (methanotrophs/01_rf_models.R) and are
# NOT reconstructed here; the faithful way to redo those is a copy-and-modify of that
# script swapping the response to asf(flux). RF variable-importance rankings are
# permutation-based and expected stable to a monotonic response transform.
# NEW file; 02_scale sourced for data only (2 outputs git-restored by wrapper).
# Output: outputs/revision/S1_rf_soil_arcsinh_report.txt
# ==============================================================================
suppressWarnings(suppressMessages(source("code/07_molecular/02_scale_dependent_gene_patterns.R")))
suppressPackageStartupMessages({ library(tidyverse); library(randomForest) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(out_path("S1_rf_soil_arcsinh_report.txt"), split = TRUE)
asf <- function(x) asinh(x / 0.1) / log(10)
line <- function() cat(strrep("-", 78), "\n")

# ---- assemble area-weighted RF data: tree_level_complete genes + dbh ----------
y <- read.csv("data/processed/integrated/merged_tree_dataset_final.csv")
rf <- tree_level_complete %>%
  left_join(y %>% dplyr::select(tree_id, dbh), by = "tree_id") %>%
  dplyr::select(CH4_flux, species, dbh, log_tree_mcra, log_pmoa, log_tree_mmox) %>%
  filter(complete.cases(.)) %>% mutate(species = factor(species))

cat("### S1 (RF). Area-weighted random forest: species + dbh + area-weighted genes\n"); line()
cat("n trees =", nrow(rf), "\n")
run_rf <- function(resp_raw) {
  d <- rf; d$y <- if (resp_raw) d$CH4_flux else asf(d$CH4_flux); d$CH4_flux <- NULL
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
sd <- read.csv("data/compiled/ddpcr_gene_abundances.csv") %>%
  filter(analysis_type == "loose", material == "Soil", target_gene %in% c("mcra_probe", "pmoa", "mmox")) %>%
  mutate(tree = tolower(gsub("[^A-Za-z0-9]", "", sub("\n.*", "", sub("^[A-Z]+_", "", sample_id))))) %>%
  group_by(tree, target_gene) %>% summarise(conc = mean(concentration_copies_per_uL, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = target_gene, values_from = conc, names_prefix = "soil_")
base <- tree_level_complete %>% mutate(tree = tolower(gsub("[^A-Za-z0-9]", "", tree_id))) %>%
  left_join(sd, by = "tree")
soil_genes <- c(mcrA = "soil_mcra_probe", pmoA = "soil_pmoa", mmoX = "soil_mmox")
for (nm in names(soil_genes)) {
  g <- soil_genes[[nm]]; if (!g %in% names(base)) { cat(sprintf("  soil %s: column missing\n", nm)); next }
  d <- base %>% filter(!is.na(.data[[g]]), !is.na(CH4_flux)) %>% mutate(lsg = log10(.data[[g]] + 1))
  for (resp in c("RAW", "ARCSINH")) {
    d$yv <- if (resp == "RAW") d$CH4_flux else asf(d$CH4_flux)
    m0 <- lm(yv ~ species + log_tree_mmox, d); m1 <- lm(yv ~ species + log_tree_mmox + lsg, d)
    dR2 <- summary(m1)$r.squared - summary(m0)$r.squared; p <- summary(m1)$coef["lsg", 4]
    cat(sprintf("  soil %-5s [%-7s] n=%d  dR2=%+.4f  p(soil)=%.3f\n", nm, resp, nrow(d), dR2, p))
  }
}
# ---- concentration RF: species + dbh + tissue-specific + soil genes ----------
cat("\n\n### S1 (RF conc). Concentration random forest: species + dbh + hw/sw genes + soil genes\n"); line()
cat("(tests the 'heartwood mcrA is top predictor' ranking; absolute variance not comparable to paper\n")
cat(" because environmental predictors from the original pipeline are not included)\n")
species_mapping <- c("ACRU"="Acer rubrum","ACSA"="Acer saccharum","BEAL"="Betula alleghaniensis","BELE"="Betula lenta",
  "BEPA"="Betula papyrifera","FAGR"="Fagus grandifolia","FRAM"="Fraxinus americana","PIST"="Pinus strobus",
  "QURU"="Quercus rubra","TSCA"="Tsuga canadensis","CAOV"="Carya ovata","KALA"="Kalmia latifolia",
  "PRSE"="Prunus serotina","QUAL"="Quercus alba","QUVE"="Quercus velutina","SAAL"="Sassafras albidum")
crf <- y %>% filter(!is.na(CH4_best.flux_125cm), !is.na(species_id)) %>%
  transmute(tree_id, species = factor(species_mapping[species_id]), dbh, CH4_flux = CH4_best.flux_125cm,
    hw_mcra=log10(ddpcr_mcra_probe_Inner_loose+1), sw_mcra=log10(ddpcr_mcra_probe_Outer_loose+1),
    hw_pmoa=log10(ddpcr_pmoa_Inner_loose+1), sw_pmoa=log10(ddpcr_pmoa_Outer_loose+1),
    hw_mmox=log10(ddpcr_mmox_Inner_loose+1), sw_mmox=log10(ddpcr_mmox_Outer_loose+1)) %>%
  mutate(tree = tolower(gsub("[^A-Za-z0-9]", "", tree_id))) %>% left_join(sd, by = "tree") %>%
  mutate(soil_mcra=log10(soil_mcra_probe+1), soil_pmoa=log10(soil_pmoa+1), soil_mmox=log10(soil_mmox+1)) %>%
  dplyr::select(CH4_flux, species, dbh, hw_mcra, sw_mcra, hw_pmoa, sw_pmoa, hw_mmox, sw_mmox, soil_mcra, soil_pmoa, soil_mmox) %>%
  filter(complete.cases(.))
cat("n trees =", nrow(crf), "\n")
for (resp in c("RAW", "ARCSINH")) {
  d <- crf; d$yv <- if (resp == "RAW") d$CH4_flux else asf(d$CH4_flux); d$CH4_flux <- NULL
  set.seed(42); m <- randomForest(yv ~ ., data = d, importance = TRUE, ntree = 1000)
  imp <- sort(importance(m)[, "%IncMSE"], decreasing = TRUE)
  cat(sprintf("\n[%s]  variance explained = %.1f%%   top predictors (%%IncMSE):\n", resp, 100*m$rsq[length(m$rsq)]))
  for (nm in names(imp)[1:6]) cat(sprintf("      %-12s %.2f\n", nm, imp[nm]))
}
sink(); cat("Wrote S1_rf_soil_arcsinh_report.txt\n")
