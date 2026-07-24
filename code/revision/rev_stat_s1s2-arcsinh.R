# ==============================================================================
# REVISION — S1/S2 tree- & species-level flux prediction, REDONE with an
# arcsinh-transformed RESPONSE (asinh(flux/0.1)/log10; cofactor matches Fig 1/2/3).
# The original workflow (manuscript_statistics.R, 01_rf_models.R) fit RAW CH4_flux;
# flux is right-skewed, so the response should be variance-stabilized. Predictors
# stay log-transformed as in the original. Reports RAW vs ARCSINH side by side.
# Covers: S1 concentration 6-gene, S1 area-weighted + AIC ladder, S2 concentration
# heartwood-mcrA, S2 area-weighted mcrA/methanotroph/ratio + AIC.
# DEFERRED (separate pass): soil-gene additions; random forests (variance-reduction
# splits are also not transform-invariant — redo, but rankings expected stable).
# NEW file; originals untouched (02_scale sourced for data only, its 2 outputs
# git-restored by the wrapper). Output: outputs/revision/S1S2_arcsinh_report.txt
# ==============================================================================
suppressWarnings(suppressMessages(source("code/05_gene_flux_analysis/02_scale_dependent_gene_patterns.R")))
suppressPackageStartupMessages({ library(tidyverse); has_car <- requireNamespace("car", quietly = TRUE) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "S1S2_arcsinh_report.txt"), split = TRUE)
asf <- function(x) asinh(x / 0.1) / log(10)
species_mapping <- c("ACRU"="Acer rubrum","ACSA"="Acer saccharum","BEAL"="Betula alleghaniensis","BELE"="Betula lenta",
  "BEPA"="Betula papyrifera","FAGR"="Fagus grandifolia","FRAM"="Fraxinus americana","PIST"="Pinus strobus",
  "QURU"="Quercus rubra","TSCA"="Tsuga canadensis","CAOV"="Carya ovata","KALA"="Kalmia latifolia",
  "PRSE"="Prunus serotina","QUAL"="Quercus alba","QUVE"="Quercus velutina","SAAL"="Sassafras albidum")
cr <- function(m, term) { s<-summary(m); c(b=unname(coef(m)[term]), se=s$coef[term,2], p=s$coef[term,4], r2=s$r.squared, adjr2=s$adj.r.squared, aic=AIC(m)) }
# species-level flux from the COMBINED 2021+2023 dataset (flux_all, sourced), raw + arcsinh
flux_sp <- flux_all %>% group_by(species, species_id) %>%
  summarise(n_flux=n(), fr=median(CH4_flux), fa=median(asf(CH4_flux)), .groups="drop") %>% filter(n_flux>=5)
line <- function() cat(strrep("-", 78), "\n")

# ============================ S1a. CONCENTRATION 6-GENE ============================
cat("\n### S1a. Concentration tissue-specific 6-gene model (species + hw/sw x 3 genes)\n"); line()
y <- read.csv("data/processed/integrated/merged_tree_dataset_final.csv")
conc <- y %>% filter(!is.na(CH4_best.flux_125cm), !is.na(species_id)) %>%
  mutate(species = species_mapping[species_id],
    log_hw_mcra=log10(ddpcr_mcra_probe_Inner_loose+1), log_sw_mcra=log10(ddpcr_mcra_probe_Outer_loose+1),
    log_hw_pmoa=log10(ddpcr_pmoa_Inner_loose+1),       log_sw_pmoa=log10(ddpcr_pmoa_Outer_loose+1),
    log_hw_mmox=log10(ddpcr_mmox_Inner_loose+1),       log_sw_mmox=log10(ddpcr_mmox_Outer_loose+1)) %>%
  filter(!is.na(species)) %>% filter(complete.cases(dplyr::select(., starts_with("log_hw_"), starts_with("log_sw_"))))
cat("n trees =", nrow(conc), " | n species =", n_distinct(conc$species), "\n")
rhs <- "species + log_hw_mcra + log_sw_mcra + log_hw_pmoa + log_sw_pmoa + log_hw_mmox + log_sw_mmox"
for (resp in c("RAW","ARCSINH")) {
  conc$yv <- if (resp=="RAW") conc$CH4_best.flux_125cm else asf(conc$CH4_best.flux_125cm)
  m <- lm(as.formula(paste("yv ~", rhs)), conc); s <- summary(m)
  cat(sprintf("\n[%s]  model R2=%.3f  adjR2=%.3f  AIC=%.1f\n", resp, s$r.squared, s$adj.r.squared, AIC(m)))
  for (g in c("log_hw_mcra","log_sw_mcra","log_hw_pmoa","log_sw_pmoa","log_hw_mmox","log_sw_mmox"))
    cat(sprintf("   %-12s b=%+.3f  SE=%.3f  p=%.4f\n", g, coef(m)[g], s$coef[g,2], s$coef[g,4]))
  if (has_car) { a<-car::Anova(m, type=2); cat("   ANOVA-II sw_mmox: F=%.3f p=%.4f | species: F=%.3f p=%.2e\n" %>%
    sprintf(a["log_sw_mmox","F value"], a["log_sw_mmox","Pr(>F)"], a["species","F value"], a["species","Pr(>F)"])) }
}

# ============================ S1b/c. AREA-WEIGHTED ============================
cat("\n\n### S1b. Area-weighted single-gene models (species + gene)\n"); line()
TL <- tree_level_complete %>% mutate(af = asf(CH4_flux))
awg <- c(mcrA="log_tree_mcra", pmoA="log_pmoa", mmoX="log_tree_mmox", methanotroph="log_methanotroph")
cat("n trees =", nrow(TL), "\n")
for (nm in names(awg)) { g<-awg[[nm]]; TL$g<-TL[[g]]; d<-TL %>% filter(is.finite(g))
  mr<-cr(lm(CH4_flux~species+g,d),"g"); ma<-cr(lm(af~species+g,d),"g")
  cat(sprintf("%-13s RAW: b=%+.3f p=%.3f R2=%.3f | ARCSINH: b=%+.3f p=%.3f R2=%.3f\n", nm, mr["b"],mr["p"],mr["r2"], ma["b"],ma["p"],ma["r2"])) }
cat("\n### S1c. Area-weighted AIC ladder (ARCSINH response)\n"); line()
mods <- list("species"="species", "species+mmoX"="species+log_tree_mmox",
  "species+mmoX+mcrA"="species+log_tree_mmox+log_tree_mcra", "species+mmoX+pmoA"="species+log_tree_mmox+log_pmoa",
  "species+mmoX+mcrA+pmoA"="species+log_tree_mmox+log_tree_mcra+log_pmoa")
for (nm in names(mods)) { m<-lm(as.formula(paste("af ~", mods[[nm]])), TL); s<-summary(m)
  cat(sprintf("%-26s R2=%.3f  adjR2=%.3f  AIC=%.1f\n", nm, s$r.squared, s$adj.r.squared, AIC(m))) }

# ============================ S2. SPECIES-LEVEL ============================
cat("\n\n### S2a. Concentration species-level: heartwood mcrA (median hw_mcra, n_trees>=5; flux from flux_all)\n"); line()
sp_hw <- y %>% filter(!is.na(species_id), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  group_by(species_id) %>% summarise(n_trees=n(), hw=log10(median(ddpcr_mcra_probe_Inner_loose)+1), .groups="drop") %>%
  filter(n_trees>=5) %>% inner_join(flux_sp, by="species_id")
cat("n species =", nrow(sp_hw), " (paper: 10, R2=0.611)\n")
for (resp in c("RAW","ARCSINH")) { yv<-if(resp=="RAW") sp_hw$fr else sp_hw$fa; ct<-cor.test(sp_hw$hw, yv)
  cat(sprintf("[%-7s]  hw_mcrA: r=%.3f  R2=%.3f  p=%.4f\n", resp, ct$estimate, ct$estimate^2, ct$p.value)) }

cat("\n### S2b. Area-weighted species-level: mcrA / methanotroph / ratio (median gene, n_trees>=5; flux from flux_all)\n"); line()
gene_sp <- tree_level_complete %>% group_by(species, species_id) %>%
  summarise(n_trees=n(), mcra=log10(median(mcrA)+1), meth=log10(median(methanotroph_total)+1), ratio=median(log_ratio), .groups="drop") %>%
  filter(n_trees>=5) %>% inner_join(flux_sp, by=c("species","species_id"))
cat("n species =", nrow(gene_sp), " (paper raw: mcrA R2=0.394, ratio R2=0.513)\n")
for (resp in c("RAW","ARCSINH")) { yv<-if(resp=="RAW") gene_sp$fr else gene_sp$fa
  cat(sprintf("\n[%-7s]\n", resp))
  for (v in c("mcra","meth","ratio")) { ct<-cor.test(gene_sp[[v]], yv); m<-lm(yv~gene_sp[[v]])
    cat(sprintf("   %-6s r=%+.3f  R2=%.3f  p=%.4f  AIC=%.2f\n", v, ct$estimate, ct$estimate^2, ct$p.value, AIC(m))) } }
sink()
cat("Wrote S1S2_arcsinh_report.txt\n")
