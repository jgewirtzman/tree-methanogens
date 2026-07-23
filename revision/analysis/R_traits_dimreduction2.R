# ==============================================================================
# REVISION — Dimension reduction, part 2: permutation-based multivariate co-structure
# ==============================================================================
# Complements PCA (structure) and PLS (prediction, which failed) with the question
# they don't answer: is there SIGNIFICANT multivariate co-structure between trait
# space and the methane-response space? Permutation tests are the right tool at
# n=10 species (valid, distribution-free, don't overfit):
#   * RDA  -- constrained ordination: how much joint response variance traits explain
#            (traits summarised as PC1-3 to avoid over-parameterising at n=10).
#   * Mantel -- trait-distance vs response-distance correlation.
#   * Procrustes/PROTEST -- do the trait and response ordinations align.
#
# NOT run (indefensible at n<=16): t-SNE / UMAP (manifold methods need hundreds of
# points -> manufacture artifacts), sparse/regularised PCA, factor analysis.
#
# Responses = manuscript canonical species vars (via _prep_species_data.R).
# NEW file. Run: Rscript revision/analysis/R_traits_dimreduction2.R
# Output: revision/outputs/traits_dimreduction2.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(vegan) })
set.seed(42)
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

# ---- response matrix (species x {mcrA, methanotroph, balance, flux}) ----------
source("revision/analysis/_prep_species_data.R")
resp <- analysis_ratio %>% transmute(species_id, flux = median_flux, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id") %>%
  drop_na()
sp10 <- resp$species_id

# ---- trait matrix (curated, complete, scaled) --------------------------------
tr <- read_csv(TRAITS, show_col_types = FALSE) %>% filter(spcode %in% sp10) %>%
  group_by(spcode) %>% slice(1) %>% ungroup()
niche <- read_csv("revision/outputs/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(spcode = names(sp_map)[match(species, sp_map)]) %>% transmute(spcode, vwc_realized = vwc_mean)
use <- c("wood_density_gcm3","bark_density_gcm3","porosity_num","gymnosperm",
         "wood_sapwood_pH","wood_heartwood_pH","try_antifungal","try_rooting_depth",
         "try_plant_height","try_plant_longevity")
X <- tr %>% select(spcode, all_of(use)) %>% left_join(niche, by = "spcode") %>%
  column_to_rownames("spcode")
X <- X[sp10, ]
for (j in names(X)) X[[j]][!is.finite(X[[j]])] <- median(X[[j]], na.rm = TRUE)
Xs <- scale(X)
Y  <- resp %>% column_to_rownames("species_id") %>% .[sp10, ] %>% scale()

# ---- trait PCs (to constrain RDA without over-parameterising) ----------------
tpca <- prcomp(Xs)
tPC  <- as.data.frame(tpca$x[, 1:3]); names(tPC) <- c("tPC1","tPC2","tPC3")

sink(file.path(out, "traits_dimreduction2.txt"))
cat("=================================================================\n")
cat("Multivariate trait <-> methane-response co-structure (permutation)\n")
cat("n =", nrow(Y), "species; responses = mcrA, methanotroph, balance, flux\n")
cat("=================================================================\n\n")

# ---- 1. RDA: responses constrained by trait PC1-3 ----------------------------
rd <- rda(Y ~ tPC1 + tPC2 + tPC3, data = tPC)
a  <- anova.cca(rd, permutations = 999)
cat("1) RDA  (responses ~ trait PC1-3):\n")
cat(sprintf("   constrained (traits explain) = %.1f%% of response variance; adjusted R2 = %.3f\n",
            100 * rd$CCA$tot.chi / rd$tot.chi, RsquareAdj(rd)$adj.r.squared))
cat(sprintf("   permutation test: F = %.2f, p = %.3f (999 perms)\n\n", a$F[1], a$`Pr(>F)`[1]))
# per-axis terms
at <- anova.cca(rd, by = "terms", permutations = 999)
cat("   by trait axis:\n")
for (i in 1:3) cat(sprintf("     tPC%d: F=%.2f p=%.3f\n", i, at$F[i], at$`Pr(>F)`[i]))

# ---- 2. Mantel: trait distance vs response distance --------------------------
mt <- mantel(dist(Xs), dist(Y), method = "spearman", permutations = 999)
cat(sprintf("\n2) Mantel (trait-dist vs response-dist): r = %.3f, p = %.3f\n", mt$statistic, mt$signif))

# ---- 3. Procrustes / PROTEST: trait vs response ordinations ------------------
rpca <- prcomp(Y)
pr <- protest(tpca$x[, 1:2], rpca$x[, 1:2], permutations = 999)
cat(sprintf("\n3) Procrustes (trait PCA vs response PCA alignment):\n"))
cat(sprintf("   correlation = %.3f, PROTEST p = %.3f\n", pr$t0, pr$signif))

cat("\n----------------------------------------------------------------\n")
cat("READING: all three ask whether trait space and methane-response space share\n")
cat("significant multivariate structure. p>0.05 => no significant joint co-structure\n")
cat("detectable at this n; consistent with the failed PLS prediction. Individual\n")
cat("a-priori links (density->methanotroph/balance; moisture->methanogen) can still\n")
cat("hold -- these omnibus tests are low-powered at n=", nrow(Y), " and test the JOINT structure.\n")
cat("NOT run: t-SNE/UMAP (need hundreds of points; artifacts at this n).\n")
sink()
cat(readLines(file.path(out, "traits_dimreduction2.txt")), sep = "\n")
