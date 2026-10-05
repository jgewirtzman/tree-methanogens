#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 25_plant-traits.R -- species traits vs the stem methane-cycling community
# (Results §9, Methods S5, Fig S25 via code/08_figures/figS28_plant-traits.R)
# ------------------------------------------------------------------------------
# Unit: the ten gene-flux species (at least five trees with genes and five flux
# measurements). Responses: species medians of area-weighted mcrA, pmoA + mmoX,
# the methanogen:methanotroph balance, and net stem flux.
#
# Traits are chosen by rules set before analysis (Methods S5), applied to the
# species-level table in the companion tree-gas-traits repository:
#   R1 data for at least 8 of the 10 species;
#   R2 at least 5 distinct values, or at least 3 species in the smallest class of a
#      categorical trait (so that no single species carries a result);
#   R3 one column per property (duplicates from other sources, ratios of included
#      traits out);
#   R4 values plausible for adult trees;
#   R5 a stated route by which the trait could affect stem methane cycling.
# Excluded under these rules: TRY wood density and estimated trunk/branch specific
# gravity (R3, duplicates of wood density; TRY gives sugar maple 0.99), trunk:branch
# ratio (R5), bark:wood density ratio (R3), porosity (R2: 2 ring-porous species),
# antifungal chemistry and decay resistance (R2: one species each), wetland indicator
# (R2: 8 of 10 tied), mycorrhiza type and growth form (R2), an unlabelled TRY column,
# TRY plant height (R4: means of 2-35 m, e.g. sugar maple 2.5 m), sapwood and
# heartwood C (R3, kept as stem C), Huber value (R1). Conifer status is not a trait
# row. Only the 16 included columns (+ gymnosperm) are vendored, by
# code/tools/vendor_traits.R; this script re-checks R1 and R2 on them.
#
# Analysis: Spearman correlation of each trait with each response across all ten
# species and, because the two conifers are extremes in trait space, across the
# eight broadleaf species alone. Benjamini-Hochberg FDR across the 64 tests of each
# set. Each correlation is also recomputed with each species left out in turn, and
# the direction of the conifer pair (P. strobus vs T. canadensis) is recorded.
# A multivariate check (RDA on trait PC1-3, permutation; PLS judged by leave-one-out
# CV) is run on all species, on broadleaf species, and on the traits complete for
# all ten (missing values median-imputed otherwise).
# Writes outputs/data/traits_included.csv, outputs/data/traits_correlations.csv,
#        outputs/audit/plant_traits.txt
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(vegan); library(pls) })
set.seed(42)

TRAITS <- "data/raw/external/tree-gas-traits/ymf_species_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",BELE="Betula lenta",
  BEPA="Betula papyrifera",FAGR="Fagus grandifolia",FRAM="Fraxinus americana",
  PIST="Pinus strobus",QURU="Quercus rubra",TSCA="Tsuga canadensis")

invisible(capture.output(suppressMessages(source("code/lib/prep_species_data.R"))))
resp <- analysis_ratio %>% transmute(species_id, flux = median_flux, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id")
sp10 <- resp$species_id
Y <- resp %>% column_to_rownames("species_id") %>% .[sp10, c("mcra", "meth", "balance", "flux")] %>% as.matrix()

traits <- tribble(~col, ~trait, ~group,
  "wood_density_gcm3",            "Wood density",             "Tissue density & bark",
  "bark_density_gcm3",            "Bark density",             "Tissue density & bark",
  "try_bark_volume_rel_wood",     "Bark volume (rel. wood)",  "Tissue density & bark",
  "wood_sapwood_pH",              "Sapwood pH",               "Wood chemistry & decay",
  "wood_heartwood_pH",            "Heartwood pH",             "Wood chemistry & decay",
  "try_wood_CN_ratio",            "Wood C:N",                 "Wood chemistry & decay",
  "try_stem_C",                   "Stem C",                   "Wood chemistry & decay",
  "try_bark_C",                   "Bark C",                   "Wood chemistry & decay",
  "try_CWD_stem_decomp_rate_k",   "Wood decay rate",          "Wood chemistry & decay",
  "vwc_realized",                 "Soil-moisture niche",      "Hydrology",
  "try_waterlogging_tolerance",   "Waterlogging tolerance",   "Hydrology",
  "try_rooting_depth",            "Rooting depth",            "Hydrology",
  "try_fine_root_diameter",       "Fine-root diameter",       "Roots & symbionts",
  "try_fine_root_tissue_density", "Fine-root tissue density", "Roots & symbionts",
  "try_ectomycorrhizal",          "Ectomycorrhizal",          "Roots & symbionts",
  "try_plant_longevity",          "Longevity",                "Life history")
tr <- read_csv(TRAITS, show_col_types = FALSE)
tr <- tr[match(sp10, tr$spcode), ]
niche <- read_csv("outputs/data/tree_species_moisture_niche.csv", show_col_types = FALSE)
tr$vwc_realized <- niche$vwc_mean[match(sp_map[tr$spcode], niche$species)]
g <- tr$gymnosperm
X <- as.data.frame(lapply(traits$col, function(c) as.numeric(tr[[c]]))); names(X) <- traits$trait; rownames(X) <- sp10

# R1, R2 re-checked on the vendored columns
inc <- traits %>% mutate(n = colSums(is.finite(as.matrix(X))),
  distinct = sapply(X, function(v) length(unique(v[is.finite(v)]))),
  smallest_class = sapply(X, function(v) { v <- v[is.finite(v)]  # categorical: ectomycorrhizal is a 0-1 share (0.98-1 = ECM)
    if (length(unique(round(v))) <= 3 && max(abs(v - round(v))) < 0.05) min(table(round(v))) else NA }))
stopifnot(all(inc$n >= 8), all(inc$distinct >= 5 | (!is.na(inc$smallest_class) & inc$smallest_class >= 3)))
write.csv(inc, out_path("traits_included.csv"), row.names = FALSE)

# ---- univariate ------------------------------------------------------------------
sc <- function(x, y) { k <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[k], y[k], method = "spearman")); c(n = sum(k), rho = unname(ct$estimate), p = ct$p.value) }
uni <- expand_grid(trait = traits$trait, response = colnames(Y)) %>% rowwise() %>% mutate(
  a = list(sc(X[[trait]], Y[, response])), b = list(sc(X[[trait]][g == 0], Y[g == 0, response])),
  n = a[["n"]], rho = a[["rho"]], p = a[["p"]], n_bl = b[["n"]], rho_bl = b[["rho"]], p_bl = b[["p"]],
  jk = list(sapply(sp10[is.finite(X[[trait]])], function(s) { k <- sp10 != s; sc(X[[trait]][k], Y[k, response])[["p"]] })),
  dropone_max_p = max(jk), dropone_worst = names(jk)[which.max(jk)],
  conifer_same_dir = { d <- X[[trait]][sp10 == "TSCA"] - X[[trait]][sp10 == "PIST"]
    e <- Y[sp10 == "TSCA", response] - Y[sp10 == "PIST", response]; sign(d * e) == sign(rho) }) %>%
  ungroup() %>% select(-a, -b, -jk) %>%
  mutate(q = p.adjust(p, "BH"), q_bl = p.adjust(p_bl, "BH"), group = traits$group[match(trait, traits$trait)]) %>%
  relocate(group, .after = trait)
write.csv(uni, out_path("traits_correlations.csv"), row.names = FALSE)
both <- uni %>% filter(p < 0.05, p_bl < 0.05)

# ---- multivariate check -------------------------------------------------------------
Xi <- X; for (j in names(Xi)) Xi[[j]][!is.finite(Xi[[j]])] <- median(Xi[[j]], na.rm = TRUE)
complete <- names(X)[colSums(is.finite(as.matrix(X))) == length(sp10)]
mv_run <- function(k, cols) {
  Xs <- scale(Xi[k, cols]); Xs <- Xs[, apply(Xs, 2, function(v) all(is.finite(v))), drop = FALSE]; Ys <- scale(Y[k, ])
  tPC <- as.data.frame(prcomp(Xs)$x[, 1:3]); rd <- rda(Ys ~ PC1 + PC2 + PC3, data = tPC); a <- anova.cca(rd, permutations = 999)
  pls_r2 <- sapply(colnames(Y), function(rn) { m <- plsr(y ~ ., data = data.frame(y = Ys[, rn], Xs, check.names = FALSE),
    ncomp = 1, validation = "LOO", scale = FALSE); drop(pls::R2(m, estimate = "CV")$val)[2] })
  c(rda_adj_r2 = RsquareAdj(rd)$adj.r.squared, rda_p = a$`Pr(>F)`[1], setNames(pls_r2, paste0("loo_r2_", colnames(Y)))) }
mv <- rbind(`all species` = mv_run(rep(TRUE, 10), names(Xi)),
            `broadleaf species` = mv_run(g == 0, names(Xi)),
            `all species, complete traits` = mv_run(rep(TRUE, 10), complete),
            `all species without P. strobus` = mv_run(sp10 != "PIST", names(Xi)))

# ---- transcript -----------------------------------------------------------------------
sink(out_path("plant_traits.txt"))
cat(sprintf("Plant traits: %d species (%d broadleaf), %d traits, %d tests per set\n\n", length(sp10), sum(g == 0), nrow(traits), nrow(uni)))
print(as.data.frame(inc), row.names = FALSE)
cat(sprintf("\nMinimum q: all species %.2f; broadleaf species %.2f\n", min(uni$q), min(uni$q_bl)))
cat(sprintf("Flux: minimum p all species %.3f; broadleaf %.3f\n", min(uni$p[uni$response == "flux"]), min(uni$p_bl[uni$response == "flux"])))
cat("\nCells with p < 0.05 in either set:\n")
print(as.data.frame(uni %>% filter(p < 0.05 | p_bl < 0.05) %>% arrange(p) %>%
  mutate(across(c(rho, rho_bl), ~ round(., 2)), across(c(p, p_bl, q, q_bl, dropone_max_p), ~ signif(., 2))) %>%
  select(trait, response, n, rho, p, q, n_bl, rho_bl, p_bl, q_bl, dropone_max_p, dropone_worst, conifer_same_dir)), row.names = FALSE)
cat("\nHeld in both sets (p < 0.05):", paste(both$trait, "-", both$response, collapse = "; "), "\n")
cat("\nMultivariate check (RDA on trait PC1-3; PLS 1 component, leave-one-out R2):\n")
print(round(mv, 3))
sink()
cat(readLines(out_path("plant_traits.txt")), sep = "\n")
