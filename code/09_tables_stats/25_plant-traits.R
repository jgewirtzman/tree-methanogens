#!/usr/bin/env Rscript
# ==============================================================================
# 25_plant-traits.R -- species traits vs the stem methane-cycling community
# (Results §9, Methods S5, Fig S25 via code/08_figures/figS28_plant-traits.R)
# ------------------------------------------------------------------------------
# Unit: the ten gene-flux species (at least five trees with genes and five flux
# measurements). Responses: species medians of area-weighted mcrA, pmoA + mmoX,
# the methanogen:methanotroph balance, and net stem flux.
#
# Traits are chosen by rules set before analysis (Methods S5), applied to the
# canonical table data/tables/species_traits.csv (code/tables/07_species_traits.R, from
# archived TRY and GWDD records; conifer status from species.csv):
#   R1 data for at least 8 of the 10 species;
#   R2 at least 5 distinct values, or at least 3 species in the smallest class of a
#      categorical trait (so that no single species carries a result). Ectomycorrhizal
#      is categorical: EcM if the dataset-weighted EcM share is >= 0.5, else AM;
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
# heartwood C (R3, kept as stem C), Huber value (R1), maximum rooting depth (R1: 6
# of 10 once USDA minimum-depth records are dropped; checked here), wood pH (no
# species-level source for the analysis species). Conifer status is not a trait row.
# R1 and R2 are applied here to the candidate columns; failures are listed in
# traits_included.csv and left out.
#
# Analysis: Spearman correlation of each trait with each response across all ten
# species and, because the two conifers are extremes in trait space, across the
# eight broadleaf species alone. Benjamini-Hochberg FDR across the tests of each set
# (traits x 4 responses). Each correlation is also recomputed with each species left
# out in turn, and the direction of the conifer pair (P. strobus vs T. canadensis) is
# recorded. An association is reported when (A) p < 0.05 across all ten and among
# broadleaf species, with the same sign; or (B) p < 0.05 among broadleaf species, the
# same sign across all ten, and the two conifers differ in the direction of the
# broadleaf trend. The criterion met is recorded (column `criterion`).
# A multivariate check (RDA on trait PC1-3, permutation; PLS judged by leave-one-out
# CV) is run on all species, on broadleaf species, and on the traits complete for
# all ten (missing values median-imputed otherwise).
# Writes outputs/data/traits_included.csv, outputs/data/traits_correlations.csv,
#        outputs/audit/plant_traits.txt
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages({ library(tidyverse); library(vegan); library(pls) })
set.seed(42)

TRAITS <- "data/tables/species_traits.csv"

invisible(capture.output(suppressMessages(source("code/lib/species_gene_data.R"))))
resp <- analysis_ratio %>% transmute(species_id, flux = median_flux, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id")
sp10 <- resp$species_id
Y <- resp %>% column_to_rownames("species_id") %>% .[sp10, c("mcra", "meth", "balance", "flux")] %>% as.matrix()

traits <- tribble(~col, ~trait, ~group, ~type,
  "wood_density_gcm3",            "Wood density",             "Tissue density & bark",  "continuous",
  "bark_density_gcm3",            "Bark density",             "Tissue density & bark",  "continuous",
  "try_bark_volume_rel_wood",     "Bark volume (rel. wood)",  "Tissue density & bark",  "continuous",
  "try_wood_CN_ratio",            "Wood C:N",                 "Wood chemistry & decay", "continuous",
  "try_stem_C",                   "Stem C",                   "Wood chemistry & decay", "continuous",
  "try_bark_C",                   "Bark C",                   "Wood chemistry & decay", "continuous",
  "try_CWD_stem_decomp_rate_k",   "Wood decay rate",          "Wood chemistry & decay", "continuous",
  "vwc_realized",                 "Soil-moisture niche",      "Hydrology",              "continuous",
  "try_waterlogging_tolerance",   "Waterlogging tolerance",   "Hydrology",              "continuous",
  "try_rooting_depth",            "Rooting depth",            "Hydrology",              "continuous",
  "try_fine_root_diameter",       "Fine-root diameter",       "Roots & symbionts",      "continuous",
  "try_fine_root_tissue_density", "Fine-root tissue density", "Roots & symbionts",      "continuous",
  "try_ectomycorrhizal",          "Ectomycorrhizal",          "Roots & symbionts",      "categorical",
  "try_plant_longevity",          "Longevity",                "Life history",           "continuous")
tr <- read_csv(TRAITS, show_col_types = FALSE) %>% rename(spcode = species_code)
tr <- tr[match(sp10, tr$spcode), ]
niche <- read_csv("outputs/data/tree_species_moisture_niche.csv", show_col_types = FALSE)
tr$vwc_realized <- niche$vwc_mean[match(SPECIES_NAMES[tr$spcode], niche$species)]
tr$try_ectomycorrhizal <- as.numeric(tr$try_ectomycorrhizal >= 0.5)   # EcM 1, AM 0
g <- as.numeric(read_table("species")$gymnosperm[match(tr$spcode, read_table("species")$species_code)])
X <- as.data.frame(lapply(traits$col, function(c) as.numeric(tr[[c]]))); names(X) <- traits$trait; rownames(X) <- sp10

# R1, R2 on the candidate columns; failures are recorded and left out
inc <- traits %>% mutate(n = colSums(is.finite(as.matrix(X))),
  distinct = sapply(X, function(v) length(unique(v[is.finite(v)]))),
  smallest_class = ifelse(type == "categorical", sapply(X, function(v) min(table(v[is.finite(v)]))), NA),
  pass_R1 = n >= 8,
  pass_R2 = ifelse(type == "categorical", smallest_class >= 3, distinct >= 5),
  included = pass_R1 & pass_R2)
write.csv(inc, out_path("traits_included.csv"), row.names = FALSE)
traits <- traits %>% filter(inc$included); X <- X[, traits$trait, drop = FALSE]

# ---- univariate ------------------------------------------------------------------
sc <- function(x, y) { k <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[k], y[k], method = "spearman")); c(n = sum(k), rho = unname(ct$estimate), p = ct$p.value) }
uni <- expand_grid(trait = traits$trait, response = colnames(Y)) %>% rowwise() %>% mutate(
  a = list(sc(X[[trait]], Y[, response])), b = list(sc(X[[trait]][g == 0], Y[g == 0, response])),
  n = a[["n"]], rho = a[["rho"]], p = a[["p"]], n_bl = b[["n"]], rho_bl = b[["rho"]], p_bl = b[["p"]],
  jk = list(sapply(sp10[is.finite(X[[trait]])], function(s) { k <- sp10 != s; sc(X[[trait]][k], Y[k, response])[["p"]] })),
  dropone_max_p = max(jk), dropone_worst = names(jk)[which.max(jk)],
  conifer_dir = { d <- X[[trait]][sp10 == "TSCA"] - X[[trait]][sp10 == "PIST"]
    e <- Y[sp10 == "TSCA", response] - Y[sp10 == "PIST", response]; sign(d * e) },
  conifer_same_dir = conifer_dir == sign(rho), conifer_same_dir_bl = conifer_dir == sign(rho_bl)) %>%
  ungroup() %>% select(-a, -b, -jk, -conifer_dir) %>%
  mutate(q = p.adjust(p, "BH"), q_bl = p.adjust(p_bl, "BH"), group = traits$group[match(trait, traits$trait)]) %>%
  relocate(group, .after = trait)
# reporting rule (Methods S5)
uni <- uni %>% mutate(same_sign = sign(rho) == sign(rho_bl),
  criterion = case_when(p < 0.05 & p_bl < 0.05 & same_sign ~ "A: all ten and broadleaf",
                        p_bl < 0.05 & same_sign & conifer_same_dir_bl ~ "B: broadleaf, same sign across ten, conifers concordant",
                        .default = NA_character_),
  reported = !is.na(criterion)) %>% select(-same_sign)
write.csv(uni, out_path("traits_correlations.csv"), row.names = FALSE)
rep_ <- uni %>% filter(reported)

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
print(as.data.frame(inc %>% select(-col)), row.names = FALSE)
cat("Left out (R1/R2):", paste(inc$trait[!inc$included], collapse = "; "), "\n")
cat(sprintf("\nMinimum q: all species %.2f; broadleaf species %.2f\n", min(uni$q), min(uni$q_bl)))
cat(sprintf("Flux: minimum p all species %.3f; broadleaf %.3f\n", min(uni$p[uni$response == "flux"]), min(uni$p_bl[uni$response == "flux"])))
cat("\nCells with p < 0.05 in either set:\n")
print(as.data.frame(uni %>% filter(p < 0.05 | p_bl < 0.05) %>% arrange(p) %>%
  mutate(across(c(rho, rho_bl), ~ round(., 2)), across(c(p, p_bl, q, q_bl, dropone_max_p), ~ signif(., 2))) %>%
  select(trait, response, n, rho, p, q, n_bl, rho_bl, p_bl, q_bl, dropone_max_p, dropone_worst, conifer_same_dir)), row.names = FALSE)
cat("\nReported associations (A: p < 0.05 in both sets, same sign; B: broadleaf p < 0.05, same sign\n",
    "across all ten, conifer pair in the direction of the broadleaf trend):\n", sep = "")
print(as.data.frame(rep_ %>% mutate(across(c(rho, rho_bl), ~ round(., 2)), across(c(p, p_bl), ~ signif(., 2))) %>%
  select(trait, response, rho, p, rho_bl, p_bl, criterion)), row.names = FALSE)
cat("\nNot reported, broadleaf p < 0.05:", with(uni %>% filter(p_bl < 0.05, !reported),
  paste(trait, "-", response, collapse = "; ")), "\n")
cat("\nMultivariate check (RDA on trait PC1-3; PLS 1 component, leave-one-out R2):\n")
print(round(mv, 3))
sink()
cat(readLines(out_path("plant_traits.txt")), sep = "\n")
