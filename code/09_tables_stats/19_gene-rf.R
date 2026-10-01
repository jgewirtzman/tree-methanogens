source("code/lib/outputs.R")
# ==============================================================================
# 19_gene-rf.R -- do gene abundances predict tree-level flux out of sample?
# ------------------------------------------------------------------------------
# Replaces the gene random forests that stood in the SI tree-level methods until 2026-10-01
# ("49.9% of variance", "46.2%", H-statistics), which no script in the pipeline
# reproduced. Those were out-of-bag or in-sample figures on ~120 trees with many
# predictors; the question that matters is narrower and testable: beyond species
# and environment, do the genes improve prediction of a held-out tree's flux?
#
# Unit: one row per tree (the 121 trees with mcrA, pmoA, mmoX and flux; the table
# behind Fig. S21 and Methods S2). Response: tree-mean flux on the arcsinh scale
# (cofactor 0.1), as in Fig. S21. Environment: DBH and soil VWC at the tree.
# Genes: log10 area-weighted mcrA, pmoA, mmoX.
# Skill: repeated 5-fold cross-validation (30 repeats), R2 = 1 - SSE/SST on the
# held-out trees, with the same folds for every model so differences are paired.
# Outputs: outputs/data/gene_rf_cv.csv, outputs/audit/gene_rf_cv.txt
# ==============================================================================
suppressWarnings(suppressMessages(source("code/07_molecular/helper_scale_dependent_gene_patterns.R")))
suppressPackageStartupMessages({ library(dplyr); library(ranger) })

asf <- function(x) asinh(x / 0.1) / log(10)
env <- read.csv("data/processed/integrated/merged_tree_dataset_final.csv") %>%
  transmute(tree_id, dbh = as.numeric(dbh), vwc = as.numeric(VWC_mean)) %>%
  distinct(tree_id, .keep_all = TRUE)
D <- tree_level_complete %>%
  transmute(tree_id, species = factor(species), y = asf(CH4_flux),
            lmcra = log10(mcrA + 1), lpmoa = log10(pmoA + 1), lmmox = log10(mmoX + 1)) %>%
  left_join(env, by = "tree_id") %>%
  filter(complete.cases(.))
stopifnot(!anyDuplicated(D$tree_id))

MODELS <- list(
  species              = c("species"),
  `species+env`        = c("species", "dbh", "vwc"),
  `species+env+genes`  = c("species", "dbh", "vwc", "lmcra", "lpmoa", "lmmox"),
  `env+genes`          = c("dbh", "vwc", "lmcra", "lpmoa", "lmmox"),
  genes                = c("lmcra", "lpmoa", "lmmox"))

cv_r2 <- function(vars, folds) {
  pred <- rep(NA_real_, nrow(D))
  for (k in unique(folds)) {
    tr <- folds != k
    m <- ranger(x = D[tr, vars, drop = FALSE], y = D$y[tr], num.trees = 500,
                min.node.size = 5, seed = 1, num.threads = 1)
    pred[!tr] <- predict(m, D[!tr, vars, drop = FALSE], num.threads = 1)$predictions
  }
  1 - sum((D$y - pred)^2) / sum((D$y - mean(D$y))^2)
}
set.seed(42); NREP <- 30
FOLDS <- replicate(NREP, sample(rep(1:5, length.out = nrow(D))), simplify = FALSE)
R <- sapply(names(MODELS), function(nm) sapply(FOLDS, function(f) cv_r2(MODELS[[nm]], f)))

S <- data.frame(model = colnames(R), n_trees = nrow(D),
                cv_r2_mean = colMeans(R), cv_r2_sd = apply(R, 2, sd))
gain <- R[, "species+env+genes"] - R[, "species+env"]
S$gain_over_species_env_mean <- NA; S$gain_frac_repeats_positive <- NA
S[S$model == "species+env+genes", c("gain_over_species_env_mean", "gain_frac_repeats_positive")] <-
  c(mean(gain), mean(gain > 0))
write.csv(S, out_path("gene_rf_cv.csv"), row.names = FALSE)

# permutation importance of the full model, fitted to all trees
full <- ranger(x = D[, MODELS[["species+env+genes"]]], y = D$y, num.trees = 1000,
               min.node.size = 5, importance = "permutation", seed = 1, num.threads = 1)
imp <- sort(full$variable.importance, decreasing = TRUE)

sink(out_path("gene_rf_cv.txt"))
cat("GENE RANDOM FORESTS: held-out skill on", nrow(D), "trees,", NREP, "x 5-fold CV (same folds for every model)\n\n")
print(S[, c("model", "cv_r2_mean", "cv_r2_sd")], row.names = FALSE, digits = 3)
cat(sprintf("\nadding genes to species+env: mean change in held-out R2 %+.3f; improved in %.0f%% of repeats\n",
            mean(gain), 100 * mean(gain > 0)))
cat(sprintf("out-of-bag R2 of the full model fitted to all trees: %.3f (optimistic; for comparison only)\n",
            full$r.squared))
cat("\npermutation importance, full model:\n"); print(round(imp, 4))
sink()
cat(readLines(out_path("gene_rf_cv.txt")), sep = "\n")
