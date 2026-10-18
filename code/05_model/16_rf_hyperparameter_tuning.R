#!/usr/bin/env Rscript
# ==============================================================================
# 16_rf_hyperparameter_tuning.R -- what tuning mtry and min.node.size is worth
# ------------------------------------------------------------------------------
# Methods S6: "tuning these settings, with selection nested inside cross-validation,
# changed grouped R2 by about ...". This script is the source of that number.
#
# For each locked forest (02_rf_models.R; design matrices X_tree / X_soil and the
# arcsinh response saved with them) it scores, under GROUPED cross-validation (whole
# trees or soil collars held out, as 04_rf_grouped_cv.R):
#   untuned  the locked setting (the forest's own mtry and min.node.size)
#   tuned    mtry and min.node.size chosen INSIDE each outer training fold by an
#            inner grouped 5-fold CV over the grid below, then scored on the outer fold
# on identical outer folds, so the difference is paired. Selection inside the outer
# loop is what makes the tuned figure honest: picking the best grid cell on the same
# folds that report skill would overstate what tuning buys.
# Also reported: each grid cell's grouped R2 on the outer folds (descriptive only).
#
# Grid: mtry 1 .. number of predictors; min.node.size 1, 3, 5, 10, 20. num.trees is
# not tuned; NUM_TREES below is used for every fit here, tuned and untuned alike, so
# the comparison is like for like (the locked forests use RF_NUM_TREES = 10,000; at
# that size a full nested run would take hours).
# Outer folds: NREP repeats x 5 folds, drawn as in 04_rf_grouped_cv.R (seed 1000 + rep).
# R2 is computed per repeat and averaged, as there.
#
# Writes outputs/data/rf_hyperparameters.csv, outputs/audit/rf_hyperparameters.txt
# ==============================================================================
source("code/lib/outputs.R"); source("code/lib/rf_spec.R")   # RF_THREADS; RUN (code/lib/run_mode.R)
suppressPackageStartupMessages(library(ranger))
load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")
stopifnot(identical(colnames(X_tree), TreeRF$forest$independent.variable.names),
          identical(colnames(as.data.frame(X_soil)), SoilRF$forest$independent.variable.names),
          nrow(X_tree) == TreeRF$num.samples, nrow(X_soil) == SoilRF$num.samples)

# run sizes from code/lib/run_mode.R (full run: 1000 trees, 10 repeats, nodes 1, 3, 5, 10, 20)
NUM_TREES <- RUN$tuning_num_trees; NREP <- RUN$tuning_repeats; NFOLD <- 5; NINNER <- 5
NODES <- RUN$tuning_node_sizes
r2 <- function(a, b) 1 - sum((a - b)^2) / sum((a - mean(a))^2)
fit_pred <- function(X, y, tr, te, mtry, node) {
  m <- ranger(x = X[tr, , drop = FALSE], y = y[tr], num.trees = NUM_TREES, mtry = mtry,
              min.node.size = node, num.threads = RF_THREADS, seed = 42)
  predict(m, X[te, , drop = FALSE], num.threads = RF_THREADS)$predictions
}
group_folds <- function(grp, seed, K) {
  ut <- sort(unique(grp), method = "radix"); set.seed(seed)
  gm <- setNames(sample(rep(seq_len(K), length.out = length(ut))), ut); unname(gm[grp])
}

tune <- function(X, y, grp, rf, label) {
  X <- as.data.frame(X); keep <- is.finite(y) & !is.na(grp)
  X <- X[keep, , drop = FALSE]; y <- y[keep]; grp <- as.character(grp[keep])
  G <- expand.grid(mtry = seq_len(ncol(X)), min.node.size = NODES)
  cur <- which(G$mtry == rf$mtry & G$min.node.size == rf$min.node.size); stopifnot(length(cur) == 1)
  r2_untuned <- r2_tuned <- numeric(NREP); r2_grid <- matrix(NA, NREP, nrow(G)); picked <- character(0)
  for (rep in seq_len(NREP)) {
    fo <- group_folds(grp, 1000 + rep, NFOLD)
    p_tuned <- p_untuned <- rep(NA_real_, length(y)); p_grid <- matrix(NA_real_, length(y), nrow(G))
    for (k in seq_len(NFOLD)) {
      tr <- which(fo != k); te <- which(fo == k)
      # inner grouped CV on the outer training part only
      fi <- group_folds(grp[tr], 5000 + 10 * rep + k, NINNER)
      inner <- vapply(seq_len(nrow(G)), function(g) {
        p <- rep(NA_real_, length(tr))
        for (j in seq_len(NINNER)) p[fi == j] <- fit_pred(X[tr, , drop = FALSE], y[tr], which(fi != j), which(fi == j),
                                                       G$mtry[g], G$min.node.size[g])
        r2(y[tr], p)
      }, 0)
      b <- which.max(inner); picked <- c(picked, sprintf("mtry %d, node %d", G$mtry[b], G$min.node.size[b]))
      for (g in seq_len(nrow(G))) p_grid[te, g] <- fit_pred(X, y, tr, te, G$mtry[g], G$min.node.size[g])
      p_tuned[te] <- p_grid[te, b]; p_untuned[te] <- p_grid[te, cur]
    }
    r2_tuned[rep] <- r2(y, p_tuned); r2_untuned[rep] <- r2(y, p_untuned)
    r2_grid[rep, ] <- apply(p_grid, 2, function(p) r2(y, p))
    cat(sprintf("  %s repeat %2d: untuned %.4f  tuned %.4f\n", label, rep, r2_untuned[rep], r2_tuned[rep]))
  }
  d <- r2_tuned - r2_untuned
  list(summary = data.frame(model = label, n_rows = length(y), n_units = length(unique(grp)),
                            mtry_locked = rf$mtry, node_locked = rf$min.node.size, num_trees = NUM_TREES,
                            n_repeats = NREP, grouped_r2_untuned = mean(r2_untuned), grouped_r2_tuned = mean(r2_tuned),
                            tuning_delta_r2 = mean(d), tuning_delta_r2_sd = sd(d),
                            repeats_tuning_better = sum(d > 0)),
       grid = data.frame(model = label, G, grouped_r2 = colMeans(r2_grid),
                         grouped_r2_sd = apply(r2_grid, 2, sd), is_locked = seq_len(nrow(G)) == cur),
       picked = sort(table(picked), decreasing = TRUE))
}

cat("Nested grouped CV of mtry x min.node.size\n")
TR <- tune(X_tree, tree_train_complete$y_flux, tree_train_complete$tree_id, TreeRF, "TreeRF")
SO <- tune(X_soil, soil_train_complete$y_flux, soil_train_complete$site_id, SoilRF, "SoilRF")

SUM <- rbind(TR$summary, SO$summary)
write.csv(merge(rbind(TR$grid, SO$grid), SUM, by = "model"), out_path("rf_hyperparameters.csv"), row.names = FALSE)

con <- file(out_path("rf_hyperparameters.txt"), "w")
wr <- function(...) cat(sprintf(...), file = con)
wr("TUNING mtry AND min.node.size, SELECTION NESTED IN GROUPED CV\n")
wr("%d repeats x %d outer folds (whole trees / collars held out); inner grouped %d-fold; %d trees per fit\n\n",
   NREP, NFOLD, NINNER, NUM_TREES)
for (R in list(TR, SO)) {
  s <- R$summary
  wr("%s  (%d rows, %d units; locked mtry %d, min.node.size %d)\n", s$model, s$n_rows, s$n_units, s$mtry_locked, s$node_locked)
  wr("  grouped R2 untuned (locked setting) %.4f\n", s$grouped_r2_untuned)
  wr("  grouped R2 tuned (nested)           %.4f\n", s$grouped_r2_tuned)
  wr("  tuning changes grouped R2 by        %+.4f (sd %.4f over repeats; tuned better in %d of %d)\n",
     s$tuning_delta_r2, s$tuning_delta_r2_sd, s$repeats_tuning_better, s$n_repeats)
  wr("  settings picked in the outer folds: %s\n", paste(sprintf("%s (%d)", names(R$picked), R$picked), collapse = "; "))
  g <- R$grid[order(-R$grid$grouped_r2), ]
  wr("  grid, grouped R2 on the outer folds (descriptive; best first):\n")
  for (i in seq_len(nrow(g))) wr("    mtry %d, node %2d  %.4f (sd %.4f)%s\n", g$mtry[i], g$min.node.size[i],
                                 g$grouped_r2[i], g$grouped_r2_sd[i], ifelse(g$is_locked[i], "  <- locked", ""))
  wr("\n")
}
close(con)
cat(readLines(out_path("rf_hyperparameters.txt")), sep = "\n")
