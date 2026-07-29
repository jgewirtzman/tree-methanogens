#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_hyperparameters.R  --  does tuning change anything that matters?
# ------------------------------------------------------------------------------
# The tree model exists to produce a STAND SUM over 8,006 inventory stems, not a
# per-tree prediction, so it is scored on two criteria at once:
#
#   R2      stratified 5-fold CV, species-stratified, repeated -- ordinary skill
#   ratio   area-weighted sum ratio, sum(observed)/sum(predicted) over the same
#           held-out folds, weighted by stem band area. This is the quantity the
#           budget actually depends on, and a regression forest is biased in it:
#           squared-error minimisation shrinks predictions toward the conditional
#           mean, so low-flux species come out high and high-flux species low.
#
# A model that is marginally worse per tree but unbiased in aggregate is the
# better model here. That trade is the entire reason the per-species calibration
# exists, so the question worth asking of min.node.size in particular is whether
# a smaller node gets `ratio` near 1 on its own and leaves the calibration a
# smaller correction to make.
#
# TUNING IS NESTED INSIDE THE CV. Selecting hyperparameters on the same folds
# used to report skill inflates the reported R2. Here the outer loop holds out a
# fold, the inner loop picks hyperparameters using only the training part, and
# the outer fold is scored with that choice -- so `nested_R2` is an honest
# estimate of what tuning buys, and is directly comparable to the untuned model.
#
# Grid (num.trees is not tuned: at 800 the OOB error has flattened and more trees
# only reduce Monte Carlo noise; max.depth is left unset because it is redundant
# with min.node.size):
#   mtry           2 .. 6      currently floor(sqrt(6)) = 2
#   min.node.size  1,3,5,10,20 currently 5
#
# Output: outputs/revision/rf_hyperparameters.csv
# ==============================================================================
suppressMessages({library(dplyr); library(ranger)})
set.seed(42)
outdir <- "outputs/revision"
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")

d <- tree_train_complete
X <- as.data.frame(X_tree); y <- d$stem_flux_corrected
k <- is.finite(y) & complete.cases(X)
X <- X[k, , drop = FALSE]; y <- y[k]
sp <- as.character(d$species_clean)[k]
# band area of the measured stem, the weight the stand sum applies to it
A <- pi * X$dbh_m * ifelse(sp == "Kalmia latifolia", 0.75, 2)

GRID <- expand.grid(mtry = 2:min(6, ncol(X)), min.node.size = c(1, 3, 5, 10, 20))
CURRENT <- list(mtry = max(1, floor(sqrt(ncol(X)))), min.node.size = 5)
cat(sprintf("predictors: %s\n", paste(names(X), collapse = ", ")))
cat(sprintf("current setting: mtry = %d, min.node.size = %d\n",
            CURRENT$mtry, CURRENT$min.node.size))
cat(sprintf("grid: %d combinations, %d measurements\n\n", nrow(GRID), length(y)))

folds <- function(strat, K) ave(seq_along(strat), strat,
                                FUN = function(i) sample(rep_len(seq_len(K), length(i))))

# held-out predictions for one hyperparameter setting
oof <- function(g, fo, K) {
  p <- numeric(length(y))
  for (f in seq_len(K)) {
    m <- ranger(x = X[fo != f, , drop = FALSE], y = y[fo != f], num.trees = 800,
                mtry = g$mtry, min.node.size = g$min.node.size, num.threads = 1)
    p[fo == f] <- predict(m, X[fo == f, , drop = FALSE], num.threads = 1)$predictions
  }
  p
}
score <- function(p) c(R2 = 1 - sum((y - p)^2)/sum((y - mean(y))^2),
                       ratio = sum(A*y)/sum(A*p))

# ---- 1. flat grid, repeated CV -- what each setting is worth ----------------
NREP <- 3; K <- 5
res <- lapply(seq_len(nrow(GRID)), function(i) {
  v <- replicate(NREP, score(oof(GRID[i, ], folds(sp, K), K)))
  data.frame(GRID[i, ], R2 = mean(v["R2", ]), R2_se = sd(v["R2", ])/sqrt(NREP),
             ratio = mean(v["ratio", ]), ratio_se = sd(v["ratio", ])/sqrt(NREP))
}) %>% bind_rows()
res$is_current <- res$mtry == CURRENT$mtry & res$min.node.size == CURRENT$min.node.size

cat("=== grid, repeated 5-fold CV stratified by species ===\n")
cat("   ratio = area-weighted sum(obs)/sum(pred); 1.000 is unbiased in aggregate\n\n")
print(as.data.frame(res %>% arrange(-R2) %>%
  transmute(mtry, min.node.size, R2 = sprintf("%.4f", R2),
            `+/-` = sprintf("%.4f", R2_se), ratio = sprintf("%.3f", ratio),
            current = ifelse(is_current, "  <- current", ""))), row.names = FALSE)

best_r2 <- res[which.max(res$R2), ]
best_bal <- res[which.min(abs(res$ratio - 1)), ]
cat(sprintf("\nbest R2      : mtry %d, node %d -> R2 %.4f, ratio %.3f\n",
            best_r2$mtry, best_r2$min.node.size, best_r2$R2, best_r2$ratio))
cat(sprintf("least biased : mtry %d, node %d -> R2 %.4f, ratio %.3f\n",
            best_bal$mtry, best_bal$min.node.size, best_bal$R2, best_bal$ratio))
cur <- res[res$is_current, ]
cat(sprintf("current      : mtry %d, node %d -> R2 %.4f, ratio %.3f\n",
            cur$mtry, cur$min.node.size, cur$R2, cur$ratio))

# ---- 2. nested CV -- what tuning is honestly worth --------------------------
# Outer fold held out; hyperparameters chosen on the inner folds only.
cat("\n=== nested CV: tuning selected inside each outer fold ===\n")
KO <- 5; fo <- folds(sp, KO)
p_nest <- numeric(length(y)); picked <- character(KO)
for (f in seq_len(KO)) {
  tr <- fo != f
  Xi <- X[tr, , drop = FALSE]; yi <- y[tr]; si <- sp[tr]
  fi <- ave(seq_along(si), si, FUN = function(i) sample(rep_len(1:5, length(i))))
  inner <- sapply(seq_len(nrow(GRID)), function(i) {
    pi_ <- numeric(length(yi))
    for (g in 1:5) {
      m <- ranger(x = Xi[fi != g, , drop = FALSE], y = yi[fi != g], num.trees = 800,
                  mtry = GRID$mtry[i], min.node.size = GRID$min.node.size[i], num.threads = 1)
      pi_[fi == g] <- predict(m, Xi[fi == g, , drop = FALSE], num.threads = 1)$predictions
    }
    1 - sum((yi - pi_)^2)/sum((yi - mean(yi))^2)
  })
  b <- GRID[which.max(inner), ]
  picked[f] <- sprintf("mtry=%d,node=%d", b$mtry, b$min.node.size)
  m <- ranger(x = Xi, y = yi, num.trees = 800, mtry = b$mtry,
              min.node.size = b$min.node.size, num.threads = 1)
  p_nest[fo == f] <- predict(m, X[fo == f, , drop = FALSE], num.threads = 1)$predictions
}
s_nest <- score(p_nest)
p_cur <- oof(CURRENT, fo, KO); s_cur <- score(p_cur)
cat(sprintf("  selected per outer fold: %s\n", paste(picked, collapse = " | ")))
cat(sprintf("  untuned (mtry %d, node %d) : R2 %.4f, ratio %.3f\n",
            CURRENT$mtry, CURRENT$min.node.size, s_cur["R2"], s_cur["ratio"]))
cat(sprintf("  nested-tuned               : R2 %.4f, ratio %.3f\n", s_nest["R2"], s_nest["ratio"]))
cat(sprintf("  => tuning is worth %+.4f R2 once selection is paid for honestly\n",
            s_nest["R2"] - s_cur["R2"]))

res$nested_R2_untuned <- s_cur["R2"]; res$nested_R2_tuned <- s_nest["R2"]
res$nested_ratio_untuned <- s_cur["ratio"]; res$nested_ratio_tuned <- s_nest["ratio"]
write.csv(res, file.path(outdir, "rf_hyperparameters.csv"), row.names = FALSE)
cat(sprintf("\nwritten: %s/rf_hyperparameters.csv\n", outdir))
