#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_predictor_selection_current.R
# ------------------------------------------------------------------------------
# Predictor selection for the CURRENT TreeRF specification. Supersedes
# rev_rf_predictor_selection.R and rev_rf_predictor_significance.R for the
# readmission question.
#
# WHY A NEW PASS. Every earlier predictor test in this repo was run against a
# specification that no longer exists:
#   * soil temperature was a 10-value step function of month; it now carries 182
#     distinct per-measurement values (paired soil, 2026-07-29), so any earlier
#     verdict on "temperature" was a verdict on month wearing temperature's name;
#   * the 2021 campaign carried a single fabricated date, so month and the monthly
#     drivers were wrong for 46% of that campaign's rows;
#   * the species ladder was not reaching the inventory.
# Those tests are void, not merely dated.
#
# WHY THE OLD SCORING WAS WRONG TOO. rev_rf_predictor_selection.R used NO held-out
# data at all (20 in-sample refits against random reference columns), and
# rev_rf_predictor_significance.R scored a binomial test on how many of 30 SEEDS
# beat a shadow -- a statement about the fitting procedure's reproducibility, not
# about the data, which drives its p-value to zero as NRUN grows.
#
# HOW THIS SCORES. Two numbers, because the model has two jobs:
#   1. grouped CV R2, folds assigned by TREE. 1,130 rows come from 478 trees, so an
#      ungrouped fold puts the same stem in train and test and inflates skill ~2x.
#   2. the SUM ratio: mean held-out prediction / mean observed. The quantity the
#      paper reports is a stand total, and a predictor can improve per-stem R2 while
#      making the aggregate worse.
# Both with a standard error across repeats, so differences can be read against
# their own noise instead of ranked on the third decimal.
#
# Output: outputs/revision/rf_predictor_selection_current.csv / .txt
# ==============================================================================
suppressPackageStartupMessages({library(ranger)})

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")

# 30 repeats -- see rev_rf_grouped_cv.R. At NREP = 5 the between-variant ranking here
# was dominated by fold noise; differences must be read against the SE, not ranked.
NREP <- 30; NFOLD <- 5
X0  <- as.data.frame(X_tree)
y   <- tree_train_complete$y_asinh
grp <- as.character(tree_train_complete$tree_id)
mo  <- tree_train_complete$month
stopifnot(identical(colnames(X0), TreeRF$forest$independent.variable.names))

# Candidate specifications. "base" is what is currently locked.
# Variants are built by DROPPING from / ADDING to whatever the locked model actually
# has. air_temp_C_mean was removed from the specification on 2026-07-30 (decision A),
# so any variant defined as setdiff(X0, "air_temp_C_mean") is now a silent no-op --
# it would duplicate `base` and be reported as a distinct comparison. Guard against
# that by constructing readmissions explicitly and asserting every variant is unique.
ADDABLE <- list(air_temp_C_mean = tree_train_complete$air_temp_C_mean,
                month           = mo)
VARIANTS <- c(
  list(base = colnames(X0)),
  setNames(lapply(colnames(X0), function(v) setdiff(colnames(X0), v)),
           paste0("drop_", colnames(X0))),
  setNames(lapply(names(ADDABLE), function(v) c(colnames(X0), v)),
           paste0("add_", names(ADDABLE))),
  list(add_both = c(colnames(X0), names(ADDABLE)))
)
if ("soil_temp_C_mean" %in% colnames(X0))
  VARIANTS$drop_temp_add_air <- c(setdiff(colnames(X0), "soil_temp_C_mean"), "air_temp_C_mean")
stopifnot(!any(duplicated(lapply(VARIANTS, sort))))

score <- function(cols) {
  X <- X0
  for (v in names(ADDABLE)) if (v %in% cols && !v %in% names(X)) X[[v]] <- ADDABLE[[v]]
  X <- X[, cols, drop = FALSE]
  stopifnot(all(sapply(X, function(z) !all(is.na(z)))))
  # hyperparameters from the FITTED object, not retyped, so this cannot silently
  # score a different configuration than the one that is locked
  mt <- max(1, floor(sqrt(ncol(X)))); nt <- TreeRF$num.trees; nn <- TreeRF$min.node.size
  r2 <- ratio <- numeric(NREP)
  for (rep in seq_len(NREP)) {
    ut <- unique(grp); set.seed(2000 + rep)
    fold <- setNames(sample(rep(seq_len(NFOLD), length.out = length(ut))), ut)[grp]
    pr <- rep(NA_real_, nrow(X))
    for (k in seq_len(NFOLD)) {
      tr <- which(fold != k); te <- which(fold == k)
      m <- ranger(x = X[tr, , drop = FALSE], y = y[tr], num.trees = nt,
                  min.node.size = nn, mtry = mt, num.threads = 1, seed = 42)
      pr[te] <- predict(m, X[te, , drop = FALSE])$predictions
    }
    ok <- is.finite(pr) & is.finite(y)
    r2[rep]    <- 1 - sum((y[ok]-pr[ok])^2)/sum((y[ok]-mean(y[ok]))^2)
    ratio[rep] <- mean(pr[ok])/mean(y[ok])
  }
  # OOB on the full data, for comparison with what gets quoted
  oob <- ranger(x = X, y = y, num.trees = nt, min.node.size = nn, mtry = mt,
                num.threads = 1, seed = 42)$r.squared
  data.frame(n_pred = ncol(X), oob_r2 = oob,
             grouped_r2 = mean(r2), grouped_r2_se = sd(r2)/sqrt(NREP),
             sum_ratio = mean(ratio), sum_ratio_se = sd(ratio)/sqrt(NREP))
}

cat(sprintf("Grouped CV by tree: %d rows, %d trees, %d repeats x %d folds\n\n",
            length(y), length(unique(grp)), NREP, NFOLD))
RES <- do.call(rbind, lapply(names(VARIANTS), function(nm) {
  cat("  scoring", nm, "...\n")
  cbind(variant = nm, score(VARIANTS[[nm]]))
}))

base <- RES[RES$variant == "base", ]
RES$d_grouped_r2 <- RES$grouped_r2 - base$grouped_r2
# is the change larger than the noise on the difference?
RES$exceeds_se <- abs(RES$d_grouped_r2) > sqrt(RES$grouped_r2_se^2 + base$grouped_r2_se^2)
RES <- RES[order(-RES$grouped_r2), ]
write.csv(RES, "outputs/data/rf_predictor_selection_current.csv", row.names = FALSE)

con <- file("outputs/audit/rf_predictor_selection_current.txt", "w")
wr <- function(...) cat(sprintf(...), file = con)
wr("PREDICTOR SELECTION, CURRENT SPECIFICATION\n")
wr("Grouped CV by tree (%d rows from %d trees), %d repeats x %d folds.\n",
   length(y), length(unique(grp)), NREP, NFOLD)
wr("sum_ratio = mean held-out prediction / mean observed; 1.000 is unbiased in the sum.\n\n")
wr("%-20s %5s %8s %10s %8s %10s %8s %6s\n", "variant","npred","OOB R2","grouped R2",
   "+/-SE","sum ratio","+/-SE","d>SE")
for (i in seq_len(nrow(RES))) with(RES[i, ],
  wr("%-20s %5d %8.4f %10.4f %8.4f %10.4f %8.4f %6s\n", variant, n_pred, oob_r2,
     grouped_r2, grouped_r2_se, sum_ratio, sum_ratio_se, ifelse(exceeds_se,"yes","no")))
wr("\nRead the 'd>SE' column before the ranking: a difference smaller than the SE on\n")
wr("the difference is not a result. The base specification is the locked model.\n")
close(con)
cat(readLines("outputs/audit/rf_predictor_selection_current.txt"), sep = "\n")
cat("\n  Written: outputs/revision/rf_predictor_selection_current.{csv,txt}\n")
