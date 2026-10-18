#!/usr/bin/env Rscript
# ==============================================================================
# 04_rf_grouped_cv.R -- held-out skill of the two forests, folds grouped by unit
# ------------------------------------------------------------------------------
# Stage B; run from the repository root. Out-of-bag R2 is not grouped: a bootstrap
# sample that omits a row still contains the same tree measured at other heights or
# in other months, so OOB scores the stem forest partly on trees it has seen. Here
# folds hold out whole trees (stem forest) or whole collars (soil forest): repeated
# 5-fold CV of the fitted specification (nothing is refitted or changed).
#
# Output: outputs/data/rf_grouped_cv.csv   (read by 06_upscale/03_budget_canonical.R)
#         outputs/audit/rf_grouped_cv.txt
# ==============================================================================
source("code/lib/rf_spec.R")   # RF_NUM_TREES, RF_THREADS
suppressPackageStartupMessages({library(ranger)})

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")

# 30 repeats, not 5. With 478 trees and a handful of high-flux stems dominating
# whichever fold receives them, the fold assignment is the largest source of variance:
# over 30 independent seeds the grouped R2 ranges 0.043-0.156. A 5-repeat SD came out
# at 0.0060 and understated the true spread (~0.025) by about four-fold, and this
# number is quoted in canonical_budget.csv and on the main-text model panel.
# 100, not 30 (2026-10-01): with a per-repeat SD of ~0.027 the 30-repeat mean has a
# standard error of ~0.005, enough to flip the second decimal between runs; at 100 it
# is ~0.003.
NREP <- 100; NFOLD <- 5
r2 <- function(a, b) 1 - sum((a - b)^2) / sum((a - mean(a))^2)

# Score one locked spec: OOB (as reported) vs grouped CV (honest), NREP x NFOLD.
grouped_cv <- function(X, y, grp, rf, label, unit) {
  X <- as.data.frame(X)
  keep <- is.finite(y) & !is.na(grp)
  X <- X[keep, , drop = FALSE]; y <- y[keep]; grp <- as.character(grp[keep])
  mt <- rf$mtry; nn <- rf$min.node.size; nt <- rf$num.trees

  per_rep <- numeric(NREP); ratio <- numeric(NREP)
  for (rep in seq_len(NREP)) {
    ut <- sort(unique(grp), method = "radix")   # fold draw independent of row order
    set.seed(1000 + rep)
    gm <- setNames(sample(rep(seq_len(NFOLD), length.out = length(ut))), ut)
    f  <- unname(gm[grp])
    pr <- rep(NA_real_, nrow(X))
    for (k in seq_len(NFOLD)) {
      tr <- which(f != k); te <- which(f == k)
      m <- ranger(x = X[tr, , drop = FALSE], y = y[tr], num.trees = nt,
                  min.node.size = nn, mtry = mt, num.threads = RF_THREADS, seed = 42)
      pr[te] <- predict(m, X[te, , drop = FALSE])$predictions
    }
    # R2 per repeat, then averaged. Averaging the PREDICTIONS across repeats and
    # scoring once inflates R2 (it averages away fold noise) -- that is what
    # 12_model_family_comparison.R:126 does, worth ~+0.010.
    per_rep[rep] <- r2(y, pr); ratio[rep] <- mean(pr) / mean(y)
  }
  data.frame(model = label, unit = unit, n_rows = length(y),
             n_units = length(unique(grp)), rows_per_unit = length(y)/length(unique(grp)),
             oob_r2 = rf$r.squared,
             grouped_r2 = mean(per_rep), grouped_r2_sd = sd(per_rep),
             inflation = rf$r.squared / mean(per_rep),
             mean_ratio = mean(ratio))
}

# X_tree / X_soil are the EXACT design matrices the locked forests were fitted on,
# saved alongside them. Use those rather than reassembling from the training frame,
# so this cannot silently score a different specification.
stopifnot(identical(colnames(X_tree), TreeRF$forest$independent.variable.names),
          identical(colnames(as.data.frame(X_soil)), SoilRF$forest$independent.variable.names),
          nrow(X_tree) == TreeRF$num.samples, nrow(X_soil) == SoilRF$num.samples)

# --- tree: group by tree_id ---------------------------------------------------
RES <- grouped_cv(X_tree, tree_train_complete$y_flux, tree_train_complete$tree_id,
                  TreeRF, "TreeRF", "tree_id")

# --- soil: group by collar ----------------------------------------------------
sc <- intersect(c("collar_id","site_id","site","collar"), names(soil_train_complete))
if (length(sc)) {
  RES <- rbind(RES, grouped_cv(as.data.frame(X_soil), soil_train_complete$y_flux,
                               soil_train_complete[[sc[1]]], SoilRF, "SoilRF", sc[1]))
} else {
  cat("NOTE: no collar identifier found in soil_train_complete; soil row omitted\n")
}

write.csv(RES, "outputs/data/rf_grouped_cv.csv", row.names = FALSE)

con <- file("outputs/audit/rf_grouped_cv.txt", "w")
wr <- function(...) cat(sprintf(...), file = con)
wr("GROUPED CROSS-VALIDATION OF THE LOCKED FORESTS\n")
wr("%d repeats x %d folds, folds assigned by sampling unit, not by row.\n", NREP, NFOLD)
wr("The models are NOT refitted or respecified; this scores what is locked.\n\n")
for (i in seq_len(nrow(RES))) with(RES[i, ], {
  wr("%s  (grouped by %s)\n", model, unit)
  wr("  %d rows from %d %ss = %.2f rows per unit\n", n_rows, n_units, unit, rows_per_unit)
  wr("  OOB R2 as reported   %.4f\n", oob_r2)
  wr("  grouped CV R2        %.4f  (sd %.4f over %d repeats)\n", grouped_r2, grouped_r2_sd, NREP)
  wr("  OOB overstates by    %.2fx\n", inflation)
  wr("  mean pred / mean obs %.4f\n\n", mean_ratio)
})
wr("Report the grouped figure, or report both with the repeated-measures structure\n")
wr("stated. OOB alone is not defensible where one tree contributes up to 13 rows.\n")
close(con)

cat(readLines("outputs/audit/rf_grouped_cv.txt"), sep = "\n")
cat("\n  Written: outputs/rf_grouped_cv.{csv,txt}\n")
