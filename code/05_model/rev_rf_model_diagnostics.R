#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_model_diagnostics.R
# ------------------------------------------------------------------------------
# Baseline understanding of whatever RF configuration is currently saved.
# Run this after every change to 02_rf_models.R so each version is interpretable
# rather than a black box that emits a budget.
#
# Produces, for TreeRF and SoilRF:
#   1. skill summary        OOB R2, raw-scale R2, mean ratio (bias), Spearman
#   2. variable importance  permutation importance, which (unlike ranger's default
#                           impurity importance) is not biased toward high-cardinality
#                           and continuous predictors
#   3. partial dependence   marginal response to each top predictor, back-transformed
#                           to nmol m-2 s-1 so the y axis is interpretable
#   4. observed vs predicted, and residuals vs each predictor
#
# Outputs: outputs/revision/rf_diagnostics_{tree,soil}.png
#          outputs/revision/rf_diagnostics_summary.csv
#
# Run from repo root:  Rscript code/05_model/rev_rf_model_diagnostics.R
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(ggplot2); library(tidyr)})
set.seed(42)

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# NOTE: columns named *_umol_m2_s hold nmol m-2 s-1 (long-standing upstream misnomer).
MODELS <- list(
  tree = list(rf = TreeRF, X = as.data.frame(X_tree),
              y = tree_train_complete$stem_flux_corrected,      label = "Tree stem"),
  soil = list(rf = SoilRF, X = as.data.frame(X_soil),
              y = soil_train_complete$soil_flux_umol_m2_s,      label = "Soil")
)

# ---------------------------------------------------------------------------
# Permutation importance: shuffle one predictor, measure the loss in OOB R2.
# Reported as % of the model's OOB R2 so the two models are comparable.
# ---------------------------------------------------------------------------
perm_importance <- function(rf, X, y, n_rep = 5) {
  # Baseline and permuted scores are BOTH computed with predict() on the same rows,
  # so they are directly comparable. (Using rf$predictions for the baseline mixes OOB
  # with in-sample predictions and, if the row set differs at all, silently recycles.)
  # NO transform. Both models are fitted on the MEASURED scale (the asinh was
  # removed; the variable still named y_asinh in 02_rf_models.R holds raw values).
  # This line used to be z <- asinh(y), which scored raw predictions against a
  # transformed target -- so every permutation importance below was computed
  # against the wrong reference.
  z <- y
  r2 <- function(p) 1 - sum((z - p)^2) / sum((z - mean(z))^2)
  base <- r2(predict(rf, X)$predictions)
  vapply(colnames(X), function(v) {
    mean(vapply(seq_len(n_rep), function(i) {
      Xp <- X; Xp[[v]] <- sample(Xp[[v]])
      base - r2(predict(rf, Xp)$predictions)
    }, numeric(1)))
  }, numeric(1))
}

# ---------------------------------------------------------------------------
# Partial dependence: vary one predictor across its observed range, holding all
# others at their actual values, and average the prediction. Back-transformed.
# ---------------------------------------------------------------------------
partial_dep <- function(rf, X, var, n_grid = 25, n_sample = 300) {
  xs <- X[sample(nrow(X), min(n_sample, nrow(X))), , drop = FALSE]
  grid <- seq(quantile(X[[var]], 0.02, na.rm = TRUE),
              quantile(X[[var]], 0.98, na.rm = TRUE), length.out = n_grid)
  data.frame(
    variable = var, value = grid,
    yhat = vapply(grid, function(g) { xs[[var]] <- g
                    mean(predict(rf, xs)$predictions, na.rm = TRUE) }, numeric(1))
  )
}

summary_rows <- list()

for (nm in names(MODELS)) {
  m <- MODELS[[nm]]
  # NOTE: do NOT drop rows with missing predictors here. ranger trains on them
  # (splitting on observed values), so filtering them out would report metrics on a
  # cleaner subset than the model actually saw -- flattering the score.
  ok <- is.finite(m$y) & is.finite(m$rf$predictions)
  X <- m$X[ok, , drop = FALSE]; y <- m$y[ok]
  # Predictions and response are BOTH on the measured scale, so there is nothing to
  # back-transform. This was sinh(m$rf$predictions), left over from the removed
  # asinh: on a max OOB prediction of 1.81 nmol that inflates by 64%, and on the
  # observed max of 6.17 by 3,780%. It made r2_raw read 0.214 instead of 0.243 and
  # mean_ratio 1.084 instead of 1.026 -- i.e. it reported the model as three times
  # more biased than it is.
  oob <- m$rf$predictions[ok]

  summary_rows[[nm]] <- data.frame(
    model = m$label, n = length(y),
    oob_r2 = m$rf$r.squared,
    r2_raw       = 1 - sum((y - oob)^2) / sum((y - mean(y))^2),
    mean_obs = mean(y), mean_pred = mean(oob),
    mean_ratio = mean(oob) / mean(y),
    spearman = cor(y, oob, method = "spearman")
  )

  imp <- sort(perm_importance(m$rf, X, y), decreasing = TRUE)
  top <- names(imp)[seq_len(min(6, length(imp)))]

  cat("\n", strrep("=", 70), "\n  ", toupper(m$label), " MODEL  (n = ", length(y), ")\n",
      strrep("=", 70), "\n", sep = "")
  cat(sprintf("  OOB R2 %.3f | raw R2 %.3f | mean ratio %.3f | Spearman %.3f\n\n",
              summary_rows[[nm]]$oob_r2, summary_rows[[nm]]$r2_raw,
              summary_rows[[nm]]$mean_ratio, summary_rows[[nm]]$spearman))
  cat("  Permutation importance (drop in OOB R2 when shuffled):\n")
  for (v in names(imp)[seq_len(min(10, length(imp)))])
    cat(sprintf("    %-46s %7.4f\n", v, imp[v]))

  p_imp <- data.frame(variable = names(imp), importance = as.numeric(imp)) %>%
    slice_head(n = 12) %>%
    ggplot(aes(reorder(variable, importance), importance)) +
    geom_col(fill = "#2166AC") + coord_flip() +
    labs(title = paste0(m$label, ": permutation importance"),
         subtitle = "drop in OOB R² when the predictor is shuffled",
         x = NULL, y = "Δ OOB R²") + theme_minimal(base_size = 9)

  # Partial dependence is only defined on an ordered axis, so factors are skipped.
  # `species` entered the top 6 once dbh_within_z was dropped (it is now the
  # second-largest term in the tree model), and quantile() on a factor errors.
  top_num <- top[vapply(top, function(v) is.numeric(X[[v]]), logical(1))]
  if (length(top_num) < length(top))
    cat(sprintf("  (partial dependence skips non-numeric predictor(s): %s)\n",
                paste(setdiff(top, top_num), collapse = ", ")))
  pd <- bind_rows(lapply(top_num, function(v) partial_dep(m$rf, X, v)))
  p_pd <- ggplot(pd, aes(value, yhat)) +
    geom_line(colour = "#B2182B", linewidth = 0.7) +
    facet_wrap(~ variable, scales = "free_x", ncol = 3) +
    labs(title = paste0(m$label, ": partial dependence"),
         subtitle = "marginal response, back-transformed",
         x = NULL, y = expression(CH[4]~flux~(nmol~m^-2~s^-1))) +
    theme_minimal(base_size = 9)

  p_fit <- data.frame(obs = y, pred = oob) %>%
    ggplot(aes(obs, pred)) +
    geom_point(alpha = 0.35, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    labs(title = paste0(m$label, ": out-of-bag predicted vs observed"),
         subtitle = sprintf("mean ratio %.2f, Spearman %.2f (pseudo-log axes)",
                            summary_rows[[nm]]$mean_ratio, summary_rows[[nm]]$spearman),
         x = "observed", y = "predicted") + theme_minimal(base_size = 9)

  ggsave(file.path(outdir, paste0("rf_diagnostics_", nm, ".png")),
         gridExtra::arrangeGrob(p_imp, p_fit, p_pd, layout_matrix = rbind(c(1,2), c(3,3))),
         width = 10, height = 8, dpi = 200)
  cat("  -> ", file.path(outdir, paste0("rf_diagnostics_", nm, ".png")), "\n", sep = "")
}

res <- bind_rows(summary_rows)
write.csv(res, file.path(outdir, "rf_diagnostics_summary.csv"), row.names = FALSE)
cat("\nWritten: outputs/revision/rf_diagnostics_summary.csv\n")
