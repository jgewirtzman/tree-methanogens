#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_transform_sensitivity.R
# ------------------------------------------------------------------------------
# Compares target-scale and back-transformation choices for the RF upscaling.
#
# WHY THIS EXISTS
#   The pipeline nominally applies an asinh transform to the flux target. But
#   asinh has a fixed knee at |x| ~ 1, so when the model ran in umol m-2 s-1
#   (median ~2.6e-5) asinh(y) was the IDENTITY to 7 significant figures -- the
#   transform did nothing and the model was fit on raw flux. Working in nmol
#   (the reported unit) makes asinh real, but then naive sinh() back-transformation
#   introduces a large Jensen bias: RF predicts the conditional MEAN on the
#   transformed scale, and sinh() of that is a median-like quantity, not the mean a
#   budget requires.
#
# WHAT IT COMPARES (features are held fixed; only the target treatment varies)
#   1. raw                  no transform
#   2. asinh + naive        asinh(y), back-transform sinh()          <- status quo
#   3. asinh + smearing     asinh(y), Duan (1983) smearing estimator
#   4. asinh(y/s) + smearing  s = MAD of training flux, then smearing
#
# DECISIVE METRIC
#   A budget needs an unbiased MEAN flux. So the headline diagnostic is
#   mean(back-transformed OOB prediction) / mean(observed) -- 1.00 is unbiased.
#   OOB R2 is reported on a common scale (raw nmol) so the rows are comparable.
#
# Inputs : outputs/models/TRAINING_DATA.RData  (X_tree, X_soil + training frames)
# Outputs: outputs/revision/rf_transform_sensitivity.csv (+ console table)
# ==============================================================================

suppressMessages({library(ranger); library(dplyr)})
set.seed(42)

ROOT <- "."   # run from repo root
load(file.path(ROOT, "outputs/models/TRAINING_DATA.RData"))

# NOTE: columns named *_umol_m2_s hold nmol m-2 s-1 (long-standing misnomer upstream).
tree_y <- tree_train_complete$stem_flux_corrected
soil_y <- soil_train_complete$soil_flux_umol_m2_s

# ------------------------------------------------------------------------------
# Duan (1983) smearing estimator.
# For a model fitted on g(y), the conditional mean on the raw scale is estimated as
# the average of g_inv(prediction + residual_i) over the empirical residuals -- this
# restores the mass that Jensen's inequality removes from a naive g_inv(prediction).
# ------------------------------------------------------------------------------
smear_backtransform <- function(pred_t, resid_t, ginv) {
  vapply(pred_t, function(p) mean(ginv(p + resid_t)), numeric(1))
}

fit_config <- function(X, y, label, transform = c("raw","asinh","asinh_scaled"),
                       backtransform = c("none","naive","smearing")) {
  transform <- match.arg(transform); backtransform <- match.arg(backtransform)

  s <- 1
  if (transform == "raw") {
    g <- function(v) v;            ginv <- function(v) v
  } else if (transform == "asinh") {
    g <- function(v) asinh(v);     ginv <- function(v) sinh(v)
  } else {
    s <- mad(y, na.rm = TRUE); if (!is.finite(s) || s <= 0) s <- sd(y, na.rm = TRUE)
    g <- function(v) asinh(v / s); ginv <- function(v) sinh(v) * s
  }

  ok <- is.finite(y) & stats::complete.cases(X)
  Xf <- as.data.frame(X)[ok, , drop = FALSE]; yf <- y[ok]
  z  <- g(yf)

  rf <- ranger(x = Xf, y = z, num.trees = 800, min.node.size = 5,
               mtry = max(1, floor(sqrt(ncol(Xf)))), num.threads = 1,
               oob.error = TRUE, seed = 42)

  oob_t <- rf$predictions                       # OOB predictions on the modelled scale
  resid <- z - oob_t

  pred_raw <- switch(backtransform,
    none     = oob_t,
    naive    = ginv(oob_t),
    smearing = smear_backtransform(oob_t, resid[is.finite(resid)], ginv)
  )

  # R2 computed on the RAW nmol scale so configurations are comparable
  ss_res <- sum((yf - pred_raw)^2); ss_tot <- sum((yf - mean(yf))^2)

  data.frame(
    dataset          = label,
    transform        = transform,
    backtransform    = backtransform,
    scale_s          = s,
    n                = length(yf),
    oob_r2_modelscale= rf$r.squared,
    r2_rawscale      = 1 - ss_res/ss_tot,
    mean_observed    = mean(yf),
    mean_predicted   = mean(pred_raw),
    mean_ratio       = mean(pred_raw)/mean(yf),      # 1.00 = unbiased  <-- decisive
    spearman         = suppressWarnings(cor(yf, pred_raw, method = "spearman")),
    stringsAsFactors = FALSE
  )
}

CONFIGS <- list(
  list("raw",          "none"),
  list("asinh",        "naive"),
  list("asinh",        "smearing"),
  list("asinh_scaled", "smearing")
)

res <- bind_rows(lapply(CONFIGS, function(cf)
         bind_rows(fit_config(X_tree, tree_y, "TREE", cf[[1]], cf[[2]]),
                   fit_config(X_soil, soil_y, "SOIL", cf[[1]], cf[[2]]))))
res <- res %>% arrange(dataset, transform, backtransform)

cat("\n", strrep("=", 92), "\n", sep = "")
cat("  RF TARGET-TRANSFORM SENSITIVITY  (features fixed; nmol m-2 s-1)\n")
cat("  mean_ratio = mean(backtransformed OOB prediction) / mean(observed);  1.000 = unbiased\n")
cat(strrep("=", 92), "\n\n", sep = "")
print(res %>%
        mutate(across(c(oob_r2_modelscale, r2_rawscale, mean_ratio, spearman), ~round(.x, 3)),
               across(c(mean_observed, mean_predicted), ~signif(.x, 4)),
               scale_s = signif(scale_s, 4)) %>%
        dplyr::select(dataset, transform, backtransform, n, oob_r2_modelscale,
                      r2_rawscale, mean_observed, mean_predicted, mean_ratio, spearman),
      row.names = FALSE)

outdir <- file.path(ROOT, "outputs/revision")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(outdir, "rf_transform_sensitivity.csv"), row.names = FALSE)
cat("\nWritten: outputs/revision/rf_transform_sensitivity.csv\n")
