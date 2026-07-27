#!/usr/bin/env Rscript
# ==============================================================================
# rev_model_family_comparison.R
# ------------------------------------------------------------------------------
# Is a random forest the right model for a STAND BUDGET?
#
# THE PROBLEM. A regression forest minimises squared error, and the minimiser of
# squared error shrinks predictions toward the conditional mean. That is correct
# behaviour for prediction and wrong behaviour for a SUM: the audit in
# rev_rf_species_bias_audit.R found low-flux species over-predicted (Tsuga 1.43x,
# Pinus strobus 1.44x) and high-flux species under-predicted (Betula
# alleghaniensis 0.83x, Fraxinus 0.79x). Because the two over-predicted conifers
# carry 54% of the band area, the stand tree term comes out 23% high.
#
# WHY A LINEAR MODEL MIGHT BE BETTER HERE, and it is a structural argument rather
# than an empirical hope: least squares with species as a factor forces the fitted
# mean to equal the observed mean WITHIN EACH LEVEL. The residuals sum to zero by
# construction, per species. That is precisely the property the budget needs and
# the one the forest cannot provide. The price is that a linear model cannot
# represent interactions or thresholds the forest finds, so it may predict
# individual stems worse.
#
# So the two model families are being asked for different things, and this script
# scores both on both:
#   RMSE / R2        can it predict an individual measurement
#   per-species bias can it get each species' MEAN right
#   stand budget     what it does to the number in the paper
#
# Candidates:
#   rf              the locked random forest
#   rf_calibrated   the forest with a per-species multiplicative correction fitted
#                   on the training fold only (so the correction is itself CV'd)
#   lm              least squares, species factor + environment, no interactions
#   lm_log          least squares on log(flux - min + eps), back-transformed with
#                   a smearing estimator (Duan) rather than naive exponentiation
#   lmm             mixed model, random species intercept -- partial pooling
#
# Output: outputs/revision/model_family_comparison.csv / .txt
# ==============================================================================
suppressPackageStartupMessages({library(dplyr); library(ranger)})
HAVE_LME4 <- requireNamespace("lme4", quietly = TRUE)
set.seed(42)
load("outputs/models/TRAINING_DATA.RData")
source("code/revision/rev_geometry.R")
CONV <- 86400 * 365.25 * 16e-6

X0 <- as.data.frame(X_tree)
PREDS <- setdiff(names(X0), "species")
d <- tree_train_complete
stopifnot(nrow(X0) == nrow(d))
keep <- is.finite(d$stem_flux_corrected) & complete.cases(X0[, PREDS, drop = FALSE])
X0 <- X0[keep, , drop = FALSE]; d <- d[keep, ]
y <- d$stem_flux_corrected
sp <- as.character(X0$species)
X0$species <- factor(sp)
cat(sprintf("rows %d | species levels %d | y: mean %.4f median %.4f max %.2f\n\n",
            nrow(X0), nlevels(X0$species), mean(y), median(y), max(y)))

FIT <- list(
  rf = function(tr, te) {
    m <- ranger(x = X0[tr, , drop = FALSE], y = y[tr], num.trees = 800,
                min.node.size = 5, mtry = max(1, floor(sqrt(ncol(X0)))),
                num.threads = 1, seed = 42)
    predict(m, X0[te, , drop = FALSE])$predictions
  },
  rf_calibrated = function(tr, te) {
    m <- ranger(x = X0[tr, , drop = FALSE], y = y[tr], num.trees = 800,
                min.node.size = 5, mtry = max(1, floor(sqrt(ncol(X0)))),
                num.threads = 1, seed = 42)
    # per-species ratio from OOB predictions on the TRAINING fold only
    cal <- data.frame(s = sp[tr], o = y[tr], p = m$predictions) %>%
      group_by(s) %>% summarise(r = ifelse(mean(p) > 0, mean(o)/mean(p), 1), .groups = "drop")
    cm <- setNames(cal$r, cal$s)
    p <- predict(m, X0[te, , drop = FALSE])$predictions
    f <- ifelse(is.na(cm[sp[te]]), 1, cm[sp[te]])
    p * pmax(pmin(as.numeric(f), 5), 0.2)      # guard against wild ratios on n=1 levels
  },
  lm = function(tr, te) {
    fo <- as.formula(paste("y ~ species +", paste(PREDS, collapse = " + ")))
    m <- lm(fo, data = cbind(X0[tr, , drop = FALSE], y = y[tr]))
    as.numeric(predict(m, newdata = X0[te, , drop = FALSE]))
  },
  lm_log = function(tr, te) {
    sh <- abs(min(y[tr])) + 0.01
    yl <- log(y[tr] + sh)
    fo <- as.formula(paste("yl ~ species +", paste(PREDS, collapse = " + ")))
    m <- lm(fo, data = cbind(X0[tr, , drop = FALSE], yl = yl))
    r <- residuals(m); smear <- mean(exp(r))          # Duan smearing, not exp(fit)
    as.numeric(exp(predict(m, newdata = X0[te, , drop = FALSE])) * smear - sh)
  })
if (HAVE_LME4) FIT$lmm <- function(tr, te) {
  fo <- as.formula(paste("y ~", paste(PREDS, collapse = " + "), "+ (1|species)"))
  m <- suppressMessages(lme4::lmer(fo, data = cbind(X0[tr, , drop = FALSE], y = y[tr])))
  as.numeric(predict(m, newdata = X0[te, , drop = FALSE], allow.new.levels = TRUE))
}

K <- 5; REP <- 3
PRED <- lapply(names(FIT), function(z) matrix(NA_real_, nrow(X0), REP))
names(PRED) <- names(FIT)
for (rep in seq_len(REP)) {
  set.seed(200 + rep)
  fold <- ave(seq_len(nrow(X0)), sp, FUN = function(i) sample(rep_len(1:K, length(i))))
  for (k in 1:K) {
    te <- which(fold == k); tr <- which(fold != k)
    for (nm in names(FIT))
      PRED[[nm]][te, rep] <- tryCatch(FIT[[nm]](tr, te), error = function(e) NA_real_)
  }
  cat(sprintf("  repeat %d done\n", rep))
}

TP <- read.csv("outputs/tables/tree_flux_predictions.csv", stringsAsFactors = FALSE)
TP <- TP[TP$in_stand, ]
area_by_sp <- TP %>% group_by(species) %>%
  summarise(A = sum(A_stem_m2), .groups = "drop")

rows <- lapply(names(FIT), function(nm) {
  p <- rowMeans(PRED[[nm]], na.rm = TRUE); ok <- is.finite(p)
  bysp <- data.frame(s = sp[ok], o = y[ok], p = p[ok]) %>% group_by(s) %>%
    summarise(obs = mean(o), pred = mean(p), n = dplyr::n(), .groups = "drop") %>%
    mutate(ratio = pred/obs)
  # area-weighted bias: how wrong is the SUM this model would give
  w <- bysp %>% left_join(area_by_sp, by = c("s" = "species")) %>%
       filter(is.finite(A))
  data.frame(model = nm, n = sum(ok),
    rmse = sqrt(mean((y[ok]-p[ok])^2)),
    r2   = 1 - sum((y[ok]-p[ok])^2)/sum((y[ok]-mean(y[ok]))^2),
    mean_abs_sp_bias = mean(abs(bysp$pred - bysp$obs)),
    max_abs_sp_ratio_dev = max(abs(log(pmax(bysp$ratio, 1e-6)))),
    area_wtd_ratio = sum(w$pred*w$A)/sum(w$obs*w$A))
})
R <- bind_rows(rows) %>% arrange(abs(area_wtd_ratio - 1))
cat("\n=== MODEL FAMILY COMPARISON (repeated 5-fold CV) ===\n")
cat("area_wtd_ratio is the one that matters for a budget: 1.00 = unbiased sum\n\n")
print(as.data.frame(R %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)

cat("\n=== per-species mean ratio (pred/obs), by model ===\n")
tab <- bind_rows(lapply(names(FIT), function(nm) {
  p <- rowMeans(PRED[[nm]], na.rm = TRUE); ok <- is.finite(p)
  data.frame(s = sp[ok], o = y[ok], p = p[ok], model = nm) }))
wide <- tab %>% group_by(model, s) %>% summarise(r = mean(p)/mean(o), n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = model, values_from = r)
print(as.data.frame(wide %>% arrange(-n) %>% mutate(across(where(is.numeric), ~round(.x, 2)))), row.names = FALSE)

dir.create("outputs/revision", showWarnings = FALSE, recursive = TRUE)
write.csv(R, "outputs/revision/model_family_comparison.csv", row.names = FALSE)
con <- file("outputs/revision/model_family_comparison.txt", "w")
cat(sprintf("MODEL FAMILY COMPARISON\nbuilt %s\nrepeated %d-fold CV x %d\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), K, REP), file = con)
utils::write.table(round(as.data.frame(R %>% select(-model)), 4), con, sep = "\t", quote = FALSE)
close(con)
cat("\nwritten: outputs/revision/model_family_comparison.{csv,txt}\n")
