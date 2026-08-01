#!/usr/bin/env Rscript
# ==============================================================================
# rev_model_family_comparison.R
# ------------------------------------------------------------------------------
# What model family is right for a stand BUDGET, as opposed to for predicting an
# individual stem?
#
# THE REQUIREMENT IS UNUSUAL. The model is fitted to ~1000 chamber measurements
# and then applied to ~8000 inventory stems at combinations of species, diameter,
# height, moisture and temperature that are mostly NOT represented in the
# training data -- and the quantity of interest is the SUM. So the model must
#   (a) capture whatever nonlinearity and interaction the data contain, and
#   (b) be unbiased in the mean over a population it was not fitted to.
# Those pull in opposite directions and most families give up one for the other.
#
# WHY THE FOREST FAILS (b). Squared-error minimisation shrinks predictions toward
# the conditional mean, so low-flux species are over-predicted and high-flux ones
# under-predicted. Because Tsuga and Pinus strobus are both low-flux and carry
# 54% of the band area, the stand sum comes out ~20% high.
#
# WHY A SCALAR PER-SPECIES CORRECTION IS NOT ENOUGH, which is the objection this
# version exists to test. Such a factor is estimated at the covariate
# distribution of the MEASUREMENTS and then applied to inventory stems that are
# systematically smaller, span all twelve months, and sit across the full
# moisture range. If the bias varies across covariate space, a single scalar per
# species transports it incorrectly. So this script also reports bias WITHIN
# quartiles of diameter, moisture and predicted flux -- if a correction is
# adequate its residual bias should be flat across all three.
#
# CANDIDATES
#   rf                  the locked forest, uncorrected
#   rf_cal_species      per-species multiplicative factor (training fold only)
#   rf_cal_linear       calibration on the PREDICTION, obs ~ a + b*pred, inverted;
#                       corrects shrinkage wherever it occurs rather than per taxon
#   rf_cal_isotonic     monotone nonparametric version of the same idea
#   gam                 mgcv, thin-plate smooths on the continuous predictors plus
#                       a species factor; penalised least squares, so residuals
#                       sum to ~zero and it is mean-unbiased LIKE lm, but it can
#                       still bend
#   lm                  least squares, the unbiased-but-blind baseline
#
# Output: outputs/data/model_family_comparison.csv / .txt
# ==============================================================================
suppressPackageStartupMessages({library(dplyr); library(ranger); library(mgcv)})
set.seed(42)
load("outputs/models/TRAINING_DATA.RData")
source("code/lib/rev_geometry.R")

X0 <- as.data.frame(X_tree)
PREDS <- setdiff(names(X0), "species")
d <- tree_train_complete
stopifnot(nrow(X0) == nrow(d))
keep <- is.finite(d$stem_flux_corrected) & complete.cases(X0[, PREDS, drop = FALSE])
X0 <- X0[keep, , drop = FALSE]; d <- d[keep, ]
y <- d$stem_flux_corrected
sp <- as.character(X0$species); X0$species <- factor(sp)
cat(sprintf("rows %d | species levels %d\n\n", nrow(X0), nlevels(X0$species)))

rf_fit <- function(tr) ranger(x = X0[tr, , drop = FALSE], y = y[tr], num.trees = 800,
  min.node.size = 5, mtry = max(1, floor(sqrt(ncol(X0)))), num.threads = 1,
  seed = 42, keep.inbag = FALSE)

FIT <- list(
  rf = function(tr, te) predict(rf_fit(tr), X0[te, , drop = FALSE])$predictions,

  rf_cal_species = function(tr, te) {
    m <- rf_fit(tr)
    cal <- data.frame(s = sp[tr], o = y[tr], p = m$predictions) %>% group_by(s) %>%
      summarise(r = ifelse(mean(p) > 0, mean(o)/mean(p), 1), .groups = "drop")
    cm <- setNames(cal$r, cal$s)
    p <- predict(m, X0[te, , drop = FALSE])$predictions
    p * pmax(pmin(as.numeric(ifelse(is.na(cm[sp[te]]), 1, cm[sp[te]])), 5), 0.2)
  },

  rf_cal_linear = function(tr, te) {
    m <- rf_fit(tr)
    cf <- coef(lm(o ~ p, data = data.frame(o = y[tr], p = m$predictions)))
    as.numeric(cf[1] + cf[2] * predict(m, X0[te, , drop = FALSE])$predictions)
  },

  rf_cal_isotonic = function(tr, te) {
    m <- rf_fit(tr)
    o <- order(m$predictions)
    ir <- isoreg(m$predictions[o], y[tr][o])
    fk <- approxfun(ir$x, ir$yf, rule = 2)
    as.numeric(fk(predict(m, X0[te, , drop = FALSE])$predictions))
  },

  gam = function(tr, te) {
    cn <- PREDS[sapply(X0[, PREDS, drop = FALSE], function(z) length(unique(z)) > 9)]
    fo <- as.formula(paste("y ~ species +",
      paste(c(sprintf("s(%s, k=5)", cn), setdiff(PREDS, cn)), collapse = " + ")))
    m <- try(gam(fo, data = cbind(X0[tr, , drop = FALSE], y = y[tr]), method = "REML"), silent = TRUE)
    if (inherits(m, "try-error")) return(rep(NA_real_, length(te)))
    as.numeric(predict(m, newdata = X0[te, , drop = FALSE]))
  },

  lm = function(tr, te) {
    fo <- as.formula(paste("y ~ species +", paste(PREDS, collapse = " + ")))
    as.numeric(predict(lm(fo, data = cbind(X0[tr, , drop = FALSE], y = y[tr])),
                       newdata = X0[te, , drop = FALSE]))
  })

K <- 5; REP <- 3
PRED <- setNames(lapply(names(FIT), function(z) matrix(NA_real_, nrow(X0), REP)), names(FIT))
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
area_by_sp <- TP[TP$in_stand, ] %>% group_by(species) %>%
  summarise(A = sum(A_stem_m2), .groups = "drop")

qb <- function(v, n = 4) cut(v, unique(quantile(v, seq(0, 1, length.out = n + 1))),
                             include.lowest = TRUE, labels = FALSE)
q_dbh <- qb(X0$dbh_m); q_moi <- qb(X0$soil_moisture_at_tree)

summ <- bind_rows(lapply(names(FIT), function(nm) {
  p <- rowMeans(PRED[[nm]], na.rm = TRUE); ok <- is.finite(p)
  bysp <- data.frame(s = sp[ok], o = y[ok], p = p[ok]) %>% group_by(s) %>%
    summarise(obs = mean(o), pred = mean(p), .groups = "drop")
  w <- bysp %>% left_join(area_by_sp, by = c("s" = "species")) %>% filter(is.finite(A))
  q_pred <- qb(p[ok])
  rq <- function(g) { r <- tapply(seq_along(p[ok]), g, function(i)
                        mean(p[ok][i])/mean(y[ok][i])); max(abs(log(pmax(r, 1e-6)))) }
  data.frame(model = nm, rmse = sqrt(mean((y[ok]-p[ok])^2)),
    r2 = 1 - sum((y[ok]-p[ok])^2)/sum((y[ok]-mean(y[ok]))^2),
    area_wtd_ratio = sum(w$pred*w$A)/sum(w$obs*w$A),
    maxdev_by_dbh = rq(q_dbh[ok]), maxdev_by_moist = rq(q_moi[ok]),
    maxdev_by_pred = rq(q_pred))
}))
R <- summ %>% arrange(abs(area_wtd_ratio - 1))
cat("\n=== MODEL FAMILY COMPARISON (repeated 5-fold CV) ===\n")
cat("area_wtd_ratio: 1.00 = unbiased SUM.  maxdev_*: worst |log ratio| across\n")
cat("quartiles of that variable -- 0 means the correction transports correctly.\n\n")
print(as.data.frame(R %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)

cat("\n=== residual bias by quartile (pred/obs), the transportability test ===\n")
for (nm in c("rf", "rf_cal_species", "rf_cal_linear", "rf_cal_isotonic", "gam")) {
  p <- rowMeans(PRED[[nm]], na.rm = TRUE); ok <- is.finite(p)
  rr <- function(g) round(tapply(seq_along(p[ok]), g, function(i)
          mean(p[ok][i])/mean(y[ok][i])), 2)
  cat(sprintf("  %-16s dbh Q1-Q4 %s | moisture Q1-Q4 %s\n", nm,
              paste(rr(q_dbh[ok]), collapse = " "), paste(rr(q_moi[ok]), collapse = " ")))
}

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(R, "outputs/data/model_family_comparison.csv", row.names = FALSE)
con <- file("outputs/audit/model_family_comparison.txt", "w")
cat(sprintf("MODEL FAMILY COMPARISON\nbuilt %s\nrepeated %d-fold CV x %d\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), K, REP), file = con)
utils::write.table(round(as.data.frame(R %>% select(-model)), 4), con, sep = "\t", quote = FALSE)
close(con)
cat("\nwritten: outputs/model_family_comparison.{csv,txt}\n")
