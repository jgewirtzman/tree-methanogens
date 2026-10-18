#!/usr/bin/env Rscript
# ==============================================================================
# 09_rf_species_fallback_loso.R -- predicting a species the model never learned
# ------------------------------------------------------------------------------
# Stage D; run from the repository root. Leave-one-species-out: hold out every
# measurement of one species, fit on the rest, predict the held-out species -- the
# situation an unmeasured inventory species presents. Three schemes:
#
#   pooled        the held-out species goes to SPECIES_OTHER
#   genus         to a genus-level factor when the genus is measured, else SPECIES_OTHER
#   species_free  species dropped from the predictors; size and environment only
#
# Output: outputs/data/species_fallback_loso.csv / .txt
# ==============================================================================
suppressPackageStartupMessages({library(ranger); library(dplyr)})
source("code/lib/run_mode.R")   # RUN: run sizes
set.seed(42)
load("outputs/models/TRAINING_DATA.RData")

# Use the model's OWN feature matrix rather than rebuilding it: X_tree already
# carries the locked 6-predictor set, including the imputation of missing
# measurement height to 125 cm (839 of 1130 rows, i.e. 74% of training data is at
# an assumed breast height -- worth knowing, but not this script's business).
X0 <- as.data.frame(X_tree)
PREDS <- setdiff(names(X0), "species")
d <- tree_train_complete
stopifnot(nrow(X0) == nrow(d))
RESP <- "stem_flux_corrected"            # raw nmol m-2 s-1, matching the locked model
keep <- is.finite(d[[RESP]]) & complete.cases(X0[, PREDS, drop = FALSE])
X0 <- X0[keep, , drop = FALSE]; d <- d[keep, ]
d$sp    <- as.character(d$species_clean)
d$genus <- sub(" .*", "", as.character(d$species))
y <- d[[RESP]]
cat(sprintf("training rows usable: %d | species levels: %d | genera: %d\n",
            nrow(d), length(unique(d$sp)), length(unique(d$genus))))
cat(sprintf("predictors: %s\n", paste(PREDS, collapse = ", ")))

# Hold out only species with enough measurements for the held-out R2 to mean
# something; everything else stays in the training pool.
HOLD <- names(which(table(d$species)[table(d$species) >= 20] > 0))
HOLD <- HOLD[!is.na(HOLD) & HOLD != ""]
cat(sprintf("species held out (n >= 20): %d -- %s\n\n", length(HOLD), paste(HOLD, collapse = ", ")))

fitpred <- function(scheme, tr, te) {
  Xtr <- X0[tr, PREDS, drop = FALSE]; Xte <- X0[te, PREDS, drop = FALSE]
  if (scheme != "species_free") {
    if (scheme == "pooled") {
      lab_tr <- d$sp[tr]
      lab_te <- rep("SPECIES_OTHER", length(te))          # unseen -> pooled bucket
    } else {                                             # genus fallback
      lab_tr <- d$genus[tr]
      seen <- unique(lab_tr)
      lab_te <- ifelse(d$genus[te] %in% seen, d$genus[te], "SPECIES_OTHER")
    }
    lv <- sort(unique(c(lab_tr, lab_te, "SPECIES_OTHER")))
    Xtr$species <- factor(lab_tr, levels = lv)
    Xte$species <- factor(lab_te, levels = lv)
  }
  m <- ranger(x = Xtr, y = y[tr], num.trees = RUN$refit_num_trees, min.node.size = 5,
              mtry = max(1, floor(sqrt(ncol(Xtr)))), num.threads = 1, seed = 42)
  predict(m, Xte)$predictions
}

SCHEMES <- c("pooled", "genus", "species_free")
res <- list()
for (sc in SCHEMES) {
  pr <- rep(NA_real_, nrow(d))
  for (s in HOLD) {
    te <- which(d$species == s); tr <- which(d$species != s)
    pr[te] <- fitpred(sc, tr, te)
  }
  ok <- !is.na(pr)
  res[[sc]] <- data.frame(scheme = sc, n = sum(ok),
    r2   = 1 - sum((y[ok]-pr[ok])^2)/sum((y[ok]-mean(y[ok]))^2),
    rmse = sqrt(mean((y[ok]-pr[ok])^2)),
    rho  = cor(y[ok], pr[ok], method = "spearman"),
    bias = mean(pr[ok]) - mean(y[ok]),
    ratio = mean(pr[ok])/mean(y[ok]))
  # per-species, so one species cannot carry the verdict
  ps <- bind_rows(lapply(HOLD, function(s) { i <- which(d$species == s)
    data.frame(scheme = sc, species = s, n = length(i),
               bias = mean(pr[i]) - mean(y[i]),
               rmse = sqrt(mean((y[i]-pr[i])^2))) }))
  res[[paste0(sc, "_by_sp")]] <- ps
}
OVER <- bind_rows(res[SCHEMES]) %>% arrange(-r2)
cat("=== LEAVE-ONE-SPECIES-OUT: predicting a species never measured ===\n")
print(as.data.frame(OVER %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)

BY <- bind_rows(res[paste0(SCHEMES, "_by_sp")])
cat("\n=== per held-out species, RMSE ===\n")
print(as.data.frame(BY %>% select(scheme, species, n, rmse) %>%
      tidyr::pivot_wider(names_from = scheme, values_from = rmse) %>%
      mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat("\nwins by species (lowest RMSE):\n")
print(table(BY %>% group_by(species) %>% slice_min(rmse, n = 1) %>% pull(scheme)))

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(BY, "outputs/data/species_fallback_loso.csv", row.names = FALSE)
con <- file("outputs/audit/species_fallback_loso.txt", "w")
cat(sprintf("SPECIES FALLBACK, LEAVE-ONE-SPECIES-OUT\nbuilt %s\nheld out: %s\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(HOLD, collapse = ", ")), file = con)
utils::write.table(round(as.data.frame(OVER %>% select(-scheme)), 4), con, sep = "\t", quote = FALSE)
close(con)
cat("\nwritten: outputs/species_fallback_loso.{csv,txt}\n")
