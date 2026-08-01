#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_species_fallback_loso.R
# ------------------------------------------------------------------------------
# How should the model predict a stem whose species it never learned?
#
# THE PROBLEM. species_clean lumps any species with fewer than 10 measurements
# into a single SPECIES_OTHER level. In this inventory that level receives 1,971
# stems -- 24.6% of the count -- but only 3.7% of the 0-2 m band area, because
# 1,898 of them are Kalmia latifolia, an understory shrub. The bucket behind that
# level holds 19 measurements: Quercus velutina 6, Carya ovata 3, Kalmia 3,
# Sassafras albidum 3, Betula populifolia 2, Carya laciniosa 1, unknown 1. So
# 1,898 mountain-laurel shrubs inherit a prediction that is 84% canopy hardwood,
# and they come out 44% HIGHER than trained species (median 0.280 vs 0.195
# nmol m-2 s-1) -- the reddest patches on the Figure 9 map.
#
# It also is not uniformly trivial in area: Quercus velutina is 19 stems but
# 2.52 m2 of basal area, 1.5% of the stand, and it IS measured (it is the felled
# black oak) -- just below the threshold that would give it its own level.
#
# THREE CANDIDATE SCHEMES, compared by leave-one-species-out: hold out every
# measurement of one species, fit on the rest, predict the held-out species. That
# is exactly the situation an unmeasured inventory species presents.
#
#   pooled        assign the held-out species to SPECIES_OTHER (current behaviour)
#   genus         assign it to a genus-level factor when the genus is measured,
#                 else SPECIES_OTHER
#   species_free  drop species from the predictor set entirely and predict from
#                 size and environment alone
#
# The species-free scheme is the honest null: for a species you never measured,
# predict from what you did measure about it -- diameter, height, moisture,
# temperature -- rather than borrowing a congener's identity.
#
# The earlier rev_rf_leave_one_species_out.R is superseded: it used the asinh
# response and feature matrix from before the model lock, and read a traits table
# from an absolute path outside this repository, so it cannot reproduce.
#
# Output: outputs/data/species_fallback_loso.csv / .txt
# ==============================================================================
suppressPackageStartupMessages({library(ranger); library(dplyr)})
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
  m <- ranger(x = Xtr, y = y[tr], num.trees = 800, min.node.size = 5,
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
