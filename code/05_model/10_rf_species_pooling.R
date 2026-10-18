#!/usr/bin/env Rscript
# ==============================================================================
# 10_rf_species_pooling.R -- how much should a rarely measured species be pooled?
# ------------------------------------------------------------------------------
# Stage D; run from the repository root. A bias-variance question, distinct from the
# leave-one-species-out test (09): here the species has records, just few, and the
# question is whether to trust them, pool them away, or shrink between the two (e.g.
# Kalmia latifolia: 3 records, carried to ~1,900 inventory stems). The pooling
# schemes are compared by repeated stratified k-fold cross-validation over
# observations, with the rare-species rows reported separately (an overall RMSE is
# dominated by the abundant species and always favours pooling):
#   overall RMSE / R2       the dataset as a whole
#   rare-species RMSE       the species the question is about
#   rare-species BIAS       signed
#   stand tree budget       the effect on the number in the paper
# The taxonomic ladder of 02_rf_models.R ("genus_then_pool") is the scheme adopted.
#
# Output: outputs/data/species_pooling_cv.csv / .txt
# ==============================================================================
suppressPackageStartupMessages({library(ranger); library(dplyr)})
set.seed(42)
load("outputs/models/TRAINING_DATA.RData")

X0 <- as.data.frame(X_tree)
PREDS <- setdiff(names(X0), "species")
d <- tree_train_complete
stopifnot(nrow(X0) == nrow(d))
keep <- is.finite(d$stem_flux_corrected) & complete.cases(X0[, PREDS, drop = FALSE])
X0 <- X0[keep, , drop = FALSE]; d <- d[keep, ]
y <- d$stem_flux_corrected
sp_raw <- as.character(d$species)
unident <- is.na(sp_raw) | grepl("^unknown", sp_raw, ignore.case = TRUE) | !grepl(" ", trimws(sp_raw))
sp_raw[unident] <- "SPECIES_OTHER"
genus <- ifelse(unident, "SPECIES_OTHER", sub(" .*", "", sp_raw))
n_by_sp <- table(sp_raw)
RARE_CUT <- 10
rare_row <- !unident & n_by_sp[sp_raw] < RARE_CUT
cat(sprintf("rows %d | species %d | rare (<%d records) rows %d across %d species\n",
            nrow(d), length(unique(sp_raw)), RARE_CUT, sum(rare_row),
            length(unique(sp_raw[rare_row]))))
cat(sprintf("rare species: %s\n\n", paste(sort(unique(sp_raw[rare_row])), collapse = ", ")))

# --- label schemes ------------------------------------------------------------
# Each returns the species label to hand the forest, given the TRAINING rows only
# (so nothing leaks from the held-out fold).
lab_lump10 <- function(tr, idx) {
  cnt <- table(sp_raw[tr]); big <- names(cnt)[cnt >= RARE_CUT]
  ifelse(sp_raw[idx] %in% big, sp_raw[idx], "SPECIES_OTHER")
}
lab_all <- function(tr, idx) {
  seen <- unique(sp_raw[tr])
  ifelse(sp_raw[idx] %in% seen, sp_raw[idx], "SPECIES_OTHER")
}
lab_genus_then_pool <- function(tr, idx) {
  cnt <- table(sp_raw[tr]); big <- names(cnt)[cnt >= RARE_CUT]
  gseen <- unique(genus[tr])
  ifelse(sp_raw[idx] %in% big, sp_raw[idx],
    ifelse(genus[idx] %in% gseen, paste0("GEN_", genus[idx]), "SPECIES_OTHER"))
}

SCHEMES <- list(
  lump10             = list(lab = lab_lump10,          node = 5),
  all_species        = list(lab = lab_all,             node = 5),
  all_species_node2  = list(lab = lab_all,             node = 2),
  genus_then_pool    = list(lab = lab_genus_then_pool, node = 5),
  shrunk_offset      = list(lab = NULL,                node = 5))

# --- shrunk offset: partial pooling, species toward genus toward grand mean ----
# Empirical-Bayes shrinkage with a fixed prior weight K: a species with n records
# keeps n/(n+K) of its own deviation and borrows the rest from its genus, and the
# genus does the same toward the grand mean. Fed to the forest as a NUMERIC
# predictor, so the forest can use it or ignore it.
shrink_offset <- function(tr, idx, K = 5) {
  gm <- mean(y[tr])
  gdf <- data.frame(g = genus[tr], r = y[tr] - gm) %>% group_by(g) %>%
         summarise(n = dplyr::n(), m = mean(r), .groups = "drop") %>%
         mutate(off = m * n/(n + K))
  gmap <- setNames(gdf$off, gdf$g)
  gof_tr <- ifelse(is.na(gmap[genus[tr]]), 0, gmap[genus[tr]])
  sdf <- data.frame(s = sp_raw[tr], r = y[tr] - gm - gof_tr) %>% group_by(s) %>%
         summarise(n = dplyr::n(), m = mean(r), .groups = "drop") %>%
         mutate(off = m * n/(n + K))
  smap <- setNames(sdf$off, sdf$s)
  g <- ifelse(is.na(gmap[genus[idx]]), 0, gmap[genus[idx]])
  s <- ifelse(is.na(smap[sp_raw[idx]]), 0, smap[sp_raw[idx]])
  as.numeric(g + s)
}

fitpred <- function(nm, tr, te) {
  cfg <- SCHEMES[[nm]]
  Xtr <- X0[tr, PREDS, drop = FALSE]; Xte <- X0[te, PREDS, drop = FALSE]
  if (nm == "shrunk_offset") {
    Xtr$taxon_offset <- shrink_offset(tr, tr); Xte$taxon_offset <- shrink_offset(tr, te)
  } else {
    a <- cfg$lab(tr, tr); b <- cfg$lab(tr, te)
    lv <- sort(unique(c(a, b, "SPECIES_OTHER")))
    Xtr$species <- factor(a, levels = lv); Xte$species <- factor(b, levels = lv)
  }
  m <- ranger(x = Xtr, y = y[tr], num.trees = 800, min.node.size = cfg$node,
              mtry = max(1, floor(sqrt(ncol(Xtr)))), num.threads = 1, seed = 42)
  predict(m, Xte)$predictions
}

# --- repeated stratified k-fold ----------------------------------------------
K_FOLD <- 5; N_REP <- 4
res <- list()
for (nm in names(SCHEMES)) {
  pr <- matrix(NA_real_, nrow(d), N_REP)
  for (rep in seq_len(N_REP)) {
    set.seed(100 + rep)
    fold <- ave(seq_len(nrow(d)), sp_raw, FUN = function(i) sample(rep_len(1:K_FOLD, length(i))))
    for (k in 1:K_FOLD) {
      te <- which(fold == k); tr <- which(fold != k)
      if (!length(te)) next
      pr[te, rep] <- fitpred(nm, tr, te)
    }
  }
  p <- rowMeans(pr, na.rm = TRUE); ok <- is.finite(p)
  rr <- ok & rare_row
  res[[nm]] <- data.frame(scheme = nm,
    rmse_all  = sqrt(mean((y[ok]-p[ok])^2)),
    r2_all    = 1 - sum((y[ok]-p[ok])^2)/sum((y[ok]-mean(y[ok]))^2),
    rmse_rare = sqrt(mean((y[rr]-p[rr])^2)),
    bias_rare = mean(p[rr]) - mean(y[rr]),
    ratio_rare = mean(p[rr])/mean(y[rr]),
    n_rare = sum(rr))
  cat(sprintf("  %-18s done\n", nm))
}
R <- bind_rows(res)
cat("\n=== REPEATED 5-FOLD CV (4 repeats), rare = species with < 10 records ===\n")
print(as.data.frame(R %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
cat("\nranking by rare-species RMSE:  ", paste(R$scheme[order(R$rmse_rare)], collapse = " < "), "\n")
cat("ranking by |rare-species bias|: ", paste(R$scheme[order(abs(R$bias_rare))], collapse = " < "), "\n")
cat("ranking by overall RMSE:        ", paste(R$scheme[order(R$rmse_all)], collapse = " < "), "\n")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(R, "outputs/data/species_pooling_cv.csv", row.names = FALSE)
con <- file("outputs/audit/species_pooling_cv.txt", "w")
cat(sprintf("SPECIES POOLING, REPEATED %d-FOLD CV x %d\nbuilt %s\nrare cut: < %d records\n\n",
            K_FOLD, N_REP, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), RARE_CUT), file = con)
utils::write.table(round(as.data.frame(R %>% select(-scheme)), 4), con, sep = "\t", quote = FALSE)
close(con)
cat("\nwritten: outputs/species_pooling_cv.{csv,txt}\n")
