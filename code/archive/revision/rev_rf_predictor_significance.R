#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_predictor_significance.R
# ------------------------------------------------------------------------------
# Three questions, answered empirically:
#
#  (1) DBH ENCODING. The submitted model uses `dbh_within_z` -- DBH standardized
#      within species -- instead of raw DBH. That is a linear-model device: it
#      deconfounds size from species when the model cannot form species x size
#      interactions on its own. A random forest can: given `species` as a split
#      variable, the tree splits on species first and then on DBH, recovering any
#      species-specific size response directly. Because the z-score is affine
#      WITHIN each species, and forests are invariant to monotone transforms of a
#      single predictor, the two encodings are information-equivalent once species
#      is in the model -- except that the z-score DESTROYS absolute size (a 10 cm
#      sapling and a 90 cm oak can both sit at z = 0). We test raw / z / both.
#
#  (2) PREDICTOR SIGNIFICANCE. Three complementary tests, none of which relies on
#      the raw importance ranking (which has no null distribution):
#        a. SHADOW FEATURES (the Boruta idea, and the generalization of the
#           "add a random column" heuristic). A single uniform random column is a
#           weak null: it has no relationship to anything and is easy to beat.
#           The stronger null is a PERMUTED COPY of each real predictor -- same
#           marginal distribution, same cardinality, same number of split points,
#           relationship to y destroyed. A predictor is confirmed only if it beats
#           the best shadow, repeatedly.
#        b. ALTMANN P-VALUES (ranger, built in). Permute y, refit, recompute
#           importance -> an exact null distribution of importance per predictor.
#           Gives a real per-predictor p-value.
#        c. DROP-COLUMN dOOB R2. Refit without each predictor; the loss in OOB R2
#           is the effect size a reader actually cares about. Repeated over seeds
#           so the noise band is visible.
#      (Janitza's fast p-value is deliberately NOT used: it is only valid in the
#      high-dimensional regime, and we have ~9 predictors.)
#
#  (3) REPORTING ORDER. OOB R2 is the headline (it is what the submitted paper
#      reported and what reviewers accepted); grouped CV is reported secondarily
#      as the stricter, upscaling-relevant check.
#
# Outputs: outputs/revision/dbh_encoding_comparison.csv
#          outputs/revision/predictor_significance.csv
#          outputs/revision/predictor_significance.png
#          outputs/revision/dbh_encoding.png
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(ggplot2); library(readxl)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cache <- file.path(outdir, "predictor_frame_cache.rds")

# ------------------------------------------------------------------------------
# Assemble the candidate frame (identical derivation to rev_rf_model_selection.R)
# ------------------------------------------------------------------------------
if (file.exists(cache)) {
  CA <- readRDS(cache); TREE <- CA$TREE; yraw <- CA$yraw; grp <- CA$grp
  cat("loaded cached predictor frame:", nrow(TREE), "rows\n")
} else {
  load("outputs/models/TRAINING_DATA.RData")
  load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
  INV <- rf_workflow_data$PLACEHOLDER_INVENTORY
  MLAT <- 111320; MLON <- 111320 * cos(41.99 * pi / 180)

  fix_dbh <- function(dbh, species) {
    d <- ifelse(!is.na(dbh) & dbh > 3, dbh / 100, dbh)
    d <- ifelse(grepl("Betula", species) & !is.na(d) & d * 100 > 200, d / 10, d)
    d <- ifelse(species == "Pinus strobus" & !is.na(d) & d * 100 > 230, d / 10, d)
    d <- ifelse(species == "Kalmia latifolia" & !is.na(d) & d * 100 > 100, d / 100,
         ifelse(species == "Kalmia latifolia" & !is.na(d) & d * 100 > 10,  d / 10, d))
    d
  }
  INV$dbh_fix <- fix_dbh(INV$dbh_m, INV$species)

  riv <- suppressMessages(read_excel("data/raw/inventory/spatial_data/River.xlsx"))
  names(riv) <- make.names(names(riv))
  rok <- !is.na(riv$Latitude) & !is.na(riv$Longitude)
  rx <- riv$Longitude[rok] * MLON; ry <- riv$Latitude[rok] * MLAT
  dist_river <- function(x, y) vapply(seq_along(x), function(i) {
    if (is.na(x[i]) || is.na(y[i])) return(NA_real_)
    min(sqrt((x[i] * MLON - rx)^2 + (y[i] * MLAT - ry)^2)) }, numeric(1))

  ELEV <- read.csv("outputs/tables/combined_elevation_data.csv", stringsAsFactors = FALSE)
  ex <- ELEV$Longitude * MLON; ey <- ELEV$Latitude * MLAT; ez <- ELEV$Elevation
  elev_at <- function(x, y, k = 6) vapply(seq_along(x), function(i) {
    if (is.na(x[i]) || is.na(y[i])) return(NA_real_)
    dd <- sqrt((x[i] * MLON - ex)^2 + (y[i] * MLAT - ey)^2)
    o <- order(dd)[seq_len(min(k, length(dd)))]
    w <- 1 / pmax(dd[o], 0.5)^2
    sum(w * ez[o]) / sum(w) }, numeric(1))

  ix <- INV$x * MLON; iy <- INV$y * MLAT; iba <- pi * (INV$dbh_fix / 2)^2
  local_ba <- function(x, y) vapply(seq_along(x), function(i) {
    if (is.na(x[i]) || is.na(y[i])) return(NA_real_)
    dd <- (x[i] * MLON - ix)^2 + (y[i] * MLAT - iy)^2
    sum(iba[dd <= 100], na.rm = TRUE) }, numeric(1))

  d <- tree_train_complete
  yraw <- d$stem_flux_corrected
  keep <- is.finite(yraw)
  d <- d[keep, ]; yraw <- yraw[keep]

  cat("deriving spatial covariates for", nrow(d), "tree measurements...\n")
  sp <- factor(d$species_clean)
  # within-species z, computed exactly as 02_rf_models.R does it
  dbh_within_z <- ave(d$dbh_m, as.character(sp), FUN = function(x)
    if (length(x) >= 5 && sd(x, na.rm = TRUE) > 0)
      (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) else rep(0, length(x)))

  TREE <- data.frame(
    species      = sp,
    dbh_m        = d$dbh_m,
    dbh_within_z = dbh_within_z,
    moisture     = d$soil_moisture_at_tree,
    soil_temp    = d$soil_temp_C_mean,
    air_temp     = d$air_temp_C_mean,
    month        = d$month,
    height_cm    = ifelse(is.na(d$measurement_height_cm), 125, d$measurement_height_cm),
    dist_river_m = dist_river(d$x, d$y),
    local_ba_m2  = local_ba(d$x, d$y),
    elevation_m  = elev_at(d$x, d$y)
  )
  grp <- d$tree_id
  saveRDS(list(TREE = TREE, yraw = yraw, grp = grp), cache)
}

r2 <- function(a, b) 1 - sum((a - b)^2) / sum((a - mean(a))^2)

# ------------------------------------------------------------------------------
# TWO CONFOUNDS THAT MUST BE CONTROLLED, or every comparison below is invalid:
#
#  (i) mtry. The usual default mtry = floor(sqrt(n_predictors)) CHANGES when you
#      add or drop a column (9 -> 3, 8 -> 2). A drop-column test using the default
#      therefore varies the REGULARIZATION at the same time as the variable set,
#      and can easily report that dropping a useful predictor *improves* OOB R2.
#      We hold mtry FIXED at MTRY for every comparison.
# (ii) rows. Complete-case filtering on different predictor subsets changes n.
#      We fix ONE row set, defined by the full candidate frame, and use it
#      everywhere so every number below is computed on identical data.
# ------------------------------------------------------------------------------
MTRY  <- 3
NTREE <- 1000

fit <- function(X, y, seed = 42, ntree = NTREE, imp = "none", mtry = MTRY) {
  X <- as.data.frame(X)
  ranger(x = X, y = y, num.trees = ntree, min.node.size = 5,
         mtry = min(mtry, ncol(X)), num.threads = 1,
         importance = imp, seed = seed)
}

grouped_r2 <- function(X, y, grp, seed = 42) {
  X <- as.data.frame(X)
  ut <- unique(grp); set.seed(seed)
  gm <- setNames(sample(rep(1:5, length.out = length(ut))), ut); f <- unname(gm[grp])
  pr <- rep(NA_real_, nrow(X))
  for (k in 1:5) {
    tr <- which(f != k); te <- which(f == k)
    m <- fit(X[tr, , drop = FALSE], y[tr], seed = seed)
    pr[te] <- predict(m, X[te, , drop = FALSE])$predictions
  }
  c(grouped_r2 = r2(y, pr), mean_ratio = mean(pr) / mean(y),
    spearman = cor(y, pr, method = "spearman"))
}

# The best specification from rev_rf_model_selection.R (config F), minus DBH,
# which is the variable under test in part 1.
F_BASE <- c("species", "moisture", "soil_temp", "month", "height_cm",
            "dist_river_m", "local_ba_m2", "elevation_m")
ALLV   <- c(F_BASE, "dbh_m", "dbh_within_z")

# NOTE ON MISSING DATA: do NOT complete-case filter. The project decision
# (REVISION_PLAN.md D0) is "no imputation of environmental covariates; NAs passed to
# ranger", and ranger splits on missingness natively. Complete-casing the full
# candidate set here costs more than half the data (1,104 -> 503 measurements,
# 224 trees) and REVERSES the part-1 ranking -- an artifact of which rows survive,
# not of the encodings. All rows are kept; every configuration sees the same 1,104.
cat(sprintf("row set: n = %d measurements, %d trees (NAs passed to ranger, no imputation)\n",
            nrow(TREE), length(unique(grp))))
cat(sprintf("mtry held fixed at %d, num.trees = %d, for every model below\n\n", MTRY, NTREE))

# ==============================================================================
# PART 1 -- does within-species DBH standardization help a random forest?
# ==============================================================================
cat("\n", strrep("=", 84), "\n", sep = "")
cat("  PART 1: DBH ENCODING  (5 seeds; OOB R2 is the headline)\n")
cat(strrep("=", 84), "\n\n", sep = "")

DBH_VARIANTS <- list(
  "raw DBH"          = c(F_BASE, "dbh_m"),
  "within-species z" = c(F_BASE, "dbh_within_z"),
  "both"             = c(F_BASE, "dbh_m", "dbh_within_z"),
  "no DBH at all"    = F_BASE
)

dbh_res <- bind_rows(lapply(names(DBH_VARIANTS), function(lab) {
  X <- TREE[, DBH_VARIANTS[[lab]], drop = FALSE]
  oob <- vapply(1:5, function(s) fit(X, yraw, seed = s)$r.squared, numeric(1))
  g <- grouped_r2(X, yraw, grp)
  data.frame(encoding = lab, n_pred = ncol(X),
             oob_r2 = mean(oob), oob_sd = sd(oob),
             grouped_r2 = g[["grouped_r2"]], mean_ratio = g[["mean_ratio"]],
             spearman = g[["spearman"]])
}))
print(as.data.frame(dbh_res), row.names = FALSE, digits = 3)
write.csv(dbh_res, file.path(outdir, "dbh_encoding_comparison.csv"), row.names = FALSE)

# How much absolute size does the z-score throw away?
cat("\n  Within-species z destroys absolute size. Trees at |z| < 0.25 span DBH:\n")
near0 <- abs(TREE$dbh_within_z) < 0.25 & !is.na(TREE$dbh_m)
cat(sprintf("    %.1f to %.1f cm  (n = %d stems at effectively the same z)\n",
            100 * min(TREE$dbh_m[near0]), 100 * max(TREE$dbh_m[near0]), sum(near0)))
cat(sprintf("    correlation(raw DBH, within-species z) = %.3f\n",
            cor(TREE$dbh_m, TREE$dbh_within_z, use = "complete.obs")))

# ==============================================================================
# PART 2 -- which predictors are significant?
# ==============================================================================
BEST <- c(F_BASE, "dbh_m")
X  <- TREE[, BEST, drop = FALSE]
yb <- yraw                            # same rows as part 1 (confound ii)

cat("\n", strrep("=", 84), "\n", sep = "")
cat("  PART 2a: SHADOW FEATURES  (Boruta-style; 30 runs)\n")
cat(strrep("=", 84), "\n\n", sep = "")

NRUN <- 30
hits <- setNames(numeric(length(BEST)), BEST)
noise_hits <- setNames(numeric(length(BEST)), BEST)   # vs a plain uniform column
for (s in 1:NRUN) {
  set.seed(1000 + s)
  SH <- X
  # a permuted copy of every real predictor: same distribution, no signal
  for (v in BEST) SH[[paste0("shadow_", v)]] <- sample(X[[v]])
  SH$noise_uniform <- runif(nrow(X))                  # the naive single-random-column test
  rf <- fit(SH, yb, seed = 1000 + s, imp = "permutation")
  imp <- rf$variable.importance
  shadow_max <- max(imp[grep("^shadow_", names(imp))])
  for (v in BEST) {
    if (imp[[v]] > shadow_max)            hits[v] <- hits[v] + 1
    if (imp[[v]] > imp[["noise_uniform"]]) noise_hits[v] <- noise_hits[v] + 1
  }
}
shadow_tab <- data.frame(
  predictor = BEST,
  beats_best_shadow = sprintf("%d/%d", hits[BEST], NRUN),
  beats_uniform_noise = sprintf("%d/%d", noise_hits[BEST], NRUN),
  shadow_p = vapply(BEST, function(v)
    binom.test(hits[[v]], NRUN, 0.5, alternative = "greater")$p.value, numeric(1)),
  row.names = NULL
)
shadow_tab$verdict <- ifelse(shadow_tab$shadow_p < 0.01, "confirmed",
                      ifelse(shadow_tab$shadow_p < 0.05, "tentative", "rejected"))
print(shadow_tab, row.names = FALSE, digits = 3)

cat("\n", strrep("=", 84), "\n", sep = "")
cat("  PART 2b: ALTMANN PERMUTATION P-VALUES  (200 y-permutations)\n")
cat(strrep("=", 84), "\n\n", sep = "")

DF <- cbind(y = yb, X)
rf_form <- ranger(y ~ ., data = DF, num.trees = 500, min.node.size = 5,
                  mtry = max(1, floor(sqrt(ncol(X)))), num.threads = 1,
                  importance = "permutation", seed = 42)
alt <- importance_pvalues(rf_form, method = "altmann", formula = y ~ ., data = DF,
                          num.permutations = 200)
alt_tab <- data.frame(predictor = rownames(alt),
                      importance = alt[, "importance"],
                      altmann_p = alt[, "pvalue"], row.names = NULL)
print(alt_tab[order(-alt_tab$importance), ], row.names = FALSE, digits = 3)

cat("\n", strrep("=", 84), "\n", sep = "")
cat("  PART 2c: DROP-COLUMN dOOB R2  (5 seeds each)\n")
cat(strrep("=", 84), "\n\n", sep = "")

full_oob <- vapply(1:5, function(s) fit(X, yb, seed = s)$r.squared, numeric(1))
drop_tab <- bind_rows(lapply(BEST, function(v) {
  o <- vapply(1:5, function(s) fit(X[, setdiff(BEST, v), drop = FALSE], yb, seed = s)$r.squared,
              numeric(1))
  data.frame(predictor = v, oob_without = mean(o),
             delta_oob = mean(full_oob) - mean(o),
             delta_sd = sqrt(var(full_oob) / 5 + var(o) / 5))
}))
drop_tab$delta_over_sd <- drop_tab$delta_oob / drop_tab$delta_sd
cat(sprintf("full model OOB R2 = %.4f (sd %.4f over 5 seeds)\n\n", mean(full_oob), sd(full_oob)))
print(as.data.frame(drop_tab[order(-drop_tab$delta_oob), ]), row.names = FALSE, digits = 3)

# ==============================================================================
# PART 2d -- BLOCK drops, because single-variable tests lie under collinearity
# ------------------------------------------------------------------------------
# soil_temp is a MONTHLY PLOT-LEVEL CONSTANT, so it is very nearly a deterministic
# function of `month`. Permutation/shadow tests split credit between such a pair and
# can reject BOTH ("if I permute soil_temp, month rescues the prediction; if I permute
# month, soil_temp rescues it"). Single-variable drops have the same blind spot.
# Dropping the whole block is the only way to see what the block is worth.
# ==============================================================================
cat("\n", strrep("=", 84), "\n", sep = "")
cat("  PART 2d: COLLINEARITY, then BLOCK drops\n")
cat(strrep("=", 84), "\n\n", sep = "")

num <- X[, vapply(X, is.numeric, logical(1)), drop = FALSE]
cm <- cor(num, use = "pairwise.complete.obs")
cm[upper.tri(cm, diag = TRUE)] <- NA
pairs <- which(abs(cm) > 0.5, arr.ind = TRUE)
if (nrow(pairs)) {
  cat("predictor pairs with |r| > 0.5:\n")
  for (i in seq_len(nrow(pairs)))
    cat(sprintf("   %-14s %-14s r = %+.3f\n", rownames(cm)[pairs[i, 1]],
                colnames(cm)[pairs[i, 2]], cm[pairs[i, 1], pairs[i, 2]]))
} else cat("no numeric pair exceeds |r| = 0.5\n")
cat(sprintf("\n   R2 of month regressed on soil_temp: %.3f  (soil_temp is a monthly constant)\n",
            summary(lm(month ~ soil_temp, data = X))$r.squared))

BLOCKS <- list(
  "temporal (soil_temp + month)"                  = c("soil_temp", "month"),
  "spatial (elev + dist_river + local_ba)"        = c("elevation_m", "dist_river_m", "local_ba_m2"),
  "tree identity (species + dbh)"                 = c("species", "dbh_m"),
  "moisture"                                      = "moisture",
  "height"                                        = "height_cm"
)
block_tab <- bind_rows(lapply(names(BLOCKS), function(b) {
  o <- vapply(1:5, function(s)
    fit(X[, setdiff(BEST, BLOCKS[[b]]), drop = FALSE], yb, seed = s)$r.squared, numeric(1))
  g <- grouped_r2(X[, setdiff(BEST, BLOCKS[[b]]), drop = FALSE], yb, grp)
  data.frame(block = b, n_dropped = length(BLOCKS[[b]]), oob_without = mean(o),
             delta_oob = mean(full_oob) - mean(o),
             delta_sd = sqrt(var(full_oob) / 5 + var(o) / 5),
             grouped_without = g[["grouped_r2"]])
}))
g_full <- grouped_r2(X, yb, grp)
cat(sprintf("\nfull model: OOB R2 = %.4f, grouped R2 = %.4f\n\n",
            mean(full_oob), g_full[["grouped_r2"]]))
print(as.data.frame(block_tab[order(-block_tab$delta_oob), ]), row.names = FALSE, digits = 3)
write.csv(block_tab, file.path(outdir, "predictor_block_drops.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
SIG <- shadow_tab %>%
  left_join(alt_tab, by = "predictor") %>%
  left_join(drop_tab, by = "predictor") %>%
  arrange(desc(delta_oob))
write.csv(SIG, file.path(outdir, "predictor_significance.csv"), row.names = FALSE)

cat("\n", strrep("=", 84), "\n", sep = "")
cat("  CONSENSUS  (agreement across all three tests)\n")
cat(strrep("=", 84), "\n\n", sep = "")
print(SIG[, c("predictor", "verdict", "shadow_p", "altmann_p", "delta_oob", "delta_over_sd")],
      row.names = FALSE, digits = 3)

# --- figures ------------------------------------------------------------------
p_dbh <- dbh_res %>%
  ggplot(aes(reorder(encoding, oob_r2), oob_r2)) +
  geom_col(fill = "#2166AC") +
  geom_errorbar(aes(ymin = oob_r2 - oob_sd, ymax = oob_r2 + oob_sd), width = .2) +
  geom_point(aes(y = grouped_r2), colour = "#B2182B", size = 2) +
  coord_flip() +
  labs(title = "DBH encoding: OOB R2 (bars, +/-1 sd over 5 seeds) and grouped CV R2 (red)",
       x = NULL, y = expression(R^2)) + theme_minimal(base_size = 9)
ggsave(file.path(outdir, "dbh_encoding.png"), p_dbh, width = 7.5, height = 3.2,
       dpi = 190, bg = "white")

p_sig <- SIG %>%
  ggplot(aes(reorder(predictor, delta_oob), delta_oob, fill = verdict)) +
  geom_col() +
  geom_errorbar(aes(ymin = delta_oob - delta_sd, ymax = delta_oob + delta_sd), width = .2) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c(confirmed = "#2166AC", tentative = "#F4A582",
                               rejected = "#B2182B")) +
  coord_flip() +
  labs(title = "Predictor significance: drop-column loss in OOB R2, coloured by shadow-feature verdict",
       x = NULL, y = expression(Delta*R[OOB]^2~"when dropped")) +
  theme_minimal(base_size = 9)
ggsave(file.path(outdir, "predictor_significance.png"), p_sig, width = 8, height = 3.6,
       dpi = 190, bg = "white")

cat("\nWrote:\n  ", file.path(outdir, "dbh_encoding_comparison.csv"),
    "\n  ", file.path(outdir, "predictor_significance.csv"),
    "\n  ", file.path(outdir, "dbh_encoding.png"),
    "\n  ", file.path(outdir, "predictor_significance.png"), "\n", sep = "")
