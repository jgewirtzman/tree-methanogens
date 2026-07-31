#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_predictor_selection.R
# ------------------------------------------------------------------------------
# Which predictors carry real signal, and which are no better than noise?
#
# METHOD: random-reference predictors ("is it better than noise?").
#
#   Three reference columns are added ONCE to the real predictor set:
#     rand_uniform  ~ U(0,1)          continuous, low information
#     rand_normal   ~ N(0,1)          continuous, different distribution
#     rand_factor   13 random levels  controls for the cardinality bias that inflates
#                                     importance for many-levelled factors like species
#   A predictor is retained if its permutation importance exceeds the LARGEST of the
#   three references, averaged over repeated fits.
#
#   WHY NOT SOMETHING FANCIER -- both alternatives were tried and both failed here:
#   * Shadow features (Boruta), one permuted copy per predictor: doubles the design
#     matrix and halved the tree model's OOB R2 (0.275 -> 0.165). With mtry = 4 of 16
#     columns, half of them noise, real predictors are often never offered to a split.
#     Every importance went negative; the test measured dilution, not signal.
#   * Response permutation (Altmann et al. 2010 PIMP): with 50 permutations the
#     smallest attainable p is 1/51 = 0.0196, so after Bonferroni across 7-11
#     predictors NOTHING can reach p < 0.05 -- it rejected soil predictors whose
#     observed importance was 5x the null 95th percentile. Arithmetically incapable of
#     confirming anything at this permutation count.
#
#   Adding three reference columns costs almost no dilution and the comparison is
#   directly interpretable.
#
# Output: outputs/revision/predictor_selection.csv
#         outputs/revision/fig_predictor_selection.png
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(ggplot2); library(readxl)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
INV <- rf_workflow_data$PLACEHOLDER_INVENTORY
MLAT <- 111320; MLON <- 111320 * cos(41.99 * pi / 180)

fix_dbh <- function(dbh, species) {
  d <- ifelse(!is.na(dbh) & dbh > 3, dbh / 100, dbh)
  d <- ifelse(grepl("Betula", species) & !is.na(d) & d*100 > 200, d/10, d)
  d <- ifelse(species == "Pinus strobus" & !is.na(d) & d*100 > 230, d/10, d)
  ifelse(species == "Kalmia latifolia" & !is.na(d) & d*100 > 100, d/100,
  ifelse(species == "Kalmia latifolia" & !is.na(d) & d*100 > 10,  d/10, d))
}
INV$dbh_fix <- fix_dbh(INV$dbh_m, INV$species)

riv <- suppressMessages(read_excel("data/raw/inventory/spatial_data/River.xlsx"))
names(riv) <- make.names(names(riv)); rok <- !is.na(riv$Latitude)
rx <- riv$Longitude[rok]*MLON; ry <- riv$Latitude[rok]*MLAT
ELEV <- read.csv("outputs/tables/combined_elevation_data.csv", stringsAsFactors = FALSE)
ex <- ELEV$Longitude*MLON; ey <- ELEV$Latitude*MLAT; ez <- ELEV$Elevation
ix <- INV$x*MLON; iy <- INV$y*MLAT; iba <- pi*(INV$dbh_fix/2)^2

dist_river <- function(x,y) vapply(seq_along(x), function(i)
  if (is.na(x[i])) NA_real_ else min(sqrt((x[i]*MLON-rx)^2 + (y[i]*MLAT-ry)^2)), numeric(1))
elev_at <- function(x,y,k=6) vapply(seq_along(x), function(i) {
  if (is.na(x[i])) return(NA_real_)
  dd <- sqrt((x[i]*MLON-ex)^2 + (y[i]*MLAT-ey)^2); o <- order(dd)[1:min(k,length(dd))]
  w <- 1/pmax(dd[o],0.5)^2; sum(w*ez[o])/sum(w) }, numeric(1))
local_ba <- function(x,y) vapply(seq_along(x), function(i) {
  if (is.na(x[i])) return(NA_real_)
  sum(iba[(x[i]*MLON-ix)^2 + (y[i]*MLAT-iy)^2 <= 100], na.rm=TRUE) }, numeric(1))

# ---------------------------------------------------------------------------
ref_test <- function(X, y, label, n_runs = 20) {
  X <- as.data.frame(X); vars <- names(X)
  imp <- matrix(NA_real_, n_runs, length(vars) + 3,
                dimnames = list(NULL, c(vars, "rand_uniform", "rand_normal", "rand_factor")))
  for (r in seq_len(n_runs)) {
    set.seed(3000 + r)
    Xr <- X
    Xr$rand_uniform <- runif(nrow(X))
    Xr$rand_normal  <- rnorm(nrow(X))
    Xr$rand_factor  <- factor(sample(1:13, nrow(X), replace = TRUE))
    rf <- ranger(x = Xr, y = y, num.trees = 500, min.node.size = 5,
                 mtry = max(1, floor(sqrt(ncol(Xr)))), num.threads = 1,
                 importance = "permutation", seed = 3000 + r)
    imp[r, names(rf$variable.importance)] <- rf$variable.importance
  }
  refs <- imp[, c("rand_uniform", "rand_normal", "rand_factor"), drop = FALSE]
  thresh_per_run <- apply(refs, 1, max, na.rm = TRUE)
  m <- colMeans(imp, na.rm = TRUE)
  beats <- vapply(vars, function(v) mean(imp[, v] > thresh_per_run, na.rm = TRUE), numeric(1))
  data.frame(model = label, predictor = vars,
             mean_importance = as.numeric(m[vars]),
             ref_threshold = mean(thresh_per_run),
             frac_runs_beating_noise = as.numeric(beats),
             decision = ifelse(beats >= 0.95, "Keep",
                        ifelse(beats >= 0.50, "Marginal", "Drop")),
             row.names = NULL) %>% arrange(desc(mean_importance)) -> out
  attr(out, "refs") <- colMeans(refs, na.rm = TRUE)
  out
}

# --- TREE -------------------------------------------------------------------
d <- tree_train_complete; y <- d$stem_flux_corrected; k <- is.finite(y)
d <- d[k, ]; y <- y[k]
TREE <- data.frame(
  species      = factor(d$species_clean),
  dbh_m        = d$dbh_m,
  dbh_within_z = d$dbh_within_z,
  moisture     = d$soil_moisture_at_tree,
  soil_temp    = d$soil_temp_C_mean,
  air_temp     = d$air_temp_C_mean,
  month        = d$month,
  height_cm    = ifelse(is.na(d$measurement_height_cm), 125, d$measurement_height_cm),
  dist_river_m = dist_river(d$x, d$y),
  local_ba_m2  = local_ba(d$x, d$y),
  elevation_m  = elev_at(d$x, d$y)
)
cat("random-reference testing: TREE (", nrow(TREE), "obs,", ncol(TREE), "candidates )\n")
rt <- ref_test(TREE, y, "Tree stem")

# --- SOIL -------------------------------------------------------------------
ds <- soil_train_complete; ys <- ds$soil_flux_umol_m2_s; ks <- is.finite(ys)
ds <- ds[ks, ]; ys <- ys[ks]; Xs0 <- as.data.frame(X_soil)[ks, ]
SOIL <- data.frame(
  moisture     = Xs0$soil_moisture_at_site,
  soil_temp    = Xs0$soil_temp_C_mean,
  air_temp     = Xs0$air_temp_C_mean,
  month        = ds$month,
  dist_river_m = dist_river(ds$x, ds$y),
  local_ba_m2  = local_ba(ds$x, ds$y),
  elevation_m  = elev_at(ds$x, ds$y)
)
cat("random-reference testing: SOIL (", nrow(SOIL), "obs,", ncol(SOIL), "candidates )\n\n")
rs <- ref_test(SOIL, ys, "Soil")

res <- bind_rows(rt, rs)
for (m in unique(res$model)) {
  cat(strrep("=", 74), "\n  ", toupper(m), " -- random-reference selection (20 runs)\n",
      strrep("=", 74), "\n", sep = "")
  q <- res[res$model == m, ]
  cat(sprintf("  %-16s %12s %12s %9s  %s\n", "predictor", "importance", "noise thr.", "beats", "decision"))
  for (i in seq_len(nrow(q)))
    cat(sprintf("  %-16s %12.4g %12.4g %8.0f%%  %s\n", q$predictor[i],
                q$mean_importance[i], q$ref_threshold[i],
                100*q$frac_runs_beating_noise[i], q$decision[i]))
  cat("\n")
}
write.csv(res, file.path(outdir, "predictor_selection.csv"), row.names = FALSE)

p <- res %>% mutate(predictor = reorder(predictor, mean_importance)) %>%
  ggplot(aes(predictor, mean_importance, fill = decision)) +
  geom_col() + coord_flip() + facet_wrap(~ model, scales = "free") +
  geom_hline(aes(yintercept = ref_threshold), linetype = "dashed", colour = "grey25") +
  scale_fill_manual(values = c(Keep = "#2166AC", Marginal = "#F4A582", Drop = "#B2182B")) +
  labs(title = "Predictor selection against random reference columns",
       subtitle = "bar = mean permutation importance; dashed = best of three random reference predictors",
       x = NULL, y = "permutation importance") +
  theme_minimal(base_size = 9)
ggsave(file.path(outdir, "fig_predictor_selection.png"), p, width = 10, height = 5, dpi = 200, bg = "white")
cat("Written: outputs/revision/predictor_selection.csv and fig_predictor_selection.png\n")
