#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_model_selection.R
# ------------------------------------------------------------------------------
# Systematic, diagnosed comparison of candidate model specifications.
#
# PRINCIPLE: a predictor may only be used if it can be evaluated at EVERY inventory
# stem / grid cell, not merely at training points. Note that soil and air temperature
# are monthly plot-level constants at prediction time -- they can explain WHEN flux is
# high but never WHERE, so they cannot shape the upscaled spatial field.
#
# Candidate predictors (all available at prediction):
#   species        factor, from the inventory
#   dbh_m          from the inventory
#   moisture       December-anchored surface, monthly affine -> varies in space and time
#   soil_temp      monthly constant
#   air_temp       monthly constant
#   month          seasonality
#   height_cm      measurement height (set to the 0-2 m band mid-point at prediction)
#   dist_river_m   derived: distance to the traced river, varies in space
#   local_ba_m2    derived: basal area within 10 m, varies in space
#   elevation_m    IDW-interpolated from 183 surveyed points (plots, wells, river,
#                  soil-moisture stations, trees) spanning the inventory; median stem
#                  is 8.6 m from its nearest measurement
#   x, y           coordinates
#
# ON TEMPERATURE: air and soil temperature come from the met tower (40,458 hourly
# records, including Tsoil). They are spatially uniform within a month but resolve the
# SEASONAL cycle -- which is essential, because the annual budget is the sum of monthly
# stand fluxes. They carry the temporal signal; moisture, elevation, distance-to-river
# and stand structure carry the spatial one. Both are needed.
#
# Every configuration is scored on BOTH out-of-bag and GROUPED CV (whole trees /
# collars held out, which is what upscaling requires), and each gets a full diagnostic
# panel: permutation importance, predicted-vs-observed, partial dependence, residuals.
#
# Outputs: outputs/revision/model_selection_summary.csv
#          outputs/revision/modelsel_<config>.png   (one panel per configuration)
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(ggplot2); library(readxl)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
INV <- rf_workflow_data$PLACEHOLDER_INVENTORY

MLAT <- 111320; MLON <- 111320 * cos(41.99 * pi / 180)

# --- DBH unit corrections (as the pipeline applies) -------------------------
fix_dbh <- function(dbh, species) {
  d <- ifelse(!is.na(dbh) & dbh > 3, dbh / 100, dbh)
  d <- ifelse(grepl("Betula", species) & !is.na(d) & d * 100 > 200, d / 10, d)
  d <- ifelse(species == "Pinus strobus" & !is.na(d) & d * 100 > 230, d / 10, d)
  d <- ifelse(species == "Kalmia latifolia" & !is.na(d) & d * 100 > 100, d / 100,
       ifelse(species == "Kalmia latifolia" & !is.na(d) & d * 100 > 10,  d / 10, d))
  d
}
INV$dbh_fix <- fix_dbh(INV$dbh_m, INV$species)

# --- derived spatial covariates ---------------------------------------------
riv <- suppressMessages(read_excel("data/raw/inventory/spatial_data/River.xlsx"))
names(riv) <- make.names(names(riv))
rok <- !is.na(riv$Latitude) & !is.na(riv$Longitude)
rx <- riv$Longitude[rok] * MLON; ry <- riv$Latitude[rok] * MLAT

dist_river <- function(x, y) {
  vapply(seq_along(x), function(i) {
    if (is.na(x[i]) || is.na(y[i])) return(NA_real_)
    min(sqrt((x[i] * MLON - rx)^2 + (y[i] * MLAT - ry)^2))
  }, numeric(1))
}
# elevation: inverse-distance-weighted from the surveyed point set
ELEV <- read.csv("outputs/tables/combined_elevation_data.csv", stringsAsFactors = FALSE)
ex <- ELEV$Longitude * MLON; ey <- ELEV$Latitude * MLAT; ez <- ELEV$Elevation
elev_at <- function(x, y, k = 6) {
  vapply(seq_along(x), function(i) {
    if (is.na(x[i]) || is.na(y[i])) return(NA_real_)
    dd <- sqrt((x[i] * MLON - ex)^2 + (y[i] * MLAT - ey)^2)
    o <- order(dd)[seq_len(min(k, length(dd)))]
    w <- 1 / pmax(dd[o], 0.5)^2
    sum(w * ez[o]) / sum(w)
  }, numeric(1))
}

ix <- INV$x * MLON; iy <- INV$y * MLAT; iba <- pi * (INV$dbh_fix / 2)^2
local_ba <- function(x, y) {
  vapply(seq_along(x), function(i) {
    if (is.na(x[i]) || is.na(y[i])) return(NA_real_)
    dd <- (x[i] * MLON - ix)^2 + (y[i] * MLAT - iy)^2
    sum(iba[dd <= 100], na.rm = TRUE)          # 10 m radius
  }, numeric(1))
}

# --- assemble the candidate training frame ----------------------------------
d <- tree_train_complete
yraw <- d$stem_flux_corrected
keep <- is.finite(yraw)
d <- d[keep, ]; yraw <- yraw[keep]
Xc <- as.data.frame(X_tree)[keep, ]

cat("deriving spatial covariates for", nrow(d), "tree measurements...\n")
TREE <- data.frame(
  species      = factor(d$species_clean),
  dbh_m        = d$dbh_m,
  moisture     = d$soil_moisture_at_tree,
  soil_temp    = d$soil_temp_C_mean,
  air_temp     = d$air_temp_C_mean,
  month        = d$month,
  height_cm    = ifelse(is.na(d$measurement_height_cm), 125, d$measurement_height_cm),
  dist_river_m = dist_river(d$x, d$y),
  local_ba_m2  = local_ba(d$x, d$y),
  elevation_m  = elev_at(d$x, d$y),
  x = d$x, y = d$y
)
cat("  dist_river non-NA:", sum(!is.na(TREE$dist_river_m)),
    " local_ba non-NA:", sum(!is.na(TREE$local_ba_m2)), "of", nrow(TREE), "\n\n")

# --- scoring ----------------------------------------------------------------
r2 <- function(a, b) 1 - sum((a - b)^2) / sum((a - mean(a))^2)

score <- function(X, y, grp, label, make_plot = TRUE) {
  X <- as.data.frame(X)
  rf <- ranger(x = X, y = y, num.trees = 800, min.node.size = 5,
               mtry = max(1, floor(sqrt(ncol(X)))), num.threads = 1,
               importance = "permutation", seed = 42)
  ut <- unique(grp); set.seed(1)
  gm <- setNames(sample(rep(1:5, length.out = length(ut))), ut); f <- unname(gm[grp])
  pr <- rep(NA_real_, nrow(X))
  for (k in 1:5) {
    tr <- which(f != k); te <- which(f == k)
    m <- ranger(x = X[tr, , drop = FALSE], y = y[tr], num.trees = 800, min.node.size = 5,
                mtry = max(1, floor(sqrt(ncol(X)))), num.threads = 1, seed = 42)
    pr[te] <- predict(m, X[te, , drop = FALSE])$predictions
  }
  s <- data.frame(config = label, n = nrow(X), n_pred = ncol(X),
                  oob_r2 = rf$r.squared, grouped_r2 = r2(y, pr),
                  mean_ratio = mean(pr) / mean(y),
                  spearman = cor(y, pr, method = "spearman"))

  if (make_plot) {
    imp <- sort(rf$variable.importance, decreasing = TRUE)
    p1 <- data.frame(v = names(imp), i = as.numeric(imp)) %>%
      ggplot(aes(reorder(v, i), i)) + geom_col(fill = "#2166AC") + coord_flip() +
      labs(title = "permutation importance", x = NULL, y = NULL) + theme_minimal(base_size = 8)
    p2 <- data.frame(obs = y, pred = pr) %>% ggplot(aes(obs, pred)) +
      geom_point(alpha = .3, size = .8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
      labs(title = sprintf("grouped CV: R2 %.3f, ratio %.2f", s$grouped_r2, s$mean_ratio),
           x = "observed", y = "predicted") + theme_minimal(base_size = 8)
    p3 <- data.frame(pred = pr, res = y - pr) %>% ggplot(aes(pred, res)) +
      geom_hline(yintercept = 0, colour = "grey40") + geom_point(alpha = .3, size = .8) +
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
      labs(title = "residual vs fitted", x = "predicted", y = "residual") +
      theme_minimal(base_size = 8)
    top <- names(imp)[seq_len(min(6, length(imp)))]
    pdp <- bind_rows(lapply(top, function(v) {
      xs <- X[sample(nrow(X), min(250, nrow(X))), , drop = FALSE]
      if (is.factor(X[[v]])) {
        lv <- levels(droplevels(X[[v]]))
        bind_rows(lapply(lv, function(l) { xs[[v]] <- factor(l, levels = levels(X[[v]]))
          data.frame(variable = v, lab = l, value = NA_real_,
                     yhat = mean(predict(rf, xs)$predictions, na.rm = TRUE)) }))
      } else {
        g <- seq(quantile(X[[v]], .02, na.rm = TRUE), quantile(X[[v]], .98, na.rm = TRUE),
                 length.out = 20)
        bind_rows(lapply(g, function(q) { xs[[v]] <- q
          data.frame(variable = v, lab = NA_character_, value = q,
                     yhat = mean(predict(rf, xs)$predictions, na.rm = TRUE)) }))
      }
    }))
    pdn <- pdp %>% filter(!is.na(value))
    p4 <- ggplot(pdn, aes(value, yhat)) + geom_line(colour = "#B2182B", linewidth = .7) +
      facet_wrap(~ variable, scales = "free_x", ncol = 3) +
      labs(title = "partial dependence", x = NULL,
           y = expression(CH[4]~(nmol~m^-2~s^-1))) + theme_minimal(base_size = 8)
    ggsave(file.path(outdir, paste0("modelsel_", gsub("[^A-Za-z0-9]+", "_", label), ".png")),
           gridExtra::arrangeGrob(p1, p2, p3, p4, layout_matrix = rbind(c(1,2,3), c(4,4,4)),
                                  top = label),
           width = 11, height = 7, dpi = 190, bg = "white")
  }
  s
}

# --- configurations ---------------------------------------------------------
h125 <- is.na(d$measurement_height_cm) | d$measurement_height_cm == 125

CONFIGS <- list(
  list(lab = "A current (34 engineered)",     X = Xc,                                    rows = rep(TRUE, nrow(d))),
  list(lab = "B minimal 6 (+height)",         X = TREE[, c("species","dbh_m","moisture","soil_temp","month","height_cm")], rows = rep(TRUE, nrow(d))),
  list(lab = "C minimal 5, 125cm only",       X = TREE[, c("species","dbh_m","moisture","soil_temp","month")],             rows = h125),
  list(lab = "D minimal 5, all heights pooled", X = TREE[, c("species","dbh_m","moisture","soil_temp","month")],           rows = rep(TRUE, nrow(d))),
  list(lab = "E full candidate set",          X = TREE,                                  rows = rep(TRUE, nrow(d))),
  list(lab = "F spatial + temps",             X = TREE[, c("species","dbh_m","moisture","soil_temp","month","height_cm","dist_river_m","local_ba_m2","elevation_m")], rows = rep(TRUE, nrow(d))),
  list(lab = "G minimal 6 + elevation",        X = TREE[, c("species","dbh_m","moisture","soil_temp","month","height_cm","elevation_m")], rows = rep(TRUE, nrow(d)))
)

res <- bind_rows(lapply(CONFIGS, function(cf) {
  r <- cf$rows
  score(cf$X[r, , drop = FALSE], yraw[r], d$tree_id[r], cf$lab)
}))

cat("\n", strrep("=", 88), "\n", sep = "")
cat("  TREE MODEL SELECTION  (raw nmol m-2 s-1; grouped CV holds out whole trees)\n")
cat(strrep("=", 88), "\n\n", sep = "")
print(as.data.frame(res), row.names = FALSE, digits = 3)
write.csv(res, file.path(outdir, "model_selection_summary.csv"), row.names = FALSE)
cat("\nDiagnostic panels written to outputs/revision/modelsel_*.png\n")
