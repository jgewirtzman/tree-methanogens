#!/usr/bin/env Rscript
# ==============================================================================
# figS24_rf-model-summary.R -- Figure S26: the stem and soil random forests
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a, b) observed against out-of-bag predicted flux (stems, soils), hexbin density,
#          1:1 line; CCC and R2 on the out-of-bag predictions, grouped-CV R2
#          (outputs/data/rf_grouped_cv.csv), n. Equal axes, so 1:1 sits at 45 degrees;
#          both on the arcsinh flux scale (display only).
#   (c, d) permutation importance on held-out data, share of the total, +/- 1 SE over
#          50 random 70/30 splits (20 permutations each). Importance is permutation,
#          not ranger's impurity importance, which favours continuous predictors.
#          Each split refits a forest with the fitted model's min.node.size and mtry
#          but 800 trees (the fitted models have 10,000). Splits, permutations and
#          trees: RUN$importance_splits, $importance_permutations, $refit_num_trees
#          (code/lib/run_mode.R; 50, 20 and 800 in the full run).
#   (e, f) partial dependence for the predictors, 2nd-98th percentile, with a rug of
#          the observed values (the tails rest on few measurements).
# Inputs: outputs/models/RF_MODELS.RData (TreeRF, SoilRF) and TRAINING_DATA.RData
#   (05_model/02_rf_models.R); outputs/data/rf_grouped_cv.csv (05_model/04_rf_grouped_cv.R)
# Output: outputs/figures/generated/figS21_rf_model_summary.png; importance shares to
#   the console
# ==============================================================================
source("code/lib/outputs.R")
suppressMessages({library(dplyr); library(ranger); library(ggplot2)
                  library(patchwork); library(tidyr)})
source("code/lib/rf_spec.R")   # RF_THREADS; RUN (code/lib/run_mode.R)
source("code/lib/fig_style.R")      # SI house style: theme_si(), si_tags(), si_save(), COL_*
source("code/lib/flux_transform.R") # flux_asinh_trans(): the flux axis of Figs 1-3 and S26
set.seed(42)
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")

# stem green / soil brown, as across the SI; hexbin counts magma on a log10 scale
TREE_COL <- COL_STEM; SOIL_COL <- COL_SOIL

Xt <- as.data.frame(X_tree); yt <- tree_train_complete$stem_flux_corrected
Xs <- as.data.frame(X_soil); ys <- soil_train_complete$y_flux
stopifnot(length(yt) == nrow(Xt), length(ys) == nrow(Xs),
          length(TreeRF$predictions) == nrow(Xt), length(SoilRF$predictions) == nrow(Xs))

# ---- permutation importance, held out ------------------------------------------
# A group of predictors is permuted with ONE shared permutation (each group is one
# predictor in the current models). 50 splits: with 10 the stem shares moved by up to
# 7 points between seeds; the split, not the forest size, is the noise.
GROUPS <- list(
  TreeRF = list(species = "species", diameter = "dbh_m",
                `soil moisture` = "soil_moisture_at_tree",
                temperature = "soil_temp_C",
                `measurement height` = "height_cm"),
  SoilRF = list(`soil moisture` = "soil_moisture_at_site",
                temperature = "soil_temp_C",
                month = "month"))
stopifnot(setequal(unlist(GROUPS$TreeRF), names(Xt)), setequal(unlist(GROUPS$SoilRF), names(Xs)))

gimp <- function(fit, X, y, groups, lab, nsplit = RUN$importance_splits, nrep = RUN$importance_permutations) {
  k <- is.finite(y) & complete.cases(X); X <- X[k, , drop = FALSE]; y <- y[k]
  R <- replicate(nsplit, {
    tr <- sample(nrow(X), floor(0.7 * nrow(X)))
    m  <- ranger(x = X[tr, , drop = FALSE], y = y[tr], num.trees = RUN$refit_num_trees,
                 min.node.size = fit$min.node.size, mtry = fit$mtry, num.threads = RF_THREADS,
                 seed = sample.int(1e6, 1))
    Xh <- X[-tr, , drop = FALSE]; yh <- y[-tr]
    base <- mean((yh - predict(m, Xh, num.threads = RF_THREADS)$predictions)^2)
    sapply(groups, function(g) mean(replicate(nrep, {
      Xp <- Xh; i <- sample(nrow(Xp))
      for (v in g) Xp[[v]] <- Xp[[v]][i]          # ONE shared permutation per group
      mean((yh - predict(m, Xp, num.threads = RF_THREADS)$predictions)^2) - base
    })))
  })
  mu <- pmax(rowMeans(R), 0)
  data.frame(model = lab, label = names(groups), imp = as.numeric(mu),
             se = apply(R, 1, sd) / sqrt(nsplit)) %>%
    mutate(share = imp / sum(imp), share_se = se / sum(imp))
}

IMP_T <- gimp(TreeRF, Xt, yt, GROUPS$TreeRF, "TreeRF")
IMP_S <- gimp(SoilRF, Xs, ys, GROUPS$SoilRF, "SoilRF")

# display names for the predictor groups, as in the partial-dependence panels
IMP_LAB <- c(species = "Species", diameter = "DBH", `soil moisture` = "Soil VWC",
             temperature = "Soil temperature", `measurement height` = "Measurement height",
             month = "Month")
imp_panel <- function(D, ttl, col) {
  stopifnot(all(D$label %in% names(IMP_LAB)))
  D <- D %>% arrange(share) %>% mutate(label = factor(IMP_LAB[label], levels = IMP_LAB[label]))
  ggplot(D, aes(share, label)) +
    geom_col(width = 0.68, fill = col) +
    geom_errorbarh(aes(xmin = pmax(0, share - share_se), xmax = share + share_se),
                   height = 0.22, linewidth = 0.35, colour = "grey25") +
    geom_text(aes(x = share + share_se, label = sprintf("%.0f%%", 100 * share)),
              hjust = -0.3, size = 3) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1.25), breaks = seq(0, 1, 0.25),
                       expand = c(0, 0)) +
    labs(title = ttl, x = "Share of total importance", y = NULL) +
    theme_si() +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = "black"))
}
pc_ <- imp_panel(IMP_T, NULL, TREE_COL)
pd_ <- imp_panel(IMP_S, NULL, SOIL_COL)

# ---- observed against out-of-bag predicted ------------------------------------
ccc <- function(o, p) {
  2 * cov(o, p) / (var(o) + var(p) + (mean(o) - mean(p))^2)
}
obsT <- data.frame(o = yt, p = TreeRF$predictions) %>% filter(is.finite(o), is.finite(p))
obsS <- data.frame(o = ys, p = SoilRF$predictions) %>% filter(is.finite(o), is.finite(p))

g <- read.csv("outputs/data/rf_grouped_cv.csv", stringsAsFactors = FALSE)
R2G <- c(g$grouped_r2[g$model == "TreeRF"], g$grouped_r2[g$model == "SoilRF"])
stopifnot(length(R2G) == 2, is.finite(R2G))

obs_panel <- function(D, o, p, ttl, col, r2g, brks) {
  r2 <- 1 - sum((D[[o]] - D[[p]])^2) / sum((D[[o]] - mean(D[[o]]))^2)
  cat(sprintf("%s: CCC = %.3f, R2 (out-of-bag) = %.3f, R2 (grouped CV) = %.3f, n = %d\n",
              ttl, ccc(D[[o]], D[[p]]), r2, r2g, nrow(D)))
  lab <- sprintf("CCC = %.2f\nR\u00b2 = %.2f (out-of-bag)\nR\u00b2 = %.2f (grouped CV)\nn = %d",
                 ccc(D[[o]], D[[p]]), r2, r2g, nrow(D))
  # equal axes: 1:1 at 45 degrees, so under-prediction of large fluxes is visible.
  # Both axes on the arcsinh flux scale (cofactor 0.1 nmol m-2 s-1, as in Figs 1-3 and
  # S26): linear near zero, where most measurements lie, and compressed in the tails,
  # so every point stays in view without the bulk shrinking into one corner. Display
  # only: the statistics above are on the untransformed fluxes.
  rng <- range(c(D[[o]], D[[p]]), finite = TRUE)
  ggplot(D, aes(.data[[p]], .data[[o]])) +
    geom_hex(bins = 40) +
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.7, alpha = 0.8) +
    # boxed so the hexes do not overprint the statistics
    annotate("label", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.06,
             label = lab, size = 2.9, lineheight = 1.05,
             fill = scales::alpha("white", 0.85), linewidth = 0) +
    scale_fill_viridis_c(name = "Count", trans = "log10", option = "magma") +
    scale_x_continuous(transform = flux_asinh_trans(), breaks = brks, labels = as.character) +
    scale_y_continuous(transform = flux_asinh_trans(), breaks = brks, labels = as.character) +
    coord_equal(xlim = rng, ylim = rng) +
    labs(title = ttl,
         x = expression(Predicted~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
         y = expression(Observed~CH[4]~flux~(nmol~m^{-2}~s^{-1}))) +
    theme_si() +
    theme(legend.position = "bottom", legend.key.width = unit(0.7, "cm"), legend.key.height = unit(0.22, "cm"),
          legend.title = element_text(size = rel(0.8), vjust = 0.9), legend.text = element_text(size = rel(0.7)),
          legend.margin = margin(0, 0, 0, 0))
}
pa_ <- obs_panel(obsT, "o", "p", NULL, TREE_COL, R2G[1], brks = c(-0.3, 0, 0.1, 0.3, 1, 3))
pb_ <- obs_panel(obsS, "o", "p", NULL, SOIL_COL, R2G[2], brks = c(-3, -1, -0.3, 0, 0.3, 1, 3))

# ---- partial dependence ------------------------------------------------------
pdp <- function(m, dat, xv, n = 40) {
  q <- quantile(dat[[xv]], c(.02, .98), na.rm = TRUE)
  g <- seq(q[1], q[2], length.out = n)
  data.frame(x = g, y = vapply(g, function(v) { D <- dat; D[[xv]] <- v
    mean(predict(m, D, num.threads = 1)$predictions) }, numeric(1)))
}
PD_LAB <- c(soil_moisture_at_tree = "Soil VWC (%)",
            soil_moisture_at_site = "Soil VWC (%)",
            soil_temp_C = "Soil temperature (\u00b0C)",
            dbh_m = "DBH (cm)",
            height_cm = "Measurement height (cm)",
            month = "Month")
# display units only: VWC m3 m-3 -> %, DBH m -> cm (the forests see the stored units)
PD_SCALE <- c(soil_moisture_at_tree = 100, soil_moisture_at_site = 100, dbh_m = 100)
pd_scale <- function(v) if (v %in% names(PD_SCALE)) PD_SCALE[[v]] else 1

pd_grid <- function(m, X, vars, ttl, col) {
  D <- bind_rows(lapply(vars, function(v)
    pdp(m, X, v) %>% mutate(x = x * pd_scale(v), var = PD_LAB[[v]])))
  D$var <- factor(D$var, levels = unname(PD_LAB[vars]))
  # rug of the observed values behind each curve
  RUG <- bind_rows(lapply(vars, function(v)
    data.frame(x = X[[v]] * pd_scale(v), var = PD_LAB[[v]])))
  RUG$var <- factor(RUG$var, levels = unname(PD_LAB[vars]))
  ggplot(D, aes(x, y)) +
    geom_rug(data = RUG, aes(x = x), inherit.aes = FALSE,
             alpha = 0.14, length = unit(0.025, "npc"), colour = col) +
    geom_line(color = col, linewidth = 0.9) +
    facet_wrap(~var, scales = "free_x", nrow = 1) +
    # month on calendar quarters; other predictors on default breaks
    scale_x_continuous(breaks = function(l) if (l[1] >= 1 && l[2] <= 12.5) c(3, 6, 9, 12)
                                            else scales::breaks_extended()(l)) +
    labs(title = ttl, x = NULL, y = lab_flux) +
    theme_si() + theme(panel.spacing.x = unit(0.8, "lines"), strip.text = element_text(size = rel(0.78)))
}
pe_ <- pd_grid(TreeRF, Xt, c("soil_moisture_at_tree", "dbh_m",
                             "height_cm", "soil_temp_C"),
               NULL, TREE_COL)
pf_ <- pd_grid(SoilRF, Xs, c("soil_moisture_at_site", "soil_temp_C", "month"),
               NULL, SOIL_COL)

# rows: observed vs predicted and importance (stems left, soils right), then partial
# dependence for stems and for soils, each across the full width
# free(): the partial-dependence y titles sit by their own axes, not aligned to (c)'s labels
fig <- (pa_ | pb_) / (pc_ | pd_) / free(pe_, type = "label", side = "l") /
  free(pf_ + plot_spacer() + plot_layout(widths = c(3, 1)), type = "label", side = "l") +
  plot_layout(heights = c(1.45, 0.8, 0.75, 0.75)) + si_tags()
si_save(out_path("figS21_rf_model_summary.png"), fig, height = 10)

cat("=== grouped permutation importance (held-out), share within model ===\n")
for (M in list(IMP_T, IMP_S)) {
  cat(sprintf("\n--- %s ---\n", M$model[1]))
  print(as.data.frame(M %>% arrange(-share) %>%
        transmute(predictor = label, share = sprintf("%.1f%%", 100 * share))),
        row.names = FALSE)
}
cat("\nwritten: outputs/figures/generated/figS21_rf_model_summary.png\n")
