#!/usr/bin/env Rscript
# ==============================================================================
# rev_figS21_rf-model-summary.R  --  SI: what the two random forests learned
# ------------------------------------------------------------------------------
# REPLACES two figures that disagreed with each other:
#   - figS15_rf_predictions.png (08_rf_publication_plots.R), which fed SI S21
#   - fig_model_findings.png (rev_fig_model_findings.R), which was generated on
#     every run but never assembled into the manuscript
# It keeps the older figure's LAYOUT -- observed-against-predicted, importance,
# partial dependence, for both models -- and the newer figure's METHOD.
#
# WHY THE METHOD HAD TO CHANGE. The superseded panel used ranger's default
# IMPURITY importance. Impurity gain counts variance reduction at splits and is
# well known to favour continuous, high-cardinality predictors; permutation
# importance measures the increase in held-out error when a predictor is
# shuffled, which is what "the model relies on this" means. The two metrics give
# opposite answers on these data:
#
#   predictor              impurity %   permutation %
#   SoilRF  soil moisture        9.0        44.0
#           air temperature      2.1        22.7
#           soil temperature    88.0        18.7
#
# The old SI panel therefore told the reader the SOIL model was a temperature
# model. It is a moisture model. That also contradicted the revision's own
# position that air temperature is a monthly mean and so is not identifiable
# separately from season.
#
# GROUPING. Soil and air temperature correlate at r = 0.97. Shuffling one leaves
# the other as a proxy, so ordinary permutation mutually suppresses both. They
# are permuted TOGETHER under one shared permutation, which asks the only
# question with an answer: how much does the model rely on temperature at all.
# Measured out-of-sample (fit on 70%, permute and score on 30%) so the baseline
# and permuted error are on the same footing.
#
# The monthly-prediction panels of the old figure are dropped: the seasonal cycle
# is already Figure 9(c), and its June/September dips are RF split-threshold
# artifacts rather than seasonality, so repeating them here invited over-reading.
#
# The old figure's DBH partial-dependence axis was labelled "DBH (within-species
# z-score)" (08_rf_publication_plots.R:892) while plotting dbh_m in metres --
# a stale label from the dropped 34-predictor specification. dbh_within_z is not
# a predictor of the frozen model.
#
# Output: outputs/figures/generated/figS21_rf_model_summary.png
# ==============================================================================
source("code/lib/outputs.R")
suppressMessages({library(dplyr); library(ranger); library(ggplot2)
                  library(patchwork); library(tidyr)})
set.seed(42)
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")

TREE_COL <- "#756BB1"; SOIL_COL <- "#2166ac"

Xt <- as.data.frame(X_tree); yt <- tree_train_complete$stem_flux_corrected
Xs <- as.data.frame(X_soil); ys <- soil_train_complete$y_asinh[seq_len(nrow(Xs))]

# ---- grouped permutation importance, held out --------------------------------
GROUPS <- list(
  TreeRF = list(species = "species", diameter = "dbh_m",
                `soil moisture` = "soil_moisture_at_tree",
                temperature = c("soil_temp_C_mean", "air_temp_C_mean"),
                `measurement height` = "height_cm"),
  SoilRF = list(`soil moisture` = "soil_moisture_at_site",
                temperature = c("soil_temp_C_mean", "air_temp_C_mean"),
                month = "month"))

gimp <- function(X, y, groups, lab, nsplit = 10, nrep = 20) {
  k <- is.finite(y) & complete.cases(X); X <- X[k, , drop = FALSE]; y <- y[k]
  R <- replicate(nsplit, {
    tr <- sample(nrow(X), floor(0.7 * nrow(X)))
    m  <- ranger(x = X[tr, , drop = FALSE], y = y[tr], num.trees = 800,
                 min.node.size = 5, mtry = max(1, floor(sqrt(ncol(X)))), num.threads = 1)
    Xh <- X[-tr, , drop = FALSE]; yh <- y[-tr]
    base <- mean((yh - predict(m, Xh, num.threads = 1)$predictions)^2)
    sapply(groups, function(g) mean(replicate(nrep, {
      Xp <- Xh; i <- sample(nrow(Xp))
      for (v in g) Xp[[v]] <- Xp[[v]][i]          # ONE shared permutation per group
      mean((yh - predict(m, Xp, num.threads = 1)$predictions)^2) - base
    })))
  })
  mu <- pmax(rowMeans(R), 0)
  data.frame(model = lab, label = names(groups), imp = as.numeric(mu),
             se = apply(R, 1, sd) / sqrt(nsplit)) %>%
    mutate(share = imp / sum(imp), share_se = se / sum(imp))
}

IMP_T <- gimp(Xt, yt, GROUPS$TreeRF, "TreeRF")
IMP_S <- gimp(Xs, ys, GROUPS$SoilRF, "SoilRF")

imp_panel <- function(D, ttl, col) {
  D <- D %>% arrange(share) %>% mutate(label = factor(label, levels = label))
  ggplot(D, aes(share, label)) +
    geom_col(width = 0.68, fill = col) +
    geom_errorbarh(aes(xmin = pmax(0, share - share_se), xmax = share + share_se),
                   height = 0.22, linewidth = 0.35, colour = "grey25") +
    geom_text(aes(x = share + share_se, label = sprintf("%.0f%%", 100 * share)),
              hjust = -0.3, size = 2.9) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1.16), expand = c(0, 0)) +
    labs(title = ttl, x = "share of total importance", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold", size = 10))
}
pc_ <- imp_panel(IMP_T, "(c) Tree feature importance", TREE_COL)
pd_ <- imp_panel(IMP_S, "(d) Soil feature importance", SOIL_COL)

# ---- observed against predicted, calibrated ----------------------------------
# Per-species calibration, UNCLAMPED, matching rev_predict_tree_flux_current.R.
ccc <- function(o, p) {
  2 * cov(o, p) / (var(o) + var(p) + (mean(o) - mean(p))^2)
}
d <- tree_train_complete
calT <- data.frame(sp = as.character(d$species_clean), o = d$stem_flux_corrected,
                   p = TreeRF$predictions) %>%
  filter(is.finite(o), is.finite(p)) %>%
  group_by(sp) %>% mutate(r = ifelse(mean(p) > 0, mean(o) / mean(p), 1)) %>%
  ungroup() %>% mutate(pc = p * r)
calS <- data.frame(o = ys, p = SoilRF$predictions) %>% filter(is.finite(o), is.finite(p))

R2G <- local({
  f <- "outputs/data/rf_grouped_cv.csv"
  if (!file.exists(f)) return(c(NA_real_, NA_real_))
  g <- read.csv(f, stringsAsFactors = FALSE)
  c(g$grouped_r2[g$model == "TreeRF"][1], g$grouped_r2[g$model == "SoilRF"][1])
})

obs_panel <- function(D, o, p, ttl, col, r2g, nlab) {
  r2 <- 1 - sum((D[[o]] - D[[p]])^2) / sum((D[[o]] - mean(D[[o]]))^2)
  sub <- sprintf("CCC %.3f | R2 %.3f out-of-bag, %.3f grouped CV | n = %d",
                 ccc(D[[o]], D[[p]]), r2, r2g, nrow(D))
  # Hexbin with a count scale, as in the figure this replaces: the response is
  # heavily zero-concentrated, so plain points hide the density that carries the
  # fit and show only the handful of large-flux outliers.
  ggplot(D, aes(.data[[p]], .data[[o]])) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
    geom_hex(bins = 34) +
    scale_fill_viridis_c(option = "magma", trans = "log10", name = "count") +
    labs(title = ttl, subtitle = sub,
         x = expression("predicted (nmol m"^-2*" s"^-1*")"),
         y = expression("observed (nmol m"^-2*" s"^-1*")")) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.42, "cm"),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
          plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 7.0, colour = "grey30"))
}
pa_ <- obs_panel(calT, "o", "pc", "(a) Tree stems",  TREE_COL, R2G[1])
pb_ <- obs_panel(calS, "o", "p",  "(b) Soil",        SOIL_COL, R2G[2])

# ---- partial dependence ------------------------------------------------------
# Air and soil temperature are shown SEPARATELY here even though importance pools
# them: the pooling exists because permutation cannot apportion a shared effect
# between collinear predictors, which is a different question from the shape of
# each marginal response. Read the two temperature curves as one driver seen
# twice, not as two independent effects.
pdp <- function(m, dat, xv, n = 40) {
  q <- quantile(dat[[xv]], c(.02, .98), na.rm = TRUE)
  g <- seq(q[1], q[2], length.out = n)
  data.frame(x = g, y = vapply(g, function(v) { D <- dat; D[[xv]] <- v
    mean(predict(m, D, num.threads = 1)$predictions) }, numeric(1)))
}
PD_LAB <- c(soil_moisture_at_tree = "soil moisture (m^3 m^-3)",
            soil_moisture_at_site = "soil moisture (m^3 m^-3)",
            soil_temp_C_mean = "soil temperature (C)",
            air_temp_C_mean  = "air temperature (C)",
            dbh_m = "DBH (m)",                       # was mislabelled a z-score
            height_cm = "measurement height (cm)",
            month = "month")

pd_grid <- function(m, X, vars, ttl, col) {
  D <- bind_rows(lapply(vars, function(v)
    pdp(m, X, v) %>% mutate(var = PD_LAB[[v]])))
  D$var <- factor(D$var, levels = unname(PD_LAB[vars]))
  ggplot(D, aes(x, y)) +
    geom_line(linewidth = 0.9, colour = col) +
    facet_wrap(~var, scales = "free_x", nrow = 2) +
    labs(title = ttl, x = NULL,
         y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 7.4),
          plot.title = element_text(face = "bold", size = 10))
}
pe_ <- pd_grid(TreeRF, Xt, c("soil_moisture_at_tree", "dbh_m",
                             "height_cm", "soil_temp_C_mean"),
               "(e) Tree partial dependence", TREE_COL)
pf_ <- pd_grid(SoilRF, Xs, c("soil_moisture_at_site", "soil_temp_C_mean",
                             "air_temp_C_mean", "month"),
               "(f) Soil partial dependence", SOIL_COL)

fig <- (pa_ | pc_ | pe_) / (pb_ | pd_ | pf_) +
  plot_layout(widths = c(0.85, 0.9, 1.15)) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))
ggsave(out_path("figS21_rf_model_summary.png"), fig,
       width = 15, height = 8.4, dpi = 200, bg = "white")

cat("=== grouped permutation importance (held-out), share within model ===\n")
for (M in list(IMP_T, IMP_S)) {
  cat(sprintf("\n--- %s ---\n", M$model[1]))
  print(as.data.frame(M %>% arrange(-share) %>%
        transmute(predictor = label, share = sprintf("%.1f%%", 100 * share))),
        row.names = FALSE)
}
cat("\nwritten: outputs/figures/generated/figS21_rf_model_summary.png\n")
