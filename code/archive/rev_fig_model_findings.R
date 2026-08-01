source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_model_findings.R  --  main-text panel: what the models actually learned
# ------------------------------------------------------------------------------
# Diagnostics belong in the SI; this is the one panel that carries a RESULT.
#
# IMPORTANCE IS PERMUTATION (incMSE), NOT IMPURITY, and the distinction is not
# cosmetic here -- the two metrics give opposite answers. Impurity gain counts
# total variance reduction at splits on a predictor and is well known to favour
# continuous, high-cardinality variables; permutation importance measures the
# increase in out-of-bag error when a predictor is shuffled, which is what
# "does the model rely on this" actually means.
#
#   predictor              impurity %   permutation %
#   TreeRF  diameter             9.0        29.6
#           diameter (within)    8.8        24.2
#           species              5.8        16.0
#           soil moisture       67.2        11.0
#   SoilRF  soil moisture        9.0        44.0
#           air temperature      2.1        22.7
#           soil temperature    88.0        18.7
#
# Impurity said the soil model was a temperature model and the tree model was a
# moisture model. Permutation says the reverse of both. Permutation also
# reproduces the earlier reference-column test in 02_rf_models.R (moisture 1.167
# here against 1.121 there), so it is the metric with two independent
# confirmations behind it.
#
# What the panel therefore reports: the SOIL model is driven by moisture, as
# expected for diffusion-limited uptake; and the TREE model is driven mostly by
# SIZE -- diameter and within-species diameter together are 54% -- with species
# third at 16% and moisture fourth at 11%.
#
#   a  predictor importance, both models, permutation (incMSE)
#   b  partial dependence on the dominant predictor, with the observed range
#   c  predicted against observed, calibrated, with the 1:1 line
#
# Output: outputs/figures/generated/fig_model_findings.png
# ==============================================================================
suppressMessages({library(dplyr); library(ranger); library(ggplot2); library(patchwork); library(tidyr)})
set.seed(42)
outdir <- "outputs"
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
FLUX_PURPLE <- "#756BB1"

nice <- c(soil_moisture_at_tree = "soil moisture", soil_moisture_at_site = "soil moisture",
          dbh_m = "diameter", dbh_within_z = "diameter (within species)",
          species = "species", height_cm = "measurement height",
          soil_temp_C_mean = "soil temperature", air_temp_C_mean = "air temperature",
          month = "month")

# GROUPED permutation importance, measured on held-out data.
#
# Ordinary permutation importance is not readable when predictors are collinear:
# shuffling one leaves the other as a proxy, so both look unimportant. Soil and
# air temperature correlate at r = 0.97 here, which splits one physical driver
# across two mutually-suppressed bars. They are therefore shuffled TOGETHER with
# one shared permutation, which keeps the group's joint structure intact and asks
# the only question that has an answer: how much does the model rely on
# temperature at all.
#
# Measured out-of-sample rather than in-sample -- fit on 70%, permute and score
# on the remaining 30% -- so baseline and permuted error are on the same footing.
# (An earlier version of this compared an OOB baseline against in-sample permuted
# predictions, which produced negative importances.)
GROUPS <- list(
  TreeRF = list(species = "species", diameter = "dbh_m",
                `soil moisture` = "soil_moisture_at_tree",
                temperature = c("soil_temp_C_mean", "air_temp_C_mean"),
                `measurement height` = "height_cm"),
  SoilRF = list(`soil moisture` = "soil_moisture_at_site",
                temperature = c("soil_temp_C_mean", "air_temp_C_mean"),
                month = "month"))

gimp <- function(X, y, groups, lab, nsplit = 10, nrep = 20) {
  k <- is.finite(y) & complete.cases(X); X <- X[k,,drop=FALSE]; y <- y[k]
  R <- replicate(nsplit, {
    tr <- sample(nrow(X), floor(0.7*nrow(X)))
    m  <- ranger(x = X[tr,,drop=FALSE], y = y[tr], num.trees = 800, min.node.size = 5,
                 mtry = max(1, floor(sqrt(ncol(X)))), num.threads = 1)
    Xh <- X[-tr,,drop=FALSE]; yh <- y[-tr]
    base <- mean((yh - predict(m, Xh, num.threads = 1)$predictions)^2)
    sapply(groups, function(g) mean(replicate(nrep, {
      Xp <- Xh; i <- sample(nrow(Xp))
      for (v in g) Xp[[v]] <- Xp[[v]][i]        # ONE shared permutation per group
      mean((yh - predict(m, Xp, num.threads = 1)$predictions)^2) - base
    })))
  })
  mu <- pmax(rowMeans(R), 0)
  data.frame(model = lab, label = names(groups), imp = as.numeric(mu),
             se = apply(R, 1, sd)/sqrt(nsplit)) %>%
    mutate(share = imp/sum(imp), share_se = se/sum(imp))
}
Xt <- as.data.frame(X_tree); yt <- tree_train_complete$stem_flux_corrected
Xs <- as.data.frame(X_soil); ys <- soil_train_complete$y_asinh[seq_len(nrow(Xs))]
IMP_T <- gimp(Xt, yt, GROUPS$TreeRF, "TreeRF")
IMP_S <- gimp(Xs, ys, GROUPS$SoilRF, "SoilRF")
TreePerm <- TreeRF; SoilPerm <- SoilRF
# "soil moisture" appears in both models, so a single shared factor cannot order
# both facets. Pad the SoilRF labels with a trailing space to make them distinct,
# then order the combined levels by model and share.
IMP <- bind_rows(IMP_T, IMP_S) %>%
  mutate(label = ifelse(model == "SoilRF", paste0(label, " "), label)) %>%
  arrange(model, share)
IMP$label <- factor(IMP$label, levels = IMP$label)

# Error bars are +/-1 SE over the 10 train/test splits. They are on the panel
# because the ranking is NOT stable: the same forest scored by ranger's OOB
# permutation instead of this held-out one moves soil moisture from 36% to 28%
# and temperature from 11% to 28%. Grouping is a small correction by comparison
# (temperature 7.8% -> 11.0%, no reordering); the estimator is what moves things.
# The panel therefore supports "no single predictor dominates, and species and
# moisture are both major terms" -- not a precise ordering.
pa <- ggplot(IMP, aes(share, label, fill = model)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_errorbarh(aes(xmin = pmax(0, share - share_se), xmax = share + share_se),
                 height = 0.22, linewidth = 0.35, colour = "grey25") +
  geom_text(aes(x = share + share_se, label = sprintf("%.0f%%", 100*share)),
            hjust = -0.3, size = 2.9) +
  facet_wrap(~model, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(TreeRF = FLUX_PURPLE, SoilRF = "#2166ac")) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.16), expand = c(0, 0)) +
  labs(title = "a   what each model relies on", x = "share of total importance", y = NULL,
       subtitle = "grouped permutation importance, held-out (+/-1 SE);\nsoil and air temperature pooled (r = 0.97)") +
  theme_bw(base_size = 9) +
  theme(panel.grid.major.y = element_blank(), strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.2, colour = "grey30"))

# ---- b: partial dependence on moisture, both models --------------------------
pdp <- function(m, dat, xv, n = 40) {
  g <- seq(quantile(dat[[xv]], .02, na.rm = TRUE), quantile(dat[[xv]], .98, na.rm = TRUE), length.out = n)
  data.frame(x = g, y = sapply(g, function(v) { D <- dat; D[[xv]] <- v
    mean(predict(m, D, num.threads = 1)$predictions) }))
}
PT <- pdp(TreeRF, Xt, "soil_moisture_at_tree") %>% mutate(model = "TreeRF")
PS <- pdp(SoilRF, Xs, "soil_moisture_at_site") %>% mutate(model = "SoilRF")
PD <- bind_rows(PT, PS)
rug <- bind_rows(data.frame(x = Xt$soil_moisture_at_tree, model = "TreeRF"),
                 data.frame(x = Xs$soil_moisture_at_site, model = "SoilRF"))
pb <- ggplot(PD, aes(x, y, colour = model)) +
  geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, alpha = 0.16, length = unit(0.02, "npc")) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 1) +
  facet_wrap(~model, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(TreeRF = FLUX_PURPLE, SoilRF = "#2166ac"), guide = "none") +
  labs(title = "b   response to soil moisture",
       subtitle = "partial dependence; ticks are the observed range",
       x = expression("soil moisture (m"^3*" m"^-3*")"),
       y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.2, colour = "grey30"))

# ---- c: predicted vs observed, calibrated -----------------------------------
# Out-of-bag predictions with the per-species calibration applied, which is what
# the budget uses. Both axes on a signed log scale: the response spans four
# orders of magnitude and includes negatives, so a linear axis shows one point.
d <- tree_train_complete
cal <- data.frame(sp = as.character(d$species_clean), o = d$stem_flux_corrected,
                  p = TreeRF$predictions) %>% filter(is.finite(o), is.finite(p)) %>%
  # NO CLAMP. This was pmax(pmin(ratio, 5), 0.2), but rev_predict_tree_flux_current.R
  # removed the clamp, so the main-text panel was showing a calibration the stand
  # budget does not use.
  group_by(sp) %>% mutate(r = ifelse(mean(p) > 0, mean(o)/mean(p), 1)) %>%
  ungroup() %>% mutate(pc = p*r)
slog <- function(x) sign(x)*log10(1 + abs(x)/0.01)
r2 <- 1 - sum((cal$o - cal$pc)^2)/sum((cal$o - mean(cal$o))^2)
# Grouped-CV skill alongside the OOB figure. OOB is not grouped by tree, and 1,130
# rows come from 478 trees, so quoting OOB alone in a main-text panel overstates the
# model by ~2x. rev_rf_grouped_cv.R writes this.
R2G <- local({
  f <- "outputs/data/rf_grouped_cv.csv"
  if (!file.exists(f)) return(NA_real_)
  g <- read.csv(f, stringsAsFactors = FALSE); g$grouped_r2[g$model == "TreeRF"][1]
})
pc <- ggplot(cal, aes(slog(pc), slog(o))) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
  geom_point(alpha = 0.3, size = 1.1, colour = FLUX_PURPLE) +
  labs(title = "c   predicted against observed",
       subtitle = sprintf("out-of-bag, calibrated; %d measurements from %d trees\nR2 %.2f out-of-bag, %.2f cross-validated by tree\nsigned log axes: flux spans 4 orders of\nmagnitude and includes uptake",
                          nrow(cal), length(unique(tree_train_complete$tree_id)), r2, R2G),
       x = "predicted", y = "observed") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.2, colour = "grey30"))

fig <- (pa | pb | pc) + plot_layout(widths = c(1, 0.95, 0.85)) +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
ggsave(out_path("fig_model_findings.png"), fig, width = 11, height = 5.2, dpi = 200, bg = "white")

cat("=== grouped permutation importance (held-out), share within model ===\n")
cat("=== TreeRF ===\n")
print(as.data.frame(IMP %>% filter(model == "TreeRF") %>% arrange(-share) %>%
      transmute(predictor = label, share = sprintf("%.1f%%", 100*share))), row.names = FALSE)
cat("\n=== SoilRF importance (share) ===\n")
print(as.data.frame(IMP %>% filter(model == "SoilRF") %>% arrange(-share) %>%
      transmute(predictor = label, share = sprintf("%.1f%%", 100*share))), row.names = FALSE)
cat(sprintf("\ncalibrated OOB R2 (tree): %.3f on %d measurements\n", r2, nrow(cal)))
cat("\nwritten: outputs/figures/generated/fig_model_findings.png\n")
