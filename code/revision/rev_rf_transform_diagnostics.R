#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_transform_diagnostics.R
# ------------------------------------------------------------------------------
# Side-by-side diagnostics for the two candidate target treatments, so the choice
# is made on behaviour rather than on a single summary number:
#
#   asinh + ratio   fit on asinh(flux); back-transform with sinh() and a single
#                   multiplicative factor mean(obs)/mean(sinh(OOB pred)) estimated
#                   from training data only
#   raw             fit on flux directly; no back-transformation at all
#
# All predictions are GROUPED 5-fold CV (whole trees / whole collars held out), which
# is what upscaling to unmeasured inventory stems requires. OOB would flatter both.
#
# Panels per model:
#   1. predicted vs observed (pseudo-log axes, 1:1 line)
#   2. residual vs fitted -- is error structure flat, or does it fan with magnitude?
#   3. observed vs predicted MEAN by observed decile -- the quantity a budget needs
#   4. calibration of the cumulative total -- how much of the true flux is recovered
#
# Output: outputs/revision/fig_transform_diagnostics.png
#         outputs/revision/transform_diagnostics_summary.csv
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(ggplot2); library(tidyr)})
set.seed(42)
load("outputs/models/TRAINING_DATA.RData")
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

grouped_pred <- function(X, yraw, grp, mode) {
  X <- as.data.frame(X); k <- is.finite(yraw)
  X <- X[k, , drop = FALSE]; yraw <- yraw[k]; grp <- grp[k]
  ut <- unique(grp); set.seed(1)
  gm <- setNames(sample(rep(1:5, length.out = length(ut))), ut); fold <- unname(gm[grp])
  pr <- rep(NA_real_, length(yraw))
  for (kk in 1:5) {
    tr <- which(fold != kk); te <- which(fold == kk)
    yy <- if (mode == "raw") yraw else asinh(yraw)
    m <- ranger(x = X[tr, , drop = FALSE], y = yy[tr], num.trees = 800, min.node.size = 5,
                mtry = max(1, floor(sqrt(ncol(X)))), num.threads = 1, seed = 42)
    p <- predict(m, X[te, , drop = FALSE])$predictions
    if (mode == "raw") pr[te] <- p
    else { ins <- sinh(m$predictions)
           f <- mean(yraw[tr]) / mean(ins[is.finite(ins)])
           pr[te] <- sinh(p) * f }
  }
  data.frame(obs = yraw, pred = pr)
}

SETS <- list(
  list(nm = "Tree stem", X = X_tree, y = tree_train_complete$stem_flux_corrected,
       g = tree_train_complete$tree_id),
  list(nm = "Soil",      X = X_soil, y = soil_train_complete$soil_flux_umol_m2_s,
       g = soil_train_complete$site_id)
)

all_df <- bind_rows(lapply(SETS, function(S)
  bind_rows(lapply(c("asinh + ratio", "raw"), function(md) {
    d <- grouped_pred(S$X, S$y, S$g, if (md == "raw") "raw" else "asinh_ratio")
    d$model <- S$nm; d$method <- md; d
  }))))

summ <- all_df %>% group_by(model, method) %>%
  summarise(n = n(),
            r2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
            spearman = cor(obs, pred, method = "spearman"),
            mean_ratio = mean(pred) / mean(obs),
            tail_ratio = { q <- abs(obs) >= quantile(abs(obs), .9)
                           mean(pred[q]) / mean(obs[q]) },
            .groups = "drop")
cat("=== grouped-CV summary ===\n"); print(as.data.frame(summ), row.names = FALSE, digits = 3)
write.csv(summ, file.path(outdir, "transform_diagnostics_summary.csv"), row.names = FALSE)

p1 <- ggplot(all_df, aes(obs, pred, colour = method)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  facet_grid(method ~ model, scales = "free") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_colour_manual(values = c("asinh + ratio" = "#B2182B", "raw" = "#2166AC"), guide = "none") +
  labs(title = "Predicted vs observed (grouped CV)", x = "observed", y = "predicted") +
  theme_minimal(base_size = 8)

p2 <- ggplot(all_df, aes(pred, obs - pred, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey30") +
  geom_point(alpha = 0.3, size = 0.8) +
  facet_grid(method ~ model, scales = "free") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_colour_manual(values = c("asinh + ratio" = "#B2182B", "raw" = "#2166AC"), guide = "none") +
  labs(title = "Residual vs fitted", x = "predicted", y = "observed - predicted") +
  theme_minimal(base_size = 8)

dec <- all_df %>% group_by(model, method) %>%
  mutate(d = ntile(obs, 10)) %>% group_by(model, method, d) %>%
  summarise(obs = mean(obs), pred = mean(pred), .groups = "drop")
p3 <- ggplot(dec, aes(obs, pred, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  geom_point(size = 1.6) + geom_line(linewidth = 0.5) +
  facet_wrap(~ model, scales = "free") +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  scale_colour_manual(values = c("asinh + ratio" = "#B2182B", "raw" = "#2166AC"), name = NULL) +
  labs(title = "Decile means: the quantity a budget depends on",
       x = "mean observed", y = "mean predicted") +
  theme_minimal(base_size = 8) + theme(legend.position = "top")

cum <- all_df %>% group_by(model, method) %>% arrange(desc(obs), .by_group = TRUE) %>%
  mutate(rank = row_number() / n(),
         cum_obs = cumsum(obs) / sum(obs), cum_pred = cumsum(pred) / sum(obs))
p4 <- ggplot(cum, aes(rank, colour = method)) +
  geom_line(aes(y = cum_pred), linewidth = 0.7) +
  geom_line(aes(y = cum_obs), colour = "grey40", linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ model, scales = "free_y") +
  scale_colour_manual(values = c("asinh + ratio" = "#B2182B", "raw" = "#2166AC"), guide = "none") +
  labs(title = "Cumulative flux recovered (grey dashed = observed)",
       x = "fraction of measurements, largest first", y = "cumulative / total observed") +
  theme_minimal(base_size = 8)

ggsave(file.path(outdir, "fig_transform_diagnostics.png"),
       gridExtra::arrangeGrob(p1, p2, p3, p4, layout_matrix = rbind(c(1,2), c(3,4))),
       width = 11, height = 8.5, dpi = 200, bg = "white")
cat("\nWritten: outputs/revision/fig_transform_diagnostics.png\n")
