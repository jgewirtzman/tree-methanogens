#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_model_findings.R  --  main-text panel: what the models actually learned
# ------------------------------------------------------------------------------
# Diagnostics belong in the SI; this is the one panel that carries a RESULT. The
# result is that the tree model is an environment model, not a species model:
# soil moisture dominates it, species is a minor term, and that single fact
# explains several things the paper otherwise has to assert separately --
#   * why the species work (ladder, per-species levels) moved the budget so little
#   * why the Figure 9 tree map should not be read as a species map
#   * why a stand-level upscaling built on this model is more sensitive to where
#     the wet ground is than to what is growing on it
#
#   a  predictor importance, both models, as a share of total impurity gain
#   b  partial dependence on the dominant predictor, with the observed range
#   c  predicted against observed, calibrated, with the 1:1 line
#
# Output: outputs/revision/fig_model_findings.png
# ==============================================================================
suppressMessages({library(dplyr); library(ranger); library(ggplot2); library(patchwork); library(tidyr)})
set.seed(42)
outdir <- "outputs/revision"
load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
FLUX_PURPLE <- "#756BB1"

nice <- c(soil_moisture_at_tree = "soil moisture", soil_moisture_at_site = "soil moisture",
          dbh_m = "diameter", dbh_within_z = "diameter (within species)",
          species = "species", height_cm = "measurement height",
          soil_temp_C_mean = "soil temperature", air_temp_C_mean = "air temperature",
          month = "month")

imp <- function(m, lab) {
  v <- m$variable.importance
  data.frame(model = lab, var = names(v), imp = as.numeric(v)) %>%
    mutate(share = imp/sum(imp), label = ifelse(is.na(nice[var]), var, nice[var]))
}
IMP <- bind_rows(imp(TreeRF, "TreeRF"), imp(SoilRF, "SoilRF")) %>%
  group_by(model) %>% arrange(model, share) %>% ungroup()
IMP$label <- factor(IMP$label, levels = unique(IMP$label[order(IMP$share)]))

pa <- ggplot(IMP, aes(share, label, fill = model)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.0f%%", 100*share)), hjust = -0.15, size = 2.9) +
  facet_wrap(~model, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(TreeRF = FLUX_PURPLE, SoilRF = "#2166ac")) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.16), expand = c(0, 0)) +
  labs(title = "a   what each model uses", x = "share of total importance", y = NULL,
       subtitle = "impurity gain, normalised within model") +
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
Xt <- as.data.frame(X_tree); Xs <- as.data.frame(X_soil)
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
  group_by(sp) %>% mutate(r = pmax(pmin(ifelse(mean(p) > 0, mean(o)/mean(p), 1), 5), 0.2)) %>%
  ungroup() %>% mutate(pc = p*r)
slog <- function(x) sign(x)*log10(1 + abs(x)/0.01)
r2 <- 1 - sum((cal$o - cal$pc)^2)/sum((cal$o - mean(cal$o))^2)
pc <- ggplot(cal, aes(slog(pc), slog(o))) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
  geom_point(alpha = 0.3, size = 1.1, colour = FLUX_PURPLE) +
  labs(title = "c   predicted against observed",
       subtitle = sprintf("out-of-bag, calibrated; %d measurements, R2 %.2f\nsigned log axes (flux spans 4 orders of magnitude and includes uptake)",
                          nrow(cal), r2),
       x = "predicted", y = "observed") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.2, colour = "grey30"))

fig <- (pa | pb | pc) + plot_layout(widths = c(1, 1, 0.9)) +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
ggsave(file.path(outdir, "fig_model_findings.png"), fig, width = 11, height = 5.2, dpi = 200, bg = "white")

cat("=== TreeRF importance (share) ===\n")
print(as.data.frame(IMP %>% filter(model == "TreeRF") %>% arrange(-share) %>%
      transmute(predictor = label, share = sprintf("%.1f%%", 100*share))), row.names = FALSE)
cat("\n=== SoilRF importance (share) ===\n")
print(as.data.frame(IMP %>% filter(model == "SoilRF") %>% arrange(-share) %>%
      transmute(predictor = label, share = sprintf("%.1f%%", 100*share))), row.names = FALSE)
cat(sprintf("\ncalibrated OOB R2 (tree): %.3f on %d measurements\n", r2, nrow(cal)))
cat("\nwritten: outputs/revision/fig_model_findings.png\n")
