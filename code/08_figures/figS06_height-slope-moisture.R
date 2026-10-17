source("code/lib/outputs.R")
# ==============================================================================
# Does the flux-decline-with-height correspond to soil moisture?
# ==============================================================================
# A decline in stem CH4 flux with height does not by itself imply a soil source.
# Tests whether the B. alleghaniensis height trend reflects wetter soil: are the
# highest B. alleghaniensis fluxes found in the wettest soil?
#
# Logic being tested: if height-decline reflects soil-derived CH4 transported up
# and attenuating with height, then a STEEPER (more negative) height slope should
# occur in WETTER soils (more soil CH4 source). We test per-tree height slope vs VWC.
#
# Outputs: outputs/figures/generated/height_slope_vs_moisture.png, height_slope_moisture_summary.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
out_dir <- "outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


source("code/lib/tables_read.R")
# 2021 survey per tree: flux by height (flux_stem.csv), species (trees.csv), landscape
# position, soil VWC at the tree (attach_soil(), code/lib/soil_fill.R)
tp <- survey_2021_trees(soil = TRUE) %>% mutate(landscape_position = plot)

# long form: height (cm) -> flux; RootCrown treated as 0 cm
hmap <- c(CH4_best.flux_RootCrowncm = 0, CH4_best.flux_50cm = 50,
          CH4_best.flux_75cm = 75, CH4_best.flux_125cm = 125, CH4_best.flux_200cm = 200)
long <- tp %>% dplyr::select(tree_id, species_id, VWC_mean, landscape_position, names(hmap)) %>%
  pivot_longer(names(hmap), names_to = "col", values_to = "flux") %>%
  mutate(height_cm = hmap[col]) %>% filter(!is.na(flux))

# per-tree height slope (need >=3 heights)
slopes <- long %>% group_by(tree_id, species_id, VWC_mean, landscape_position) %>%
  filter(n_distinct(height_cm) >= 3) %>%
  summarise(height_slope = 100 * coef(lm(flux ~ height_cm))[2],   # per m, as Fig. 2b
            flux_base = flux[which.min(height_cm)][1], .groups = "drop") %>%
  mutate(species = species_name(species_id)) %>% filter(!is.na(VWC_mean))

# ---- Correlations: per-tree slope vs VWC -------------------------------------
ct_all  <- cor.test(slopes$height_slope, slopes$VWC_mean, method = "pearson")
ct_sp_r <- cor.test(rank(slopes$height_slope), rank(slopes$VWC_mean))  # spearman-ish

# species-level: median slope vs median VWC
sp_lvl <- slopes %>% group_by(species) %>%
  filter(n() >= 3) %>%
  summarise(n = n(), med_slope = median(height_slope), med_vwc = median(VWC_mean), .groups = "drop")
ct_species <- if (nrow(sp_lvl) > 3) cor.test(sp_lvl$med_slope, sp_lvl$med_vwc) else NULL

# ---- Figures -----------------------------------------------------------------
p1 <- ggplot(slopes, aes(VWC_mean, height_slope)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(aes(color = species == "Betula alleghaniensis"), alpha = .7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  scale_color_manual(values = c(`FALSE` = "gray50", `TRUE` = "firebrick"),
                     labels = c("Other species", expression(italic("B. alleghaniensis"))), name = NULL) +
  labs(x = "Soil VWC (%)", y = expression("Height slope (nmol m"^-2*" s"^-1*" per m)")) +
  theme_bw(base_size = 12)

# Birch check: B. alleghaniensis flux (base height) vs VWC
birch <- slopes %>% filter(species == "Betula alleghaniensis")
p2 <- ggplot(birch, aes(VWC_mean, flux_base)) +
  geom_point(size = 2.5, color = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  labs(x = "Soil VWC (%)", y = expression(CH[4]~"flux at the lowest height (nmol m"^-2*" s"^-1*")")) +
  theme_bw(base_size = 12)

# statistics (r, p, n for each panel) are in the SI caption, from the summary below
ggsave(out_path("height_slope_vs_moisture.png"),
       (p1 | p2) + patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
         theme(plot.tag = element_text(size = 14, face = "bold")),
       width = 12, height = 5, dpi = 300)

# ---- Summary -----------------------------------------------------------------
sink(out_path("height_slope_moisture_summary.txt"))
cat("=================================================================\n")
cat("HEIGHT-DECLINE SLOPE vs SOIL MOISTURE\n")
cat("=================================================================\n\n")
cat(sprintf("Trees with >=3 heights AND VWC: %d\n\n", nrow(slopes)))
cat("Hypothesis: if height-decline = soil-transported CH4, wetter soils (more soil\n")
cat("source) should show STEEPER (more negative) height slopes -> NEGATIVE r between\n")
cat("height_slope and VWC.\n\n")
cat(sprintf("Per-tree:  Pearson r(height_slope, VWC) = %.3f, p = %.4f\n",
            ct_all$estimate, ct_all$p.value))
if (!is.null(ct_species))
  cat(sprintf("Species:   Pearson r(median slope, median VWC) = %.3f, p = %.4f (n=%d species)\n",
              ct_species$estimate, ct_species$p.value, nrow(sp_lvl)))
cat("\nSpecies-level slope vs VWC:\n"); print(as.data.frame(sp_lvl), row.names = FALSE)
cat("\nB. alleghaniensis wettest-soil check:\n")
if (nrow(birch) > 2) {
  bc <- cor.test(birch$VWC_mean, birch$flux_base)
  cat(sprintf("  r(base flux, VWC) = %.2f, p = %.3f, n = %d\n", bc$estimate, bc$p.value, nrow(birch)))
  cat("  -> if positive & significant, highest birch fluxes ARE in wettest soil (supports\n")
  cat("     transport); if not, this check does not support a soil source.\n")
} else cat("  n too small for a birch-only correlation.\n")
cat("\nINTERPRETATION: if the per-tree slope-vs-VWC correlation is weak, height-decline\n")
cat("alone does not establish a soil source; the cross-species height-effect vs soil-\n")
cat("methanogen correlation (r=-0.76) is the stronger evidence.\n")
sink()
cat(readLines(out_path("height_slope_moisture_summary.txt")), sep="\n")
