# ==============================================================================
# REVISION — Does the flux-decline-with-height correspond to soil moisture?
# ==============================================================================
# Referee 3's logical-gap challenge: "Reduction in stem methane flux with height
# does not necessarily imply a soil-source... A trend with height for B.
# alleghaniensis is used as evidence that wetter soil dominates... Are the highest
# B. alleghaniensis fluxes found in the wettest soil? It doesn't appear so from
# Figure 3."
#
# Logic being tested: if height-decline reflects soil-derived CH4 transported up
# and attenuating with height, then a STEEPER (more negative) height slope should
# occur in WETTER soils (more soil CH4 source). We test per-tree height slope vs VWC.
#
# NEW file; edits nothing. Run: Rscript code/revision/R_height_slope_vs_moisture.R
# Outputs: outputs/revision/height_slope_vs_moisture.png, height_slope_moisture_summary.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "outputs/revision"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

tp <- read_csv("data/compiled/tree_properties.csv", show_col_types = FALSE)

# long form: height (cm) -> flux; RootCrown treated as 0 cm
hmap <- c(CH4_best.flux_RootCrowncm = 0, CH4_best.flux_50cm = 50,
          CH4_best.flux_75cm = 75, CH4_best.flux_125cm = 125, CH4_best.flux_200cm = 200)
long <- tp %>% dplyr::select(tree_id, species_id, VWC_mean, landscape_position, names(hmap)) %>%
  pivot_longer(names(hmap), names_to = "col", values_to = "flux") %>%
  mutate(height_cm = hmap[col]) %>% filter(!is.na(flux))

# per-tree height slope (need >=3 heights)
slopes <- long %>% group_by(tree_id, species_id, VWC_mean, landscape_position) %>%
  filter(n_distinct(height_cm) >= 3) %>%
  summarise(height_slope = coef(lm(flux ~ height_cm))[2],
            flux_base = flux[which.min(height_cm)][1], .groups = "drop") %>%
  mutate(species = sp_map[species_id]) %>% filter(!is.na(VWC_mean))

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
                     labels = c("other", "B. alleghaniensis"), name = NULL) +
  annotate("text", -Inf, Inf, hjust = -0.05, vjust = 1.4, size = 3,
           label = sprintf("r = %.2f, %s, n = %d", ct_all$estimate,
                           ifelse(ct_all$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", ct_all$p.value)), nrow(slopes))) +
  labs(x = "Soil VWC (%)", y = "Height slope (flux change per cm)") + theme_bw(base_size = 10)

# R3's specific birch check: birch flux (base height) vs VWC
birch <- slopes %>% filter(species == "Betula alleghaniensis")
p2 <- ggplot(birch, aes(VWC_mean, flux_base)) +
  geom_point(size = 2.5, color = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  annotate("text", -Inf, Inf, hjust = -0.05, vjust = 1.4, size = 3,
           label = if (nrow(birch) > 2) sprintf("r = %.2f, p = %.3f, n = %d",
                    cor(birch$VWC_mean, birch$flux_base), cor.test(birch$VWC_mean, birch$flux_base)$p.value, nrow(birch)) else "n too small") +
  labs(x = "Soil VWC (%)", y = "Base-height CH4 flux") + theme_bw(base_size = 10)

if (requireNamespace("gridExtra", quietly = TRUE))
  ggsave(file.path(out_dir, "height_slope_vs_moisture.png"),
         gridExtra::arrangeGrob(p1, p2, nrow = 1), width = 12, height = 5, dpi = 150)

# ---- Summary -----------------------------------------------------------------
sink(file.path(out_dir, "height_slope_moisture_summary.txt"))
cat("=================================================================\n")
cat("HEIGHT-DECLINE SLOPE vs SOIL MOISTURE (Referee 3's logical gap)\n")
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
cat("\nB. alleghaniensis wettest-soil check (R3):\n")
if (nrow(birch) > 2) {
  bc <- cor.test(birch$VWC_mean, birch$flux_base)
  cat(sprintf("  r(base flux, VWC) = %.2f, p = %.3f, n = %d\n", bc$estimate, bc$p.value, nrow(birch)))
  cat("  -> if positive & significant, highest birch fluxes ARE in wettest soil (supports\n")
  cat("     transport); if not, R3's skepticism stands and the claim needs softening.\n")
} else cat("  n too small for a birch-only correlation.\n")
cat("\nHONEST READ: report whatever the sign/significance actually is. If the per-tree\n")
cat("slope-vs-VWC correlation is weak, concede R3's point that height-decline alone does\n")
cat("not establish a soil source, and lean on the cross-species height-effect vs soil-\n")
cat("methanogen correlation (r=-0.76, manuscript L302-304) instead.\n")
sink()
cat(readLines(file.path(out_dir,"height_slope_moisture_summary.txt")), sep="\n")
