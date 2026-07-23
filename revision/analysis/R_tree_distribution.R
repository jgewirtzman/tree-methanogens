# ==============================================================================
# REVISION — Tree species distribution: DBH x landscape table + moisture niches
# ==============================================================================
# Serves: R3 L286 (table of #trees/species + DBH mean+/-sd by location),
# R3 L125 (species distribution across the transect), and the editor + R3 + Mark
# plant-framing ask (species niche preferences, e.g. yellow birch in wetter soils).
#
# NEW file; edits nothing.
# Inputs: data/compiled/tree_properties.csv (measured trees: species, dbh, VWC,
#         landscape_position), data/compiled/forest_inventory.csv (full census).
# Run: Rscript revision/analysis/R_tree_distribution.R
# Outputs (revision/outputs/):
#   tree_dbh_landscape_table.csv     R3 L286 table
#   tree_species_moisture_niche.csv  median VWC per species (niche summary)
#   tree_distribution_panels.png     DBH-by-landscape + species-VWC niche figure
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "revision/outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")
gymno <- c("Pinus strobus","Tsuga canadensis")

tp <- read_csv("data/compiled/tree_properties.csv", show_col_types = FALSE) %>%
  filter(!is.na(species_id)) %>%
  mutate(species = sp_map[species_id],
         group = if_else(species %in% gymno, "Gymnosperm", "Angiosperm"),
         landscape_position = na_if(landscape_position, "Lab")) %>%
  filter(!is.na(species))

# ---- R3 L286: species x landscape count + DBH mean+/-sd -----------------------
dbh_tab <- tp %>% filter(!is.na(landscape_position)) %>%
  group_by(species, group, landscape_position) %>%
  summarise(n = n(),
            dbh_mean = mean(dbh, na.rm = TRUE), dbh_sd = sd(dbh, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(dbh_summary = ifelse(is.na(dbh_sd), sprintf("%.1f", dbh_mean),
                              sprintf("%.1f +/- %.1f", dbh_mean, dbh_sd))) %>%
  arrange(group, species, landscape_position)
write.csv(dbh_tab, file.path(out_dir, "tree_dbh_landscape_table.csv"), row.names = FALSE)

# wide version (species rows, landscape columns) — closer to a manuscript table
dbh_wide <- dbh_tab %>%
  mutate(cell = sprintf("n=%d, DBH %s", n, dbh_summary)) %>%
  dplyr::select(species, group, landscape_position, cell) %>%
  pivot_wider(names_from = landscape_position, values_from = cell)
write.csv(dbh_wide, file.path(out_dir, "tree_dbh_landscape_table_wide.csv"), row.names = FALSE)

# ---- Species moisture niche (VWC) --------------------------------------------
niche <- tp %>% filter(!is.na(VWC_mean)) %>%
  group_by(species, group) %>%
  summarise(n = n(), vwc_median = median(VWC_mean), vwc_mean = mean(VWC_mean),
            vwc_q25 = quantile(VWC_mean, .25), vwc_q75 = quantile(VWC_mean, .75),
            .groups = "drop") %>%
  arrange(desc(vwc_median))
write.csv(niche, file.path(out_dir, "tree_species_moisture_niche.csv"), row.names = FALSE)

# ---- Figures -----------------------------------------------------------------
sp_order <- niche %>% arrange(vwc_median) %>% pull(species)
tp_v <- tp %>% filter(!is.na(VWC_mean)) %>% mutate(species = factor(species, levels = sp_order))

p_niche <- ggplot(tp_v, aes(VWC_mean, species, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .5) +
  geom_jitter(height = .15, size = 1, alpha = .5) +
  scale_color_manual(values = c(Angiosperm = "#1b7837", Gymnosperm = "#762a83")) +
  labs(title = "Species soil-moisture niches (ordered by median VWC)",
       subtitle = "Wetter-preferring species at top (e.g. yellow birch); gymnosperms on drier ridges",
       x = "Soil VWC (%)", y = NULL, color = NULL) +
  theme_bw(base_size = 10) + theme(axis.text.y = element_text(face = "italic"))

comp <- tp %>% filter(!is.na(landscape_position)) %>%
  count(landscape_position, group, species)
p_comp <- ggplot(comp, aes(landscape_position, n, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c(Angiosperm = "#1b7837", Gymnosperm = "#762a83")) +
  labs(title = "Sampled-tree composition by landscape position",
       x = NULL, y = "n trees", fill = NULL) + theme_bw(base_size = 10)

p_dbh <- ggplot(tp %>% filter(!is.na(dbh), !is.na(landscape_position)),
                aes(landscape_position, dbh, color = group)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 1, alpha = .4) +
  scale_color_manual(values = c(Angiosperm = "#1b7837", Gymnosperm = "#762a83")) +
  labs(title = "DBH by landscape position", x = NULL, y = "DBH (cm)", color = NULL) +
  theme_bw(base_size = 10)

if (requireNamespace("gridExtra", quietly = TRUE)) {
  g <- gridExtra::arrangeGrob(p_niche, gridExtra::arrangeGrob(p_comp, p_dbh, ncol = 1), ncol = 2,
                              widths = c(1.3, 1))
  ggsave(file.path(out_dir, "tree_distribution_panels.png"), g, width = 13, height = 7, dpi = 140)
}

# ---- Forest-inventory context (full census) ----------------------------------
fi <- read_csv("data/compiled/forest_inventory.csv", show_col_types = FALSE)
fi_summary <- fi %>% mutate(species = sp_map[species_code]) %>%
  group_by(species_code) %>%
  summarise(n_stems = n(), dbh_mean_mm = mean(dbh_mm, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(n_stems))
write.csv(fi_summary, file.path(out_dir, "forest_inventory_species_summary.csv"), row.names = FALSE)

# ---- Console summary ---------------------------------------------------------
cat("=== R3 L286 species x landscape x DBH table (head) ===\n")
print(as.data.frame(head(dbh_tab, 12)), row.names = FALSE)
cat("\n=== Species moisture niches (VWC median, wettest first) ===\n")
print(as.data.frame(niche %>% dplyr::select(species, group, n, vwc_median)), row.names = FALSE)
cat("\nWrote tables + tree_distribution_panels.png to", out_dir, "\n")
