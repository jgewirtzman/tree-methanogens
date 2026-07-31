# ==============================================================================
# REVISION — Trait x (gas / microbe / flux) correlation HEATMAP (Editor + R3 #3.1)
# ==============================================================================
# A DESCRIPTIVE reframe artifact, not a significance test. Shows how plant traits
# structure the manuscript's own methane-cycling quantities, columns ordered along
# the microbial chain:  methanogens (mcrA) -> methanotrophs -> balance -> stem flux.
# Cell = Spearman rho (species level). Faint * marks raw p<0.05 (context only;
# we do NOT claim FDR-corrected significance at n=10 -- see the finding note).
#
# NOTE: internal O2/CH4 columns were DROPPED -- in the YMF gas set 61% of stem CH4
# samples are near-atmospheric (<3 ppm) and O2 medians are near-ambient, so species
# medians are unreliable and the simple O2-diffusion-barrier mechanism is NOT
# supported (wood density -> O2 deficit rho=-0.43 wrong sign; O2 deficit ->
# methanotroph rho=+0.12 n.s.). Density->methanotroph suppression is a robust
# association with an OPEN mechanism, not an O2 story. See finding note.
#
# NEW file. Run: Rscript code/revision/R_traits_heatmap.R
# Outputs: outputs/revision/traits_heatmap.png, traits_heatmap_matrix.csv
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

# ---- responses (manuscript canonical microbial + flux quantities) ------------
source("code/revision/rev_prep_species_data.R")
resp <- analysis_ratio %>% transmute(species_id, `Stem CH4 flux` = median_flux,
                                      `Balance (mcrA:MOB)` = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, `mcrA (methanogen)` = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, `Methanotroph`      = log10(value + 1)), by = "species_id")
sp10 <- resp$species_id
tr_raw <- read_csv(TRAITS, show_col_types = FALSE) %>% filter(spcode %in% sp10)

# ---- trait panel + our realized moisture niche -------------------------------
keep <- c("wood_density_gcm3","bark_density_gcm3","bark_wood_ratio","porosity_num","gymnosperm",
          "wood_sapwood_pH","wood_heartwood_pH","try_wood_CN_ratio","try_stem_C","try_antifungal",
          "try_CWD_stem_decomp_rate_k","try_rooting_depth","try_fine_root_diameter",
          "try_fine_root_tissue_density","try_plant_height","try_plant_longevity")
trait_lab <- c(wood_density_gcm3="Wood density", bark_density_gcm3="Bark density",
  bark_wood_ratio="Bark:wood ratio", porosity_num="Wood porosity", gymnosperm="Gymnosperm",
  wood_sapwood_pH="Sapwood pH", wood_heartwood_pH="Heartwood pH", try_wood_CN_ratio="Wood C:N",
  try_stem_C="Stem C", try_antifungal="Antifungal chem.", try_CWD_stem_decomp_rate_k="CWD decomp rate",
  try_rooting_depth="Rooting depth", try_fine_root_diameter="Fine-root diameter",
  try_fine_root_tissue_density="Fine-root tissue density", try_plant_height="Plant height",
  try_plant_longevity="Plant longevity", vwc_realized="Realized soil moisture (VWC)")
# trait category (row grouping)
cat_of <- c(`Wood density`="Structure",`Bark density`="Structure",`Bark:wood ratio`="Structure",
  `Wood porosity`="Structure",`Gymnosperm`="Structure",`Sapwood pH`="Chemistry",`Heartwood pH`="Chemistry",
  `Wood C:N`="Chemistry",`Stem C`="Chemistry",`Antifungal chem.`="Chemistry",`CWD decomp rate`="Chemistry",
  `Rooting depth`="Roots",`Fine-root diameter`="Roots",`Fine-root tissue density`="Roots",
  `Plant height`="Whole-plant",`Plant longevity`="Whole-plant",`Realized soil moisture (VWC)`="Moisture")

traits <- tr_raw %>% group_by(species_id = spcode) %>% slice(1) %>% ungroup() %>%
  select(species_id, any_of(keep))
niche <- read_csv("outputs/revision/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(species_id = names(sp_map)[match(species, sp_map)]) %>% transmute(species_id, vwc_realized = vwc_mean)
dat <- resp %>% left_join(traits, by = "species_id") %>% left_join(niche, by = "species_id")

pred_cols <- c(keep, "vwc_realized")
resp_cols <- c("mcrA (methanogen)","Methanotroph","Balance (mcrA:MOB)","Stem CH4 flux")

# ---- correlation matrix ------------------------------------------------------
grid <- expand_grid(trait = pred_cols, response = resp_cols) %>%
  rowwise() %>%
  mutate(x = list(dat[[trait]]), y = list(dat[[response]])) %>%
  mutate(ok = list(is.finite(unlist(x)) & is.finite(unlist(y)))) %>%
  mutate(n = sum(unlist(ok)),
         rho = ifelse(n >= 6, suppressWarnings(cor(unlist(x)[unlist(ok)], unlist(y)[unlist(ok)], method = "spearman")), NA_real_),
         p   = ifelse(n >= 6, suppressWarnings(cor.test(unlist(x)[unlist(ok)], unlist(y)[unlist(ok)], method = "spearman")$p.value), NA_real_)) %>%
  ungroup() %>%
  mutate(trait_label = trait_lab[trait], category = cat_of[trait_label],
         star = ifelse(!is.na(p) & p < 0.05, "*", ""))

write.csv(grid %>% select(category, trait_label, response, n, rho, p),
          file.path(out, "traits_heatmap_matrix.csv"), row.names = FALSE)

# order rows by category then by correlation with flux
row_order <- grid %>% filter(response == "Stem CH4 flux") %>%
  arrange(factor(category, levels=c("Structure","Chemistry","Roots","Moisture","Whole-plant")), rho) %>%
  pull(trait_label)
grid <- grid %>% mutate(trait_label = factor(trait_label, levels = row_order),
                        response = factor(response, levels = resp_cols),
                        category = factor(category, levels=c("Structure","Chemistry","Roots","Moisture","Whole-plant")))

p <- ggplot(grid, aes(response, trait_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f%s", rho, star)), size = 2.6) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Plant traits structure the stem methane-cycling community (n<=10 species)",
       caption = paste("Cell = Spearman rho (species level). Columns follow the microbial chain:",
                       "methanogens -> methanotrophs -> balance -> net flux.\n",
                       "* raw p<0.05, EXPLORATORY: nothing survives FDR at n<=10 -- read as descriptive structure, not hypothesis tests.")) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x.top = element_text(angle = 15, hjust = 0, size = 8.5),
        axis.text.y = element_text(size = 8),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
        strip.placement = "outside", panel.grid = element_blank(),
        plot.title = element_text(size = 12, margin = margin(b = 22)),
        plot.caption = element_text(size = 7.5, hjust = 0, margin = margin(t = 10)),
        plot.margin = margin(t = 10, r = 45, b = 6, l = 6))
ggsave(file.path(out, "traits_heatmap.png"), p, width = 9, height = 7.4, dpi = 300, bg = "white")
cat("Wrote", file.path(out, "traits_heatmap.png"), "\n")
cat("Rows =", length(unique(grid$trait_label)), "traits (5 categories); columns = 6 responses along the causal chain.\n")
