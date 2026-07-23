# ==============================================================================
# REVISION — Procrustes superimposition plot (trait ordination vs response ordination)
# ==============================================================================
# Visualises the Procrustes/PROTEST result (correlation 0.56, p=0.10). Each species
# appears twice: at its position in TRAIT space (open circle) and its position in
# RESPONSE (methane) space after Procrustes rotation onto trait space (arrow head).
# Short arrow = the species' traits and its methane responses agree; long arrow =
# the species breaks the joint pattern. Arrow length = Procrustes residual.
#
# NEW file. Run: Rscript revision/analysis/R_traits_procrustes_plot.R
# Outputs: revision/outputs/traits_procrustes.png, traits_procrustes_residuals.csv
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(vegan); library(ggrepel) })
set.seed(42)
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

# ---- same matrices as R_traits_dimreduction2.R -------------------------------
source("revision/analysis/_prep_species_data.R")
resp <- analysis_ratio %>% transmute(species_id, flux = median_flux, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id") %>%
  drop_na()
sp10 <- resp$species_id
tr <- read_csv(TRAITS, show_col_types = FALSE) %>% filter(spcode %in% sp10) %>%
  group_by(spcode) %>% slice(1) %>% ungroup()
niche <- read_csv("revision/outputs/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(spcode = names(sp_map)[match(species, sp_map)]) %>% transmute(spcode, vwc_realized = vwc_mean)
use <- c("wood_density_gcm3","bark_density_gcm3","porosity_num","gymnosperm",
         "wood_sapwood_pH","wood_heartwood_pH","try_antifungal","try_rooting_depth",
         "try_plant_height","try_plant_longevity")
X <- tr %>% select(spcode, all_of(use)) %>% left_join(niche, by = "spcode") %>% column_to_rownames("spcode")
X <- X[sp10, ]; for (j in names(X)) X[[j]][!is.finite(X[[j]])] <- median(X[[j]], na.rm = TRUE)
Xs <- scale(X); Y <- resp %>% column_to_rownames("species_id") %>% .[sp10, ] %>% scale()

tpca <- prcomp(Xs); rpca <- prcomp(Y)
pr <- procrustes(tpca$x[, 1:2], rpca$x[, 1:2], symmetric = TRUE)

# ---- assemble plotting frame -------------------------------------------------
Xt <- as.data.frame(pr$X);    names(Xt) <- c("x1","x2")          # trait config (target)
Yt <- as.data.frame(pr$Yrot); names(Yt) <- c("y1","y2")          # response config, rotated
df <- tibble(species = sp10, tx = Xt$x1, ty = Xt$x2, rx = Yt$y1, ry = Yt$y2) %>%
  mutate(resid = sqrt((tx-rx)^2 + (ty-ry)^2)) %>% arrange(desc(resid))
write.csv(df %>% transmute(species, procrustes_residual = round(resid, 3)),
          file.path(out, "traits_procrustes_residuals.csv"), row.names = FALSE)

p <- ggplot(df) +
  geom_segment(aes(tx, ty, xend = rx, yend = ry, color = resid),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.8) +
  geom_point(aes(tx, ty), shape = 21, fill = "white", color = "grey40", size = 2.5) +   # trait position
  geom_point(aes(rx, ry), size = 1.6, color = "grey20") +                                # response position (arrowhead)
  geom_text_repel(aes((tx+rx)/2, (ty+ry)/2, label = species), size = 2.7, max.overlaps = 20) +
  scale_color_viridis_c(option = "C", direction = -1, name = "Procrustes\nresidual\n(misfit)") +
  labs(x = "Procrustes axis 1", y = "Procrustes axis 2",
       title = "Procrustes: trait ordination vs methane-response ordination (n=10)",
       subtitle = "Open circle = species in TRAIT space; arrowhead = same species in RESPONSE space (rotated onto traits).\nShort arrow = traits & methane responses agree; long arrow = species breaks the joint pattern.  PROTEST r=0.56, p=0.10.") +
  coord_equal() + theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 7.5))
ggsave(file.path(out, "traits_procrustes.png"), p, width = 8, height = 7.2, dpi = 300, bg = "white")

cat("Wrote", file.path(out, "traits_procrustes.png"), "\n\n")
cat("Procrustes residuals per species (largest = worst fit to the joint pattern):\n")
print(as.data.frame(df %>% transmute(species, residual = round(resid, 3))), row.names = FALSE)
pt <- protest(tpca$x[, 1:2], rpca$x[, 1:2], permutations = 999)
cat(sprintf("\nPROTEST correlation = %.3f, p = %.3f (matches dimreduction2)\n", pt$t0, pt$signif))
