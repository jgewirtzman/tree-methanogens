# ==============================================================================
# REVISION — Dimension reduction: trait axes vs the stem methane cycle (Editor+R3)
# ==============================================================================
# With ~14 collinear traits and n=10 species, univariate scanning is the wrong
# tool. Here we (1) collapse traits into a few interpretable AXES via PCA and ask
# which axis aligns with the methane-cycling responses, and (2) run SUPERVISED
# PLS (traits -> response) with leave-one-out CV -- the standard approach when
# predictors outnumber samples -- to see whether ANY trait combination predicts
# flux / the microbial balance out-of-sample.
#
# Responses = manuscript canonical species vars (via _prep_species_data.R;
# reproduces ratio->flux R2=0.513). Honest about n=10: PCA is descriptive; PLS is
# judged ONLY by cross-validated R2 (in-sample R2 at n=10 is meaningless).
#
# NEW file. Run: Rscript code/revision/R_traits_dimreduction.R
# Outputs: traits_pca_biplot.png, traits_pca_loadings.csv, traits_pc_response_cor.csv,
#          traits_pls_cv.txt, traits_dimreduction_summary.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(pls); library(ggrepel) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

# ---- responses ---------------------------------------------------------------
source("code/lib/rev_prep_species_data.R")
resp <- analysis_ratio %>% transmute(species_id, flux = median_flux, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id")
sp10 <- resp$species_id
resp_lab <- c(flux="Stem CH4 flux", balance="mcrA:MOB balance", mcra="mcrA (methanogen)", meth="Methanotroph")

# ---- trait matrix (curated, median-imputed to complete) + our VWC ------------
tr <- read_csv(TRAITS, show_col_types = FALSE) %>% filter(spcode %in% sp10) %>%
  group_by(spcode) %>% slice(1) %>% ungroup()
niche <- read_csv("outputs/revision/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(species_id = names(sp_map)[match(species, sp_map)]) %>% transmute(species_id, vwc_realized = vwc_mean)
use <- c("wood_density_gcm3","bark_density_gcm3","bark_wood_ratio","porosity_num","gymnosperm",
         "wood_sapwood_pH","wood_heartwood_pH","try_wood_CN_ratio","try_antifungal",
         "try_CWD_stem_decomp_rate_k","try_rooting_depth","try_fine_root_diameter",
         "try_plant_height","try_plant_longevity")
short <- c(wood_density_gcm3="WoodDens", bark_density_gcm3="BarkDens", bark_wood_ratio="Bark:Wood",
  porosity_num="Porosity", gymnosperm="Gymnosperm", wood_sapwood_pH="SapwoodpH", wood_heartwood_pH="HeartwoodpH",
  try_wood_CN_ratio="WoodC:N", try_antifungal="Antifungal", try_CWD_stem_decomp_rate_k="DecompRate",
  try_rooting_depth="RootDepth", try_fine_root_diameter="FineRootDiam", try_plant_height="Height",
  try_plant_longevity="Longevity", vwc_realized="RealizedVWC")

X <- tr %>% select(species_id = spcode, all_of(use)) %>% left_join(niche, by = "species_id")
Xm <- X %>% column_to_rownames("species_id") %>% as.data.frame()
Xm <- Xm[sp10, ]                                   # align order to responses
# median-impute the few missing cells, then scale
for (j in names(Xm)) Xm[[j]][!is.finite(Xm[[j]])] <- median(Xm[[j]], na.rm = TRUE)
colnames(Xm) <- short[colnames(Xm)]
Xs <- scale(Xm)

# ============================ 1. TRAIT PCA ====================================
pca <- prcomp(Xs, center = FALSE, scale. = FALSE)
ve <- (pca$sdev^2) / sum(pca$sdev^2)
scores <- as.data.frame(pca$x[, 1:3]); scores$species_id <- rownames(scores)
load <- as.data.frame(pca$rotation[, 1:3]) %>% rownames_to_column("trait")
write.csv(cbind(load, var_explained_note = ""), file.path(out, "traits_pca_loadings.csv"), row.names = FALSE)

# correlate PC axes with responses (few tests)
pcr_rows <- list()
sc <- left_join(scores, resp, by = "species_id")
for (pc in c("PC1","PC2","PC3")) for (rn in names(resp_lab)) {
  ct <- suppressWarnings(cor.test(sc[[pc]], sc[[rn]], method = "spearman"))
  pcr_rows[[length(pcr_rows)+1]] <- tibble(PC = pc, response = resp_lab[[rn]],
    rho = unname(ct$estimate), p = ct$p.value)
}
pc_resp <- bind_rows(pcr_rows)
write.csv(pc_resp, file.path(out, "traits_pc_response_cor.csv"), row.names = FALSE)

# ---- biplot: species (colored by flux) + trait loadings + response arrows ----
mult <- 2.6
tl <- load %>% transmute(trait, x = PC1 * mult, y = PC2 * mult)
# project each response as a supplementary arrow (correlation of response with PC1/PC2)
rp <- map_dfr(names(resp_lab), function(rn) {
  tibble(response = resp_lab[[rn]],
         x = cor(sc$PC1, sc[[rn]], use = "complete.obs") * mult,
         y = cor(sc$PC2, sc[[rn]], use = "complete.obs") * mult) })
sc <- left_join(sc, resp, by = "species_id", suffix = c("", ".y"))
bip <- ggplot() +
  geom_hline(yintercept = 0, color = "grey85") + geom_vline(xintercept = 0, color = "grey85") +
  geom_segment(data = tl, aes(0, 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.15, "cm")), color = "grey55") +
  geom_text_repel(data = tl, aes(x, y, label = trait), color = "grey40", size = 2.7, max.overlaps = 20) +
  geom_segment(data = rp, aes(0, 0, xend = x, yend = y), color = "#b2182b",
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.9) +
  geom_text_repel(data = rp, aes(x, y, label = response), color = "#b2182b",
                  fontface = "bold", size = 3, max.overlaps = 20) +
  geom_point(data = sc, aes(PC1, PC2, fill = flux), shape = 21, size = 4) +
  geom_text_repel(data = sc, aes(PC1, PC2, label = species_id), size = 2.6, max.overlaps = 20) +
  scale_fill_viridis_c(option = "C", name = "Stem CH4\nflux") +
  labs(x = sprintf("PC1 (%.0f%%)", 100*ve[1]), y = sprintf("PC2 (%.0f%%)", 100*ve[2]),
       title = "Trait-space PCA with methane responses projected (red)",
       subtitle = "Grey = plant traits; red = methane-cycling responses (as supplementary correlation arrows); points = species (n=10)") +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 7.5))
ggsave(file.path(out, "traits_pca_biplot.png"), bip, width = 8.5, height = 7, dpi = 300, bg = "white")

# ============================ 2. SUPERVISED PLS ===============================
# LOO-CV PLS for each response; report cross-validated R2 (in-sample is meaningless at n=10)
sink(file.path(out, "traits_pls_cv.txt"))
cat("PLS regression traits -> response, leave-one-out CV (n=10).\n")
cat("Judge ONLY by CV R2 (>0 means predicts out-of-sample). Up to 3 components.\n\n")
pls_summary <- list()
for (rn in names(resp_lab)) {
  df <- data.frame(y = resp[[rn]], Xs); df <- df[is.finite(df$y), ]
  m <- plsr(y ~ ., data = df, ncomp = 3, validation = "LOO", scale = FALSE)
  r2cv <- drop(pls::R2(m, estimate = "CV")$val)[-1]          # drop intercept-only
  best <- which.max(r2cv);
  cat(sprintf("%-18s CV R2 by ncomp: %s  -> best = %d comp, CV R2 = %.3f\n",
              resp_lab[[rn]], paste(sprintf("%.2f", r2cv), collapse = " / "), best, r2cv[best]))
  # trait loadings on component 1 (which trait combo defines the predictive axis)
  ld <- sort(loadings(m)[,1], decreasing = TRUE)
  cat("   comp-1 trait loadings (top +):", paste(sprintf("%s=%.2f", names(head(ld,4)), head(ld,4)), collapse=", "), "\n")
  cat("   comp-1 trait loadings (top -):", paste(sprintf("%s=%.2f", names(tail(ld,4)), tail(ld,4)), collapse=", "), "\n\n")
  pls_summary[[rn]] <- tibble(response = resp_lab[[rn]], best_ncomp = best, cv_r2 = r2cv[best])
}
sink()

# ============================ summary =========================================
sink(file.path(out, "traits_dimreduction_summary.txt"))
cat("=================================================================\n")
cat("DIMENSION REDUCTION -- trait axes vs stem methane cycle (Editor+R3)\n")
cat("=================================================================\n")
cat("n=10 species; 15 traits (incl. roots + realized VWC), scaled.\n\n")
cat("PCA variance explained: ", paste(sprintf("PC%d=%.0f%%", 1:5, 100*ve[1:5]), collapse="  "), "\n\n")
cat("PC1 top loadings:\n"); print(round(sort(pca$rotation[,1]),2))
cat("\nPC2 top loadings:\n"); print(round(sort(pca$rotation[,2]),2))
cat("\nPC axis vs response (Spearman):\n")
print(as.data.frame(pc_resp %>% mutate(rho=round(rho,2), p=round(p,3))), row.names = FALSE)
cat("\nSupervised PLS cross-validated R2 (the honest predictive test):\n")
print(as.data.frame(bind_rows(pls_summary) %>% mutate(cv_r2 = round(cv_r2,3))), row.names = FALSE)
cat("\nReading: CV R2 > 0 = a trait axis predicts that response out-of-sample.\n")
cat("CAVEAT: n=10; PCA descriptive, PLS judged by CV R2 only. ddPCR x10 pending;\n")
cat("balance is a ratio -> unaffected.\n")
sink()
cat(readLines(file.path(out,"traits_dimreduction_summary.txt")), sep="\n")
cat("\n\n---- PLS CV ----\n"); cat(readLines(file.path(out,"traits_pls_cv.txt")), sep="\n")
