# ==============================================================================
# REVISION — EXPANDED plant-trait analysis of stem CH4 flux & microbial balance
# (Editor + R3 #3.1; strengthens R_traits_flux_microbial_bridge.R)
# ==============================================================================
# Builds on the bridge: (1) uses the FULL set of well-covered traits (>=8 of the
# 10 flux-species) from ../tree-gas-traits -- now including ROOTS (rooting depth,
# fine-root diameter/tissue density), MOISTURE AFFINITY (waterlogging tolerance,
# wetland indicator), wood C:N, antifungal chemistry, decomposition rate -- PLUS
# (2) our OWN realized soil-moisture niche per species (vwc_mean from
# tree_species_moisture_niche.csv), and (3) reports predictor collinearity and a
# cautious best 1-2 predictor model with LOO-CV. Responses are the manuscript's
# canonical species vars (via _prep_species_data.R; ratio->flux R2=0.513).
#
# n=10 species -> UNDERPOWERED for true multivariable models. We therefore lead
# with the univariate atlas + collinearity, and treat any multi-predictor model
# as explicitly exploratory (LOO-CV to check it isn't just overfitting).
#
# NEW file. Run: Rscript revision/analysis/R_traits_multitrait.R
# Outputs: traits_multitrait_atlas.csv, traits_multitrait_collinearity.csv,
#          traits_multitrait_models.txt, traits_multitrait_atlas.png, _summary.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

# ---- canonical species responses ---------------------------------------------
source("revision/analysis/_prep_species_data.R")
resp <- analysis_ratio %>%
  transmute(species_id, flux = median_flux, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id")
sp10 <- resp$species_id

# ---- full well-covered trait panel + our realized moisture niche -------------
# (traits with >=8 of 10 flux-species non-NA; drop the gas cols themselves)
keep_traits <- c(
  "porosity_num","gymnosperm","wood_density_gcm3","bark_density_gcm3","bark_wood_ratio",
  "try_bark_volume_rel_wood","wood_sapwood_pH","wood_heartwood_pH","try_wood_CN_ratio",
  "try_stem_C","try_antifungal","try_CWD_stem_decomp_rate_k","try_rooting_depth",
  "try_fine_root_diameter","try_fine_root_tissue_density","try_waterlogging_tolerance",
  "try_wetland_indicator","try_plant_longevity","try_plant_height")
trait_lab <- c(porosity_num="Wood porosity", gymnosperm="Gymnosperm",
  wood_density_gcm3="Wood density", bark_density_gcm3="Bark density",
  bark_wood_ratio="Bark:wood ratio", try_bark_volume_rel_wood="Bark volume rel. wood",
  wood_sapwood_pH="Sapwood pH", wood_heartwood_pH="Heartwood pH",
  try_wood_CN_ratio="Wood C:N", try_stem_C="Stem C", try_antifungal="Antifungal chem.",
  try_CWD_stem_decomp_rate_k="CWD decomp rate", try_rooting_depth="Rooting depth",
  try_fine_root_diameter="Fine-root diameter", try_fine_root_tissue_density="Fine-root tissue density",
  try_waterlogging_tolerance="Waterlogging tol. (trait)", try_wetland_indicator="Wetland indicator (trait)",
  try_plant_longevity="Plant longevity", try_plant_height="Plant height",
  vwc_realized="Realized soil moisture (OUR VWC)")

traits <- read_csv(TRAITS, show_col_types = FALSE) %>%
  filter(spcode %in% sp10) %>% group_by(spcode) %>% slice(1) %>% ungroup() %>%
  select(species_id = spcode, any_of(keep_traits))

# our realized moisture niche (per species, from our data)
niche <- read_csv("revision/outputs/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(species_id = names(sp_map)[match(species, sp_map)]) %>%
  transmute(species_id, vwc_realized = vwc_mean)

dat <- resp %>% left_join(traits, by = "species_id") %>% left_join(niche, by = "species_id")
pred_cols <- c(setdiff(names(traits), "species_id"), "vwc_realized")

# ---- univariate atlas: every trait x every response --------------------------
responses <- c(flux="Stem CH4 flux", balance="Methanogen:methanotroph balance",
               mcra="mcrA abundance", meth="Methanotroph abundance")
rows <- list()
for (tn in pred_cols) for (rn in names(responses)) {
  x <- dat[[tn]]; y <- dat[[rn]]; ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 6 || sd(x[ok]) == 0) next
  s <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  p <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  rows[[length(rows)+1]] <- tibble(trait = tn, trait_label = trait_lab[[tn]],
    response = rn, response_label = responses[[rn]], n = sum(ok),
    rho = unname(s$estimate), p_spearman = s$p.value,
    pearson_r = unname(p$estimate), p_pearson = p$p.value)
}
atlas <- bind_rows(rows) %>% arrange(response, p_spearman)
write.csv(atlas, file.path(out, "traits_multitrait_atlas.csv"), row.names = FALSE)

# ---- predictor collinearity (are the "hits" one axis?) -----------------------
pm <- dat[, pred_cols]
cc <- cor(pm, use = "pairwise.complete.obs", method = "spearman")
write.csv(round(cc, 2), file.path(out, "traits_multitrait_collinearity.csv"))

# ---- cautious model: best 1- and 2-predictor per response, with LOO-CV -------
loo_r2 <- function(df, form) {           # leave-one-out CV R2
  n <- nrow(df); pr <- numeric(n)
  for (i in seq_len(n)) {
    m <- try(lm(form, df[-i, ]), silent = TRUE)
    if (inherits(m, "try-error")) { pr[i] <- NA; next }
    pr[i] <- predict(m, df[i, ])
  }
  y <- df[[all.vars(form)[1]]]
  1 - sum((y - pr)^2, na.rm = TRUE) / sum((y - mean(y))^2)
}
sink(file.path(out, "traits_multitrait_models.txt"))
cat("Best-subset models (n=10 species; EXPLORATORY -- LOO-CV guards overfit).\n")
cat("in-sample adjR2 will flatter; trust LOO-CV R2.\n\n")
for (rn in names(responses)) {
  cat("================", responses[[rn]], "================\n")
  d2 <- dat[is.finite(dat[[rn]]), ]
  cand <- pred_cols[sapply(pred_cols, function(p) sum(is.finite(d2[[p]])) >= 9)]
  # 1-predictor
  best1 <- map_dfr(cand, function(p) {
    dd <- d2[is.finite(d2[[p]]), c(rn, p)]; names(dd) <- c("y","x")
    m <- lm(y ~ x, dd)
    tibble(p, adjR2 = summary(m)$adj.r.squared, loo = loo_r2(dd, y ~ x))
  }) %>% arrange(desc(loo))
  cat("-- top single predictors (by LOO-CV R2):\n")
  print(as.data.frame(head(best1, 5) %>% mutate(across(c(adjR2,loo), ~round(.,3)))), row.names = FALSE)
  # 2-predictor (top-4 singles, all pairs), require both finite
  top4 <- head(best1$p, 4)
  pairs <- combn(top4, 2, simplify = FALSE)
  best2 <- map_dfr(pairs, function(pp) {
    dd <- d2[is.finite(d2[[pp[1]]]) & is.finite(d2[[pp[2]]]), c(rn, pp)]
    names(dd) <- c("y","x1","x2"); m <- lm(y ~ x1 + x2, dd)
    tibble(pair = paste(pp, collapse=" + "), adjR2 = summary(m)$adj.r.squared, loo = loo_r2(dd, y ~ x1 + x2))
  }) %>% arrange(desc(loo))
  cat("-- top 2-predictor models (by LOO-CV R2):\n")
  print(as.data.frame(head(best2, 3) %>% mutate(across(c(adjR2,loo), ~round(.,3)))), row.names = FALSE)
  cat("\n")
}
sink()

# ---- atlas figure (forest plot, faceted by response) -------------------------
atl <- atlas %>% mutate(sig = p_spearman < 0.05,
  trait_label = fct_reorder(trait_label, rho))
p <- ggplot(atl, aes(rho, trait_label, color = sig)) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_point(aes(size = n)) +
  scale_color_manual(values = c(`FALSE`="grey65", `TRUE`="#b2182b"), guide = "none") +
  scale_size(range = c(1, 3), guide = "none") +
  facet_wrap(~response_label, ncol = 4) +
  labs(x = "Spearman rho (species level)", y = NULL,
       title = "Plant traits vs stem CH4 flux & microbial balance (red = p<0.05)",
       subtitle = "n<=10 species; incl. roots, trait & realized moisture, wood chemistry") +
  theme_bw(base_size = 8) + theme(axis.text.y = element_text(size = 7))
ggsave(file.path(out, "traits_multitrait_atlas.png"), p, width = 12, height = 5.5, dpi = 300, bg = "white")

# ---- summary -----------------------------------------------------------------
sink(file.path(out, "traits_multitrait_summary.txt"))
cat("=================================================================\n")
cat("EXPANDED trait analysis -- flux & microbial balance (Editor + R3)\n")
cat("=================================================================\n")
cat("Predictors:", length(pred_cols), "traits (incl. roots, trait+realized moisture,\n")
cat("wood chemistry). Responses = manuscript canonical species vars. n<=10 species.\n\n")
cat("SIGNIFICANT (p_spearman<0.05):\n")
sig <- atlas %>% filter(p_spearman < 0.05) %>%
  transmute(response = response_label, trait = trait_label, n, rho = round(rho,2), p = round(p_spearman,3))
print(as.data.frame(sig), row.names = FALSE)
cat("\nNULLS worth reporting (answer R3 directly):\n")
for (t in c("try_rooting_depth","try_wetland_indicator","try_waterlogging_tolerance","vwc_realized","porosity_num")) {
  sub <- atlas %>% filter(trait == t, response == "flux")
  if (nrow(sub)) cat(sprintf("  %-28s vs flux: rho=%.2f p=%.2f\n", trait_lab[[t]], sub$rho, sub$p_spearman))
}
cat("\nSee traits_multitrait_models.txt for LOO-CV best-subset (exploratory, n=10).\n")
cat("See traits_multitrait_collinearity.csv for predictor inter-correlations.\n")
cat("CAVEAT: absolute ddPCR loads subject to pending x10; ratios/correlations unaffected.\n")
sink()
cat(readLines(file.path(out, "traits_multitrait_summary.txt")), sep = "\n")
