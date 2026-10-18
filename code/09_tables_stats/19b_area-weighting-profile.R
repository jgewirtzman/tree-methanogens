# ==============================================================================
# 19b_area-weighting-profile.R -- does the radial profile of the area weighting matter?
# (Methods S4)
# ------------------------------------------------------------------------------
# Tree-level genes are averaged over the breast-height cross-section with the sapwood
# value in the outer 5 cm, the heartwood value within min(5, R - 5) cm of the pith and
# a linear transition between (area_weighted_gene(), code/lib/species_gene_data.R, the
# one copy of the formula). Here the transition is replaced by:
#   s_moderate / s_steep - logistic in the transition fraction, k = 6 / k = 20
#   step_max             - heartwood value everywhere inside the outer 5 cm
#   step_min             - heartwood value only within min(5, R - 5) cm of the pith
# and the gene-flux results recomputed exactly as reported:
#   species level (ten species, >= 5 trees and >= 5 fluxes): lm(median_flux ~
#     median log ratio), lm(median_flux ~ log10(median value + 1)) for mcrA and
#     pmoA + mmoX (R2, p);
#   tree level: R2 of log10 mcrA alone; added R2 and p of each gene over species
#     (lm(CH4_flux ~ species + gene));
#   the three gene-based trait associations named in Results (25_plant-traits.R):
#     Spearman rho and p across all ten and the eight broadleaf species, and whether
#     the Methods S5 reporting rule still holds.
# Output: outputs/data/area_weighting_profile_sensitivity.csv (printed by
#   20_manuscript_statistics.R)
# ==============================================================================
source("code/lib/outputs.R")
invisible(capture.output(suppressMessages(source("code/lib/species_gene_data.R"))))

PROFILES <- c("linear", "s_moderate", "s_steep", "step_max", "step_min")

TREE_GENES <- c(mcrA = "log_tree_mcra", mmoX = "log_tree_mmox",
                methanotroph = "log_methanotroph", ratio = "log_ratio")
lm_fit <- function(y, x) { s <- summary(lm(y ~ x)); c(R2 = s$r.squared, p = coef(s)[2, 4]) }
added <- function(d, g) {
  m0 <- lm(CH4_flux ~ species, d); m1 <- lm(as.formula(paste("CH4_flux ~ species +", g)), d)
  c(dR2 = summary(m1)$r.squared - summary(m0)$r.squared, p = coef(summary(m1))[g, 4])
}

# Trait values for the named associations (columns as in 25_plant-traits.R)
tr <- read_csv("data/tables/species_traits.csv", show_col_types = FALSE)
niche <- read_csv("outputs/data/tree_species_moisture_niche.csv", show_col_types = FALSE)
gym <- read_table("species") %>% transmute(species_id = species_code, gym = as.numeric(gymnosperm))
NAMED <- tribble(~trait, ~col, ~response,
  "Bark density", "bark_density_gcm3", "balance",
  "Soil-moisture niche", "vwc_realized", "mcra",
  "Longevity", "try_plant_longevity", "meth")
trait_test <- function(tb, trait_col, response) {
  d <- tb$analysis_ratio %>% transmute(species_id, species, balance = median_log_ratio) %>%
    left_join(tb$analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
    left_join(tb$analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id") %>%
    left_join(gym, by = "species_id") %>%
    mutate(x = if (trait_col == "vwc_realized") niche$vwc_mean[match(species, niche$species)]
               else as.numeric(tr[[trait_col]][match(species_id, tr$species_code)]),
           y = .data[[response]]) %>% filter(is.finite(x), is.finite(y))
  sp <- function(k) { ct <- suppressWarnings(cor.test(d$x[k], d$y[k], method = "spearman"))
                      c(unname(ct$estimate), ct$p.value) }
  a <- sp(rep(TRUE, nrow(d))); b <- sp(d$gym == 0)
  con <- with(d, sign((x[species_id == "TSCA"] - x[species_id == "PIST"]) *
                      (y[species_id == "TSCA"] - y[species_id == "PIST"])))
  rule <- (a[2] < 0.05 & b[2] < 0.05 & sign(a[1]) == sign(b[1])) |
          (b[2] < 0.05 & sign(a[1]) == sign(b[1]) & con == sign(b[1]))
  c(rho = a[1], p = a[2], rho_bl = b[1], p_bl = b[2], rule = as.numeric(rule))
}

res <- lapply(PROFILES, function(pr) {
  tb <- build_gene_tables(pr); tl <- tb$tree_level_complete
  r <- lm_fit(tb$analysis_ratio$median_flux, tb$analysis_ratio$median_log_ratio)
  m <- lm_fit(tb$analysis_mcra$median_flux, log10(tb$analysis_mcra$value + 1))
  t <- lm_fit(tb$analysis_meth$median_flux, log10(tb$analysis_meth$value + 1))
  row <- tibble(profile = pr, n_species = nrow(tb$analysis_ratio), n_trees = nrow(tl),
    sp_ratio_R2 = r[["R2"]], sp_ratio_p = r[["p"]], sp_mcrA_R2 = m[["R2"]], sp_mcrA_p = m[["p"]],
    sp_methanotroph_R2 = t[["R2"]], sp_methanotroph_p = t[["p"]],
    tree_mcrA_R2 = summary(lm(CH4_flux ~ log_tree_mcra, tl))$r.squared)
  for (g in names(TREE_GENES)) {
    a <- added(tl, TREE_GENES[[g]])
    row[[paste0("tree_", g, "_addedR2")]] <- a[["dR2"]]; row[[paste0("tree_", g, "_p")]] <- a[["p"]]
  }
  for (i in seq_len(nrow(NAMED))) {
    v <- trait_test(tb, NAMED$col[i], NAMED$response[i])
    key <- paste0("trait_", gsub("[^a-z]+", "_", tolower(NAMED$trait[i])), "_", NAMED$response[i])
    row[[paste0(key, "_rho")]] <- v[["rho"]]; row[[paste0(key, "_p")]] <- v[["p"]]
    row[[paste0(key, "_rho_bl")]] <- v[["rho_bl"]]; row[[paste0(key, "_p_bl")]] <- v[["p_bl"]]
    row[[paste0(key, "_reported")]] <- as.logical(v[["rule"]])
  }
  row
}) %>% bind_rows()

write.csv(res, out_path("area_weighting_profile_sensitivity.csv"), row.names = FALSE)
print(as.data.frame(res %>% select(profile:tree_mcrA_R2) %>% mutate(across(where(is.numeric), ~ signif(., 3)))),
      row.names = FALSE)
