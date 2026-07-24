# ==============================================================================
# REVISION — Plant traits -> stem CH4 FLUX and -> microbial BALANCE (Editor + R3 #3.1)
# ==============================================================================
# The load-bearing revision concern (editor De Kauwe; R3): the paper is
# "plant-light / microbe-heavy" and never says WHAT ABOUT A SPECIES drives stem
# CH4 flux or the microbial community. A sibling project (../tree-gas-traits) has
# already assembled per-species wood/plant traits (TRY + GWDD wood density + BAAD
# heartwood fraction + curated wood pH + USDA porosity) and shown they predict
# internal stem CH4 CONCENTRATION (porosity R^2=0.41; phylo lambda=0.97).
#
# BUT that analysis used gas *concentration*, not this manuscript's stem *flux*
# or its methanogen:methanotroph *balance*. This script builds the missing bridge:
# it joins those per-species traits (by 4-letter species code) to THIS repo's
# per-tree stem CH4 flux and ddPCR gene abundances, and asks whether plant traits
# predict (a) stem CH4 flux and (b) the methanogen:methanotroph balance -- i.e.
# the "traits -> microbial balance -> flux" narrative the editor asked for.
#
# Species-level analysis (traits are species-level), mirroring the manuscript's
# species framing. Uses the manuscript's OWN canonical species variables via
# _prep_species_data.R (area-weighted genes across heartwood+sapwood; pooled
# 2021+2023 flux; species with >=5 trees & >=5 flux), which reproduces the
# as-reviewed ratio->flux R2 = 0.513 -- so the trait correlations attach to the
# exact quantities the paper reports, not an ad-hoc reconstruction.
# Spearman rho is primary (flux is heavy-tailed); Pearson r reported alongside.
#
# NEW file. Reads ../tree-gas-traits (read-only) + this repo's compiled data.
# Run: Rscript code/revision/R_traits_flux_microbial_bridge.R
# Outputs: outputs/revision/traits_bridge_correlations.csv,
#          traits_bridge_species_table.csv, traits_bridge_summary.txt,
#          traits_bridge_scatter.png
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(scales) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"

# ---- per-species trait panel (one row per species) ---------------------------
trait_cols <- c("porosity_num","gymnosperm","wood_density_gcm3","bark_density_gcm3",
                "wood_sapwood_pH","wood_heartwood_pH","baad_heartwood_fraction_mean")
traits <- read_csv(TRAITS, show_col_types = FALSE) %>%
  select(species_id = spcode, all_of(trait_cols)) %>% distinct()

trait_labels <- c(porosity_num = "Wood porosity (tracheid<diffuse<ring)",
  gymnosperm = "Gymnosperm (vs angiosperm)", wood_density_gcm3 = "Wood density",
  bark_density_gcm3 = "Bark density", wood_sapwood_pH = "Sapwood pH",
  wood_heartwood_pH = "Heartwood pH", baad_heartwood_fraction_mean = "Heartwood fraction")

# ---- canonical species-level responses (manuscript's own construction) -------
source("code/revision/rev_prep_species_data.R")   # analysis_ratio/mcra/meth, flux_by_species
sp <- analysis_ratio %>%
  transmute(species_id, n_trees, flux = median_flux, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id") %>%
  left_join(traits, by = "species_id")

# signed pseudo-log flux for Pearson (flux spans orders of magnitude incl. small neg)
psl <- function(x) asinh(x / (2 * 0.005))    # arcsinh, ~linear below sigma, ~log above
sp$flux_psl <- psl(sp$flux)

# ---- correlate each trait vs each response -----------------------------------
responses <- list(
  flux    = list(col = "flux",    plog = "flux",     label = "Stem CH4 flux (species median)"),
  balance = list(col = "balance", plog = "balance",  label = "Methanogen:methanotroph ratio (log10, area-wtd)"),
  mcra    = list(col = "mcra",    plog = "mcra",     label = "mcrA abundance (log10, area-wtd)"),
  meth    = list(col = "meth",    plog = "meth",     label = "Methanotroph abundance (log10, area-wtd)"))

rows <- list()
for (tn in trait_cols) for (rn in names(responses)) {
  r <- responses[[rn]]
  x <- sp[[tn]]; y <- sp[[r$col]]; yp <- sp[[r$plog]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) next
  sp_r <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  pe_r <- suppressWarnings(cor.test(x[ok], yp[ok & is.finite(yp)], method = "pearson"))
  rows[[length(rows)+1]] <- tibble(
    trait = tn, trait_label = trait_labels[[tn]], response = rn, response_label = r$label,
    n = sum(ok), rho = unname(sp_r$estimate), p_spearman = sp_r$p.value,
    pearson_r = unname(pe_r$estimate), p_pearson = pe_r$p.value)
}
corr <- bind_rows(rows) %>% arrange(response, p_spearman)
write.csv(corr, file.path(out, "traits_bridge_correlations.csv"), row.names = FALSE)
write.csv(sp,   file.path(out, "traits_bridge_species_table.csv"), row.names = FALSE)

# ---- figure: scatter of the strongest trait relationships --------------------
top <- corr %>% group_by(response) %>% slice_min(p_spearman, n = 2) %>% ungroup()
plt <- top %>% rowwise() %>% mutate(data = list(tibble(
  x = sp[[trait]], y = sp[[responses[[response]]$col]], sp = sp$species_id))) %>%
  ungroup() %>% mutate(panel = paste0(response_label, "\n~ ", trait_label,
    "  (rho=", round(rho,2), ", p=", signif(p_spearman,2), ")")) %>%
  select(panel, data) %>% unnest(data) %>% filter(is.finite(x) & is.finite(y))

p <- ggplot(plt, aes(x, y)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey60", linewidth = 0.6, formula = y~x) +
  geom_text(aes(label = sp), size = 2.6) +
  facet_wrap(~panel, scales = "free", ncol = 2) +
  labs(x = "Trait value", y = "Species median response",
       title = "Plant traits vs stem CH4 flux and microbial balance (species level)") +
  theme_bw(base_size = 9) + theme(strip.text = element_text(size = 7))
ggsave(file.path(out, "traits_bridge_scatter.png"), p, width = 9, height = 7, dpi = 300, bg = "white")

# ---- summary -----------------------------------------------------------------
sink(file.path(out, "traits_bridge_summary.txt"))
cat("==================================================================\n")
cat("Plant traits -> stem CH4 FLUX and -> microbial BALANCE (Editor + R3)\n")
cat("==================================================================\n")
cat("Species-level, manuscript's own variables (via _prep_species_data.R).\n")
cat("Traits from ../tree-gas-traits (TRY+GWDD+BAAD+wood pH+USDA porosity).\n")
cat("Balance = median log10((mcrA+1)/(methanotroph_total+1)), area-weighted.\n\n")
# VALIDITY: reproduce the manuscript's ratio->flux EXACTLY. The paper regresses
# RAW median_flux on median_log_ratio (Pearson r=0.717, R2=0.513). (Correlating a
# pseudo-log-transformed flux instead gives ~0.47 -- a transform artifact, not a
# data discrepancy; raw flux is the paper's convention and reproduces 0.513.)
vok <- is.finite(sp$balance) & is.finite(sp$flux)
vp  <- suppressWarnings(cor.test(sp$balance[vok], sp$flux[vok], method = "pearson"))
vs  <- suppressWarnings(cor.test(sp$balance[vok], sp$flux[vok], method = "spearman"))
cat(sprintf("VALIDITY CHECK (reproduces manuscript ratio->flux): balance vs RAW flux, n=%d,\n  Pearson r=%.3f R2=%.3f (paper: r=0.717, R2=0.513) | Spearman rho=%.2f p=%.3f\n\n",
            sum(vok), vp$estimate, vp$estimate^2, vs$estimate, vs$p.value))
cat("Species with data:", nrow(sp), "| trait coverage varies by trait.\n\n")
for (rn in names(responses)) {
  cat("---", responses[[rn]]$label, "---\n")
  sub <- corr %>% filter(response == rn) %>%
    transmute(trait = trait_label, n, rho = round(rho,2), p_spearman = round(p_spearman,3),
              pearson_r = round(pearson_r,2))
  print(as.data.frame(sub), row.names = FALSE); cat("\n")
}
cat("Reading: rho>0 means the trait tracks HIGHER flux / more methanogen-dominated.\n")
cat("CAVEAT: absolute ddPCR loads subject to the pending x10 issue; the BALANCE\n")
cat("(a ratio) and all trait correlations here are UNAFFECTED by that.\n")
sink()
cat(readLines(file.path(out, "traits_bridge_summary.txt")), sep = "\n")
