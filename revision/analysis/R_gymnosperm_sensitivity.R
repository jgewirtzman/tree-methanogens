# ==============================================================================
# REVISION — Are the trait effects real WITHIN angiosperms, or just the
# gymnosperm/angiosperm split? (gymnosperms = Procrustes misfits)
# ==============================================================================
# The 2 gymnosperms (PIST, TSCA) are the biggest Procrustes misfits and anchor the
# low-density end of trait PC1. Question: do the trait->response associations
# (density->methanotroph/balance; moisture->methanogen) reflect a CONTINUOUS trait
# gradient that also holds among the 8-13 angiosperms, or are they driven by the 2
# gymnosperms (i.e. really a conifer-vs-hardwood categorical difference)?
#
# Approach: recompute each key species-level relationship (a) ALL species and
# (b) ANGIOSPERMS ONLY, and describe where the 2 gymnosperms sit on each response.
# n=2 gymnosperms cannot be modelled as their own group -- only flagged/described.
#
# NEW file. Run: Rscript revision/analysis/R_gymnosperm_sensitivity.R
# Output: revision/outputs/gymnosperm_sensitivity.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

# ---- per-species responses (use WIDE set: microbe from all wood-sampled spp) --
source("revision/analysis/_prep_species_data.R")
micro <- tree_level_complete %>% group_by(species_id) %>%
  summarise(meth = median(log_meth), mcra = median(log_mcra), balance = median(log_ratio),
            n_tree = n(), .groups = "drop")
flux <- flux_by_species %>% transmute(species_id, flux = median_flux)
tr <- read_csv(TRAITS, show_col_types = FALSE) %>% group_by(spcode) %>% slice(1) %>% ungroup() %>%
  transmute(species_id = spcode, wood_density = wood_density_gcm3, bark_density = bark_density_gcm3,
            porosity_type, gymnosperm)
niche <- read_csv("revision/outputs/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(species_id = names(sp_map)[match(species, sp_map)]) %>% transmute(species_id, vwc = vwc_mean)

sp <- micro %>% full_join(flux, by = "species_id") %>%
  left_join(tr, by = "species_id") %>% left_join(niche, by = "species_id") %>%
  mutate(clade = if_else(gymnosperm == 1, "Gymnosperm", "Angiosperm"))

cmp <- function(yvar, xvar) {
  f <- function(d) { d <- d %>% filter(is.finite(.data[[yvar]]), is.finite(.data[[xvar]]))
    if (nrow(d) < 5) return(c(n = nrow(d), rho = NA, p = NA))
    s <- suppressWarnings(cor.test(d[[xvar]], d[[yvar]], method = "spearman"))
    c(n = nrow(d), rho = unname(s$estimate), p = s$p.value) }
  all <- f(sp); ang <- f(sp %>% filter(clade == "Angiosperm"))
  sprintf("  %-9s ~ %-13s | ALL: n=%d rho=%+.2f p=%.3f  |  ANGIOSPERM-ONLY: n=%d rho=%+.2f p=%.3f",
          yvar, xvar, all["n"], all["rho"], all["p"], ang["n"], ang["rho"], ang["p"])
}

sink(file.path(out, "gymnosperm_sensitivity.txt"))
cat("=================================================================\n")
cat("Trait effects: ALL species vs ANGIOSPERM-ONLY (drop 2 gymnosperms)\n")
cat("=================================================================\n")
cat("If a relationship HOLDS angiosperm-only -> continuous trait effect (real).\n")
cat("If it VANISHES -> it was really the conifer/hardwood categorical split.\n\n")
cat(cmp("meth",    "wood_density"), "\n")
cat(cmp("meth",    "bark_density"), "\n")
cat(cmp("balance", "bark_density"), "\n")
cat(cmp("balance", "wood_density"), "\n")
cat(cmp("mcra",    "vwc"),          "\n")
cat(cmp("flux",    "bark_density"), "\n")
cat(cmp("flux",    "wood_density"), "\n\n")

# ---- where do the gymnosperms sit on each response? (z vs angiosperm mean) ----
cat("Where the 2 gymnosperms sit on each response (z-score vs angiosperm mean):\n")
zof <- function(v) { a <- sp %>% filter(clade == "Angiosperm") %>% pull(v)
  g <- sp %>% filter(clade == "Gymnosperm"); (g[[v]] - mean(a, na.rm=T)) / sd(a, na.rm=T) }
gz <- tibble(species = sp$species_id[sp$clade=="Gymnosperm"])
for (v in c("meth","mcra","balance","flux","wood_density","bark_density","vwc"))
  gz[[v]] <- round(zof(v), 2)
print(as.data.frame(gz), row.names = FALSE)
cat("\n(porosity_type: gymnosperms are 'tracheid'; angiosperms are diffuse-/ring-porous.)\n")

# ---- does clade ALONE separate the responses? --------------------------------
cat("\nClade means (Angiosperm vs Gymnosperm) on each response:\n")
print(as.data.frame(sp %>% group_by(clade) %>%
  summarise(across(c(meth, mcra, balance, flux), ~round(median(., na.rm=TRUE), 2)),
            n = n(), .groups = "drop")), row.names = FALSE)

cat("\nREADING: n=2 gymnosperms -> cannot model as a group; interpret descriptively.\n")
cat("Whether trait effects survive angiosperm-only is the key diagnostic above.\n")
sink()
cat(readLines(file.path(out, "gymnosperm_sensitivity.txt")), sep = "\n")
