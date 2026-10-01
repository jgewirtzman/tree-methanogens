source("code/lib/outputs.R")
# ==============================================================================
# 07_mass-basis-sensitivity.R
# ------------------------------------------------------------------------------
# Referee 2 #1: gene copies per gram on a dry vs wet basis. Wood was freeze-dried
# before weighing (dry basis); soil was weighed field-moist (fresh basis). How much
# does the choice of a common basis change the results?
#
#   wood dry -> fresh: exact, from each core's measured water content
#                      (tree_properties: *_moisture_fresh_percent, per tree)
#   soil fresh -> dry: estimated only, from the black-oak soils' gravimetric water
#                      content (no per-sample soil GWC exists)
#
# Tests: heartwood:soil mcrA ratio under each basis; heartwood species ranking,
# dry vs fresh; paired heartwood:sapwood contrast. Ratios of two genes from one
# extract (e.g. mcrA:methanotroph) are basis-invariant and are not re-tested.
#
# Result (2026-09-30): heartwood exceeds soil 14-30x on every basis (never two
# orders of magnitude); species ranking rho = 0.72 dry vs fresh, because core water
# content spans 10-74%; heartwood:sapwood unchanged. Recommendation: dry basis.
#
# NEW file. Output: outputs/audit/mass_basis_sensitivity.txt
# ==============================================================================
suppressMessages(library(dplyr))
norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
L    <- function(x) log10(x + 1)

tp <- read.csv("data/compiled/tree_properties.csv", check.names = FALSE) %>%
  transmute(key = norm(tree_id), mc_in = inner_moisture_fresh_percent / 100,
            mc_out = outer_moisture_fresh_percent / 100) %>%
  group_by(key) %>% summarise(mc_in = mean(mc_in, na.rm = TRUE), mc_out = mean(mc_out, na.rm = TRUE), .groups = "drop")

d <- read.csv("data/compiled/ddpcr_gene_abundances.csv") %>%
  filter(analysis_type == "loose", sample_mass_mg > 0, target_gene %in% c("mcra_probe", "mcra")) %>%
  group_by(sample_id) %>% filter(if (any(target_gene == "mcra_probe")) target_gene == "mcra_probe" else TRUE) %>% ungroup() %>%
  mutate(id = trimws(sub("\n.*", "", sample_id)), key = norm(sub("^[A-Z]{4}_", "", id)), sp = substr(id, 1, 4),
         cpg = concentration_copies_per_uL * 75 / sample_mass_mg * 1000,   # as 04_harmonize_all_data.R
         comp = case_when(material == "Wood" & core_type == "Inner" ~ "Heartwood",
                          material == "Wood" & core_type == "Outer" ~ "Sapwood",
                          material == "Soil" ~ paste0("Soil_", core_type)))

w <- d %>% filter(comp %in% c("Heartwood", "Sapwood")) %>% left_join(tp, by = "key") %>%
  mutate(mc = ifelse(comp == "Heartwood", mc_in, mc_out)) %>% filter(is.finite(mc)) %>%
  mutate(fresh = cpg * (1 - mc))

# soil GWC from the black-oak soils, as in 06_copies-per-gram.R
bo  <- read.csv("data/processed/molecular/black_oak/bo_soil_moisture.csv", check.names = FALSE)
g   <- suppressWarnings(as.numeric(bo$GWC))
gwc <- c(Soil_Organic = median(g[grepl("Organic", bo$`Sample Name`)], na.rm = TRUE),
         Soil_Mineral = median(g[grepl("Mineral", bo$`Sample Name`)], na.rm = TRUE))

sink(out_path("mass_basis_sensitivity.txt"))
cat("MASS BASIS SENSITIVITY (R2 #1)\n\n")
cat(sprintf("wood samples with a per-tree water content: %d\n", nrow(w)))
cat("dry -> fresh factor (1 - water fraction):\n")
print(w %>% group_by(comp) %>% summarise(n = n(), median = round(median(1 - mc), 2),
        p10 = round(quantile(1 - mc, .1), 2), p90 = round(quantile(1 - mc, .9), 2), .groups = "drop") %>%
      as.data.frame(), row.names = FALSE)

cat("\nheartwood : soil, median mcrA copies per g (zeros included)\n")
hw <- w %>% filter(comp == "Heartwood")
for (k in names(gwc)) { sf <- median(d$cpg[d$comp %in% k])
  cat(sprintf("  %-13s wood dry / soil fresh %4.0fx | all fresh %4.0fx | all dry (assumed GWC %.2f) %4.0fx\n",
              k, median(hw$cpg) / sf, median(hw$fresh) / sf, gwc[k], median(hw$cpg) / (sf * (1 + gwc[k])))) }

s <- hw %>% group_by(sp) %>% filter(n() >= 5) %>%
  summarise(n = n(), water = round(median(mc), 2), log_dry = round(median(L(cpg)), 2),
            log_fresh = round(median(L(fresh)), 2), .groups = "drop") %>%
  mutate(rank_dry = rank(-log_dry), rank_fresh = rank(-log_fresh)) %>% arrange(rank_dry)
cat(sprintf("\nheartwood species medians (>=5 samples); Spearman rho dry vs fresh = %.2f\n",
            cor(s$log_dry, s$log_fresh, method = "spearman")))
print(as.data.frame(s), row.names = FALSE)

p <- w %>% group_by(key) %>% filter(all(c("Heartwood", "Sapwood") %in% comp)) %>%
  summarise(dry = L(cpg[comp == "Heartwood"][1]) - L(cpg[comp == "Sapwood"][1]),
            fresh = L(fresh[comp == "Heartwood"][1]) - L(fresh[comp == "Sapwood"][1]), .groups = "drop")
cat(sprintf("\npaired heartwood - sapwood, median log10 difference (n = %d trees): dry %.2f, fresh %.2f\n",
            nrow(p), median(p$dry), median(p$fresh)))
sink()
cat(readLines(out_path("mass_basis_sensitivity.txt")), sep = "\n")
