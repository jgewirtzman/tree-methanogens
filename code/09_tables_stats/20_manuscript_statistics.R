# ==============================================================================
# Manuscript Statistics — All Quantitative Results
# ==============================================================================
# Purpose: Recomputes every quantitative result reported in the manuscript
#   from underlying data. Organized by manuscript section for easy
#   cross-referencing. No figures are generated — console output only.
#
# Usage: Open tree-methanogens.Rproj in RStudio, then:
#   source("code/09_tables_stats/20_manuscript_statistics.R")
#
# Output: Saves all statistics to outputs/manuscript_statistics.txt
#   and prints to console. Ends with a machine-readable summary for
#   easy diff between revisions.
#
# Run time: ~2-3 minutes (phyloseq rarefaction is the bottleneck)
# ==============================================================================

# ==============================================================================
# SECTION 0: SETUP
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(lme4)
library(car)
library(ranger)
library(phyloseq)
library(microeco)

# --- Set up output file ---
output_file <- "outputs/audit/manuscript_statistics.txt"
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
sink(output_file, split = TRUE)  # write to file AND console

cat("\n")
cat(strrep("*", 72), "\n")
cat("  MANUSCRIPT STATISTICS — Tree Microbiomes & Methane\n")
cat("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("*", 72), "\n")

# --- Shared constants ---
NMOL_TO_UG <- 57.744  # nmol m-2 s-1 -> ug CH4 m-2 hr-1

source("code/lib/tables_read.R")   # data/tables/ readers; SPECIES_NAMES = species.csv code -> binomial


# --- Shared utilities ---
section_header <- function(title) {
  cat("\n", strrep("=", 72), "\n")
  cat(" ", title, "\n")
  cat(strrep("=", 72), "\n\n")
}

sub_header <- function(title) {
  cat("---", title, "---\n")
}

stat <- function(label, value, unit = "") {
  if (is.numeric(value)) {
    cat(sprintf("  %-55s %s %s\n", label, format(round(value, 4), nsmall = 4), unit))
  } else {
    cat(sprintf("  %-55s %s %s\n", label, as.character(value), unit))
  }
}

detection_rate <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  100 * sum(x > 0) / length(x)
}

# --- Summary collector for footer ---
summary_stats <- list()
record <- function(name, value) {
  summary_stats[[name]] <<- value
}

# --- Load ALL data upfront ---
cat("\nLoading data...\n")

# The monthly survey as Figure 1 plots it (fig01_temporal-flux.R): every monthly stem
# deployment in data/tables/flux_stem.csv (recovered untagged and dead stems included) and
# every soil deployment in flux_soil.csv. Measurements are excluded only for procedural
# failure: the chamber screen in 02_qc_c0_screen.R (two soil chambers not at ambient at
# closure; flux_soil.qc_pass). The r2 >= 0.7 and -100..200 filters that stood here were
# removed 2026-10-01: an r2 screen is not a detection criterion and inflated the mean stem
# flux 144% (Methods S4). Landscape position: collars.csv (soil), trees.csv (stems).
.collars <- read_table("collars")
semirigid_soil <- read_table("flux_soil") %>% filter(campaign == "monthly_2020_2021", qc_pass) %>%
  transmute(UniqueID = measurement_id, Date = date, Plot.Tag = collar_id,
            Plot.letter = .collars$landscape_position[match(collar_id, .collars$collar_id)],
            CH4_best.flux = CH4_flux_nmol_m2_s)
semirigid_tree <- read_table("flux_stem") %>% filter(campaign == "monthly_2020_2021") %>%
  transmute(UniqueID = measurement_id, Date = date, Plot.Tag = tree_id,
            Plot.Letter = trees_table()$landscape_position[match(tree_id, trees_table()$tree_id)],
            CH4_best.flux = CH4_flux_nmol_m2_s)
stopifnot(!anyNA(semirigid_soil$Plot.letter), !anyNA(semirigid_tree$Plot.Letter))
source("code/lib/soil_collars.R")
moisture_data  <- collar_soil_readings()   # data/tables/soil_readings.csv, one row per collar and date

# 2021 survey per tree (genes, gas, cores, soil holes, flux by height; species and DBH from
# trees.csv) and the 2023 cross-species deployments (species from trees.csv)
ymf2021 <- survey_2021_trees(soil = TRUE)
ymf2023 <- flux_stem_campaign("2023_cross_species")

# CANONICAL BUDGET is authoritative for every stand-level number (2026-07-30).
# This script is named "manuscript statistics" and was quoting the LEGACY
# ANNUAL_SUMMARY.csv / MONTHLY_FLUXES.csv, which come from the old upscaling block
# inside 02_rf_models.R -- built on a 7,010-stem inventory that omits the entire bytag
# survey, with a two-rung species ladder that buckets all 1,899 Kalmia into
# SPECIES_OTHER, and annualised with 30.4 days/month. It reported tree 2.725 and soil
# -543.43 against the canonical 4.378 and -422.16: 38% and 29% apart. Those legacy
# files are still read below for the monthly SHAPE only; every scalar now comes from
# canonical_budget.csv.
monthly_fluxes <- read.csv("outputs/tables/MONTHLY_FLUXES.csv")
annual_summary <- read.csv("outputs/tables/ANNUAL_SUMMARY.csv")
CANON <- local({
  f <- "outputs/data/canonical_budget.csv"
  if (!file.exists(f))
    stop("20_manuscript_statistics.R needs ", f,
         " -- run code/06_upscale/03_budget_canonical.R first")
  b <- read.csv(f, stringsAsFactors = FALSE)
  setNames(b$value, b$quantity)
})
canon <- function(q) { if (!q %in% names(CANON)) stop("canonical_budget.csv lacks ", q); unname(CANON[[q]]) }
load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")  # tree_train_complete, X_tree, X_soil, soil_train_complete

cat("  All data loaded.\n")


# ==============================================================================
# SECTION 1: TEMPORAL FLUX PATTERNS (Figure 1)
# ==============================================================================

section_header("SECTION 1: TEMPORAL FLUX PATTERNS (Figure 1)")

# Quality filter (matches 06_soil_tree_timeseries.R)
soil_qc <- semirigid_soil %>%
  filter(!is.na(CH4_best.flux)) %>%
  mutate(Plot_Type = case_when(
    Plot.letter %in% c("WD", "WS") ~ "W",
    Plot.letter == "I" ~ "I",
    Plot.letter == "U" ~ "U",
    TRUE ~ Plot.letter
  ))

tree_qc <- semirigid_tree %>%
  filter(!is.na(CH4_best.flux)) %>%
  mutate(
    Plot_Type = case_when(
      Plot.Letter %in% c("WD", "WS") ~ "W",
      Plot.Letter == "I" ~ "I",
      Plot.Letter == "U" ~ "U",
      TRUE ~ Plot.Letter
    ),
    CH4_flux = CH4_best.flux
  )

# Sample sizes
stat("Total tree flux measurements", nrow(tree_qc))
stat("Unique trees in temporal survey", length(unique(tree_qc$Plot.Tag)))
stat("Total soil measurements", nrow(soil_qc))
stat("Unique soil plots", length(unique(soil_qc$Plot.Tag)))

# Tree fluxes by landscape position
sub_header("Tree fluxes (nmol m-2 s-1)")
for (pos in c("U", "I", "W")) {
  d <- tree_qc %>% filter(Plot_Type == pos)
  pos_label <- c(U = "Upland", I = "Intermediate", W = "Transitional wetland")[pos]
  mean_val <- mean(d$CH4_flux, na.rm = TRUE)
  cat(sprintf("\n  %s (n=%d):\n", pos_label, nrow(d)))
  stat("    Mean", mean_val, "nmol m-2 s-1")
  stat("    Mean (ug CH4 m-2 hr-1)", mean_val * NMOL_TO_UG)
  stat("    Range", paste(round(min(d$CH4_flux, na.rm = TRUE) * NMOL_TO_UG, 1),
                          "to",
                          round(max(d$CH4_flux, na.rm = TRUE) * NMOL_TO_UG, 1)), "ug CH4 m-2 hr-1")
  stat("    % positive", round(100 * sum(d$CH4_flux > 0, na.rm = TRUE) / nrow(d), 1), "%")
  record(paste0(tolower(pos_label), "_tree_mean_nmol"), mean_val)
}

# Soil fluxes by landscape position
sub_header("Soil fluxes (nmol m-2 s-1)")
for (pos in c("U", "I", "W")) {
  d <- soil_qc %>% filter(Plot_Type == pos)
  pos_label <- c(U = "Upland", I = "Intermediate", W = "Transitional wetland")[pos]
  mean_val <- mean(d$CH4_best.flux, na.rm = TRUE)
  cat(sprintf("\n  %s (n=%d):\n", pos_label, nrow(d)))
  stat("    Mean", mean_val, "nmol m-2 s-1")
  stat("    Mean (ug CH4 m-2 hr-1)", mean_val * NMOL_TO_UG)
  stat("    Range", paste(round(min(d$CH4_best.flux, na.rm = TRUE) * NMOL_TO_UG, 1),
                          "to",
                          round(max(d$CH4_best.flux, na.rm = TRUE) * NMOL_TO_UG, 1)), "ug CH4 m-2 hr-1")
  if (pos == "W") {
    stat("    Median", median(d$CH4_best.flux, na.rm = TRUE), "nmol m-2 s-1")
    stat("    Median (ug CH4 m-2 hr-1)", median(d$CH4_best.flux, na.rm = TRUE) * NMOL_TO_UG)
    stat("    % negative", round(100 * sum(d$CH4_best.flux < 0, na.rm = TRUE) / nrow(d), 1), "%")
  }
  record(paste0(tolower(pos_label), "_soil_mean_nmol"), mean_val)
}

# One-sample t-test: upland tree flux > 0
upland_tree <- tree_qc %>% filter(Plot_Type == "U")
tt <- t.test(upland_tree$CH4_flux, mu = 0, alternative = "greater")
sub_header("One-sample t-test: Upland tree flux > 0")
stat("t-statistic", tt$statistic)
stat("p-value", tt$p.value)
record("upland_tree_ttest_p", tt$p.value)

# VWC by position
sub_header("Volumetric Water Content by position")
moisture_data <- moisture_data %>%
  mutate(VWC = soil_vwc_pct,
         Plot_Type = case_when(
    plot_letter == "U" ~ "U",
    plot_letter == "I" ~ "I",
    plot_letter %in% c("WD", "WS") ~ "W",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Plot_Type), !is.na(VWC))

for (pos in c("U", "I", "W")) {
  d <- moisture_data %>% filter(Plot_Type == pos)
  pos_label <- c(U = "Upland", I = "Intermediate", W = "Transitional wetland")[pos]
  cat(sprintf("  %s: %.1f +/- %.1f%% (n=%d)\n",
              pos_label, mean(d$VWC, na.rm = TRUE), sd(d$VWC, na.rm = TRUE), nrow(d)))
}

# Soil temperature range: the collar readings (soil_readings.csv), as Figure 1 plots them
stat("Soil temperature range (collar readings)",
     paste(round(min(moisture_data$soil_temp_C, na.rm = TRUE), 1), "to",
           round(max(moisture_data$soil_temp_C, na.rm = TRUE), 1)), "deg C")


# ==============================================================================
# SECTION 2: HEIGHT EFFECTS (Figure 2)
# ==============================================================================

section_header("SECTION 2: HEIGHT EFFECTS (Figure 2)")

# Prepare data (matches 01_height_effects.R): 2021 multi-height deployments, species and
# landscape position (plot) from trees.csv; the root crown has no height_cm
height_data <- flux_stem_campaign("2021_multiheight") %>%
  transmute(UniqueID = measurement_id, tree_id, species = species_code,
            plot = trees_table()$landscape_position[match(tree_id, trees_table()$tree_id)],
            height_numeric = height_cm, CH4_best.flux = CH4_flux_nmol_m2_s) %>%
  filter(!is.na(CH4_best.flux), !is.na(height_numeric)) %>%
  mutate(
    height_m = height_numeric / 100,
    tree_unique = paste(plot, tree_id, sep = "_"),
    species_latin = SPECIES_NAMES[species]
  ) %>%
  filter(species != "KALA")

stat("Total measurements in height analysis", nrow(height_data))
stat("Total individual trees", length(unique(height_data$tree_unique)))
stat("Total species analyzed", length(unique(height_data$species)))

# Per-species LME models
test_height <- function(sp, data) {
  d <- data %>% filter(species == sp)
  n_multi <- d %>% group_by(tree_unique) %>%
    summarise(nh = n_distinct(height_numeric)) %>% filter(nh > 1) %>% nrow()
  if (n_multi < 3) return(data.frame(species = sp, coef = NA, se = NA, p = NA, n = nrow(d)))
  tryCatch({
    m <- lmer(CH4_best.flux ~ height_numeric + (1 | tree_unique), data = d)
    ct <- summary(m)$coefficients
    data.frame(species = sp, coef = ct["height_numeric", "Estimate"],
               se = ct["height_numeric", "Std. Error"],
               p = 2 * (1 - pnorm(abs(ct["height_numeric", "t value"]))),
               n = nrow(d))
  }, error = function(e) data.frame(species = sp, coef = NA, se = NA, p = NA, n = nrow(d)))
}

sp_list <- sort(unique(height_data$species))
height_results <- do.call(rbind, lapply(sp_list, test_height, height_data))

n_tested <- sum(!is.na(height_results$p))
n_sig <- sum(height_results$p < 0.05, na.rm = TRUE)
stat("Species with significant height effects", paste(n_sig, "of", n_tested))
record("height_nsig_of_ntested", paste0(n_sig, "/", n_tested))

sub_header("Significant species (nmol m-2 s-1 per METER)")
sig_sp <- height_results %>% filter(p < 0.05)
for (i in seq_len(nrow(sig_sp))) {
  cat(sprintf("  %s: coef = %.3f, p = %.4f\n",
              SPECIES_NAMES[sig_sp$species[i]],
              sig_sp$coef[i] * 100, sig_sp$p[i]))  # *100 converts /cm to /m
}

# Jackknife robustness
sub_header("Jackknife robustness (leave-one-tree-out)")
for (sp in c("BEAL", "TSCA", "ACSA")) {
  d <- height_data %>% filter(species == sp)
  trees <- unique(d$tree_unique)
  if (length(trees) < 4) next
  orig_m <- lmer(CH4_best.flux ~ height_numeric + (1 | tree_unique), data = d)
  orig_coef <- summary(orig_m)$coefficients["height_numeric", "Estimate"]
  jk <- data.frame()
  for (tr in trees) {
    sub <- d %>% filter(tree_unique != tr)
    tryCatch({
      m <- lmer(CH4_best.flux ~ height_numeric + (1 | tree_unique), data = sub)
      ct <- summary(m)$coefficients
      jk <- rbind(jk, data.frame(coef = ct["height_numeric", "Estimate"],
                                  p = 2 * (1 - pnorm(abs(ct["height_numeric", "t value"])))))
    }, error = function(e) {})
  }
  if (nrow(jk) > 0) {
    prop_sig <- round(100 * sum(jk$p < 0.05) / nrow(jk), 0)
    cat(sprintf("  %s: %d%% iterations significant\n", SPECIES_NAMES[sp], prop_sig))
  }
  # Within-species: proportion of trees with negative slopes
  tree_slopes <- d %>% group_by(tree_unique) %>% filter(n() >= 2) %>%
    summarise(slope = coef(lm(CH4_best.flux ~ height_numeric))[2], .groups = "drop")
  prop_neg <- round(100 * sum(tree_slopes$slope < 0) / nrow(tree_slopes), 0)
  cat(sprintf("  %s: %d/%d trees with negative slopes (%d%%)\n",
              SPECIES_NAMES[sp], sum(tree_slopes$slope < 0), nrow(tree_slopes), prop_neg))
}

# B. alleghaniensis flux at specific heights
beal <- height_data %>% filter(species == "BEAL")
for (h in c(50, 200)) {
  d <- beal %>% filter(height_numeric == h)
  if (nrow(d) > 0) {
    cat(sprintf("  BEAL flux at %d cm: %.2f +/- %.2f nmol m-2 s-1 (n=%d)\n",
                h, mean(d$CH4_best.flux), sd(d$CH4_best.flux) / sqrt(nrow(d)), nrow(d)))
  }
}

# ACSA outlier analysis
acsa <- height_data %>% filter(species == "ACSA")
acsa_200 <- acsa %>% filter(height_numeric == 200) %>% arrange(desc(CH4_best.flux))
if (nrow(acsa_200) > 0) {
  outlier_tree <- acsa_200$tree_unique[1]
  stat("ACSA outlier flux at 2m", acsa_200$CH4_best.flux[1], "nmol m-2 s-1")
  acsa_no_outlier <- acsa %>% filter(tree_unique != outlier_tree)
  m <- lmer(CH4_best.flux ~ height_numeric + (1 | tree_unique), data = acsa_no_outlier)
  ct <- summary(m)$coefficients
  p_no_outlier <- 2 * (1 - pnorm(abs(ct["height_numeric", "t value"])))
  cat(sprintf("  ACSA after outlier removal: coef = %.3f per m, p = %.3f\n",
              ct["height_numeric", "Estimate"] * 100, p_no_outlier))
}

# Soil condition correlations
sub_header("Soil condition correlations")
soil_by_sp <- ymf2021 %>%
  filter(species_id %in% unique(height_data$species)) %>%
  group_by(species_id) %>%
  summarise(mean_VWC = mean(VWC_mean, na.rm = TRUE),
            mean_mcra_mineral = mean(ddpcr_mcra_probe_Mineral_loose, na.rm = TRUE),
            .groups = "drop") %>%
  rename(species = species_id)

sp_env <- height_results %>%
  left_join(soil_by_sp, by = "species") %>%
  filter(!is.na(coef), !is.na(mean_VWC))

sp_env$log_mcra <- log10(sp_env$mean_mcra_mineral + 1)

vwc_cor <- cor.test(sp_env$coef, sp_env$mean_VWC)
mcra_cor <- cor.test(sp_env$coef, sp_env$log_mcra)
cat(sprintf("  Height coef vs VWC: r = %.2f, p = %.4f\n", vwc_cor$estimate, vwc_cor$p.value))
cat(sprintf("  Height coef vs log(mcrA+1): r = %.2f, p = %.4f\n", mcra_cor$estimate, mcra_cor$p.value))

vwc_mcra <- ymf2021 %>%
  filter(!is.na(VWC_mean), !is.na(ddpcr_mcra_probe_Mineral_loose)) %>%
  mutate(log_mcra = log10(ddpcr_mcra_probe_Mineral_loose + 1))
vwc_mcra_cor <- cor.test(vwc_mcra$VWC_mean, vwc_mcra$log_mcra)
cat(sprintf("  VWC vs log(mcrA+1) across trees: r = %.2f, p = %.4f (n=%d)\n",
            vwc_mcra_cor$estimate, vwc_mcra_cor$p.value, nrow(vwc_mcra)))

# Soil conditions for BEAL specifically
beal_env <- sp_env %>% filter(species == "BEAL")
if (nrow(beal_env) > 0) {
  stat("BEAL soil VWC", beal_env$mean_VWC, "%")
  stat("BEAL soil mcrA (log10)", beal_env$log_mcra)
}


# ==============================================================================
# SECTION 3: VARIANCE PARTITIONING (Figure 3)
# ==============================================================================

section_header("SECTION 3: VARIANCE PARTITIONING (Figure 3)")

# The partition and species means reported in Results §3 come from the Figure 3
# script (one growing-season breast-height value per tree, 472 trees, inverse-variance
# weighted), printed to outputs/logs/fig03_variance-partition.txt. The earlier
# per-measurement version that stood here (476 rows, a different model) disagreed
# with the figure and the text, so it was removed (2026-10-05) rather than kept as a
# second, conflicting set of numbers.
cat("  See outputs/logs/fig03_variance-partition.txt (Figure 3 script) for the partition\n")
cat("  and species means quoted in Results §3.\n")


# ==============================================================================
# SECTION 4: GENE ABUNDANCE ACROSS COMPARTMENTS (Figure 4)
# ==============================================================================

section_header("SECTION 4: GENE ABUNDANCE ACROSS COMPARTMENTS (Figure 4)")

gene_stats <- function(values, label) {
  values <- values[!is.na(values)]
  n_total <- length(values)
  n_pos <- sum(values > 0)
  pos_vals <- values[values > 0]
  cat(sprintf("\n  %s (n=%d):\n", label, n_total))
  stat("    Detection rate", round(100 * n_pos / n_total, 1), "%")
  if (length(pos_vals) > 0) {
    stat("    Mean (log10, among positive)", round(mean(log10(pos_vals)), 2))
    stat("    Median (log10, among positive)", round(median(log10(pos_vals)), 2))
    stat("    Max (log10)", round(max(log10(pos_vals)), 2))
  }
}

sub_header("mcrA (probe, loose) by compartment")
gene_stats(ymf2021$ddpcr_mcra_probe_Inner_loose, "Heartwood")
gene_stats(ymf2021$ddpcr_mcra_probe_Outer_loose, "Sapwood")
gene_stats(ymf2021$ddpcr_mcra_probe_Mineral_loose, "Mineral soil")
gene_stats(ymf2021$ddpcr_mcra_probe_Organic_loose, "Organic soil")

# Share of samples at or above the assay's limit of detection: 3 copies per reaction
# (6/6 replicates; Arnold et al. 2024), i.e. 0.15 copies per uL of a 20 uL partitioned reaction.
# Judged per sample on the reaction scale, so it does not depend on sample mass; in recovered
# copies per gram (scaled to the whole lysate) it is ~3,600 for 100 mg of wood and ~1,900 for 250 mg of soil.
source("code/lib/ddpcr_constants.R")
LOD_CONC <- 3 / 20
sub_header("mcrA (probe, loose): samples at or above the limit of detection (3 copies per reaction)")
# every 2021 probe-mcrA run that passed QC (data/tables/ddpcr.csv; position and mass from samples.csv)
cg <- read_table("ddpcr") %>% filter(gene == "mcrA", assay == "probe", qc_pass, !is.na(copies_per_uL_rxn)) %>%
  inner_join(read_table("samples") %>% filter(campaign == "2021_survey") %>%
               transmute(sample_id, core_type = c(inner = "Inner", outer = "Outer", organic = "Organic",
                                                  mineral = "Mineral")[position], sample_mass_mg = extraction_mass_mg),
             by = "sample_id")
for (ct in c("Inner", "Outer", "Mineral", "Organic")) {
  v <- cg$copies_per_uL_rxn[cg$core_type == ct]; m <- cg$sample_mass_mg[cg$core_type == ct]
  stat(sprintf("    %s: %% at or above (n = %d)", ct, length(v)), round(100 * mean(v >= LOD_CONC), 1), "%")
  stat(sprintf("    %s: limit in copies/g at median mass %.0f mg", ct, median(m, na.rm = TRUE)),
       round(ddpcr_copies_per_g(LOD_CONC, median(m, na.rm = TRUE), ct)))
  record(sprintf("mcra_%s_pct_above_lod", tolower(ct)), round(100 * mean(v >= LOD_CONC), 1))
}

# Record key values
hw_mcra <- ymf2021$ddpcr_mcra_probe_Inner_loose[!is.na(ymf2021$ddpcr_mcra_probe_Inner_loose)]
record("mcra_heartwood_detection_pct", round(100 * sum(hw_mcra > 0) / length(hw_mcra), 1))

# Wilcoxon heartwood vs. sapwood
hw <- ymf2021$ddpcr_mcra_probe_Inner_loose
sw <- ymf2021$ddpcr_mcra_probe_Outer_loose
both <- !is.na(hw) & !is.na(sw)
if (sum(both) > 0) {
  wt <- wilcox.test(hw[both], sw[both])
  sub_header("Wilcoxon: heartwood vs sapwood mcrA")
  stat("W statistic", wt$statistic)
  stat("p-value", wt$p.value)
}

# Sample sizes per compartment
sub_header("Sample sizes per compartment")
stat("Heartwood (Inner)", sum(!is.na(ymf2021$ddpcr_mcra_probe_Inner_loose)))
stat("Sapwood (Outer)", sum(!is.na(ymf2021$ddpcr_mcra_probe_Outer_loose)))
stat("Mineral soil", sum(!is.na(ymf2021$ddpcr_mcra_probe_Mineral_loose)))
stat("Organic soil", sum(!is.na(ymf2021$ddpcr_mcra_probe_Organic_loose)))

# Methanotrophs by compartment
sub_header("Methanotroph abundance (pmoA + mmoX) by compartment")
for (loc in c("Inner", "Outer", "Mineral", "Organic")) {
  pmoa_col <- paste0("ddpcr_pmoa_", loc, "_loose")
  mmox_col <- paste0("ddpcr_mmox_", loc, "_loose")
  if (pmoa_col %in% names(ymf2021) && mmox_col %in% names(ymf2021)) {
    p <- ymf2021[[pmoa_col]]
    m <- ymf2021[[mmox_col]]
    total <- p + m          # both genes measured; NA where either is missing
    pos <- total[!is.na(total) & total > 0]
    label <- c(Inner = "Heartwood", Outer = "Sapwood",
               Mineral = "Mineral soil", Organic = "Organic soil")[loc]
    if (length(pos) > 0)
      cat(sprintf("  %s: mean log10 = %.2f (n=%d positive of %d total)\n",
                  label, mean(log10(pos)), length(pos), sum(!is.na(total))))
  }
}

# pmoA and mmoX detection rates separately
sub_header("pmoA and mmoX detection rates (separately)")
for (gene in c("pmoa", "mmox")) {
  for (loc in c("Inner", "Outer", "Mineral", "Organic")) {
    col <- paste0("ddpcr_", gene, "_", loc, "_loose")
    if (col %in% names(ymf2021)) {
      vals <- ymf2021[[col]][!is.na(ymf2021[[col]])]
      label_loc <- c(Inner = "HW", Outer = "SW", Mineral = "Min", Organic = "Org")[loc]
      cat(sprintf("  %s %s: %.1f%% detected (n=%d)\n",
                  toupper(gene), label_loc, detection_rate(vals), length(vals)))
    }
  }
}


# ==============================================================================
# SECTION 5: 16S COMMUNITY COMPOSITION (Figure 5)
# ==============================================================================

section_header("SECTION 5: 16S COMMUNITY COMPOSITION (Figure 5)")

cat("Building phyloseq object...\n")

# Replicate phyloseq build from 08c_combined_methane_cycling_composition.R

# 16S phyloseq assembly: single definition in code/lib/build_phyloseq.R.
# All nine copies of this ~45-line block produced an identical object
# (71,511 taxa / 587 samples / 6,427,747 counts) before consolidation.
source("code/lib/build_phyloseq.R")
.ps16   <- build_16s_phyloseq()
no_mito <- .ps16$ps
phylo_tree <- .ps16$tree
ddpcr      <- .ps16$ddpcr
taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))

set.seed(46814)
ps.rare <- rarefy_even_depth(no_mito, sample.size = 3500)
colnames(tax_table(ps.rare)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
ps.ra <- transform_sample_counts(ps.rare, function(x) x / sum(x) * 100)

# Filter samples
samp_df <- data.frame(sample_data(ps.ra), stringsAsFactors = FALSE)
samp_df$compartment <- case_when(
  samp_df$core_type == "Inner" ~ "Heartwood",
  samp_df$core_type == "Outer" ~ "Sapwood",
  samp_df$core_type == "Mineral" ~ "Mineral Soil",
  samp_df$core_type == "Organic" ~ "Organic Soil",
  TRUE ~ NA_character_
)
samp_df <- samp_df %>% filter(!is.na(compartment), !is.na(species.x), species.x != "")
sample_data(ps.ra) <- sample_data(samp_df)
ps.filt <- prune_samples(!is.na(sample_data(ps.ra)$compartment), ps.ra)

otu_df <- as.data.frame(otu_table(ps.filt))
tax_df <- as.data.frame(tax_table(ps.filt))
samp_meta <- data.frame(sample_data(ps.filt))

# Methanogen families (code/lib/methanogen_families.R)
source("code/lib/methanogen_families.R")
methanogen_families <- METHANOGEN_FAMILIES

methanogen_asvs <- rownames(tax_df)[tax_df$Family %in% methanogen_families]
methanogen_asvs <- intersect(methanogen_asvs, rownames(otu_df))
cat(sprintf("  Found %d methanogen ASVs in rarefied OTU table\n", length(methanogen_asvs)))
methanogen_pct <- colSums(otu_df[methanogen_asvs, , drop = FALSE])

samp_meta$methanogen_pct <- methanogen_pct[rownames(samp_meta)]

# Methanogen relative abundance by compartment
sub_header("Methanogen relative abundance (% of community)")
for (comp in c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil")) {
  d <- samp_meta %>% filter(compartment == comp)
  vals <- d$methanogen_pct
  cat(sprintf("  %s: %.2f +/- %.2f%% (range: %.2f to %.2f%%, n=%d)\n",
              comp, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
              min(vals, na.rm = TRUE), max(vals, na.rm = TRUE), length(vals)))
}

# ANOVA for species effect
sub_header("ANOVA for species effect on methanogen abundance")
hw_samples <- samp_meta %>% filter(compartment == "Heartwood")
anova_16s <- aov(methanogen_pct ~ species.x, data = hw_samples)
anova_p_16s <- summary(anova_16s)[[1]]["species.x", "Pr(>F)"]
cat(sprintf("  16S ANOVA p-value: %.4f\n", anova_p_16s))

# A. saccharum heartwood
acsa_hw <- hw_samples %>% filter(species.x == "ACSA")
cat(sprintf("  A. saccharum heartwood methanogen: %.1f +/- %.1f%% (n=%d)\n",
            mean(acsa_hw$methanogen_pct), sd(acsa_hw$methanogen_pct), nrow(acsa_hw)))

# Top methanogen families by abundance (VERIFY TAXA NAMES)
sub_header("Top methanogen families by mean relative abundance (heartwood)")
hw_family_means <- data.frame()
for (fam in methanogen_families) {
  asvs <- rownames(tax_df)[tax_df$Family == fam]
  asvs <- intersect(asvs, rownames(otu_df))
  if (length(asvs) > 0) {
    pcts <- colSums(otu_df[asvs, rownames(hw_samples), drop = FALSE])
    hw_family_means <- rbind(hw_family_means,
                              data.frame(Family = fam, mean_pct = mean(pcts),
                                         sd_pct = sd(pcts), max_pct = max(pcts)))
  }
}
hw_family_means <- hw_family_means %>% arrange(desc(mean_pct))
for (i in seq_len(nrow(hw_family_means))) {
  cat(sprintf("  %s: %.2f +/- %.2f%% (max %.2f%%)\n",
              hw_family_means$Family[i], hw_family_means$mean_pct[i],
              hw_family_means$sd_pct[i], hw_family_means$max_pct[i]))
}

# Sapwood family abundances
sub_header("Key methanogen families in sapwood")
sw_samples <- samp_meta %>% filter(compartment == "Sapwood")
for (fam in c("Methanobacteriaceae", "Methanomassiliicoccaceae")) {
  asvs <- rownames(tax_df)[tax_df$Family == fam]
  asvs <- intersect(asvs, rownames(otu_df))
  if (length(asvs) > 0) {
    pcts <- colSums(otu_df[asvs, rownames(sw_samples), drop = FALSE])
    cat(sprintf("  %s: range %.2f to %.2f%%\n", fam, min(pcts), max(pcts)))
  }
}

# Tests behind Results §5 (heartwood vs sapwood families; A. saccharum against other species)
sub_header("Methanogen families: heartwood vs sapwood (Wilcoxon) and soil means")
fam_pct <- function(fam, ids) {
  asvs <- intersect(rownames(tax_df)[tax_df$Family == fam], rownames(otu_df))
  colSums(otu_df[asvs, ids, drop = FALSE])
}
for (fam in c("Methanobacteriaceae", "Methanomassiliicoccaceae")) {
  w <- wilcox.test(fam_pct(fam, rownames(hw_samples)), fam_pct(fam, rownames(sw_samples)))
  cat(sprintf("  %s heartwood vs sapwood: Wilcoxon p = %.2g\n", fam, w$p.value))
  for (comp in c("Mineral Soil", "Organic Soil")) {
    v <- fam_pct(fam, rownames(samp_meta)[samp_meta$compartment == comp])
    cat(sprintf("  %s, %s: %.3f +/- %.2f%%\n", fam, comp, mean(v), sd(v)))
  }
}
sub_header("A. saccharum heartwood against other species (Wilcoxon)")
w16 <- wilcox.test(hw_samples$methanogen_pct[hw_samples$species.x == "ACSA"],
                   hw_samples$methanogen_pct[hw_samples$species.x != "ACSA"])
cat(sprintf("  16S methanogen share: p = %.2g\n", w16$p.value))
hw_mcra_sp <- ymf2021 %>% filter(!is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  transmute(acsa = species_id == "ACSA", lg = log10(ddpcr_mcra_probe_Inner_loose + 1),
            pos = ddpcr_mcra_probe_Inner_loose > 0)
wm <- wilcox.test(lg ~ acsa, data = hw_mcra_sp)
acsa_pos <- hw_mcra_sp %>% filter(acsa, pos)
cat(sprintf("  heartwood mcrA, log10 among positive: %.1f +/- %.1f (SD, n = %d); vs other species (log10(x+1), all samples): p = %.2g\n",
            mean(acsa_pos$lg), sd(acsa_pos$lg), nrow(acsa_pos), wm$p.value))
pos_only <- hw_mcra_sp %>% filter(pos)
cat(sprintf("  heartwood mcrA among positive, vs other species: Welch t p = %.2g; Wilcoxon p = %.2g\n",
            t.test(lg ~ acsa, data = pos_only)$p.value, wilcox.test(lg ~ acsa, data = pos_only)$p.value))

# Soil taxa (VERIFY Bathyarchaeia dominance)
sub_header("Top soil archaeal taxa (Mineral + Organic)")
soil_samples <- samp_meta %>% filter(compartment %in% c("Mineral Soil", "Organic Soil"))
all_classes <- tax_df$Class
archaea_asvs <- rownames(tax_df)[tax_df$Kingdom == "Archaea"]
archaea_asvs <- intersect(archaea_asvs, rownames(otu_df))
if (length(archaea_asvs) > 0) {
  archaea_class_pcts <- data.frame()
  for (cls in unique(tax_df[archaea_asvs, "Class"])) {
    if (is.na(cls)) next
    asvs <- archaea_asvs[tax_df[archaea_asvs, "Class"] == cls]
    asvs <- intersect(asvs, rownames(otu_df))
    if (length(asvs) == 0) next
    pcts <- colSums(otu_df[asvs, rownames(soil_samples), drop = FALSE])
    archaea_class_pcts <- rbind(archaea_class_pcts,
                                 data.frame(Class = cls, mean_pct = mean(pcts)))
  }
  archaea_class_pcts <- archaea_class_pcts %>% arrange(desc(mean_pct))
  for (i in 1:min(5, nrow(archaea_class_pcts))) {
    cat(sprintf("  %s: %.2f%%\n", archaea_class_pcts$Class[i], archaea_class_pcts$mean_pct[i]))
  }
}

# Families correlated with mcrA (VERIFY taxa names)
sub_header("Top families correlated with mcrA (heartwood)")
hw_with_mcra <- hw_samples %>% filter(!is.na(mcra_probe_strict))
if (nrow(hw_with_mcra) > 5) {
  all_families <- unique(tax_df$Family)
  all_families <- all_families[!is.na(all_families) & all_families != ""]
  fam_cors <- data.frame()
  for (fam in all_families) {
    asvs <- rownames(tax_df)[tax_df$Family == fam]
    asvs <- intersect(asvs, rownames(otu_df))
    if (length(asvs) > 0) {
      pcts <- colSums(otu_df[asvs, rownames(hw_with_mcra), drop = FALSE])
      sd_val <- sd(pcts, na.rm = TRUE)
      if (!is.na(sd_val) && sd_val > 0) {
        ct <- cor.test(pcts, log10(hw_with_mcra$mcra_probe_strict + 1))
        fam_cors <- rbind(fam_cors, data.frame(Family = fam, r = ct$estimate, p = ct$p.value))
      }
    }
  }
  fam_cors <- fam_cors %>% filter(p < 0.05) %>% arrange(desc(abs(r)))
  cat("  Top 10 positively correlated families:\n")
  top_pos <- fam_cors %>% filter(r > 0) %>% head(10)
  for (i in seq_len(nrow(top_pos))) {
    cat(sprintf("    %s: r = %.3f, p = %.4f\n", top_pos$Family[i], top_pos$r[i], top_pos$p[i]))
  }
}

# Methanotroph composition (VERIFY taxonomy)
sub_header("Methanotroph families detected (known/putative)")
# One rule for every count: the shared classifier (genus-aware; a mixed-family ASV is
# Putative only when its genus is unresolved). This section previously classified inline
# and disagreed with Figure 5 and the clade census (2026-10-02).
source("code/lib/load_methanotroph_definitions.R")
mt_defs <- load_methanotroph_defs()
mt_status <- classify_methanotrophs(tax_df, mt_defs)
known_asvs   <- intersect(rownames(tax_df)[mt_status %in% "Known"], rownames(otu_df))
putative_only <- intersect(rownames(tax_df)[mt_status %in% "Putative"], rownames(otu_df))
{
  tt <- read.csv("outputs/data/methanotroph_tier_totals.csv")
  cat("  Tier totals over the analysed samples (03_clade-census.R):\n")
  for (r in seq_len(nrow(tt))) cat(sprintf("    %-17s %4d ASVs; heartwood %.3f%%, sapwood %.3f%%, mineral %.3f%%, organic %.3f%%\n",
      tt$tier[r], tt$ASVs_analysed[r], tt$Heartwood_pct[r], tt$Sapwood_pct[r], tt$Mineral_pct[r], tt$Organic_pct[r]))
}

cat(sprintf("  Known methanotroph ASVs: %d\n", length(known_asvs)))
cat(sprintf("  Putative methanotroph ASVs: %d\n", length(putative_only)))

# List families detected
all_mt_asvs <- unique(c(known_asvs, putative_only))
mt_families_detected <- sort(unique(tax_df[all_mt_asvs, "Family"]))
mt_families_detected <- mt_families_detected[!is.na(mt_families_detected)]
cat("  Detected methanotroph families:\n")
for (fam in mt_families_detected) {
  asvs <- intersect(all_mt_asvs, rownames(tax_df)[tax_df$Family == fam])
  asvs <- intersect(asvs, rownames(otu_df))
  if (length(asvs) == 0) next
  pcts <- colSums(otu_df[asvs, rownames(hw_samples), drop = FALSE])
  cat(sprintf("    %s: mean %.3f%% in heartwood (n_ASVs=%d)\n", fam, mean(pcts), length(asvs)))
}


# ==============================================================================
# SECTION 6: FAPROTAX / PICRUSt (Figures S3-S5)
# ==============================================================================

section_header("SECTION 6: FAPROTAX / PICRUSt (Figures S3-S5)")

# FAPROTAX — computed via microeco (same pipeline as figS10_faprotax.R)
sub_header("FAPROTAX functional predictions (via microeco)")

# Prepare phyloseq for microeco: need QIIME-style taxonomy prefixes
ps.wood <- subset_samples(ps.filt, material == "Wood")
ps.wood.fap <- ps.wood
tax_table(ps.wood.fap)[, 1] <- paste0("k__", tax_table(ps.wood.fap)[, 1])
tax_table(ps.wood.fap)[, 2] <- paste0("p__", tax_table(ps.wood.fap)[, 2])
tax_table(ps.wood.fap)[, 3] <- paste0("c__", tax_table(ps.wood.fap)[, 3])
tax_table(ps.wood.fap)[, 4] <- paste0("o__", tax_table(ps.wood.fap)[, 4])
tax_table(ps.wood.fap)[, 5] <- paste0("f__", tax_table(ps.wood.fap)[, 5])
tax_table(ps.wood.fap)[, 6] <- paste0("g__", tax_table(ps.wood.fap)[, 6])
tax_table(ps.wood.fap)[, 7] <- paste0("s__", tax_table(ps.wood.fap)[, 7])

# Panel (a): Merge by core_type (Inner = Heartwood, Outer = Sapwood)
meco_all <- microtable$new(
  otu_table = as.data.frame(otu_table(ps.wood.fap)),
  tax_table = noquote(as.data.frame(tax_table(ps.wood.fap))),
  sample_table = as.data.frame(as.matrix(sample_data(ps.wood.fap))),
  phylo_tree = phy_tree(ps.wood.fap)
)
meco_core <- clone(meco_all)$merge_samples("core_type")
cat("Running FAPROTAX for core type grouping...\n")
t_func_core <- trans_func$new(meco_core)
t_func_core$cal_spe_func(prok_database = "FAPROTAX")
t_func_core$cal_spe_func_perc(abundance_weighted = TRUE, dec = 2)
perc_core <- t_func_core$res_spe_func_perc

# Report compartment-level percentages for key functions
key_funcs <- c("fermentation", "dark_hydrogen_oxidation", "methylotrophy",
               "methanogenesis", "methanogenesis_by_CO2_reduction_with_H2",
               "methanogenesis_by_reduction_of_methyl_compounds_with_H2",
               "methanol_oxidation", "methanotrophy", "aerobic_chemoheterotrophy")
for (func in key_funcs) {
  if (func %in% colnames(perc_core)) {
    hw_val <- perc_core["Inner", func]
    sw_val <- perc_core["Outer", func]
    cat(sprintf("  %s: HW=%.2f%%, SW=%.2f%%\n", func, hw_val, sw_val))
  }
}

# Panel (b): Heartwood only, per species
ps.hw.fap <- subset_samples(ps.wood.fap, core_type == "Inner")
meco_hw <- microtable$new(
  otu_table = as.data.frame(otu_table(ps.hw.fap)),
  tax_table = noquote(as.data.frame(tax_table(ps.hw.fap))),
  sample_table = as.data.frame(as.matrix(sample_data(ps.hw.fap))),
  phylo_tree = phy_tree(ps.hw.fap)
)
meco_species <- clone(meco_hw)$merge_samples("species.x")
cat("Running FAPROTAX for species grouping (heartwood)...\n")
t_func_sp <- trans_func$new(meco_species)
t_func_sp$cal_spe_func(prok_database = "FAPROTAX")
t_func_sp$cal_spe_func_perc(abundance_weighted = TRUE, dec = 2)
perc_species <- t_func_sp$res_spe_func_perc

# Species-level heartwood means for key metabolisms
sub_header("Species-level heartwood FAPROTAX means")
panel_b_funcs <- c("fermentation", "dark_hydrogen_oxidation", "methanogenesis",
                    "methylotrophy", "methanotrophy", "methanol_oxidation",
                    "aerobic_chemoheterotrophy")
for (func in panel_b_funcs) {
  if (func %in% colnames(perc_species)) {
    vals <- perc_species[, func]
    cat(sprintf("  %s: mean=%.2f%%, range=%.2f-%.2f%%\n",
                func, mean(vals), min(vals), max(vals)))
  }
}

# A. saccharum heartwood peaks
if ("ACSA" %in% rownames(perc_species)) {
  sub_header("A. saccharum heartwood FAPROTAX peaks")
  for (func in c("fermentation", "dark_hydrogen_oxidation", "methylotrophy",
                  "methanogenesis", "aerobic_chemoheterotrophy")) {
    if (func %in% colnames(perc_species)) {
      cat(sprintf("  %s: %.2f%%\n", func, perc_species["ACSA", func]))
    }
  }
}

# All FAPROTAX functions >1% mean in heartwood
sub_header("All FAPROTAX functions >1% mean in heartwood")
func_means <- colMeans(perc_species)
func_means <- sort(func_means[func_means > 1], decreasing = TRUE)
for (fn in names(func_means)) {
  cat(sprintf("  %s: %.2f%%\n", fn, func_means[fn]))
}

# PICRUSt pathway verification
sub_header("PICRUSt pathway associations (VERIFY pathway names)")
picrust_mcra <- "outputs/data/picrust/pathway_associations_mcra_all.csv"
{
  pvals_mcra <- read.csv(picrust_mcra)
  sig_mcra <- pvals_mcra %>% filter(FDR < 0.01) %>% arrange(FDR)
  cat(sprintf("  Significant mcrA-associated pathways (FDR < 0.01): %d\n", nrow(sig_mcra)))
  cat("  Top 15 pathways:\n")
  for (i in 1:min(15, nrow(sig_mcra))) {
    direction <- ifelse(sig_mcra$t[i] > 0, "+", "-")
    cat(sprintf("    %s %s (FDR=%.4f, t=%.2f)\n",
                direction, sig_mcra$pathway[i], sig_mcra$FDR[i], sig_mcra$t[i]))
  }
}

picrust_pmoa <- "outputs/data/picrust/pathway_associations_pmoa.csv"
{
  pvals_pmoa <- read.csv(picrust_pmoa)
  sig_pmoa <- pvals_pmoa %>% filter(FDR < 0.01) %>% arrange(FDR)
  cat(sprintf("\n  Significant pmoA-associated pathways (FDR < 0.01): %d\n", nrow(sig_pmoa)))
  cat("  Top 10 pathways:\n")
  for (i in 1:min(10, nrow(sig_pmoa))) {
    direction <- ifelse(sig_pmoa$t[i] > 0, "+", "-")
    cat(sprintf("    %s %s (FDR=%.4f, t=%.2f)\n",
                direction, sig_pmoa$pathway[i], sig_pmoa$FDR[i], sig_pmoa$t[i]))
  }
}


# ==============================================================================
# SECTION 6b: PICRUSt PATHWAY-GENE HEATMAP (Figure 6)
# ==============================================================================

section_header("SECTION 6b: PICRUSt PATHWAY-GENE HEATMAP (Figure 6)")

# Figure 6 is the mcrA-pathway heatmap from figS11-13_picrust-heatmaps.R
# Uses the no-mcrA-ASV pipeline: pathway abundances reconstructed after removing
# methanogen ASV contributions, then tested against mcrA via LMER.
# FDR threshold 0.001 with mcrA-ASV contribution < 10%.

picrust_no_mcra_file <- "outputs/data/picrust/pathway_associations_mcra_no_mcra_otus.csv"
picrust_combined_file <- "outputs/data/picrust/pathway_associations_combined.csv"

{
  pvals_no_mcra <- read.csv(picrust_no_mcra_file, stringsAsFactors = FALSE)
  stat("Total pathways tested (no-mcrA-ASV pipeline)", nrow(pvals_no_mcra))

  sig_001 <- pvals_no_mcra %>% filter(!is.na(FDR), FDR < 0.001)
  stat("Significant at FDR < 0.001", nrow(sig_001))

  sig_01 <- pvals_no_mcra %>% filter(!is.na(FDR), FDR < 0.01)
  stat("Significant at FDR < 0.01", nrow(sig_01))

  # Apply mcrA-ASV contribution filter (< 10%)
  {
    combined_pvals <- read.csv(picrust_combined_file, stringsAsFactors = FALSE)
    contrib_lookup <- setNames(combined_pvals$mean_percent_from_mcra, combined_pvals$pathway)
    sig_001$gene_contrib <- contrib_lookup[sig_001$pathway]
    stopifnot(!anyNA(sig_001$gene_contrib))   # every tested pathway has a contribution
    sig_filtered <- sig_001 %>% filter(gene_contrib < 0.10)
    stat("After mcrA-ASV contribution < 10% filter", nrow(sig_filtered))
    record("fig6_n_pathways_heatmap", nrow(sig_filtered))

    # List top pathways (positive and negative associations)
    sub_header("Top positively associated pathways (Fig 6)")
    pos_paths <- sig_filtered %>% filter(t > 0) %>% arrange(FDR)
    for (i in 1:min(10, nrow(pos_paths))) {
      desc <- ifelse(!is.na(pos_paths$description[i]) & pos_paths$description[i] != "",
                     pos_paths$description[i], pos_paths$pathway[i])
      cat(sprintf("    %s (FDR=%.2e, t=%.2f)\n", desc, pos_paths$FDR[i], pos_paths$t[i]))
    }

    sub_header("Top negatively associated pathways (Fig 6)")
    neg_paths <- sig_filtered %>% filter(t < 0) %>% arrange(FDR)
    for (i in 1:min(10, nrow(neg_paths))) {
      desc <- ifelse(!is.na(neg_paths$description[i]) & neg_paths$description[i] != "",
                     neg_paths$description[i], neg_paths$pathway[i])
      cat(sprintf("    %s (FDR=%.2e, t=%.2f)\n", desc, neg_paths$FDR[i], neg_paths$t[i]))
    }
  }

  # Verify pathway names mentioned in manuscript
  sub_header("VERIFY: Manuscript-cited pathway names in Fig 6 data")
  manuscript_pathways <- c("reductive acetyl-CoA", "methanogenesis", "homolactic fermentation",
                           "mannan degradation", "L-arginine biosynthesis",
                           "purine nucleotide degradation", "heme biosynthesis",
                           "peptidoglycan", "teichoic acid")
  for (mp in manuscript_pathways) {
    matches <- pvals_no_mcra %>%
      filter(grepl(mp, pathway, ignore.case = TRUE) |
             grepl(mp, description, ignore.case = TRUE))
    if (nrow(matches) > 0) {
      best <- matches %>% arrange(FDR) %>% slice(1)
      cat(sprintf("  FOUND: '%s' -> %s (FDR=%.4f, t=%.2f)\n",
                  mp, best$pathway, best$FDR, best$t))
    } else {
      cat(sprintf("  NOT FOUND: '%s'\n", mp))
    }
  }
}


# ==============================================================================
# SECTION 7: FELLED BLACK OAK PROFILES (Figure 7)
# ==============================================================================

section_header("SECTION 7: FELLED BLACK OAK PROFILES (Figure 7)")

# Felled-oak data from the canonical tables: trunk gas (internal_gas.csv, calibrated in
# code/tables/32_internal_gas.R), mcrA (ddpcr.csv + samples.csv), flux (flux_stem.csv)
# Three panels: internal CH4, mcrA abundance, CH4 flux along stem height
{
  int_gas_bo <- read_table("internal_gas") %>% filter(campaign == "2022_felled_oak") %>%
    transmute(Tree.Height = height_m, CH4_conc = CH4_ppm, O2_conc = O2_ppm)

  sub_header("Internal CH4 concentrations")
  stat("Number of internal gas measurements", nrow(int_gas_bo))
  stat("Height range", paste(min(int_gas_bo$Tree.Height, na.rm = TRUE), "to",
                              max(int_gas_bo$Tree.Height, na.rm = TRUE)), "m")
  stat("Max internal CH4", round(max(int_gas_bo$CH4_conc, na.rm = TRUE), 0), "ppm")
  stat("Mean internal CH4", round(mean(int_gas_bo$CH4_conc, na.rm = TRUE), 0), "ppm")
  stat("Median internal CH4", round(median(int_gas_bo$CH4_conc, na.rm = TRUE), 0), "ppm")
  record("black_oak_max_ch4_ppm", max(int_gas_bo$CH4_conc, na.rm = TRUE))

  # CH4 at mid-stem heights (4-6 m)
  mid_gas <- int_gas_bo %>% filter(Tree.Height >= 4 & Tree.Height <= 6)
  if (nrow(mid_gas) > 0) {
    stat("CH4 at 4-6 m (mean)", round(mean(mid_gas$CH4_conc, na.rm = TRUE), 0), "ppm")
    stat("CH4 at 4-6 m (max)", round(max(mid_gas$CH4_conc, na.rm = TRUE), 0), "ppm")
  }

  # O2 range
  o2_pct <- int_gas_bo$O2_conc / 1e6 * 100
  stat("O2 range", paste(round(min(o2_pct, na.rm = TRUE), 1), "to",
                          round(max(o2_pct, na.rm = TRUE), 1)), "%")
  stat("O2 mean", round(mean(o2_pct, na.rm = TRUE), 1), "%")
  record("black_oak_o2_range_pct", paste(round(min(o2_pct, na.rm = TRUE), 1),
                                          round(max(o2_pct, na.rm = TRUE), 1)))

  # --- mcrA from black oak cores ---
  sub_header("Black oak mcrA abundance")
  .bo_smp <- read_table("samples") %>% filter(campaign == "2022_felled_oak")
  mcra_wood <- read_table("ddpcr") %>% filter(gene == "mcrA", qc_pass, sample_id %in% .bo_smp$sample_id) %>%
    transmute(sample_id, Conc = copies_per_uL_rxn,
              Height_cm = .bo_smp$height_cm[match(sample_id, .bo_smp$sample_id)],
              Component = c(heartwood = "Heartwood", sapwood = "Sapwood")[.bo_smp$position[match(sample_id, .bo_smp$sample_id)]]) %>%
    filter(Component %in% c("Heartwood", "Sapwood"))

  stat("Total mcrA measurements (wood)", nrow(mcra_wood))
  stat("Heartwood samples", sum(mcra_wood$Component == "Heartwood"))
  stat("Sapwood samples", sum(mcra_wood$Component == "Sapwood"))

  for (comp in c("Heartwood", "Sapwood")) {
    d <- mcra_wood %>% filter(Component == comp)
    pos <- d$Conc[d$Conc > 0]
    cat(sprintf("  %s: detection = %.0f%%, among positive: mean = %.1f, max = %.1f copies/uL\n",
                comp, 100 * length(pos) / nrow(d),
                ifelse(length(pos) > 0, mean(pos), NA),
                ifelse(length(pos) > 0, max(pos), NA)))
    if (length(pos) > 0) {
      cat(sprintf("    Range (log10): %.1f to %.1f (%.0f orders of magnitude)\n",
                  min(log10(pos)), max(log10(pos)),
                  max(log10(pos)) - min(log10(pos))))
    }
  }

  # Overall range in positive wood samples
  all_pos <- mcra_wood$Conc[mcra_wood$Conc > 0]
  if (length(all_pos) > 1) {
    stat("mcrA range across all wood (orders of magnitude)",
         round(max(log10(all_pos)) - min(log10(all_pos)), 1))
  }
  record("black_oak_mcra_orders_magnitude",
         round(max(log10(all_pos)) - min(log10(all_pos)), 1))

  # Detection at mid-stem (4-6 m)
  mid_mcra <- mcra_wood %>% filter(Height_cm >= 400 & Height_cm <= 600)
  if (nrow(mid_mcra) > 0) {
    det_mid <- round(100 * sum(mid_mcra$Conc > 0) / nrow(mid_mcra), 0)
    stat("mcrA detection at 4-6 m height", det_mid, "%")
    record("black_oak_mcra_detection_mid_pct", det_mid)
  }

  # --- Flux data ---
  sub_header("Black oak CH4 flux")
  flux_bo <- flux_stem_campaign("2022_felled_oak") %>%
    transmute(Height_m = height_cm / 100, CH4_best.flux = CH4_flux_nmol_m2_s) %>%
    filter(!is.na(Height_m))                      # the basal seam has no height

  stat("Total flux measurements", nrow(flux_bo))
  stat("Mean flux", round(mean(flux_bo$CH4_best.flux, na.rm = TRUE), 4), "nmol m-2 s-1")
  stat("Flux range", paste(round(min(flux_bo$CH4_best.flux, na.rm = TRUE), 4), "to",
                            round(max(flux_bo$CH4_best.flux, na.rm = TRUE), 4)), "nmol m-2 s-1")

  # Flux by height zone
  flux_mid <- flux_bo %>% filter(Height_m >= 4 & Height_m <= 6)
  flux_base <- flux_bo %>% filter(Height_m < 2)
  flux_top <- flux_bo %>% filter(Height_m > 6)

  if (nrow(flux_mid) > 0) {
    stat("Flux at 4-6 m (mean)", round(mean(flux_mid$CH4_best.flux), 4), "nmol m-2 s-1")
    stat("Flux at 4-6 m (max)", round(max(flux_mid$CH4_best.flux), 4), "nmol m-2 s-1")
  }
  if (nrow(flux_base) > 0) {
    stat("Flux at 0-2 m (mean)", round(mean(flux_base$CH4_best.flux), 4), "nmol m-2 s-1")
  }
  if (nrow(flux_top) > 0) {
    stat("Flux at >6 m (mean)", round(mean(flux_top$CH4_best.flux), 4), "nmol m-2 s-1")
  }

  # Merge flux with internal CH4 to show relationship
  flux_gas_bo <- flux_bo %>%
    left_join(
      int_gas_bo %>% group_by(Tree.Height) %>%
        summarize(ch4_internal = mean(CH4_conc, na.rm = TRUE), .groups = "drop"),
      by = c("Height_m" = "Tree.Height")
    )

  if (sum(!is.na(flux_gas_bo$ch4_internal)) >= 3) {
    cor_flux_ch4 <- cor.test(flux_gas_bo$CH4_best.flux, flux_gas_bo$ch4_internal,
                              use = "complete.obs")
    sub_header("Flux vs internal CH4 correlation (black oak)")
    stat("Pearson r", round(cor_flux_ch4$estimate, 3))
    stat("p-value", round(cor_flux_ch4$p.value, 4))
    stat("n measurements", sum(!is.na(flux_gas_bo$ch4_internal) & !is.na(flux_gas_bo$CH4_best.flux)))
  }

}


# ==============================================================================
# SECTION 8: INDIVIDUAL GENE-FLUX (Figure S7, S10)
# ==============================================================================

section_header("SECTION 8: INDIVIDUAL GENE-FLUX RELATIONSHIPS (Figs S7, S10)")

# Area-weighted genes per tree and species-level flux: code/lib/species_gene_data.R,
# the prep every gene-flux analysis uses (helper_scale_dependent_gene_patterns.R, 17).
# Sourced into its own environment so this script's ymf2021 (with soil, above) stays.
SG <- new.env()
sys.source("code/lib/species_gene_data.R", envir = SG)
tree_level <- SG$tree_level_complete %>%
  rename(log_mcra = log_tree_mcra, log_mmox = log_tree_mmox, ratio = ratio_mcra_methanotroph) %>%
  relocate(log_pmoa, .before = log_mmox)
flux_by_sp <- SG$flux_by_species %>% dplyr::select(species, species_id, n_flux, median_flux)

stat("Trees with complete genes AND flux", nrow(tree_level))

sub_header("Individual tree-level R2 (species-controlled lm)")
genes <- list(
  mcrA = "log_mcra", pmoA = "log_pmoa", mmoX = "log_mmox",
  methanotroph = "log_methanotroph", ratio = "log_ratio"
)
for (nm in names(genes)) {
  fml <- as.formula(paste("CH4_flux ~ species +", genes[[nm]]))
  m <- lm(fml, data = tree_level)
  r2 <- summary(m)$r.squared
  p_gene <- summary(m)$coefficients[genes[[nm]], "Pr(>|t|)"]
  cat(sprintf("  %s: R2 = %.4f, p(gene) = %.4f\n", nm, r2, p_gene))
}

# --- Concentration-based models (S1 methods: heartwood/sapwood separately) ---
sub_header("Concentration-based LMs (tissue-specific, Figure S7)")
conc_cols <- c("ddpcr_mcra_probe_Inner_loose", "ddpcr_mcra_probe_Outer_loose",
               "ddpcr_pmoa_Inner_loose", "ddpcr_pmoa_Outer_loose",
               "ddpcr_mmox_Inner_loose", "ddpcr_mmox_Outer_loose")
conc_available <- all(conc_cols %in% names(ymf2021))

if (conc_available) {
  conc_data_s7 <- ymf2021 %>%
    filter(!is.na(CH4_best.flux_125cm), !is.na(species_id)) %>%
    mutate(
      species = SPECIES_NAMES[species_id],
      log_hw_mcra = log10(ddpcr_mcra_probe_Inner_loose + 1),
      log_sw_mcra = log10(ddpcr_mcra_probe_Outer_loose + 1),
      log_hw_pmoa = log10(ddpcr_pmoa_Inner_loose + 1),
      log_sw_pmoa = log10(ddpcr_pmoa_Outer_loose + 1),
      log_hw_mmox = log10(ddpcr_mmox_Inner_loose + 1),
      log_sw_mmox = log10(ddpcr_mmox_Outer_loose + 1)
    ) %>%
    filter(!is.na(species))

  # Full 6-gene concentration model
  conc_complete <- conc_data_s7 %>%
    filter(complete.cases(dplyr::select(., starts_with("log_hw_"), starts_with("log_sw_"))))
  stat("Trees with complete concentration data", nrow(conc_complete))

  if (nrow(conc_complete) > 20) {
    m_conc_full <- lm(CH4_best.flux_125cm ~ species + log_hw_mcra + log_sw_mcra +
                        log_hw_pmoa + log_sw_pmoa + log_hw_mmox + log_sw_mmox,
                      data = conc_complete)
    s <- summary(m_conc_full)
    stat("Full concentration model R2", round(s$r.squared, 3))
    stat("Full concentration model adj-R2", round(s$adj.r.squared, 3))
    stat("Full concentration model AIC", round(AIC(m_conc_full), 1))
    record("conc_full_r2", s$r.squared)

    # Extract sapwood mmoX coefficient
    if ("log_sw_mmox" %in% rownames(s$coefficients)) {
      coef_sw_mmox <- s$coefficients["log_sw_mmox", ]
      cat(sprintf("  Sapwood mmoX: beta = %.3f, SE = %.3f, p = %.4f\n",
                  coef_sw_mmox[1], coef_sw_mmox[2], coef_sw_mmox[4]))
    }

    # Type II ANOVA
    if (requireNamespace("car", quietly = TRUE)) {
      a2 <- car::Anova(m_conc_full, type = 2)
      cat("  Type II ANOVA p-values:\n")
      for (term in rownames(a2)) {
        if (term != "Residuals") {
          cat(sprintf("    %s: F = %.3f, p = %.4f\n", term, a2[term, "F value"], a2[term, "Pr(>F)"]))
        }
      }
    }
  }

  # Best area-weighted model: species + mmoX
  m_sp_mmox <- lm(CH4_flux ~ species + log_mmox,
                   data = tree_level %>% filter(!is.na(CH4_flux)))
  s_sp_mmox <- summary(m_sp_mmox)
  sub_header("Best area-weighted model: species + mmoX")
  stat("R2", round(s_sp_mmox$r.squared, 3))
  stat("Adj R2", round(s_sp_mmox$adj.r.squared, 3))
  stat("AIC", round(AIC(m_sp_mmox), 1))
  if ("log_mmox" %in% rownames(s_sp_mmox$coefficients)) {
    cat(sprintf("  mmoX: beta = %.3f, SE = %.3f, p = %.4f\n",
                s_sp_mmox$coefficients["log_mmox", 1],
                s_sp_mmox$coefficients["log_mmox", 2],
                s_sp_mmox$coefficients["log_mmox", 4]))
  }
  record("aw_sp_mmox_r2", s_sp_mmox$r.squared)
}

# --- Concentration-based species-level models (S8) ---
sub_header("Concentration-based species-level models (Figure S8)")
if (conc_available) {
  # Species medians for heartwood mcrA
  sp_hw_mcra <- ymf2021 %>%
    filter(!is.na(species_id), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
    group_by(species_id) %>%
    summarise(n_trees = n(),
              median_hw_mcra = median(ddpcr_mcra_probe_Inner_loose, na.rm = TRUE),
              .groups = "drop") %>%
    filter(n_trees >= 5)

  # Merge with species-level flux
  sp_hw_flux <- sp_hw_mcra %>%
    inner_join(flux_by_sp %>% filter(n_flux >= 5), by = c("species_id"))

  if (nrow(sp_hw_flux) >= 5) {
    cor_hw <- cor.test(log10(sp_hw_flux$median_hw_mcra + 1), sp_hw_flux$median_flux)
    cat(sprintf("  Heartwood mcrA: R2 = %.3f, r = %.3f, p = %.4f (n=%d species)\n",
                cor_hw$estimate^2, cor_hw$estimate, cor_hw$p.value, nrow(sp_hw_flux)))
    record("species_hw_mcra_flux_r2", cor_hw$estimate^2)
  }
}


# ==============================================================================
# SECTION 9: SPECIES-LEVEL GENE-FLUX (Figure 8, S8, S9, S10 right, S13)
# ==============================================================================

section_header("SECTION 9: SPECIES-LEVEL GENE-FLUX (Fig 8, S8, S9, S10, S13)")

# Flux by species (combined 2021 + 2023): flux_by_sp, from code/lib/species_gene_data.R (Section 8)

# Species-level gene summaries
sp_mcra <- tree_level %>%
  group_by(species, species_id) %>%
  summarise(n_trees = n(), median_mcra = median(mcrA), .groups = "drop") %>%
  inner_join(flux_by_sp, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

sp_methanotroph <- tree_level %>%
  group_by(species, species_id) %>%
  summarise(n_trees = n(), median_mt = median(methanotroph_total), .groups = "drop") %>%
  inner_join(flux_by_sp, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

sp_ratio <- tree_level %>%
  group_by(species, species_id) %>%
  summarise(n_trees = n(), median_ratio = median(ratio),
            median_log_ratio = median(log_ratio), .groups = "drop") %>%
  inner_join(flux_by_sp, by = c("species", "species_id")) %>%
  filter(n_trees >= 5, n_flux >= 5)

# Individual tree area-weighted mcrA vs flux (no species control)
ind_cor <- cor.test(log10(tree_level$mcrA + 1), tree_level$CH4_flux)
sub_header("Individual-level mcrA vs flux (no species control)")
stat("R2", ind_cor$estimate^2)
stat("n trees", nrow(tree_level))

# Species-level correlations
sub_header("Species-level correlations (Pearson)")
cor_mcra <- cor.test(log10(sp_mcra$median_mcra + 1), sp_mcra$median_flux)
cat(sprintf("  mcrA: R2 = %.3f, r = %.3f, p = %.4f (n=%d species)\n",
            cor_mcra$estimate^2, cor_mcra$estimate, cor_mcra$p.value, nrow(sp_mcra)))
record("species_mcra_flux_r2", cor_mcra$estimate^2)
record("species_mcra_flux_p", cor_mcra$p.value)

cor_mt <- cor.test(log10(sp_methanotroph$median_mt + 1), sp_methanotroph$median_flux)
cat(sprintf("  methanotroph: R2 = %.3f, r = %.3f, p = %.4f (n=%d species)\n",
            cor_mt$estimate^2, cor_mt$estimate, cor_mt$p.value, nrow(sp_methanotroph)))

cor_ratio <- cor.test(sp_ratio$median_log_ratio, sp_ratio$median_flux)
cat(sprintf("  ratio: R2 = %.3f, r = %.3f, p = %.4f (n=%d species)\n",
            cor_ratio$estimate^2, cor_ratio$estimate, cor_ratio$p.value, nrow(sp_ratio)))
record("species_ratio_flux_r2", cor_ratio$estimate^2)

# AIC comparison
m_mcra <- lm(median_flux ~ log10(median_mcra + 1), data = sp_mcra)
m_ratio <- lm(median_flux ~ median_log_ratio, data = sp_ratio)
sub_header("AIC comparison")
stat("mcrA model AIC", AIC(m_mcra))
stat("Ratio model AIC", AIC(m_ratio))

# mcrA vs methanotroph correlation (S13)
sub_header("mcrA vs methanotroph (Figure S13)")
tree_mcra_mt <- cor.test(tree_level$log_mcra, tree_level$log_methanotroph)
cat(sprintf("  Tree-level: r = %.3f, p = %.4f (n=%d)\n",
            tree_mcra_mt$estimate, tree_mcra_mt$p.value, nrow(tree_level)))

sp_mcra_mt <- inner_join(
  tree_level %>% group_by(species) %>%
    summarise(n_trees = n(), med_mcra = median(mcrA), med_mt = median(methanotroph_total), .groups = "drop"),
  flux_by_sp %>% dplyr::select(species, n_flux), by = "species"
) %>% filter(n_trees >= 5, n_flux >= 5)   # the ten gene-flux species (as Fig 8 and Fig S24b)

if (nrow(sp_mcra_mt) >= 3) {
  sp_cor <- cor.test(log10(sp_mcra_mt$med_mcra + 1), log10(sp_mcra_mt$med_mt + 1))
  cat(sprintf("  Species-level: r = %.3f, p = %.4f (n=%d)\n",
              sp_cor$estimate, sp_cor$p.value, nrow(sp_mcra_mt)))
}


# ==============================================================================
# SECTION 10: RANDOM FOREST & UPSCALING (Figure 9, S9)
# ==============================================================================

section_header("SECTION 10: RANDOM FOREST & UPSCALING (Fig 9, S9)")

# OOB R²
sub_header("Random Forest OOB performance")
if (exists("TreeRF")) {
  tree_r2 <- TreeRF$r.squared
  stat("Tree stem RF OOB R2", tree_r2)
  record("rf_tree_oob_r2", tree_r2)
}
if (exists("SoilRF")) {
  soil_r2 <- SoilRF$r.squared
  stat("Soil RF OOB R2", soil_r2)
  record("rf_soil_oob_r2", soil_r2)
}

# Monthly fluxes
sub_header("Monthly predicted fluxes")
# Phi_*_umol_m2_s already holds nmol m-2 s-1 -- the name is a documented misnomer.
# The x1000 that stood here reported the tree mean as 4.96 nmol m-2 s-1 instead of
# 0.00496 and soil as -1078 instead of -1.08.
# canonical_monthly.csv (code/06_upscale) is the source of the monthly numbers in
# Results §10 and Fig 9c; the legacy MONTHLY_FLUXES.csv used here before disagreed.
canon_monthly <- read.csv("outputs/data/canonical_monthly.csv")
tree_monthly_nmol <- canon_monthly$tree_nmol_m2_s
soil_monthly_nmol <- canon_monthly$soil_nmol_m2_s

stat("Tree stem annual mean", mean(tree_monthly_nmol), "nmol m-2 s-1")
stat("Tree stem monthly range", paste(round(min(tree_monthly_nmol), 4), "to",
                                       round(max(tree_monthly_nmol), 4)), "nmol m-2 s-1")
stat("Soil annual mean", mean(soil_monthly_nmol), "nmol m-2 s-1")
stat("Soil monthly range", paste(round(min(soil_monthly_nmol), 4), "to",
                                  round(max(soil_monthly_nmol), 4)), "nmol m-2 s-1")
stat("Per-unit-area offset", paste0(round(abs(mean(tree_monthly_nmol) / mean(soil_monthly_nmol)) * 100, 1), "%"))

# Annual totals
sub_header("Annual totals")
stat("Tree emissions", round(canon("tree_measured_mg_m2_yr"), 3), "mg CH4 m-2 yr-1")
stat("Soil uptake", round(canon("soil_mg_m2_yr"), 3), "mg CH4 m-2 yr-1")
stat("Net budget", round(canon("net_measured_mg_m2_yr"), 3), "mg CH4 m-2 yr-1")
record("annual_tree_mg_m2_yr", canon("tree_measured_mg_m2_yr"))
record("annual_soil_mg_m2_yr", canon("soil_mg_m2_yr"))
record("net_measured_mg_m2_yr", canon("net_measured_mg_m2_yr"))
cat(sprintf("  [legacy ANNUAL_SUMMARY.csv, NOT used: tree %.3f, soil %.3f -- see note at the top]\n",
            annual_summary$annual_tree_mg_m2, annual_summary$annual_soil_mg_m2))

# Surface area and offset
sub_header("Stem surface area and offset")
# Load inventory for geometry
# Read, not recomputed. This block used PLACEHOLDER_INVENTORY (7,010 stems, no bytag
# survey), a hardcoded plot_area_ha <- 10.2 -- the lat/lon BOUNDING BOX of a plot
# rotated ~10 degrees from north, 2.7x the censused stand -- and pi*DBH*2 for every
# stem including the 1,904 Kalmia whose band is 0.75 m. It was also inside a tryCatch
# whose guard ("DBH" %in% names(inv)) is FALSE for that object, so it silently emitted
# nothing at all.
stat("Total stem surface area (0-2m)", round(canon("stem_area_0_2m_m2"), 0), "m2")
stat("Stem area as % of censused stand", round(100*canon("stem_area_index"), 2), "%")
stat("Censused stand area", round(canon("stand_area_m2")/1e4, 2), "ha")
stat("Stems in stand", canon("n_stems"))

tree_ann <- canon("tree_measured_mg_m2_yr")
soil_ann <- canon("soil_mg_m2_yr")
offset_pct <- abs(tree_ann / soil_ann) * 100
stat("Plot-level offset", round(offset_pct, 2), "%")
record("plot_offset_pct", offset_pct)

# WAI extrapolation
# WHOLE-SURFACE SCALING: read the grid, do not re-derive it here.
# What stood here multiplied a PER-M2-OF-GROUND tree flux by a woody area index, which
# is per m2 of woody surface per m2 of ground -- the two are not composable that way,
# so that calculation was removed. The scaling is
# now a 240-combination grid (6 flux forms x 5 WAI x 4 branch x 2 bole) produced by
# 07_scaling_full_grid.R, which holds the measured 0-2 m band fixed and lets only the
# area above 2 m absorb the WAI choice. Quote the RANGE; the named pair (flux declines /
# holds above 2 m) is two cells of 240.
local({
  f <- "outputs/data/scaling_full_grid.csv"; h <- "outputs/data/scaling_headline.csv"
  G <- read.csv(f, stringsAsFactors = FALSE); H <- read.csv(h, stringsAsFactors = FALSE)
  stat("Whole-surface range, all 240 combinations",
       sprintf("%.1f to %.1f", min(G$total_mg), max(G$total_mg)), "mg CH4 m-2 yr-1")
  stat("As % of soil uptake",
       sprintf("%.1f to %.1f", min(G$pct_of_soil), max(G$pct_of_soil)), "%")
  stat("Net budget across the grid",
       sprintf("%.1f to %.1f", min(G$net_mg), max(G$net_mg)), "mg CH4 m-2 yr-1")
  stat("Sink in", sprintf("%d of %d combinations", sum(G$net_mg < 0), nrow(G)))
  Gp <- G[G$flux != "linear_bounded_median", ]
  stat(sprintf("Range excluding the uptake form (%d)", nrow(Gp)),
       sprintf("%.1f to %.1f", min(Gp$total_mg), max(Gp$total_mg)), "mg CH4 m-2 yr-1")
  for (k in seq_len(nrow(H))) {
    stat(sprintf("Whole surface, flux %s above 2 m (%s)",
                 ifelse(H$case[k] == "decline", "declines", "holds"), H$flux[k]),
         sprintf("%.2f mg (%.1f%% of soil, %.0f%% extrapolated; net %.1f)",
                 H$total_mg[k], H$pct_of_soil[k], H$pct_extrapolated[k], H$net_mg[k]))
    record(sprintf("whole_surface_%s_mg_m2_yr", H$case[k]), H$total_mg[k])
  }
  record("whole_surface_lo_mg_m2_yr", min(G$total_mg))
  record("whole_surface_hi_mg_m2_yr", max(G$total_mg))
})


# --- Scaling sensitivities (06_upscale/10 and 11) ---
sub_header("Sensitivity: forests with air temperature (05_model/15_rf_air_temperature_test.R)")
SA <- read.csv("outputs/data/rf_air_temperature_test.csv", stringsAsFactors = FALSE)
for (i in which(SA$section == "skill"))
  cat(sprintf("  %-13s OOB R2 %.3f | grouped CV R2 %.3f (sd %.3f)\n", SA$spec[i], SA$oob_r2[i], SA$grouped_r2[i], SA$grouped_r2_sd[i]))
for (i in which(SA$section == "budget"))
  cat(sprintf("  %-28s %-16s soil %8.2f | tree band %.3f mg m-2 yr-1 | tree %.2f%% of soil\n", SA$spec[i], SA$drivers[i],
              SA$soil_mg_m2_yr[i], SA$tree_band_mg_m2_yr[i], SA$tree_pct_of_soil[i]))
sub_header("Scaling sensitivity: daily vs monthly drivers (06_upscale/11_sens_daily_drivers.R)")
SY <- read.csv("outputs/data/sens_daily_drivers.csv", stringsAsFactors = FALSE)
for (i in seq_len(nrow(SY)))
  cat(sprintf("  %-46s monthly %.4f | daily %.4f %s | daily/monthly %.3f (%d days)\n", SY$term[i],
              SY$monthly_drivers[i], SY$daily_drivers[i], SY$units[i], SY$ratio_daily_to_monthly[i], SY$n_days[i]))
record("sens_daily_soil_ratio", SY$ratio_daily_to_monthly[1]); record("sens_daily_tree_ratio", SY$ratio_daily_to_monthly[2])

# ==============================================================================
# SECTION 11: ADDITIONAL METHANOTROPH RESULTS (S10, S12)


# ==============================================================================

section_header("SECTION 11: ADDITIONAL METHANOTROPH RESULTS (S10, S12)")
cat("  (Statistics from existing analyses not yet in manuscript)\n\n")

# pmoA vs mmoX correlation
sub_header("pmoA vs mmoX correlation (wood samples)")
both_mt <- tree_level %>% filter(!is.na(pmoA), !is.na(mmoX), pmoA > 0, mmoX > 0)
if (nrow(both_mt) > 3) {
  cor_pm <- cor.test(log10(both_mt$pmoA), log10(both_mt$mmoX))
  cat(sprintf("  Pearson r = %.3f, R2 = %.3f, p = %.4f (n=%d trees)\n",
              cor_pm$estimate, cor_pm$estimate^2, cor_pm$p.value, nrow(both_mt)))
  m_pm <- lm(log10(mmoX) ~ log10(pmoA), data = both_mt)
  cat(sprintf("  Slope = %.3f\n", coef(m_pm)[2]))
}

# pmoA:mmoX ratio statistics
sub_header("pmoA:mmoX ratio by compartment")
for (loc in c("Inner", "Outer")) {
  pmoa_col <- paste0("ddpcr_pmoa_", loc, "_loose")
  mmox_col <- paste0("ddpcr_mmox_", loc, "_loose")
  if (pmoa_col %in% names(ymf2021) && mmox_col %in% names(ymf2021)) {
    p <- ymf2021[[pmoa_col]]
    m <- ymf2021[[mmox_col]]
    both_pos <- !is.na(p) & !is.na(m) & p > 0 & m > 0
    if (sum(both_pos) > 0) {
      ratio_vals <- log10(p[both_pos] / m[both_pos])
      label <- c(Inner = "Heartwood", Outer = "Sapwood")[loc]
      cat(sprintf("  %s: mean log10(pmoA/mmoX) = %.2f, median = %.2f, SD = %.2f (n=%d)\n",
                  label, mean(ratio_vals), median(ratio_vals), sd(ratio_vals), sum(both_pos)))
    }
  }
}

# mmoX as best individual predictor
sub_header("Model comparison: which gene best predicts flux?")
m_sp_only <- lm(CH4_flux ~ species, data = tree_level)
r2_sp <- summary(m_sp_only)$r.squared
for (g in c("log_mcra", "log_pmoa", "log_mmox", "log_methanotroph", "log_ratio")) {
  fml <- as.formula(paste("CH4_flux ~ species +", g))
  m <- lm(fml, data = tree_level)
  delta_r2 <- summary(m)$r.squared - r2_sp
  p_gene <- summary(m)$coefficients[g, "Pr(>|t|)"]
  sig <- ifelse(p_gene < 0.05, "*", "")
  cat(sprintf("  Species + %s: delta-R2 = %.4f, p = %.4f %s\n", g, delta_r2, p_gene, sig))
}

# Soil genes adding no value. Complete cases: trees without the soil gene are left out of
# both models (until 2026-10-16 a missing soil gene entered as 0 copies).
sub_header("Soil genes beyond tree mmoX")
if ("log_mmox" %in% names(tree_level)) {
  for (soil_gene in c("ddpcr_mcra_probe_Mineral_loose", "ddpcr_pmoa_Mineral_loose", "ddpcr_mmox_Mineral_loose")) {
    if (soil_gene %in% names(ymf2021)) {
      tree_w_soil <- tree_level %>%
        left_join(ymf2021 %>% dplyr::select(tree_id, !!soil_gene), by = "tree_id") %>%
        filter(!is.na(.data[[soil_gene]]))
      sg <- tree_w_soil[[soil_gene]]
      tree_w_soil$soil_gene_log <- log10(ifelse(sg <= 0, 1, sg))
      m_base <- lm(CH4_flux ~ species + log_mmox, data = tree_w_soil)
      m_add <- lm(CH4_flux ~ species + log_mmox + soil_gene_log, data = tree_w_soil)
      delta <- summary(m_add)$r.squared - summary(m_base)$r.squared
      p_soil <- summary(m_add)$coefficients["soil_gene_log", "Pr(>|t|)"]
      short_name <- sub("ddpcr_", "", sub("_Mineral_loose", "", soil_gene))
      cat(sprintf("  + soil %s: delta-R2 = %.4f, p = %.4f (n = %d of %d trees)\n", short_name, delta, p_soil,
                  nrow(tree_w_soil), nrow(tree_level)))
    }
  }
}


# ==============================================================================
# SECTION 12: SUPPLEMENTARY FIGURE STATISTICS
# ==============================================================================

section_header("SECTION 12: SUPPLEMENTARY FIGURE STATISTICS")

# --- Methanotroph relative abundance (supporting Figure 5) ---
sub_header("Methanotroph relative abundances by compartment (Figure 5 support)")
# Already have ps.filt, otu_df, tax_df, samp_meta from Section 5
if (exists("ps.filt") && exists("mt_defs")) {

  # Known + putative methanotroph ASVs (reuse from Section 5)
  all_mt_asvs <- unique(c(known_asvs, putative_only))

  for (comp in c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil")) {
    samps <- samp_meta %>% filter(compartment == comp) %>% rownames()
    if (length(samps) > 0 && length(all_mt_asvs) > 0) {
      mt_pcts <- colSums(otu_df[intersect(all_mt_asvs, rownames(otu_df)), samps, drop = FALSE])
      cat(sprintf("  %s: methanotroph mean = %.3f%%, range = %.3f-%.3f%% (n=%d)\n",
                  comp, mean(mt_pcts), min(mt_pcts), max(mt_pcts), length(samps)))
    }
  }

  # Per-species methanotroph abundance (heartwood only)
  sub_header("  Per-species methanotroph abundance (heartwood)")
  hw_samps <- samp_meta %>% filter(compartment == "Heartwood", !is.na(species.x))
  mt_pct_vec <- colSums(otu_df[intersect(all_mt_asvs, rownames(otu_df)),
                                rownames(hw_samps), drop = FALSE])
  hw_samps$mt_pct <- mt_pct_vec[rownames(hw_samps)]
  sp_mt_means <- hw_samps %>%
    group_by(species.x) %>%
    summarise(mean_mt = mean(mt_pct, na.rm = TRUE), n = n(), .groups = "drop") %>%
    arrange(desc(mean_mt))
  for (i in 1:min(5, nrow(sp_mt_means))) {
    cat(sprintf("    %s: %.3f%% (n=%d)\n",
                SPECIES_NAMES[sp_mt_means$species.x[i]],
                sp_mt_means$mean_mt[i], sp_mt_means$n[i]))
  }
} else {
  cat("  [Requires phyloseq objects from Section 5]\n")
}

# --- Figure S3: FAPROTAX (already in Section 6, just cross-reference) ---
sub_header("Figure S3: FAPROTAX (see Section 6 above)")

# --- PICRUSt pathway statistics (supporting Figures S4-S5) ---
sub_header("PICRUSt pathway statistics (Figures S4-S5)")
picrust_mcra_all <- "outputs/data/picrust/pathway_associations_mcra_all.csv"
{
  pa_all <- read.csv(picrust_mcra_all, stringsAsFactors = FALSE)
  sig_all <- pa_all %>% filter(!is.na(FDR), FDR < 0.01)
  stat("Total pathways tested (all ASVs)", nrow(pa_all))
  stat("Significant at FDR < 0.01 (all ASVs)", nrow(sig_all))

  # Top positive (heartwood-enriched) and negative (sapwood-enriched)
  pos_all <- sig_all %>% filter(t > 0) %>% arrange(FDR)
  neg_all <- sig_all %>% filter(t < 0) %>% arrange(FDR)
  cat(sprintf("  Positively associated: %d, Negatively associated: %d\n",
              nrow(pos_all), nrow(neg_all)))
}

sub_header("Figure S5: pmoA pathway associations")
picrust_pmoa_no <- "outputs/data/picrust/pathway_associations_pmoa_no_pmoa_otus.csv"
{
  pa_pmoa <- read.csv(picrust_pmoa_no, stringsAsFactors = FALSE)
  sig_pmoa_no <- pa_pmoa %>% filter(!is.na(FDR), FDR < 0.01)
  stat("Significant pmoA pathways (no-pmoA ASVs, FDR < 0.01)", nrow(sig_pmoa_no))
  pos_pmoa <- sig_pmoa_no %>% filter(t > 0)
  neg_pmoa <- sig_pmoa_no %>% filter(t < 0)
  cat(sprintf("  Positive: %d, Negative: %d\n", nrow(pos_pmoa), nrow(neg_pmoa)))
}

# --- Figure S4: PICRUSt mcrA no-mcrA-ASV heatmap (expanded version of Fig 6) ---
sub_header("Figure S4: PICRUSt mcrA (no-mcrA ASVs, full FDR < 0.01) heatmap")
{
  cat(sprintf("  See pathway counts above. Full FDR<0.01 set: %d pathways\n", nrow(sig_all)))
}

# --- Figure S5: PICRUSt pmoA heatmap ---
sub_header("Figure S5: PICRUSt pmoA heatmap")
cat("  See pmoA pathway counts above.\n")

# --- Figure S7: Internal gas beeswarm by species ---
sub_header("Figure S7: Internal gas concentrations by species")
gas_data <- ymf2021 %>%
  filter(!is.na(CH4_concentration), !is.na(species_id)) %>%
  mutate(species_latin = SPECIES_NAMES[species_id]) %>%
  filter(!is.na(species_latin))

stat("Total trees with internal gas", nrow(gas_data))
stat("Number of species", n_distinct(gas_data$species_latin))

# Overall CH4 stats
stat("Overall CH4 mean", round(mean(gas_data$CH4_concentration), 0), "ppm")
stat("Overall CH4 median", round(median(gas_data$CH4_concentration), 0), "ppm")
stat("Overall CH4 range", paste(round(min(gas_data$CH4_concentration), 0), "to",
                                 round(max(gas_data$CH4_concentration), 0)), "ppm")
# Share of trees above thresholds (Discussion: comparison with Wang et al. 2026, who report
# >= 1,000 ppm in ~5% of upland stems)
for (thr in c(10, 1000)) {
  stat(sprintf("Trees with internal CH4 >= %g ppm", thr),
       sprintf("%d of %d (%.0f%%)", sum(gas_data$CH4_concentration >= thr), nrow(gas_data),
               100 * mean(gas_data$CH4_concentration >= thr)), "")
  record(sprintf("internal_ch4_pct_ge_%g", thr), 100 * mean(gas_data$CH4_concentration >= thr))
}

# Species with highest/lowest
sp_ch4 <- gas_data %>%
  group_by(species_latin) %>%
  summarise(mean_ch4 = mean(CH4_concentration), n = n(), .groups = "drop") %>%
  arrange(desc(mean_ch4))
cat("  Top 3 species by mean CH4:\n")
for (i in 1:min(3, nrow(sp_ch4))) {
  cat(sprintf("    %s: %.0f ppm (n=%d)\n", sp_ch4$species_latin[i], sp_ch4$mean_ch4[i], sp_ch4$n[i]))
}

# O2 stats
if ("O2_concentration" %in% names(ymf2021)) {
  o2_data_s6 <- ymf2021 %>% filter(!is.na(O2_concentration))
  stat("O2 range", paste(round(min(o2_data_s6$O2_concentration / 10000, na.rm = TRUE), 1), "to",
                          round(max(o2_data_s6$O2_concentration / 10000, na.rm = TRUE), 1)), "%")
}

# --- Figure S8: Internal gas profiles (4 panels) ---
sub_header("Figure S8: Internal gas correlation panels")

# Panel A: CH4 vs CO2
ch4_co2 <- ymf2021 %>% filter(!is.na(CH4_concentration), !is.na(CO2_concentration))
if (nrow(ch4_co2) > 3) {
  cor_ch4_co2 <- cor.test(log10(ch4_co2$CH4_concentration + 1),
                           log10(ch4_co2$CO2_concentration + 1))
  m_ch4_co2 <- lm(log10(CH4_concentration + 1) ~ log10(CO2_concentration + 1), data = ch4_co2)
  cat(sprintf("  (a) CH4 vs CO2: R2 = %.3f, p = %s (n=%d)\n",
              summary(m_ch4_co2)$r.squared,
              ifelse(cor_ch4_co2$p.value < 0.001, "< 0.001", sprintf("%.3f", cor_ch4_co2$p.value)),
              nrow(ch4_co2)))
  record("internal_ch4_co2_r2", summary(m_ch4_co2)$r.squared)
}

# Panel B: CH4 vs O2
ch4_o2 <- ymf2021 %>% filter(!is.na(CH4_concentration), !is.na(O2_concentration))
if (nrow(ch4_o2) > 3) {
  cor_ch4_o2 <- cor.test(log10(ch4_o2$CH4_concentration + 1),
                          log10(ch4_o2$O2_concentration + 1))
  m_ch4_o2 <- lm(log10(CH4_concentration + 1) ~ log10(O2_concentration + 1), data = ch4_o2)
  cat(sprintf("  (b) CH4 vs O2: R2 = %.3f, p = %s (n=%d)\n",
              summary(m_ch4_o2)$r.squared,
              ifelse(cor_ch4_o2$p.value < 0.001, "< 0.001", sprintf("%.3f", cor_ch4_o2$p.value)),
              nrow(ch4_o2)))
  record("internal_ch4_o2_r2", summary(m_ch4_o2)$r.squared)
}

# Panel C: CH4 concentration vs heartwood mcrA
ch4_mcra <- ymf2021 %>%
  filter(!is.na(CH4_concentration), !is.na(ddpcr_mcra_probe_Inner_loose)) %>%
  mutate(mcra_copies = ddpcr_mcra_probe_Inner_loose + 1)
if (nrow(ch4_mcra) > 3) {
  cor_ch4_mcra <- cor.test(log10(ch4_mcra$mcra_copies),
                            log10(ch4_mcra$CH4_concentration + 1))
  m_ch4_mcra <- lm(log10(CH4_concentration + 1) ~ log10(mcra_copies), data = ch4_mcra)
  cat(sprintf("  (c) CH4 conc vs mcrA: R2 = %.3f, p = %s (n=%d)\n",
              summary(m_ch4_mcra)$r.squared,
              ifelse(cor_ch4_mcra$p.value < 0.001, "< 0.001", sprintf("%.3f", cor_ch4_mcra$p.value)),
              nrow(ch4_mcra)))
  record("internal_ch4_mcra_r2", summary(m_ch4_mcra)$r.squared)
}

# Panel D: CH4 flux vs CH4 concentration at each height
sub_header("  Panel D: Flux vs internal CH4 by height")
for (h_col in c("CH4_best.flux_50cm", "CH4_best.flux_125cm", "CH4_best.flux_200cm")) {
  h_label <- gsub("CH4_best.flux_", "", gsub("cm", " cm", h_col))
  flux_ch4 <- ymf2021 %>%
    filter(!is.na(CH4_concentration), !is.na(.data[[h_col]])) %>%
    mutate(flux_val = .data[[h_col]])
  if (nrow(flux_ch4) > 3) {
    m_fc <- lm(log10(abs(flux_val) + 0.01) ~ log10(CH4_concentration + 1), data = flux_ch4)
    cor_fc <- cor.test(log10(flux_ch4$CH4_concentration + 1),
                        log10(abs(flux_ch4$flux_val) + 0.01))
    cat(sprintf("    %s: R2 = %.3f, p = %s (n=%d)\n",
                h_label, summary(m_fc)$r.squared,
                ifelse(cor_fc$p.value < 0.001, "< 0.001", sprintf("%.3f", cor_fc$p.value)),
                nrow(flux_ch4)))
  }
}

# --- Figure S9: δ13CH4 isotope analysis ---
sub_header("Figure S9: d13CH4 isotope analysis")
# isotopes.csv: whole-tree set flagged by the rule in code/lib/isotope_samples.R; species from trees.csv
{
  internal_iso <- isotope_whole_tree() %>% dplyr::select(SampleName, species, d13CH4, ch4_ppm)
  atm_iso      <- isotope_atmosphere() %>% dplyr::mutate(species = "Atmosphere")
  iso_s8       <- dplyr::bind_rows(internal_iso, atm_iso)

  stat("Total isotope samples", nrow(iso_s8))
  stat("Internal tree samples", nrow(internal_iso))
  stat("Atmosphere samples", nrow(atm_iso))
  stat("Number of species", n_distinct(internal_iso$species))

  stat("Internal d13CH4 mean", round(mean(internal_iso$d13CH4), 1), "permil VPDB")
  stat("Internal d13CH4 median", round(median(internal_iso$d13CH4), 1), "permil VPDB")
  stat("Internal d13CH4 range", paste(round(min(internal_iso$d13CH4), 1), "to",
                                       round(max(internal_iso$d13CH4), 1)), "permil VPDB")
  stat("Internal d13CH4 SD", round(sd(internal_iso$d13CH4), 1), "permil")
  n_below70 <- sum(internal_iso$d13CH4 <= -70)
  stat("Samples <= -70 permil", sprintf("%d of %d (%.1f%%)",
       n_below70, nrow(internal_iso), 100 * n_below70 / nrow(internal_iso)))
  stat("Atmosphere d13CH4 mean", round(mean(atm_iso$d13CH4), 1), "permil VPDB")

  # Fraction in each pathway range
  n_hydro <- sum(internal_iso$d13CH4 >= -110 & internal_iso$d13CH4 <= -60)
  n_aceto <- sum(internal_iso$d13CH4 >= -65 & internal_iso$d13CH4 <= -50)
  n_methyl <- sum(internal_iso$d13CH4 >= -70 & internal_iso$d13CH4 <= -50)
  cat(sprintf("  In hydrogenotrophic range (-110 to -60): %d (%.0f%%)\n",
              n_hydro, 100 * n_hydro / nrow(internal_iso)))
  cat(sprintf("  In acetoclastic range (-65 to -50): %d (%.0f%%)\n",
              n_aceto, 100 * n_aceto / nrow(internal_iso)))

  # Keeling plot intercept
  keeling <- internal_iso %>% filter(ch4_ppm > 5) %>% mutate(inv_ch4 = 1 / ch4_ppm)
  if (nrow(keeling) > 3) {
    keeling_lm <- lm(d13CH4 ~ inv_ch4, data = keeling)
    stat("Keeling intercept (OLS, CH4 > 5 ppm; legacy check)", round(coef(keeling_lm)[1], 1), "permil")
    stat("Keeling R2", round(summary(keeling_lm)$r.squared, 3))
    stat("Canonical source estimates", "see outputs/data/ISOTOPES_summary.csv (MM robust Keeling)")
  }
}

# --- Figure S15: RF model diagnostics (expanded) ---
sub_header("Figure S15: RF model diagnostics (expanded)")
if (exists("TreeRF") && exists("SoilRF")) {
  # CCC requires DescTools; compute manually if not available
  # Already nmol m-2 s-1. The x1000 that stood here was a umol-era conversion and made
  # every RMSE printed below 1000x too large, on lines labelled "nmol m-2 s-1".
  tree_obs <- tree_train_complete$stem_flux_corrected
  tree_pred <- tree_train_complete$pred_flux
  soil_obs <- soil_train_complete$soil_flux_umol_m2_s
  soil_pred <- soil_train_complete$pred_flux

  # Manual CCC: 2 * r * sx * sy / (sx^2 + sy^2 + (mx - my)^2)
  ccc_manual <- function(x, y) {
    mx <- mean(x, na.rm = TRUE); my <- mean(y, na.rm = TRUE)
    sx <- sd(x, na.rm = TRUE); sy <- sd(y, na.rm = TRUE)
    r <- cor(x, y, use = "complete.obs")
    2 * r * sx * sy / (sx^2 + sy^2 + (mx - my)^2)
  }

  tree_ccc <- ccc_manual(tree_obs, tree_pred)
  soil_ccc <- ccc_manual(soil_obs, soil_pred)
  tree_rmse <- sqrt(mean((tree_pred - tree_obs)^2, na.rm = TRUE))
  soil_rmse <- sqrt(mean((soil_pred - soil_obs)^2, na.rm = TRUE))

  cat(sprintf("  Tree: CCC = %.3f, R2 = %.3f, RMSE = %.1f nmol m-2 s-1, n = %d\n",
              tree_ccc, cor(tree_obs, tree_pred)^2, tree_rmse, length(tree_obs)))
  cat(sprintf("  Soil: CCC = %.3f, R2 = %.3f, RMSE = %.1f nmol m-2 s-1, n = %d\n",
              soil_ccc, cor(soil_obs, soil_pred)^2, soil_rmse, length(soil_obs)))
  record("rf_tree_ccc", tree_ccc)
  record("rf_soil_ccc", soil_ccc)

  # Feature importance rankings
  sub_header("  Tree RF top features")
  tree_imp_df <- data.frame(feature = names(importance(TreeRF)),
                             imp = as.numeric(importance(TreeRF))) %>%
    arrange(desc(imp))
  for (i in 1:min(8, nrow(tree_imp_df))) {
    cat(sprintf("    %d. %s (%.2e)\n", i, tree_imp_df$feature[i], tree_imp_df$imp[i]))
  }

  sub_header("  Soil RF top features")
  soil_imp_df <- data.frame(feature = names(importance(SoilRF)),
                              imp = as.numeric(importance(SoilRF))) %>%
    arrange(desc(imp))
  for (i in 1:min(8, nrow(soil_imp_df))) {
    cat(sprintf("    %d. %s (%.2e)\n", i, soil_imp_df$feature[i], soil_imp_df$imp[i]))
  }
}

# --- Figure S11: already covered in Sections 8 + 11 ---
sub_header("Figure S11: Scale-dependent genes (see Sections 8-9)")
cat("  Individual gene-flux R2: see Section 8\n")
cat("  Species-level correlations: see Section 9\n")

# --- Figure S13: Tree radial mcrA cross-sections ---
sub_header("Figure S13: Tree radial mcrA cross-sections")
# mcrA by species (inner vs outer)
sp_mcra_comps <- ymf2021 %>%
  filter(!is.na(species_id)) %>%
  group_by(species_id) %>%
  summarise(
    n = n(),
    mean_inner = mean(ddpcr_mcra_probe_Inner_loose, na.rm = TRUE),
    mean_outer = mean(ddpcr_mcra_probe_Outer_loose, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n >= 3) %>%
  mutate(species_latin = SPECIES_NAMES[species_id],
         inner_outer_ratio = mean_inner / pmax(mean_outer, 1))

stat("Species with >= 3 trees", nrow(sp_mcra_comps))
cat("  Species with highest inner/outer mcrA ratio:\n")
sp_mcra_comps <- sp_mcra_comps %>% arrange(desc(inner_outer_ratio))
for (i in 1:min(5, nrow(sp_mcra_comps))) {
  cat(sprintf("    %s: inner=%.0f, outer=%.0f, ratio=%.1f\n",
              sp_mcra_comps$species_latin[i],
              sp_mcra_comps$mean_inner[i], sp_mcra_comps$mean_outer[i],
              sp_mcra_comps$inner_outer_ratio[i]))
}

# --- Figure S10: pmoA vs mmoX patterns (already in Section 11) ---
sub_header("Figure S10: pmoA/mmoX patterns (see Section 11)")
# Add S12-specific stats from manuscript text
# pmoA detection, mmoX detection, co-occurrence
pmoa_inner <- ymf2021$ddpcr_pmoa_Inner_loose
mmox_inner <- ymf2021$ddpcr_mmox_Inner_loose
pmoa_outer <- ymf2021$ddpcr_pmoa_Outer_loose
mmox_outer <- ymf2021$ddpcr_mmox_Outer_loose

# Total wood samples with either gene
all_pmoa <- c(pmoa_inner[!is.na(pmoa_inner)], pmoa_outer[!is.na(pmoa_outer)])
all_mmox <- c(mmox_inner[!is.na(mmox_inner)], mmox_outer[!is.na(mmox_outer)])
stat("pmoA detected in wood samples", sum(all_pmoa > 0))
stat("mmoX detected in wood samples", sum(all_mmox > 0))

# Co-occurrence (paired samples)
for (loc in c("Inner", "Outer")) {
  p_col <- paste0("ddpcr_pmoa_", loc, "_loose")
  m_col <- paste0("ddpcr_mmox_", loc, "_loose")
  if (p_col %in% names(ymf2021) && m_col %in% names(ymf2021)) {
    both <- ymf2021 %>% filter(!is.na(.data[[p_col]]), !is.na(.data[[m_col]]))
    either <- both %>% filter(.data[[p_col]] > 0 | .data[[m_col]] > 0)
    co_occur <- both %>% filter(.data[[p_col]] > 0 & .data[[m_col]] > 0)
    if (nrow(either) > 0) {
      cat(sprintf("  %s co-occurrence: %d/%d (%.0f%%) of samples with either gene\n",
                  c(Inner = "Heartwood", Outer = "Sapwood")[loc],
                  nrow(co_occur), nrow(either), 100 * nrow(co_occur) / nrow(either)))
    }
    # pmoA:mmoX ratio (when both positive)
    if (nrow(co_occur) > 0) {
      ratios <- co_occur[[p_col]] / co_occur[[m_col]]
      cat(sprintf("    pmoA:mmoX median ratio = %.1f-fold (log10 = %.2f +/- %.2f)\n",
                  median(ratios), median(log10(ratios)), sd(log10(ratios))))
    }
  }
}

# HW vs SW ratio test
hw_ratios <- ymf2021 %>%
  filter(!is.na(ddpcr_pmoa_Inner_loose), !is.na(ddpcr_mmox_Inner_loose),
         ddpcr_pmoa_Inner_loose > 0, ddpcr_mmox_Inner_loose > 0) %>%
  mutate(ratio = log10(ddpcr_pmoa_Inner_loose / ddpcr_mmox_Inner_loose))
sw_ratios <- ymf2021 %>%
  filter(!is.na(ddpcr_pmoa_Outer_loose), !is.na(ddpcr_mmox_Outer_loose),
         ddpcr_pmoa_Outer_loose > 0, ddpcr_mmox_Outer_loose > 0) %>%
  mutate(ratio = log10(ddpcr_pmoa_Outer_loose / ddpcr_mmox_Outer_loose))
if (nrow(hw_ratios) > 3 && nrow(sw_ratios) > 3) {
  wt_ratio <- wilcox.test(hw_ratios$ratio, sw_ratios$ratio)
  cat(sprintf("  HW vs SW ratio test: p = %.3f\n", wt_ratio$p.value))
}

# pmoA:mmoX ratio vs methanogen abundance
if (nrow(hw_ratios) > 3) {
  hw_with_mcra <- hw_ratios %>%
    left_join(ymf2021 %>% dplyr::select(tree_id, mcra = ddpcr_mcra_probe_Inner_loose),
              by = "tree_id") %>%
    filter(!is.na(mcra), mcra > 0)
  if (nrow(hw_with_mcra) > 3) {
    cor_ratio_mcra <- cor.test(hw_with_mcra$ratio, log10(hw_with_mcra$mcra))
    cat(sprintf("  Ratio vs mcrA: r = %.2f, p = %.3f\n",
                cor_ratio_mcra$estimate, cor_ratio_mcra$p.value))
  }
}

# --- Figure S14: mcrA vs methanotroph independence (already in Section 9) ---
sub_header("Figure S14: mcrA vs methanotroph (see Section 9)")
cat("  Tree-level and species-level correlations reported above.\n")

# --- Figure S12: Black oak methanome heatmap ---
sub_header("Figure S12: Black oak tissue methanome heatmap")
# This is from 10_black_oak_methanome_heatmap.R — qualitative, but key result:
# Methanobacteriaceae only in heartwood, no methanogens in soil
cat("  Key qualitative results (from manuscript text):\n")
cat("  - Only Methanobacterium (Methanobacteriaceae) detected as methanogen\n")
cat("  - Found exclusively in heartwood and sapwood\n")
cat("  - No methanogens in bark, branches, foliage, or soil\n")
cat("  - Methanotrophs (Methylocella) in heartwood, sapwood, bark, foliage\n")
cat("  - Methyloferula in sapwood, bark, leaf litter\n")
cat("  - Methylorosula in heartwood, sapwood, leaf litter\n")

# --- Figure S2: Taxonomy x pmoA heatmap ---
sub_header("Figure S2: Family-level 16S x pmoA associations")
picrust_pmoa_full <- "outputs/data/picrust/pathway_associations_pmoa.csv"
{
  pa_pmoa_f <- read.csv(picrust_pmoa_full, stringsAsFactors = FALSE)
  sig_pmoa_f <- pa_pmoa_f %>% filter(!is.na(FDR), FDR < 0.01)
  stat("Total pathways tested (pmoA)", nrow(pa_pmoa_f))
  stat("Significant at FDR < 0.01 (pmoA)", nrow(sig_pmoa_f))
}

# Also check taxonomy heatmap for S6
taxonomy_mcra_file <- "outputs/data/picrust/pathway_associations_mcra_no_mcra_otus.csv"
sub_header("Figure S6: Family-level 16S x mcrA associations")
if (exists("fam_cors") && nrow(fam_cors) > 0) {
  stat("Families significantly correlated with mcrA (p<0.05)", nrow(fam_cors))
  stat("Positively correlated", sum(fam_cors$r > 0))
  stat("Negatively correlated", sum(fam_cors$r < 0))
} else {
  cat("  [See Section 5 family-mcrA correlations]\n")
}


# ==============================================================================
# FOOTER: CHANGE DETECTION SUMMARY
# ==============================================================================

section_header("CHANGE DETECTION SUMMARY")
cat("  (Machine-readable key-value pairs for diffing between runs)\n\n")
for (nm in sort(names(summary_stats))) {
  val <- summary_stats[[nm]]
  if (is.numeric(val)) {
    cat(sprintf("  %-45s %.6f\n", nm, val))
  } else {
    cat(sprintf("  %-45s %s\n", nm, as.character(val)))
  }
}

cat("\n", strrep("*", 72), "\n")
cat("  DONE — ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("*", 72), "\n")

# --- Close output file ---
sink()
cat(sprintf("\nResults saved to: %s\n", normalizePath(output_file)))
