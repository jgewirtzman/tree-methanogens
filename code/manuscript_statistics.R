# ==============================================================================
# Manuscript Statistics — All Quantitative Results
# ==============================================================================
# Purpose: Recomputes every quantitative result reported in the manuscript
#   from underlying data. Organized by manuscript section for easy
#   cross-referencing. No figures are generated — console output only.
#
# Usage: Open tree-methanogens.Rproj in RStudio, then:
#   source("code/manuscript_statistics.R")
#
# Output: Prints all statistics to console. Ends with a machine-readable
#   summary for easy diff between revisions.
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

cat("\n")
cat(strrep("*", 72), "\n")
cat("  MANUSCRIPT STATISTICS — Tree Microbiomes & Methane\n")
cat("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(strrep("*", 72), "\n")

# --- Shared constants ---
NMOL_TO_UG <- 57.744  # nmol m-2 s-1 -> ug CH4 m-2 hr-1

species_mapping <- c(
  "ACRU" = "Acer rubrum", "ACSA" = "Acer saccharum",
  "BEAL" = "Betula alleghaniensis", "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera", "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana", "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra", "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata", "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina", "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina", "SAAL" = "Sassafras albidum"
)

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

area_weighted_gene <- function(gene_inner, gene_outer, dbh_cm) {
  R <- dbh_cm / 2
  if (!is.finite(R) || R <= 0) return(NA_real_)
  r1 <- min(5, R)
  r2 <- max(R - 5, r1)
  S <- r1^2 + r1 * r2 + r2^2
  gene_outer + (gene_inner - gene_outer) * (S / (3 * R^2))
}

# --- Summary collector for footer ---
summary_stats <- list()
record <- function(name, value) {
  summary_stats[[name]] <<- value
}

# --- Load ALL data upfront ---
cat("\nLoading data...\n")

semirigid_tree <- read.csv("data/processed/flux/semirigid_tree_final_complete_dataset.csv")
semirigid_soil <- read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv")
moisture_data  <- read.csv("data/raw/field_data/ipad_data/Cleaned data/soilmoisture_total.csv",
                           fileEncoding = "UTF-8-BOM")

ymf2021 <- read.csv("data/processed/integrated/merged_tree_dataset_final.csv")
ymf2023 <- read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv")

aux      <- read.csv("data/processed/flux/goflux_auxfile.csv")
ch4_flux <- read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv")

monthly_fluxes <- read.csv("outputs/tables/MONTHLY_FLUXES.csv")
annual_summary <- read.csv("outputs/tables/ANNUAL_SUMMARY.csv")
load("outputs/models/RF_MODELS.RData")

cat("  All data loaded.\n")


# ==============================================================================
# SECTION 1: TEMPORAL FLUX PATTERNS (Figure 1)
# ==============================================================================

section_header("SECTION 1: TEMPORAL FLUX PATTERNS (Figure 1)")

# Quality filter (matches 06_soil_tree_timeseries.R)
soil_qc <- semirigid_soil %>%
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2 >= 0.7 | CH4_HM.r2 >= 0.7),
         CH4_best.flux >= -100 & CH4_best.flux <= 200) %>%
  mutate(Plot_Type = case_when(
    Plot.letter %in% c("WD", "WS") ~ "W",
    Plot.letter == "I" ~ "I",
    Plot.letter == "U" ~ "U",
    TRUE ~ Plot.letter
  ))

tree_qc <- semirigid_tree %>%
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2.x >= 0.7 | CH4_HM.r2.x >= 0.7),
         CH4_best.flux.x >= -100 & CH4_best.flux.x <= 200) %>%
  mutate(
    Plot_Type = case_when(
      Plot.Letter %in% c("WD", "WS") ~ "W",
      Plot.Letter == "I" ~ "I",
      Plot.Letter == "U" ~ "U",
      TRUE ~ Plot.Letter
    ),
    CH4_flux = CH4_best.flux.x
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
moisture_data$Date <- as.Date(moisture_data$Date)
moisture_data <- moisture_data %>%
  mutate(Plot_Type = case_when(
    grepl("^U", Plot) | Plot == "U" ~ "U",
    grepl("^I", Plot) | Plot == "I" ~ "I",
    grepl("^W", Plot) | Plot %in% c("WD", "WS", "W") ~ "W",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Plot_Type), !is.na(VWC))

for (pos in c("U", "I", "W")) {
  d <- moisture_data %>% filter(Plot_Type == pos)
  pos_label <- c(U = "Upland", I = "Intermediate", W = "Transitional wetland")[pos]
  cat(sprintf("  %s: %.1f +/- %.1f%% (n=%d)\n",
              pos_label, mean(d$VWC, na.rm = TRUE), sd(d$VWC, na.rm = TRUE), nrow(d)))
}

# Soil temperature range
if ("SoilTemp" %in% names(semirigid_soil) || "Soil.Temp" %in% names(semirigid_soil)) {
  temp_col <- intersect(c("SoilTemp", "Soil.Temp", "Soil_Temp"), names(semirigid_soil))[1]
  if (!is.na(temp_col)) {
    temps <- semirigid_soil[[temp_col]]
    stat("Soil temperature range", paste(round(min(temps, na.rm = TRUE), 1),
                                          "to", round(max(temps, na.rm = TRUE), 1)), "deg C")
  }
}


# ==============================================================================
# SECTION 2: HEIGHT EFFECTS (Figure 2)
# ==============================================================================

section_header("SECTION 2: HEIGHT EFFECTS (Figure 2)")

# Prepare data (matches 04_height_effect_analysis.R)
height_data <- merge(ch4_flux, aux[, c("UniqueID", "measurement_height", "tree_id",
                                        "species", "plot")],
                     by = "UniqueID", all.x = TRUE)
names(height_data)[names(height_data) == "best.flux"] <- "CH4_best.flux"

height_data <- height_data %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux), !is.na(height_numeric)) %>%
  mutate(
    species = ifelse(species == "TSLA", "TSCA", species),
    height_m = height_numeric / 100,
    tree_unique = paste(plot, tree_id, sep = "_"),
    species_latin = species_mapping[species]
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
              species_mapping[sig_sp$species[i]],
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
    cat(sprintf("  %s: %d%% iterations significant\n", species_mapping[sp], prop_sig))
  }
  # Within-species: proportion of trees with negative slopes
  tree_slopes <- d %>% group_by(tree_unique) %>% filter(n() >= 2) %>%
    summarise(slope = coef(lm(CH4_best.flux ~ height_numeric))[2], .groups = "drop")
  prop_neg <- round(100 * sum(tree_slopes$slope < 0) / nrow(tree_slopes), 0)
  cat(sprintf("  %s: %d/%d trees with negative slopes (%d%%)\n",
              species_mapping[sp], sum(tree_slopes$slope < 0), nrow(tree_slopes), prop_neg))
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

# Combine datasets (matches 04_variance_partition.R)
data_2023 <- ymf2023 %>%
  select(Species.Code, DBH = DBH..cm., Air_temp = air_temp_C,
         Soil_temp = Soil.Temp....C., VWC = vwc_mean, CH4_flux = CH4_best.flux) %>%
  mutate(Year = "2023", Species_Latin = species_mapping[Species.Code]) %>%
  drop_na(CH4_flux)

data_2021 <- ymf2021 %>%
  select(Species.Code = species_id, DBH = dbh, Air_temp = Temp_Air_125cm,
         Soil_temp = SoilTemp_mean, VWC = VWC_mean, CH4_flux = CH4_best.flux_125cm) %>%
  mutate(Year = "2021", Species_Latin = species_mapping[Species.Code]) %>%
  drop_na(CH4_flux) %>% filter(!is.nan(CH4_flux))

combined <- bind_rows(data_2023, data_2021) %>%
  filter(!is.na(Species_Latin)) %>%
  group_by(Species_Latin) %>% filter(n() > 3) %>% ungroup()

stat("Total observations", nrow(combined))
stat("2021 observations", sum(combined$Year == "2021"))
stat("2023 observations", sum(combined$Year == "2023"))
stat("Number of species", n_distinct(combined$Species_Latin))

# Standardize
for (v in c("DBH", "Air_temp", "Soil_temp", "VWC")) {
  if (v %in% names(combined))
    combined[[paste0(v, "_std")]] <- scale(combined[[v]])[, 1]
}

# Models
m_env <- lm(CH4_flux ~ DBH_std + Air_temp_std + Soil_temp_std + VWC_std, data = combined)
m_sp  <- lm(CH4_flux ~ Species_Latin, data = combined)
m_full <- lm(CH4_flux ~ DBH_std + Air_temp_std + Soil_temp_std + VWC_std + Species_Latin, data = combined)
m_int  <- lm(CH4_flux ~ (DBH_std + Air_temp_std + Soil_temp_std + VWC_std) * Species_Latin, data = combined)

r2_env <- summary(m_env)$r.squared
r2_sp  <- summary(m_sp)$r.squared
r2_full <- summary(m_full)$r.squared
r2_int  <- summary(m_int)$r.squared

env_pct <- max(0, r2_full - r2_sp) * 100
sp_pct  <- max(0, r2_full - r2_env) * 100
int_pct <- max(0, r2_int - r2_full) * 100
unexp_pct <- (1 - r2_int) * 100

sub_header("Variance partitioning (Method 2)")
stat("Environment unique", env_pct, "%")
stat("Species unique", sp_pct, "%")
stat("Env x Species interaction", int_pct, "%")
stat("Unexplained", unexp_pct, "%")
record("variance_species_pct", sp_pct)
record("variance_interaction_pct", int_pct)
record("variance_unexplained_pct", unexp_pct)

# Top species mean flux
sub_header("Species mean fluxes (nmol m-2 s-1)")
sp_means <- combined %>%
  group_by(Species_Latin) %>%
  summarise(mean_flux = mean(CH4_flux), se_flux = sd(CH4_flux) / sqrt(n()),
            n = n(), .groups = "drop") %>%
  arrange(desc(mean_flux))

for (i in 1:min(5, nrow(sp_means))) {
  cat(sprintf("  %s: %.3f +/- %.3f (n=%d)\n",
              sp_means$Species_Latin[i], sp_means$mean_flux[i],
              sp_means$se_flux[i], sp_means$n[i]))
}


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
    stat("    % exceeding 500 copies/g", round(100 * sum(values > 500) / n_total, 1), "%")
  }
}

sub_header("mcrA (probe, loose) by compartment")
gene_stats(ymf2021$ddpcr_mcra_probe_Inner_loose, "Heartwood")
gene_stats(ymf2021$ddpcr_mcra_probe_Outer_loose, "Sapwood")
gene_stats(ymf2021$ddpcr_mcra_probe_Mineral_loose, "Mineral soil")
gene_stats(ymf2021$ddpcr_mcra_probe_Organic_loose, "Organic soil")

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
    total <- ifelse(!is.na(p) & !is.na(m), p + m,
                    ifelse(!is.na(p), p, ifelse(!is.na(m), m, NA)))
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
ddpcr_raw <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv")
water <- read.csv("data/raw/tree_cores/Tree_Core_Sectioning_Data.csv")
qpcr_16s <- read.csv("data/raw/16s/16s_w_metadata.csv")
qpcr_16s <- subset(qpcr_16s, qpcr_16s$Sample.ID != "None")
qpcr_16s <- qpcr_16s[, c(3, 4, 6)]

water$seq_id <- toupper(water$seq_id)
ddpcr_raw <- merge(ddpcr_raw, water, by = c("core_type", "seq_id"), all.x = TRUE)
ddpcr_raw <- merge(ddpcr_raw, qpcr_16s, by = c("core_type", "seq_id"), all.x = TRUE)

otu_tab <- read.delim("data/raw/16s/OTU_table.txt", header = TRUE, row.names = 1)
bastard_tax <- otu_tab[, 590:596]
bastard_tax[bastard_tax == ""] <- NA
tax_tab_pre <- tax_table(bastard_tax)
taxa_names(tax_tab_pre) <- sub("sp", "seq", taxa_names(tax_tab_pre))
otu_tab_corr <- otu_tab[, 1:589]
otu_table_pre <- otu_table(otu_tab_corr, taxa_are_rows = TRUE)

phylo_tree <- read_tree("data/raw/16s/unrooted_tree.nwk")
samp_data <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", row.names = 1)
samp_data$RowName <- row.names(samp_data)
samp_data$seq_id <- sub("prime", "'", samp_data$seq_id)
samp_data$seq_id <- sub("star", "*", samp_data$seq_id)
samp_data$seq_id <- sub("HM", "H", samp_data$seq_id)
samp_data$seq_id[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "RO104"
samp_data$core_type[samp_data$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"] <- "Inner"

samp_data_merged <- merge(ddpcr_raw, samp_data, by = c("seq_id", "core_type"), all.y = TRUE)
dups <- which(duplicated(samp_data_merged$RowName))
samp_data_merged <- samp_data_merged[-c(dups), ]
row.names(samp_data_merged) <- samp_data_merged$RowName

raw_ps <- phyloseq(tax_tab_pre, otu_table_pre, phylo_tree, sample_data(samp_data_merged))

pop_taxa <- function(physeq, badTaxa) {
  allTaxa <- taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  prune_taxa(allTaxa, physeq)
}

mitochondria <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 5] == "Mitochondria")]
chloroplast <- rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[, 4] == "Chloroplast")]
no_mito <- pop_taxa(raw_ps, c(mitochondria, chloroplast))
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

# Methanogen families
methanogen_families <- c(
  "Methanobacteriaceae", "Methanocellaceae", "Methanomassiliicoccaceae",
  "Methanomicrobiaceae", "Methanosaetaceae", "Methanosarcinaceae",
  "Methanospirillaceae", "Methanotrichaceae", "Methanomethylophilaceae",
  "Methanoperedenaceae", "Methanofastidiosaceae"
)

methanogen_asvs <- rownames(tax_df)[tax_df$Family %in% methanogen_families]
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
  if (length(asvs) > 0) {
    pcts <- colSums(otu_df[asvs, rownames(sw_samples), drop = FALSE])
    cat(sprintf("  %s: range %.2f to %.2f%%\n", fam, min(pcts), max(pcts)))
  }
}

# Soil taxa (VERIFY Bathyarchaeia dominance)
sub_header("Top soil archaeal taxa (Mineral + Organic)")
soil_samples <- samp_meta %>% filter(compartment %in% c("Mineral Soil", "Organic Soil"))
all_classes <- tax_df$Class
archaea_asvs <- rownames(tax_df)[tax_df$Kingdom == "Archaea"]
if (length(archaea_asvs) > 0) {
  archaea_class_pcts <- data.frame()
  for (cls in unique(tax_df[archaea_asvs, "Class"])) {
    if (is.na(cls)) next
    asvs <- archaea_asvs[tax_df[archaea_asvs, "Class"] == cls]
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
    if (length(asvs) > 0) {
      pcts <- colSums(otu_df[asvs, rownames(hw_with_mcra), drop = FALSE])
      if (sd(pcts) > 0) {
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
mt_defs <- read.csv("data/processed/molecular/methanotroph_definitions.csv")
known_families <- mt_defs %>% filter(Taxon_rank == "Family", Include_known == "YES") %>% pull(Taxon)
known_genera <- mt_defs %>% filter(Taxon_rank == "Genus", Include_known == "YES") %>% pull(Taxon)
putative_families <- mt_defs %>% filter(Taxon_rank == "Family", Include_putative == "YES") %>% pull(Taxon)

# Known methanotrophs
known_asvs <- c(
  rownames(tax_df)[tax_df$Family %in% known_families],
  rownames(tax_df)[tax_df$Genus %in% known_genera]
)
known_asvs <- unique(known_asvs)

# Putative (family-level only for unresolved genera)
putative_only <- rownames(tax_df)[
  tax_df$Family %in% putative_families &
  !(tax_df$Genus %in% known_genera) &
  (is.na(tax_df$Genus) | tax_df$Genus == "")
]

cat(sprintf("  Known methanotroph ASVs: %d\n", length(known_asvs)))
cat(sprintf("  Putative methanotroph ASVs: %d\n", length(putative_only)))

# List families detected
mt_families_detected <- sort(unique(tax_df[c(known_asvs, putative_only), "Family"]))
cat("  Detected methanotroph families:\n")
for (fam in mt_families_detected) {
  asvs <- intersect(c(known_asvs, putative_only), rownames(tax_df)[tax_df$Family == fam])
  pcts <- colSums(otu_df[asvs, rownames(hw_samples), drop = FALSE])
  cat(sprintf("    %s: mean %.3f%% in heartwood (n_ASVs=%d)\n", fam, mean(pcts), length(asvs)))
}


# ==============================================================================
# SECTION 6: FAPROTAX / PICRUSt (Figures S2-S5)
# ==============================================================================

section_header("SECTION 6: FAPROTAX / PICRUSt (Figures S2-S5)")

# FAPROTAX
sub_header("FAPROTAX functional predictions")
faprotax_file <- "data/raw/16s/faprotax_output/functional_table.tsv"
if (file.exists(faprotax_file)) {
  fap <- read.delim(faprotax_file, row.names = 1)
  # Match to rarefied samples
  common_samps <- intersect(colnames(fap), rownames(samp_meta))
  fap_sub <- fap[, common_samps]
  fap_ra <- sweep(fap_sub, 2, colSums(otu_df[, common_samps]), FUN = "/") * 100

  for (func in c("fermentation", "dark_hydrogen_oxidation", "methylotrophy",
                  "methanogenesis", "methane_oxidation", "methanotrophy")) {
    if (func %in% rownames(fap_ra)) {
      for (comp in c("Heartwood", "Sapwood")) {
        samps <- samp_meta %>% filter(compartment == comp) %>% rownames()
        samps <- intersect(samps, colnames(fap_ra))
        if (length(samps) > 0) {
          vals <- as.numeric(fap_ra[func, samps])
          cat(sprintf("  %s in %s: %.2f +/- %.2f%%\n",
                      func, comp, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE)))
        }
      }
      # Wilcoxon test hw vs sw
      hw_s <- samp_meta %>% filter(compartment == "Heartwood") %>% rownames()
      sw_s <- samp_meta %>% filter(compartment == "Sapwood") %>% rownames()
      hw_s <- intersect(hw_s, colnames(fap_ra))
      sw_s <- intersect(sw_s, colnames(fap_ra))
      if (length(hw_s) > 3 && length(sw_s) > 3) {
        wt <- wilcox.test(as.numeric(fap_ra[func, hw_s]), as.numeric(fap_ra[func, sw_s]))
        cat(sprintf("    HW vs SW p-value: %.4f\n", wt$p.value))
      }
    }
  }

  # A. saccharum peaks
  acsa_hw_s <- samp_meta %>% filter(compartment == "Heartwood", species.x == "ACSA") %>% rownames()
  acsa_hw_s <- intersect(acsa_hw_s, colnames(fap_ra))
  if (length(acsa_hw_s) > 0) {
    sub_header("A. saccharum heartwood FAPROTAX peaks")
    for (func in c("dark_hydrogen_oxidation", "methylotrophy")) {
      if (func %in% rownames(fap_ra)) {
        vals <- as.numeric(fap_ra[func, acsa_hw_s])
        cat(sprintf("  %s: %.2f%%\n", func, mean(vals, na.rm = TRUE)))
      }
    }
  }

  # VERIFY function names: list all >1% in heartwood
  sub_header("All FAPROTAX functions >1% mean in heartwood")
  hw_all <- samp_meta %>% filter(compartment == "Heartwood") %>% rownames()
  hw_all <- intersect(hw_all, colnames(fap_ra))
  func_means <- rowMeans(fap_ra[, hw_all, drop = FALSE], na.rm = TRUE)
  func_means <- sort(func_means[func_means > 1], decreasing = TRUE)
  for (fn in names(func_means)) {
    cat(sprintf("  %s: %.2f%%\n", fn, func_means[fn]))
  }
} else {
  cat("  [SKIPPED] FAPROTAX output file not found.\n")
}

# PICRUSt pathway verification
sub_header("PICRUSt pathway associations (VERIFY pathway names)")
picrust_mcra <- "data/processed/molecular/picrust/pathway_associations_mcra_all.csv"
if (file.exists(picrust_mcra)) {
  pvals_mcra <- read.csv(picrust_mcra)
  sig_mcra <- pvals_mcra %>% filter(fdr < 0.01) %>% arrange(fdr)
  cat(sprintf("  Significant mcrA-associated pathways (FDR < 0.01): %d\n", nrow(sig_mcra)))
  cat("  Top 15 pathways:\n")
  for (i in 1:min(15, nrow(sig_mcra))) {
    direction <- ifelse(sig_mcra$t_value[i] > 0, "+", "-")
    cat(sprintf("    %s %s (FDR=%.4f, t=%.2f)\n",
                direction, sig_mcra$pathway[i], sig_mcra$fdr[i], sig_mcra$t_value[i]))
  }
}

picrust_pmoa <- "data/processed/molecular/picrust/pathway_associations_pmoa.csv"
if (file.exists(picrust_pmoa)) {
  pvals_pmoa <- read.csv(picrust_pmoa)
  sig_pmoa <- pvals_pmoa %>% filter(fdr < 0.01) %>% arrange(fdr)
  cat(sprintf("\n  Significant pmoA-associated pathways (FDR < 0.01): %d\n", nrow(sig_pmoa)))
  cat("  Top 10 pathways:\n")
  for (i in 1:min(10, nrow(sig_pmoa))) {
    direction <- ifelse(sig_pmoa$t_value[i] > 0, "+", "-")
    cat(sprintf("    %s %s (FDR=%.4f, t=%.2f)\n",
                direction, sig_pmoa$pathway[i], sig_pmoa$fdr[i], sig_pmoa$t_value[i]))
  }
}


# ==============================================================================
# SECTION 7: INDIVIDUAL GENE-FLUX (Figure S10 left)
# ==============================================================================

section_header("SECTION 7: INDIVIDUAL GENE-FLUX RELATIONSHIPS")

# Prepare area-weighted gene data (matches 02_scale_dependent_gene_patterns.R)
prepare_long_all_genes <- function(df) {
  df %>%
    dplyr::select(tree_id, starts_with("ddpcr_")) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    pivot_longer(cols = starts_with("ddpcr_"), names_to = "measurement_type", values_to = "gene_copies") %>%
    filter(!is.na(gene_copies)) %>%
    separate(measurement_type, into = c("method", "gene", "part1", "part2", "part3"),
             sep = "_", extra = "merge", fill = "right") %>%
    mutate(
      is_probe = (part1 == "probe"),
      location = if_else(is_probe, part2, part1),
      stringency = if_else(is_probe, part3, part2),
      location = str_remove(location, "probe"),
      location = case_when(location %in% c("Inner", "inner") ~ "Inner",
                           location %in% c("Outer", "outer") ~ "Outer", TRUE ~ location),
      gene = case_when(gene == "mcra" ~ "mcrA", gene == "mmox" ~ "mmoX",
                       gene == "pmoa" ~ "pmoA", TRUE ~ gene),
      sample_type = case_when(location == "Inner" ~ "Heartwood",
                              location == "Outer" ~ "Sapwood", TRUE ~ NA_character_)
    ) %>%
    filter(stringency == "loose", !is.na(sample_type), sample_type %in% c("Heartwood", "Sapwood"))
}

long_all <- prepare_long_all_genes(ymf2021)

tree_genes_weighted <- long_all %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, species_id, dbh), by = "tree_id") %>%
  filter(is.finite(dbh)) %>%
  group_by(tree_id, species_id, dbh, gene) %>%
  summarise(gene_inner = mean(gene_copies[sample_type == "Heartwood"], na.rm = TRUE),
            gene_outer = mean(gene_copies[sample_type == "Sapwood"], na.rm = TRUE),
            .groups = "drop") %>%
  filter(is.finite(gene_inner), is.finite(gene_outer)) %>%
  mutate(gene_area_weighted = mapply(area_weighted_gene, gene_inner, gene_outer, dbh),
         species = species_mapping[species_id]) %>%
  dplyr::select(tree_id, species_id, species, gene, gene_area_weighted) %>%
  pivot_wider(names_from = gene, values_from = gene_area_weighted) %>%
  mutate(methanotroph_total = case_when(
    is.na(pmoA) & !is.na(mmoX) ~ mmoX,
    !is.na(pmoA) & is.na(mmoX) ~ pmoA,
    !is.na(pmoA) & !is.na(mmoX) ~ pmoA + mmoX,
    TRUE ~ NA_real_
  ))

tree_genes_complete <- tree_genes_weighted %>%
  filter(!is.na(mcrA), !is.na(pmoA), !is.na(mmoX), !is.na(methanotroph_total))

tree_level <- tree_genes_complete %>%
  left_join(ymf2021 %>% dplyr::select(tree_id, CH4_flux = CH4_best.flux_125cm), by = "tree_id") %>%
  filter(!is.na(CH4_flux), !is.nan(CH4_flux)) %>%
  mutate(
    log_mcra = log10(mcrA + 1), log_pmoa = log10(pmoA + 1),
    log_mmox = log10(mmoX + 1), log_methanotroph = log10(methanotroph_total + 1),
    ratio = (mcrA + 1) / (methanotroph_total + 1), log_ratio = log10((mcrA + 1) / (methanotroph_total + 1))
  )

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


# ==============================================================================
# SECTION 8: SPECIES-LEVEL GENE-FLUX (Figure 8, S10 right, S13)
# ==============================================================================

section_header("SECTION 8: SPECIES-LEVEL GENE-FLUX RELATIONSHIPS")

# Flux by species (combined 2021 + 2023)
flux_all <- bind_rows(
  ymf2023 %>% dplyr::select(species_id = Species.Code, CH4_flux = CH4_best.flux) %>%
    mutate(species = species_mapping[species_id]),
  ymf2021 %>% dplyr::select(species_id, CH4_flux = CH4_best.flux_125cm) %>%
    mutate(species = species_mapping[species_id])
) %>% filter(!is.na(CH4_flux), !is.nan(CH4_flux), !is.na(species))

flux_by_sp <- flux_all %>%
  group_by(species, species_id) %>%
  summarise(n_flux = n(), median_flux = median(CH4_flux), .groups = "drop")

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
    summarise(med_mcra = median(mcrA), med_mt = median(methanotroph_total), .groups = "drop"),
  flux_by_sp %>% dplyr::select(species, n_flux), by = "species"
) %>% filter(n_flux >= 5)

if (nrow(sp_mcra_mt) >= 3) {
  sp_cor <- cor.test(log10(sp_mcra_mt$med_mcra + 1), log10(sp_mcra_mt$med_mt + 1))
  cat(sprintf("  Species-level: r = %.3f, p = %.4f (n=%d)\n",
              sp_cor$estimate, sp_cor$p.value, nrow(sp_mcra_mt)))
}


# ==============================================================================
# SECTION 9: RANDOM FOREST & UPSCALING (Figure 9)
# ==============================================================================

section_header("SECTION 9: RANDOM FOREST & UPSCALING (Figure 9)")

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
# Convert from umol to nmol: multiply by 1000
tree_monthly_nmol <- monthly_fluxes$Phi_tree_umol_m2_s * 1000
soil_monthly_nmol <- monthly_fluxes$Phi_soil_umol_m2_s * 1000

stat("Tree stem annual mean", mean(tree_monthly_nmol), "nmol m-2 s-1")
stat("Tree stem monthly range", paste(round(min(tree_monthly_nmol), 4), "to",
                                       round(max(tree_monthly_nmol), 4)), "nmol m-2 s-1")
stat("Soil annual mean", mean(soil_monthly_nmol), "nmol m-2 s-1")
stat("Soil monthly range", paste(round(min(soil_monthly_nmol), 4), "to",
                                  round(max(soil_monthly_nmol), 4)), "nmol m-2 s-1")
stat("Per-unit-area offset", paste0(round(abs(mean(tree_monthly_nmol) / mean(soil_monthly_nmol)) * 100, 1), "%"))

# Annual totals
sub_header("Annual totals")
stat("Tree emissions", annual_summary$annual_tree_mg_m2, "mg CH4 m-2 yr-1")
stat("Soil uptake", annual_summary$annual_soil_mg_m2, "mg CH4 m-2 yr-1")
record("annual_tree_mg_m2_yr", annual_summary$annual_tree_mg_m2)
record("annual_soil_mg_m2_yr", annual_summary$annual_soil_mg_m2)

# Surface area and offset
sub_header("Stem surface area and offset")
# Load inventory for geometry
tryCatch({
  load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
  inv <- rf_workflow_data$PLACEHOLDER_INVENTORY
  if (!is.null(inv) && "DBH" %in% names(inv)) {
    # Compute lateral surface area to 2m height
    inv$stem_area_m2 <- pi * (inv$DBH / 100) * 2  # pi * diameter * height
    total_stem_area <- sum(inv$stem_area_m2, na.rm = TRUE)
    plot_area_ha <- 10.2
    plot_area_m2 <- plot_area_ha * 10000
    stem_frac <- total_stem_area / plot_area_m2 * 100
    stat("Total stem surface area (0-2m)", round(total_stem_area, 0), "m2")
    stat("Stem area as % of plot", round(stem_frac, 2), "%")
  }
}, error = function(e) {
  cat("  [Using pre-computed values from annual summary]\n")
})

# From annual summary
tree_ann <- annual_summary$annual_tree_mg_m2
soil_ann <- annual_summary$annual_soil_mg_m2
offset_pct <- abs(tree_ann / soil_ann) * 100
stat("Plot-level offset", round(offset_pct, 2), "%")
record("plot_offset_pct", offset_pct)

# WAI extrapolation
WAI <- 3.07
# Per-stem-area flux = mean tree nmol * conversion to annual mg
mean_tree_nmol <- mean(tree_monthly_nmol)
# Annualize: nmol m-2 s-1 * 16e-9 g/nmol * 1e6 ug/g * 1e-3 mg/ug * 86400 s/d * 365.25 d/yr
wai_annual <- mean_tree_nmol * 16e-9 * 1e6 * 1e-3 * 86400 * 365.25 * WAI
stat("WAI extrapolation (annual)", round(wai_annual, 1), "mg CH4 m-2 yr-1")
wai_offset <- abs(wai_annual / soil_ann) * 100
stat("WAI offset of soil uptake", round(wai_offset, 1), "%")
record("wai_tree_mg_m2_yr", wai_annual)
record("wai_offset_pct", wai_offset)


# ==============================================================================
# SECTION 10: ADDITIONAL METHANOTROPH RESULTS (NEW)
# ==============================================================================

section_header("SECTION 10: ADDITIONAL METHANOTROPH RESULTS")
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

# Soil genes adding no value
sub_header("Soil genes beyond tree mmoX")
if ("log_mmox" %in% names(tree_level)) {
  m_base <- lm(CH4_flux ~ species + log_mmox, data = tree_level)
  r2_base <- summary(m_base)$r.squared
  for (soil_gene in c("ddpcr_mcra_probe_Mineral_loose", "ddpcr_pmoa_Mineral_loose", "ddpcr_mmox_Mineral_loose")) {
    if (soil_gene %in% names(ymf2021)) {
      tree_w_soil <- tree_level %>%
        left_join(ymf2021 %>% dplyr::select(tree_id, !!soil_gene), by = "tree_id")
      sg <- tree_w_soil[[soil_gene]]
      tree_w_soil$soil_gene_log <- log10(ifelse(is.na(sg) | sg <= 0, 1, sg))
      m_add <- lm(CH4_flux ~ species + log_mmox + soil_gene_log, data = tree_w_soil)
      delta <- summary(m_add)$r.squared - r2_base
      p_soil <- tryCatch(summary(m_add)$coefficients["soil_gene_log", "Pr(>|t|)"], error = function(e) NA)
      short_name <- sub("ddpcr_", "", sub("_Mineral_loose", "", soil_gene))
      cat(sprintf("  + soil %s: delta-R2 = %.4f, p = %.4f\n", short_name, delta, ifelse(is.na(p_soil), 1, p_soil)))
    }
  }
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
