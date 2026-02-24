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
if (file.exists("outputs/models/TRAINING_DATA.RData")) {
  load("outputs/models/TRAINING_DATA.RData")  # tree_train_complete, X_tree, X_soil, soil_train_complete
}

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
# SECTION 6b: PICRUSt PATHWAY-GENE HEATMAP (Figure 6)
# ==============================================================================

section_header("SECTION 6b: PICRUSt PATHWAY-GENE HEATMAP (Figure 6)")

# Figure 6 is the mcrA-pathway heatmap from 12b_picrust_pathway_heatmap.R
# Uses the no-mcrA-OTU pipeline: pathway abundances reconstructed after removing
# methanogen OTU contributions, then tested against mcrA via LMER.
# FDR threshold 1e-4 with mcrA-OTU contribution < 50%.

picrust_no_mcra_file <- "data/processed/molecular/picrust/pathway_associations_mcra_no_mcra_otus.csv"
picrust_combined_file <- "data/processed/molecular/picrust/pathway_associations_combined.csv"

if (file.exists(picrust_no_mcra_file)) {
  pvals_no_mcra <- read.csv(picrust_no_mcra_file, stringsAsFactors = FALSE)
  stat("Total pathways tested (no-mcrA-OTU pipeline)", nrow(pvals_no_mcra))

  sig_1e4 <- pvals_no_mcra %>% filter(!is.na(FDR), FDR < 1e-4)
  stat("Significant at FDR < 1e-4", nrow(sig_1e4))

  sig_001 <- pvals_no_mcra %>% filter(!is.na(FDR), FDR < 0.01)
  stat("Significant at FDR < 0.01", nrow(sig_001))

  # Apply mcrA-OTU contribution filter (< 50%)
  if (file.exists(picrust_combined_file)) {
    combined_pvals <- read.csv(picrust_combined_file, stringsAsFactors = FALSE)
    contrib_lookup <- setNames(combined_pvals$mean_percent_from_mcra, combined_pvals$pathway)
    sig_1e4$gene_contrib <- contrib_lookup[sig_1e4$pathway]
    sig_1e4$gene_contrib[is.na(sig_1e4$gene_contrib)] <- 0
    sig_filtered <- sig_1e4 %>% filter(gene_contrib < 0.50)
    stat("After mcrA-OTU contribution < 50% filter", nrow(sig_filtered))
    record("fig6_n_pathways_heatmap", nrow(sig_filtered))

    # List top pathways (positive and negative associations)
    sub_header("Top positively associated pathways (Fig 6)")
    pos_paths <- sig_filtered %>% filter(t_value > 0) %>% arrange(FDR)
    for (i in 1:min(10, nrow(pos_paths))) {
      desc <- ifelse(!is.na(pos_paths$description[i]) & pos_paths$description[i] != "",
                     pos_paths$description[i], pos_paths$pathway[i])
      cat(sprintf("    %s (FDR=%.2e, t=%.2f)\n", desc, pos_paths$FDR[i], pos_paths$t_value[i]))
    }

    sub_header("Top negatively associated pathways (Fig 6)")
    neg_paths <- sig_filtered %>% filter(t_value < 0) %>% arrange(FDR)
    for (i in 1:min(10, nrow(neg_paths))) {
      desc <- ifelse(!is.na(neg_paths$description[i]) & neg_paths$description[i] != "",
                     neg_paths$description[i], neg_paths$pathway[i])
      cat(sprintf("    %s (FDR=%.2e, t=%.2f)\n", desc, neg_paths$FDR[i], neg_paths$t_value[i]))
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
                  mp, best$pathway, best$FDR, best$t_value))
    } else {
      cat(sprintf("  NOT FOUND: '%s'\n", mp))
    }
  }
} else {
  cat("  [SKIPPED] No-mcrA-OTU pathway associations file not found.\n")
}


# ==============================================================================
# SECTION 7: FELLED BLACK OAK PROFILES (Figure 7)
# ==============================================================================

section_header("SECTION 7: FELLED BLACK OAK PROFILES (Figure 7)")

# Replicates data processing from 09_felled_oak_profiles.R
# Three panels: internal CH4, mcrA abundance, CH4 flux along stem height

# --- Load and process GC data ---
gc_file <- "data/raw/field_data/black_oak/quve_gc_data.csv"
o2_std_file <- "data/raw/field_data/black_oak/o2_standards.csv"
ghg_std_file <- "data/raw/field_data/black_oak/ghg_standards.csv"
mcra_bo_file <- "data/raw/ddpcr/black_oak_mcrA.csv"
flux_bo_file <- "data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv"

if (all(file.exists(gc_file, o2_std_file, ghg_std_file, mcra_bo_file, flux_bo_file))) {

  O2_standards_bo <- read.csv(o2_std_file, fileEncoding = "UTF-8-BOM")
  GHG_standards_bo <- read.csv(ghg_std_file, fileEncoding = "UTF-8-BOM")
  GC_data_bo <- read.csv(gc_file, fileEncoding = "UTF-8-BOM")

  colnames(GC_data_bo) <- make.names(colnames(GC_data_bo))
  colnames(O2_standards_bo) <- make.names(colnames(O2_standards_bo))
  colnames(GHG_standards_bo) <- make.names(colnames(GHG_standards_bo))

  GC_data_bo$O2.Area[which(GC_data_bo$Sample.Type == "N2")] <- NA

  # O2 calibration
  O2_data_bo <- O2_standards_bo %>%
    left_join(GC_data_bo, by = c("Sample" = "Sample.Type"))
  O2_data_bo$O2.Area <- as.numeric(O2_data_bo$O2.Area)
  O2_data_bo$O2_ppm <- as.numeric(gsub(",", "", O2_data_bo$X.O2...ppm.))

  low_range_o2_stds <- c("Outdoor Air 2", "Outdoor Air 3", "Outdoor Air 4",
                          "Oxygen Standard 1", "Oxygen Standard 2",
                          "Oxygen Standard 3", "Oxygen Standard 4")
  high_range_o2_stds <- c("Outdoor Air 5", "Outdoor Air 6", "Outdoor Air 7",
                           "Oxygen Standard 5")

  low_o2 <- O2_data_bo %>% filter(Sample %in% low_range_o2_stds)
  all_o2 <- bind_rows(low_o2, O2_data_bo %>% filter(Sample %in% high_range_o2_stds))

  O2_low_curve_bo <- lm(O2_ppm ~ O2.Area, data = low_o2)
  O2_high_curve_bo <- lm(O2_ppm ~ O2.Area, data = all_o2)
  max_O2_low_area_bo <- max(low_o2$O2.Area, na.rm = TRUE)

  # GHG calibration
  GHG_data_bo <- GHG_standards_bo %>%
    left_join(GC_data_bo, by = c("Sample" = "Sample.Type"))
  GHG_data_bo$CH4.Area <- as.numeric(GHG_data_bo$CH4.Area)
  GHG_data_bo$CH4_ppm <- as.numeric(gsub(",", "", GHG_data_bo$X.CH4...ppm.))
  GHG_data_bo <- GHG_data_bo %>%
    filter(!is.na(CH4.Area) & !is.na(CH4_ppm))

  low_range_ghg_stds <- c("N2", "SB1", "SB2", "SB3")
  low_ghg <- GHG_data_bo %>% filter(Sample %in% low_range_ghg_stds)
  CH4_low_curve_bo <- lm(CH4_ppm ~ CH4.Area, data = low_ghg)
  CH4_high_curve_bo <- lm(CH4_ppm ~ CH4.Area, data = GHG_data_bo)
  max_CH4_low_area_bo <- 1000

  # Apply calibration
  GC_data_bo <- GC_data_bo %>%
    mutate(
      O2_conc = if_else(
        O2.Area <= max_O2_low_area_bo,
        predict(O2_low_curve_bo, newdata = data.frame(O2.Area = O2.Area)),
        predict(O2_high_curve_bo, newdata = data.frame(O2.Area = O2.Area))
      ),
      CH4_conc = if_else(
        CH4.Area <= max_CH4_low_area_bo,
        predict(CH4_low_curve_bo, newdata = data.frame(CH4.Area = CH4.Area)),
        predict(CH4_high_curve_bo, newdata = data.frame(CH4.Area = CH4.Area))
      )
    )

  # Internal gas samples
  int_gas_bo <- GC_data_bo %>% filter(!is.na(Lab.ID), Tree.Tissue == "Trunk Gas")
  int_gas_bo$Tree.Height <- as.numeric(int_gas_bo$Tree.Height)

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
  mcra_bo <- readr::read_csv(mcra_bo_file, show_col_types = FALSE)
  conc_col_bo <- grep("Conc", colnames(mcra_bo), value = TRUE)[1]
  height_col_bo <- grep("Height", colnames(mcra_bo), value = TRUE)[1]
  mcra_bo <- mcra_bo %>% rename(Height_cm = !!height_col_bo, Conc = !!conc_col_bo)
  mcra_wood <- mcra_bo %>% filter(Component %in% c("Heartwood", "Sapwood"))

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
  flux_bo <- read.csv(flux_bo_file)
  flux_bo$Height_m <- suppressWarnings(as.numeric(flux_bo$Height_m))
  flux_bo <- flux_bo %>% filter(!is.na(Height_m), !is.na(CH4_best.flux))

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

} else {
  cat("  [SKIPPED] Black oak data files not found.\n")
  missing <- c(gc_file, o2_std_file, ghg_std_file, mcra_bo_file, flux_bo_file)
  cat("  Missing:", paste(missing[!file.exists(missing)], collapse = ", "), "\n")
}


# ==============================================================================
# SECTION 8: INDIVIDUAL GENE-FLUX (Figure S7, S10)
# ==============================================================================

section_header("SECTION 8: INDIVIDUAL GENE-FLUX RELATIONSHIPS (Figs S7, S10)")

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
      species = species_mapping[species_id],
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
  m_sp_mmox <- lm(CH4_best.flux_125cm ~ species + log_mmox,
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
# SECTION 12: SUPPLEMENTARY FIGURE STATISTICS
# ==============================================================================

section_header("SECTION 12: SUPPLEMENTARY FIGURE STATISTICS")

# --- Figure S1: Methanotroph relative abundance (from 08c) ---
sub_header("Figure S1: Methanotroph relative abundances by compartment")
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
  sp_mt_means <- hw_samps %>%
    mutate(mt_pct = colSums(otu_df[intersect(all_mt_asvs, rownames(otu_df)),
                                    rownames(hw_samps), drop = FALSE])) %>%
    group_by(species.x) %>%
    summarise(mean_mt = mean(mt_pct), n = n(), .groups = "drop") %>%
    arrange(desc(mean_mt))
  for (i in 1:min(5, nrow(sp_mt_means))) {
    cat(sprintf("    %s: %.3f%% (n=%d)\n",
                species_mapping[sp_mt_means$species.x[i]],
                sp_mt_means$mean_mt[i], sp_mt_means$n[i]))
  }
} else {
  cat("  [Requires phyloseq objects from Section 5]\n")
}

# --- Figure S2: FAPROTAX (already in Section 6, just cross-reference) ---
sub_header("Figure S2: FAPROTAX (see Section 6 above)")

# --- Figure S3: PICRUSt inner vs outer + methanotroph associations ---
sub_header("Figure S3a: PICRUSt inner vs outer core pathways")
picrust_mcra_all <- "data/processed/molecular/picrust/pathway_associations_mcra_all.csv"
if (file.exists(picrust_mcra_all)) {
  pa_all <- read.csv(picrust_mcra_all, stringsAsFactors = FALSE)
  sig_all <- pa_all %>% filter(!is.na(FDR), FDR < 0.01)
  stat("Total pathways tested (all OTUs)", nrow(pa_all))
  stat("Significant at FDR < 0.01 (all OTUs)", nrow(sig_all))

  # Top positive (heartwood-enriched) and negative (sapwood-enriched)
  pos_all <- sig_all %>% filter(t_value > 0) %>% arrange(FDR)
  neg_all <- sig_all %>% filter(t_value < 0) %>% arrange(FDR)
  cat(sprintf("  Positively associated: %d, Negatively associated: %d\n",
              nrow(pos_all), nrow(neg_all)))
}

sub_header("Figure S3b: pmoA pathway associations")
picrust_pmoa_no <- "data/processed/molecular/picrust/pathway_associations_pmoa_no_pmoa_otus.csv"
if (file.exists(picrust_pmoa_no)) {
  pa_pmoa <- read.csv(picrust_pmoa_no, stringsAsFactors = FALSE)
  sig_pmoa_no <- pa_pmoa %>% filter(!is.na(FDR), FDR < 0.01)
  stat("Significant pmoA pathways (no-pmoA OTUs, FDR < 0.01)", nrow(sig_pmoa_no))
  pos_pmoa <- sig_pmoa_no %>% filter(t_value > 0)
  neg_pmoa <- sig_pmoa_no %>% filter(t_value < 0)
  cat(sprintf("  Positive: %d, Negative: %d\n", nrow(pos_pmoa), nrow(neg_pmoa)))
} else {
  pmoa_alt <- "data/processed/molecular/picrust/pathway_associations_pmoa.csv"
  if (file.exists(pmoa_alt)) {
    pa_pmoa <- read.csv(pmoa_alt, stringsAsFactors = FALSE)
    sig_pmoa_full <- pa_pmoa %>% filter(!is.na(FDR), FDR < 0.01)
    stat("Significant pmoA pathways (FDR < 0.01)", nrow(sig_pmoa_full))
  }
}

# --- Figure S4: PICRUSt mcrA all-OTU heatmap (expanded version of Fig 6) ---
sub_header("Figure S4: PICRUSt mcrA (all OTUs) heatmap")
if (file.exists(picrust_mcra_all)) {
  cat(sprintf("  See S3a counts above. Full FDR<0.01 set: %d pathways\n", nrow(sig_all)))
}

# --- Figure S5: PICRUSt pmoA heatmap ---
sub_header("Figure S5: PICRUSt pmoA heatmap")
cat("  See S3b counts above.\n")

# --- Figure S6: Internal gas beeswarm by species ---
sub_header("Figure S6: Internal gas concentrations by species")
gas_data <- ymf2021 %>%
  filter(!is.na(CH4_concentration), !is.na(species_id)) %>%
  mutate(species_latin = species_mapping[species_id]) %>%
  filter(!is.na(species_latin))

stat("Total trees with internal gas", nrow(gas_data))
stat("Number of species", n_distinct(gas_data$species_latin))

# Overall CH4 stats
stat("Overall CH4 mean", round(mean(gas_data$CH4_concentration), 0), "ppm")
stat("Overall CH4 median", round(median(gas_data$CH4_concentration), 0), "ppm")
stat("Overall CH4 range", paste(round(min(gas_data$CH4_concentration), 0), "to",
                                 round(max(gas_data$CH4_concentration), 0)), "ppm")

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

# --- Figure S7: Internal gas profiles (4 panels) ---
sub_header("Figure S7: Internal gas correlation panels")

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

# --- Figure S8: δ13CH4 isotope analysis ---
sub_header("Figure S8: d13CH4 isotope analysis")
pic_files <- c("data/raw/internal_gas/picarro/20251128_211030_results.csv",
               "data/raw/internal_gas/picarro/20251128_213226_results.csv",
               "data/raw/internal_gas/picarro/20251128_215521_results.csv")
if (all(file.exists(pic_files))) {
  pic_all <- bind_rows(lapply(pic_files, function(f) readr::read_csv(f, show_col_types = FALSE)))

  # Map species
  ddpcr_meta <- read.csv("data/raw/ddpcr/ddPCR_meta_all_data.csv")
  sp_map_iso <- ddpcr_meta %>% distinct(seq_id, species) %>% rename(SampleName = seq_id)
  pic_all <- pic_all %>% left_join(sp_map_iso, by = "SampleName") %>%
    mutate(species = ifelse(grepl("^Amb", SampleName), "Atmosphere", species)) %>%
    filter(!is.na(species))

  # Extract key columns
  iso_s8 <- pic_all %>%
    dplyr::select(SampleName, species,
                  d13CH4 = HR_Delta_iCH4_Raw_mean,
                  ch4_ppm = HR_12CH4_dry_mean) %>%
    filter(!is.na(d13CH4), ch4_ppm >= 1.5) %>%
    filter(!(species == "Atmosphere" & ch4_ppm > 5))

  internal_iso <- iso_s8 %>% filter(species != "Atmosphere")
  atm_iso <- iso_s8 %>% filter(species == "Atmosphere")

  stat("Total isotope samples", nrow(iso_s8))
  stat("Internal tree samples", nrow(internal_iso))
  stat("Atmosphere samples", nrow(atm_iso))
  stat("Number of species", n_distinct(internal_iso$species))

  stat("Internal d13CH4 mean", round(mean(internal_iso$d13CH4), 1), "permil VPDB")
  stat("Internal d13CH4 range", paste(round(min(internal_iso$d13CH4), 1), "to",
                                       round(max(internal_iso$d13CH4), 1)), "permil VPDB")
  stat("Internal d13CH4 SD", round(sd(internal_iso$d13CH4), 1), "permil")
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
    stat("Keeling intercept", round(coef(keeling_lm)[1], 1), "permil")
    stat("Keeling R2", round(summary(keeling_lm)$r.squared, 3))
  }
} else {
  cat("  [SKIPPED] Picarro isotope files not found.\n")
}

# --- Figure S9: RF model diagnostics (expanded) ---
sub_header("Figure S9: RF model diagnostics (expanded)")
if (exists("TreeRF") && exists("SoilRF")) {
  # CCC requires DescTools; compute manually if not available
  tree_obs <- tree_train_complete$stem_flux_corrected * 1000
  tree_pred <- tree_train_complete$pred_flux * 1000
  soil_obs <- soil_train_complete$soil_flux_umol_m2_s * 1000
  soil_pred <- soil_train_complete$pred_flux * 1000

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
    cat(sprintf("    %d. %s (%.1f)\n", i, tree_imp_df$feature[i], tree_imp_df$imp[i]))
  }

  sub_header("  Soil RF top features")
  soil_imp_df <- data.frame(feature = names(importance(SoilRF)),
                              imp = as.numeric(importance(SoilRF))) %>%
    arrange(desc(imp))
  for (i in 1:min(8, nrow(soil_imp_df))) {
    cat(sprintf("    %d. %s (%.1f)\n", i, soil_imp_df$feature[i], soil_imp_df$imp[i]))
  }
}

# --- Figure S10: already covered in Sections 8 + 11 ---
sub_header("Figure S10: Scale-dependent genes (see Sections 8-9)")
cat("  Individual gene-flux R2: see Section 8\n")
cat("  Species-level correlations: see Section 9\n")

# --- Figure S11: Tree radial mcrA cross-sections ---
sub_header("Figure S11: Tree radial mcrA cross-sections")
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
  mutate(species_latin = species_mapping[species_id],
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

# --- Figure S12: pmoA vs mmoX patterns (already in Section 11) ---
sub_header("Figure S12: pmoA/mmoX patterns (see Section 11)")
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

# --- Figure S13: mcrA vs methanotroph independence (already in Section 9) ---
sub_header("Figure S13: mcrA vs methanotroph (see Section 9)")
cat("  Tree-level and species-level correlations reported above.\n")

# --- Figure S14: Black oak methanome heatmap ---
sub_header("Figure S14: Black oak tissue methanome heatmap")
# This is from 10_black_oak_methanome_heatmap.R — qualitative, but key result:
# Methanobacteriaceae only in heartwood, no methanogens in soil
cat("  Key qualitative results (from manuscript text):\n")
cat("  - Only Methanobacterium (Methanobacteriaceae) detected as methanogen\n")
cat("  - Found exclusively in heartwood and sapwood\n")
cat("  - No methanogens in bark, branches, foliage, or soil\n")
cat("  - Methanotrophs (Methylocella) in heartwood, sapwood, bark, foliage\n")
cat("  - Methyloferula in sapwood, bark, leaf litter\n")
cat("  - Methylorosula in heartwood, sapwood, leaf litter\n")

# --- Figure S15: Taxonomy x pmoA heatmap ---
sub_header("Figure S15: Family-level 16S x pmoA associations")
picrust_pmoa_full <- "data/processed/molecular/picrust/pathway_associations_pmoa.csv"
if (file.exists(picrust_pmoa_full)) {
  pa_pmoa_f <- read.csv(picrust_pmoa_full, stringsAsFactors = FALSE)
  sig_pmoa_f <- pa_pmoa_f %>% filter(!is.na(FDR), FDR < 0.01)
  stat("Total pathways tested (pmoA)", nrow(pa_pmoa_f))
  stat("Significant at FDR < 0.01 (pmoA)", nrow(sig_pmoa_f))
}

# Also check taxonomy heatmap for S3
taxonomy_mcra_file <- "data/processed/molecular/picrust/pathway_associations_mcra_no_mcra_otus.csv"
sub_header("Figure S3: Family-level 16S x mcrA associations")
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
