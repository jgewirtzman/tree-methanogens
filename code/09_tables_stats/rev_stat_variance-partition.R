# ==============================================================================
# REVISION — Variance partition review: reconcile 82.9% vs 91.3% + mixed-model ICC
# ==============================================================================
# Jon's question: is the new (~91%) partition different from the ~82.9% already
# in the manuscript? Which is correct?
#
# ANSWER (see printed output):
#   * Manuscript 82.9% comes from code/08_figures/04_variance_partition.R
#     "Method 2": an OLS hierarchical partition of a model that includes
#     ENVIRONMENT (DBH, air/soil temp, VWC) + SPECIES + their INTERACTIONS.
#     Unexplained = 1 - R2(interaction model) = 82.9%; Species-UNIQUE = 5.3%.
#     It is NOT a mixed / random-effects model despite the "mixed" shorthand.
#   * The ~91% earlier was a simpler lm(CH4 ~ species) residual (species-ONLY),
#     i.e. 1 - R2(species-only). Different quantity, not a contradiction.
#   * The cleanest single number for "share of individual flux variance that is
#     between species" is the mixed-model ICC = var(species) / total var.
#
# This script reproduces the manuscript Method-2 numbers AND reports the ICC on
# raw and arcsinh(pseudo-log) flux. NEW file; edits nothing.
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(car); library(lme4) })
out_dir <- "outputs/revision"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ymf2023 <- read.csv("data/processed/flux/methanogen_tree_flux_complete_dataset.csv")
ymf2021 <- read.csv("data/processed/integrated/merged_tree_dataset_final.csv")
species_mapping <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",FAGR="Fagus grandifolia",FRAM="Fraxinus americana",
  PIST="Pinus strobus",QURU="Quercus rubra",TSCA="Tsuga canadensis",CAOV="Carya ovata",
  KALA="Kalmia latifolia",PRSE="Prunus serotina",QUAL="Quercus alba",QUVE="Quercus velutina",
  SAAL="Sassafras albidum")
psl <- function(x) asinh(x / 0.2) / log(10)   # arcsinh, manuscript pseudolog10_individual

# --- rebuild combined_data exactly as 04_variance_partition.R ------------------
d23 <- ymf2023 %>% dplyr::select(Species.Code, DBH=DBH..cm., Air_temp=air_temp_C,
        Stem_temp=stem_temp_C, Soil_temp=Soil.Temp....C., VWC=vwc_mean, CH4_flux=CH4_best.flux) %>%
  mutate(Year="2023", Species_Latin=species_mapping[Species.Code]) %>% drop_na(CH4_flux)
d21 <- ymf2021 %>% dplyr::select(Species.Code=species_id, DBH=dbh, Air_temp=Temp_Air_125cm,
        Soil_temp=SoilTemp_mean, VWC=VWC_mean, CH4_flux=CH4_best.flux_125cm) %>%
  mutate(Year="2021", Species_Latin=species_mapping[Species.Code], Stem_temp=NA) %>%
  drop_na(CH4_flux) %>% filter(!is.nan(CH4_flux))
combined <- bind_rows(d23, d21) %>% filter(!is.na(Species_Latin)) %>%
  group_by(Species_Latin) %>% filter(n() > 3) %>% ungroup()
for (v in c("DBH","Air_temp","Stem_temp","Soil_temp","VWC"))
  combined[[paste0(v,"_std")]] <- as.numeric(scale(combined[[v]]))

# --- Manuscript Method 2 (OLS hierarchical partition, incl. interactions) ------
m_env  <- lm(CH4_flux ~ DBH_std+Air_temp_std+Soil_temp_std+VWC_std, combined)
m_sp   <- lm(CH4_flux ~ Species_Latin, combined)
m_full <- lm(CH4_flux ~ DBH_std+Air_temp_std+Soil_temp_std+VWC_std+Species_Latin, combined)
m_int  <- lm(CH4_flux ~ (DBH_std+Air_temp_std+Soil_temp_std+VWC_std)*Species_Latin, combined)
r2 <- c(env=summary(m_env)$r.squared, sp=summary(m_sp)$r.squared,
        full=summary(m_full)$r.squared, int=summary(m_int)$r.squared)
r2 <- unname(r2); names(r2) <- c("env","sp","full","int")
method2 <- c(
  Environment_unique     = unname(max(0, r2["full"]-r2["sp"])*100),
  Species_unique         = unname(max(0, r2["full"]-r2["env"])*100),
  EnvSpecies_interaction = unname(max(0, r2["int"]-r2["full"])*100),
  Unexplained            = unname((1-r2["int"])*100))

# --- Mixed-model ICC (the clean "between-species share") -----------------------
icc <- function(y) {
  m <- lmer(y ~ 1 + (1|Species_Latin), data = combined,
            control = lmerControl(calc.derivs = FALSE))
  vc <- as.data.frame(VarCorr(m)); v_sp <- vc$vcov[vc$grp=="Species_Latin"]; v_res <- vc$vcov[vc$grp=="Residual"]
  c(between_species = v_sp, within = v_res, ICC_pct = 100*v_sp/(v_sp+v_res))
}
icc_raw <- icc(combined$CH4_flux)
icc_psl <- icc(psl(combined$CH4_flux))

# --- also the simple species-only lm residual (the earlier ~91%) ---------------
species_only_unexpl <- (1 - summary(m_sp)$r.squared) * 100

sink(file.path(out_dir, "variance_partition_review.txt"))
cat("=====================================================================\n")
cat("VARIANCE PARTITION REVIEW — what each number means, and what's correct\n")
cat("=====================================================================\n\n")
cat(sprintf("combined n = %d (2021=%d, 2023=%d), species = %d\n\n",
            nrow(combined), sum(combined$Year=="2021"), sum(combined$Year=="2023"),
            n_distinct(combined$Species_Latin)))
cat("Model R-squared values:\n")
cat(sprintf("  Environment only      : %.1f%%\n", r2["env"]*100))
cat(sprintf("  Species only          : %.1f%%\n", r2["sp"]*100))
cat(sprintf("  Full additive         : %.1f%%\n", r2["full"]*100))
cat(sprintf("  With interactions     : %.1f%%\n\n", r2["int"]*100))
cat("(1) MANUSCRIPT 'Method 2' (OLS hierarchical, env+species+interaction):\n")
cat(sprintf("      Environment (unique)      : %.1f%%\n", method2["Environment_unique"]))
cat(sprintf("      Species (unique)          : %.1f%%   <- the manuscript's 5.3%%\n", method2["Species_unique"]))
cat(sprintf("      Env x Species interaction : %.1f%%\n", method2["EnvSpecies_interaction"]))
cat(sprintf("      Unexplained               : %.1f%%   <- the manuscript's 82.9%%\n\n", method2["Unexplained"]))
cat("(2) SPECIES-ONLY lm residual (the earlier ~91%):\n")
cat(sprintf("      1 - R2(species-only)      : %.1f%%\n", species_only_unexpl))
cat("      -> Not a contradiction: it excludes environment & interaction terms,\n")
cat("         so more variance is left 'unexplained' than in Method 2.\n\n")
cat("(3) MIXED-MODEL ICC  (clean between-species share of individual flux var):\n")
cat(sprintf("      raw flux      : ICC = %.1f%%  (between=%.4g, within=%.4g)\n",
            icc_raw["ICC_pct"], icc_raw["between_species"], icc_raw["within"]))
cat(sprintf("      arcsinh flux  : ICC = %.1f%%  (between=%.4g, within=%.4g)\n\n",
            icc_psl["ICC_pct"], icc_psl["between_species"], icc_psl["within"]))
cat("WHAT'S CORRECT / RECOMMENDATION:\n")
cat("  * All three are internally correct; they answer different questions.\n")
cat("  * Keep the manuscript's Method-2 numbers (5.3% species-unique, 82.9%\n")
cat("    unexplained) as the published partition, but call it what it is: an\n")
cat("    OLS hierarchical partition, NOT a mixed model. Caveat: the interaction\n")
cat("    R2 is inflated by many species x env parameters, so 'unexplained=82.9%'\n")
cat("    is a slight lower bound on the true residual.\n")
cat("  * For the aggregation argument, the ICC is the cleaner statistic: only\n")
cat(sprintf("    ~%.0f%% (raw) / ~%.0f%% (arcsinh) of INDIVIDUAL flux variance is between\n",
            icc_raw["ICC_pct"], icc_psl["ICC_pct"]))
cat("    species; the large within-species remainder is exactly the point-\n")
cat("    sampling heterogeneity that motivates species-level aggregation.\n")
sink()
cat(readLines(file.path(out_dir,"variance_partition_review.txt")), sep="\n")
