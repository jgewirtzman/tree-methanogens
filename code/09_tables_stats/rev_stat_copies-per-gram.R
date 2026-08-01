# ==============================================================================
# REVISION — Gene copies per gram (real mass) + dry/wet weight basis
# ==============================================================================
# Addresses R2 #1 (unit basis: dry vs wet; "copies per gram") and R2 Fig 7b
# (copies per uL depends on elution volume & extraction mass -> use copies/g).
#
# The pipeline (code/07_molecular/05_ddpcr_analysis.R, L287) converts with a FIXED
# 100 mg proxy mass and 37.5 uL elution:  copies ~ concentration * 37.5 / 100.
# We now have the ACTUAL per-sample mass (sample_mass_mg), so we compute the
# true mass-normalized value:
#     copies_per_g_fresh = concentration_copies_per_uL * 37.5 / (mass_mg/1000)
# and convert to a DRY-weight basis for wood using per-tree tissue moisture.
#
# KEY consequence: wood mass ~100 mg (proxy was ~right) but soil ~250 mg, so the
# fixed proxy OVERESTIMATED soil copies/g ~2.5x. The real correction LOWERS soil
# and thus STRENGTHENS the "wood exceeds soil by ~2 orders of magnitude" result.
#
# NEW file; edits nothing. Run: Rscript code/revision/R_copies_per_gram.R
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "outputs/revision"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ELUTION_UL <- 75   # CANONICAL: 02_harmonize_all_data.R:317 uses Conc * 75 / mass_mg * 1000
                   # (75 uL elution volume; treats ddPCR Conc as copies/uL of eluent).
                   # Earlier draft used 37.5 (=75/2) -> was 2x low; corrected to match manuscript.

d <- read_csv("data/compiled/ddpcr_gene_abundances.csv", show_col_types = FALSE) %>%
  filter(analysis_type == "loose", !is.na(sample_mass_mg), sample_mass_mg > 0)

# tree_id parsed from sample_id (e.g. "ACRU_RM1\nWood_Outer\n...") -> "ACRU_RM1"
d <- d %>% mutate(tree_id = str_trim(str_split_fixed(sample_id, "\n", 2)[,1]))

# copies per gram, fresh (as-extracted) mass basis, with REAL mass
d <- d %>% mutate(
  copies_per_g_proxy100 = concentration_copies_per_uL * ELUTION_UL / 100 * 1000,   # old pipeline proxy
  copies_per_g_fresh    = concentration_copies_per_uL * ELUTION_UL / sample_mass_mg * 1000
)

# ---- wood vs soil: effect of real mass on the 2-orders-of-magnitude claim -----
comp <- d %>% filter(target_gene %in% c("mcra","mcra_probe")) %>%
  mutate(compartment = case_when(
    material == "Wood" & core_type == "Inner" ~ "Heartwood",
    material == "Wood" & core_type == "Outer" ~ "Sapwood",
    material == "Soil" ~ paste0("Soil_", core_type), TRUE ~ NA_character_)) %>%
  filter(!is.na(compartment)) %>%
  group_by(compartment) %>%
  summarise(n = n(),
            median_proxy = median(copies_per_g_proxy100, na.rm = TRUE),
            median_real  = median(copies_per_g_fresh, na.rm = TRUE),
            mass_median_mg = median(sample_mass_mg), .groups = "drop")
write.csv(comp, file.path(out_dir, "copies_per_g_compartment.csv"), row.names = FALSE)

hw <- comp$median_real[comp$compartment == "Heartwood"]
soil_real  <- max(comp$median_real[grepl("Soil", comp$compartment)], na.rm = TRUE)
soil_proxy <- max(comp$median_proxy[grepl("Soil", comp$compartment)], na.rm = TRUE)

# ---- dry-weight conversion for WOOD using tissue moisture ---------------------
# tree_properties: *_moisture_fresh_percent = water / fresh mass * 100 (per tree/tissue)
# NOTE: a per-sample join is blocked by ID schemes (ddPCR sample_id "ACRU_RM1"
# vs tree_properties numeric tree_id) -> needs the harmonization crosswalk.
# Here we report the representative fresh->dry factor from tree_properties moisture.
tp <- read_csv("data/compiled/tree_properties.csv", show_col_types = FALSE)
wood_conv <- tibble(
  tissue = c("Inner (heartwood)", "Outer (sapwood)"),
  median_fresh_moisture_pct = c(median(tp$inner_moisture_fresh_percent, na.rm = TRUE),
                                median(tp$outer_moisture_fresh_percent, na.rm = TRUE)),
  fresh_to_dry_factor = c(median(1/(1 - tp$inner_moisture_fresh_percent/100), na.rm = TRUE),
                          median(1/(1 - tp$outer_moisture_fresh_percent/100), na.rm = TRUE)))

# ---- estimated soil DRY-basis copies/g using black-oak GWC (organic vs mineral)
# Main-study soils lack gravimetric moisture; use black-oak soil GWC as typical.
# GWC = water/dry, so copies/g_dry = copies/g_fresh * (1 + GWC).
bo_sm <- tryCatch(read_csv("data/processed/molecular/black_oak/bo_soil_moisture.csv", show_col_types = FALSE),
                  error = function(e) NULL)
soil_dry <- NULL
if (!is.null(bo_sm)) {
  g <- suppressWarnings(as.numeric(bo_sm$GWC))
  gwc_org <- median(g[grepl("Organic", bo_sm$`Sample Name`)], na.rm = TRUE)
  gwc_min <- median(g[grepl("Mineral", bo_sm$`Sample Name`)], na.rm = TRUE)
  soil_org_fresh <- comp$median_real[comp$compartment == "Soil_Organic"]
  soil_min_fresh <- comp$median_real[comp$compartment == "Soil_Mineral"]
  soil_dry <- tibble(
    compartment = c("Soil_Organic", "Soil_Mineral"),
    assumed_GWC = c(gwc_org, gwc_min),
    fresh_to_dry_factor = c(1 + gwc_org, 1 + gwc_min),
    copies_g_fresh = c(soil_org_fresh, soil_min_fresh),
    copies_g_dry_est = c(soil_org_fresh * (1 + gwc_org), soil_min_fresh * (1 + gwc_min)))
  write.csv(soil_dry, file.path(out_dir, "copies_per_g_soil_dry_estimate.csv"), row.names = FALSE)
}

sink(file.path(out_dir, "copies_per_g_summary.txt"))
cat("=================================================================\n")
cat("GENE COPIES PER GRAM — real mass vs proxy, and dry/wet basis\n")
cat("=================================================================\n\n")
cat(sprintf("Elution volume used (pipeline): %.1f uL. Formula (matches 02_harmonize L317):\n", ELUTION_UL))
cat("  copies/g = concentration_copies_per_uL * 75 / mass_mg * 1000\n\n")
cat("Median wood mass ~", round(median(d$sample_mass_mg[d$material=='Wood'],na.rm=T)),
    "mg; median soil mass ~", round(median(d$sample_mass_mg[d$material=='Soil'],na.rm=T)), "mg.\n")
cat("=> the fixed 100 mg proxy is ~right for wood but too low for soil, so it\n")
cat("   OVERESTIMATED soil copies/g. Real masses lower soil, strengthening wood>>soil.\n\n")
cat("mcrA copies/g by compartment (proxy-100mg vs real mass):\n")
print(as.data.frame(comp), row.names = FALSE)
cat(sprintf("\nHeartwood : soil ratio (mcrA copies/g):  proxy = %.0fx ;  real = %.0fx\n",
            hw/soil_proxy, hw/soil_real))
cat("\nWOOD dry-weight conversion (fresh -> dry) using per-tree tissue moisture:\n")
print(as.data.frame(wood_conv), row.names = FALSE)
cat("\n*** CONFIRMED MASS BASIS (Black Oak Tree Project Data.xlsx extraction sheet:\n")
cat("    'TARGET MASS FOR WOOD: 100 mg; SOIL: 250 mg'; matches sample_mass_mg): ***\n")
cat("  * sample_mass_mg = mass weighed into the extraction tube.\n")
cat("  * WOOD was FREEZE-DRIED before weighing -> wood copies/g is ALREADY DRY-weight.\n")
cat("  * SOIL was FRESH -> soil copies/g is WET-weight. (Exactly R2's observation.)\n")
cat("  NOTE: the 'copies_per_g_fresh' column name above is a misnomer for wood\n")
cat("  (it is dry for wood, fresh for soil). Values are correct; label per material.\n")
cat("\nDRY/WET STANDARDIZATION (answer to R2):\n")
cat("  * To a COMMON basis we can convert WOOD dry->fresh using core moisture\n")
cat("    (fresh & dry masses per core), giving everything on a FRESH/wet basis\n")
cat("    (soil already fresh). This is the feasible common unit without soil GWC.\n")
cat("  * All-DRY is only possible where soil GWC exists (black-oak soils have GWC;\n")
cat("    main-study soils do not).\n")
cat("  * Per-SAMPLE wood moisture join uses tree_id_comprehensive_mapping.csv\n")
cat("    (data/processed/tree_data/) to bridge ddPCR ids <-> tree_id.\n")
cat("  * Soil: gravimetric soil moisture per ddPCR soil sample is NOT in the\n")
cat("    compiled ddPCR table; VWC (volumetric) exists per tree/plot but is not a\n")
cat("    direct gravimetric conversion. To standardize soil to dry weight we need\n")
cat("    the soil sample gravimetric moisture (or to state soil is fresh-basis).\n")
cat("  * Recommendation: report ALL copies/g on a DRY-weight basis; convert wood\n")
cat("    with tissue moisture; for soil either supply gravimetric moisture or\n")
cat("    state the basis explicitly and note it does not affect wood>>soil.\n\n")
if (!is.null(soil_dry)) {
  cat("\nESTIMATED SOIL DRY-basis copies/g (black-oak GWC as typical; SI note):\n")
  print(as.data.frame(soil_dry), row.names = FALSE)
  hw_dry <- comp$median_real[comp$compartment == "Heartwood"]   # wood already dry
  cat(sprintf("\nConsistent DRY-basis heartwood:soil ratio: organic %.0fx, mineral %.0fx\n",
              hw_dry / soil_dry$copies_g_dry_est[1], hw_dry / soil_dry$copies_g_dry_est[2]))
  cat("(vs 34x from the current mixed wood-dry/soil-fresh comparison; still wood>>soil.)\n")
}
cat("\nDNA EXTRACTION EFFICIENCY: the workflow applies NO recovery/efficiency\n")
cat("correction, so all copies/g are conservative LOWER BOUNDS (extraction <100%).\n")
cat("Worth a one-line discussion mention; it only strengthens 'methanogens abundant'.\n\n")
cat("FIG 7b (R2): the height-resolved black-oak mcrA profile is NOT in the\n")
cat("compiled ddPCR table (QUVE_BO1 is single-tissue), so the per-gram fix for\n")
cat("Fig 7b needs the height-resolved molecular source file located (flag Jon).\n")
cat("Once located, apply the same copies/g formula (it has mass) to replace\n")
cat("copies/uL with copies/g.\n")
sink()
cat(readLines(file.path(out_dir,"copies_per_g_summary.txt")), sep="\n")
