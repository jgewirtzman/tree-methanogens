# ==============================================================================
# 06_copies-per-gram.R -- gene copies per gram of sample (Fig 7b basis)
# ------------------------------------------------------------------------------
# Stage D; run from the repository root. Copies per uL depend on the elution volume
# and the extracted mass, so abundances are reported per gram, from each sample's own
# extraction mass (samples.extraction_mass_mg) and code/lib/ddpcr_constants.R:
#     copies g-1 = Conc x (25/2.5) x elution (75 wood, 100 soil) x (800/250) / mass (g)
# Also reported: the same with a fixed 100 mg proxy mass (soil extractions were
# ~250 mg, so the proxy overstates soil copies per g ~2.5x), and wood on a dry basis
# from each tree's tissue moisture.
# Outputs: outputs/data/copies_per_g_compartment.csv, outputs/audit/copies_per_g_summary.txt
# ==============================================================================
source("code/lib/outputs.R")

suppressPackageStartupMessages({ library(tidyverse) })
source("code/lib/ddpcr_constants.R")   # copies g-1 = Conc x (25/2.5) x elution (75 wood, 100 soil) x (800/250) / mass (g), as code/tables/21_ddpcr.R

# per DNA sample (data/tables/ddpcr.csv + samples.csv via dna_samples(); mean over the
# sample's runs that passed QC): mcrA probe assay copies per uL extract and per g
source("code/lib/tables_read.R")
d <- dna_samples("2021_survey") %>%
  filter(!is.na(extraction_mass_mg), extraction_mass_mg > 0, !is.na(mcra_probe_ul)) %>%
  mutate(sample_mass_mg = extraction_mass_mg,
         copies_per_g_proxy100 = extract_copies_per_g(mcra_probe_ul, 100, material),   # 100 mg proxy mass
         copies_per_g_fresh    = mcra_probe_g)                                         # real mass (ddpcr.copies_per_g)

# ---- wood vs soil: effect of real mass on the 2-orders-of-magnitude claim -----
comp <- d %>%   # probe assay only (EvaGreen mcrA gives false positives in wood; Fig. S4)
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
write.csv(comp, out_path("copies_per_g_compartment.csv"), row.names = FALSE)

hw <- comp$median_real[comp$compartment == "Heartwood"]
soil_real  <- max(comp$median_real[grepl("Soil", comp$compartment)], na.rm = TRUE)
soil_proxy <- max(comp$median_proxy[grepl("Soil", comp$compartment)], na.rm = TRUE)

# ---- dry-weight conversion for WOOD using tissue moisture ---------------------
# wood_cores.csv moisture_fresh_pct = water / fresh mass * 100, as the per-tree mean of each
# radial position (survey_2021_trees(): inner_/outer_moisture_fresh_percent).
# Here we report the representative fresh->dry factor (07_mass-basis-sensitivity.R joins per sample).
tp <- survey_2021_trees()
wood_conv <- tibble(
  tissue = c("Inner (heartwood)", "Outer (sapwood)"),
  median_fresh_moisture_pct = c(median(tp$inner_moisture_fresh_percent, na.rm = TRUE),
                                median(tp$outer_moisture_fresh_percent, na.rm = TRUE)),
  fresh_to_dry_factor = c(median(1/(1 - tp$inner_moisture_fresh_percent/100), na.rm = TRUE),
                          median(1/(1 - tp$outer_moisture_fresh_percent/100), na.rm = TRUE)))

# ---- estimated soil DRY-basis copies/g using black-oak GWC (organic vs mineral)
# Main-study soils lack gravimetric moisture; use black-oak soil GWC as typical (a labelled
# estimate). The felled-oak soil GWC sheet is not in data/tables/, so it is read raw.
# GWC = water/dry, so copies/g_dry = copies/g_fresh * (1 + GWC).
bo_sm <- read_csv("data/raw/field_data/black_oak/bo_soil_moisture.csv", show_col_types = FALSE)
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
write.csv(soil_dry, out_path("copies_per_g_soil_dry_estimate.csv"), row.names = FALSE)


sink(out_path("copies_per_g_summary.txt"))
cat("=================================================================\n")
cat("GENE COPIES PER GRAM — real mass vs proxy, and dry/wet basis\n")
cat("=================================================================\n\n")
cat(sprintf("Elution %.0f uL (wood) / %.0f uL (soil); template %.1f uL in a %.0f uL reaction. Formula (code/lib/ddpcr_constants.R):\n", DDPCR_ELUTION_UL[["Wood"]], DDPCR_ELUTION_UL[["Soil"]], DDPCR_TEMPLATE_UL, DDPCR_REACTION_UL))
cat("  copies/g = concentration_copies_per_uL * (25/2.5) * elution * (800/250) / mass_mg * 1000\n\n")
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
cat("  * SOIL was FRESH -> soil copies/g is WET-weight.\n")
cat("  NOTE: the 'copies_per_g_fresh' column name above is a misnomer for wood\n")
cat("  (it is dry for wood, fresh for soil). Values are correct; label per material.\n")
cat("\nDRY/WET STANDARDIZATION:\n")
cat("  * To a COMMON basis we can convert WOOD dry->fresh using core moisture\n")
cat("    (fresh & dry masses per core), giving everything on a FRESH/wet basis\n")
cat("    (soil already fresh). This is the feasible common unit without soil GWC.\n")
cat("  * All-DRY is only possible where soil GWC exists (black-oak soils have GWC;\n")
cat("    main-study soils do not).\n")
cat("  * Per-SAMPLE wood moisture join: samples.csv tree_id <-> wood_cores.csv tree_id\n")
cat("    (07_mass-basis-sensitivity.R).\n")
cat("  * Soil: gravimetric soil moisture per ddPCR soil sample is NOT in the\n")
cat("    ddPCR or samples tables; VWC (volumetric) exists per tree/plot but is not a\n")
cat("    direct gravimetric conversion. To standardize soil to dry weight we need\n")
cat("    the soil sample gravimetric moisture (or to state soil is fresh-basis).\n")
cat("  * Recommendation: report ALL copies/g on a DRY-weight basis; convert wood\n")
cat("    with tissue moisture; for soil either supply gravimetric moisture or\n")
cat("    state the basis explicitly and note it does not affect wood>>soil.\n\n")
cat("\nESTIMATED SOIL DRY-basis copies/g (black-oak GWC as typical; SI note):\n")
print(as.data.frame(soil_dry), row.names = FALSE)
hw_dry <- comp$median_real[comp$compartment == "Heartwood"]   # wood already dry
cat(sprintf("\nConsistent DRY-basis heartwood:soil ratio: organic %.0fx, mineral %.0fx\n",
            hw_dry / soil_dry$copies_g_dry_est[1], hw_dry / soil_dry$copies_g_dry_est[2]))
cat("(vs 34x from the current mixed wood-dry/soil-fresh comparison; still wood>>soil.)\n")

cat("\nDNA EXTRACTION EFFICIENCY: the workflow applies NO recovery/efficiency\n")
cat("correction, so all copies/g are conservative LOWER BOUNDS (extraction <100%).\n")
cat("\n")
cat("FIG 7: the height-resolved felled-oak mcrA profile is in data/tables/ddpcr.csv\n")
cat("(plate black_oak_mcrA) with masses in samples.csv; fig07 plots it as copies/g.\n")
sink()
cat(readLines(out_path("copies_per_g_summary.txt")), sep="\n")
