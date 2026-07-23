# ==============================================================================
# REVISION — Scaling scenario: wood vs soil standing methanogen GENE POOL per m2
# ==============================================================================
# Jon's question: heartwood is ~25x more concentrated (mcrA copies/g) than soil,
# but a forest has vastly more soil mass than heartwood mass. So per unit GROUND
# AREA, how do the standing methanogen gene pools compare? A sensitivity scenario.
#
# pool (copies/m2) = copies/g_dry * dry_mass (g/m2), for heartwood and for soil.
# STRICT CAVEAT: gene copies != activity != flux. Soil methanogens' CH4 is largely
# oxidized before emission (soils here are net sinks); wood methanogens sit near
# the emitting surface. This scenario is about STANDING genetic potential only.
#
# NEW file; edits nothing. Run: Rscript revision/analysis/R_scaling_wood_soil_pool.R
# Outputs: revision/outputs/scaling_wood_soil_pool.txt, scaling_wood_soil_pool.csv
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "revision/outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- copies/g DRY from our own outputs ---------------------------------------
comp <- read_csv(file.path(out_dir, "copies_per_g_compartment.csv"), show_col_types = FALSE)
soil_dry <- read_csv(file.path(out_dir, "copies_per_g_soil_dry_estimate.csv"), show_col_types = FALSE)
cpg_heartwood <- comp$median_real[comp$compartment == "Heartwood"]      # wood already dry
cpg_soil_org  <- soil_dry$copies_g_dry_est[soil_dry$compartment == "Soil_Organic"]
cpg_soil_min  <- soil_dry$copies_g_dry_est[soil_dry$compartment == "Soil_Mineral"]

# ---- WOOD (heartwood) dry mass per m2 -----------------------------------------
# NOTE: data/compiled/forest_inventory.csv is a SAPLING subplot (median DBH 2.7 cm,
# only 7 stems >=10 cm) -- it CANNOT give stand biomass. So we use an assumed
# above-ground biomass (AGB) range from the temperate-hardwood literature instead.
# (If the full ForestGEO 8-ha census is located, swap it in for a site-specific AGB.)
AGB_MGHA <- 200        # Mg/ha, mid temperate hardwood (range 100-300 explored below)
agb_kg_m2 <- AGB_MGHA * 1e6 / 1e4 / 1e3                     # Mg/ha -> kg/m2  (= AGB_MGHA/10)
STEM_FRAC <- 0.65      # stem wood as fraction of AGB
HW_FRAC   <- 0.30      # heartwood as fraction of stem wood (mid estimate)
heartwood_g_m2 <- agb_kg_m2 * STEM_FRAC * HW_FRAC * 1000

# ---- SOIL dry mass per m2 from measured depths + literature bulk density ------
tp <- read_csv("data/compiled/tree_properties.csv", show_col_types = FALSE)
org_depth_m <- median(tp$OrganicDepth_mean, na.rm = TRUE) / 100
min_depth_m <- median(tp$MineralDepth_mean, na.rm = TRUE) / 100
BD_ORG <- 0.20; BD_MIN <- 1.20                              # g/cm3 (kg/L); literature mid
soil_org_g_m2 <- org_depth_m * BD_ORG * 1e6                 # depth(m)*BD(g/cm3)*1e6 = g/m2
soil_min_g_m2 <- min_depth_m * BD_MIN * 1e6

# ---- pools -------------------------------------------------------------------
pool_heartwood <- cpg_heartwood * heartwood_g_m2
pool_soil <- cpg_soil_org * soil_org_g_m2 + cpg_soil_min * soil_min_g_m2

res <- tibble(
  compartment = c("Heartwood", "Soil (organic+mineral)"),
  copies_per_g_dry = c(cpg_heartwood, NA),
  dry_mass_g_m2 = c(heartwood_g_m2, soil_org_g_m2 + soil_min_g_m2),
  pool_copies_m2 = c(pool_heartwood, pool_soil))
write.csv(res, file.path(out_dir, "scaling_wood_soil_pool.csv"), row.names = FALSE)

# ---- sensitivity over AGB, HW fraction, and mineral bulk density -------------
sens <- expand_grid(agb_mgha = c(100,200,300), hw_frac = c(0.20,0.30,0.45), bd_min = c(1.0,1.2,1.4)) %>%
  mutate(hw_g = (agb_mgha/10)*STEM_FRAC*hw_frac*1000,
         ratio_wood_to_soil = (cpg_heartwood*hw_g) /
           (cpg_soil_org*soil_org_g_m2 + cpg_soil_min*min_depth_m*bd_min*1e6))

sink(file.path(out_dir, "scaling_wood_soil_pool.txt"))
cat("=================================================================\n")
cat("SCALING SCENARIO: standing methanogen gene pool, wood vs soil, per m2\n")
cat("=================================================================\n\n")
cat("Inputs (all clearly assumptions where noted):\n")
cat(sprintf("  copies/g DRY: heartwood %.0f ; soil organic %.1f ; soil mineral %.1f\n",
            cpg_heartwood, cpg_soil_org, cpg_soil_min))
cat(sprintf("  AGB = %d Mg/ha = %.1f kg/m2 [ASSUMED, temperate-hardwood literature; range 100-300 below]\n",
            AGB_MGHA, agb_kg_m2))
cat(sprintf("  heartwood mass = %.0f g/m2 (AGB x stem %.2f x heartwood %.2f) [ASSUMED fractions]\n",
            heartwood_g_m2, STEM_FRAC, HW_FRAC))
cat(sprintf("  soil depths: organic %.0f cm, mineral %.0f cm (measured); BD org %.1f, min %.1f g/cm3 [ASSUMED]\n",
            org_depth_m*100, min_depth_m*100, BD_ORG, BD_MIN))
cat(sprintf("  soil dry mass = %.0f g/m2\n\n", soil_org_g_m2 + soil_min_g_m2))
cat("STANDING GENE POOLS (copies mcrA per m2 ground):\n")
print(as.data.frame(res %>% mutate(pool_copies_m2 = signif(pool_copies_m2,3))), row.names = FALSE)
cat(sprintf("\n  Wood : Soil pool ratio (mid scenario, AGB=%d Mg/ha) = %.2f\n", AGB_MGHA, pool_heartwood/pool_soil))
cat(sprintf("  Heartwood is ~%.0fx more concentrated per gram, but soil has ~%.0fx more\n",
            cpg_heartwood/mean(c(cpg_soil_org,cpg_soil_min)), (soil_org_g_m2+soil_min_g_m2)/heartwood_g_m2))
cat("  dry mass. Net: standing pools are the SAME ORDER of magnitude -- the per-gram\n")
cat("  wood enrichment substantially compensates for wood's much smaller mass.\n\n")
cat("Sensitivity — wood:soil pool ratio across AGB x heartwood-fraction x soil BD:\n")
print(as.data.frame(sens %>% transmute(agb_mgha, hw_frac, bd_min, ratio = round(ratio_wood_to_soil,2))), row.names = FALSE)
cat(sprintf("\n  Range: wood:soil pool ratio = %.2f to %.2f (mostly 0.2-0.8) -> wood standing\n",
            min(sens$ratio_wood_to_soil), max(sens$ratio_wood_to_soil)))
cat("  methanogen pool is ~20-80%% of the soil pool: comparable, not negligible.\n")
cat("\nCAVEATS (report prominently): (1) gene copies are NOT activity or flux;\n")
cat("(2) soil methanogen CH4 is largely oxidized before emission (soils net sinks\n")
cat("here) whereas wood methanogens sit near the emitting bark surface, so equal\n")
cat("POOLS do NOT imply equal emissions; (3) copies/g are extraction lower bounds;\n")
cat("(4) soil dry masses use black-oak GWC + literature bulk densities; (5) shallow\n")
cat("soil only (~20 cm) -- deeper soil adds mass but little methanogen relevance here.\n")
cat("This is a bounding CONTEXT scenario, not a budget.\n")
sink()
cat(readLines(file.path(out_dir,"scaling_wood_soil_pool.txt")), sep="\n")
