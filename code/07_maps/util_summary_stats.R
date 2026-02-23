# ==============================================================================
# Spatial Summary Statistics
# ==============================================================================
# Purpose: Summary statistics for spatial data.
# ==============================================================================

cat("\n=== DIRECT FLUX COMPARISON ===\n\n")

# Constants
WAI_temperate <- 3.07  # Gauci et al. 2024

# =============================================================================
# STEP 1: Per-unit-area fluxes (the actual measurements)
# =============================================================================

# These are your predictions per m² of stem surface or soil surface
flux_per_stem_surface <- monthly_results_nmol$Phi_tree_nmol_m2_s / stem_area_fraction
flux_per_soil_surface <- monthly_results_nmol$Phi_soil_nmol_m2_s

cat("PER-UNIT-AREA FLUXES (nmol m⁻² s⁻¹):\n")
cat("Month  | Stem Surface | Soil Surface | Ratio\n")
cat("-------|--------------|--------------|-------\n")
for (i in 1:12) {
  ratio <- abs(flux_per_stem_surface[i] / flux_per_soil_surface[i])
  cat(sprintf("%-6s | %12.3f | %12.2f | %.3f\n",
              month.abb[i], flux_per_stem_surface[i], 
              flux_per_soil_surface[i], ratio))
}

mean_stem <- mean(flux_per_stem_surface)
mean_soil <- mean(flux_per_soil_surface)
ratio_means <- mean_stem / abs(mean_soil)

cat("\nAnnual means:\n")
cat("  Stem: ", sprintf("%.3f nmol m⁻² s⁻¹\n", mean_stem))
cat("  Soil: ", sprintf("%.2f nmol m⁻² s⁻¹\n", mean_soil))
cat("  Ratio: ", sprintf("%.3f (%.1f%%)\n", ratio_means, 100 * ratio_means))

# =============================================================================
# STEP 2: Scale to WAI = 3.07 for plot-level fluxes
# =============================================================================

cat("\n=== SCALING TO WAI = 3.07 ===\n\n")

# Plot-level flux = per-stem-area flux × WAI
flux_plot_tree_WAI <- flux_per_stem_surface * WAI_temperate
flux_plot_soil <- flux_per_soil_surface  # unchanged

cat("PLOT-LEVEL FLUXES with WAI = 3.07 (nmol m⁻² s⁻¹):\n")
cat("Month  | Tree | Soil\n")
cat("-------|------|------\n")
for (i in 1:12) {
  cat(sprintf("%-6s | %4.2f | %5.2f\n",
              month.abb[i], flux_plot_tree_WAI[i], flux_plot_soil[i]))
}

# =============================================================================
# STEP 3: Annual budgets
# =============================================================================

cat("\n=== ANNUAL BUDGETS ===\n\n")

# Conversion: nmol m⁻² s⁻¹ to mg m⁻² yr⁻¹
# nmol/m²/s × (1 μmol/1000 nmol) × (86400 s/day) × (30.4 days/month) × (12 months/yr) × (16 g/mol) × (0.001 mg/g)
conv <- (1/1000) * 86400 * 30.4 * 12 * 16 * 0.001

annual_tree_WAI <- sum(flux_plot_tree_WAI) * conv / 12  # divide by 12 because we sum 12 months
annual_soil <- sum(flux_plot_soil) * conv / 12

offset_WAI <- abs(annual_tree_WAI / annual_soil) * 100

cat("With WAI = 3.07:\n")
cat("  Tree emissions: ", round(annual_tree_WAI, 1), " mg CH₄ m⁻² yr⁻¹\n")
cat("  Soil uptake: ", round(annual_soil, 1), " mg CH₄ m⁻² yr⁻¹\n")
cat("  Offset: ", round(offset_WAI, 1), "%\n\n")

# Summer months
summer_idx <- 6:8
summer_tree_WAI <- flux_plot_tree_WAI[summer_idx]
summer_soil <- flux_plot_soil[summer_idx]

cat("Summer (Jun-Aug) fluxes with WAI = 3.07:\n")
cat("  Tree: ", sprintf("%.2f-%.2f nmol m⁻² s⁻¹\n", 
                        min(summer_tree_WAI), max(summer_tree_WAI)))
cat("  Soil: ", sprintf("%.2f-%.2f nmol m⁻² s⁻¹\n",
                        min(summer_soil), max(summer_soil)))

# =============================================================================
# FORMATTED TEXT
# =============================================================================

cat("\n=== MANUSCRIPT TEXT ===\n\n")

cat("On a per-unit-area basis, tree stem emissions averaged ",
    sprintf("%.2f nmol m⁻² s⁻¹", mean_stem),
    " (range: ", sprintf("%.2f-%.2f", min(flux_per_stem_surface), max(flux_per_stem_surface)),
    "), while soil uptake averaged ",
    sprintf("%.2f nmol m⁻² s⁻¹", abs(mean_soil)),
    " (range: ", sprintf("%.2f-%.2f", abs(max(flux_per_soil_surface)), abs(min(flux_per_soil_surface))),
    "). Per unit area, stem emissions represented ",
    sprintf("%.1f%%", 100 * ratio_means),
    " of soil uptake, with peak emissions in July (",
    sprintf("%.2f nmol m⁻² s⁻¹", max(flux_per_stem_surface)), ").\n\n")

cat("Scaling to a temperate forest with WAI = 3.07 (Gauci et al., 2024), ",
    "tree emissions would contribute ", 
    sprintf("%.2f-%.2f nmol m⁻² s⁻¹", min(flux_plot_tree_WAI), max(flux_plot_tree_WAI)),
    " at the plot level (summer: ",
    sprintf("%.2f-%.2f", min(summer_tree_WAI), max(summer_tree_WAI)),
    " nmol m⁻² s⁻¹), totaling ",
    round(annual_tree_WAI, 0), " mg CH₄ m⁻² yr⁻¹ annually ",
    "and offsetting ", round(offset_WAI, 0), "% of soil uptake (",
    round(abs(annual_soil), 0), " mg CH₄ m⁻² yr⁻¹).\n\n")

cat("=== COMPLETE ===\n")
















# =============================================================================
# VERIFY ALL MANUSCRIPT VALUES
# =============================================================================

cat("\n=== VERIFICATION OF ALL MANUSCRIPT VALUES ===\n\n")

# =============================================================================
# 1. MODEL PERFORMANCE
# =============================================================================

cat("1. MODEL PERFORMANCE:\n")
cat("   Tree model R²: ", round(TreeRF$r.squared, 2), "\n")
cat("   Soil model R²: ", round(SoilRF$r.squared, 2), "\n\n")

# =============================================================================
# 2. PER-UNIT-AREA FLUXES
# =============================================================================

cat("2. PER-UNIT-AREA FLUXES:\n\n")

# Calculate per-unit-area fluxes
stem_area_fraction <- total_stem_surface / PLOT_AREA
flux_per_stem_m2 <- monthly_results_nmol$Phi_tree_nmol_m2_s / stem_area_fraction
flux_per_soil_m2 <- monthly_results_nmol$Phi_soil_nmol_m2_s

cat("   Stem flux (nmol m⁻² s⁻¹):\n")
cat("     Mean: ", round(mean(flux_per_stem_m2), 2), "\n")
cat("     Min:  ", round(min(flux_per_stem_m2), 2), "\n")
cat("     Max:  ", round(max(flux_per_stem_m2), 2), "\n")
cat("     Range: ", round(min(flux_per_stem_m2), 2), "-", 
    round(max(flux_per_stem_m2), 2), "\n\n")

cat("   Soil flux (nmol m⁻² s⁻¹):\n")
cat("     Mean: ", round(mean(flux_per_soil_m2), 2), "\n")
cat("     Min:  ", round(min(flux_per_soil_m2), 2), "\n")
cat("     Max:  ", round(max(flux_per_soil_m2), 2), "\n")
cat("     Range: ", round(abs(max(flux_per_soil_m2)), 2), "-", 
    round(abs(min(flux_per_soil_m2)), 2), "\n\n")

# =============================================================================
# 3. RATIO
# =============================================================================

cat("3. STEM/SOIL RATIO:\n")
ratio <- mean(flux_per_stem_m2) / abs(mean(flux_per_soil_m2))
cat("   ", round(ratio, 3), " = ", round(100 * ratio, 1), "%\n\n")

# =============================================================================
# 4. PEAK IN JULY
# =============================================================================

cat("4. PEAK EMISSIONS:\n")
july_stem <- flux_per_stem_m2[7]  # July is month 7
cat("   July stem flux: ", round(july_stem, 2), " nmol m⁻² s⁻¹\n")
cat("   Verify it's the max: ", max(flux_per_stem_m2) == july_stem, "\n\n")

# =============================================================================
# 5. WAI SCALING
# =============================================================================

cat("5. WAI SCALING (WAI = 3.07):\n\n")

WAI <- 3.07
flux_plot_tree_WAI <- flux_per_stem_m2 * WAI

cat("   Plot-level tree flux (nmol m⁻² s⁻¹):\n")
cat("     Mean: ", round(mean(flux_plot_tree_WAI), 2), "\n")
cat("     Min:  ", round(min(flux_plot_tree_WAI), 2), "\n")
cat("     Max:  ", round(max(flux_plot_tree_WAI), 2), "\n\n")

# =============================================================================
# 6. ANNUAL BUDGETS
# =============================================================================

cat("6. ANNUAL BUDGETS:\n\n")

# Conversion factor
# nmol m⁻² s⁻¹ to mg m⁻² yr⁻¹:
# × (1 µmol / 1000 nmol) × (86400 s/day) × (30.4 day/month) × (12 month/yr) × (16 g/mol) × (0.001 mg/g)
conv <- (1/1000) * 86400 * 30.4 * 12 * 16 * 0.001

cat("   Conversion factor: ", conv, "\n\n")

# Annual tree (WAI scaled)
annual_tree_WAI <- mean(flux_plot_tree_WAI) * conv
cat("   Tree emissions (WAI = 3.07):\n")
cat("     Mean flux: ", round(mean(flux_plot_tree_WAI), 3), " nmol m⁻² s⁻¹\n")
cat("     × ", conv, "\n")
cat("     = ", round(annual_tree_WAI, 1), " mg CH₄ m⁻² yr⁻¹\n\n")

# Annual soil
annual_soil <- mean(flux_per_soil_m2) * conv
cat("   Soil uptake:\n")
cat("     Mean flux: ", round(mean(flux_per_soil_m2), 3), " nmol m⁻² s⁻¹\n")
cat("     × ", conv, "\n")
cat("     = ", round(annual_soil, 1), " mg CH₄ m⁻² yr⁻¹\n\n")

# =============================================================================
# 7. OFFSET PERCENTAGE
# =============================================================================

cat("7. OFFSET PERCENTAGE:\n")
offset <- abs(annual_tree_WAI / annual_soil) * 100
cat("   ", round(abs(annual_tree_WAI), 1), " / ", round(abs(annual_soil), 1), 
    " = ", round(offset, 1), "%\n\n")

# =============================================================================
# 8. SUMMARY TABLE
# =============================================================================

cat("=== SUMMARY TABLE FOR MANUSCRIPT ===\n\n")
cat("Model Performance:\n")
cat("  Tree R² = ", round(TreeRF$r.squared, 2), "\n")
cat("  Soil R² = ", round(SoilRF$r.squared, 2), "\n\n")

cat("Per-unit-area fluxes:\n")
cat("  Stem: ", round(mean(flux_per_stem_m2), 2), " nmol m⁻² s⁻¹ (range: ",
    round(min(flux_per_stem_m2), 2), "-", round(max(flux_per_stem_m2), 2), ")\n")
cat("  Soil: ", round(abs(mean(flux_per_soil_m2)), 2), " nmol m⁻² s⁻¹ (range: ",
    round(abs(max(flux_per_soil_m2)), 2), "-", round(abs(min(flux_per_soil_m2)), 2), ")\n")
cat("  Ratio: ", round(100 * ratio, 1), "%\n")
cat("  July peak: ", round(july_stem, 2), " nmol m⁻² s⁻¹\n\n")

cat("WAI scaling (WAI = 3.07):\n")
cat("  Plot-level mean: ", round(mean(flux_plot_tree_WAI), 2), " nmol m⁻² s⁻¹\n")
cat("  Annual tree: ", round(annual_tree_WAI, 0), " mg CH₄ m⁻² yr⁻¹\n")
cat("  Annual soil: ", round(abs(annual_soil), 0), " mg CH₄ m⁻² yr⁻¹\n")
cat("  Offset: ", round(offset, 0), "%\n\n")

cat("=== COMPLETE ===\n")















# =============================================================================
# SENSITIVITY ANALYSIS: VERTICAL PROFILE ASSUMPTIONS
# =============================================================================

cat("\n=== SENSITIVITY TO VERTICAL PROFILE ASSUMPTIONS ===\n\n")

WAI <- 3.07
mean_stem_flux <- mean(flux_per_stem_m2)  # 0.073 nmol m⁻² s⁻¹
conv <- (1/1000) * 86400 * 30.4 * 12 * 16 * 0.001
annual_soil <- mean(flux_per_soil_m2) * conv

cat("Per-unit-area stem flux (measured up to 2m height): ", round(mean_stem_flux, 3), " nmol m⁻² s⁻¹\n\n")

# SCENARIO 1: Zero emissions above 2m (lower bound)
cat("SCENARIO 1 (Lower bound): Zero emissions above 2m\n")
cat("Only measured section contributes. Assuming measured section (0-2m)\n")
cat("represents ~15-25% of total woody surface area:\n\n")

for (frac in c(0.15, 0.20, 0.25)) {
  WAI_eff <- WAI * frac
  flux_plot <- mean_stem_flux * WAI_eff
  annual <- flux_plot * conv
  offset <- 100 * abs(annual / annual_soil)
  cat(sprintf("  %d%% of WAI: %.0f mg m⁻² yr⁻¹ (%.1f%% offset)\n", 
              round(100*frac), annual, offset))
}

# SCENARIO 2: Constant emissions across all woody area (upper bound)
cat("\nSCENARIO 2 (Upper bound): Constant emissions to full canopy height\n")
flux_upper <- mean_stem_flux * WAI
annual_upper <- flux_upper * conv
offset_upper <- 100 * abs(annual_upper / annual_soil)

cat("  Full WAI: ", WAI, "\n")
cat("  Plot-level flux: ", round(flux_upper, 2), " nmol m⁻² s⁻¹\n")
cat("  Annual: ", round(annual_upper, 0), " mg CH₄ m⁻² yr⁻¹\n")
cat("  Offset: ", round(offset_upper, 1), "%\n\n")

# Use 20% for lower bound in text
frac_lower <- 0.20
WAI_lower <- WAI * frac_lower
flux_lower <- mean_stem_flux * WAI_lower
annual_lower <- flux_lower * conv
offset_lower <- 100 * abs(annual_lower / annual_soil)

cat("=== FORMATTED TEXT ===\n\n")

cat("As a scaling exercise, we applied our per-unit-area stem flux measurements\n")
cat("(0.07 nmol m⁻² s⁻¹, measured up to 2 m height) to the Woody Area Index of 3.07\n")
cat("reported for temperate forests by Gauci et al. (2024). Because we lack empirical\n")
cat("woody area data for our site and vertical flux profiles above 2 m, we present\n")
cat("bounding scenarios. Assuming zero emissions above 2 m yields ", round(annual_lower, 0), 
    " mg CH₄ m⁻² yr⁻¹\n")
cat("(", round(offset_lower, 1), "% offset), while assuming constant per-area emissions across all woody\n")
cat("surfaces yields ", round(annual_upper, 0), " mg CH₄ m⁻² yr⁻¹ (", round(offset_upper, 1), 
    "% offset) relative to soil uptake\n")
cat("(", round(abs(annual_soil), 0), " mg CH₄ m⁻² yr⁻¹). These estimates are sensitive to both assumptions about\n")
cat("vertical flux patterns and site-specific woody area. While we observed no evidence\n")
cat("of declining fluxes within our measurement range and detected a flux hotspot at\n")
cat("~4 m height in one tree measured to 10 m, insufficient data exist to constrain\n")
cat("vertical profiles or woody area for this site. This scaling exercise illustrates\n")
cat("potential tree contributions under various assumptions but should not be interpreted\n")
cat("as a plot-level methane budget.\n\n")

cat("=== COMPLETE ===\n")











# =============================================================================
# ACTUAL SITE MEASUREMENTS - LATERAL CYLINDER SURFACE TO 2M
# =============================================================================

cat("\n=== ACTUAL SITE FLUX ESTIMATES ===\n\n")

# Site characteristics
cat("SITE CHARACTERISTICS:\n")
cat("Plot area: ", round(PLOT_AREA, 0), " m² (", round(PLOT_AREA/10000, 2), " ha)\n")
cat("Number of trees: ", nrow(INVENTORY), "\n")
cat("Lateral stem surface area (0-2m height): ", round(total_stem_surface, 0), " m²\n")
cat("  (calculated as cylinder: circumference × height)\n")
cat("Woody surface fraction: ", round(stem_area_fraction, 4), "\n\n")

# Per-unit-area fluxes
flux_per_stem_m2 <- monthly_results_nmol$Phi_tree_nmol_m2_s / stem_area_fraction
flux_per_soil_m2 <- monthly_results_nmol$Phi_soil_nmol_m2_s

mean_stem <- mean(flux_per_stem_m2)
mean_soil <- mean(flux_per_soil_m2)

cat("PER-UNIT-AREA FLUXES:\n")
cat("Stem (0-2m lateral surface): ", round(mean_stem, 3), " nmol m⁻² s⁻¹\n")
cat("Soil: ", round(abs(mean_soil), 2), " nmol m⁻² s⁻¹\n")
cat("Ratio: ", round(100 * mean_stem / abs(mean_soil), 1), "%\n\n")

# Plot-level fluxes (actual site)
flux_plot_tree <- monthly_results_nmol$Phi_tree_nmol_m2_s
flux_plot_soil <- monthly_results_nmol$Phi_soil_nmol_m2_s

cat("PLOT-LEVEL FLUXES (actual site):\n")
cat("Tree (0-2m section): ", round(mean(flux_plot_tree), 4), " nmol m⁻² s⁻¹\n")
cat("Soil: ", round(mean(flux_plot_soil), 2), " nmol m⁻² s⁻¹\n\n")

# Annual budgets
conv <- (1/1000) * 86400 * 30.4 * 12 * 16 * 0.001

annual_tree_site <- mean(flux_plot_tree) * conv
annual_soil_site <- mean(flux_plot_soil) * conv
offset_site <- 100 * abs(annual_tree_site / annual_soil_site)

cat("ANNUAL BUDGETS (actual site, 0-2m stems only):\n")
cat("Tree emissions: ", round(annual_tree_site, 1), " mg CH₄ m⁻² yr⁻¹\n")
cat("Soil uptake: ", round(abs(annual_soil_site), 0), " mg CH₄ m⁻² yr⁻¹\n")
cat("Offset: ", round(offset_site, 2), "%\n\n")

# Scaling scenario for comparison
WAI <- 3.07
flux_plot_WAI <- mean_stem * WAI
annual_WAI <- flux_plot_WAI * conv
offset_WAI <- 100 * abs(annual_WAI / annual_soil_site)

cat("COMPARISON: Scaled to WAI = 3.07 (constant emissions assumption):\n")
cat("Tree emissions: ", round(annual_WAI, 0), " mg CH₄ m⁻² yr⁻¹\n")
cat("Offset: ", round(offset_WAI, 1), "%\n\n")

# =============================================================================
# FORMATTED TEXT
# =============================================================================

cat("=== MANUSCRIPT TEXT ===\n\n")

cat("On a per-unit-area basis, tree stem emissions averaged ", round(mean_stem, 2), 
    " nmol m⁻² s⁻¹\n")
cat("across the year (monthly means ranged from ", round(min(flux_per_stem_m2), 2), "-",
    round(max(flux_per_stem_m2), 2), " nmol m⁻² s⁻¹) while\n")
cat("soil uptake averaged ", round(abs(mean_soil), 2), " nmol m⁻² s⁻¹ (monthly means: ",
    round(abs(max(flux_per_soil_m2)), 2), "-", round(abs(min(flux_per_soil_m2)), 2), 
    " nmol m⁻² s⁻¹).\n")
cat("Per unit area, stem emissions represented ", round(100 * mean_stem / abs(mean_soil), 1),
    "% of soil uptake, with peak\n")
cat("emissions in July (", round(max(flux_per_stem_m2), 2), " nmol m⁻² s⁻¹).\n\n")

cat("At our study site, based on lateral stem surface area of cylinders to 2 m height\n")
cat("(3,555 m² across 10.2 ha, representing ", round(100*stem_area_fraction, 2), 
    "% of plot area), tree emissions\n")
cat("contributed ", round(mean(flux_plot_tree), 3), " nmol m⁻² s⁻¹ on a plot-area basis, totaling\n")
cat(round(annual_tree_site, 1), " mg CH₄ m⁻² yr⁻¹ and offsetting ", 
    round(offset_site, 2), "% of soil uptake\n")
cat("(", round(abs(annual_soil_site), 0), " mg CH₄ m⁻² yr⁻¹). As a scaling exercise,\n")
cat("applying our per-unit-area stem flux to the Woody Area Index of 3.07 reported\n")
cat("for temperate forests by Gauci et al. (2024)—assuming constant emissions across\n")
cat("all woody surfaces—would yield ", round(annual_WAI, 0), " mg CH₄ m⁻² yr⁻¹, offsetting\n")
cat(round(offset_WAI, 0), "% of soil uptake. These estimates are based on measurements\n")
cat("restricted to 0-2 m stem height and lack empirical data on vertical flux profiles\n")
cat("or site-specific woody area above 2 m. While we observed no evidence of declining\n")
cat("fluxes within our measurement range and detected a flux hotspot at ~4 m height in\n")
cat("one tree measured to 10 m, insufficient data exist to constrain these uncertainties.\n\n")

cat("=== COMPLETE ===\n")