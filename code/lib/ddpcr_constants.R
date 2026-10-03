# ==============================================================================
# ddpcr_constants.R -- converting a QX200 concentration to gene copies per gram
# ------------------------------------------------------------------------------
# The QX200 reports Conc(copies/µL) per µL of the 1x reaction (droplet volume
# 0.85 nL; "Copies/20µLWell" is Conc x 20, software bookkeeping only). Each reaction
# held 2.5 µL of DNA extract in 25 µL (12.5 µL 2x supermix, primers, probe, water;
# Arnold et al. 2024, SI Table 1). Wood extracts were eluted in 75 µL (Arnold et al.
# 2024) and soil extracts in 100 µL (Arnold et al. 2025, Nature, same extracts), so
#
#   copies g-1 = Conc x (25 / 2.5) x elution (75 wood, 100 soil) / mass (g)
#
# Values are recovered copies, not corrected for extraction recovery (~20% for wood,
# measured against the whole spike; it includes processing ~250 of the 800 µL of
# lysate). Check: the methods paper's spike recoveries (22% wood, 30% liquid control)
# reproduce from its Dryad data only with the x10 reaction dilution included.
#
# Sourced by: 03_merge/04_harmonize_all_data.R, 08_figures/fig07_decay-methanogenesis.R,
#   09_tables_stats/06_copies-per-gram.R, 09_tables_stats/07_mass-basis-sensitivity.R
# ==============================================================================
DDPCR_REACTION_UL <- 25     # total reaction volume (µL)
DDPCR_TEMPLATE_UL <- 2.5    # DNA extract per reaction (µL)
DDPCR_ELUTION_UL  <- c(Wood = 75, Soil = 100)   # elution volume (µL) by material
# material: "Wood"/"Soil", or a compartment ("Inner", "Outer", "Mineral", "Organic")
ddpcr_elution_ul <- function(material) {
  soil <- material %in% c("Soil", "Mineral", "Organic")
  stopifnot(all(soil | material %in% c("Wood", "Inner", "Outer", "Heartwood", "Sapwood")))
  unname(ifelse(soil, DDPCR_ELUTION_UL[["Soil"]], DDPCR_ELUTION_UL[["Wood"]]))
}
# extract concentration (copies per µL of eluate, e.g. facility qPCR) -> copies per g
extract_copies_per_g <- function(conc_per_ul_extract, mass_mg, material)
  conc_per_ul_extract * ddpcr_elution_ul(material) / (mass_mg / 1000)
ddpcr_copies_per_g <- function(conc_per_ul_reaction, mass_mg, material)
  extract_copies_per_g(conc_per_ul_reaction * (DDPCR_REACTION_UL / DDPCR_TEMPLATE_UL), mass_mg, material)
