# ==============================================================================
# ddpcr_constants.R -- converting a QX200 concentration to gene copies per gram
# ------------------------------------------------------------------------------
# The QX200 reports Conc(copies/µL) per µL of the 1x reaction (droplet volume
# 0.85 nL; "Copies/20µLWell" is Conc x 20, software bookkeeping only). Each reaction
# held 2.5 µL of DNA extract in 25 µL (12.5 µL 2x supermix, primers, probe, water;
# Arnold et al. 2024, SI Table 1). Wood extracts were eluted in 75 µL (Arnold et al.
# 2024) and soil extracts in 100 µL (Arnold et al. 2025, Nature, same extracts), so
#
#   copies g-1 = Conc x (25 / 2.5) x elution (75 wood, 100 soil) x (800 / 250) / mass (g)
#
# The last factor scales the eluate back to the whole lysate: samples were lysed in
# 800 µL and 250 µL of cleared lysate was taken for cleanup (method-study Dryad
# records; Methods S3). 250 µL was the target and transfers never exceeded it
# (127-238 µL taken), so the correction is a lower bound. Values are not corrected
# for extraction losses (release from wood, cleanup: together ~80% kept in the method
# study) or for freeze-drying and grinding. Check: the methods paper's spike
# recoveries (22% wood, 30% liquid control) reproduce from its Dryad data only with
# the x10 reaction dilution included; the 30% liquid control is ~250/800.
#
# Sourced by: 03_merge/04_harmonize_all_data.R, 08_figures/fig07_decay-methanogenesis.R,
#   09_tables_stats/06_copies-per-gram.R, 09_tables_stats/07_mass-basis-sensitivity.R
# ==============================================================================
DDPCR_REACTION_UL <- 25     # total reaction volume (µL)
DDPCR_TEMPLATE_UL <- 2.5    # DNA extract per reaction (µL)
DDPCR_ELUTION_UL  <- c(Wood = 75, Soil = 100)   # elution volume (µL) by material
DDPCR_LYSATE_UL    <- 800                         # lysis buffer per sample (µL)
DDPCR_PROCESSED_UL <- c(Wood = 250, Soil = 250)   # cleared lysate taken to cleanup (µL); confirm with Wyatt
# material: "Wood"/"Soil", or a compartment ("Inner", "Outer", "Mineral", "Organic")
ddpcr_elution_ul <- function(material) {
  soil <- material %in% c("Soil", "Mineral", "Organic")
  stopifnot(all(soil | material %in% c("Wood", "Inner", "Outer", "Heartwood", "Sapwood")))
  unname(ifelse(soil, DDPCR_ELUTION_UL[["Soil"]], DDPCR_ELUTION_UL[["Wood"]]))
}
ddpcr_lysate_factor <- function(material) {
  soil <- material %in% c("Soil", "Mineral", "Organic")
  unname(DDPCR_LYSATE_UL / ifelse(soil, DDPCR_PROCESSED_UL[["Soil"]], DDPCR_PROCESSED_UL[["Wood"]]))
}
# extract concentration (copies per µL of eluate, e.g. facility qPCR) -> copies per g
extract_copies_per_g <- function(conc_per_ul_extract, mass_mg, material)
  conc_per_ul_extract * ddpcr_elution_ul(material) * ddpcr_lysate_factor(material) / (mass_mg / 1000)
ddpcr_copies_per_g <- function(conc_per_ul_reaction, mass_mg, material)
  extract_copies_per_g(conc_per_ul_reaction * (DDPCR_REACTION_UL / DDPCR_TEMPLATE_UL), mass_mg, material)
# The per-sample tables (ddPCR_meta_all_data.csv; tree_data_methanogen_group.csv, the
# *_loose / *_strict gene columns) hold "Copies/20µLWell" = Conc x 20, NOT Conc. A 20 µL
# well holds 20 x 2.5/25 = 2 µL of extract, so copies per µL of extract = well / 2.
# The tree-microbiome (Nature 2025) code used exactly this: well x (75/2) / mass.
DDPCR_WELL_UL <- 20
ddpcr_copies_per_g_from_well <- function(copies_per_well, mass_mg, material)
  ddpcr_copies_per_g(copies_per_well / DDPCR_WELL_UL, mass_mg, material)
