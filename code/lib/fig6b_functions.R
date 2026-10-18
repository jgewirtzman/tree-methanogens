# fig6b_functions.R -- the FAPROTAX functions of Figure 6b: the one list, sourced by
# 07_molecular/01_faprotax_dump_HW_SW.R (heartwood vs sapwood tests) and
# 08_figures/fig06_hydrogenotrophy.R (the panel).
# Inclusion rule: energy and redox functions with mean relative abundance >= 0.1% in a
# compartment. Left out: host-association, pathogen and phototrophy functions, and the
# redundant parent annotations (generic chemoheterotrophy; generic, CO2 and methyl
# methanogenesis). Hydrogenotrophic methanogenesis stands for the methanogenesis routes.
FIG6B_FUNCTIONS <- data.frame(
  func = c("hydrogenotrophic_methanogenesis", "dark_hydrogen_oxidation", "fermentation",
           "anaerobic_chemoheterotrophy", "methylotrophy", "nitrate_reduction",
           "aerobic_chemoheterotrophy", "methanol_oxidation", "methanotrophy", "cellulolysis"),
  label = c("H₂-trophic methanogenesis", "Dark hydrogen oxidation", "Fermentation",
            "Anaerobic chemoheterotrophy", "Methylotrophy", "Nitrate reduction",
            "Aerobic chemoheterotrophy", "Methanol oxidation", "Methanotrophy", "Cellulolysis"),
  stringsAsFactors = FALSE)
