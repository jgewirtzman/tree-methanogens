# ==============================================================================
# REVISION — Report pmoA and mmoX SEPARATELY (Referee 2 #1c)
# ==============================================================================
# R2: "pmoA and mmoX abundance should be reported separately as many methanotrophs
# encode both." We currently combine them as methanotroph_total. Separated, a clear
# biological pattern emerges (which is R2's point): pmoA dominates in soil, while in
# wood -- especially heartwood -- pmoA and mmoX are comparable, consistent with a
# larger relative role for mmoX-bearing / pMMO-lacking methanotrophs (e.g.
# Methylocella, Methyloferula) in wood.
#
# NOTE: absolute copies/g are subject to the pending x10 conversion question
# (see ddpcr_controls_LoD.md); the pmoA:mmoX RATIO and all relative patterns here
# are UNAFFECTED by that.
#
# NEW file. Run: Rscript revision/analysis/R_pmoa_mmox_separate.R
# Outputs: revision/outputs/pmoa_mmox_by_compartment.csv, _by_species.csv, _summary.txt
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

d <- read_csv("data/compiled/ddpcr_gene_abundances.csv", show_col_types = FALSE) %>%
  filter(analysis_type == "loose", target_gene %in% c("pmoa","mmox")) %>%
  mutate(compartment = case_when(
    material=="Wood" & core_type=="Inner" ~ "Heartwood",
    material=="Wood" & core_type=="Outer" ~ "Sapwood",
    material=="Soil" ~ paste0("Soil_", core_type), TRUE ~ NA_character_)) %>%
  filter(!is.na(compartment))

# ---- by compartment ----------------------------------------------------------
comp_tab <- d %>% group_by(compartment, target_gene) %>%
  summarise(n = n(), pct_detect = round(100*mean(positives > 0, na.rm = TRUE)),
            median_conc = median(concentration_copies_per_uL, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = target_gene, values_from = c(n, pct_detect, median_conc)) %>%
  mutate(pmoA_to_mmoX_ratio = round(median_conc_pmoa / median_conc_mmox, 2))
write.csv(comp_tab, file.path(out, "pmoa_mmox_by_compartment.csv"), row.names = FALSE)

# ---- by species (wood only) --------------------------------------------------
sp_tab <- d %>% filter(material == "Wood", !is.na(species)) %>%
  mutate(species_name = sp_map[species]) %>% filter(!is.na(species_name)) %>%
  group_by(species_name, target_gene) %>%
  summarise(median_conc = median(concentration_copies_per_uL, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = target_gene, values_from = median_conc, names_prefix = "median_") %>%
  mutate(pmoA_to_mmoX = round(median_pmoa / median_mmox, 2)) %>% arrange(desc(pmoA_to_mmoX))
write.csv(sp_tab, file.path(out, "pmoa_mmox_by_species.csv"), row.names = FALSE)

sink(file.path(out, "pmoa_mmox_summary.txt"))
cat("=================================================================\n")
cat("pmoA and mmoX reported SEPARATELY (Referee 2 #1c)\n")
cat("=================================================================\n\n")
cat("By compartment (median concentration, copies/uL; detection %):\n")
print(as.data.frame(comp_tab %>% transmute(compartment,
  pmoA_n = n_pmoa, pmoA_detect = pct_detect_pmoa, pmoA_median = round(median_conc_pmoa,3),
  mmoX_n = n_mmox, mmoX_detect = pct_detect_mmox, mmoX_median = round(median_conc_mmox,3),
  pmoA_to_mmoX_ratio)), row.names = FALSE)
cat("\nKEY PATTERN (R2's point that combining obscures):\n")
cat(" * SOIL: pmoA >> mmoX (ratio ~", comp_tab$pmoA_to_mmoX_ratio[grepl("Soil_Min",comp_tab$compartment)],
    "in mineral) -- pmoA-bearing (pMMO) methanotrophs dominate soil.\n")
cat(" * HEARTWOOD: pmoA ~ mmoX (ratio ~", comp_tab$pmoA_to_mmoX_ratio[comp_tab$compartment=="Heartwood"],
    ") -- mmoX-bearing / pMMO-lacking lineages relatively more important in wood.\n")
cat(" * This is consistent with the sMMO/copper interpretation and the mmoX",
    "predictive value already noted (Methods S1, L695-700).\n\n")
cat("By species (wood, pmoA:mmoX ratio):\n")
print(as.data.frame(sp_tab %>% transmute(species = species_name,
  pmoA = round(median_pmoa,3), mmoX = round(median_mmox,3), pmoA_to_mmoX)), row.names = FALSE)
cat("\nCAVEAT: absolute copies subject to pending x10 conversion (ddpcr_controls_LoD.md);\n")
cat("the pmoA:mmoX ratio and relative patterns here are unaffected by that.\n")
sink()
cat(readLines(file.path(out,"pmoa_mmox_summary.txt")), sep="\n")
