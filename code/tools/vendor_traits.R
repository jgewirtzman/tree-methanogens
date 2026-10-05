#!/usr/bin/env Rscript
# ==============================================================================
# vendor_traits.R -- build data/raw/external/tree-gas-traits/ymf_species_traits.csv
#
# Fig S25 and code/09_tables_stats/25_plant-traits.R relate plant traits to methane
# cycling. The trait table is maintained in the companion tree-gas-traits repository;
# this copies the species-level columns they use into this repository so they build
# without that checkout.
#
# Run from the repository root, with tree-gas-traits checked out alongside:
#   Rscript code/tools/vendor_traits.R
#
# Only the columns the analysis reads are copied, one row per species: the 16 traits
# that pass the inclusion rules in 25_plant-traits.R (Methods S5) and the gymnosperm
# indicator. The script checks that each column is constant within a species before
# collapsing. (Before 2026-10-04 the copied set was the 16 columns of the earlier
# heatmap; the shared columns are unchanged.)
# ==============================================================================
suppressMessages(library(dplyr)); suppressMessages(library(readr))

SRC <- Sys.getenv("TREE_GAS_TRAITS",
         "../tree-gas-traits/data/clean/species_traits.csv")
keep <- c("gymnosperm", "wood_density_gcm3","bark_density_gcm3","try_bark_volume_rel_wood",
          "wood_sapwood_pH","wood_heartwood_pH","try_wood_CN_ratio","try_stem_C","try_bark_C",
          "try_CWD_stem_decomp_rate_k","try_waterlogging_tolerance","try_rooting_depth",
          "try_fine_root_diameter","try_fine_root_tissue_density","try_ectomycorrhizal",
          "try_plant_longevity")

x <- read_csv(SRC, show_col_types = FALSE)
# the 16 species sampled at Yale-Myers (the source table also covers other sites)
ymf <- c("ACRU","ACSA","BEAL","BELE","BEPA","CAOV","FAGR","FRAM","KALA","PIST","PRSE","QUAL","QURU","QUVE","SAAL","TSCA")
x <- x %>% filter(spcode %in% ymf)
cat("source:", nrow(x), "rows x", ncol(x), "cols\n")

# One row per species; verify the columns are species-constant before collapsing.
chk <- x %>% group_by(spcode) %>%
  summarise(across(any_of(keep), ~ n_distinct(.x[!is.na(.x)])), .groups = "drop")
nonconst <- chk %>% select(-spcode) %>% summarise(across(everything(), ~ max(.x, na.rm = TRUE)))
bad <- names(nonconst)[unlist(nonconst) > 1]
cat("columns varying WITHIN a species:", if (length(bad)) paste(bad, collapse = ", ") else "none", "\n")

out <- x %>% group_by(spcode) %>% slice(1) %>% ungroup() %>%
  select(spcode, any_of(keep)) %>% arrange(spcode)

dir.create("data/raw/external/tree-gas-traits", showWarnings = FALSE, recursive = TRUE)
write_csv(out, "data/raw/external/tree-gas-traits/ymf_species_traits.csv")
cat("wrote data/raw/external/tree-gas-traits/ymf_species_traits.csv:",
    nrow(out), "species x", ncol(out), "cols,",
    file.size("data/raw/external/tree-gas-traits/ymf_species_traits.csv"), "bytes\n")
cat("species:", paste(out$spcode, collapse = ", "), "\n")
