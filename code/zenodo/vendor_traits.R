#!/usr/bin/env Rscript
# ==============================================================================
# vendor_traits.R -- build data/processed/traits/ymf_species_traits.csv
#
# Figure S17 relates plant traits to methane cycling. The trait table is
# maintained in the companion tree-gas-traits repository; this copies the
# species-level columns the figure uses into this repository so it builds
# without that checkout.
#
# Run from the repository root, with tree-gas-traits checked out alongside:
#   Rscript code/zenodo/vendor_traits.R
#
# Only the 16 columns Figure S17 reads are copied, one row per species. The
# script checks that each column is constant within a species before collapsing.
# ==============================================================================
suppressMessages(library(dplyr)); suppressMessages(library(readr))

SRC <- Sys.getenv("TREE_GAS_TRAITS",
         "../tree-gas-traits/data/processed/ymf_with_traits.csv")
keep <- c("wood_density_gcm3","bark_density_gcm3","bark_wood_ratio","porosity_num","gymnosperm",
          "wood_sapwood_pH","wood_heartwood_pH","try_wood_CN_ratio","try_stem_C","try_antifungal",
          "try_CWD_stem_decomp_rate_k","try_rooting_depth","try_fine_root_diameter",
          "try_fine_root_tissue_density","try_plant_height","try_plant_longevity")

x <- read_csv(SRC, show_col_types = FALSE)
cat("source:", nrow(x), "rows x", ncol(x), "cols\n")

# One row per species: the figure does group_by(spcode) %>% slice(1), so these
# columns are already species-constant. Verify that before collapsing.
chk <- x %>% group_by(spcode) %>%
  summarise(across(any_of(keep), ~ n_distinct(.x[!is.na(.x)])), .groups = "drop")
nonconst <- chk %>% select(-spcode) %>% summarise(across(everything(), ~ max(.x, na.rm = TRUE)))
bad <- names(nonconst)[unlist(nonconst) > 1]
cat("columns varying WITHIN a species:", if (length(bad)) paste(bad, collapse = ", ") else "none", "\n")

out <- x %>% group_by(spcode) %>% slice(1) %>% ungroup() %>%
  select(spcode, any_of(keep)) %>% arrange(spcode)

dir.create("data/processed/traits", showWarnings = FALSE, recursive = TRUE)
write_csv(out, "data/processed/traits/ymf_species_traits.csv")
cat("wrote data/processed/traits/ymf_species_traits.csv:",
    nrow(out), "species x", ncol(out), "cols,",
    file.size("data/processed/traits/ymf_species_traits.csv"), "bytes\n")
cat("species:", paste(out$spcode, collapse = ", "), "\n")
