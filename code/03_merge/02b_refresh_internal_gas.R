# ==============================================================================
# 02b_refresh_internal_gas.R -- swap the recalibrated internal-gas columns into
# the merged tree dataset
# ------------------------------------------------------------------------------
# WHY THIS EXISTS. 02_harmonize_all_data.R cannot currently be rerun:
# it reads data/processed/flux/methanogen_tree_flux_complete_dataset.csv as the
# 2021 multi-height flux table, but 02_goflux_trees_2021.R and
# 02_goflux_trees_2023.R both write to that filename, and the 2023 run
# (Sep 2025) overwrote the 2021 output. merged_tree_dataset_final.csv on disk
# predates the overwrite and is otherwise correct.
#
# This script changes ONLY the four gas columns (CO2/CH4/N2O/O2_concentration),
# rebuilt exactly as 02_harmonize_all_data.R builds them (same tree-ID
# standardiser, same join key), after 03_process_internal_gas.R was
# recalibrated (2026-09-30). Every other column is checked to be unchanged.
# Retire it once the filename collision is fixed and 02_harmonize reruns.
#
# Run from code/03_merge/:  Rscript 02b_refresh_internal_gas.R
# Reads:  data/processed/internal_gas/sample_data_only.csv
#         data/processed/tree_data/tree_id_comprehensive_mapping.csv
#         data/processed/integrated/merged_tree_dataset_final.csv
# Writes: data/processed/integrated/merged_tree_dataset_final.csv (gas columns)
# ==============================================================================
suppressMessages({ library(readr); library(dplyr) })

# --- tree-ID standardiser: evaluated from 02_harmonize_all_data.R itself, so the
#     two scripts cannot drift apart
src  <- parse("02_harmonize_all_data.R")
want <- c("mapping", "create_comprehensive_lookup", "tree_lookup", "standardize_tree_id")
for (e in src) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% want) eval(e)
}
stopifnot(exists("standardize_tree_id"))

GAS <- c("CO2_concentration", "CH4_concentration", "N2O_concentration", "O2_concentration")
gas <- read_csv("../../data/processed/internal_gas/sample_data_only.csv", show_col_types = FALSE) %>%
  transmute(tree_id = standardize_tree_id(Tree.ID), across(all_of(GAS))) %>%
  filter(!is.na(tree_id), tree_id != "untagged")
stopifnot(!anyDuplicated(gas$tree_id))

f   <- "../../data/processed/integrated/merged_tree_dataset_final.csv"
old <- read_csv(f, show_col_types = FALSE, guess_max = 1e5)
new <- old %>% select(-all_of(GAS)) %>% left_join(gas, by = "tree_id") %>% select(all_of(names(old)))

# trees with gas but absent from the merged table would have been added as new
# rows by the full_join in 02_harmonize; they must not exist
missing <- setdiff(gas$tree_id, old$tree_id)
if (length(missing)) stop("gas trees absent from the merged table: ", paste(missing, collapse = ", "))
other <- setdiff(names(old), GAS)
stopifnot(identical(as.data.frame(old[other]), as.data.frame(new[other])))
stopifnot(sum(!is.na(new$CH4_concentration)) == sum(!is.na(old$CH4_concentration)))

write_csv(new, f)
cat(sprintf("Refreshed %s for %d trees (%d rows). CH4 exact zeros: before %d, after %d\n",
            paste(GAS, collapse = ", "), sum(!is.na(new$CH4_concentration)), nrow(new),
            sum(old$CH4_concentration == 0, na.rm = TRUE), sum(new$CH4_concentration == 0, na.rm = TRUE)))
