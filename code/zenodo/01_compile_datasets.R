# ==============================================================================
# 01_compile_datasets.R -- the archived datasets in data/compiled/
# ------------------------------------------------------------------------------
# Runs at the end of stage B (after the canonical tables are built and validated and
# the models are fitted), so every later stage and every reader of the archive sees the
# same data.
#
# The archive IS the canonical tables: every data/tables/*.csv is copied byte for byte
# (one copy of each data item; columns, units and keys are documented in
# code/tables/dictionary*.csv, from which 03_generate_data_readme.R writes the README).
# Nothing is joined, recomputed or filled here. Added alongside them:
#   methanotroph_definitions.csv                    the methanotroph taxon list the 16S
#     analyses use (code/lib/methanotroph_definitions.csv)
# Results computed later are archived separately by 02_compile_results.R.
#
# Replaces (2026-10-16) the archive built from data/processed/: flux_stem.csv and
# flux_soil.csv with model covariates joined by campaign + flux value and a hard-coded
# qc_pass; tree_properties.csv (the merged per-tree file, with its own species and DBH
# rules: now trees.csv + wood_cores.csv + soil_readings.csv + internal_gas.csv);
# ddpcr_gene_abundances.csv (now ddpcr.csv + samples.csv); otu_table_16S.csv,
# taxonomy_key_16S.csv, sample_metadata_16S.csv (now asv_counts_16S*.csv,
# asv_taxonomy_16S*.csv, samples_16S.csv, qpcr.csv); environmental_timeseries.csv (the
# raw logger file in logger-clock time; now weather.csv, true time, screened);
# black_oak_experiment.csv and black_oak_its_load.csv (the felled oak is a campaign in
# flux_stem, internal_gas, samples, ddpcr and qpcr); forest_inventory_stems.csv (now
# inventory_stems.csv; census-only attributes such as status stay in the census sheet in
# data/raw/inventory/).
# ==============================================================================
out_dir <- "data/compiled"
tables_dir <- "data/tables"
dir.create(out_dir, showWarnings = FALSE)

superseded <- c("tree_properties.csv", "ddpcr_gene_abundances.csv", "otu_table_16S.csv",
                "taxonomy_key_16S.csv", "sample_metadata_16S.csv", "environmental_timeseries.csv",
                "black_oak_experiment.csv", "black_oak_its_load.csv", "forest_inventory_stems.csv",
                "semirigid_chamber_flux.csv", "static_chamber_flux.csv", "height_chamber_flux.csv",
                "flux_measurements_tree.csv", "flux_measurements_soil.csv", "picrust_pathway_associations.csv")
unlink(file.path(out_dir, superseded))          # superseded files must not linger in the archive

src <- c(list.files(tables_dir, "[.]csv$", full.names = TRUE),
         "code/lib/methanotroph_definitions.csv")
miss <- src[!file.exists(src)]
if (length(miss)) stop("missing input(s) -- run stage B first:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)
ok <- file.copy(src, file.path(out_dir, basename(src)), overwrite = TRUE, copy.date = FALSE)
stopifnot(all(ok))
same <- unname(tools::md5sum(src)) == unname(tools::md5sum(file.path(out_dir, basename(src))))
stopifnot(all(same))
for (f in basename(src)) {
  n <- length(count.fields(file.path(out_dir, f), sep = ",", quote = "\"")) - 1L
  cat(sprintf("  %-34s %7d rows\n", f, n))
}
stale <- setdiff(list.files(out_dir, "[.]csv$"), basename(src))
# the archive mirrors data/tables exactly: anything else is a leftover and is removed
if (length(stale)) { cat("\nremoved from", out_dir, "(no longer a table):", paste(stale, collapse = ", "), "\n"); unlink(file.path(out_dir, stale)) }
cat(sprintf("\n%d datasets in %s/ (%d canonical tables)\n", length(src), out_dir,
            sum(dirname(src) == tables_dir)))
