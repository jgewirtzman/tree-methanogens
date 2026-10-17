# ==============================================================================
# 24_picrust.R -- canonical tables: PICRUSt2 pathway abundances per 16S library
# ------------------------------------------------------------------------------
# picrust_pathways       library_id x pathway: summed PICRUSt2 contributions of the
#                        library's non-organelle ASVs (abundance), and the same without
#                        the ASVs predicted to carry mcrA or pmoA
# picrust_pathway_info   pathway, description
# Built in stage A by code/01_import/06_picrust_pathways.R from the PICRUSt2 predictions
# (data/raw/picrust/, the large Zenodo download) into data/processed/molecular/picrust/.
# This step checks them against the 16S library table and copies them byte for byte
# (no value is re-parsed or rewritten).
# Writes data/tables/{picrust_pathways, picrust_pathway_info}.csv
# ==============================================================================
SRC <- "data/processed/molecular/picrust"; OUT <- "data/tables"
f <- file.path(SRC, c("picrust_pathways.csv", "picrust_pathway_info.csv"))
if (!all(file.exists(f)))
  stop("missing ", paste(f[!file.exists(f)], collapse = ", "),
       ": from the Zenodo processed archive, or rebuilt by code/01_import/06_picrust_pathways.R",
       " (stage A, needs data/raw/picrust/)", call. = FALSE)
P <- utils::read.csv(f[1], colClasses = c("character", "character", "numeric", "numeric", "numeric"))
I <- utils::read.csv(f[2], colClasses = "character")
L <- utils::read.csv(file.path(OUT, "samples_16S.csv"), colClasses = "character", na.strings = "")
stopifnot(identical(names(P), c("library_id", "pathway", "abundance", "abundance_no_mcra_asvs", "abundance_no_pmoa_asvs")),
          identical(names(I), c("pathway", "description")),
          all(P$library_id %in% L$library_id[L$run == "main"]), all(P$pathway %in% I$pathway),
          !anyDuplicated(P[, c("library_id", "pathway")]), !anyDuplicated(I$pathway),
          !anyNA(P[, 3:5]), all(P$abundance_no_mcra_asvs <= P$abundance * (1 + 1e-12)),
          all(P$abundance_no_pmoa_asvs <= P$abundance * (1 + 1e-12)))
stopifnot(all(file.copy(f, file.path(OUT, basename(f)), overwrite = TRUE, copy.date = FALSE)))
lt <- L$library_type[match(unique(P$library_id), L$library_id)]
cat(sprintf("picrust_pathways: %d rows, %d libraries (%s) x %d pathways; picrust_pathway_info: %d pathways\n",
            nrow(P), length(lt), paste(names(table(lt)), table(lt), collapse = ", "), length(unique(P$pathway)), nrow(I)))
