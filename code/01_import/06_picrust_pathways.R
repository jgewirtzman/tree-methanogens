# ==============================================================================
# 06_picrust_pathways.R -- PICRUSt2 MetaCyc pathway abundances per 16S library
# ------------------------------------------------------------------------------
# Stage A; workdir: repository root. Turns the PICRUSt2 predictions
# (data/raw/picrust/, the 7.9 GB tree-methanogens_picrust.zip) into the one small
# input every PICRUSt2 analysis and figure reads. Nothing analytical happens here:
# the contributions are summed, no sample is selected and nothing is normalised.
#
# Rule. A pathway's abundance in a library is the sum of its per-ASV contributions
# (path_abun_contrib.tsv taxon_function_abun) over the library's ASVs, leaving out
# organelle ASVs (is_organelle(), code/lib/build_phyloseq.R: SILVA family
# Mitochondria or order Chloroplast, the rule of every 16S analysis). Two further
# versions leave out, in addition, the ASVs PICRUSt2 predicts to carry
#   mcrA (KEGG K00399, KO_predicted.tsv > 0)  -> abundance_no_mcra_asvs
#   pmoA (KEGG K10944, KO_predicted.tsv > 0)  -> abundance_no_pmoa_asvs
# so an association can be tested without the gene carriers' own contribution.
#
# Reads   data/raw/picrust/path_abun_contrib.tsv     per-ASV pathway contributions
#         data/raw/picrust/KO_predicted.tsv          KEGG orthologs per ASV
#         data/raw/picrust/path_abun_unstrat_descrip.tsv  pathway descriptions
#         data/raw/16s/OTU_table.txt                 SILVA ranks per ASV (organelles)
# Writes  data/processed/molecular/picrust/picrust_pathways.csv
#             library_id x pathway (rows with a non-zero abundance in any version)
#         data/processed/molecular/picrust/picrust_pathway_info.csv   pathway, description
#         outputs/audit/picrust_pathways_build.txt
# code/tables/24_picrust.R checks them against samples_16S and makes them the canonical
# tables data/tables/picrust_pathways.csv and picrust_pathway_info.csv (as the flux fits
# here in stage A become flux_stem.csv in stage B). They are in the Zenodo processed
# archive, so stages B-G never need data/raw/picrust/.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(readr) })
source("code/lib/build_phyloseq.R")          # is_organelle()
source("code/lib/outputs.R")
RAW <- "data/raw/picrust"
need <- file.path(RAW, c("path_abun_contrib.tsv", "KO_predicted.tsv", "path_abun_unstrat_descrip.tsv"))
if (!all(file.exists(need)))
  stop("PICRUSt2 predictions not found: ", paste(need[!file.exists(need)], collapse = ", "),
       "\n  Download tree-methanogens_picrust.zip from the Zenodo record (README.md) and unzip it",
       " at the repository root (it unpacks into data/raw/picrust/).", call. = FALSE)

# ---- organelle ASVs (main-run taxonomy, as data/tables/asv_taxonomy_16S.csv) -----
RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
O <- fread("data/raw/16s/OTU_table.txt", sep = "\t", colClasses = "character",
           select = c(1, 591:597), na.strings = "")
stopifnot(identical(names(O)[-1], RANKS))
setnames(O, 1, "feature_id")
organelle <- O$feature_id[is_organelle(O)]

# ---- gene carriers (PICRUSt2 KEGG predictions) -------------------------------------
KO <- fread(file.path(RAW, "KO_predicted.tsv"), select = c("sequence", "K00399", "K10944"))
mcra_asvs <- KO$sequence[KO$K00399 > 0]
pmoa_asvs <- KO$sequence[KO$K10944 > 0]

# ---- per-library pathway abundance ----------------------------------------------------
C <- fread(file.path(RAW, "path_abun_contrib.tsv"),
           select = c("sample", "function", "taxon", "taxon_function_abun"))
n_all <- nrow(C)
C <- C[!taxon %in% organelle]
n_org <- n_all - nrow(C)
stopifnot(all(unique(C$taxon) %in% O$feature_id))
setnames(C, "function", "pathway")
# grouped sums in file order (the order fixes the floating-point summation)
A  <- C[, .(abundance = sum(taxon_function_abun)), by = .(sample, pathway)]
NM <- C[!taxon %in% mcra_asvs, .(abundance_no_mcra_asvs = sum(taxon_function_abun)), by = .(sample, pathway)]
NP <- C[!taxon %in% pmoa_asvs, .(abundance_no_pmoa_asvs = sum(taxon_function_abun)), by = .(sample, pathway)]
rm(C); invisible(gc())
P <- merge(merge(A, NM, by = c("sample", "pathway"), all = TRUE), NP, by = c("sample", "pathway"), all = TRUE)
for (v in c("abundance", "abundance_no_mcra_asvs", "abundance_no_pmoa_asvs")) set(P, which(is.na(P[[v]])), v, 0)
P <- P[abundance > 0 | abundance_no_mcra_asvs > 0 | abundance_no_pmoa_asvs > 0]
setnames(P, "sample", "library_id")
setorder(P, library_id, pathway)
stopifnot(!anyDuplicated(P[, .(library_id, pathway)]),
          all(P$abundance_no_mcra_asvs <= P$abundance * (1 + 1e-12)),
          all(P$abundance_no_pmoa_asvs <= P$abundance * (1 + 1e-12)))

# ---- pathway descriptions ----------------------------------------------------------
U <- fread(file.path(RAW, "path_abun_unstrat_descrip.tsv"), select = c("pathway", "description"),
           quote = "")
stopifnot(!anyDuplicated(U$pathway), all(P$pathway %in% U$pathway))
INFO <- U[pathway %in% P$pathway][order(pathway)]

OUT <- "data/processed/molecular/picrust"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
write_csv(P, file.path(OUT, "picrust_pathways.csv"), na = "")
write_csv(INFO, file.path(OUT, "picrust_pathway_info.csv"), na = "")
# what is written reads back to the same numbers (to the last bit or two: decimal text
# round-trips doubles to ~1e-16 relative)
back <- utils::read.csv(file.path(OUT, "picrust_pathways.csv"))
for (v in c("abundance", "abundance_no_mcra_asvs", "abundance_no_pmoa_asvs"))
  stopifnot(identical(back$library_id, P$library_id), max(abs(back[[v]] - P[[v]]) / pmax(P[[v]], 1)) < 1e-14)

sink(out_path("picrust_pathways_build.txt", "audit"))
cat(sprintf("path_abun_contrib.tsv: %d contribution rows; %d from %d organelle ASVs removed (%.1f%%)\n",
            n_all, n_org, length(organelle), 100 * n_org / n_all))
cat(sprintf("gene carriers (KO_predicted.tsv > 0): mcrA K00399 %d ASVs; pmoA K10944 %d ASVs\n",
            length(mcra_asvs), length(pmoa_asvs)))
cat(sprintf("picrust_pathways.csv: %d rows, %d libraries x %d pathways\n",
            nrow(P), uniqueN(P$library_id), uniqueN(P$pathway)))
sink()
cat(readLines(out_path("picrust_pathways_build.txt", "audit")), sep = "\n")
