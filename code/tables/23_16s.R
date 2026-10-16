#!/usr/bin/env Rscript
source("code/lib/outputs.R")
# ==============================================================================
# 23_16s.R -- 16S amplicon tables (data/tables/)
# ------------------------------------------------------------------------------
# samples_16S              one row per sequenced library, both runs
# asv_counts_16S           main run: ASV x library read counts (wide; feature_id +
#                          one column per library_id), as data/raw/16s/OTU_table.txt
# asv_taxonomy_16S         main run: SILVA 138 ranks per ASV (OTU_table.txt columns
#                          Kingdom..Species, as code/lib/build_phyloseq.R uses them)
# asv_counts_16S_felled_oak / asv_taxonomy_16S_felled_oak
#                          the felled-oak run (data/raw/16s/black_oak/OTU_table.txt),
#                          denoised separately: its ASV ids (seq1, seq2, ...) reuse the
#                          main run's ids for DIFFERENT sequences. Never join the two runs
#                          on feature_id.
#
# Library -> sample (main run): data/raw/16s/tree_16s_mapping_dada2_corrected.txt
# (589 libraries = the 589 count columns), with the seq_id patches of
# code/lib/build_phyloseq.R: prime -> ', star -> *, HM -> H, and library
# 206.B01.16S.S (fastq 206_B01_16S_S3) is RO104 Inner. A library carrying a lab ID
# is matched to the sample with that lab ID; otherwise by tree (NCBI spelling) +
# compartment (+ lab-ID suffix, e.g. YB9Inner44c). library_type: mock (Zymo),
# other (libraries labelled with no survey sample: muck, aquatic, microbial goo,
# inside rot, "WMAVan", "Unknown ACSA SM1 or ACRU RMS1", and 202.C02 -- sample_id
# NA), sample (the rest). No extraction or PCR blanks are among the 589 libraries.
# Felled-oak run: library names are <sample>.16S.S<n>; no fastq map or run date in
# the repository.
# Writes data/tables/{samples_16S, asv_counts_16S, asv_taxonomy_16S,
#        asv_counts_16S_felled_oak, asv_taxonomy_16S_felled_oak}.csv,
#        outputs/audit/samples_16S_build.txt
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/tree_ids.R")
OUT <- "data/tables"
SAMPLES <- read.csv(file.path(OUT, "samples.csv"), stringsAsFactors = FALSE)
RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

to_sample_id <- function(x) {   # tree spelling + compartment [+ lab ID] -> samples.sample_id
  m <- str_match(x, "^(.*?)(Inner|Outer|Mineral|Organic)([0-9]+[a-e])?$")
  tree <- gsub("star", "*", gsub("prime", "'", gsub("slash", "/", m[, 2])))
  tid <- ifelse(tree == "RO7/8", "RO7slash8", ncbi_tree_id(tree))          # soil between RO7 and RO8
  y <- paste0(tid, m[, 3], coalesce(m[, 4], ""))
  ifelse(!is.na(m[, 1]) & y %in% SAMPLES$sample_id, y, NA_character_)
}

# ==============================================================================
# 1. MAIN RUN
# ==============================================================================
M <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", check.names = FALSE, colClasses = "character")
# seq_id patches (code/lib/build_phyloseq.R)
M$seq_id <- sub("prime", "'", M$seq_id); M$seq_id <- sub("star", "*", M$seq_id); M$seq_id <- sub("HM", "H", M$seq_id)
fix <- M$ForwardFastqFile == "206_B01_16S_S3_R1_001.fastq"
stopifnot(sum(fix) == 1); M$seq_id[fix] <- "RO104"; M$core_type[fix] <- "Inner"

S21 <- SAMPLES %>% filter(campaign == "2021_survey")
L1 <- M %>%
  mutate(library_id = `#SampleID`,
         by_lab = S21$sample_id[match(lab_id, S21$lab_id)],
         by_name = to_sample_id(paste0(seq_id, core_type)),
         sample_id = coalesce(by_lab, by_name),
         library_type = case_when(grepl("^Zymo", library_id) ~ "mock",
                                  is.na(sample_id) ~ "other",
                                  TRUE ~ "sample"),
         label = sub("\\.16S\\..*$", "", library_id),
         date_seq = as.Date(substr(Date_Seq, 1, 10))) %>%
  transmute(library_id, run = "main", sample_id, library_type, label, fastq_forward = ForwardFastqFile,
            fastq_reverse = ReverseFastqFile, date_seq,
            note = case_when(fix ~ "fastq 206_B01_16S_S3 reassigned to RO104 Inner (code/lib/build_phyloseq.R)",
                             !is.na(by_lab) & !is.na(by_name) & by_lab != by_name ~ paste("lab ID gives", by_lab, "; name gives", by_name),
                             TRUE ~ NA_character_))
stopifnot(nrow(L1) == 589, !anyDuplicated(L1$library_id))
lab_name_conflict <- L1 %>% filter(grepl("^lab ID gives", note))

# counts and taxonomy
O <- read_tsv("data/raw/16s/OTU_table.txt", col_types = cols(.default = "c"), progress = FALSE, name_repair = "minimal")
names(O)[1] <- "feature_id"
stopifnot(ncol(O) == 597, identical(names(O)[591:597], RANKS), setequal(names(O)[2:590], L1$library_id),
          !anyDuplicated(O$feature_id))
CNT <- O[, 1:590] %>% mutate(across(-feature_id, as.integer))
stopifnot(!anyNA(CNT))
TAX <- O[, c(1, 591:597)] %>% rename_with(tolower, all_of(RANKS)) %>% mutate(across(-feature_id, ~ na_if(.x, "")))
# taxonomy_table.txt (the SILVA 138 assignment as strings) agrees, "none" = blank
TX <- read_tsv("data/raw/16s/taxonomy_table.txt", col_types = cols(.default = "c"), progress = FALSE, name_repair = "minimal")
tx <- str_split_fixed(TX$Taxon, ";\\s*", 7); tx[tx %in% c("none", "")] <- NA
tax_mismatch <- sum(TX$`Feature ID` != TAX$feature_id) +
  sum(sapply(1:7, function(j) sum(!(is.na(tx[, j]) & is.na(TAX[[j + 1]])) & (is.na(tx[, j]) | is.na(TAX[[j + 1]]) | tx[, j] != TAX[[j + 1]]))))

# ==============================================================================
# 2. FELLED-OAK RUN
# ==============================================================================
OB <- read_tsv("data/raw/16s/black_oak/OTU_table.txt", col_types = cols(.default = "c"), progress = FALSE, name_repair = "minimal")
names(OB)[1] <- "feature_id"
nb <- ncol(OB) - 7
stopifnot(identical(names(OB)[(nb + 1):ncol(OB)], RANKS), !anyDuplicated(OB$feature_id))
CNTB <- OB[, 1:nb] %>% mutate(across(-feature_id, as.integer))
TAXB <- OB[, c(1, (nb + 1):ncol(OB))] %>% rename_with(tolower, all_of(RANKS)) %>% mutate(across(-feature_id, ~ na_if(.x, "")))
L2 <- tibble(library_id = names(OB)[2:nb]) %>%
  mutate(sample_id = sub("\\.16S\\..*$", "", library_id), run = "felled_oak", library_type = "sample", label = sample_id,
         fastq_forward = NA_character_, fastq_reverse = NA_character_, date_seq = as.Date(NA), note = NA_character_)
stopifnot(all(L2$sample_id %in% SAMPLES$sample_id[SAMPLES$campaign == "2022_felled_oak"]))

# ==============================================================================
# 3. WRITE
# ==============================================================================
LIB <- bind_rows(L1, L2) %>% select(library_id, run, sample_id, library_type, label, fastq_forward, fastq_reverse, date_seq, note) %>%
  arrange(run, library_id)
stopifnot(!anyDuplicated(LIB$library_id))
write_csv(LIB, file.path(OUT, "samples_16S.csv"), na = "")
write_csv(CNT, file.path(OUT, "asv_counts_16S.csv"), na = "")
write_csv(TAX, file.path(OUT, "asv_taxonomy_16S.csv"), na = "")
write_csv(CNTB, file.path(OUT, "asv_counts_16S_felled_oak.csv"), na = "")
write_csv(TAXB, file.path(OUT, "asv_taxonomy_16S_felled_oak.csv"), na = "")

# ==============================================================================
# 4. AUDIT
# ==============================================================================
reads <- colSums(as.matrix(CNT[, -1])); readsB <- colSums(as.matrix(CNTB[, -1]))
LIB$reads <- c(reads, readsB)[LIB$library_id]
shared <- intersect(CNT$feature_id, CNTB$feature_id)
sink(out_path("samples_16S_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("samples_16S %d libraries | main run: %d ASVs x %d libraries, %d reads | felled oak: %d ASVs x %d libraries, %d reads\n",
            nrow(LIB), nrow(CNT), ncol(CNT) - 1, sum(reads), nrow(CNTB), ncol(CNTB) - 1, sum(readsB)))
print(as.data.frame(LIB %>% count(run, library_type, date_seq)), row.names = FALSE)
cat("\nLibraries by sample material (run x material x position):\n")
print(as.data.frame(LIB %>% left_join(SAMPLES %>% select(sample_id, material, position), by = "sample_id") %>%
                      count(run, material, position)), row.names = FALSE)
cat("\nLibraries of type other / mock (sample_id NA):\n")
print(as.data.frame(LIB %>% filter(library_type != "sample") %>% select(library_id, library_type, label, reads)), row.names = FALSE)
cat(sprintf("\nSamples sequenced more than once: %d\n", sum(duplicated(na.omit(LIB$sample_id)))))
cat(sprintf("2021 samples without a library: %d of %d\n", sum(!S21$sample_id %in% LIB$sample_id), nrow(S21)))
cat(sprintf("Library lab ID and name disagree: %d\n", nrow(lab_name_conflict)))
if (nrow(lab_name_conflict)) print(as.data.frame(lab_name_conflict %>% select(library_id, note)), row.names = FALSE)
cat(sprintf("Taxonomy: OTU_table.txt ranks vs taxonomy_table.txt strings: %d differences\n", tax_mismatch))
cat(sprintf("feature_id values present in both runs: %d (different sequences; the runs were denoised separately)\n", length(shared)))
cat("\nNote: read.delim(check.names = TRUE) in code/lib/build_phyloseq.R renames the count columns 202.C02.16S.S and\n",
    "206.B01.16S.S to X202... / X206..., so phyloseq keeps 587 libraries and drops these two -- including the RO104 Inner\n",
    "library that the fastq patch reassigns.\n", sep = "")
options(op); sink()
cat(readLines(out_path("samples_16S_build.txt", "audit")), sep = "\n")
