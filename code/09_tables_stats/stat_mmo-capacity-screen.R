source("code/lib/outputs.R")
# ==============================================================================
# stat_mmo-capacity-screen.R
# ------------------------------------------------------------------------------
# Genome evidence for methane monooxygenase (MMO) capacity in the resolved genera
# that sit inside mixed methanotroph families. Supports the capacity rule set
# (manuscript OPEN_DECISIONS, 2026-09-30): a taxon counts as a methanotroph if it
# carries pMMO or sMMO genes.
#
#   1. per genus: NCBI Assembly count, and NCBI Protein records named as copper-MMO
#      (pMMO/AMO) or soluble di-iron MMO family members
#   2. each named hit aligned (BLOSUM62, local) to M. capsulatus PmoA (Q607G3) and
#      M. trichosporium MmoX (P27353); a hit counts as a true subunit at >=40% id over
#      >=150 aa to PmoA, or >=55% id over >=300 aa to MmoX. Names alone are unreliable:
#      RefSeq's "ammonia monooxygenase" also labels an unrelated membrane protein.
#   3. in silico match of the mmoX primers 536f/898r (Fuse et al. 1998) to any
#      di-iron hydroxylase found, with true mmoX as the positive control
#
# LIVE: queries NCBI E-utilities, so results drift as NCBI grows. The evidence the
# definitions file cites is the committed snapshot code/lib/mmo_capacity_screen_2026-09-30.csv.
# Not part of run_all.R (network). Needs the Bioconductor package pwalign.
#
# Output: outputs/data/mmo_capacity_screen.csv
#
# STATUS 2026-09-30: executed with pwalign 1.2.0 (Bioconductor 3.20); reproduces every
# verdict in the snapshot and the primer mismatches (0/0 true mmoX; 9/8 AYO81047).
# Identities are pwalign PID1 (gaps count against identity), so they read lower than the
# Biopython figures quoted in the snapshot; the verdicts do not change.
# ==============================================================================
suppressMessages({library(jsonlite); library(Biostrings); library(pwalign)})
E  <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
eu <- function(path, ...) { Sys.sleep(0.4); paste(readLines(url(paste0(E, path, "?",
        paste(names(list(...)), vapply(lapply(list(...), as.character), URLencode, "", reserved = TRUE), sep = "=", collapse = "&"))),
        warn = FALSE), collapse = "\n") }
esearch <- function(db, term, retmax = 200) fromJSON(eu("/esearch.fcgi", db = db, term = term,
                                                         retmode = "json", retmax = retmax))$esearchresult
fasta_aa <- function(ids) { f <- tempfile(); writeLines(eu("/efetch.fcgi", db = "protein",
                              id = paste(ids, collapse = ","), rettype = "fasta", retmode = "text"), f); readAAStringSet(f) }
uniprot  <- function(acc) { f <- tempfile(); download.file(sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", acc), f, quiet = TRUE); readAAStringSet(f)[[1]] }

REF <- list(PmoA = uniprot("Q607G3"), MmoX = uniprot("P27353"))
ident <- function(q, r) {
  a <- pairwiseAlignment(q, r, substitutionMatrix = "BLOSUM62", type = "local",
                         gapOpening = 10, gapExtension = 0.5)
  c(pid = pid(a, type = "PID1"), len = nchar(pattern(a)))
}
CUMMO <- '("methane monooxygenase"[Title] OR "ammonia monooxygenase"[Title] OR pmoA[Gene Name] OR amoA[Gene Name])'
SDIMO <- '("aromatic/alkene monooxygenase hydroxylase"[Title] OR "soluble methane monooxygenase"[Title] OR mmoX[Gene Name])'
GENERA <- c("Methylocapsa","Methylocella",                 # positive controls
            "Roseiarcus","Bosea","Microvirga","Methylovirgula","Methylorosula","Rhodoblastus",
            "Psychroglaciecola","Lichenibacterium","Methylobacterium","Methylorubrum")

rows <- lapply(GENERA, function(g) {
  asm <- as.integer(esearch("assembly", sprintf("%s[Organism]", g))$count)
  ids <- unique(c(esearch("protein", sprintf("%s[Organism] AND %s", g, CUMMO))$idlist,
                  esearch("protein", sprintf("%s[Organism] AND %s", g, SDIMO))$idlist))
  best <- c(PmoA = 0, MmoX = 0); true_hit <- FALSE
  if (length(ids)) {
    s <- fasta_aa(head(ids, 60))
    for (k in seq_along(s)) {
      p <- ident(s[[k]], REF$PmoA); m <- ident(s[[k]], REF$MmoX)
      # identity only counts over a meaningful alignment: very short local alignments
      # (e.g. 12 aa) reach high % identity by chance and would mislead
      if (p["len"] >= 50) best["PmoA"] <- max(best["PmoA"], p["pid"])
      if (m["len"] >= 50) best["MmoX"] <- max(best["MmoX"], m["pid"])
      if ((p["pid"] >= 40 && p["len"] >= 150) || (m["pid"] >= 55 && m["len"] >= 300)) true_hit <- TRUE
    }
  }
  data.frame(genus = g, assemblies = asm, named_hits = length(ids),
             best_pid_PmoA_50aa = round(best["PmoA"], 1), best_pid_MmoX_50aa = round(best["MmoX"], 1),
             verified_mmo = if (asm == 0) "untested (no genome)" else if (true_hit) "yes" else "no")
})
R <- do.call(rbind, rows); rownames(R) <- NULL
write.csv(R, out_path("mmo_capacity_screen.csv"), row.names = FALSE)
print(R, row.names = FALSE)

# ---- mmoX primer match to the Methylobacterium di-iron hydroxylase ------------
primer_mm <- function(primer, target) {
  hits <- lapply(list(DNAString(primer), reverseComplement(DNAString(primer))), function(p)
    min(vapply(seq_len(nchar(target) - nchar(p) + 1), function(i)
      sum(strsplit(as.character(p), "")[[1]] != strsplit(as.character(subseq(target, i, i + nchar(p) - 1)), "")[[1]]), 0)))
  min(unlist(hits))
}
cds <- function(protein_acc) { f <- tempfile(); writeLines(eu("/efetch.fcgi", db = "protein", id = protein_acc,
                                 rettype = "fasta_cds_na", retmode = "text"), f); readDNAStringSet(f)[[1]] }
gb  <- tempfile(); writeLines(eu("/efetch.fcgi", db = "nuccore", id = "X55394", rettype = "fasta", retmode = "text"), gb)
mmo_operon <- readDNAStringSet(gb)[[1]]              # contains M. trichosporium mmoX (positive control)
P <- c(`536f` = "CGCTGTGGAAGGGCATGAAGCG", `898r` = "GCTCGACCTTGAACTTGGAGCC")
cat("\nmmoX primer mismatches (best site, either strand):\n")
for (nm in names(P)) cat(sprintf("  %s  M. trichosporium mmo operon: %d   M. brachiatum AYO81047 CDS: %d\n",
                                 nm, primer_mm(P[[nm]], mmo_operon), primer_mm(P[[nm]], cds("AYO81047"))))
