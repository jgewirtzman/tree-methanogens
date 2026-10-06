#!/usr/bin/env Rscript
# ==============================================================================
# 23_si-table-workbook.R -- the data-heavy SI tables as one Excel workbook
# ------------------------------------------------------------------------------
# Tables S2 (methane-cycling taxa and their classification) and S6 (diameter of
# measured trees by campaign, location, species and status) are long reference
# tables, so they are supplied as a workbook a reader can sort and filter rather
# than as pages of print. Each sheet carries its own legend above the table. The
# short tables (S1, S3–S5, S7) are formatted in the SI tables document.
#
# Reads:  outputs/data/known_putative_taxa_table.csv (04_known-putative-table.R)
#         outputs/data/dbh_by_species_campaign.csv   (08_dbh_by_species_campaign.R)
# Writes: outputs/tables/SI_Tables_S2_S6.xlsx
# ==============================================================================
suppressPackageStartupMessages({ library(openxlsx) })
source("code/lib/outputs.R")

S2 <- read.csv(out_path("known_putative_taxa_table.csv"), stringsAsFactors = FALSE, check.names = FALSE)
S2$note <- sub("^Capacity rule [0-9-]+: *", "", S2$note)            # internal dating, not for readers
S2$classification <- sub("^Methanotroph_", "Methanotroph, ", S2$classification)
S2[is.na(S2)] <- ""
names(S2) <- c("Class", "Family", "Taxon", "Level", "ASVs", "Share of all 16S reads, pooled (%)", "Source", "Note")

S4 <- read.csv(out_path("dbh_by_species_campaign.csv"), stringsAsFactors = FALSE, check.names = FALSE)
names(S4) <- c("Campaign", "Location", "Species", "Status", "Trees", "DBH (cm), mean ± SD")

LEG <- list(
  `Table S2` = c("Table S2. Methane-cycling taxa in the 16S data and their classification.",
    paste("Methanotrophs are classified by methane monooxygenase capacity. Known: genera that carry particulate or soluble",
          "methane monooxygenase, including the NC10 genus Candidatus Methylomirabilis. Putative: an upper bound -- members of",
          "methanotroph-containing families unresolved at genus level, and genera in those families whose capacity",
          "is untested (no sequenced genome) or found in only some genomes (Methylovirgula, Rhodoblastus; Table S3). Genera inside those families with no",
          "annotated methane monooxygenase among their NCBI protein records are listed in Table S3 and not counted.",
          "ASVs: number of amplicon sequence variants in the full, unrarefied 16S table (all libraries); share: the taxon's reads as a percentage of all 16S reads pooled across libraries.",
          "Counts differ from Methods S3, which counts ASVs present in survey samples after rarefaction to 3,500 reads.")),
  `Table S6` = c("Table S6. Diameter at breast height of measured trees by campaign, location, species and status.",
    "Trees: number of individual trees; DBH: mean ± standard deviation across those trees (a single value where n = 1)."))

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", border = "bottom")
ls <- createStyle(wrapText = TRUE, valign = "top")
for (nm in names(LEG)) {
  d <- if (nm == "Table S2") S2 else S4
  addWorksheet(wb, nm)
  writeData(wb, nm, LEG[[nm]][1], startRow = 1); addStyle(wb, nm, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  writeData(wb, nm, LEG[[nm]][2], startRow = 2)
  mergeCells(wb, nm, cols = 1:ncol(d), rows = 2); addStyle(wb, nm, ls, rows = 2, cols = 1)
  setRowHeights(wb, nm, rows = 2, heights = 75)
  writeData(wb, nm, d, startRow = 4, headerStyle = hs)
  setColWidths(wb, nm, cols = 1:ncol(d), widths = "auto")
  freezePane(wb, nm, firstActiveRow = 5)
}
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
OUT_XLSX <- "outputs/tables/SI_Tables_S2_S6.xlsx"
saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)

# openxlsx 4.2.8 writes relationships to drawing parts it never adds to the archive
# (xl/drawings/drawing1.xml, vmlDrawing1.vml), so Excel offers to repair the file and
# strict readers fail. This workbook has no drawings: drop the dangling references.
fix_dangling_drawings <- function(path) {
  td <- tempfile("xlsx_"); dir.create(td); unzip(path, exdir = td)
  have <- function(rel, base) file.exists(file.path(td, "xl", base, rel))
  for (rf in list.files(file.path(td, "xl/worksheets/_rels"), full.names = TRUE)) {
    x <- paste(readLines(rf, warn = FALSE), collapse = "")
    rels <- regmatches(x, gregexpr("<Relationship [^>]*/>", x))[[1]]
    keep <- vapply(rels, function(r) {
      tgt <- sub('.*Target="([^"]+)".*', "\\1", r)
      !grepl("drawing", tgt, ignore.case = TRUE) || file.exists(file.path(td, "xl/worksheets", tgt))
    }, logical(1))
    for (r in rels[!keep]) {
      id <- sub('.*Id="([^"]+)".*', "\\1", r)
      x <- sub(r, "", x, fixed = TRUE)
      sh <- file.path(td, "xl/worksheets", sub("\\.rels$", "", basename(rf)))
      y <- paste(readLines(sh, warn = FALSE), collapse = "")
      y <- gsub(sprintf('<(drawing|legacyDrawing) r:id="%s"/>', id), "", y)
      writeLines(y, sh, sep = "")
    }
    writeLines(x, rf, sep = "")
  }
  ct <- file.path(td, "[Content_Types].xml"); z <- paste(readLines(ct, warn = FALSE), collapse = "")
  ov <- regmatches(z, gregexpr('<Override PartName="/xl/drawings/[^"]+"[^>]*/>', z))[[1]]
  for (o in ov) if (!file.exists(file.path(td, sub('.*PartName="/([^"]+)".*', "\\1", o)))) z <- sub(o, "", z, fixed = TRUE)
  writeLines(z, ct, sep = "")
  files <- list.files(td, recursive = TRUE, all.files = TRUE)
  files <- c("[Content_Types].xml", setdiff(files, "[Content_Types].xml"))
  abs <- file.path(normalizePath(dirname(path)), basename(path)); unlink(abs)
  zip::zip(abs, files = files, root = td, mode = "mirror")
  unlink(td, recursive = TRUE)
}
fix_dangling_drawings(OUT_XLSX)
cat(sprintf("Wrote outputs/tables/SI_Tables_S2_S6.xlsx (Table S2: %d rows; Table S6: %d rows)\n", nrow(S2), nrow(S4)))
