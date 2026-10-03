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
names(S2) <- c("Class", "Family", "Taxon", "Level", "ASVs", "Mean relative abundance (%)", "Source", "Note")

S4 <- read.csv(out_path("dbh_by_species_campaign.csv"), stringsAsFactors = FALSE, check.names = FALSE)
names(S4) <- c("Campaign", "Location", "Species", "Status", "Trees", "DBH (cm), mean ± SD")

LEG <- list(
  `Table S2` = c("Table S2. Methane-cycling taxa in the 16S data and their classification.",
    paste("Methanotrophs are classified by methane monooxygenase capacity. Known: genera that carry particulate or soluble",
          "methane monooxygenase, including the NC10 genus Candidatus Methylomirabilis. Putative: members of methanotroph-",
          "containing families unresolved at genus level; an upper bound only. Genera inside those families with genome",
          "evidence of no methane monooxygenase (Table S3) are listed in Table S3 and not counted.",
          "ASVs: number of amplicon sequence variants; mean relative abundance across all samples.")),
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
saveWorkbook(wb, "outputs/tables/SI_Tables_S2_S6.xlsx", overwrite = TRUE)
cat(sprintf("Wrote outputs/tables/SI_Tables_S2_S6.xlsx (Table S2: %d rows; Table S6: %d rows)\n", nrow(S2), nrow(S4)))
