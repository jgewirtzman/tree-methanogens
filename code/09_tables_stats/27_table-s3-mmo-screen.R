# ==============================================================================
# 27_table-s3-mmo-screen.R -- Table S3: genome screen for methane monooxygenase
# ------------------------------------------------------------------------------
# Stage D; run from the repository root. The screen itself queries NCBI
# (code/tools/mmo_capacity_screen.R), so its 2026-09-30 run is committed and read
# here; this script adds the classification each genus takes in Table S2.
#
# Inputs
#   code/lib/mmo_capacity_screen_table_S3_2026-09-30.csv   the tool's output
#   code/lib/mmo_capacity_screen_2026-09-30.csv            capacity verdict per genus
#   code/lib/methanotroph_definitions.csv                   Table S2 tiers
# Output
#   outputs/data/table_S3_mmo_capacity_screen.csv
#
# Columns
#   assemblies, named_hits      NCBI assemblies; protein records named as MMO subunits
#   best_pid_PmoA_50aa,         best identity (%) to M. capsulatus PmoA (Q607G3) /
#   best_pid_MmoX_50aa          M. trichosporium MmoX (P27353) among alignments of
#                               >= 50 amino acids
#   subunit_passing_threshold   a record passes >= 40% identity over >= 150 aa to PmoA
#                               or >= 55% over >= 300 aa to MmoX (the caption's rule)
#   capacity_verdict            the verdict after checking the passing records
#                               (positive control; no; yes, one genome; yes,
#                               strain-specific; untested, no genome)
#   table_S2_tier               Known, Putative, or not counted (methanotroph_definitions.csv)
# A passing record is not a verified methane monooxygenase on its own: Rhodoblastus
# (one genome of unconfirmed source) and Methylovirgula (one strain) are Putative.
# ==============================================================================
SCREEN <- read.csv("code/lib/mmo_capacity_screen_table_S3_2026-09-30.csv", stringsAsFactors = FALSE)
VERDICT <- read.csv("code/lib/mmo_capacity_screen_2026-09-30.csv", stringsAsFactors = FALSE)
DEF <- read.csv("code/lib/methanotroph_definitions.csv", stringsAsFactors = FALSE)

# genera recorded under another name in methanotroph_definitions.csv
DEF_TAXON <- c(Lichenibacterium = "1174-901-12", Methylobacterium = "Methylobacterium-Methylorubrum",
               Methylorubrum = "Methylobacterium-Methylorubrum")
taxon <- ifelse(SCREEN$genus %in% names(DEF_TAXON), DEF_TAXON[SCREEN$genus], SCREEN$genus)
d <- DEF[DEF$Taxon_rank == "Genus", ][match(taxon, DEF$Taxon[DEF$Taxon_rank == "Genus"]), ]
stopifnot(!anyNA(d$Taxon), setequal(SCREEN$genus, VERDICT$genus))

verdict <- VERDICT$verified_mmo_capacity[match(SCREEN$genus, VERDICT$genus)]
OUT <- data.frame(
  genus = SCREEN$genus, assemblies = SCREEN$assemblies, named_hits = SCREEN$named_hits,
  best_pid_PmoA_50aa = SCREEN$best_pid_PmoA_50aa, best_pid_MmoX_50aa = SCREEN$best_pid_MmoX_50aa,
  subunit_passing_threshold = SCREEN$verified_mmo,
  capacity_verdict = sub("^untested$", "untested (no genome)", verdict),
  table_S2_tier = ifelse(d$Include_known == "YES", "Known", ifelse(d$Include_putative == "YES", "Putative", "not counted")))
# a genus with a passing record is counted (Known or Putative); one without is not counted
stopifnot(all(OUT$table_S2_tier[OUT$subunit_passing_threshold == "yes"] != "not counted"),
          all(OUT$table_S2_tier[OUT$subunit_passing_threshold == "no"] == "not counted"))

source("code/lib/outputs.R")
write.csv(OUT, out_path("table_S3_mmo_capacity_screen.csv"), row.names = FALSE)
print(OUT[, c("genus", "subunit_passing_threshold", "capacity_verdict", "table_S2_tier")], row.names = FALSE)
cat("\nWritten: outputs/data/table_S3_mmo_capacity_screen.csv\n")
