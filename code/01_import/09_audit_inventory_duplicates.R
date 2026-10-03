# ==============================================================================
# 09_audit_inventory_duplicates.R -- are any trees counted in both the 2018 and 2019
# censuses?
# ------------------------------------------------------------------------------
# The two censuses split 14 quadrats along the east edge of the plot. A 2018 stem
# with a 2019 stem of the same species nearby and a similar diameter looks like a
# double count -- but in dense sapling patches such pairs occur by chance, and the
# 2019 census itself records multi-stem sprouts metres apart. So each rule is scored
# against a baseline: how often does a 2019 stem match ANOTHER 2019 stem by it?
# A rule is evidence of duplication only where 2018->2019 matches exceed that baseline.
#
# Result 2026-09-30: they never do. Among stems >= 5 cm (baseline 0.7%), 0 of 384
# 2018 stems match a 2019 stem; small-stem matches (4%) are well below baseline (42%).
# 08_inventory_build.R therefore drops nothing across censuses.
#
# Run from the repo root after 08_inventory_build.R:
#   Rscript code/01_import/09_audit_inventory_duplicates.R
# Writes outputs/audit/inventory_duplicate_test.txt
# ==============================================================================
source("code/lib/outputs.R")
suppressMessages(library(dplyr))
inv <- read.csv("outputs/tables/inventory_stems.csv", stringsAsFactors = FALSE) %>% filter(located)
B <- inv %>% filter(source == "bytag") %>% transmute(sp = species_code, dbh = dbh_cm, PX, PY)
F <- inv %>% filter(source == "fg19", PX >= 95, PY >= 95) %>%   # the block the 2018 census covers (+5 m)
  transmute(sp = species_code, dbh = dbh_cm, PX, PY)

hit <- function(px, py, sp, dbh, P, r, self = FALSE) {
  if (dbh < r$min_dbh) return(FALSE)
  d <- sqrt((P$PX - px)^2 + (P$PY - py)^2); j <- which(d < r$dmax & (!self | d > 0))
  any(P$sp[j] == sp & abs(P$dbh[j] - dbh) <= r$rel * dbh)
}
RULES <- list(list(dmax = 1, min_dbh = 0, rel = .15), list(dmax = .5, min_dbh = 0, rel = .10),
              list(dmax = 1, min_dbh = 5, rel = .15), list(dmax = .5, min_dbh = 5, rel = .10),
              list(dmax = 1, min_dbh = 10, rel = .15))
sink(out_path("inventory_duplicate_test.txt"))
cat("CROSS-CENSUS DUPLICATE TEST (code/01_import/09_audit_inventory_duplicates.R)\n\n")
for (r in RULES) {
  m18 <- mapply(hit, B$PX, B$PY, B$sp, B$dbh, MoreArgs = list(P = F, r = r))
  m19 <- mapply(hit, F$PX, F$PY, F$sp, F$dbh, MoreArgs = list(P = F, r = r, self = TRUE))
  n18 <- sum(B$dbh >= r$min_dbh); n19 <- sum(F$dbh >= r$min_dbh)
  cat(sprintf("<%.1f m, DBH >= %2.0f cm, within %2.0f%%: 2018->2019 %3d of %4d (%4.1f%%) | 2019->2019 baseline %4.1f%%\n",
              r$dmax, r$min_dbh, 100 * r$rel, sum(m18), n18, 100 * sum(m18) / n18, 100 * sum(m19) / n19))
}
cat("\nDuplication is indicated only where the 2018->2019 rate exceeds the baseline.\n")
sink()
cat(readLines(out_path("inventory_duplicate_test.txt")), sep = "\n")
