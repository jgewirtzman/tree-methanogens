# ==============================================================================
# tree_ids.R -- spelling rules for tree names (sourced by code/tables/)
# ------------------------------------------------------------------------------
# 2021 survey trees (and the 2024 heartwood/sapwood resample) are named in the
# field ("AB 3/5", "AB6'", "WA103*"). Their tree_id is the NCBI spelling: upper
# case, no spaces, special characters spelled out -- "/" -> slash, prime -> prime,
# "*" -> star ("?" -> query, for the two 2024 names recorded with a question mark).
#   AB 3/5 -> AB3slash5    AB6' -> AB6prime    WA103* -> WA103star
# Two further rules: hemlock names written HM<n> are H<n> (HM7 = H7; the 16S map
# and the 2024 gas run use HM), and leading zeros are dropped (WP06 -> WP6).
# Whole-name corrections (typos, nicknames) are in code/lib/tree_name_aliases.csv.
# ==============================================================================

TREE_ALIASES_FILE <- "code/lib/tree_name_aliases.csv"

ncbi_tree_id <- function(x, aliases = utils::read.csv(TREE_ALIASES_FILE, stringsAsFactors = FALSE)) {
  y <- toupper(gsub("\\s+", "", trimws(as.character(x))))
  y <- gsub("['’`]", "prime", y)
  y <- gsub("PRIME", "prime", y)
  y <- gsub("\\*", "star", y)
  y <- sub("STAR$", "star", y)
  y <- gsub("/", "slash", y)
  y <- gsub("\\?+", "query", y)
  y <- sub("^([A-Z]+)0+([1-9])", "\\1\\2", y)
  y <- sub("^HM(?=[0-9XZ])", "H", y, perl = TRUE)
  hit <- y %in% aliases$spelling
  y[hit] <- aliases$tree_id[match(y[hit], aliases$spelling)]
  y[is.na(x) | trimws(x) == ""] <- NA_character_
  y
}

# Species from the field-name prefix (2021 and 2024 surveys); the fallback when a
# record carries no species code.
FIELD_PREFIX_SPECIES <- c(AB = "FAGR", BB = "BELE", BC = "PRSE", BO = "QUVE", H = "TSCA",
                          ML = "KALA", PB = "BEPA", RM = "ACRU", RO = "QURU", SA = "SAAL",
                          SH = "CAOV", SM = "ACSA", WA = "FRAM", WO = "QUAL", WP = "PIST",
                          YB = "BEAL", CUL8R = "BEAL")
prefix_species <- function(id) {
  p <- ifelse(id %in% names(FIELD_PREFIX_SPECIES), id,
       ifelse(substr(id, 1, 2) %in% names(FIELD_PREFIX_SPECIES), substr(id, 1, 2),
       ifelse(grepl("^H", id), "H", NA_character_)))
  unname(FIELD_PREFIX_SPECIES[p])
}

# Species-code corrections for data-entry typos (not species).
SPECIES_CODE_FIXES <- c(TSLA = "TSCA")
fix_species_code <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA")] <- NA_character_
  ifelse(x %in% names(SPECIES_CODE_FIXES), SPECIES_CODE_FIXES[x], x)
}
