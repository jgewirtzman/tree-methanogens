# ==============================================================================
# Methanotroph Definitions Loader
# ==============================================================================
# Purpose: Provides shared functions for loading the canonical methanotroph
#   definitions CSV and classifying ASVs/OTUs as Known, Putative, or
#   non-methanotroph. Used by figure scripts 08b, 08c, and 10.
#
# Definitions file: data/processed/molecular/methanotroph_definitions.csv
#   Curated from Knief (2015) with SILVA 138 taxonomy mapping.
#   Includes Known/Putative/Conditional flags per taxon.
#
# Usage:
#   source("code/00_harmonization/load_methanotroph_definitions.R")
#   mt_defs <- load_methanotroph_defs()
#   tax_df$mt_status <- classify_methanotrophs(tax_df, mt_defs)
#   tax_df$mt_family <- assign_display_family(tax_df, mt_defs)
# ==============================================================================

# ------------------------------------------------------------------------------
# load_methanotroph_defs: Read the canonical CSV
# ------------------------------------------------------------------------------
load_methanotroph_defs <- function(
    path = "data/processed/molecular/methanotroph_definitions.csv") {
  mt_defs <- read.csv(path, stringsAsFactors = FALSE)
  mt_defs$Include_known    <- toupper(trimws(mt_defs$Include_known))
  mt_defs$Include_putative <- toupper(trimws(mt_defs$Include_putative))
  mt_defs$Taxon_rank       <- trimws(mt_defs$Taxon_rank)
  mt_defs$Taxon            <- trimws(mt_defs$Taxon)
  return(mt_defs)
}

# ------------------------------------------------------------------------------
# Helper: check if a value is "unresolved" (NA, empty, "none", etc.)
# ------------------------------------------------------------------------------
.is_unresolved <- function(x) {
  is.na(x) | x == "" | tolower(x) == "none" | tolower(x) == "unclassified"
}

# ------------------------------------------------------------------------------
# classify_methanotrophs: Return "Known", "Putative", or NA for each row
# ------------------------------------------------------------------------------
#' Classification is hierarchical (first match wins, no double-counting):
#'   1. Genus-level Known  -> "Known"
#'   2. Family-level Known -> "Known"
#'   3. Genus-level Putative (Include_putative == "YES") -> "Putative"
#'   4. Family-level Putative -> "Putative"  ** only if genus is unresolved **
#'      (Mixed families like Beijerinckiaceae/Methylocystaceae contain many
#'       non-methanotrophic genera. An ASV in these families is only "Putative"
#'       if its genus is unresolved — a resolved genus that isn't in the Known
#'       genus whitelist is NOT a methanotroph.)
#'   5. Conditional (genus/family/phylum with CONDITIONAL) -> "Putative"
#'      (only if include_conditional = TRUE)
#'   6. Everything else -> NA
#'
#' @param tax_df  Data frame with columns "Family", "Genus", and optionally "Phylum".
#' @param mt_defs Data frame from load_methanotroph_defs().
#' @param include_conditional Logical. If TRUE, CONDITIONAL rows are included
#'        as "Putative" (e.g., Methylomirabilaceae, Methylomirabilota).
#' @return Character vector: "Known", "Putative", or NA.
classify_methanotrophs <- function(tax_df, mt_defs,
                                   include_conditional = FALSE) {
  n <- nrow(tax_df)
  status <- rep(NA_character_, n)

  genus_defs <- mt_defs[mt_defs$Taxon_rank == "Genus", ]
  fam_defs   <- mt_defs[mt_defs$Taxon_rank == "Family", ]
  phy_defs   <- mt_defs[mt_defs$Taxon_rank == "Phylum", ]

  has_genus  <- !.is_unresolved(tax_df$Genus)
  has_phylum <- "Phylum" %in% colnames(tax_df)

  # --- Known genera ---
  known_genera <- genus_defs$Taxon[genus_defs$Include_known == "YES"]
  if (length(known_genera) > 0) {
    idx <- has_genus & tax_df$Genus %in% known_genera & is.na(status)
    status[idx] <- "Known"
  }

  # --- Known families (exclusive methanotroph families) ---
  known_fams <- fam_defs$Taxon[fam_defs$Include_known == "YES"]
  if (length(known_fams) > 0) {
    idx <- tax_df$Family %in% known_fams & is.na(status)
    status[idx] <- "Known"
  }

  # --- Putative genera ---
  put_genera <- genus_defs$Taxon[genus_defs$Include_known != "YES" &
                                   genus_defs$Include_putative == "YES"]
  if (length(put_genera) > 0) {
    idx <- has_genus & tax_df$Genus %in% put_genera & is.na(status)
    status[idx] <- "Putative"
  }

  # --- Putative families (mixed families: Beijerinckiaceae, Methylocystaceae, etc.) ---
  # CRITICAL: Only mark as Putative when genus is UNRESOLVED. A resolved genus
  # that isn't in the Known genus whitelist (e.g., Roseiarcus, 1174-901-12,
  # Methylobacterium-Methylorubrum) is a non-methanotroph, not putative.
  put_fams <- fam_defs$Taxon[fam_defs$Include_known != "YES" &
                               fam_defs$Include_putative == "YES"]
  if (length(put_fams) > 0) {
    idx <- tax_df$Family %in% put_fams & !has_genus & is.na(status)
    status[idx] <- "Putative"
  }

  # --- Conditional (treated as Putative if requested) ---
  if (include_conditional) {
    # Conditional genera
    cond_genera <- genus_defs$Taxon[genus_defs$Include_putative == "CONDITIONAL"]
    if (length(cond_genera) > 0) {
      idx <- has_genus & tax_df$Genus %in% cond_genera & is.na(status)
      status[idx] <- "Putative"
    }
    # Conditional families (same genus-unresolved rule as putative families)
    cond_fams <- fam_defs$Taxon[fam_defs$Include_putative == "CONDITIONAL"]
    if (length(cond_fams) > 0) {
      idx <- tax_df$Family %in% cond_fams & !has_genus & is.na(status)
      status[idx] <- "Putative"
    }
    # Conditional phyla
    if (has_phylum && nrow(phy_defs) > 0) {
      cond_phyla <- phy_defs$Taxon[phy_defs$Include_putative == "CONDITIONAL"]
      if (length(cond_phyla) > 0) {
        idx <- tax_df$Phylum %in% cond_phyla & is.na(status)
        status[idx] <- "Putative"
      }
    }
  }

  return(status)
}

# ------------------------------------------------------------------------------
# identify_methanotrophs: Convenience wrapper — returns logical vector
# ------------------------------------------------------------------------------
#' @param include_putative If TRUE, both Known and Putative are TRUE.
#'        If FALSE, only Known are TRUE.
identify_methanotrophs <- function(tax_df, mt_defs,
                                   include_putative = FALSE,
                                   include_conditional = FALSE) {
  status <- classify_methanotrophs(tax_df, mt_defs,
                                   include_conditional = include_conditional)
  if (include_putative) {
    return(!is.na(status))
  } else {
    return(!is.na(status) & status == "Known")
  }
}

# ------------------------------------------------------------------------------
# assign_display_family: Return display family for each methanotroph row
# ------------------------------------------------------------------------------
#' Returns the ASV's Family as the display name for matched rows.
#' Non-methanotroph rows get NA.
assign_display_family <- function(tax_df, mt_defs,
                                  include_putative = FALSE,
                                  include_conditional = FALSE) {
  is_mt <- identify_methanotrophs(tax_df, mt_defs,
                                  include_putative = include_putative,
                                  include_conditional = include_conditional)
  display <- rep(NA_character_, nrow(tax_df))
  display[is_mt] <- tax_df$Family[is_mt]
  return(display)
}
