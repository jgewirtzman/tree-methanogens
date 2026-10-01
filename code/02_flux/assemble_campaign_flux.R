# ==============================================================================
# assemble_campaign_flux.R -- build one campaign's tree-flux table from its
# goFlux results. Non-interactive and re-runnable.
# ------------------------------------------------------------------------------
# WHY THIS EXISTS (2026-09-30). The three tree goFlux scripts (semirigid monthly,
# 2021 multi-height, 2023 cross-species) wrote their intermediates AND their final
# table to shared filenames. Whichever ran last won:
#   - methanogen_tree_flux_complete_dataset.csv held 2021 data until the 2023 run
#     overwrote it (Sep 2025). 04_harmonize_all_data.R expected 2021 and broke;
#     fourteen other readers expected 2023 and happened to be right.
#   - CO2_flux_lgr_results.csv held 2023 while its CH4/best siblings held 2021.
# Every per-campaign file now carries the campaign in its name, and the final
# assembly step lives here, so it can be rerun without repeating the interactive
# window selection (click.peak2) in the goFlux scripts.
#
# Usage (from the repo root, or sourced at the end of a goFlux script):
#   Rscript code/02_flux/assemble_campaign_flux.R 2021_multiheight
#   source("code/02_flux/assemble_campaign_flux.R"); assemble_campaign("2023_cross_species")
#
# Reads:  data/processed/flux/<aux file>
#         data/processed/flux/{CO2,CH4}_best_flux_lgr_results_<campaign>.csv
# Writes: data/processed/flux/<campaign output> (see CAMPAIGNS)
# ==============================================================================
suppressMessages({ library(readr); library(dplyr) })

CAMPAIGNS <- list(
  `2021_multiheight`   = list(aux = "goflux_auxfile.csv",
                              out = "tree_flux_2021_multiheight.csv"),
  `2023_cross_species` = list(aux = "ymf2023_goflux_auxfile.csv",
                              out = "tree_flux_2023_cross_species.csv"),
  semirigid_tree       = list(aux = "auxfile_goFlux_with_weather_formatted_datetime.csv",
                              out = "semirigid_tree_final_complete_dataset.csv")
)

flux_dir <- function() {
  for (d in c("data/processed/flux", "../../data/processed/flux", "../../../data/processed/flux"))
    if (dir.exists(d)) return(d)
  stop("cannot find data/processed/flux from ", getwd(), call. = FALSE)
}

# Path of a per-campaign goFlux intermediate, e.g. results_path("CH4_best_flux", "2021_multiheight")
results_path <- function(stem, campaign, dir = flux_dir())
  file.path(dir, sprintf("%s_lgr_results_%s.csv", stem, campaign))

assemble_campaign <- function(campaign, dir = flux_dir()) {
  cfg <- CAMPAIGNS[[campaign]]
  if (is.null(cfg)) stop("unknown campaign: ", campaign, call. = FALSE)
  aux <- read_csv(file.path(dir, cfg$aux), show_col_types = FALSE)
  out <- aux
  for (gas in c("CO2", "CH4")) {
    f <- results_path(paste0(gas, "_best_flux"), campaign, dir)
    if (!file.exists(f)) stop("missing ", f, "; run the campaign's goFlux script first", call. = FALSE)
    res <- read_csv(f, show_col_types = FALSE) %>% rename_with(~ paste0(gas, "_", .x), -UniqueID)
    stopifnot(!anyDuplicated(res$UniqueID))
    out <- left_join(out, res, by = "UniqueID")
  }
  stopifnot(nrow(out) == nrow(aux))   # a left join on a unique key cannot add or drop rows
  write_csv(out, file.path(dir, cfg$out))
  cat(sprintf("assemble_campaign(%s): %d rows (%d with a CH4 flux) -> %s\n", campaign,
              nrow(out), sum(!is.na(out$CH4_best.flux)), cfg$out))
  invisible(out)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("usage: Rscript code/02_flux/assemble_campaign_flux.R <campaign>", call. = FALSE)
  for (a in args) assemble_campaign(a)
}
