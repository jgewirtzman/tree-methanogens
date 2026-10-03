# ==============================================================================
# rebuild_closure_windows.R -- recover hand-picked closure windows that were never saved
# ------------------------------------------------------------------------------
# One-time recovery tool (run 2026-10-01). The closure windows picked by hand with
# click.peak2() for the monthly-survey stems and the 2023 cross-species campaign were
# never written to disk, so those fluxes could not be re-fitted from the raw traces.
# Every stored fit keeps its slope, intercept (C0) and number of points, and the raw
# analyzer traces are archived. At one reading per second the window is the stretch of
# trace, within the deployment's observation window, whose linear fit reproduces the
# stored slope and C0 with the stored number of points.
#
# Result on 2026-10-01: every window recovered to machine precision (2023: 338/338;
# monthly: 369/369), and re-fitting reproduces every stored flux when the stored
# linear/non-linear model choice is kept (below).
#
# Writes, in the format click.peak2() produces (columns, flag, Etime and the corrected
# start/end times), keeping only the rows inside each picked window:
#   data/raw/flux_windows/lgr_manual_identification_semirigid_tree.csv
#   data/raw/flux_windows/lgr_manual_identification_2023_cross_species.csv
# and the stored model choices for every campaign:
#   data/raw/flux_windows/flux_model_choices.csv
# (written to data/processed/flux/ on 2026-10-01; moved to data/raw/ the same day)
# The goFlux fitting scripts load these instead of opening the picker.
#
# Usage (repository root):  Rscript code/tools/rebuild_closure_windows.R [monthly|2023]
# ==============================================================================
suppressPackageStartupMessages({ library(goFlux); library(dplyr); library(readr) })
FD <- "data/processed/flux"; OUT <- "data/raw/flux_windows"
CFG <- list(
  monthly = list(raw = "data/raw/lgr/semirigid_2020-2021", aux = "auxfile_goFlux_with_weather.csv", ol = 600, sh = 300,
                 stored = "semirigid_tree_final_complete_dataset.csv", sfx = ".x",
                 out = "lgr_manual_identification_semirigid_tree.csv"),
  `2023`  = list(raw = "data/raw/lgr/static_2023", aux = "ymf2023_goflux_auxfile.csv", ol = 300, sh = 1200,
                 stored = "tree_flux_2023_cross_species.csv", sfx = "",
                 out = "lgr_manual_identification_2023_cross_species.csv"))
camps <- commandArgs(trailingOnly = TRUE); if (!length(camps)) camps <- names(CFG)

import_lgr <- function(dp) {   # same import as the goFlux scripts (zips + loose .txt, UGGA)
  zips <- list.files(dp, recursive = TRUE, pattern = "\\.zip$", full.names = TRUE)
  tx <- tempfile(); dir.create(tx); for (z in zips) try(unzip(z, exdir = tx, overwrite = TRUE), silent = TRUE)
  ex <- list.files(dp, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE); ex <- ex[file.size(ex) > 0]
  cd <- tempfile(); dir.create(cd); file.copy(ex, cd)
  file.copy(list.files(tx, recursive = TRUE, pattern = "\\.txt$", full.names = TRUE), cd)
  on.exit(unlink(c(tx, cd), recursive = TRUE))
  import2RData(path = cd, instrument = "UGGA", date.format = "mdy", timezone = "UTC",
               keep_all = FALSE, prec = c(0.35, 0.9, 200), merge = TRUE)
}

for (cp in camps) {
  cfg <- CFG[[cp]]; cat(sprintf("\n== %s ==\n", cp))
  lgr <- import_lgr(cfg$raw)
  aux <- read_csv(file.path(FD, cfg$aux), show_col_types = FALSE) %>% mutate(start.time = as.POSIXct(start.time, tz = "UTC"))
  ow  <- obs.win(inputfile = lgr, auxfile = aux, gastype = "CO2dry_ppm", obs.length = cfg$ol, shoulder = cfg$sh)
  names(ow) <- vapply(ow, function(x) as.character(unique(x$UniqueID)[1]), "")
  st <- read.csv(file.path(FD, cfg$stored), check.names = FALSE)
  for (v in c("CH4_nb.obs", "CH4_LM.slope", "CH4_LM.C0")) st[[v]] <- suppressWarnings(as.numeric(st[[paste0(v, cfg$sfx)]]))
  ids <- st$UniqueID[!is.na(st$CH4_nb.obs) & !is.na(st$CH4_LM.slope) & !is.na(st$CH4_LM.C0) & st$UniqueID %in% names(ow)]
  out <- list(); score <- c()
  for (id in ids) {
    w <- ow[[id]] %>% arrange(POSIX.time); r <- st[match(id, st$UniqueID), ]
    nb <- as.integer(r$CH4_nb.obs); y <- w$CH4dry_ppb; tt <- as.numeric(w$POSIX.time)
    if (nb > length(y)) next
    best <- c(s = NA_real_, sc = Inf)
    for (s in seq_len(length(y) - nb + 1)) {
      i <- s:(s + nb - 1); if (anyNA(y[i])) next
      b <- unname(coef(lm.fit(cbind(1, tt[i] - tt[s]), y[i]))); if (anyNA(b)) next
      sc <- abs(b[2] - r$CH4_LM.slope) / max(abs(r$CH4_LM.slope), 1e-6) + abs(b[1] - r$CH4_LM.C0) / max(abs(r$CH4_LM.C0), 1)
      if (sc < best[["sc"]]) best <- c(s = s, sc = sc)
    }
    i <- best[["s"]]:(best[["s"]] + nb - 1)
    out[[id]] <- w %>% mutate(start.time_corr = w$POSIX.time[best[["s"]]], end.time_corr = w$POSIX.time[max(i)],
                              obs.length_corr = as.numeric(end.time_corr - start.time_corr, units = "secs"),
                              Etime = as.numeric(POSIX.time - start.time_corr, units = "secs"),
                              flag = as.integer(seq_len(n()) %in% i))
    score[id] <- best[["sc"]]
  }
  stopifnot(all(score < 1e-9))                   # every window must reproduce its stored fit exactly
  # only the picked window is kept (flag == 1): goFlux fits those rows alone, and the full
  # observation windows made the 2023 file 0.5 GB
  write_csv(bind_rows(out) %>% filter(flag == 1), file.path(OUT, cfg$out))
  cat(sprintf("%d of %d stored fits: windows recovered (max mismatch %.1e) -> %s\n",
              length(out), length(ids), max(score), cfg$out))
}

# ---- stored model choices (LM or HM) for every campaign ----------------------
# goFlux 0.2.0 (renv.lock) converges the Hutchinson-Mosier fit on some deployments
# where the version used in 2025 did not, and so would choose differently for ~9%
# of deployments. The choice made at the time is part of the record.
pick <- function(f, id = "UniqueID", ch4 = "CH4_model", co2 = "CO2_model", campaign) {
  d <- read.csv(file.path(FD, f), check.names = FALSE)
  bind_rows(if (ch4 %in% names(d)) data.frame(UniqueID = d[[id]], gas = "CH4", model = d[[ch4]]),
            if (co2 %in% names(d)) data.frame(UniqueID = d[[id]], gas = "CO2", model = d[[co2]])) %>%
    filter(!is.na(UniqueID), !is.na(model)) %>% mutate(campaign = campaign) }
MC <- bind_rows(
  pick("semirigid_tree_final_complete_dataset.csv", ch4 = "CH4_model.x", campaign = "monthly_2020_2021"),
  pick("semirigid_tree_final_complete_dataset_soil.csv", campaign = "monthly_2020_2021_soil"),
  pick("tree_flux_2021_multiheight.csv", campaign = "2021_multiheight"),
  pick("tree_flux_2023_cross_species.csv", campaign = "2023_cross_species")) %>% distinct()
stopifnot(!anyDuplicated(MC[, c("UniqueID", "gas")]))
write_csv(MC, file.path(OUT, "flux_model_choices.csv"))
cat(sprintf("\nflux_model_choices.csv: %d choices (%s)\n", nrow(MC),
            paste(names(table(MC$campaign)), table(MC$campaign), collapse = "; ")))
