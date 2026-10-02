#!/usr/bin/env Rscript
# ==============================================================================
# 22_chamber-excess.R -- does apparent stem uptake track a chamber opening above ambient?
# ------------------------------------------------------------------------------
# SI Methods S2 argues that apparent stem uptake is measurement resolution, not a sink.
# One line of that argument: in the 2023 cross-species survey, a chamber whose starting
# concentration sits above the air around it decays toward ambient, which a linear fit
# reports as uptake. If so, detected uptakes should open further above local ambient
# than other measurements, and detected emissions should not.
#
# Local ambient = 10th percentile of the analyser trace in the 180 s either side of the
# enclosure (the lowest air the analyser saw while not enclosed; a low percentile, not a
# mean, because back-to-back enclosures leave elevated tails). Excess = the median of the
# first 10 enclosed readings minus local ambient. Diagnostic only: nothing is removed.
#
# Rebuilt 2026-10-02 from archive/revision/rev_mdf_09_local_ambient.R (July 2026), which
# used the pre-28 cm3 flux table and its own detection limit; fluxes and detection classes
# now come from outputs/data/flux_FINAL.csv (02_flux/03_precision_and_detection.R).
#
# Reads:  data/processed/flux/ymf2023_goflux_auxfile.csv, data/raw/lgr/static_2023/,
#         outputs/data/flux_FINAL.csv
# Writes: outputs/data/chamber_excess.csv (per measurement),
#         outputs/data/chamber_excess_summary.csv (the statistics SI Methods S2 quotes)
# ==============================================================================
suppressMessages({ library(dplyr) })
source("code/lib/outputs.R")
BUF <- 180; AMB_Q <- 0.10; PRE_EXCURSION_S <- 900

a <- read.csv("data/processed/flux/ymf2023_goflux_auxfile.csv", stringsAsFactors = FALSE, check.names = TRUE)
a$st <- as.POSIXct(a$start.time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
a$ln <- suppressWarnings(as.numeric(a$obs.length))
a <- a[is.finite(a$ln) & !is.na(a$st), ]
a$day <- format(a$st, "%Y-%m-%d")

# One day's analyser traces. Most days hold each file twice -- extracted (inside a
# directory named *.txt) and zipped -- so search recursively and de-duplicate.
read_day <- function(day) {
  dir <- file.path("data/raw/lgr/static_2023", day)
  if (!dir.exists(dir)) return(NULL)
  fs <- list.files(dir, pattern = "_f\\d+\\.txt$", full.names = TRUE, recursive = TRUE)
  fi <- file.info(fs); fs <- fs[!is.na(fi$isdir) & !fi$isdir & fi$size > 50000]
  if (length(fs)) fs <- fs[!duplicated(basename(fs))]
  zs <- list.files(dir, pattern = "_f\\d+\\.txt\\.zip$", full.names = TRUE)
  zs <- zs[!(sub("\\.zip$", "", basename(zs)) %in% basename(fs))]
  srcs <- c(as.list(fs), lapply(zs, function(z) list(zip = z)))
  out <- lapply(srcs, function(f) {
    con <- if (is.list(f)) unz(f$zip, utils::unzip(f$zip, list = TRUE)$Name[1]) else f
    d <- tryCatch(read.csv(con, skip = 1, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    names(d) <- trimws(names(d))
    tc <- grep("^Time$", names(d), value = TRUE)[1]; ch <- grep("^\\[CH4\\]d_ppm$", names(d), value = TRUE)[1]
    if (is.na(tc) || is.na(ch)) return(NULL)
    data.frame(tm = as.POSIXct(trimws(d[[tc]]), format = "%m/%d/%Y %H:%M:%S", tz = "UTC"),
               CH4 = suppressWarnings(as.numeric(d[[ch]])) * 1000)
  })
  out <- bind_rows(out[!vapply(out, is.null, logical(1))])
  if (!nrow(out)) return(NULL)
  out[is.finite(out$CH4) & !is.na(out$tm), ] %>% arrange(tm)
}

res <- list()
for (dy in sort(unique(a$day))) {
  tr <- read_day(dy); if (is.null(tr)) next
  ad <- a[a$day == dy, ] %>% arrange(st)
  ad$prev_end <- c(as.POSIXct(NA), head(ad$st + ad$ln, -1))
  res[[dy]] <- bind_rows(lapply(seq_len(nrow(ad)), function(i) {
    s <- ad$st[i]; e <- s + ad$ln[i]
    win <- tr[tr$tm >= s - BUF & tr$tm <= e + BUF, ]; inm <- win$tm >= s & win$tm <= e
    if (sum(inm) < 10 || sum(!inm) < 20) return(NULL)
    amb <- as.numeric(quantile(win$CH4[!inm], AMB_Q, na.rm = TRUE))
    pre <- tr$CH4[tr$tm >= s - PRE_EXCURSION_S & tr$tm < s]
    data.frame(UniqueID = ad$UniqueID[i], day = dy,
               C0_obs = median(head(win$CH4[inm], 10), na.rm = TRUE), local_ambient = amb,
               gap_prev_s = as.numeric(difftime(s, ad$prev_end[i], units = "secs")),
               max_pre15min = if (length(pre)) max(pre, na.rm = TRUE) else NA_real_)
  }))
}
L <- bind_rows(res); L$excess <- L$C0_obs - L$local_ambient

F <- read.csv(out_path("flux_FINAL.csv"), stringsAsFactors = FALSE)
F <- F[F$camp == "Cross-species" & F$type == "stem", c("UniqueID", "best.flux", "MDF", "class")]
L <- inner_join(L, F, by = "UniqueID") %>% filter(is.finite(best.flux))
write.csv(L, out_path("chamber_excess.csv"), row.names = FALSE)
cat(sprintf("2023 cross-species: %d of %d stem measurements matched to raw traces\n\n", nrow(L), nrow(F)))

med <- function(cl) median(L$excess[L$class == cl]); n <- function(cl) sum(L$class == cl)
cat(sprintf("  %-16s %4s %12s %18s\n", "class", "n", "med excess", "IQR (ppb)"))
for (cl in c("uptake", "below detection", "emission")) {
  x <- L$excess[L$class == cl]
  cat(sprintf("  %-16s %4d %12.1f %8.1f to %6.1f\n", cl, length(x), median(x), quantile(x, .25), quantile(x, .75)))
}
w_up   <- wilcox.test(L$excess[L$class == "uptake"], L$excess[L$class != "uptake"])
w_emis <- wilcox.test(L$excess[L$class == "emission"], L$excess[L$class == "below detection"])
rho    <- cor.test(L$excess, L$best.flux, method = "spearman", exact = FALSE)
rho_g  <- cor.test(L$excess, L$gap_prev_s, method = "spearman", exact = FALSE)
cat(sprintf("\n  detected uptake vs all others, Wilcoxon p = %.2g\n", w_up$p.value))
cat(sprintf("  detected emission vs below detection, Wilcoxon p = %.2g\n", w_emis$p.value))
cat(sprintf("  Spearman rho(excess, flux) = %.3f, p = %.2g\n", rho$estimate, rho$p.value))
cat(sprintf("  noise check: rho(excess, gap to previous enclosure) = %.3f, p = %.2g\n", rho_g$estimate, rho_g$p.value))

m <- L[which.min(L$best.flux), ]
cat(sprintf("\n  most negative: %s, %.3f nmol m-2 s-1; excess %.0f ppb; max CH4 in the 15 min before %.0f ppb (%.2fx local ambient)\n",
            m$UniqueID, m$best.flux, m$excess, m$max_pre15min, m$max_pre15min / m$local_ambient))

S <- data.frame(quantity = c("n_matched", "n_cross_species", "n_uptake", "n_below", "n_emission",
                             "median_excess_uptake_ppb", "median_excess_below_ppb", "median_excess_emission_ppb",
                             "wilcox_p_uptake_vs_rest", "wilcox_p_emission_vs_below", "spearman_rho_excess_flux",
                             "spearman_p_excess_flux", "spearman_rho_excess_gap", "spearman_p_excess_gap",
                             "most_negative_flux", "most_negative_excess_ppb", "most_negative_pre15min_max_ppb",
                             "most_negative_pre15min_ratio"),
                value = c(nrow(L), nrow(F), n("uptake"), n("below detection"), n("emission"),
                          med("uptake"), med("below detection"), med("emission"),
                          w_up$p.value, w_emis$p.value, rho$estimate, rho$p.value, rho_g$estimate, rho_g$p.value,
                          m$best.flux, m$excess, m$max_pre15min, m$max_pre15min / m$local_ambient))
write.csv(S, out_path("chamber_excess_summary.csv"), row.names = FALSE)
cat("\n  Written: outputs/data/chamber_excess.csv, chamber_excess_summary.csv\n")
