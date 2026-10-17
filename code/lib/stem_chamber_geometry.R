# stem_chamber_geometry.R -- semi-rigid stem chamber area and volume (monthly survey
# 2020-21). One formula and one rule, used by 02_flux/semirigid/01_calc_chamber_dims.R
# (the formula), 02_flux/semirigid/02_join_flux_geometry.R (tagged stems) and
# 02_flux/rescue/01_untagged_auxfile.R (untagged stems).
#
# FORMULA (Siegenthaler et al. 2016, with the wraparound correction): sector fraction
#   K = min(1, L / (pi * (D + 2T))); Vc = K * pi/4 * ((D + 2T)^2 - D^2) * H - Vwedges;
#   Sc = K * pi * D * H - Swedges (both floored at 0). D stem diameter, T sleeve
#   thickness, L chamber length, H chamber height (cm).
#
# RULE for one deployment (geometry_source says which):
#   sheet              the chamber sheet of that month (data/raw/flux_chamber_corrections/
#                      Chamber vol corr/Tree_fluxes_<month>.xlsx) has a row for that stem
#                      and the chamber used that day: that row's geometry (the formula on
#                      the row's own D, T, L, H and wedges).
#   sheet (chamber not recorded on the flux sheet)
#                      the flux sheet names no numbered chamber ("na", "Baby chamber") and
#                      the chamber sheet of that month has a row for the stem: that row.
#   computed           otherwise: the formula on the stem's diameter (the median of its
#                      chamber-sheet D values; stems on no chamber sheet: the monthly
#                      tree list D_stem) and the dimensions of the chamber used:
#                        chambers 1-6: the dimensions the chamber sheets give it (its
#                          most frequent T, L, H, wedge row);
#                        chambers 7, 8, 9 (used from March 2021, after the last chamber
#                          sheet): their dimensions were not found in the raw sheets, the
#                          field files or the notes; they take the mean of the dimensions
#                          of chambers 1-6 (Jon, 2026-10-16: a minor difference);
#                        chamber not recorded and no sheet row: the mean of the
#                          dimensions of the numbered chambers used that month.
# Months with a chamber sheet: June, July, August, September, October, December 2020.
SEMIRIGID_SHEET_DIR <- "data/raw/flux_chamber_corrections/Chamber vol corr"
SEMIRIGID_SHEET_MONTHS <- c(june = "2020-06", july = "2020-07", august = "2020-08", september = "2020-09",
                            october = "2020-10", december = "2020-12")
SEMIRIGID_MEAN_CHAMBERS <- c("7", "8", "9")      # dimensions not found: mean of chambers 1-6

semirigid_geometry <- function(D, T, L, H, Swedges, Vwedges) {
  Dext <- D + 2 * T
  K <- pmin(L / (pi * Dext), 1)
  Swedges[!is.finite(Swedges)] <- 0; Vwedges[!is.finite(Vwedges)] <- 0
  Vc <- K * (pi / 4) * (Dext^2 - D^2) * H - Vwedges
  Sc <- K * pi * D * H - Swedges
  data.frame(Sc_cm2 = ifelse(is.finite(Sc), pmax(Sc, 0), NA_real_), Vc_cm3 = ifelse(is.finite(Vc), pmax(Vc, 0), NA_real_),
             K = K, wraparound = is.finite(L / (pi * Dext)) & L / (pi * Dext) > 1)
}

# every row of the chamber sheets: stem, month, chamber, dimensions and its geometry
semirigid_sheet_rows <- function(root = ".") {
  num <- function(x) suppressWarnings(readr::parse_number(as.character(x)))
  out <- do.call(rbind, lapply(names(SEMIRIGID_SHEET_MONTHS), function(m) {
    f <- file.path(root, SEMIRIGID_SHEET_DIR, paste0("Tree_fluxes_", m, ".xlsx"))
    x <- as.data.frame(readxl::read_excel(f, sheet = 1, col_types = "text"))
    x <- x[!apply(is.na(x) | x == "", 1, all), , drop = FALSE]
    col <- function(p) x[[grep(p, names(x), ignore.case = TRUE)[1]]]
    data.frame(month = unname(SEMIRIGID_SHEET_MONTHS[m]), sheet_row = seq_len(nrow(x)),
               sample = trimws(sub("\\.0$", "", as.character(x$Sample))), chamber = trimws(as.character(x$`Chamber used`)),
               D = num(col("^D\\s*stem")), T = num(col("^T\\s*\\(")), L = num(col("^Chamber\\s*length")),
               H = num(col("^Chamber\\s*height")), Swedges = num(col("^Swedges")), Vwedges = num(col("^Two\\s*V\\s*wedges")),
               stringsAsFactors = FALSE)
  }))
  g <- semirigid_geometry(out$D, out$T, out$L, out$H, out$Swedges, out$Vwedges)
  cbind(out, g[, c("Sc_cm2", "Vc_cm3")])
}

# dimensions of each chamber (1-6: its most frequent row on the sheets; 7-9: the mean of 1-6)
semirigid_chamber_dims <- function(rows) {
  r <- rows[rows$chamber %in% as.character(1:6) & is.finite(rows$T) & is.finite(rows$L) & is.finite(rows$H), ]
  k <- paste(r$chamber, r$T, r$L, r$H, r$Swedges, r$Vwedges)
  d <- do.call(rbind, lapply(split(seq_len(nrow(r)), r$chamber), function(i) {
    tab <- sort(table(k[i]), decreasing = TRUE); j <- i[match(names(tab)[1], k[i])]
    data.frame(chamber = r$chamber[j], T = r$T[j], L = r$L[j], H = r$H[j], Swedges = r$Swedges[j], Vwedges = r$Vwedges[j],
               n_rows = length(i), n_modal = as.integer(tab[1]), stringsAsFactors = FALSE)
  }))
  mn <- data.frame(T = mean(d$T), L = mean(d$L), H = mean(d$H), Swedges = mean(d$Swedges), Vwedges = mean(d$Vwedges))
  rbind(d, do.call(rbind, lapply(SEMIRIGID_MEAN_CHAMBERS, function(ch)
    data.frame(chamber = ch, mn, n_rows = 0L, n_modal = 0L))))
}

# geometry for deployments: sample (stem label as on the chamber sheets), month (YYYY-MM),
# chamber (as on the flux sheet), D_fallback (stem diameter where the stem is on no sheet),
# month_chambers (list: month -> chambers named on the flux sheet that month)
semirigid_deployment_geometry <- function(sample, month, chamber, D_fallback = NA_real_, rows, dims, month_chambers) {
  n <- length(sample); D_fallback <- rep_len(D_fallback, n)
  chamber <- sub("\\.0$", "", trimws(as.character(chamber)))
  numbered <- chamber %in% dims$chamber
  Dmed <- tapply(rows$D, rows$sample, median, na.rm = TRUE)
  out <- data.frame(Sc_cm2 = rep(NA_real_, n), Vc_cm3 = NA_real_, D_stem_cm = NA_real_, chamber_dims = NA_character_,
                    geometry_source = NA_character_, stringsAsFactors = FALSE)
  for (i in seq_len(n)) {
    s <- rows[rows$sample == sample[i] & rows$month == month[i], ]
    hit <- if (numbered[i]) s[s$chamber == chamber[i], ] else s
    if (nrow(hit) >= 1) {
      hit <- hit[1, ]
      out[i, c("Sc_cm2", "Vc_cm3", "D_stem_cm")] <- hit[, c("Sc_cm2", "Vc_cm3", "D")]
      out$chamber_dims[i] <- hit$chamber
      out$geometry_source[i] <- if (numbered[i]) "sheet" else "sheet (chamber not recorded on the flux sheet)"
      next
    }
    D <- if (sample[i] %in% names(Dmed)) unname(Dmed[sample[i]]) else D_fallback[i]
    if (numbered[i]) { dm <- dims[dims$chamber == chamber[i], ]; lab <- chamber[i] }
    else {
      used <- unique(month_chambers[[month[i]]]); used <- used[used %in% dims$chamber]
      dm <- dims[dims$chamber %in% used, ]
      dm <- data.frame(T = mean(dm$T), L = mean(dm$L), H = mean(dm$H), Swedges = mean(dm$Swedges), Vwedges = mean(dm$Vwedges))
      lab <- paste0("mean of ", paste(sort(used), collapse = "/"))
    }
    g <- semirigid_geometry(D, dm$T, dm$L, dm$H, dm$Swedges, dm$Vwedges)
    out$Sc_cm2[i] <- g$Sc_cm2; out$Vc_cm3[i] <- g$Vc_cm3; out$D_stem_cm[i] <- D; out$chamber_dims[i] <- lab
    out$geometry_source[i] <- paste0("computed (", if (chamber[i] %in% SEMIRIGID_MEAN_CHAMBERS) "chamber dimensions not found: mean of chambers 1-6"
                                     else if (numbered[i]) "chamber sheet dimensions" else "chamber not recorded: mean of the chambers used that month", ")")
  }
  out
}
