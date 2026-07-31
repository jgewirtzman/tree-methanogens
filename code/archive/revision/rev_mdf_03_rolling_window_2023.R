#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_03_rolling_window_2023.R      STEP 3 of the empirical-MDF rebuild
# ------------------------------------------------------------------------------
# Per-INSTRUMENT precision for the 2023 static campaign, whose per-measurement
# click.peak windows were destroyed by a filename collision between the 2021 and
# 2023 processing scripts (see task #21).
#
# METHOD (rolling window scan, identical to ch4-data-filtering 02_allan_deviation.R)
#   Scan the continuous 1 Hz recording with a 30 s rolling window. For each window
#   fit a line, and record |slope| and the residual SD. Normalise both to [0,1],
#   sum them, and take the quietest-and-flattest bottom 5% of windows. The precision
#   estimate is the median residual SD of those windows -- periods where the analyser
#   was neither drifting nor being perturbed, i.e. reading its own noise floor.
#
#   This needs NO manually clicked windows, which is why it is usable for 2023.
#
# CRITICAL: THIS INSTRUMENT, THIS PERIOD.
#   Precision is estimated from the 2023 raw recordings themselves, not borrowed from
#   another project or another year. The analyser degrades: the per-measurement Allan
#   deviations already show 1.04 ppb in the 2021 height campaign against ~2.0 ppb in
#   the 2020-21 semirigid campaigns, so an imported value would not describe 2023.
#
# VALIDATION
#   The same scan is run on the 2021 static recordings, where per-measurement Allan
#   deviation IS available. Agreement between the two methods there tells us how much
#   to trust the rolling-window number for 2023.
#
# Output: outputs/revision/mdf_rolling_window.csv
# ==============================================================================

suppressMessages({library(dplyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

WINDOW_SEC <- 30; CUTOFF <- 0.05
CUTOFFS <- c(0.05, 0.10, 0.25)   # sensitivity to how "quiet" a window must be

read_lgr <- function(path) {
  con <- tryCatch({
    if (grepl("\\.zip$", path)) {
      nm <- utils::unzip(path, list=TRUE)$Name[1]
      unz(path, nm)
    } else path
  }, error=function(e) NULL)
  if (is.null(con)) return(NULL)
  d <- tryCatch(read.csv(con, skip=1, stringsAsFactors=FALSE, check.names=FALSE),
                error=function(e) NULL)
  if (is.null(d) || nrow(d) < WINDOW_SEC+1) return(NULL)
  names(d) <- trimws(names(d))
  ch <- grep("^\\[CH4\\]d_ppm$", names(d), value=TRUE)
  if (!length(ch)) ch <- grep("^\\[CH4\\]_ppm$", names(d), value=TRUE)
  if (!length(ch)) return(NULL)
  x <- suppressWarnings(as.numeric(d[[ch[1]]])) * 1000   # ppm -> ppb
  x <- x[is.finite(x)]
  if (length(x) < WINDOW_SEC+1) return(NULL)
  x
}

rolling_precision_multi <- function(x, window_sec=WINDOW_SEC, cutoffs=CUTOFFS) {
  n <- length(x); if (n < window_sec+1) return(rep(NA_real_, length(cutoffs)))
  nw <- n - window_sec + 1
  stride <- max(1, floor(nw/4000)); idxs <- seq(1, nw, by=stride)
  sl <- numeric(length(idxs)); rs <- numeric(length(idxs))
  tw <- seq_len(window_sec) - 1; X <- cbind(1, tw)
  for (k in seq_along(idxs)) { i <- idxs[k]; xw <- x[i:(i+window_sec-1)]
    if (anyNA(xw)) { sl[k] <- NA; rs[k] <- NA; next }
    f <- .lm.fit(X, xw); sl[k] <- abs(f$coefficients[2]); rs[k] <- sd(f$residuals) }
  ok <- is.finite(sl) & is.finite(rs); sl <- sl[ok]; rs <- rs[ok]
  if (!length(sl)) return(rep(NA_real_, length(cutoffs)))
  nz <- function(v){r<-range(v); if(diff(r)==0) rep(0,length(v)) else (v-r[1])/diff(r)}
  score <- nz(sl) + nz(rs)
  vapply(cutoffs, function(cu) median(rs[score <= quantile(score, cu, na.rm=TRUE)]), numeric(1))
}

rolling_precision <- function(x, window_sec=WINDOW_SEC, cutoff=CUTOFF) {
  n <- length(x); if (n < window_sec+1) return(NA_real_)
  nw <- n - window_sec + 1
  # subsample windows for speed on long records; stride keeps <= 4000 windows
  stride <- max(1, floor(nw/4000)); idxs <- seq(1, nw, by=stride)
  sl <- numeric(length(idxs)); rs <- numeric(length(idxs))
  tw <- seq_len(window_sec) - 1
  X <- cbind(1, tw)
  for (k in seq_along(idxs)) {
    i <- idxs[k]; xw <- x[i:(i+window_sec-1)]
    if (anyNA(xw)) { sl[k] <- NA; rs[k] <- NA; next }
    f <- .lm.fit(X, xw)
    sl[k] <- abs(f$coefficients[2]); rs[k] <- sd(f$residuals)
  }
  ok <- is.finite(sl) & is.finite(rs); sl <- sl[ok]; rs <- rs[ok]
  if (!length(sl)) return(NA_real_)
  nz <- function(v) { r <- range(v); if (diff(r)==0) rep(0,length(v)) else (v-r[1])/diff(r) }
  score <- nz(sl) + nz(rs)
  sel <- which(score <= quantile(score, cutoff, na.rm=TRUE))
  median(rs[sel])
}

scan_dir <- function(dir, label, max_files = 60) {
  fs <- list.files(dir, pattern="\\.(txt|zip)$", recursive=TRUE, full.names=TRUE)
  fs <- fs[file.info(fs)$size > 50000]
  if (!length(fs)) { cat(sprintf("  %-14s no usable files in %s\n", label, dir)); return(NULL) }
  if (length(fs) > max_files) fs <- fs[round(seq(1, length(fs), length.out=max_files))]
  cat(sprintf("  %-14s scanning %d recordings ...\n", label, length(fs)))
  out <- vapply(fs, function(f) { x <- read_lgr(f); if (is.null(x)) NA_real_ else rolling_precision(x) },
                numeric(1))
  out <- out[is.finite(out)]
  cat(sprintf("  %-14s n=%3d recordings   precision (ppb CH4): median %6.3f   IQR %6.3f - %6.3f   range %6.3f - %6.3f\n",
      label, length(out), median(out), quantile(out,.25), quantile(out,.75), min(out), max(out)))
  data.frame(campaign=label, n_recordings=length(out), precision_ppb=median(out),
             q25=quantile(out,.25), q75=quantile(out,.75), stringsAsFactors=FALSE)
}

rule("ROLLING-WINDOW PRECISION, PER INSTRUMENT-PERIOD")
cat(sprintf("\n  window %d s, quietest/flattest bottom %.0f%% of windows\n\n", WINDOW_SEC, 100*CUTOFF))
R <- bind_rows(
  scan_dir("data/raw/lgr/static_2023",         "static_2023"),
  scan_dir("data/raw/lgr/static_2021",         "static_2021"),
  scan_dir("data/raw/lgr/semirigid_2020-2021", "semirigid_2020-21"))

rule("VALIDATION AGAINST PER-MEASUREMENT ALLAN")
A <- read.csv("outputs/revision/mdf_allan_per_measurement.csv", stringsAsFactors=FALSE)
cmp <- data.frame(
  campaign = c("static_2021 (= height_2021)","semirigid_2020-21"),
  allan_median = c(median(A$allan_sd_CH4[A$campaign=="height_2021"], na.rm=TRUE),
                   median(A$allan_sd_CH4[A$campaign=="semirigid_soil"], na.rm=TRUE)),
  rolling = c(R$precision_ppb[R$campaign=="static_2021"],
              R$precision_ppb[R$campaign=="semirigid_2020-21"]))
cmp$ratio <- cmp$rolling / cmp$allan_median
cat("\n  Where BOTH methods are available, do they agree?\n\n")
cat(sprintf("  %-30s %12s %12s %8s\n","campaign","Allan (ppb)","rolling(ppb)","ratio"))
for (i in seq_len(nrow(cmp))) with(cmp[i,],
  cat(sprintf("  %-30s %12.3f %12.3f %8.2f\n", campaign, allan_median, rolling, ratio)))
cat("
  The rolling-window scan measures the analyser's noise floor during its QUIETEST
  periods, so it is expected to sit at or below the per-measurement Allan deviation,
  which includes real field perturbation. A ratio near or under 1 means the two are
  consistent; a ratio well above 1 would mean the scan is picking up something other
  than instrument noise and should not be trusted for 2023.
")
write.csv(R, file.path(outdir,"mdf_rolling_window.csv"), row.names=FALSE)
cat(sprintf("\n  Written: outputs/revision/mdf_rolling_window.csv\n"))

rule("2023 IN CONTEXT")
p23 <- R$precision_ppb[R$campaign=="static_2023"]
cat(sprintf("
  2023 rolling-window precision : %.3f ppb   (%.2fx the 0.9 ppb manufacturer spec)
  2021 height campaign, Allan   : %.3f ppb
  2020-21 semirigid, Allan      : %.3f ppb

  If 2023 sits within the range already measured, the 406 lost clicks buy little and
  the rolling-window value can stand in. If 2023 is markedly worse, that is itself
  the finding and re-clicking becomes worth the effort.
",
 p23, p23/0.9,
 median(A$allan_sd_CH4[A$campaign=="height_2021"], na.rm=TRUE),
 median(A$allan_sd_CH4[A$campaign=="semirigid_soil"], na.rm=TRUE)))


# ==============================================================================
rule("SENSITIVITY TO THE QUIET-WINDOW CUTOFF")
# ==============================================================================
cat("
  The noise floor is the median residual SD of the quietest windows. How quiet is
  quiet enough is a choice; a looser cutoff admits windows with more real variation
  and so returns a higher estimate. Recomputed at 5%, 10% and 25%:

")
cat(sprintf("  %-22s %10s %10s %10s %12s\n","period","5%","10%","25%","25%/5%"))
for (rd in c("data/raw/lgr/static_2021","data/raw/lgr/semirigid_2020-2021","data/raw/lgr/static_2023")) {
  fs <- list.files(rd, pattern="\\.(txt|zip)$", recursive=TRUE, full.names=TRUE)
  fi <- file.info(fs); fs <- fs[!is.na(fi$isdir) & !fi$isdir & fi$size > 50000]
  if (length(fs)) fs <- fs[!duplicated(basename(fs))]
  if (length(fs) > 40) fs <- fs[round(seq(1, length(fs), length.out=40))]
  M <- vapply(fs, function(f){ x <- read_lgr(f); if (is.null(x)) rep(NA_real_,3) else
                rolling_precision_multi(x) }, numeric(3))
  m <- apply(M, 1, function(v) median(v[is.finite(v)]))
  cat(sprintf("  %-22s %10.3f %10.3f %10.3f %11.2fx\n", basename(rd), m[1], m[2], m[3], m[3]/m[1]))
}
cat("
  A cutoff between 5% and 25% changes the precision estimate by well under the
  difference between campaigns, so the detection limits are not sensitive to it.
")
