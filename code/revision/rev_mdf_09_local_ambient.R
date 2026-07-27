#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_09_local_ambient.R      STEP 9
# ------------------------------------------------------------------------------
# IS THE CHAMBER STARTING ABOVE LOCAL AMBIENT?
#
# The campaign-level C0 screen (task #16) compares each chamber's starting
# concentration to the campaign median. That caught two grossly contaminated soil
# chambers but MISSED the clearest stem artefact we have found:
#
#   20230623_974_091900   C0 2116 ppb, campaign median 2036 -> only 1.04x, passes
#   but local ambient in the minutes either side was 1990 ppb, so the chamber opened
#   ~126 ppb HIGH, and the "flux" of -0.391 (r2 0.964) is the decay back to ambient.
#   A CH4 spike to 3796 ppb occurred minutes before.
#
# So the diagnostic quantity is elevation above LOCAL ambient, not above a campaign
# constant.
#
# WHY THIS WILL BE NOISY (anticipated, and measured below)
#   - measurements run back-to-back, so the window before one enclosure often contains
#     the tail of the previous one, which may be genuinely elevated by a real flux
#   - moving the analyser and fitting the chamber perturbs the intake
#   - the gap between consecutive measurements varies
#   A low PERCENTILE of the surrounding air is therefore used rather than a mean, and
#   both sides of the measurement are pooled: the aim is the lowest concentration the
#   analyser saw while not enclosed, which is the best available proxy for ambient.
#   The gap to the previous measurement is reported so the noise can be seen.
#
# This is diagnostic, not a filter. Nothing is removed here.
#
# Output: outputs/revision/local_ambient_check.csv / .txt
#         outputs/revision/fig_local_ambient.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

BUF <- 180      # s of trace to inspect either side
AMB_Q <- 0.10   # percentile of surrounding air taken as "ambient"

a <- read.csv("data/processed/flux/ymf2023_goflux_auxfile.csv", stringsAsFactors=FALSE, check.names=TRUE)
a$st <- as.POSIXct(a$start.time, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
a$ln <- suppressWarnings(as.numeric(a$obs.length))
a <- a[is.finite(a$ln) & !is.na(a$st), ]
a$day <- format(a$st, "%Y-%m-%d")
cat(sprintf("2023 auxfile: %d measurements across %d days\n", nrow(a), dplyr::n_distinct(a$day)))

read_day <- function(day) {
  dir <- file.path("data/raw/lgr/static_2023", day)
  if (!dir.exists(dir)) return(NULL)
  # On most days the "*.txt" entries are DIRECTORIES containing the extracted file one
  # level down, with a sibling "*.txt.zip" holding the same content. A non-recursive
  # listing therefore matches directories, which a size filter then discards -- which
  # is why an earlier version found traces for only 16% of measurements. Search
  # recursively for real files and de-duplicate on basename.
  fs <- list.files(dir, pattern="_f\\d+\\.txt$", full.names=TRUE, recursive=TRUE)
  fi <- file.info(fs)
  fs <- fs[!is.na(fi$isdir) & !fi$isdir & fi$size > 50000]
  if (length(fs)) fs <- fs[!duplicated(basename(fs))]
  # fall back to the zip archives where no extracted copy exists
  zs <- list.files(dir, pattern="_f\\d+\\.txt\\.zip$", full.names=TRUE)
  zs <- zs[!(sub("\\.zip$","",basename(zs)) %in% basename(fs))]
  if (!length(fs) && !length(zs)) return(NULL)
  srcs <- c(as.list(fs), lapply(zs, function(z) list(zip=z)))
  out <- lapply(srcs, function(f) {
    con <- if (is.list(f)) { nm <- utils::unzip(f$zip, list=TRUE)$Name[1]; unz(f$zip, nm) } else f
    d <- tryCatch(read.csv(con, skip=1, stringsAsFactors=FALSE, check.names=FALSE),
                  error=function(e) NULL)
    if (is.null(d)) return(NULL)
    names(d) <- trimws(names(d))
    tc <- grep("^Time$", names(d), value=TRUE)[1]
    ch <- grep("^\\[CH4\\]d_ppm$", names(d), value=TRUE)[1]
    if (is.na(tc) || is.na(ch)) return(NULL)
    tm <- as.POSIXct(trimws(d[[tc]]), format="%m/%d/%Y %H:%M:%S", tz="UTC")
    data.frame(tm=tm, CH4=suppressWarnings(as.numeric(d[[ch]]))*1000)
  })
  out <- bind_rows(out[!vapply(out, is.null, logical(1))])
  out[is.finite(out$CH4) & !is.na(out$tm), ]
}

res <- list()
for (dy in sort(unique(a$day))) {
  tr <- read_day(dy)
  if (is.null(tr) || !nrow(tr)) { cat(sprintf("  %s  no trace\n", dy)); next }
  ad <- a[a$day==dy, ] %>% arrange(st)
  ad$prev_end <- c(as.POSIXct(NA), head(ad$st + ad$ln, -1))
  r <- lapply(seq_len(nrow(ad)), function(i) {
    s <- ad$st[i]; e <- s + ad$ln[i]
    win <- tr[tr$tm >= s-BUF & tr$tm <= e+BUF, ]
    inm <- win$tm >= s & win$tm <= e
    if (sum(inm) < 10 || sum(!inm) < 20) return(NULL)
    amb <- as.numeric(quantile(win$CH4[!inm], AMB_Q, na.rm=TRUE))
    c0  <- median(head(win$CH4[inm], 10), na.rm=TRUE)
    data.frame(UniqueID=ad$UniqueID[i], day=dy,
               C0_obs=c0, local_ambient=amb, excess=c0-amb,
               n_in=sum(inm), n_out=sum(!inm),
               gap_prev_s=as.numeric(difftime(s, ad$prev_end[i], units="secs")))
  })
  r <- bind_rows(r[!vapply(r, is.null, logical(1))])
  res[[dy]] <- r
  cat(sprintf("  %s  %3d/%3d measurements matched to trace\n", dy, nrow(r), nrow(ad)))
}
L <- bind_rows(res)
cat(sprintf("\nMatched %d of %d 2023 measurements to raw traces\n", nrow(L), nrow(a)))

# --- join flux + detection ------------------------------------------------------
ALL <- read.csv("outputs/revision/mdf_allan_all_campaigns.csv", stringsAsFactors=FALSE)
G <- read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE)
sig <- 1.298/0.86
G$MDF95 <- sig*qnorm(.975)/G$nb.obs*G$flux.term
L <- L %>% inner_join(G %>% transmute(UniqueID, flux=best.flux, r2=LM.r2, MDF95, C0_goflux=C0),
                      by="UniqueID") %>% filter(is.finite(flux))
L$class <- ifelse(L$flux < -L$MDF95, "sig.negative",
           ifelse(L$flux >  L$MDF95, "sig.positive", "below MDF"))
write.csv(L, file.path(outdir,"local_ambient_check.csv"), row.names=FALSE)

rule("EXCESS ABOVE LOCAL AMBIENT, BY DETECTION CLASS")
cat(sprintf("\n  ambient = %.0fth percentile of the %d s either side of the enclosure\n\n",
            100*AMB_Q, BUF))
cat(sprintf("  %-14s %6s %10s %10s %10s %10s\n","class","n","med excess","IQR lo","IQR hi","% > 50 ppb"))
for (cl in c("sig.negative","below MDF","sig.positive")) {
  d <- L[L$class==cl,]
  if (!nrow(d)) next
  cat(sprintf("  %-14s %6d %10.1f %10.1f %10.1f %9.0f%%\n", cl, nrow(d),
      median(d$excess), quantile(d$excess,.25), quantile(d$excess,.75),
      100*mean(d$excess > 50)))
}
w <- suppressWarnings(wilcox.test(L$excess[L$class=="sig.negative"],
                                  L$excess[L$class!="sig.negative"]))
cat(sprintf("\n  sig.negative vs rest, Wilcoxon p = %.4g\n", w$p.value))
cat(sprintf("  Spearman rho(excess, flux) = %.3f (p = %.3g)\n",
    suppressWarnings(cor.test(L$excess, L$flux, method="spearman")$estimate),
    suppressWarnings(cor.test(L$excess, L$flux, method="spearman")$p.value)))

rule("HOW NOISY IS THE METRIC? (the concern raised)")
cat(sprintf("
  gap to previous measurement: median %.0f s, IQR %.0f - %.0f s, %.0f%% under 60 s
  excess vs gap, Spearman rho = %.3f (p = %.3g)
    -> a strong negative rho would mean the metric is mostly measuring how recently
       the previous enclosure ended, i.e. carryover from back-to-back sampling rather
       than a property of this measurement.
",
 median(L$gap_prev_s,na.rm=TRUE), quantile(L$gap_prev_s,.25,na.rm=TRUE),
 quantile(L$gap_prev_s,.75,na.rm=TRUE), 100*mean(L$gap_prev_s<60,na.rm=TRUE),
 suppressWarnings(cor.test(L$excess,L$gap_prev_s,method="spearman")$estimate),
 suppressWarnings(cor.test(L$excess,L$gap_prev_s,method="spearman")$p.value)))

rule("THE 32 SIGNIFICANT NEGATIVES, RANKED BY EXCESS")
n <- L[L$class=="sig.negative",] %>% arrange(desc(excess))
cat(sprintf("\n  %-28s %9s %9s %9s %8s %7s %8s\n","UniqueID","flux","C0 obs","ambient","excess","r2","gap s"))
for (i in seq_len(nrow(n))) with(n[i,],
  cat(sprintf("  %-28s %9.4f %9.1f %9.1f %8.1f %7.3f %8.0f\n",
      substr(UniqueID,1,28), flux, C0_obs, local_ambient, excess, r2, gap_prev_s)))

p <- ggplot(L, aes(excess, flux, colour=class)) +
  geom_hline(yintercept=0, linetype="22", colour="grey50") +
  geom_vline(xintercept=0, linetype="22", colour="grey50") +
  geom_point(size=1.2, alpha=.75) +
  scale_y_continuous(trans=scales::pseudo_log_trans(sigma=.01, base=10)) +
  scale_colour_manual(values=c(`sig.negative`="#B2182B", `below MDF`="grey65",
                               `sig.positive`="#2166AC")) +
  labs(title="Chamber excess above local ambient vs measured flux (2023 stems)",
       subtitle="ambient = 10th percentile of the 3 min either side of each enclosure",
       x="C0 - local ambient (ppb)", y=expression(CH[4]~flux~(nmol~m^-2~s^-1))) +
  theme_bw(base_size=9)
ggsave(file.path(outdir,"fig_local_ambient.png"), p, width=8, height=5, dpi=200, bg="white")
cat("\n  Written: outputs/revision/fig_local_ambient.png\n")
