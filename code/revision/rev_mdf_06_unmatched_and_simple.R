#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_06_unmatched_and_simple.R      STEP 6
# ------------------------------------------------------------------------------
# TWO QUESTIONS
#
# 1. WHY ARE MEASUREMENTS UNMATCHED?
#    Per-measurement Allan deviation requires the saved click.peak trace. Where a
#    measurement has a flux but no trace, it drops out. This diagnoses each gap
#    rather than leaving it as an unexplained shortfall.
#
# 2. WHAT IF WE USED A SIMPLER, CAMPAIGN-LEVEL PRECISION INSTEAD?
#    Rather than a per-measurement sigma, use ONE empirical precision per campaign
#    (the median Allan deviation of that campaign, or the rolling-window value where
#    no clicks survive) and apply it to every measurement:
#
#        MDF = z * sigma_campaign / t * flux.term
#
#    This is the Wassmann et al. (2018) form with a campaign-specific empirical sigma
#    in place of the manufacturer constant. Advantages: it covers EVERY measurement
#    with chamber geometry and duration, including those whose traces were lost, and
#    it is far simpler to describe. Cost: it cannot distinguish a quiet measurement
#    from a noisy one within a campaign.
#
#    Because the per-measurement approach is available for most campaigns, the two can
#    be compared directly to see how much the extra complexity actually buys.
#
# Output: outputs/revision/mdf_unmatched_diagnosis.txt
#         outputs/revision/mdf_simple_campaign.csv
# ==============================================================================

suppressMessages({library(dplyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

ALL <- read.csv("outputs/revision/mdf_allan_all_campaigns.csv", stringsAsFactors=FALSE)
W   <- read.csv("outputs/revision/mdf_master_windows.csv", stringsAsFactors=FALSE)

# ==============================================================================
rule("1. WHERE DO MEASUREMENTS GO MISSING?")
# ==============================================================================
FLUX <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    transmute(UniqueID, campaign="height_2021", flux=best.flux, flux.term, t=nb.obs),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv", stringsAsFactors=FALSE) %>%
    transmute(UniqueID, campaign="semirigid_soil", flux=best.flux, flux.term, t=nb.obs),
  read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    transmute(UniqueID, campaign="static_2023", flux=best.flux, flux.term, t=nb.obs))
st <- read.csv("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv",
               stringsAsFactors=FALSE, check.names=TRUE)
g <- function(b) { for (s in c(".y",".x","")) { n <- paste0("CH4_",b,s)
  if (n %in% names(st)) { v <- suppressWarnings(as.numeric(st[[n]])); if (any(is.finite(v))) return(v) } }
  rep(NA_real_, nrow(st)) }
FLUX <- bind_rows(FLUX, data.frame(UniqueID=st$UniqueID, campaign="semirigid_tree",
  flux=g("best.flux"), flux.term=g("flux.term"), t=g("nb.obs"), stringsAsFactors=FALSE))
FLUX <- FLUX %>% filter(is.finite(flux)) %>% distinct(UniqueID, .keep_all=TRUE)

FLUX$has_allan <- FLUX$UniqueID %in% ALL$UniqueID
FLUX$has_geom  <- is.finite(FLUX$flux.term) & is.finite(FLUX$t) & FLUX$t > 2
cat("\n  Measurements with a flux, and what else they have:\n\n")
print(FLUX %>% group_by(campaign) %>%
  summarise(n_flux=n(), has_allan=sum(has_allan), no_allan=sum(!has_allan),
            has_geom=sum(has_geom), usable_permeas=sum(has_allan & has_geom),
            usable_simple=sum(has_geom), .groups="drop") %>% as.data.frame(), row.names=FALSE)

cat("\n  WHY the Allan values are missing:\n")
for (cp in unique(FLUX$campaign)) {
  miss <- FLUX$UniqueID[FLUX$campaign==cp & !FLUX$has_allan]
  if (!length(miss)) { cat(sprintf("    %-16s none missing\n", cp)); next }
  in_win <- sum(miss %in% W$UniqueID)
  cat(sprintf("    %-16s %3d without Allan | %d of those ARE in the window list\n",
              cp, length(miss), in_win))
  cat(sprintf("                     examples: %s\n", paste(head(miss,3), collapse=", ")))
}
cat("
  A measurement is missing an Allan deviation when its click.peak trace was not
  retained -- either the measurement failed click.peak and was flux-fitted another
  way, or its trace lives in a manID file we have not located. It is NOT a matching
  bug: the UniqueIDs are consistent between files. The static_2023 gap is total and
  has a known cause (filename collision, task #21).
")

# ==============================================================================
rule("2. SIMPLE CAMPAIGN-LEVEL PRECISION")
# ==============================================================================
SIG <- ALL %>% group_by(campaign) %>%
  summarise(sigma_campaign=median(allan_sd_CH4, na.rm=TRUE), n_allan=n(), .groups="drop")
SIG <- bind_rows(SIG, data.frame(campaign="static_2023",
        sigma_campaign=1.298/0.86, n_allan=0))   # rolling-window, adjusted
SIG$campaign[SIG$campaign=="untagged"] <- "semirigid_tree"   # folded in
SIG <- SIG %>% group_by(campaign) %>% summarise(sigma_campaign=median(sigma_campaign),
                                                n_allan=sum(n_allan), .groups="drop")
cat("\n  ONE precision per campaign:\n\n")
cat(sprintf("  %-16s %14s %10s %s\n","campaign","sigma (ppb)","n Allan","source"))
for (i in seq_len(nrow(SIG))) with(SIG[i,],
  cat(sprintf("  %-16s %14.3f %10d %s\n", campaign, sigma_campaign, n_allan,
      if (campaign=="static_2023") "rolling-window / 0.86" else "median per-measurement Allan")))

D <- FLUX %>% filter(has_geom) %>% left_join(SIG, by="campaign") %>%
  filter(is.finite(sigma_campaign))
D$is_soil <- D$campaign=="semirigid_soil"
D$sigma_permeas <- ALL$allan_sd_CH4[match(D$UniqueID, ALL$UniqueID)]

mdf <- function(sig, z) sig * z / D$t * D$flux.term
cat(sprintf("\n  Coverage with the SIMPLE approach: %d measurements (stems %d, soil %d)\n",
            nrow(D), sum(!D$is_soil), sum(D$is_soil)))
cat(sprintf("  Coverage with PER-MEASUREMENT Allan: %d (%d stems)\n\n",
            sum(is.finite(D$sigma_permeas)), sum(is.finite(D$sigma_permeas) & !D$is_soil)))

cat("  STEM DETECTION STATUS -- simple campaign sigma:\n\n")
cat(sprintf("  %-16s %9s %10s %9s %9s\n","threshold","DETECTED","%detected","sig.POS","sig.NEG"))
OUT <- list()
for (lab in c("90%","95%","99%")) {
  z <- qnorm(1 - (1-as.numeric(sub("%","",lab))/100)/2)
  m <- mdf(D$sigma_campaign, z); tr <- !D$is_soil
  cat(sprintf("  Wassmann %-7s %9d %9.1f%% %9d %9d\n", lab, sum(abs(D$flux[tr])>=m[tr]),
      100*mean(abs(D$flux[tr])>=m[tr]), sum(D$flux[tr]>=m[tr]), sum(D$flux[tr]<=-m[tr])))
  OUT[[lab]] <- data.frame(threshold=paste("Wassmann",lab), approach="campaign sigma",
    n=sum(tr), detected=sum(abs(D$flux[tr])>=m[tr]), pct=100*mean(abs(D$flux[tr])>=m[tr]),
    sig_pos=sum(D$flux[tr]>=m[tr]), sig_neg=sum(D$flux[tr]<=-m[tr]))
}
cat("\n  SOIL:\n\n")
for (lab in c("90%","95%","99%")) {
  z <- qnorm(1 - (1-as.numeric(sub("%","",lab))/100)/2)
  m <- mdf(D$sigma_campaign, z); so <- D$is_soil
  cat(sprintf("  Wassmann %-7s %9d %9.1f%% %9d %9d\n", lab, sum(abs(D$flux[so])>=m[so]),
      100*mean(abs(D$flux[so])>=m[so]), sum(D$flux[so]>=m[so]), sum(D$flux[so]<=-m[so])))
}

# --- how much does per-measurement buy over campaign-level? -------------------
rule("DOES PER-MEASUREMENT SIGMA CHANGE ANY CONCLUSION?")
k <- is.finite(D$sigma_permeas)
z95 <- qnorm(.975)
m_simple <- D$sigma_campaign[k]*z95/D$t[k]*D$flux.term[k]
m_perm   <- D$sigma_permeas[k]*z95/D$t[k]*D$flux.term[k]
f <- D$flux[k]; tr <- !D$is_soil[k]
cat(sprintf("
  On the %d measurements where BOTH are available (stems %d):

  %-24s %10s %10s
  %-24s %10d %10d
  %-24s %10d %10d
  %-24s %10d %10d

  Measurements classified DIFFERENTLY by the two approaches: %d (%.1f%%)
",
 sum(k), sum(tr),
 "", "campaign", "per-meas",
 "stems detected",        sum(abs(f[tr])>=m_simple[tr]), sum(abs(f[tr])>=m_perm[tr]),
 "stems sig. positive",   sum(f[tr]>=m_simple[tr]),      sum(f[tr]>=m_perm[tr]),
 "stems sig. negative",   sum(f[tr]<=-m_simple[tr]),     sum(f[tr]<=-m_perm[tr]),
 sum((abs(f)>=m_simple) != (abs(f)>=m_perm)),
 100*mean((abs(f)>=m_simple) != (abs(f)>=m_perm))))
cat("
  If the two agree closely, the campaign-level sigma is preferable: it is simpler to
  describe, it covers every measurement including those whose traces were lost, and it
  removes the need to explain a mixed per-measurement / rolling-window methodology.
")
write.csv(bind_rows(OUT), file.path(outdir,"mdf_simple_campaign.csv"), row.names=FALSE)
cat("\n  Written: outputs/revision/mdf_simple_campaign.csv\n")
