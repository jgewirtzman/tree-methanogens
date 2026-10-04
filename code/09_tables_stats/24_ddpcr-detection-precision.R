#!/usr/bin/env Rscript
# ==============================================================================
# 24_ddpcr-detection-precision.R -- mcrA detection and counting precision by threshold
# (table for SI Methods S3, §8)
# ------------------------------------------------------------------------------
# Digital PCR counting error is Poisson: CV of a concentration ~ 1/sqrt(positive
# droplets). For each threshold, reports the positive droplets it corresponds to, the
# CV, the same error on the log10 scale, the share of samples at or above it by
# compartment, and copies per gram at the median accepted droplets and sample mass.
# The limit of detection (3 copies per reaction, 0.15 copies/µL) is judged on
# concentration, as in the main text; only ~2/3 of a 20 µL well is read, so it is about
# 2 positive droplets. 8 positive droplets (CV 35%) is the qPCR-style limit of
# quantification of Forootan et al. 2017. Probe mcrA, loose thresholds, one well per sample.
# Writes outputs/data/ddpcr_detection_precision.csv
# ==============================================================================
suppressMessages(library(dplyr))
source("code/lib/ddpcr_constants.R")
source("code/lib/outputs.R")

DROPLET_NL <- 0.85
THRESHOLDS <- c(1, 4, 8, 16)
cg <- read.csv("data/compiled/ddpcr_gene_abundances.csv") %>%
  filter(analysis_type == "loose", target_gene == "mcra_probe", !is.na(concentration_copies_per_uL))
stopifnot(!any(duplicated(cg$sample_id)))

LOD_CONC <- 3 / 20   # copies per µL reaction
rows <- c(list(list(label = "single droplet", k = 1)), list(list(label = "limit of detection", k = NA)),
          lapply(THRESHOLDS[-1], function(k) list(label = if (k == 8) "limit of quantification" else "", k = k)))
res <- do.call(rbind, lapply(rows, function(r) {
  do.call(rbind, lapply(c("Inner", "Outer", "Mineral", "Organic"), function(ct) {
    d <- cg[cg$core_type == ct, ]; read_ul <- median(d$accepted_droplets) * DROPLET_NL / 1000
    if (is.na(r$k)) {          # LoD on concentration; droplets = expected positives in the read volume
      conc <- LOD_CONC; k <- conc * read_ul; above <- d$concentration_copies_per_uL >= conc
    } else {
      k <- r$k; conc <- -log(1 - k / median(d$accepted_droplets)) / (DROPLET_NL / 1000); above <- d$positives >= k
    }
    data.frame(threshold = r$label, positive_droplets = round(k, 1), cv_pct = round(100 / sqrt(k)),
               log10_error = round(log10(1 + 1 / sqrt(k)), 2), core_type = ct, n = nrow(d),
               pct_at_or_above = round(100 * mean(above), 1),
               copies_per_g = round(ddpcr_copies_per_g(conc, median(d$sample_mass_mg, na.rm = TRUE), ct), -2))
  }))
}))
write.csv(res, out_path("ddpcr_detection_precision.csv"), row.names = FALSE)
print(res)
