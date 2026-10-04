#!/usr/bin/env Rscript
# ==============================================================================
# 24_ddpcr-detection-precision.R -- mcrA detection and counting precision by threshold
# (table for SI Methods S3, §8)
# ------------------------------------------------------------------------------
# Digital PCR counting error is Poisson: CV of a concentration ~ 1/sqrt(positive
# droplets). For each threshold, reports the positive droplets it corresponds to, the
# CV, the same error on the log10 scale, and copies per gram at reference conditions
# (16,000 accepted droplets; 100 mg wood, 250 mg soil; as the Methods S3 §7 example).
# The S3 table shows these assay columns only. The shares of survey samples at or
# above each threshold (probe mcrA, loose, one well per sample) are kept here for
# reference; the main text reports only detection and the limit of detection.
# The limit of detection (3 copies per reaction, 0.15 copies/µL) is judged on
# concentration, as in the main text; only ~2/3 of a 20 µL well is read, so it is about
# 2 positive droplets. 8 positive droplets (CV 35%) is the qPCR-style limit of
# quantification of Forootan et al. 2017.
# Writes outputs/data/ddpcr_detection_precision.csv
# ==============================================================================
suppressMessages(library(dplyr))
source("code/lib/ddpcr_constants.R")
source("code/lib/outputs.R")

DROPLET_NL <- 0.85
REF_DROPLETS <- 16000
REF_MASS_MG <- c(Wood = 100, Soil = 250)
LOD_CONC <- 3 / 20   # copies per µL reaction
read_ul <- REF_DROPLETS * DROPLET_NL / 1000
conc_at <- function(k) -log(1 - k / REF_DROPLETS) / (DROPLET_NL / 1000)

cg <- read.csv("data/compiled/ddpcr_gene_abundances.csv") %>%
  filter(analysis_type == "loose", target_gene == "mcra_probe", !is.na(concentration_copies_per_uL))
stopifnot(!any(duplicated(cg$sample_id)))

rows <- data.frame(threshold = c("single droplet", "limit of detection", "", "limit of quantification", ""),
                   k = c(1, NA, 4, 8, 16))
res <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
  k <- rows$k[i]; lod <- is.na(k)
  conc <- if (lod) LOD_CONC else conc_at(k)
  droplets <- if (lod) conc * read_ul else k
  share <- sapply(c(Inner = "Inner", Outer = "Outer", Mineral = "Mineral", Organic = "Organic"), function(ct) {
    d <- cg[cg$core_type == ct, ]
    round(100 * mean(if (lod) d$concentration_copies_per_uL >= LOD_CONC else d$positives >= k), 1)
  })
  data.frame(threshold = rows$threshold[i], positive_droplets = round(droplets, 1),
             cv_pct = round(100 / sqrt(droplets)), log10_error = round(log10(1 + 1 / sqrt(droplets)), 2),
             copies_per_g_wood_100mg = round(ddpcr_copies_per_g(conc, REF_MASS_MG[["Wood"]], "Wood"), -1),
             copies_per_g_soil_250mg = round(ddpcr_copies_per_g(conc, REF_MASS_MG[["Soil"]], "Soil"), -1),
             pct_heartwood = share[["Inner"]], pct_sapwood = share[["Outer"]],
             pct_mineral = share[["Mineral"]], pct_organic = share[["Organic"]])
}))
write.csv(res, out_path("ddpcr_detection_precision.csv"), row.names = FALSE)
print(res)
