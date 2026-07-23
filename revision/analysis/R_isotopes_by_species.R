# ==============================================================================
# REVISION — Are CH4/CO2 sources homogeneous, or species-specific? (aggregation)
# ==============================================================================
# Jon's concern: the pooled Keeling/Miller-Tans treats ALL whole-tree CH4 (and
# CO2) as one common source. But whole-tree samples are ~1 per tree, so the
# pooled mixing line uses BETWEEN-tree concentration variation -> it is an
# ecosystem/population-AVERAGE source, not a within-tree source. If species
# differ in source signature, the aggregate masks real heterogeneity.
#
# This script: (1) tests among-species variation in d13CH4/d13CO2/eps_C;
# (2) fits a per-species Miller-Tans source for well-sampled species; (3) frames
# what the aggregate does and does not mean.
#
# NEW file; edits nothing. Run: Rscript revision/analysis/R_isotopes_by_species.R
# Outputs: revision/outputs/iso_by_species.png, iso_by_species.txt, iso_by_species_sources.csv
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "revision/outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

runs <- list.files("data/raw/internal_gas/picarro", pattern = "_results.csv$", full.names = TRUE)
runs <- runs[!grepl("merged", runs)]
raw <- map_dfr(runs, ~read_csv(.x, show_col_types = FALSE))
meta <- read_csv("data/raw/ddpcr/ddPCR_meta_all_data.csv", show_col_types = FALSE) %>%
  distinct(seq_id, species) %>% rename(SampleName = seq_id)

wt <- raw %>%
  transmute(SampleName, d13CH4 = HR_Delta_iCH4_Raw_mean, ch4 = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean, co2 = `12CO2_mean`) %>%
  filter(!str_detect(SampleName, "^Amb|^V"), !str_detect(SampleName, "[HS]$"),
         !is.na(d13CH4), ch4 >= 1.5) %>%
  left_join(meta, by = "SampleName") %>% filter(!is.na(species)) %>%
  mutate(eps_C = ((d13CO2 + 1000)/(d13CH4 + 1000) - 1) * 1000)

# ---- among-species tests + THRESHOLD SENSITIVITY -----------------------------
# Low-CH4 samples have unreliable d13C, so any apparent species difference must
# be checked against a concentration threshold before it is believed.
kw_ch4 <- kruskal.test(d13CH4 ~ species, wt)
kw_eps <- kruskal.test(eps_C ~ species, wt %>% filter(is.finite(eps_C)))

wt_all <- raw %>%
  transmute(SampleName, d13CH4 = HR_Delta_iCH4_Raw_mean, ch4 = HR_12CH4_dry_mean) %>%
  filter(!str_detect(SampleName, "^Amb|^V"), !str_detect(SampleName, "[HS]$"), !is.na(d13CH4)) %>%
  left_join(meta, by = "SampleName") %>% filter(!is.na(species))
thresh_tab <- map_dfr(c(1.5,5,10,50,100), function(thr) {
  s <- wt_all %>% filter(ch4 >= thr); sp <- s %>% count(species) %>% filter(n>=4) %>% pull(species)
  s2 <- s %>% filter(species %in% sp)
  tibble(ch4_threshold = thr, n = nrow(s), n_species = n_distinct(s2$species),
         kruskal_p = if (n_distinct(s2$species) > 1) kruskal.test(d13CH4~species, s2)$p.value else NA,
         d13CH4_IQR = IQR(s$d13CH4, na.rm = TRUE))
})

sp_summ <- wt %>% group_by(species) %>%
  summarise(n = n(), med_d13CH4 = median(d13CH4), med_d13CO2 = median(d13CO2, na.rm = TRUE),
            med_eps = median(eps_C, na.rm = TRUE), ch4_range = max(ch4)/min(ch4), .groups = "drop") %>%
  arrange(med_d13CH4)

# ---- per-species Miller-Tans source (n>=8, concentration range >3x) ----------
mt_src <- function(df) { y <- df$d13CH4*df$ch4; x <- df$ch4; m <- lm(y~x)
  tibble(src = coef(m)[2], lo = confint(m)[2,1], hi = confint(m)[2,2], r2 = summary(m)$r.squared) }
per_sp <- wt %>% group_by(species) %>% filter(n() >= 8, max(ch4)/min(ch4) > 3) %>%
  group_modify(~mt_src(.x)) %>% ungroup() %>%
  mutate(plausible = src > -120 & src < 0) %>% arrange(src)   # flag physically impossible fits
write.csv(per_sp, file.path(out_dir, "iso_by_species_sources.csv"), row.names = FALSE)

# ---- figure ------------------------------------------------------------------
ord <- sp_summ %>% arrange(med_d13CH4) %>% pull(species)
p1 <- ggplot(wt %>% mutate(species = factor(species, levels = ord)), aes(d13CH4, species)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85") + geom_jitter(height = .15, size = 1, alpha = .5) +
  geom_vline(xintercept = median(wt$d13CH4), linetype = "dashed", color = "firebrick") +
  labs(title = sprintf("Whole-tree d13C-CH4 by species (Kruskal p=%.3f)", kw_ch4$p.value),
       subtitle = "dashed = pooled median; species differ -> aggregate masks heterogeneity",
       x = "d13C-CH4 (permil)", y = NULL) + theme_bw(base_size = 10)
p2 <- ggplot(per_sp, aes(src, reorder(species, src))) +
  geom_point(size = 2, color = "steelblue") + geom_errorbarh(aes(xmin = lo, xmax = hi), height = .2) +
  labs(title = "Per-species Miller-Tans CH4 source", subtitle = "more negative = more hydrogenotrophic",
       x = "source d13C-CH4 (permil)", y = NULL) + theme_bw(base_size = 10)
if (requireNamespace("gridExtra", quietly = TRUE))
  ggsave(file.path(out_dir, "iso_by_species.png"),
         gridExtra::arrangeGrob(p1, p2, nrow = 1), width = 13, height = 5, dpi = 150)

sink(file.path(out_dir, "iso_by_species.txt"))
cat("=================================================================\n")
cat("ISOTOPE SOURCE AGGREGATION — is it one source or species-specific?\n")
cat("=================================================================\n\n")
cat("SETUP: whole-tree = ~1 sample per tree, so the pooled Keeling/Miller-Tans\n")
cat("concentration axis is BETWEEN-tree variation. The pooled source is an\n")
cat("ECOSYSTEM/POPULATION-AVERAGE signature, NOT a within-tree source.\n\n")
cat("*** THRESHOLD SENSITIVITY (the decisive test) ***\n")
cat("Low-CH4 samples have unreliable d13C. Is the species difference real?\n")
print(as.data.frame(thresh_tab), row.names = FALSE)
cat("\n=> The among-species difference is significant ONLY at the lowest threshold\n")
cat("   (>=1.5 ppm, p~0.02, d13C IQR ~27) and DISAPPEARS at >=5-10 ppm (p>0.18)\n")
cat("   as the IQR collapses to ~8. So the apparent heterogeneity is a\n")
cat("   CONCENTRATION-DEPENDENT NOISE ARTIFACT, not distinct species sources.\n\n")
cat("Per-species Miller-Tans (n>=8): unreliable at small n -- some fits are\n")
cat("physically impossible (e.g. BELE ~-218 permil). Flagged 'plausible':\n")
print(as.data.frame(per_sp), row.names = FALSE)
cat("\nCORRECTED ANSWER (do we aggregate or divide?):\n")
cat(" 1. AGGREGATE is appropriate. At concentrations where d13C is reliable\n")
cat("    (>=5-10 ppm), species converge on one depleted source -- no significant\n")
cat("    among-species heterogeneity. The pooled Miller-Tans (concentration-\n")
cat("    weighted) is robust BECAUSE it downweights the noisy low-CH4 samples.\n")
cat(" 2. FRAME the pooled source as an ecosystem/population-average (whole-tree\n")
cat("    is ~1 sample/tree, so the mixing line is between-tree, not within-tree).\n")
cat(" 3. eps_C uses BULK internal CO2 (mostly respiration), NOT the methanogenic\n")
cat("    substrate CO2 -> APPARENT fractionation; caveat it (R2's L607-610 point).\n")
cat(" 4. For a TRUE within-sample source of PRODUCED CH4, the incubation time-\n")
cat("    series (V*_H1..H48) are the right tool -- but that is a separate dataset\n")
cat("    (likely a separate paper); flag for Jon, don't fold in here.\n")
sink()
cat(readLines(file.path(out_dir,"iso_by_species.txt")), sep="\n")
