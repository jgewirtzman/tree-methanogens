# ==============================================================================
# REVISION — Explore heartwood vs sapwood isotopes (ALL 4 runs, NO conc filter)
# ==============================================================================
# Jon: look at the full tissue set (incl. run 205833 the original omitted) and
# drop the >=1.5 ppm CH4 filter, to see the heartwood-vs-sapwood picture before
# deciding whether to keep or drop the tissue dataset.
#
# NEW file; edits nothing. Run: Rscript code/revision/R_isotopes_hwsw_explore.R
# Outputs: outputs/revision/iso_hwsw_explore.png, iso_hwsw_explore.txt, iso_hwsw_explore_paired.csv
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "outputs/revision"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

runs <- list.files("data/raw/internal_gas/picarro", pattern = "_results.csv$", full.names = TRUE)
runs <- runs[!grepl("merged", runs)]
raw <- map_dfr(runs, ~read_csv(.x, show_col_types = FALSE) %>% mutate(run = basename(.x)))

d <- raw %>%
  transmute(SampleName, run,
            d13CH4 = HR_Delta_iCH4_Raw_mean, ch4 = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean,    co2 = `12CO2_mean`) %>%
  mutate(tissue = case_when(str_detect(SampleName, "H$") ~ "Heartwood",
                            str_detect(SampleName, "S$") ~ "Sapwood", TRUE ~ NA_character_),
         base = str_remove(SampleName, "[HS]$")) %>%
  filter(!is.na(tissue), !is.na(d13CH4))

# pair HW/SW by base tree (mean if replicated), NO concentration filter
pair_at <- function(thr) {
  d %>% filter(ch4 >= thr) %>%
    group_by(base, tissue) %>%
    summarise(d13CH4 = mean(d13CH4), d13CO2 = mean(d13CO2, na.rm = TRUE), ch4 = mean(ch4), .groups = "drop") %>%
    pivot_wider(names_from = tissue, values_from = c(d13CH4, d13CO2, ch4)) %>%
    filter(!is.na(d13CH4_Heartwood), !is.na(d13CH4_Sapwood))
}
p0 <- pair_at(0)     # NO filter
p1 <- pair_at(1.5)   # manuscript filter, for comparison

p0 <- p0 %>% mutate(enrich_CH4 = d13CH4_Sapwood - d13CH4_Heartwood,
                    enrich_CO2 = d13CO2_Sapwood - d13CO2_Heartwood)
write.csv(p0, file.path(out_dir, "iso_hwsw_explore_paired.csv"), row.names = FALSE)

wtest <- function(x) if (sum(is.finite(x)) > 2) wilcox.test(x)$p.value else NA
tt <- function(a,b) if (sum(is.finite(a-b)) > 2) t.test(a, b, paired = TRUE)$p.value else NA

# figure: paired CH4 (no filter) + enrichment hist + CO2 paired
long0 <- p0 %>% transmute(base, Heartwood = d13CH4_Heartwood, Sapwood = d13CH4_Sapwood) %>%
  pivot_longer(-base, names_to = "tissue", values_to = "d13CH4") %>%
  mutate(tissue = factor(tissue, levels = c("Heartwood","Sapwood")))
pA <- ggplot(long0, aes(tissue, d13CH4, group = base)) +
  geom_line(alpha = .35, color = "gray50") + geom_point(aes(color = tissue), size = 2) +
  scale_color_manual(values = c(Heartwood = "#8c510a", Sapwood = "#2166ac"), guide = "none") +
  labs(title = sprintf("Paired d13C-CH4 HW->SW (NO conc filter, n=%d pairs)", nrow(p0)),
       subtitle = sprintf("median SW-HW = %+.1f permil; paired t p=%.3f",
                          median(p0$enrich_CH4, na.rm=T), tt(p0$d13CH4_Sapwood, p0$d13CH4_Heartwood)),
       x = NULL, y = "d13C-CH4 (permil)") + theme_bw(base_size = 10)
pB <- ggplot(p0, aes(enrich_CH4)) + geom_histogram(bins = 20, fill = "gray40", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = median(p0$enrich_CH4, na.rm=T), color = "firebrick") +
  labs(title = "Sapwood - Heartwood d13C-CH4", subtitle = "positive = sapwood enriched (oxidation signature)",
       x = "SW - HW d13C-CH4 (permil)", y = "n pairs") + theme_bw(base_size = 10)
pC <- ggplot(p0, aes(d13CO2_Heartwood, d13CO2_Sapwood)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(alpha = .7) +
  labs(title = "d13C-CO2: heartwood vs sapwood", subtitle = "points above 1:1 = sapwood CO2 more enriched",
       x = "Heartwood d13C-CO2", y = "Sapwood d13C-CO2") + theme_bw(base_size = 10)
if (requireNamespace("gridExtra", quietly = TRUE))
  ggsave(file.path(out_dir, "iso_hwsw_explore.png"),
         gridExtra::arrangeGrob(pA, pB, pC, nrow = 1), width = 14, height = 5, dpi = 150)

# summary
sink(file.path(out_dir, "iso_hwsw_explore.txt"))
cat("===============================================================\n")
cat("HEARTWOOD vs SAPWOOD ISOTOPES — all 4 runs, NO concentration filter\n")
cat("===============================================================\n\n")
cat(sprintf("Total tissue samples (HW+SW, d13CH4 present): %d\n", nrow(d)))
cat(sprintf("  Heartwood: %d, Sapwood: %d\n", sum(d$tissue=="Heartwood"), sum(d$tissue=="Sapwood")))
cat(sprintf("Paired trees (both tissues): NO filter = %d ; >=1.5 ppm = %d\n\n", nrow(p0), nrow(p1)))
cat("d13C-CH4 (permil):\n")
cat(sprintf("  Heartwood median = %.1f ; Sapwood median = %.1f\n",
            median(p0$d13CH4_Heartwood,na.rm=T), median(p0$d13CH4_Sapwood,na.rm=T)))
cat(sprintf("  SW-HW enrichment: median %+.1f ; paired t p=%.3f ; Wilcoxon(SW-HW) p=%.3f\n",
            median(p0$enrich_CH4,na.rm=T), tt(p0$d13CH4_Sapwood,p0$d13CH4_Heartwood), wtest(p0$enrich_CH4)))
cat(sprintf("  n pairs SW more enriched than HW: %d / %d\n\n",
            sum(p0$enrich_CH4 > 0, na.rm=T), sum(is.finite(p0$enrich_CH4))))
cat("d13C-CO2 (permil):\n")
cat(sprintf("  Heartwood median = %.1f ; Sapwood median = %.1f\n",
            median(p0$d13CO2_Heartwood,na.rm=T), median(p0$d13CO2_Sapwood,na.rm=T)))
cat(sprintf("  SW-HW enrichment: median %+.1f ; paired t p=%.3f\n\n",
            median(p0$enrich_CO2,na.rm=T), tt(p0$d13CO2_Sapwood,p0$d13CO2_Heartwood)))
cat("CAUTION: without the conc filter, many tissue samples are near detection\n")
cat("(HW median ~1.3 ppm, SW ~0.9 ppm CH4), so raw deltas are noisier. Read the\n")
cat("sign/consistency, not precise magnitudes. Decide keep-vs-drop from this.\n")
sink()
cat(readLines(file.path(out_dir,"iso_hwsw_explore.txt")), sep="\n")
