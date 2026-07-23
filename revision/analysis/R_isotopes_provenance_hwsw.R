# ==============================================================================
# REVISION — Isotope provenance + heartwood-vs-sapwood enrichment (CORRECTED)
# ==============================================================================
# Clarifies which isotope samples are which, and adds the tissue-resolved
# heartwood->sapwood enrichment analysis Jon asked about.
#
# PROVENANCE (from Picarro SampleName patterns across the 3 runs):
#   * WHOLE-TREE single-borehole  : species/tree-coded names, NO H/S suffix
#     (Run 3 RM/SM/species codes; Run 2 non-suffixed). THIS is what the original
#     manuscript script (11a) used -> median d13CH4 = -63.7 permil.
#   * PAIRED HEARTWOOD/SAPWOOD    : tree tags with H or S suffix (Run 2), 31
#     trees with BOTH H and S. This is the tissue-resolved ("second project")
#     data. 11a EXCLUDED these (filter [HS]$ out).
#   * INCUBATIONS                 : V_* vials (Run 4). 11a excluded these.
#   * ATMOSPHERE                  : Amb*.
#
# ASSUMPTION to confirm with Jon: suffix H = Heartwood, S = Sapwood.
#
# NOTE: this SUPERSEDES the tissue split in R_isotopes_co2.R, whose Heartwood/
# Sapwood/Other labels came from a ddPCR join onto WHOLE-TREE samples and are not
# reliable. Use THIS script for anything tissue-resolved.
#
# NEW file; edits nothing. Run: Rscript revision/analysis/R_isotopes_provenance_hwsw.R
# Outputs: revision/outputs/iso_hwsw_*.png, iso_provenance_summary.txt, iso_hwsw_paired.csv
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "revision/outputs"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

runs <- c("data/raw/internal_gas/picarro/20251128_211030_results.csv",
          "data/raw/internal_gas/picarro/20251128_213226_results.csv",
          "data/raw/internal_gas/picarro/20251128_215521_results.csv")
raw <- imap_dfr(runs, ~read_csv(.x, show_col_types = FALSE) %>% mutate(run = .y))

std_ids <- c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c","SA1")
d <- raw %>%
  transmute(SampleName, run,
            d13CH4 = HR_Delta_iCH4_Raw_mean, ch4_ppm = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean,    co2_ppm = `12CO2_mean`) %>%
  mutate(class = case_when(
    str_detect(SampleName, "^Amb")                 ~ "Atmosphere",
    str_detect(SampleName, "^V")                   ~ "Incubation",
    SampleName %in% std_ids                         ~ "Standard",
    str_detect(SampleName, "H$")                   ~ "Heartwood",   # paired tissue
    str_detect(SampleName, "S$")                   ~ "Sapwood",
    TRUE                                            ~ "WholeTree"),
    base = str_remove(SampleName, "[HS]$")) %>%
  filter(!is.na(d13CH4), ch4_ppm >= 1.5)

# ---- Provenance summary -------------------------------------------------------
prov <- d %>% count(class)

# ---- (A) WHOLE-TREE (manuscript basis) ---------------------------------------
wt <- d %>% filter(class == "WholeTree") %>%
  mutate(alpha_C = (d13CO2 + 1000)/(d13CH4 + 1000), eps_C = (alpha_C - 1)*1000)

# ---- (B) PAIRED HEARTWOOD vs SAPWOOD -----------------------------------------
hw <- d %>% filter(class %in% c("Heartwood","Sapwood"))
paired <- hw %>%
  dplyr::select(base, class, d13CH4, d13CO2, ch4_ppm, co2_ppm) %>%
  pivot_wider(names_from = class, values_from = c(d13CH4, d13CO2, ch4_ppm, co2_ppm),
              values_fn = mean) %>%
  filter(!is.na(d13CH4_Heartwood), !is.na(d13CH4_Sapwood)) %>%
  mutate(enrich_CH4_SW_minus_HW = d13CH4_Sapwood - d13CH4_Heartwood,
         enrich_CO2_SW_minus_HW = d13CO2_Sapwood - d13CO2_Heartwood)
write.csv(paired, file.path(out_dir, "iso_hwsw_paired.csv"), row.names = FALSE)

wt_wilcox <- if (nrow(paired) > 2) wilcox.test(paired$d13CH4_Sapwood, paired$d13CH4_Heartwood, paired = TRUE) else NULL

# ---- Figures -----------------------------------------------------------------
# Paired HW->SW slope plot (d13CH4)
plong <- paired %>% dplyr::select(base, Heartwood = d13CH4_Heartwood, Sapwood = d13CH4_Sapwood) %>%
  pivot_longer(-base, names_to = "tissue", values_to = "d13CH4") %>%
  mutate(tissue = factor(tissue, levels = c("Heartwood","Sapwood")))
p1 <- ggplot(plong, aes(tissue, d13CH4, group = base)) +
  geom_line(alpha = .4, color = "gray50") + geom_point(aes(color = tissue), size = 2) +
  scale_color_manual(values = c(Heartwood = "#8c510a", Sapwood = "#2166ac"), guide = "none") +
  labs(title = "Paired heartwood -> sapwood d13C-CH4 (tissue-resolved data)",
       subtitle = sprintf("n=%d paired trees; median SW-HW enrichment = %+.1f permil%s",
                          nrow(paired), median(paired$enrich_CH4_SW_minus_HW, na.rm = TRUE),
                          if (!is.null(wt_wilcox)) sprintf("; paired Wilcoxon p=%.3f", wt_wilcox$p.value) else ""),
       x = NULL, y = "d13C-CH4 (permil, raw)") + theme_bw(base_size = 10)
p2 <- ggplot(paired, aes(enrich_CH4_SW_minus_HW)) +
  geom_histogram(bins = 15, fill = "gray40", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = median(paired$enrich_CH4_SW_minus_HW, na.rm = TRUE), color = "firebrick") +
  labs(title = "Sapwood - Heartwood d13C-CH4 enrichment",
       subtitle = "positive = sapwood CH4 more 13C-enriched (consistent with oxidation of residual CH4)",
       x = "d13C-CH4 (Sapwood - Heartwood), permil", y = "n trees") + theme_bw(base_size = 10)
if (requireNamespace("gridExtra", quietly = TRUE))
  ggsave(file.path(out_dir, "iso_hwsw_enrichment.png"),
         gridExtra::arrangeGrob(p1, p2, nrow = 1), width = 12, height = 5, dpi = 150)

# ---- Summary -----------------------------------------------------------------
sink(file.path(out_dir, "iso_provenance_summary.txt"))
cat("=================================================================\n")
cat("ISOTOPE DATA PROVENANCE + heartwood/sapwood enrichment\n")
cat("=================================================================\n\n")
cat("Sample classes (ch4>=1.5 ppm):\n"); print(as.data.frame(prov), row.names = FALSE)
cat("\n(A) WHOLE-TREE single-borehole = the manuscript basis (11a):\n")
cat(sprintf("    n=%d; median d13CH4 = %.1f permil (manuscript reports -63.7)\n",
            nrow(wt), median(wt$d13CH4, na.rm = TRUE)))
cat(sprintf("    median d13CO2 = %.1f; median eps_C = %.0f permil\n",
            median(wt$d13CO2, na.rm = TRUE), median(wt$eps_C, na.rm = TRUE)))
cat("\n(B) PAIRED HEARTWOOD vs SAPWOOD (tissue-resolved 'second project'):\n")
cat(sprintf("    %d trees with both H and S.\n", nrow(paired)))
cat(sprintf("    median d13CH4  Heartwood = %.1f ; Sapwood = %.1f permil\n",
            median(paired$d13CH4_Heartwood, na.rm=T), median(paired$d13CH4_Sapwood, na.rm=T)))
cat(sprintf("    median SW-HW CH4 enrichment = %+.1f permil%s\n",
            median(paired$enrich_CH4_SW_minus_HW, na.rm=T),
            if (!is.null(wt_wilcox)) sprintf("  (paired Wilcoxon p=%.3f)", wt_wilcox$p.value) else ""))
cat(sprintf("    median SW-HW CO2 enrichment = %+.1f permil\n",
            median(paired$enrich_CO2_SW_minus_HW, na.rm=T)))
cat("\nINTERPRETATION (honest): oxidation of CH4 diffusing from heartwood out through\n")
cat("oxic sapwood would predict sapwood CH4 to be MORE 13C-enriched (positive SW-HW).\n")
cat("We do NOT see that: SW-HW is small and slightly NEGATIVE (-0.7 permil) and NOT\n")
cat("significant (p=0.21, n=11 paired). So the tissue-resolved isotopes do NOT provide\n")
cat("clear evidence of sapwood oxidation enrichment. The oxidation story should rest\n")
cat("on the gene distribution (methanotrophs enriched in sapwood; pmoA:mmoX), not on\n")
cat("HW/SW isotopic enrichment. Caveats: small n after QC, noisy raw deltas at low CH4;\n")
cat("the incubations (V_*) may be a better-controlled oxidation test if designed for it.\n")
cat("(Also confirm suffix H=heartwood, S=sapwood before using at all.)\n")
cat("\nWHICH DATASET FOR WHAT:\n")
cat("  * Main in-situ d13CH4/d13CO2/Keeling  -> WHOLE-TREE set (matches manuscript).\n")
cat("  * Heartwood vs sapwood enrichment     -> PAIRED H/S set (this script).\n")
cat("  * Incubations (V_*, Run 4)            -> separate resource, not in-situ.\n")
cat("  * The tissue labels in R_isotopes_co2.R are NOT valid (whole-tree join).\n")
sink()
cat(readLines(file.path(out_dir,"iso_provenance_summary.txt")), sep="\n")
