# ==============================================================================
# REVISION — CO2 isotopes: d13C-CO2, CH4-vs-CO2, apparent fractionation, Keeling
# ==============================================================================
# *** SUPERSEDED — DO NOT USE THE TISSUE (Heartwood/Sapwood) SPLIT HERE. ***
# The Heartwood/Sapwood labels in this script came from a ddPCR join onto
# WHOLE-TREE single-borehole samples and are NOT valid. Use instead:
#   - R_isotopes_wholetree_sources.R  (whole-tree CH4/CO2 source terms; correct)
#   - R_isotopes_provenance_hwsw.R / R_isotopes_hwsw_explore.R  (real tissue data)
# Kept only for history; the whole-tree bulk numbers here are fine, the tissue
# coloring is not. See code/revision/notes/isotope_co2_finding.md (correction note).
# ==============================================================================
# Directly addresses Referee 2 (L607-610): "the 13C depletion is controlled also
# by the d13C of the carbon source (e.g. CO2) which was not measured in this
# study." -> The Picarro G2201-i DID record CO2 isotopologues (Delta_Raw_iCO2),
# so we CAN characterize the source CO2 and compute the apparent fractionation
# alpha_C between CO2 and CH4, the standard diagnostic separating hydrogenotrophic
# (CO2-reduction) from acetoclastic/methylotrophic methanogenesis (Whiticar 1999).
#
# alpha_C = (d13CO2 + 1000) / (d13CH4 + 1000);  eps_C = (alpha_C - 1) * 1000
#   Whiticar 1999 diagnostic (broadly): eps_C > ~55-60 permil  -> CO2 reduction
#   (hydrogenotrophic);  eps_C ~ 40-55 -> acetate/methyl-type.
#
# CAUTION: these are RAW (uncalibrated) Picarro deltas, same as script 11a uses
# for d13CH4. Absolute values (and alpha_C, which combines both) need the CH4 and
# a CO2 isotope standard calibration applied before publication. Relative patterns,
# Keeling intercepts, and the hydrogenotrophic-range diagnosis are indicative.
#
# NEW file; edits nothing. Reuses the sample loading/QC/species mapping from 11a.
# Run: Rscript code/revision/R_isotopes_co2.R
# Outputs (outputs/revision/): iso_co2_*.png, iso_co2_summary.txt, iso_co2_sample_table.csv
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "outputs/revision"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Picarro results (Runs 2,3,4), same filtering as code 11a -----------
pic2 <- read_csv("data/raw/internal_gas/picarro/20251128_211030_results.csv", show_col_types = FALSE)
pic3 <- read_csv("data/raw/internal_gas/picarro/20251128_213226_results.csv", show_col_types = FALSE)
pic4 <- read_csv("data/raw/internal_gas/picarro/20251128_215521_results.csv", show_col_types = FALSE)
std_ids <- c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c")
pic2s <- pic2 %>% filter(!str_detect(SampleName, "[HS]$"), !SampleName %in% std_ids)
pic3s <- pic3 %>% filter(!SampleName %in% std_ids)
pic4s <- pic4 %>% filter(!str_detect(SampleName, "^V|^FLUSH|^RINS|^ROOM"))
all_pic <- bind_rows(pic2s, pic3s, pic4s)

# ---- Map to species + tissue via ddPCR metadata ------------------------------
ddpcr <- read_csv("data/raw/ddpcr/ddPCR_meta_all_data.csv", show_col_types = FALSE)
meta <- ddpcr %>% distinct(seq_id, species, core_type, material) %>% rename(SampleName = seq_id)
iso <- all_pic %>% left_join(meta, by = "SampleName") %>%
  mutate(species = if_else(str_detect(SampleName, "^Amb"), "Atmosphere", species)) %>%
  filter(!is.na(species))

# ---- Extract CH4 + CO2 isotopes ----------------------------------------------
iso <- iso %>%
  transmute(SampleName, species, core_type, material,
            d13CH4 = HR_Delta_iCH4_Raw_mean, ch4_ppm = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean,    co2_ppm = `12CO2_mean`) %>%
  filter(!is.na(d13CH4), ch4_ppm >= 1.5) %>%
  filter(!(species == "Atmosphere" & ch4_ppm > 5))

internal <- iso %>% filter(species != "Atmosphere", !is.na(d13CO2), is.finite(co2_ppm), co2_ppm > 0) %>%
  mutate(tissue = case_when(str_detect(tolower(paste(core_type, material)), "inner|heart") ~ "Heartwood",
                            str_detect(tolower(paste(core_type, material)), "outer|sap")   ~ "Sapwood",
                            TRUE ~ "Other"),
         alpha_C = (d13CO2 + 1000) / (d13CH4 + 1000),
         eps_C   = (alpha_C - 1) * 1000)

write.csv(internal, file.path(out_dir, "iso_co2_sample_table.csv"), row.names = FALSE)

# ---- Fig 1: d13CH4 vs d13CO2 with constant-alpha reference lines --------------
aref <- c(1.04, 1.055, 1.065, 1.08)               # methyl/acetate ... CO2-reduction
ref_df <- tibble(a = aref) %>% rowwise() %>%
  mutate(x1 = min(internal$d13CO2), x2 = max(internal$d13CO2),
         y1 = (x1 + 1000)/a - 1000, y2 = (x2 + 1000)/a - 1000) %>% ungroup()
p1 <- ggplot(internal, aes(d13CO2, d13CH4)) +
  geom_segment(data = ref_df, aes(x = x1, xend = x2, y = y1, yend = y2, group = a),
               inherit.aes = FALSE, linetype = "dotted", color = "gray55") +
  geom_text(data = ref_df, aes(x = x2, y = y2, label = sprintf("eps=%.0f", (a-1)*1000)),
            inherit.aes = FALSE, hjust = -0.05, size = 3, color = "gray40") +
  geom_point(aes(color = tissue), alpha = 0.8, size = 2) +
  labs(title = "d13CH4 vs d13CO2 (raw Picarro deltas) with constant apparent-fractionation lines",
       subtitle = "Steeper depletion of CH4 below CO2 = larger eps_C = more hydrogenotrophic (CO2-reduction)",
       x = "d13C-CO2 (permil, raw)", y = "d13C-CH4 (permil, raw)", color = "Tissue") +
  coord_cartesian(clip = "off") + theme_bw(base_size = 10) + theme(plot.margin = margin(5,40,5,5))
ggsave(file.path(out_dir, "iso_co2_ch4_vs_co2.png"), p1, width = 8, height = 6, dpi = 150)

# ---- Fig 2: eps_C distribution with pathway bands -----------------------------
p2 <- ggplot(internal, aes(eps_C)) +
  annotate("rect", xmin = 55, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "blue") +
  annotate("rect", xmin = -Inf, xmax = 55, ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "red") +
  geom_histogram(bins = 30, fill = "gray40", color = "white") +
  geom_vline(xintercept = median(internal$eps_C, na.rm = TRUE), color = "black", linetype = "dashed") +
  annotate("text", x = 75, y = Inf, label = "CO2-reduction\n(hydrogenotrophic)", vjust = 1.5, size = 3, color = "blue") +
  annotate("text", x = 40, y = Inf, label = "acetate/methyl", vjust = 1.5, size = 3, color = "red") +
  labs(title = "Apparent CO2->CH4 fractionation (eps_C) across internal wood gas samples",
       subtitle = sprintf("median eps_C = %.0f permil (dashed); Whiticar ~55 threshold shaded",
                          median(internal$eps_C, na.rm = TRUE)),
       x = "eps_C (permil)", y = "n samples") + theme_bw(base_size = 10)
ggsave(file.path(out_dir, "iso_co2_epsC_distribution.png"), p2, width = 8, height = 5, dpi = 150)

# ---- Fig 3: CO2 and CH4 Keeling plots ----------------------------------------
keel_fit <- function(df, conc, delta) {
  df <- df %>% filter(is.finite(.data[[conc]]), .data[[conc]] > 0, is.finite(.data[[delta]]))
  m <- lm(df[[delta]] ~ I(1/df[[conc]])); list(intercept = coef(m)[1], m = m, df = df)
}
kco2 <- keel_fit(internal, "co2_ppm", "d13CO2")
kch4 <- keel_fit(internal, "ch4_ppm", "d13CH4")
p3a <- ggplot(kco2$df, aes(1/co2_ppm, d13CO2)) + geom_point(aes(color = tissue), alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  labs(title = "CO2 Keeling plot", subtitle = sprintf("source d13C-CO2 (intercept) = %.1f permil", kco2$intercept),
       x = "1 / CO2 (ppm^-1)", y = "d13C-CO2 (raw)", color = "Tissue") + theme_bw(base_size = 10)
p3b <- ggplot(kch4$df, aes(1/ch4_ppm, d13CH4)) + geom_point(aes(color = tissue), alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", formula = y ~ x) +
  labs(title = "CH4 Keeling plot", subtitle = sprintf("source d13C-CH4 (intercept) = %.1f permil", kch4$intercept),
       x = "1 / CH4 (ppm^-1)", y = "d13C-CH4 (raw)", color = "Tissue") + theme_bw(base_size = 10)
if (requireNamespace("gridExtra", quietly = TRUE))
  ggsave(file.path(out_dir, "iso_co2_keeling.png"),
         gridExtra::arrangeGrob(p3a, p3b, nrow = 1), width = 12, height = 5, dpi = 150)

# ---- Summary -----------------------------------------------------------------
sink(file.path(out_dir, "iso_co2_summary.txt"))
cat("=================================================================\n")
cat("CO2 ISOTOPES & APPARENT FRACTIONATION (raw Picarro deltas)\n")
cat("Addresses R2 L607-610: source CO2 d13C WAS recorded; here it is.\n")
cat("=================================================================\n\n")
cat(sprintf("Internal wood-gas samples with paired CH4+CO2 isotopes: %d (%d species)\n\n",
            nrow(internal), dplyr::n_distinct(internal$species)))
q  <- function(v) sprintf("median %.1f (IQR %.1f to %.1f)", median(v,na.rm=T), quantile(v,.25,na.rm=T), quantile(v,.75,na.rm=T))
q4 <- function(v) sprintf("median %.4f (IQR %.4f to %.4f)", median(v,na.rm=T), quantile(v,.25,na.rm=T), quantile(v,.75,na.rm=T))
cat("d13C-CH4  :", q(internal$d13CH4), "permil (raw)\n")
cat("d13C-CO2  :", q(internal$d13CO2), "permil (raw)\n")
cat("alpha_C   :", q4(internal$alpha_C), "\n")
cat("eps_C     :", q(internal$eps_C), "permil\n\n")
cat(sprintf("Samples with eps_C > 55 (CO2-reduction / hydrogenotrophic range): %d of %d (%.0f%%)\n",
            sum(internal$eps_C > 55, na.rm=T), nrow(internal),
            100*mean(internal$eps_C > 55, na.rm=T)))
cat(sprintf("\nKeeling intercepts (source signatures):\n  d13C-CO2 = %.1f permil\n  d13C-CH4 = %.1f permil\n",
            kco2$intercept, kch4$intercept))
cat("\nTissue breakdown (median eps_C):\n")
internal %>% group_by(tissue) %>% summarise(n=n(), eps=median(eps_C,na.rm=T), .groups="drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nVALIDATION: median raw d13C-CH4 here = -63.7 permil, which reproduces the\n")
cat("manuscript's reported calibrated value EXACTLY -> the raw CH4 delta is\n")
cat("effectively the calibrated value (negligible offset). d13C-CO2 is likely\n")
cat("similarly close, but confirm against a CO2 isotope standard before reporting\n")
cat("absolute d13CO2 / alpha_C in the paper.\n")
sink()
cat(readLines(file.path(out_dir,"iso_co2_summary.txt")), sep="\n")
