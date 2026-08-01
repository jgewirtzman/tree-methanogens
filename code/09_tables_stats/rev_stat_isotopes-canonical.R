# ==============================================================================
# ISOTOPES — FINAL settled workflow (NPH-MS-2026-56441 revision)
# ==============================================================================
# Supersedes the exploratory isotope scripts. Self-contained. Produces:
#   outputs/revision/ISOTOPES_methods.md      methods text
#   outputs/revision/ISOTOPES_results.md      results text (numbers filled in)
#   outputs/revision/ISOTOPES_fig1_source.png Keeling + convergence panels
#   outputs/revision/ISOTOPES_fig2_crossplot.png  d13CH4 vs d13CO2 + alpha_C isolines
#   outputs/revision/ISOTOPES_sample_table.csv    per-sample data
#
# SETTLED DECISIONS (from the revision thread):
#  * Use WHOLE-TREE single-borehole samples (manuscript basis). Exclude
#    calibration standards, atmosphere (Amb), incubations (V*), and the paired
#    heartwood/sapwood tissue set (dropped: no significant HW->SW enrichment).
#  * CH4 source: Keeling intercept (robust here) + two-endmember atmosphere
#    mixing correction. Report Miller-Tans too but flag it is leverage-sensitive.
#  * eps_C computed WITHIN-tree (each sample's own CH4+CO2); report bulk and
#    source-corrected. Report the CO2-CH4 coupling check (bulk CO2 is respiration-
#    dominated -> eps_C is an APPARENT, corroborating fractionation).
#  * Single delta-13C only (no dD-CH4): a stated limitation.
#
# Run: Rscript code/revision/ISOTOPES_final.R
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(gridExtra) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- 1. LOAD + QC ------------------------------------------------------------
runs <- list.files("data/raw/internal_gas/picarro", pattern = "_results.csv$", full.names = TRUE)
runs <- runs[!grepl("merged", runs)]
raw  <- map_dfr(runs, ~read_csv(.x, show_col_types = FALSE))
STANDARDS <- c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c","SA1")

d <- raw %>%
  transmute(SampleName,
            d13CH4 = HR_Delta_iCH4_Raw_mean, ch4_ppm = HR_12CH4_dry_mean,
            d13CO2 = Delta_Raw_iCO2_mean,     co2_ppm = `12CO2_mean`) %>%
  filter(!str_detect(SampleName, "^Amb"),        # atmosphere
         !str_detect(SampleName, "^V"),          # incubations (unrelated)
         !str_detect(SampleName, "[HS]$"),       # paired tissue set (dropped)
         !SampleName %in% STANDARDS,             # calibration standards
         !is.na(d13CH4), ch4_ppm >= 1.5)         # manuscript isotope QC floor
wt <- d %>% filter(!is.na(d13CO2), co2_ppm > 0) %>%
  mutate(alpha_C = (d13CO2 + 1000)/(d13CH4 + 1000), eps_C = (alpha_C - 1)*1000,
         invCH4 = 1/ch4_ppm)
write.csv(wt, file.path(out, "ISOTOPES_sample_table.csv"), row.names = FALSE)

# ---- 2. BULK DESCRIPTIVES ----------------------------------------------------
med <- function(x) median(x, na.rm = TRUE)
iqr <- function(x) paste(round(quantile(x, c(.25,.75), na.rm = TRUE),1), collapse = " to ")
bulk_ch4 <- med(d$d13CH4); bulk_co2 <- med(wt$d13CO2); bulk_eps <- med(wt$eps_C)

# ---- 3. CH4 SOURCE -----------------------------------------------------------
# (a) Keeling: d13CH4 ~ 1/CH4, intercept = source (robust here)
mk <- lm(d13CH4 ~ invCH4, data = d %>% mutate(invCH4 = 1/ch4_ppm))
keel_src <- unname(coef(mk)[1]); keel_ci <- confint(mk)[1,]
# leverage check: drop the k highest-CH4 points
kdrop <- d %>% arrange(desc(ch4_ppm)) %>% { sapply(c(0,3,5), function(k)
  unname(coef(lm(d13CH4 ~ I(1/ch4_ppm), data = slice(., (k+1):n())))[1])) }
# (b) Miller-Tans (report + flag leverage)
mt_src <- unname(coef(lm((d13CH4*ch4_ppm) ~ ch4_ppm, data = d))[2])
mt_drop <- d %>% arrange(desc(ch4_ppm)) %>% { sapply(c(0,3,5), function(k)
  unname(coef(lm((d13CH4*ch4_ppm) ~ ch4_ppm, data = slice(., (k+1):n())))[2])) }
# (c) two-endmember atmosphere mixing correction (per sample), endmember-robust
C_ATM <- 1.9; D_ATM <- -47
atmcorr <- function(datm) d %>% mutate(tc = (ch4_ppm*d13CH4 - C_ATM*datm)/(ch4_ppm - C_ATM)) %>%
  filter(ch4_ppm - C_ATM > 0.5) %>% pull(tc) %>% med()
atm_src <- atmcorr(D_ATM); atm_sens <- sapply(c(-47,-40,-34.8), atmcorr)

# ---- 4. eps_C (within-tree) + source-value eps_C + coupling -------------------
eps_bulk <- bulk_eps
eps_from <- function(dch4, dco2) ((dco2 + 1000)/(dch4 + 1000) - 1)*1000
eps_src_keel <- eps_from(keel_src, bulk_co2)          # source d13CH4 + median CO2
eps_src_atm  <- eps_from(atm_src,  bulk_co2)
# coupling: is bulk CO2 tied to CH4 (substrate-product) or decoupled?
cpl_all <- cor.test(wt$d13CH4, wt$d13CO2)
hi <- wt %>% filter(ch4_ppm >= 10); cpl_hi <- cor.test(hi$d13CH4, hi$d13CO2)
# among-species homogeneity (is one pooled estimate ok?) at reliable conc
meta_sp <- read_csv("data/raw/ddpcr/ddPCR_meta_all_data.csv", show_col_types = FALSE) %>%
  distinct(seq_id, species) %>% rename(SampleName = seq_id)
hi_join <- hi %>% left_join(meta_sp, by = "SampleName")
kw_matched <- sum(!is.na(hi_join$species)); kw_total <- nrow(hi)   # join-count check (reviewer #4)
sp10 <- hi_join %>% filter(!is.na(species)) %>% group_by(species) %>% filter(n() >= 4) %>% ungroup()
kw_ngrp <- n_distinct(sp10$species); kw_n <- nrow(sp10)
kw_p <- if (kw_ngrp > 1) kruskal.test(d13CH4 ~ species, sp10)$p.value else NA

# ---- 5. FIGURES --------------------------------------------------------------
# Fig 1: Keeling + convergence
p_keel <- ggplot(d, aes(1/ch4_ppm, d13CH4)) + geom_point(alpha=.5, size=1.3) +
  geom_smooth(method="lm", se=TRUE, color="black", formula=y~x) +
  annotate("point", x=0, y=keel_src, color="firebrick", size=3) +
  labs(title="(a) CH4 Keeling plot", subtitle=sprintf("source intercept = %.0f permil [%.0f, %.0f]", keel_src, keel_ci[1], keel_ci[2]),
       x="1 / CH4 (ppm^-1)", y=expression(delta^13*C-CH[4]~"(permil)")) + theme_bw(base_size=10)
p_conv_ch4 <- ggplot(d, aes(ch4_ppm, d13CH4)) + geom_point(alpha=.5, size=1.3) + scale_x_log10() +
  geom_smooth(method="loess", se=TRUE, color="grey30", formula=y~x) +
  geom_hline(yintercept=atm_src, color="firebrick", linetype="dashed") +
  labs(title="(b) CH4 convergence", subtitle=sprintf("plateau ~ %.0f (atmosphere-corrected)", atm_src),
       x="internal CH4 (ppm, log)", y=expression(delta^13*C-CH[4]~"(permil)")) + theme_bw(base_size=10)
p_conv_eps <- ggplot(wt, aes(ch4_ppm, eps_C)) + geom_point(alpha=.5, size=1.3) + scale_x_log10() +
  geom_smooth(method="loess", se=TRUE, color="steelblue", formula=y~x) +
  coord_cartesian(ylim=c(-50,120)) +
  labs(title="(c) apparent fractionation convergence", subtitle=sprintf("plateau eps_C ~ %.0f-%.0f permil", eps_src_atm, eps_src_keel),
       x="internal CH4 (ppm, log)", y=expression(epsilon[C]~"(permil)")) + theme_bw(base_size=10)
ggsave(file.path(out,"ISOTOPES_fig1_source.png"), arrangeGrob(p_keel, p_conv_ch4, p_conv_eps, nrow=1),
       width=13.5, height=4.3, dpi=150)

# Fig 2: d13CH4 vs d13CO2 crossplot with constant-alpha (eps) isolines
aref <- c(40,55,65,80)  # eps permil -> alpha = 1+eps/1000
iso <- tibble(eps=aref) %>% rowwise() %>%
  mutate(x1=min(wt$d13CO2), x2=max(wt$d13CO2),
         y1=(x1+1000)/(1+eps/1000)-1000, y2=(x2+1000)/(1+eps/1000)-1000) %>% ungroup()
p_cross <- ggplot(wt, aes(d13CO2, d13CH4)) +
  geom_segment(data=iso, aes(x=x1,xend=x2,y=y1,yend=y2,group=eps), inherit.aes=FALSE,
               linetype="dotted", color="grey55") +
  geom_text(data=iso, aes(x=x2,y=y2,label=paste0("eps=",eps)), inherit.aes=FALSE, hjust=-.05, size=3, color="grey40") +
  geom_point(aes(color=ch4_ppm>=10), alpha=.75, size=1.8) +
  scale_color_manual(values=c(`FALSE`="grey70",`TRUE`="black"), labels=c("<10 ppm (atmos/oxid-biased)",">=10 ppm (reliable)"), name=NULL) +
  coord_cartesian(clip="off") + theme_bw(base_size=10) + theme(plot.margin=margin(5,45,5,5), legend.position="bottom") +
  labs(title="d13C-CH4 vs d13C-CO2 with apparent-fractionation (eps_C) isolines",
       subtitle=sprintf("CO2-CH4 coupling r=%.2f (>=10 ppm), p=%.3f -> weak: bulk CO2 mostly respiration", cpl_hi$estimate, cpl_hi$p.value),
       x=expression(delta^13*C-CO[2]~"(permil)"), y=expression(delta^13*C-CH[4]~"(permil)"))
ggsave(file.path(out,"ISOTOPES_fig2_crossplot.png"), p_cross, width=8, height=6, dpi=150)

# ---- 6. METHODS + RESULTS markdown ------------------------------------------
methods <- c(
"## Methods — stable carbon isotopes of internal CH4 and CO2","",
"Internal stem gas d13C-CH4 and d13C-CO2 were measured by cavity ring-down",
"spectroscopy (Picarro G2201-i; precision <1.15 permil at 10 ppm CH4). Analyses",
"used whole-tree single-borehole samples; calibration standards, atmospheric",
"references, incubation vials, and separately-cored heartwood/sapwood samples",
"were excluded, and samples with CH4 < 1.5 ppm were removed (unreliable d13C).",
"",
"The CH4 source signature was estimated two ways: (i) a Keeling model",
"(d13C-CH4 vs 1/[CH4]; intercept = source; 95% CI parametric, so it assumes",
"error-free 1/[CH4] -- the classic Keeling caveat -- but the intercept is stable",
"to dropping the highest-concentration points), and (ii) a two-endmember mixing",
"correction removing an atmospheric background (1.9 ppm, -47 permil) from each",
"sample in isotope-ratio space, applied only where CH4 exceeds background by a",
"margin (operative enrichment floor ~2.4 ppm). A Miller-Tans regression was also",
"computed; because its slope is carried by the few highest-CH4 samples, we treat",
"it as determined by (not robust to) those points and use the atmosphere",
"correction as the tie-breaker between the Keeling and Miller-Tans estimates.",
"",
"The apparent carbon isotope fractionation between co-occurring CO2 and CH4 was",
"computed within each sample as alpha_C = (d13C-CO2 + 1000)/(d13C-CH4 + 1000),",
"eps_C = (alpha_C - 1) x 1000, and summarized across trees (bulk) and from the",
"source-corrected CH4 value. Because bulk internal CO2 integrates respiration and",
"is not the isolated methanogenic substrate, eps_C is interpreted as an apparent,",
"corroborating fractionation; we report the CO2-CH4 isotopic coupling to bound",
"this. Analyses use d13C only (no dD-CH4), which limits definitive pathway",
"discrimination.")

res <- c(
"## Results — isotopes","",
sprintf("Whole-tree samples: n=%d with d13C-CH4 (n=%d also with paired d13C-CO2).", nrow(d), nrow(wt)),"",
"### Composition",
sprintf("- d13C-CH4 (bulk median) = %.1f permil (IQR %s)", bulk_ch4, iqr(d$d13CH4)),
sprintf("- d13C-CO2 (bulk median) = %.1f permil (IQR %s)", bulk_co2, iqr(wt$d13CO2)),"",
"### CH4 source signature",
sprintf("- Keeling intercept = %.0f permil [95%% CI %.0f, %.0f]; robust (drop top 0/3/5 high-CH4 pts: %s)",
        keel_src, keel_ci[1], keel_ci[2], paste(round(kdrop), collapse="/")),
sprintf("- Atmosphere mixing-corrected = %.0f permil (endmember-robust: %s for d_atm -47/-40/-35)",
        atm_src, paste(round(atm_sens), collapse="/")),
sprintf("- Miller-Tans = %.0f permil; its slope is carried by the few highest-CH4 points", mt_src),
sprintf("  (drop top 0/3/5 -> %s), so it is determined by, not robust to, those points.", paste(round(mt_drop), collapse="/")),
"  Keeling and Miller-Tans weight opposite ends of the concentration range; we take",
"  the atmosphere correction (-70) as the independent tie-breaker between them.",
sprintf("- Interpretation: tree-produced source (~%.0f to %.0f permil) is more depleted than the", atm_src, keel_src),
"  bulk measured value because measured gas is diluted with atmospheric CH4.","",
"### Apparent CO2-CH4 fractionation (eps_C)",
sprintf("- bulk within-tree median eps_C = %.0f permil", eps_bulk),
sprintf("- from source-corrected CH4: eps_C = %.0f (atmosphere-corrected) to %.0f (Keeling) permil", eps_src_atm, eps_src_keel),
sprintf("  [= f(source d13CH4, median d13CO2 %.0f); this ASSUMES methanogenic CO2 ~ bulk CO2]", bulk_co2),
sprintf("- CO2-CH4 coupling: r=%.2f overall (p=%.2f); r=%.2f at >=10 ppm (p=%.3f) -> WEAK.",
        cpl_all$estimate, cpl_all$p.value, cpl_hi$estimate, cpl_hi$p.value),
"  Bulk CO2 is largely respiration-dominated, so the methanogenic-CO2 ~ bulk-CO2 assumption",
"  is questionable and eps_C is an apparent (corroborating) indicator with unknown-direction",
"  bias, NOT a measured fractionation.",
sprintf("- Among-species d13C-CH4: no detectable difference at reliable CH4 (Kruskal p=%.2f), but", kw_p),
sprintf("  only %d species-groups / %d samples entered the test (of %d/%d matched to species),",
        kw_ngrp, kw_n, kw_matched, kw_total),
"  i.e. modest power -> consistent with, not proof of, a single pooled source.","",
"### Bottom line",
sprintf("Depleted d13C-CH4 (bulk %.1f; tree-produced ~%.0f to %.0f permil), a within-tree", bulk_ch4, atm_src, keel_src),
sprintf("apparent CO2-CH4 fractionation of ~%.0f-%.0f permil (at/into the CO2-reduction range),", eps_src_atm, eps_src_keel),
"and hydrogenotrophic Methanobacteriaceae dominance are together consistent with",
"substantial hydrogenotrophic methanogenesis. The taxonomy is the primary evidence;",
"eps_C is consistent-with, not diagnostic (respiration-dominated bulk CO2, d13C only,",
"and possible oxidation shifting the apparent source).")

writeLines(methods, file.path(out, "ISOTOPES_methods.md"))
writeLines(res,     file.path(out, "ISOTOPES_results.md"))
cat(paste(res, collapse = "\n"), "\n\nWrote: ISOTOPES_methods.md, ISOTOPES_results.md,",
    "ISOTOPES_fig1_source.png, ISOTOPES_fig2_crossplot.png, ISOTOPES_sample_table.csv\n")
