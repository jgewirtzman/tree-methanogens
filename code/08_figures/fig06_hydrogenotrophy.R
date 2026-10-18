# ==============================================================================
# fig06_hydrogenotrophy.R -- Figure 6: hydrogenotrophic, syntrophy-fed methanogenesis
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root, after figS14_taxonomy-mcra.R.
#   (a) taxa: 16S families associated with mcrA (p < 0.01; filled if FDR < 0.05; the table of Figure S10),
#       grouped as methanogen / associate or syntroph / other fermenter
#   (b) function: FAPROTAX energy and redox functions, log2 heartwood / sapwood; asterisks,
#       paired Wilcoxon (heartwood vs sapwood of the same tree), BH FDR
#   (c) pathways: PICRUSt2 MetaCyc pathways associated with mcrA (FDR < 0.05, < 10% of
#       the pathway from mcrA-carrying ASVs), by functional category
#   (d) isotopes: d13C-CH4 of whole-stem gas (density and samples sized by CH4) with
#       the Keeling source and the ranges of the three methanogenic pathways
# Colours: red = methanogen / mcrA-associated / C1; blue = aerobic / anti-correlated;
# gold = syntrophic fermenters / fermentation; grey = other; purple = isotopes.
# Inputs
#   outputs/tables/family_mcra_associations_strict.csv       figS14_taxonomy-mcra.R
#   outputs/data/FAPROTAX_HW_SW_tests.csv                    07_molecular/01_faprotax_dump_HW_SW.R
#   outputs/data/picrust/pathway_associations_mcra_no_mcra_otus.csv,
#     pathway_associations_combined.csv                      07_molecular/02_picrust_pathway_associations.R
#   outputs/data/ISOTOPES_sample_table.csv, ISOTOPES_summary.csv  09_tables_stats/12_isotopes-canonical.R
# Output: outputs/figures/generated/fig_hydrogenotrophy.png
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

# ---- palette: red/blue diverging (as in the map figure) + orange + grey ------
mlab <- function(b) sub("^-", "\u2212", format(b, trim = TRUE, drop0trailing = TRUE))   # true minus on axes
RED  <- "#b2182b"   # methanogen (a) / C1 (c) / mcrA-associated, positive (b)
BLU  <- "#2166ac"   # aerobic (c) / anti-correlated, negative (b)
SYN  <- "#d4a017"   # syntrophic fermenters (a) / fermentation (c)  (mustard/gold)
OTH  <- "#969696"   # other / neutral
METH <- "#762a83"   # isotope panel (d) methanogen theme (unchanged)
DENS_F <- "#d7c6e0" # isotope density fill (light purple)

th <- theme_bw(base_size = 13) +
  theme(
        axis.title = element_text(size = 12), axis.text = element_text(size = 11.5),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.45, "cm"), legend.text = element_text(size = 10.5),
        legend.title = element_text(size = 10.5))

# ---- (a) TAXA: family-mcrA associations (Figure S10), p < 0.01; filled = FDR < 0.05 ----
# Three evidence-graded groups:
#   Methanogen (red) | Associate / syntroph (gold) | Other (grey)
source("code/lib/methanogen_families.R"); methanogen_fams <- METHANOGEN_FAMILIES
syntroph_fams   <- c("Christensenellaceae","Dysgonomonadaceae")  # demonstrated / community-level H2 associates
fam_group <- function(f) ifelse(f %in% methanogen_fams, "Methanogen",
                          ifelse(f %in% syntroph_fams,  "Associate / syntroph", "Other"))
fa <- read.csv("outputs/tables/family_mcra_associations_strict.csv") %>%
  mutate(t = as.numeric(t)) %>% filter(is.finite(t), p < 0.01) %>% arrange(desc(t)) %>%
  mutate(group = fam_group(family),
         group = factor(group, levels = c("Methanogen","Associate / syntroph","Other")),
         fdr = factor(ifelse(FDR < 0.05, "FDR < 0.05", "p < 0.01 only"), levels = c("FDR < 0.05", "p < 0.01 only")),
         family = fct_reorder(family, t))
XLAB_AC <- expression("Association with "*italic(mcrA)*" (t-statistic)")   # panels (a) and (c)
pa <- ggplot(fa, aes(t, family, color = group)) +
  geom_segment(aes(x = 0, xend = t, yend = family), color = "grey78", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(aes(shape = fdr), size = 3.6, stroke = 1.1, fill = "white") +
  scale_shape_manual(values = c(`FDR < 0.05` = 16, `p < 0.01 only` = 21), name = NULL) +
  scale_color_manual(values = c(Methanogen = RED, `Associate / syntroph` = SYN, Other = OTH), name = NULL) +
  guides(color = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1)) +
  scale_x_continuous(labels = mlab) +
  labs(x = XLAB_AC, y = NULL) +
  th +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 9.5), legend.key.size = unit(0.35, "cm"),
        legend.background = element_blank(),
        axis.text.y = element_text(face = "italic", size = 10))

# ---- (b) FUNCTION (FAPROTAX, heartwood / sapwood) ------------------------------
# The functions are listed once in code/lib/fig6b_functions.R (inclusion rule there).
# Methanogen ASVs are removed from every function except the methanogenesis row:
# FAPROTAX files Methanobacteriaceae under dark H2 oxidation and Methanomassiliicoccaceae
# under methylotrophy, which would recount the methanogens as independent functions.
# Bars: log2 of the unrounded pooled shares; asterisks: paired Wilcoxon per tree, BH FDR
# (the test is described in the caption, not on the panel)
# (01_faprotax_dump_HW_SW.R). Non-significant bars carry no mark.
source("code/lib/fig6b_functions.R")
fap <- FIG6B_FUNCTIONS %>% left_join(read.csv("outputs/data/FAPROTAX_HW_SW_tests.csv"), by = "func") %>%
  mutate(lr = log2_ratio, fn = fct_reorder(label, lr),
         stars = case_when(paired_FDR < 0.001 ~ "***", paired_FDR < 0.01 ~ "**", paired_FDR < 0.05 ~ "*", TRUE ~ ""))
stopifnot(all(is.finite(fap$lr)), !anyNA(fap$paired_FDR))
pb <- ggplot(fap, aes(lr, fn, fill = lr)) +
  geom_col() + geom_vline(xintercept = 0, color = "grey50") +
  geom_text(aes(x = lr + sign(lr) * 0.12, label = stars, hjust = ifelse(lr > 0, 0, 1)),
            vjust = 0.75, size = 4.6, color = "grey20") +
  scale_fill_gradient2(low = BLU, mid = "grey93", high = RED, midpoint = 0, guide = "none") +
  scale_x_continuous(labels = mlab, expand = expansion(mult = c(0.14, 0.16))) +
  labs(x = expression(log[2]*"(heartwood / sapwood)"), y = NULL) + th

# ---- (c) PATHWAY: a-priori FUNCTIONAL categories --------------------------------
# Not top-t. The pathways associated with mcrA (FDR < 0.05) with < 10% of their
# abundance from mcrA-carrying ASVs are classified by MetaCyc function (membership
# independent of the association); housekeeping is dropped and the four pathways with
# the largest |t| are shown per category.
ABLU <- "#92c5de"  # light blue: sulfur/nitrogen metabolism (negative, aerobic/oxidant-associated)
p6 <- read.csv("outputs/data/picrust/pathway_associations_mcra_no_mcra_otus.csv")
cb <- read.csv("outputs/data/picrust/pathway_associations_combined.csv")
contrib <- setNames(cb$mean_percent_from_mcra, cb$pathway)
sig <- p6 %>% filter(!is.na(FDR), FDR < 0.05) %>%
  mutate(gc = ifelse(is.na(contrib[pathway]), 0, contrib[pathway])) %>% filter(gc < 0.10)
classify_pw <- function(id, desc) { x <- tolower(paste(id, desc))
  # A pathway goes in a category only when its MetaCyc function places it there
  # unambiguously: TCA variants of anaerobes and autotrophs (Fd-dependent, reductive, acetate-
  # producer, Helicobacter) and menaquinone (an anaerobic-respiration quinone) are not counted
  # as aerobic respiration; carbohydrate breakdown is labelled as such, not as fermentation.
  if (grepl("tca cycle (v|vi|vii|viii) |reductive tca|menaquin", x)) return("other")
  if (grepl("acetyl coenzyme a|calvin|rump|formaldehyde|carbon fix", x)) return("C1 / carbon fixation")
  if (grepl("ferment|glycolysis|glycogen|mannan|starch|xylan|cellulose|propanoate|butanoate|glutamate degrad|lactate|pyruvate", x)) return("Fermentation and carbohydrate breakdown")
  if (grepl("aerobic respiration|tca cycle|glyoxylate cycle|glyoxylate bypass|ubiquin|cytochrome", x)) return("TCA cycle and aerobic respiration")
  if (grepl("sulfate|sulfur|thiosulf|sulfhydr|cysteine|nitrate|nitrite|denitrif|nitrogen", x)) return("Sulfur / nitrogen metabolism")
  return("other") }  # lipid/fatty-acid biosynthesis -> 'other' (in SI table, not displayed)
DISP <- c("C1 / carbon fixation","Fermentation and carbohydrate breakdown","TCA cycle and aerobic respiration","Sulfur / nitrogen metabolism")
pwc <- sig %>% rowwise() %>% mutate(cat = classify_pw(pathway, description)) %>% ungroup() %>%
  filter(cat %in% DISP, !grepl("anaerobic", tolower(description))) %>%
  group_by(cat) %>% arrange(desc(abs(t))) %>% slice(1:4) %>% ungroup() %>%
  mutate(cat = factor(cat, levels = DISP),
         lab = description %>%
           str_replace_all("superpathway of |superpathay of ", "") %>%
           str_replace_all("biosynthesis", "biosynth.") %>%
           str_replace_all("assimilation", "assim.") %>%
           str_replace_all("\\(2-oxoglutarate:ferredoxin oxidoreductase\\)", "(Fd oxidoreductase)") %>%
           str_replace_all("\\(2-oxoglutarate decarboxylase\\)", "(2-OG decarboxylase)") %>%
           str_replace_all("\\(by sulfhydrylation\\)", "(sulfhydrylation)") %>%
           str_wrap(24),
         lab = fct_reorder(lab, t))
pc <- ggplot(pwc, aes(t, lab, color = cat)) +
  guides(color = guide_legend(ncol = 2)) +
  scale_x_continuous(labels = mlab) +
  geom_segment(aes(x = 0, xend = t, yend = lab), color = "grey78", linewidth = 0.6) +
  geom_point(size = 3.2) + geom_vline(xintercept = 0, color = "grey50") +
  scale_color_manual(values = c(`C1 / carbon fixation`=RED, `Fermentation and carbohydrate breakdown`=SYN,
                                `TCA cycle and aerobic respiration`=BLU, `Sulfur / nitrogen metabolism`=ABLU), name = NULL) +
  labs(x = XLAB_AC, y = NULL) + th +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(),
        legend.text = element_text(size = 8), axis.text.y = element_text(size = 8.3, lineheight = 0.85))

# ---- (d) ISOTOPES: d13C-CH4 raincloud + Keeling source; points sized by CH4 ----
XLO <- -118; XHI <- 28
iso <- read.csv("outputs/data/ISOTOPES_sample_table.csv")   # already screened (code/lib/isotope_samples.R)
dens <- density(iso$d13CH4); dd <- data.frame(x = dens$x, y = dens$y/max(dens$y)*0.55) %>% filter(x >= XLO, x <= XHI)
set.seed(42); iso$jy <- runif(nrow(iso), -0.30, -0.06); iso$lc <- log10(iso$ch4_ppm)
.iso_sum <- with(read.csv("outputs/data/ISOTOPES_summary.csv"), setNames(value, quantity))   # 12_isotopes-canonical.R
keel <- .iso_sum[["keeling_source"]]; keel_lo <- .iso_sum[["keeling_ci_lo"]]; keel_hi <- .iso_sum[["keeling_ci_hi"]]; by <- -0.52; bs <- 0.15; bc <- 0.035
brk <- function(x1,x2,y,lab) list(
  annotate("segment", x=x1, xend=x2, y=y, yend=y, color="gray30", linewidth=0.9),
  annotate("segment", x=x1, xend=x1, y=y, yend=y+bc, color="gray30", linewidth=0.5),
  annotate("segment", x=x2, xend=x2, y=y, yend=y+bc, color="gray30", linewidth=0.5),
  annotate("text", x=(x1+x2)/2, y=y-0.05, label=lab, color="gray30", size=3, fontface="italic"))
pd <- ggplot() +
  geom_vline(xintercept=-47, linetype="dashed", color="gray30", linewidth=0.5) +
  annotate("text", x=-47, y=0.55, label="Atmosphere", color="gray30", size=3, fontface="italic", hjust=-0.08) +
  geom_ribbon(data=dd, aes(x=x, ymin=0, ymax=y), fill=DENS_F, color=METH, alpha=0.7, linewidth=0.7) +
  # extremes beyond the axis are drawn off-panel otherwise; they are in all statistics (legend)
  geom_point(data=subset(iso, d13CH4 >= XLO & d13CH4 <= XHI), aes(x=d13CH4, y=jy, size=lc), color=METH, fill=DENS_F, shape=21, alpha=0.6, stroke=0.3) +
  scale_size_continuous(range=c(0.6, 3.4), name=expression("CH"[4]*" conc. (ppm)"),
                        breaks=c(1,2,3,4), labels=c("10","100","1000","10k"),
                        guide=guide_legend(override.aes=list(alpha=0.85))) +
  annotate("segment", x=keel_lo, xend=keel_hi, y=0.72, yend=0.72, color=METH, linewidth=0.9) +
  annotate("segment", x=c(keel_lo,keel_hi), xend=c(keel_lo,keel_hi), y=0.69, yend=0.75, color=METH, linewidth=0.6) +
  annotate("point", x=keel, y=0.72, color=METH, size=3) +
  annotate("text", x=keel, y=0.84, label="Keeling source (95% CI)", color=METH, size=3.2, fontface="bold", hjust=0.5) +
  brk(-110,-60, by, "Hydrogenotrophic") + brk(-65,-50, by-bs, "Acetoclastic") + brk(-70,-50, by-2*bs, "Methylotrophic") +
  scale_x_continuous(breaks=seq(-120,20,20), labels = mlab) + coord_cartesian(xlim=c(XLO,XHI), ylim=c(-0.92, 0.92), clip="off") +
  labs(x=expression(delta^13*"C-CH"[4]*" (‰ VPDB)"), y=NULL) +
  th + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
             legend.position="bottom", legend.margin=margin(0,0,0,0),
             legend.background=element_blank(),
             plot.margin=margin(6,10,42,6), panel.grid=element_blank())

fig <- wrap_plots(pa, pb, pc, free(pd), ncol = 2) + plot_layout(widths = c(1, 1), heights = c(0.8, 1.45)) +   # equal panel widths
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & theme(plot.tag = element_text(face = "bold", size = 14))
ggsave(out_path("fig_hydrogenotrophy.png"), fig, width = 10.5, height = 9.4, dpi = 300, bg = "white")   # near print size
cat("Wrote fig_hydrogenotrophy.png\n")
cat("(a) families:", paste(as.character(fa$family), collapse=", "), "\n")
cat("(c) categories:", paste(levels(droplevels(pwc$cat)), collapse=" | "), "\n")
