# ==============================================================================
# REVISION — Synthesis figure: hydrogenotrophic, syntrophy-fed methanogenesis.
# Reproduces existing figures' EXACT data/filters (S6 strict, FAPROTAX/S3, Fig6, S9),
# restyled into 4 panels with a SINGLE semantic colour scheme used consistently:
#   purple = methanogenesis/methanogen/hydrogenotrophic   orange = fermentation/H2 source
#   teal   = C1 metabolism / methylotrophy                amber  = H2 oxidation
#   blue   = aerobic / methanotrophy                      grey   = other
# NEW file. Run: Rscript revision/analysis/R_fig_hydrogenotrophy.R -> fig_hydrogenotrophy.png
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- palette: red/blue diverging (as in the map figure) + orange + grey ------
RED  <- "#b2182b"   # methanogen (a) / C1 (c) / mcrA-associated, positive (b)
BLU  <- "#2166ac"   # aerobic (c) / anti-correlated, negative (b)
SYN  <- "#d4a017"   # syntrophic fermenters (a) / fermentation (c)  (mustard/gold)
OTH  <- "#969696"   # other / neutral
METH <- "#762a83"   # isotope panel (d) methanogen theme (unchanged)
DENS_F <- "#d7c6e0" # isotope density fill (light purple)

th <- theme_bw(base_size = 13) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12), axis.text = element_text(size = 11.5),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.45, "cm"), legend.text = element_text(size = 10.5),
        legend.title = element_text(size = 10.5))

# ---- (a) TAXA (S6 STRICT family-mcrA associations, FDR < 0.05) ----------------
# Three evidence-graded groups (see revision/text/panel_a_taxa_citations.md):
#   Methanogen (red) | Associate / syntroph (gold) | Other fermenter (grey)
methanogen_fams <- c("Methanobacteriaceae","Methanomassiliicoccaceae","Methanosarcinaceae","Methanocellaceae")
syntroph_fams   <- c("Christensenellaceae","Dysgonomonadaceae")  # demonstrated / community-level H2 associates
fam_group <- function(f) ifelse(f %in% methanogen_fams, "Methanogen",
                          ifelse(f %in% syntroph_fams,  "Associate / syntroph", "Other fermenter"))
fa <- read.csv("outputs/tables/family_mcra_associations_strict.csv") %>%
  mutate(t = as.numeric(t)) %>% filter(is.finite(t), FDR < 0.05) %>% arrange(desc(t)) %>%
  mutate(group = fam_group(family),
         group = factor(group, levels = c("Methanogen","Associate / syntroph","Other fermenter")),
         family = fct_reorder(family, t))
pa <- ggplot(fa, aes(t, family, color = group)) +
  geom_segment(aes(x = 0, xend = t, yend = family), color = "grey78", linewidth = 0.6) +
  geom_point(size = 3.6) +
  scale_color_manual(values = c(Methanogen = RED, `Associate / syntroph` = SYN, `Other fermenter` = OTH), name = NULL) +
  labs(x = "Co-occurrence with mcrA (t-statistic; FDR < 0.05)", y = NULL) +
  ggtitle(expression(bold("a  Taxa: methanogens + fermentative associates"))) + th +
  theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
        legend.background = element_rect(fill="white", color="grey80"),
        axis.text.y = element_text(face = "italic", size = 10))

# ---- (b) FUNCTION (FAPROTAX HW/SW) — documented energy/redox-metabolism set ----
# Values reproduced from the FAPROTAX pipeline (revision/outputs/FAPROTAX_all_functions_HW_SW.csv).
# Inclusion rule: energy/redox functions with mean rel. abundance >= 0.1% in a compartment;
# host-association/pathogen/phototrophy and redundant parent annotations (generic chemoheterotrophy,
# generic/CO2/methyl methanogenesis) excluded; hydrogenotrophic_methanogenesis shown as the
# representative methanogenesis route. See panel_b_faprotax_selection.md.
fap_all <- read.csv("revision/outputs/FAPROTAX_all_functions_HW_SW.csv")
fap_sel <- tribble(~func, ~disp,
  "hydrogenotrophic_methanogenesis", "Hydrogenotrophic methanogenesis",
  "dark_hydrogen_oxidation",         "Dark hydrogen oxidation",
  "fermentation",                    "Fermentation",
  "anaerobic_chemoheterotrophy",     "Anaerobic chemoheterotrophy",
  "methylotrophy",                   "Methylotrophy",
  "nitrate_reduction",               "Nitrate reduction",
  "aerobic_chemoheterotrophy",       "Aerobic chemoheterotrophy",
  "methanol_oxidation",              "Methanol oxidation",
  "methanotrophy",                   "Methanotrophy",
  "cellulolysis",                    "Cellulolysis")
fap <- fap_sel %>% left_join(fap_all, by = "func") %>%
  mutate(lr = log2(HW/SW), fn = fct_reorder(disp, lr))
pb <- ggplot(fap, aes(lr, fn, fill = lr)) +
  geom_col() + geom_vline(xintercept = 0, color = "grey50") +
  scale_fill_gradient2(low = BLU, mid = "grey93", high = RED, midpoint = 0, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.05))) +
  labs(x = expression(log[2]*"(heartwood / sapwood)"), y = NULL) + th

# ---- (c) PATHWAY (Fig6 filter gc<0.10): a-priori FUNCTIONAL categories --------
# Not top-t. Classify the full FDR<0.01, low-mcrA-contribution set by MetaCyc
# function (membership independent of the association), drop housekeeping, show
# representative pathways per mechanistic category. See panel_c_pathway_categories.md
ABLU <- "#92c5de"  # light blue: sulfur/nitrogen metabolism (negative, aerobic/oxidant-associated)
p6 <- read.csv("data/processed/molecular/picrust/pathway_associations_mcra_no_mcra_otus.csv")
cb <- read.csv("data/processed/molecular/picrust/pathway_associations_combined.csv")
contrib <- setNames(cb$mean_percent_from_mcra, cb$pathway)
sig <- p6 %>% filter(!is.na(FDR), FDR < 0.05) %>%
  mutate(gc = ifelse(is.na(contrib[pathway]), 0, contrib[pathway])) %>% filter(gc < 0.10)
classify_pw <- function(id, desc) { x <- tolower(paste(id, desc))
  if (grepl("acetyl coenzyme a|calvin|rump|formaldehyde|carbon fix", x)) return("C1 / carbon fixation")
  if (grepl("ferment|glycolysis|glycogen|mannan|starch|xylan|cellulose|propanoate|butanoate|glutamate degrad|lactate|pyruvate", x)) return("Fermentation / carbohydrate")
  if (grepl("aerobic respiration|tca cycle|glyoxylate|ubiquin|cytochrome|menaquin", x)) return("Aerobic respiration")
  if (grepl("sulfate|sulfur|thiosulf|sulfhydr|cysteine|nitrate|nitrite|denitrif|nitrogen", x)) return("Sulfur / nitrogen metabolism")
  return("other") }  # lipid/fatty-acid biosynthesis -> 'other' (in SI table, not displayed)
DISP <- c("C1 / carbon fixation","Fermentation / carbohydrate","Aerobic respiration","Sulfur / nitrogen metabolism")
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
  geom_segment(aes(x = 0, xend = t, yend = lab), color = "grey78", linewidth = 0.6) +
  geom_point(size = 3.2) + geom_vline(xintercept = 0, color = "grey50") +
  scale_color_manual(values = c(`C1 / carbon fixation`=RED, `Fermentation / carbohydrate`=SYN,
                                `Aerobic respiration`=BLU, `Sulfur / nitrogen metabolism`=ABLU), name = NULL) +
  labs(x = "Association with mcrA (t-statistic)", y = NULL,
       title = "c  Pathway (PICRUSt2): anaerobic C-metabolism tracks mcrA") + th +
  theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
        legend.background = element_rect(fill="white", color="grey80"),
        legend.text = element_text(size = 8), axis.text.y = element_text(size = 8.3, lineheight = 0.85))

# ---- (d) ISOTOPES (Fig-S9 raincloud) + Keeling source; points sized by CH4 ----
XLO <- -118; XHI <- 28
iso <- read.csv("revision/outputs/ISOTOPES_sample_table.csv") %>%
  filter(is.finite(d13CH4), d13CH4 > -115, d13CH4 < 25, ch4_ppm >= 1.5)
dens <- density(iso$d13CH4); dd <- data.frame(x = dens$x, y = dens$y/max(dens$y)*0.55) %>% filter(x >= XLO, x <= XHI)
set.seed(42); iso$jy <- runif(nrow(iso), -0.30, -0.06); iso$lc <- log10(iso$ch4_ppm)
keel <- -79; keel_lo <- -92; keel_hi <- -66; by <- -0.52; bs <- 0.15; bc <- 0.035
brk <- function(x1,x2,y,lab) list(
  annotate("segment", x=x1, xend=x2, y=y, yend=y, color="gray30", linewidth=0.9),
  annotate("segment", x=x1, xend=x1, y=y, yend=y+bc, color="gray30", linewidth=0.5),
  annotate("segment", x=x2, xend=x2, y=y, yend=y+bc, color="gray30", linewidth=0.5),
  annotate("text", x=(x1+x2)/2, y=y-0.05, label=lab, color="gray30", size=3, fontface="italic"))
pd <- ggplot() +
  geom_vline(xintercept=-47, linetype="dashed", color="gray30", linewidth=0.5) +
  annotate("text", x=-47, y=0.55, label="Atmosphere", color="gray30", size=3, fontface="italic", hjust=-0.08) +
  geom_ribbon(data=dd, aes(x=x, ymin=0, ymax=y), fill=DENS_F, color=METH, alpha=0.7, linewidth=0.7) +
  geom_point(data=iso, aes(x=d13CH4, y=jy, size=lc), color=METH, fill=DENS_F, shape=21, alpha=0.6, stroke=0.3) +
  scale_size_continuous(range=c(0.6, 3.4), name=expression("CH"[4]*" conc. (ppm)"),
                        breaks=c(1,2,3,4), labels=c("10","100","1000","10k"),
                        guide=guide_legend(override.aes=list(alpha=0.85))) +
  annotate("segment", x=keel_lo, xend=keel_hi, y=0.72, yend=0.72, color=METH, linewidth=0.9) +
  annotate("segment", x=c(keel_lo,keel_hi), xend=c(keel_lo,keel_hi), y=0.69, yend=0.75, color=METH, linewidth=0.6) +
  annotate("point", x=keel, y=0.72, color=METH, size=3) +
  annotate("text", x=keel, y=0.84, label="Keeling source (95% CI)", color=METH, size=3.2, fontface="bold", hjust=0.5) +
  brk(-110,-60, by, "Hydrogenotrophic") + brk(-65,-50, by-bs, "Acetoclastic") + brk(-70,-50, by-2*bs, "Methylotrophic") +
  scale_x_continuous(breaks=seq(-120,20,20)) + coord_cartesian(xlim=c(XLO,XHI), ylim=c(-0.92, 0.92), clip="off") +
  labs(x=expression(delta^13*"C-CH"[4]*" (per mil VPDB)"), y=NULL,
       title=expression(bold("d  Isotopes: "*delta^13*"C-CH"[4]*" depleted; source hydrogenotrophic"))) +
  th + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
             legend.position=c(0.99,0.05), legend.justification=c(1,0),
             legend.background=element_rect(fill="white",color="grey80"),
             plot.margin=margin(6,10,42,6), panel.grid=element_blank())

fig <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & theme(plot.tag = element_text(face = "bold", size = 18))
ggsave(file.path(out, "fig_hydrogenotrophy.png"), fig, width = 13.5, height = 11, dpi = 300, bg = "white")
cat("Wrote fig_hydrogenotrophy.png (semantic palette; base_size 11)\n")
cat("(a) families:", paste(as.character(fa$family), collapse=", "), "\n")
cat("(c) categories:", paste(levels(droplevels(pwc$cat)), collapse=" | "), "\n")
