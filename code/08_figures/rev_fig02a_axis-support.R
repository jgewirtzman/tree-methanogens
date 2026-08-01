# ==============================================================================
# REVISION — Fig 2A with a CONSISTENT x-axis (Referee 2 #6)
# ==============================================================================
# R2: Fig 2A x-axis is "non-intuitive... neither a standard linear nor log scale,
# and it is inconsistent across panels (e.g. Quercus rubra looks linear, Carya
# ovata does not); outliers (Acer rubrum, Quercus rubra) distort the panels; and
# the dashed line is undefined."
#
# DIAGNOSIS (from code/02_flux/static/04_height_effect_analysis.R):
#   The original panel (a) drew each species on its OWN free LINEAR x-axis whose
#   three ticks were literally that species' (min, midpoint, max). So (i) every
#   panel had a different scale, (ii) a single outlier maximum stretched the whole
#   panel and squashed the rest near zero, and (iii) the ticks looked neither
#   linear-standard nor log. The dashed line is the zero-net-flux reference.
#
# FIX: one SHARED axis for all panels. Because flux spans ~3 orders of magnitude
#   (-0.028 to 6.2 nmol m-2 s-1) AND includes small negatives (net uptake), we use
#   a single signed pseudo-log (symlog / arcsinh-like) transform -- the same
#   family of transform already adopted for the skewed regressions in this revision
#   (scales::pseudo_log_trans). Identical breaks on every panel => directly
#   comparable across species, outliers no longer distort, negatives/zeros handled.
#   Dashed line is drawn and DEFINED (zero net flux).
#
# NEW file (does not edit the pipeline). Reuses the same processed inputs as
#   04_height_effect_analysis.R.
# Run: Rscript code/revision/R_fig2a_consistent_axis.R
# Output: outputs/revision/fig2a_consistent_axis.png (drop-in replacement panel a)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales); library(gghalves)
source("code/lib/outputs.R")
})
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

species_mapping <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",FAGR="Fagus grandifolia",FRAM="Fraxinus americana",
  PIST="Pinus strobus",QURU="Quercus rubra",TSCA="Tsuga canadensis",CAOV="Carya ovata",
  KALA="Kalmia latifolia",PRSE="Prunus serotina",QUAL="Quercus alba",QUVE="Quercus velutina",
  SAAL="Sassafras albidum")

# ---- load exactly as the pipeline does ---------------------------------------
aux <- read.csv("data/processed/flux/goflux_auxfile.csv")
ch4 <- read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv")
d <- merge(ch4, aux[, c("UniqueID","measurement_height","tree_id","species","plot")],
           by = "UniqueID", all.x = TRUE)
names(d)[names(d) == "best.flux"] <- "CH4_best.flux"

d <- d %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux), !is.na(height_numeric)) %>%
  mutate(species = ifelse(species == "TSLA", "TSCA", species),
         height_m = height_numeric / 100) %>%
  filter(species != "KALA") %>%
  mutate(tree_unique = paste(plot, tree_id, sep = "_"),
         species_latin = species_mapping[species])

counts <- d %>% group_by(species, species_latin) %>%
  summarise(n_trees = n_distinct(tree_unique), .groups = "drop") %>%
  mutate(species_label = paste0(species_latin, "\n(n=", n_trees, ")"))
d <- d %>% left_join(counts[, c("species","species_label")], by = "species") %>%
  mutate(species_label = factor(species_label, levels = sort(unique(species_label))))

# ---- shared signed pseudo-log axis -------------------------------------------
# sigma sets the near-zero linear window; below ~sigma the scale is ~linear
# (so small values and negatives are readable), above it ~log10 (tames outliers).
sig <- 0.005
trans_symlog <- pseudo_log_trans(sigma = sig, base = 10)
brks   <- c(-0.01, 0, 0.01, 0.1, 1)
lims   <- range(c(d$CH4_best.flux, brks), na.rm = TRUE)
lab_fun <- function(x) ifelse(x == 0, "0",
                       ifelse(abs(x) < 0.1, sprintf("%.2f", x), sprintf("%.1f", x)))

boxplot_color <- "#756BB1"

p <- ggplot(d, aes(x = factor(height_m), y = CH4_best.flux)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40",
             linewidth = 0.5, alpha = 0.9) +
  geom_half_boxplot(alpha = 0.6, outlier.shape = NA, fill = boxplot_color,
                    color = "gray30", side = "r") +
  geom_jitter(alpha = 0.4, size = 0.6, color = "gray20",
              position = position_jitter(width = 0.15, height = 0)) +
  coord_flip() +
  scale_y_continuous(trans = trans_symlog, breaks = brks, limits = lims,
                     labels = lab_fun, expand = expansion(mult = c(0.03, 0.05))) +
  facet_wrap(~species_label, ncol = 5) +
  labs(x = "Height (m)",
       y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*"),  signed pseudo-log scale")) +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(size = 7, face = "italic"),
        axis.text = element_text(size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
        panel.border = element_rect(color = "gray40", fill = NA, linewidth = 0.4),
        panel.spacing = unit(0.4, "lines"))

ggsave(out_path("fig2a_consistent_axis.png"), p,
       width = 10, height = 6.2, dpi = 300, bg = "white")

cat("Wrote", out_path("fig2a_consistent_axis.png"), "\n")
cat("Shared signed pseudo-log axis (sigma =", sig, "), identical breaks on every panel:\n  ",
    paste(brks, collapse = ", "), "\n")
cat("Dashed line = zero net flux (left = net uptake, right = net emission).\n")
cat("Flux range:", round(lims[1], 3), "to", round(lims[2], 2), "nmol m-2 s-1.\n")
