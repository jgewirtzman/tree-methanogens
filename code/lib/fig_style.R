# ==============================================================================
# fig_style.R -- shared style for the Supporting Information figures. SOURCED.
# ------------------------------------------------------------------------------
# One export width, one base font size, one tag style and one set of colours, taken
# from the main figures so the SI matches them:
#   stem / soil          Fig. 1 (fig01_temporal-flux.R)
#   emission / uptake    Fig. 9 (fig09_budget.R)
#   compartments         Fig. 4 (fig04_gene-abundance.R)
# Figures are saved SI_WIDTH inches wide at 300 dpi, so placed at page width
# (about 6.5 in) text prints at about 0.8 x SI_BASE points, as in the main figures.
# ==============================================================================
suppressPackageStartupMessages(library(ggplot2))

SI_WIDTH <- 8          # inches; every SI figure
SI_DPI   <- 300
SI_BASE  <- 11         # base font size (pt) at SI_WIDTH

COL_STEM   <- "#228B22"
COL_SOIL   <- "#8B4513"
COL_EMIT   <- "#b2182b"
COL_UPTAKE <- "#2166ac"
COL_NS     <- "grey60"   # not detected / not significant

COMP_COLS <- c(Heartwood = "#a6611a", Sapwood = "#dfc27d",
               `Mineral soil` = "#80cdc1", `Organic soil` = "#018571")
COMP_LEVELS <- names(COMP_COLS)

# broadleaf / conifer (Fig. S6 and wherever species are grouped)
GROUP_COLS <- c(Broadleaf = "#1b7837", Conifer = "#762a83")
CONIFERS <- c("Pinus strobus", "Tsuga canadensis")
species_group <- function(species) ifelse(species %in% CONIFERS, "Conifer", "Broadleaf")

theme_si <- function(base_size = SI_BASE) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = "grey60"),
          strip.text = element_text(size = rel(0.9)),
          legend.key.size = unit(0.9, "lines"),
          plot.title = element_text(size = rel(1), face = "plain"),
          plot.tag = element_text(size = SI_BASE + 1, face = "bold"))
}

# panel tags "(a)", "(b)", ... for patchwork assemblies:  p + si_tags()
# (the tag style is in theme_si(); a plot_annotation() cannot carry an `& theme()`)
si_tags <- function() {
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
}

# save at the SI width; height in inches
si_save <- function(filename, plot, height, width = SI_WIDTH, ...) {
  ggsave(filename, plot, width = width, height = height, dpi = SI_DPI, bg = "white", ...)
}

# axis label helpers (plotmath)
lab_flux  <- expression(CH[4] ~ flux ~ (nmol ~ m^{-2} ~ s^{-1}))
lab_log10 <- function(gene) bquote(log[10] ~ italic(.(gene)))
