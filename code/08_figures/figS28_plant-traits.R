# ==============================================================================
# figS28_plant-traits.R -- Figure S25: plant traits and the stem methane-cycling community
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) Spearman correlations of the rule-selected traits with mcrA, pmoA + mmoX, the
#       methanogen:methanotroph ratio and net stem flux, across the gene-flux species
#       (left) and the broadleaf species alone (right); * p < 0.05
#   (b) the three associations named in the Results (trait-response pairs typed in
#       `sel` below), conifers as open symbols; dashed line in every panel, linear fit to
#       the eight broadleaf species (it spans their range only)
# Correlations and the trait inclusion rules come from
# code/09_tables_stats/25_plant-traits.R; descriptive, n = 10 species.
# Inputs: outputs/data/traits_correlations.csv, traits_included.csv (25_plant-traits.R);
#   outputs/data/tree_species_moisture_niche.csv (09_tables_stats/09_tree-distribution.R);
#   data/tables/species_traits.csv (code/tables/07_species_traits.R), species.csv
#   (gymnosperm); species gene medians from code/lib/species_gene_data.R
# Output: outputs/figures/generated/plant_traits.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(ggrepel) })
source("code/lib/outputs.R")
source("code/lib/tables_read.R")
source("code/lib/fig_style.R")   # SI house style: theme_si(), si_tags(), si_save()

S <- read_csv("outputs/data/traits_correlations.csv", show_col_types = FALSE)
inc <- read_csv("outputs/data/traits_included.csv", show_col_types = FALSE) %>% filter(included)
invisible(capture.output(suppressMessages(source("code/lib/species_gene_data.R"))))
resp <- analysis_ratio %>% transmute(species_id, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id")
SP <- read_table("species")
tr <- read_csv("data/tables/species_traits.csv", show_col_types = FALSE) %>% rename(spcode = species_code) %>%
  mutate(gymnosperm = as.numeric(SP$gymnosperm[match(spcode, SP$species_code)]))
niche <- read_csv("outputs/data/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(spcode = names(SPECIES_NAMES)[match(species, SPECIES_NAMES)]) %>% select(spcode, vwc_realized = vwc_mean)
dat <- resp %>% left_join(tr, by = c("species_id" = "spcode")) %>% left_join(niche, by = c("species_id" = "spcode"))

# ---- (a) twin heatmaps --------------------------------------------------------------
# column headers in plotmath: italic gene name, subscript in CH4
RESP <- c("mcra", "meth", "balance", "flux")
rlab <- expression(italic(mcrA), atop(Methano*"-", trophs), Ratio, atop(Stem~CH[4], flux))
gord <- unique(inc$group)
d <- S %>% mutate(group = factor(group, levels = gord), trait = factor(trait, levels = rev(inc$trait)),
                  response = factor(response, levels = RESP))
hm <- function(est, p, title, strip = TRUE) {
  dd <- d %>% mutate(e = {{ est }}, lab = paste0(sprintf("%.2f", e), if_else({{ p }} < 0.05, "*", "")))
  ggplot(dd, aes(response, trait, fill = e)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 2.7) +
    scale_fill_gradient2(low = COL_UPTAKE, mid = "white", high = COL_EMIT, limits = c(-1, 1), name = "Spearman ρ") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(group = label_wrap_gen(13))) +
    scale_x_discrete(position = "top", labels = setNames(as.list(rlab), RESP)) +
    labs(x = NULL, y = NULL, title = title) + theme_si() +
    theme(strip.text.y.left = if (strip) element_text(angle = 0, hjust = 1, face = "bold", size = rel(0.85)) else element_blank(),
          strip.background = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(),
          strip.placement = "outside", axis.text.x.top = element_text(colour = "black", vjust = 0),
          axis.text.y = element_text(colour = "black"), panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, margin = margin(b = 4))) }
h1 <- hm(rho, p, "All ten species") + theme(legend.position = "none")
h2 <- hm(rho_bl, p_bl, "Eight broadleaf species", strip = FALSE) + theme(axis.text.y = element_blank())

# ---- (b) the associations named in the text ----------------------------------------------
sel <- tribble(~trait, ~col, ~response,
  "Bark density",        "bark_density_gcm3",   "balance",
  "Soil-moisture niche", "vwc_realized",        "mcra",
  "Longevity",           "try_plant_longevity", "meth")
xlab <- list(`Bark density` = expression(Bark~density~(g~cm^{-3})),
             `Soil-moisture niche` = "Soil-moisture niche, VWC (%)",
             `Longevity` = "Longevity (years)")
ylab <- list(balance = expression(Methanogen*":"*methanotroph~ratio~(log[10])),
             mcra = expression(italic(mcrA)~(log[10]~copies~g^{-1})),
             meth = expression(Methanotrophs~(log[10]~copies~g^{-1})))
sp <- pmap(sel, function(trait, col, response) {
  s <- S %>% filter(trait == !!trait, response == !!response)
  dd <- dat %>% transmute(x = .data[[col]], y = .data[[response]], grp = if_else(gymnosperm == 1, "Conifer", "Broadleaf"),
                          lab = unname(SPECIES_ABBREV[species_id])) %>% filter(is.finite(x))
  ggplot(dd, aes(x, y)) +
    geom_smooth(data = filter(dd, grp == "Broadleaf"), method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(shape = grp), size = 2, fill = "white", stroke = 0.7) +
    geom_text_repel(aes(label = lab), size = 2.6, fontface = "italic", seed = 1, box.padding = 0.4,
                    point.padding = 0.15, force = 3, max.overlaps = Inf, min.segment.length = 0.3,
                    segment.size = 0.2, segment.colour = "grey55") +
    scale_shape_manual(values = c(Broadleaf = 16, Conifer = 21), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.14, 0.10))) +   # room for labels at both edges
    scale_y_continuous(expand = expansion(mult = c(0.16, 0.12))) +
    labs(x = xlab[[trait]], y = ylab[[response]], subtitle = sprintf("All species: ρ = %.2f (p = %.3f)\nBroadleaf: ρ = %.2f (p = %.3f)",
                                                           s$rho, s$p, s$rho_bl, s$p_bl)) +
    theme_si() + theme(plot.subtitle = element_text(size = rel(0.8), lineheight = 1.05),
                       axis.title = element_text(size = rel(0.9))) })

top <- (h1 | h2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.key.width = unit(1.6, "cm"), legend.key.height = unit(0.3, "cm"),
        legend.title = element_text(vjust = 0.85))
bot <- wrap_plots(sp, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# the wrapped halves do not carry theme_si(), so the tag style is set here as well
fig <- wrap_elements(top) / wrap_elements(bot) + plot_layout(heights = c(1.5, 1)) + si_tags() &
  theme(plot.tag = element_text(size = SI_BASE + 1, face = "bold"))
si_save(out_path("plant_traits.png"), fig, height = 9.8)
cat("Wrote plant_traits.png\n")
