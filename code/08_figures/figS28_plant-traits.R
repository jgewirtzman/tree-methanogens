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
rlab <- c(mcra = "mcrA", meth = "Methano-\ntrophs", balance = "Ratio", flux = "Stem CH\u2084\nflux")
gord <- unique(inc$group)
d <- S %>% mutate(group = factor(group, levels = gord), trait = factor(trait, levels = rev(inc$trait)),
                  response = factor(rlab[response], levels = rlab))
hm <- function(est, p, title, strip = TRUE) {
  dd <- d %>% mutate(e = {{ est }}, lab = paste0(sprintf("%.2f", e), if_else({{ p }} < 0.05, "*", "")))
  ggplot(dd, aes(response, trait, fill = e)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 2.5) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "Spearman ρ") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") + scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = title) + theme_minimal(base_size = 9) +
    theme(strip.text.y.left = if (strip) element_text(angle = 0, hjust = 1, face = "bold") else element_blank(),
          strip.placement = "outside", axis.text.x.top = element_text(size = 8.5, lineheight = 0.9, face = c("italic", "plain", "plain", "plain")), panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 4))) }
h1 <- hm(rho, p, "All ten species") + theme(legend.position = "none")
h2 <- hm(rho_bl, p_bl, "Eight broadleaf species", strip = FALSE) + theme(axis.text.y = element_blank())

# ---- (b) the associations named in the text ----------------------------------------------
sel <- tribble(~trait, ~col, ~response, ~xlab,
  "Bark density",        "bark_density_gcm3",   "balance", "Bark density (g cm⁻³)",
  "Soil-moisture niche", "vwc_realized",        "mcra",    "Soil-moisture niche (% VWC)",
  "Longevity",           "try_plant_longevity", "meth",    "Longevity (years)")
ylab <- c(balance = "Methanogen:methanotroph ratio (log₁₀)", mcra = "mcrA (log₁₀ copies g⁻¹)",
          meth = "Methanotrophs (log₁₀ copies g⁻¹)")
sp <- pmap(sel, function(trait, col, response, xlab) {
  s <- S %>% filter(trait == !!trait, response == !!response)
  dd <- dat %>% transmute(x = .data[[col]], y = .data[[response]], grp = if_else(gymnosperm == 1, "Conifer", "Broadleaf"),
                          lab = unname(SPECIES_ABBREV[species_id])) %>% filter(is.finite(x))
  ggplot(dd, aes(x, y)) +
    geom_smooth(data = filter(dd, grp == "Broadleaf"), method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(shape = grp), size = 2, fill = "white", stroke = 0.7) +
    geom_text_repel(aes(label = lab), size = 2.1, fontface = "italic", seed = 1, box.padding = 0.4,
                    point.padding = 0.15, force = 3, max.overlaps = Inf, min.segment.length = 0.3,
                    segment.size = 0.2, segment.colour = "grey55") +
    scale_shape_manual(values = c(Broadleaf = 16, Conifer = 21), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.14, 0.10))) +   # room for labels at both edges
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
    labs(x = xlab, y = ylab[[response]], subtitle = sprintf("all: ρ = %.2f (p = %.3f)\nbroadleaf: ρ = %.2f (p = %.3f)",
                                                           s$rho, s$p, s$rho_bl, s$p_bl)) +
    theme_bw(base_size = 8.5) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7.5)) })

top <- (h1 | h2) + plot_layout(guides = "collect") & theme(legend.position = "right")
bot <- wrap_plots(sp, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
fig <- wrap_elements(top) / wrap_elements(bot) + plot_layout(heights = c(1.55, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(plot.tag = element_text(face = "bold"))
ggsave(out_path("plant_traits.png"), fig, width = 10, height = 11, dpi = 300, bg = "white")
cat("Wrote plant_traits.png\n")
