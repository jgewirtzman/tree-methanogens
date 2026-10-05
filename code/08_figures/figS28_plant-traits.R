source("code/lib/outputs.R")
# ==============================================================================
# Fig S25: plant traits and the stem methane-cycling community
#   (a) Spearman correlations of the 16 rule-selected traits with mcrA, pmoA + mmoX,
#       the methanogen:methanotroph balance and net flux, across all ten gene-flux
#       species (left) and the eight broadleaf species alone (right); * p < 0.05.
#   (b) the three associations named in Results §9, conifers as open symbols.
# Correlations come from code/09_tables_stats/25_plant-traits.R (traits_correlations.csv),
# which also documents the trait inclusion rules. Descriptive, n = 10 species.
# (File name kept from the earlier numbering; it draws Figure S25.)
# Output: outputs/figures/generated/plant_traits.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(ggrepel) })

S <- read_csv("outputs/data/traits_correlations.csv", show_col_types = FALSE)
inc <- read_csv("outputs/data/traits_included.csv", show_col_types = FALSE)
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",BELE="Betula lenta",
  BEPA="Betula papyrifera",FAGR="Fagus grandifolia",FRAM="Fraxinus americana",
  PIST="Pinus strobus",QURU="Quercus rubra",TSCA="Tsuga canadensis")
invisible(capture.output(suppressMessages(source("code/lib/prep_species_data.R"))))
resp <- analysis_ratio %>% transmute(species_id, balance = median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, mcra = log10(value + 1)), by = "species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, meth = log10(value + 1)), by = "species_id")
tr <- read_csv("data/raw/external/tree-gas-traits/ymf_species_traits.csv", show_col_types = FALSE)
niche <- read_csv("outputs/data/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(spcode = names(sp_map)[match(species, sp_map)]) %>% select(spcode, vwc_realized = vwc_mean)
dat <- resp %>% left_join(tr, by = c("species_id" = "spcode")) %>% left_join(niche, by = c("species_id" = "spcode"))

# ---- (a) twin heatmaps --------------------------------------------------------------
rlab <- c(mcra = "mcrA", meth = "Methano-\ntrophs", balance = "Balance", flux = "Stem CH4\nflux")
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
          strip.placement = "outside", axis.text.x.top = element_text(size = 8.5, lineheight = 0.9), panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 4))) }
h1 <- hm(rho, p, "All ten species") + theme(legend.position = "none")
h2 <- hm(rho_bl, p_bl, "Eight broadleaf species", strip = FALSE) + theme(axis.text.y = element_blank())

# ---- (b) the associations named in the text ----------------------------------------------
abbr <- function(id) vapply(strsplit(sp_map[id], " "), function(x) paste0(substr(x[1], 1, 1), ". ", x[2]), "")
sel <- tribble(~trait, ~col, ~response, ~xlab,
  "Bark density",        "bark_density_gcm3",   "balance", "Bark density (g cm⁻³)",
  "Soil-moisture niche", "vwc_realized",        "mcra",    "Soil-moisture niche (% VWC)",
  "Longevity",           "try_plant_longevity", "meth",    "Longevity (years)")
ylab <- c(balance = "Balance (log₁₀ mcrA:methanotroph)", mcra = "mcrA (log₁₀ copies g⁻¹)",
          meth = "Methanotrophs (log₁₀ copies g⁻¹)")
sp <- pmap(sel, function(trait, col, response, xlab) {
  s <- S %>% filter(trait == !!trait, response == !!response)
  dd <- dat %>% transmute(x = .data[[col]], y = .data[[response]], grp = if_else(gymnosperm == 1, "Conifer", "Broadleaf"),
                          lab = abbr(species_id)) %>% filter(is.finite(x))
  ggplot(dd, aes(x, y)) +
    geom_smooth(data = filter(dd, grp == "Broadleaf"), method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(shape = grp), size = 2, fill = "white", stroke = 0.7) +
    geom_text_repel(aes(label = lab), size = 2.1, fontface = "italic", seed = 1, box.padding = 0.25) +
    scale_shape_manual(values = c(Broadleaf = 16, Conifer = 21), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
    labs(x = xlab, y = ylab[[response]], subtitle = sprintf("all: ρ = %.2f (p = %.3f)\nbroadleaf: ρ = %.2f (p = %.3f)",
                                                           s$rho, s$p, s$rho_bl, s$p_bl)) +
    theme_bw(base_size = 8.5) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 7.5)) })

top <- (h1 | h2) + plot_layout(guides = "collect") & theme(legend.position = "right")
bot <- wrap_plots(sp, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
fig <- wrap_elements(top) / wrap_elements(bot) + plot_layout(heights = c(1.55, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(plot.tag = element_text(face = "bold"))
ggsave(out_path("plant_traits.png"), fig, width = 10, height = 11, dpi = 300, bg = "white")
cat("Wrote plant_traits.png\n")
