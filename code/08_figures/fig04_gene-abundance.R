source("code/lib/outputs.R")
# ==============================================================================
# Fig 4: methanogen (mcrA) and methanotroph (pmoA + mmoX) gene abundance
#   (a) by species and compartment, for the species in the gene-flux analyses
#       (Results §9, Fig 8: at least five trees with genes and flux; the set is taken
#       from code/07_molecular/helper_scale_dependent_gene_patterns.R)
#   (b) mcrA against pmoA + mmoX for every sample, with a 1:1 line
# Also writes the all-species version with pmoA and mmoX separate (SI figure).
#
# UNITS: the ddpcr_*_loose columns are ALREADY copies g^-1, converted once in
# code/03_merge/04_harmonize_all_data.R with code/lib/ddpcr_constants.R
# (Conc x (25/2.5) x elution (75 wood, 100 soil) x (800/250) / mass). Basis: DRY for wood
# (freeze-dried cores); soil uses fresh sample mass.
# Points in (a) are means of log10 among positive samples (as in the text) +- SE;
# hollow points mark compartments where fewer than half of samples were positive.
# Outputs: outputs/figures/generated/fig4_final.png, figS_gene-abundance-species.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(cowplot); library(patchwork); library(ggside) })

mf <- read_csv("data/processed/integrated/merged_tree_dataset_final.csv", show_col_types = FALSE)
source("code/08_figures/helper_species_barplots.R")   # species_mapping
invisible(capture.output(suppressMessages(source("code/07_molecular/helper_scale_dependent_gene_patterns.R"))))
flux_species <- analysis_mcra$species                     # the gene-flux species set

comp <- c(Inner = "Heartwood", Outer = "Sapwood", Mineral = "Mineral soil", Organic = "Organic soil")
cols <- c(Heartwood = "#a6611a", Sapwood = "#dfc27d", `Mineral soil` = "#80cdc1", `Organic soil` = "#018571")
wide <- map_dfr(names(comp), function(l) {
  p <- mf[[paste0("ddpcr_pmoa_", l, "_loose")]]; m <- mf[[paste0("ddpcr_mmox_", l, "_loose")]]
  tibble(tree_id = mf$tree_id, species = unname(species_mapping[mf$species_id]), compartment = comp[[l]],
         mcrA = mf[[paste0("ddpcr_mcra_probe_", l, "_loose")]], pmoA = p, mmoX = m,
         mt = ifelse(is.na(p) & is.na(m), NA, coalesce(p, 0) + coalesce(m, 0)))
}) %>% filter(!is.na(species)) %>% mutate(compartment = factor(compartment, levels = names(cols)))

gene_labels <- c(mcrA = "italic(mcrA)", mt = "italic(pmoA)+italic(mmoX)", pmoA = "italic(pmoA)", mmoX = "italic(mmoX)")
species_dots <- function(dat, genes) {
  long <- dat %>% pivot_longer(all_of(genes), names_to = "gene", values_to = "copies") %>% filter(!is.na(copies))
  s <- long %>% group_by(species, compartment, gene) %>%
    summarise(n = n(), npos = sum(copies > 0),
              mean = if (npos > 0) mean(log10(copies[copies > 0])) else NA_real_,
              se = if (npos > 1) sd(log10(copies[copies > 0])) / sqrt(npos) else NA_real_, .groups = "drop") %>%
    mutate(detected = npos / n >= 0.5)
  trees <- long %>% distinct(species, tree_id) %>% count(species, name = "trees")
  s <- s %>% left_join(trees, by = "species") %>% mutate(label = paste0("italic('", species, "')~(", trees, ")"))
  ord <- s %>% filter(gene == "mcrA", compartment == "Heartwood") %>% arrange(mean) %>% pull(label)
  s %>% mutate(label = factor(label, levels = unique(c(ord, label))),
               gene = factor(gene, levels = genes, labels = paste0(gene_labels[genes], "~(copies~g^-1)")))
}
dot_plot <- function(s) {
  # each compartment gets a fixed vertical offset within its species row, so points and bars share it
  off <- c(Heartwood = 0.27, Sapwood = 0.09, `Mineral soil` = -0.09, `Organic soil` = -0.27)
  labs_y <- levels(s$label)
  s <- s %>% mutate(yy = as.numeric(label) + off[as.character(compartment)])
  ggplot(s, aes(mean, yy, colour = compartment)) +
    geom_errorbar(aes(xmin = mean - se, xmax = mean + se), orientation = "y", width = 0, linewidth = 0.5, na.rm = TRUE) +
    geom_point(aes(shape = detected), size = 2.4, stroke = 0.9, fill = "white", na.rm = TRUE) +
    facet_grid(. ~ gene, scales = "free_x", labeller = label_parsed) +
    scale_y_continuous(breaks = seq_along(labs_y), labels = parse(text = labs_y),
                       limits = c(0.55, length(labs_y) + 0.45), expand = c(0, 0)) +
    scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), name = NULL,
                       labels = c(`TRUE` = "detected in at least half of samples", `FALSE` = "detected in fewer than half")) +
    scale_x_continuous(labels = function(b) parse(text = paste0("10^", b)), breaks = 1:9) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "grey92"),
          strip.background = element_blank(), legend.position = "none")
}
legend_of <- function(p) get_plot_component(
  p + theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 10),
            legend.spacing.y = unit(0, "pt")) +
    guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 3)), shape = guide_legend(order = 2)),
  "guide-box-bottom")

# ---------- (a) gene-flux species ----------
pa <- dot_plot(species_dots(filter(wide, species %in% flux_species), c("mcrA", "mt")))

# ---------- (b) all samples: non-detects of mcrA in their own strip, 1:1 line ----------
sc <- wide %>% filter(!is.na(mcrA), !is.na(mt), mt > 0)          # one pmoA + mmoX non-detect left out
lx <- log10(sc$mcrA[sc$mcrA > 0]); ly <- log10(sc$mt)
ndx <- floor(min(lx)) - 0.8
set.seed(1)
sc <- sc %>% mutate(x = if_else(mcrA > 0, log10(mcrA), ndx), y = log10(mt),
                    xj = if_else(mcrA > 0, x, x + runif(n(), -0.18, 0.18)))
cen <- sc %>% group_by(compartment) %>%
  summarise(mx = mean(log10(mcrA[mcrA > 0])), sx = sd(log10(mcrA[mcrA > 0])), my = mean(y), sy = sd(y), .groups = "drop")
lo <- max(min(lx), min(ly)); hi <- min(max(lx), max(ly))
pb <- ggplot(sc, aes(xj, y, colour = compartment)) +
  annotate("segment", x = lo, xend = hi, y = lo, yend = hi, linetype = "dashed", colour = "grey40") +
  annotate("text", x = hi, y = hi, label = "1:1", hjust = -0.2, vjust = 1.4, size = 3.3, colour = "grey30") +
  geom_vline(xintercept = ndx + 0.45, colour = "grey75", linewidth = 0.3) +
  geom_point(size = 1.8, alpha = 0.5) +
  geom_errorbar(data = cen, aes(x = mx, ymin = my - sy, ymax = my + sy), inherit.aes = FALSE,
                width = 0, linewidth = 0.9, colour = "black") +
  geom_errorbar(data = cen, aes(y = my, xmin = mx - sx, xmax = mx + sx), orientation = "y", inherit.aes = FALSE,
                width = 0, linewidth = 0.9, colour = "black") +
  geom_point(data = cen, aes(mx, my, fill = compartment), inherit.aes = FALSE, shape = 21, size = 4.2, stroke = 1.1) +
  geom_xsidedensity(data = filter(sc, mcrA > 0), aes(x = x, y = after_stat(density), fill = compartment),
                    inherit.aes = FALSE, alpha = 0.5, colour = NA, adjust = 1.5) +
  geom_ysidedensity(aes(y = y, x = after_stat(density), fill = compartment),
                    inherit.aes = FALSE, alpha = 0.5, colour = NA, adjust = 1.5) +
  scale_x_continuous(breaks = c(ndx, seq(ceiling(min(lx)), floor(max(lx)), 1)),
                     labels = function(b) parse(text = ifelse(b == ndx, "n.d.", paste0("10^", b))),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(breaks = seq(ceiling(min(ly)), floor(max(ly)), 1),
                     labels = function(b) parse(text = paste0("10^", b)), expand = expansion(mult = 0.03)) +
  scale_colour_manual(values = cols, guide = "none") + scale_fill_manual(values = cols, guide = "none") +
  labs(x = expression(italic(mcrA)~"(copies g"^-1*")"), y = expression(italic(pmoA)+italic(mmoX)~"(copies g"^-1*")")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 0.72,
        ggside.panel.scale = 0.22, ggside.axis.text = element_blank(), ggside.axis.ticks = element_blank(),
        ggside.panel.border = element_blank(), ggside.panel.grid = element_blank(), ggside.panel.background = element_blank())

right <- plot_grid(pb, legend_of(pa), ncol = 1, rel_heights = c(1, 0.32), labels = c("(b)", ""),
                   label_size = 13, label_fontface = "bold")
fig <- plot_grid(pa, right, ncol = 2, rel_widths = c(1.15, 1), labels = c("(a)", ""), label_size = 13, label_fontface = "bold") +
  theme(plot.background = element_rect(fill = "white", colour = NA))
OUT_FIG4 <- if (exists("OUT_FIG4")) OUT_FIG4 else out_path("fig4_final.png")
ggsave(OUT_FIG4, fig, width = 13, height = 5.6, dpi = 300, bg = "white")

# ---------- SI: all species, mcrA, pmoA and mmoX separately ----------
ps <- dot_plot(species_dots(wide, c("mcrA", "pmoA", "mmoX")))
si <- plot_grid(ps, legend_of(ps), ncol = 1, rel_heights = c(1, 0.12)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
OUT_SI <- if (exists("OUT_SI")) OUT_SI else out_path("figS_gene-abundance-species.png")
ggsave(OUT_SI, si, width = 12, height = 8, dpi = 300, bg = "white")
cat("Wrote Fig 4 and the all-species SI figure\n")
