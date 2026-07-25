# ==============================================================================
# REVISION — NEW main figure: plot CH4 budget (maps + seasonal + net waterfall)
# Redesign: color = PROCESS (sink cool / source warm), consistent across panels.
# ==============================================================================
# Design principles (Jon):
#   * color encodes CH4 direction: SINK = cool (teal/blue), SOURCE = warm (orange/red),
#     used consistently in maps AND budget so soil reads sink, tree reads source.
#   * spatial -> temporal -> net reading order.
#   * net flux made explicit via a waterfall (0 -> soil -> tree offset -> net).
#   * units: maps + seasonal = nmol m-2 s-1 (rate); budget = mg m-2 yr-1 (integrated).
#   * honesty: bounding scenario flagged unconstrained; RF OOB R2 noted.
#
# NEW file. Run: Rscript code/revision/R_fig_component_budget.R
# Output: outputs/revision/fig_budget_maps.png
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(scales) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- palette: RdBu family matching the paper's red/blue (sink = blue #2166ac,
#      source = red #b2182b), consistent with the other main figures -----------
BLUES <- c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")  # ~0 -> deep uptake (sink)
REDS  <- c("#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b")  # ~0 -> deep emission (source)
SINK <- "#2166ac"; SRC <- "#b2182b"; SRC_LT <- "#f4a582"; NETC <- "grey40"

soil_map <- read.csv("outputs/tables/soil_flux_extended_annual.csv")   # x,y,mean_flux_nmol
tree_pts <- read.csv("outputs/tables/tree_flux_predictions.csv")       # x,y,dbh_m,flux_nmol_m2_s
mf       <- read.csv("outputs/tables/MONTHLY_FLUXES.csv")

# budget numbers (manuscript-consistent), mg m-2 yr-1
soil_ann <- -904.30; tree_meas <- 1.25; tree_scen <- 114
off <- function(x) 100 * abs(x)/abs(soil_ann)
rf_tree_r2 <- 0.15; rf_soil_r2 <- 0.28; n_trees <- nrow(tree_pts)

# shared map extent so a & b align exactly
xl <- range(c(soil_map$x, tree_pts$x)); yl <- range(c(soil_map$y, tree_pts$y))
# lat/long graticule labels (coords are decimal degrees; lon W, lat N)
deg_lon <- function(x) paste0(formatC(abs(x), format = "f", digits = 3), " W")
deg_lat <- function(y) paste0(formatC(y,      format = "f", digits = 3), " N")
map_scales <- list(
  scale_x_continuous(breaks = scales::breaks_pretty(3), labels = deg_lon),
  scale_y_continuous(breaks = scales::breaks_pretty(3), labels = deg_lat))
th_map <- theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 5.5, color = "grey45"), axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4),
        legend.position = "bottom", legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.9, "cm"), legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

# ---- (a) soil map (sink = blue; per m2 GROUND) -------------------------------
pa <- ggplot(soil_map, aes(x, y, fill = mean_flux_nmol)) + geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(colours = rev(BLUES), name = expression("Soil CH"[4]*"  (nmol m"^-2*" s"^-1*", per m"^2*" ground)"),
                      limits = c(min(soil_map$mean_flux_nmol), 0),
                      guide = guide_colorbar(title.position = "top")) +
  map_scales +
  coord_quickmap(xlim = xl, ylim = c(yl[1], yl[2] + 0.07*diff(yl)), expand = FALSE) + th_map  # top strip for N arrow

# ---- (b) tree map (source = red; per m2 WOODY SURFACE; quantile-spread color) -
tv <- tree_pts$flux_nmol_m2_s
pb <- ggplot(tree_pts %>% arrange(flux_nmol_m2_s), aes(x, y, color = flux_nmol_m2_s, size = dbh_m)) +
  geom_point(alpha = 0.9) +
  scale_color_gradientn(colours = REDS, name = expression("Tree CH"[4]*"  (nmol m"^-2*" s"^-1*", per m"^2*" woody surface)"),
                        values = scales::rescale(quantile(tv, probs = seq(0, 1, length.out = length(REDS)))),
                        guide = guide_colorbar(title.position = "top")) +
  scale_size_continuous(range = c(0.15, 1.9), guide = "none") +
  map_scales +
  coord_quickmap(xlim = xl, ylim = c(yl[1], yl[2] + 0.07*diff(yl)), expand = FALSE) + th_map  # top strip for N arrow

# ---- add scale bar (50 m) + traditional north arrow to both maps -------------
m_per_deg_lon <- 111320 * cos(mean(yl) * pi/180)
sb <- 50 / m_per_deg_lon                                  # 50 m in degrees lon
dx <- diff(xl); dy <- diff(yl)
add_ctx <- function(p) {
  # scale bar: bottom-left corner (white gap)
  x0 <- xl[1] + 0.05*dx; y0 <- yl[1] + 0.04*dy
  # traditional two-tone compass needle, in the white header strip ABOVE the raster
  ax <- xl[2] - 0.08*dx; apex <- yl[2] + 0.048*dy; base <- yl[2] + 0.008*dy; w <- 0.022*dx
  left  <- data.frame(x = c(ax, ax - w, ax),            y = c(apex, base, base + 0.028*dy))
  right <- data.frame(x = c(ax, ax + w, ax),            y = c(apex, base, base + 0.028*dy))
  p +
    annotate("segment", x = x0, xend = x0 + sb, y = y0, yend = y0, linewidth = 1, color = "grey10") +
    annotate("text", x = x0 + sb/2, y = y0, label = "50 m", vjust = -0.7, size = 2.3, color = "grey10") +
    annotate("polygon", x = left$x,  y = left$y,  fill = "grey10", color = "grey10", linewidth = 0.2) +
    annotate("polygon", x = right$x, y = right$y, fill = "white",  color = "grey10", linewidth = 0.2) +
    annotate("text", x = ax, y = apex, label = "N", vjust = -0.35, size = 2.6, fontface = "bold", color = "grey10")
}
pa <- add_ctx(pa); pb <- add_ctx(pb)

# ---- (c) seasonal (4 seasons): soil (blue) + tree (red), faceted own scales ---
seas_of <- function(m) c("Winter","Winter","Spring","Spring","Spring","Summer","Summer","Summer",
                         "Fall","Fall","Fall","Winter")[m]
seas <- mf %>% mutate(season = factor(seas_of(month), levels = c("Winter","Spring","Summer","Fall"))) %>%
  group_by(season) %>%
  summarise(Soil = mean(Phi_soil_umol_m2_s)*1000, Tree = mean(Phi_tree_umol_m2_s)*1000, .groups="drop") %>%
  pivot_longer(c(Soil, Tree), names_to = "src", values_to = "nmol") %>%
  mutate(src = factor(src, levels = c("Soil","Tree")))
pc <- ggplot(seas, aes(season, nmol, fill = src)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  facet_wrap(~src, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_fill_manual(values = c(Soil = SINK, Tree = SRC), guide = "none") +
  labs(x = NULL, y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*", per m"^2*" ground)")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 8.5, color = "black"),
        panel.grid.minor = element_blank())

# ---- (d) net budget waterfall (mg m-2 yr-1) ----------------------------------
# steps: baseline 0 -> soil -> tree(measured) -> tree(scenario) -> NET(total)
wf <- tibble(
  step  = factor(c("Soil\nuptake","Tree to 2 m\n(measured)","Tree full woody\nsurface (scenario)","Net\nbudget"),
                 levels = c("Soil\nuptake","Tree to 2 m\n(measured)","Tree full woody\nsurface (scenario)","Net\nbudget")),
  ymin  = c(soil_ann, soil_ann, soil_ann + tree_meas, 0),
  ymax  = c(0, soil_ann + tree_meas, soil_ann + tree_meas + (tree_scen - tree_meas),
            soil_ann + tree_scen),
  fill  = c("sink","src","src_lt","net"),
  # scenario bar (whole woody surface, constant-flux extrapolation) drawn dashed = hypothetical
  lty   = c("solid","solid","dashed","solid"))
netmeas <- soil_ann + tree_meas    # net if only measured tree offset
# foliage: an unmeasured term of unknown MAGNITUDE and SIGN (illustrative half-height only)
fol_h <- tree_scen
fol   <- tibble(x = 5)
FOLIAGE_STYLE <- "dashed"   # "dashed" (hard box; preferred) or "gradient" (alpha fades with |flux| — too noisy)
if (FOLIAGE_STYLE == "dashed") {
  fol_layer <- list(geom_rect(data = fol, aes(xmin = x-0.42, xmax = x+0.42, ymin = -fol_h, ymax = fol_h),
                    fill = "grey88", color = "grey45", linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE))
} else {                       # opaque at zero -> increasingly transparent away from zero (both signs)
  .bs   <- seq(-fol_h, fol_h, length.out = 121)
  fgrad <- tibble(x = 5, ymin = head(.bs, -1), ymax = tail(.bs, -1))
  fgrad$a <- 1 - abs((fgrad$ymin + fgrad$ymax) / 2) / fol_h
  fol_layer <- list(
    geom_rect(data = fgrad, aes(xmin = x-0.42, xmax = x+0.42, ymin = ymin, ymax = ymax, alpha = a),
              fill = "grey35", inherit.aes = FALSE),
    scale_alpha_identity())
}
pd <- ggplot(wf) +
  geom_rect(aes(xmin = as.numeric(step)-0.42, xmax = as.numeric(step)+0.42,
                ymin = ymin, ymax = ymax, fill = fill, linetype = lty),
            color = "grey30", linewidth = 0.35) +
  scale_linetype_identity() +
  geom_segment(aes(x = as.numeric(step)+0.42, xend = as.numeric(step)+1-0.42, y = ymax, yend = ymax),
               data = wf[1:3,], color = "grey55", linetype = "dashed", linewidth = 0.3) +
  # foliage — bar straddling zero (unknown sign); style set by FOLIAGE_STYLE above
  fol_layer +
  annotate("text", x = 5, y = 0, label = "?", size = 4.4, fontface = "bold",
           color = if (FOLIAGE_STYLE == "gradient") "white" else "grey35") +
  geom_hline(yintercept = 0, color = "grey60") +
  scale_fill_manual(values = c(sink = SINK, src = SRC, src_lt = SRC_LT, net = NETC), guide = "none") +
  scale_x_continuous(breaks = 1:5, labels = c(levels(wf$step), "Foliage\n(unknown)"), limits = c(0.4, 5.6)) +
  annotate("text", x = 1, y = soil_ann,             label = "-904", vjust = 1.5, size = 2.9) +
  annotate("text", x = 2, y = soil_ann + tree_meas, label = "+1.3", vjust = -0.8, size = 2.9) +
  annotate("text", x = 3, y = soil_ann + tree_scen, label = "+113", vjust = -0.8, size = 2.9) +
  annotate("text", x = 4, y = soil_ann + tree_scen, label = "-790", vjust = 1.5, size = 2.9) +
  labs(x = NULL, y = expression("CH"[4]*" (mg m"^-2*" yr"^-1*", per m"^2*" ground)")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(size = 7.5), panel.grid.minor = element_blank())

# ---- assemble (net-sink / RF-R2 / basis notes go in the figure CAPTION) -------
fig <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(out, "fig_budget_maps.png"), fig, width = 10.5, height = 9, dpi = 300, bg = "white")
cat("Wrote fig_budget_maps.png\n")
cat(sprintf("Soil %.0f | tree measured %.2f (%.2f%%) | tree scenario %d (~%.0f%%) | net -903 to %.0f mg/m2/yr\n",
            soil_ann, tree_meas, off(tree_meas), tree_scen, off(tree_scen), soil_ann+tree_scen))
