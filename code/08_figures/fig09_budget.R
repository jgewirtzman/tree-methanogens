source("code/lib/outputs.R")
# ==============================================================================
# Main figure: plot CH4 budget (maps + seasonal + net waterfall)
# Redesign: color = PROCESS (sink cool / source warm), consistent across panels.
# ==============================================================================
# Design principles:
#   * color encodes CH4 direction: SINK = cool (teal/blue), SOURCE = warm (orange/red),
#     used consistently in maps AND budget so soil reads sink, tree reads source.
#   * spatial -> temporal -> net reading order.
#   * net flux made explicit via a waterfall (0 -> soil -> tree offset -> net).
#   * units: maps + seasonal = nmol m-2 s-1 (rate); budget = mg m-2 yr-1 (integrated).
#   * honesty: whole-surface cases drawn dashed (extrapolated), with the grid range.
#
# Output: outputs/figures/generated/fig_budget_maps.png
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(scales) })
out <- "outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- palette: RdBu family matching the paper's red/blue (sink = blue #2166ac,
#      source = red #b2182b), consistent with the other main figures -----------
BLUES <- c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")  # ~0 -> deep uptake (sink)
REDS  <- c("#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b")  # ~0 -> deep emission (source)
SINK <- "#2166ac"; SRC <- "#b2182b"; SRC_LT <- "#f4a582"; NETC <- "grey40"

# ---- inputs ------------------------------------------------------------------
# ALL FOUR PANELS NOW DRAW ON THE SAME SOURCES. Previously the map panel drew the
# interpolated soil surface while the seasonal and budget panels took soil from
# MONTHLY_FLUXES.csv -- two independent estimates of one quantity, 4.7% apart --
# and the budget panel stacked a 5.41 mg measured band (canonical) under an 85 mg
# whole-surface bar built on a 5.79 mg band (grid). Soil is the interpolated
# surface everywhere; the tree term is the canonical band.
source("code/lib/geometry.R")
soil_map <- read.csv("outputs/tables/soil_surface_annual.csv")     # clipped to the stand
tree_pts <- read.csv("outputs/tables/tree_flux_predictions.csv")
tree_pts <- tree_pts[tree_pts$in_stand & tree_pts$located, ]       # mappable subset
soil_map$mean_flux_nmol <- soil_map$flux_nmol_m2_s

# ---- budget numbers: READ, never hardcoded ----------------------------------
# All values come from 03_budget_canonical.R, which computes them from the locked
# models and the inventory. Run that script first if any input has changed.
B <- read.csv("outputs/data/canonical_budget.csv", stringsAsFactors = FALSE)
val <- function(q) { v <- B$value[B$quantity == q]
  if (!length(v)) stop("canonical_budget.csv is missing: ", q); v }
soil_ann  <- val("soil_mg_m2_yr")
tree_meas <- val("tree_measured_mg_m2_yr")

# whole-woody-surface scenarios, from the sensitivity grid
GR <- read.csv("outputs/data/scaling_full_grid.csv", stringsAsFactors = FALSE)
# THE WHOLE-SURFACE BARS ARE THE NAMED PAIR (scaling_headline.csv, defined in
# 07_scaling_full_grid.R): flux declines above 2 m (exponential fitted to the
# measured three-height means) and flux holds at each stem's 2 m value (constant).
# One whisker shows the grid range. It used to be drawn on a single "illustrative"
# bar, which read as a central estimate; the pair plus the range does not.
HLF <- read.csv("outputs/data/scaling_headline.csv", stringsAsFactors = FALSE)
stopifnot(identical(HLF$case, c("decline", "hold")))
tree_dec <- HLF$total_mg[HLF$case == "decline"]
tree_hold <- HLF$total_mg[HLF$case == "hold"]
# whisker: the grid excluding the uptake form, which is contradicted by the climbed
# tree; the full range including it is printed below for the caption
GRp <- GR[GR$flux != "linear_bounded_median", ]
tree_lo <- min(GRp$total_mg); tree_hi <- max(GRp$total_mg)
tree_grid_lo <- min(GR$total_mg); tree_grid_hi <- max(GR$total_mg)

off <- function(x) 100 * abs(x)/abs(soil_ann)
rf_tree_r2 <- val("tree_r2_oob"); rf_soil_r2 <- val("soil_r2_oob")
n_trees <- nrow(tree_pts)
cat(sprintf("read canonical: soil %.2f | tree measured %.2f | decline %.1f | hold %.1f | grid %.1f-%.1f (all %.1f-%.1f)\n",
            soil_ann, tree_meas, tree_dec, tree_hold, tree_lo, tree_hi, tree_grid_lo, tree_grid_hi))

# shared map extent so a & b align exactly
# MAPPED IN PLOT-LOCAL METRES, not lon/lat. The moisture grid is regular in PX/PY
# but the plot is rotated ~10 degrees from north, so in lon/lat its cells are not
# axis-aligned: geom_raster cannot draw them at all, and geom_tile leaves ragged
# overhangs that made the staircase boundary look like detached chunks. In PX/PY
# the grid is axis-aligned, the rectilinear stand boundary renders exactly, and
# the axes read directly in metres over a 200 m plot -- more useful than a lon/lat
# graticule quoted to three decimals. North is drawn rotated by the plot's own
# rotation angle instead.
xl <- c(0, PLOT_SIDE_M); yl <- c(0, PLOT_SIDE_M)
map_scales <- list(
  scale_x_continuous(breaks = seq(0, 200, 50), expand = c(0, 0)),
  scale_y_continuous(breaks = seq(0, 200, 50), expand = c(0, 0)))
STAND <- stand_ring_local()
th_map <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9, color = "grey30"),
        axis.title = element_text(size = 10, color = "grey20"),
        panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4),
        legend.position = "bottom", legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(1.1, "cm"), legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# ---- (a) soil map (sink = blue; per m2 GROUND) -------------------------------
# geom_TILE, not geom_raster: the grid is regular in plot-local metres but the
# plot is rotated ~10 degrees from north, so in lon/lat the cells are not
# axis-aligned and geom_raster cannot represent them.
pa <- ggplot(soil_map, aes(PX, PY, fill = mean_flux_nmol)) +
  geom_raster() +
  # DIVERGING, centred on zero. The previous scale ran from the minimum to 0, so
  # every cell with a positive flux fell outside it and rendered as grey missing
  # data -- which silently erased the wet zones that are net CH4 SOURCES. They are
  # 17.5% of the stand and reach +0.67 nmol m-2 s-1, and they sit where the soil is
  # wettest (mean VWC 45.9 against 14.1 in the sink cells). Showing them also makes
  # the panel consistent with this figure's own convention: sink cool, source warm.
  scale_fill_gradientn(colours = c(rev(BLUES), REDS[-1]),
                      name = expression("Soil CH"[4]*" flux (nmol m"^-2*" ground s"^-1*")"),
                      limits = max(abs(soil_map$mean_flux_nmol))*c(-1, 1),
                      labels = function(b) sub("^-", "\u2212", format(b, trim = TRUE, drop0trailing = TRUE)),
                      guide = guide_colorbar(title.position = "top")) +
  map_scales + labs(x = "metres", y = "metres") +
  geom_path(data = STAND, aes(PX, PY), inherit.aes = FALSE,
            colour = "grey25", linewidth = 0.45) +
  coord_equal(xlim = c(-6, 206), ylim = c(-6, 206), expand = FALSE) + th_map

# ---- (b) tree map (source = red; per m2 WOODY SURFACE; quantile-spread color) -
tv <- tree_pts$flux_nmol_m2_s
# LINEAR colour scale. This was previously rescaled to equal quantiles of the
# prediction distribution, which is pathological here: RF predictions are clumped
# (30% of stems sit within 1% of one value), so two of the five colour stops
# landed 1% apart and the ramp jumped from pale to mid-red at the median. The map
# read as uniformly hot when most stems sit at the low-middle of the range.
# ASINH colour scaling (2026-07-30). UNITS ARE UNCHANGED -- the legend still reads
# nmol m-2 s-1 and the breaks are round numbers in those units; only the mapping from
# value to colour is transformed. The distribution is strongly right-skewed: median
# 0.102 against a maximum of 0.758, so on a linear ramp the median sat at 13% of the
# range, four fifths of the stems shared the palest two colours, and the spatial
# structure in the bulk was invisible. asinh rather than log because it is well behaved
# as the value approaches zero (stems run down to 0.018) and needs no offset, and it is
# the transform already used for flux axes elsewhere in this figure set.
asinh_col <- scales::trans_new(
  "asinh_flux",
  transform = function(x) asinh(x / 0.1),
  inverse   = function(x) 0.1 * sinh(x))
# Breaks chosen for even spacing ON THE TRANSFORMED scale, anchored at zero. On
# asinh(x/0.1) these sit at 0, 0.48, 0.88, 1.44, 2.10 and 2.71, i.e. gaps of
# 0.48/0.40/0.56/0.65/0.62 -- near-uniform. The old 0.02 stop crowded the bottom of the
# bar (0.20 on the transformed scale, a third of the way to the next tick) and told the
# reader nothing the zero anchor does not.
col_brks <- c(0, 0.05, 0.1, 0.2, 0.4, 0.75)
col_brks <- col_brks[col_brks <= max(tv) * 1.02]

pb <- ggplot(tree_pts %>% arrange(flux_nmol_m2_s), aes(PX, PY, color = flux_nmol_m2_s, size = dbh_m)) +
  geom_point(alpha = 0.9) +
  scale_color_gradientn(colours = REDS, name = expression("Tree CH"[4]*" flux, 0–2 m mean (nmol m"^-2*" bark s"^-1*")"),
                        limits = c(0, max(tv)), trans = asinh_col,
                        breaks = col_brks, labels = col_brks,
                        guide = guide_colorbar(title.position = "top")) +
  scale_size_continuous(range = c(0.15, 1.9), guide = "none") +
  map_scales + labs(x = "metres", y = "metres") +
  geom_path(data = STAND, aes(PX, PY), inherit.aes = FALSE,
            colour = "grey25", linewidth = 0.45) +
  coord_equal(xlim = c(-6, 206), ylim = c(-6, 206), expand = FALSE) + th_map

# ---- north arrow, rotated by the plot's own rotation ------------------------
# The axes are metres, so no scale bar is needed. The plot is rotated from north,
# so the needle is drawn at that angle rather than pointing up the page.
.rot <- local({ e <- new.env()
  load("data/processed/integrated/rf_workflow_input_data_with_2023.RData", envir = e)
  e$rf_workflow_data$rotation_angle })
# Placed in the BOTTOM-RIGHT, inside the uncensused quadrat staircase (PX 140-200,
# PY 0-60), which carries no stems and no soil cells -- genuinely empty space rather
# than somewhere the marker has to compete with data.
#
# The panels used to reserve headroom to ylim 224 purely for a top-right arrow; with
# the marker moved into the staircase that band was empty, so the maps are now square
# at ylim 206 and the stand fills the panel.
# It previously sat at (186, 214) with the "N" at 1.5x the needle length above the
# tip, i.e. y ~ 230, while coord_equal clipped at ylim 224 -- so the letter was drawn
# outside the panel and silently cropped. That is why the arrow had no N.
#
# Style: a conventional shaft-and-head arrow rather than the half-filled compass
# needle, with the N above the tip and a white halo so it reads over any background.
add_ctx <- function(p) {
  cx <- 178; cy <- 22; L <- 13                 # arrow centre and half-length, metres
  rr <- function(dx, dy) c(cx + dx*cos(-.rot) - dy*sin(-.rot),
                           cy + dx*sin(-.rot) + dy*cos(-.rot))
  tail <- rr(0, -L); tip <- rr(0, L); lab <- rr(0, L * 1.62)
  p +
    annotate("segment", x = tail[1], y = tail[2], xend = tip[1], yend = tip[2],
             colour = "grey15", linewidth = 0.55,
             arrow = grid::arrow(length = unit(0.055, "npc"), type = "closed",
                                 angle = 22)) +
    annotate("label", x = lab[1], y = lab[2], label = "N",
             size = 2.9, fontface = "bold", colour = "grey15",
             fill = "white", alpha = 0.75, label.size = 0, label.padding = unit(0.06, "lines"))
}
pa <- add_ctx(pa); pb <- add_ctx(pb)

# ---- (c) MONTHLY series: soil (blue) + tree (red), points joined by lines --------
# Reads the canonical monthly series, which uses the height-integrated 0-2 m band for
# the tree term. MONTHLY_FLUXES.csv carries the superseded single-value tree scaling.
CM <- read.csv("outputs/data/canonical_monthly.csv", stringsAsFactors = FALSE)
# The tree series is MEAN FLUX AT BREAST HEIGHT, per m2 of woody surface -- a
# measured quantity, directly comparable to the map in panel (b). It used to be
# the 0-2 m band's contribution per m2 of GROUND, which exists only because we
# summed stem area and divided by plot area; placed beside soil it invited a
# like-for-like reading when it is diluted by that denominator. Per-ground-area
# belongs in the budget panel, and that is where it now lives exclusively.
mon <- CM %>%
  transmute(month,
            `Soil (per m² ground)` = soil_nmol_m2_s,
            `Tree at 1.3 m (per m² bark)` = tree_bh_nmol_m2_s) %>%
  pivot_longer(-month, names_to = "src", values_to = "nmol") %>%
  mutate(src = factor(src, levels = c("Soil (per m² ground)",
                                      "Tree at 1.3 m (per m² bark)")))
pc <- ggplot(mon, aes(month, nmol, colour = src)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~src, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_colour_manual(values = setNames(c(SINK, SRC), levels(mon$src)), guide = "none") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(labels = function(b) sub("^-", "\u2212", format(b, trim = TRUE, drop0trailing = TRUE))) +
  labs(x = NULL, y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 9, color = "black"),
        panel.grid.minor = element_blank())

# ---- (d) net budget waterfall (mg m-2 yr-1) ----------------------------------
# steps: 0 -> soil -> stems 0-2 m -> whole surface (two alternative cases, each
# starting from the top of the 0-2 m bar) -> net (one half-bar per case)
STEPS <- c("Soil\nuptake", "Stems,\n0–2 m", "Whole\nsurface:\ndeclines\nabove 2 m",
           "Whole\nsurface:\nholds\nabove 2 m", "Net\nbudget", "Foliage\n(unknown)")
base <- soil_ann + tree_meas
HW <- 0.42                      # bar half-width
wf <- tibble(
  x    = c(1, 2, 3, 4, 5 - HW/2, 5 + HW/2),
  hw   = c(HW, HW, HW, HW, HW/2 - 0.02, HW/2 - 0.02),
  ymin = c(soil_ann, soil_ann, base, base, soil_ann + tree_dec, soil_ann + tree_hold),
  ymax = c(0, base, soil_ann + tree_dec, soil_ann + tree_hold, 0, 0),
  fill = c("sink", "src", "src_lt", "src_lt", "net", "net"),
  # whole-surface bars are extrapolations: drawn dashed
  lty  = c("solid", "solid", "dashed", "dashed", "solid", "solid"))
# foliage: an unmeasured term of unknown MAGNITUDE and SIGN; the box is drawn at the
# larger (hold) whole-surface magnitude, purely to signal scale
fol_h <- max(tree_dec, tree_hold)
fol   <- tibble(x = 6)
lab1 <- function(v) sprintf(if (abs(v) < 20) "+%.1f" else "+%.0f", v)
neg  <- function(v) sub("-", "\u2212", sprintf("%.0f", v))
pd <- ggplot(wf) +
  geom_rect(aes(xmin = x - hw, xmax = x + hw, ymin = ymin, ymax = ymax,
                fill = fill, linetype = lty), color = "grey30", linewidth = 0.35) +
  scale_linetype_identity() +
  # connectors: soil -> 0-2 m -> both whole-surface cases
  annotate("segment", x = 1 + HW, xend = 2 - HW, y = soil_ann, yend = soil_ann,
           color = "grey55", linetype = "dashed", linewidth = 0.3) +
  annotate("segment", x = 2 + HW, xend = 4 - HW, y = base, yend = base,
           color = "grey55", linetype = "dashed", linewidth = 0.3) +
  # grid range (200 combinations excluding the uptake form), between the two cases
  annotate("errorbar", x = 3.5, ymin = soil_ann + tree_lo, ymax = soil_ann + tree_hi,
           width = 0.1, colour = "grey25", linewidth = 0.45) +
  annotate("text", x = 3.5, y = soil_ann + tree_hi, label = "grid\nrange",
           vjust = -0.25, size = 2.9, colour = "grey30", lineheight = 0.85) +
  geom_rect(data = fol, aes(xmin = x - HW, xmax = x + HW, ymin = -fol_h, ymax = fol_h),
            fill = "grey88", color = "grey45", linetype = "dashed", linewidth = 0.4,
            inherit.aes = FALSE) +
  annotate("text", x = 6, y = 0, label = "?", size = 4.4, fontface = "bold", color = "grey35") +
  geom_hline(yintercept = 0, color = "grey60") +
  scale_fill_manual(values = c(sink = SINK, src = SRC, src_lt = SRC_LT, net = NETC), guide = "none") +
  scale_x_continuous(breaks = 1:6, labels = STEPS, limits = c(0.45, 6.55), expand = c(0, 0)) +
  scale_y_continuous(labels = function(b) sub("^-", "\u2212", format(b, trim = TRUE, drop0trailing = TRUE))) +
  # labels computed from the variables above -- never hardcode, they drift.
  # Whole-surface labels are TOTALS (they include the 0-2 m band).
  annotate("text", x = 1, y = soil_ann, label = neg(soil_ann),
           vjust = -0.7, size = 3.4, colour = "white", fontface = "bold") +
  annotate("text", x = 2, y = base, label = sprintf("+%.1f", tree_meas), vjust = -0.8, size = 3.4) +
  annotate("text", x = 3, y = soil_ann + tree_dec, label = lab1(tree_dec), vjust = -0.6, size = 3.4) +
  annotate("text", x = 4, y = soil_ann + tree_hold, label = lab1(tree_hold), vjust = -0.6, size = 3.4) +
  # net: one half-bar per case, named and valued inside (too narrow for a level label)
  annotate("text", x = wf$x[5:6], y = (soil_ann + c(tree_dec, tree_hold)) / 2,
           label = sprintf("%s  %s", c("declines", "holds"), neg(soil_ann + c(tree_dec, tree_hold))),
           angle = 90, size = 3.0, colour = "white", fontface = "bold") +
  labs(x = NULL, y = expression("CH"[4]*" (mg m"^-2*" yr"^-1*", per m"^2*" ground)")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 8, lineheight = 0.9),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

# ---- assemble (net-sink / RF-R2 / basis notes go in the figure CAPTION) -------
# panel (d) carries six labelled columns, so the bottom row gives it more width
fig <- (pa | pb) / ((pc | pd) + plot_layout(widths = c(0.85, 1.15))) + plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(out_path("fig_budget_maps.png"), fig, width = 9.2, height = 8, dpi = 300, bg = "white")   # near print size
cat("Wrote fig_budget_maps.png\n")
for (k in seq_len(nrow(HLF)))
  cat(sprintf("%-8s %s | %s | %s | %s : total %.1f (%.1f%% of soil, %.0f%% extrapolated), net %.0f\n",
              HLF$case[k], HLF$flux[k], HLF$WAI[k], HLF$branch[k], HLF$bole[k], HLF$total_mg[k],
              HLF$pct_of_soil[k], HLF$pct_extrapolated[k], HLF$net_mg[k]))
cat(sprintf("Soil %.1f | tree measured %.2f (%.2f%% of soil)\n", soil_ann, tree_meas, off(tree_meas)))
cat(sprintf("grid-range whisker %.1f-%.1f (%d combinations excluding the uptake form)\n",
            tree_lo, tree_hi, nrow(GRp)))
cat(sprintf("full %d-combination range including it: %.1f-%.1f\n",
            nrow(GR), tree_grid_lo, tree_grid_hi))
cat(sprintf("Net budget %.0f (measured only); %.0f (declines) to %.0f (holds) mg CH4 m-2 yr-1 -- sink throughout\n",
            soil_ann + tree_meas, soil_ann + tree_dec, soil_ann + tree_hold))
