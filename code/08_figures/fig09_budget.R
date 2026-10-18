# ==============================================================================
# fig09_budget.R -- Figure 9: the stand CH4 budget
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) annual soil flux, the interpolated surface over the stand (per m2 ground)
#   (b) predicted stem flux, 0-2 m mean, per mapped stem (per m2 bark)
#   (c) monthly soil flux and stem flux at breast height
#   (d) the annual budget as a waterfall (mg CH4 m-2 ground yr-1): soil, stems 0-2 m,
#       the whole woody surface (decline case; whisker = the scaling grid without the
#       uptake form), net, and foliage (unknown)
# Colour encodes direction throughout: sink = blue, source = red (RdBu).
# Every number is read from the upscaling outputs, none typed.
# Inputs
#   outputs/tables/soil_surface_annual.csv, tree_flux_predictions.csv (06_upscale/)
#   outputs/data/canonical_budget.csv     06_upscale/03_budget_canonical.R
#   outputs/data/canonical_monthly.csv    06_upscale/03_budget_canonical.R
#   outputs/data/scaling_full_grid.csv, scaling_headline.csv  06_upscale/07_scaling_full_grid.R
#   data/tables/site.csv (plot rotation); code/lib/geometry.R (stand outline)
# Output: outputs/figures/generated/fig_budget_maps.png; the budget numbers to the console
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(scales) })

# ---- palette: RdBu family matching the paper's red/blue (sink = blue #2166ac,
#      source = red #b2182b), consistent with the other main figures -----------
BLUES <- c("#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")  # ~0 -> deep uptake (sink)
REDS  <- c("#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b")  # ~0 -> deep emission (source)
SINK <- "#2166ac"; SRC <- "#b2182b"; SRC_LT <- "#f4a582"; NETC <- "grey40"

# ---- inputs ------------------------------------------------------------------
# All four panels draw on the same sources: soil is the interpolated surface
# everywhere; the stem term is the canonical band.
source("code/lib/geometry.R")
source("code/lib/flux_transform.R")
soil_map <- read.csv("outputs/tables/soil_surface_annual.csv")     # clipped to the stand
tree_pts <- read.csv("outputs/tables/tree_flux_predictions.csv")
tree_pts <- tree_pts[tree_pts$in_stand & tree_pts$located, ]       # mappable subset
soil_map$mean_flux_nmol <- soil_map$flux_nmol_m2_s

# ---- budget numbers: READ, never hardcoded (06_upscale/03_budget_canonical.R) --
B <- read.csv("outputs/data/canonical_budget.csv", stringsAsFactors = FALSE)
val <- function(q) { v <- B$value[B$quantity == q]
  if (!length(v)) stop("canonical_budget.csv is missing: ", q); v }
soil_ann  <- val("soil_mg_m2_yr")
tree_meas <- val("tree_measured_mg_m2_yr")

# whole-woody-surface scenarios, from the sensitivity grid
GR <- read.csv("outputs/data/scaling_full_grid.csv", stringsAsFactors = FALSE)
# THE NAMED PAIR (scaling_headline.csv, defined in 07_scaling_full_grid.R): flux
# declines above 2 m (exponential fitted to the measured three-height means) or
# holds at each stem's 2 m value (constant). Panel (d) draws the decline case as the
# bar, with the grid range as its whisker; the hold case is given in the caption.
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
cat(sprintf("read canonical: soil %.2f | tree measured %.2f | decline %.1f | hold %.1f | grid %.1f-%.1f (all %.1f-%.1f)\n",
            soil_ann, tree_meas, tree_dec, tree_hold, tree_lo, tree_hi, tree_grid_lo, tree_grid_hi))

# Maps in plot-local metres (PX/PY), not lon/lat: the grid is regular in PX/PY but the
# plot is rotated ~10 degrees from north, so only in PX/PY are the cells axis-aligned
# (geom_raster) and the rectilinear stand boundary exact. North is drawn rotated by
# the plot's own rotation. Shared extent so (a) and (b) align.
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
pa <- ggplot(soil_map, aes(PX, PY, fill = mean_flux_nmol)) +
  geom_raster() +
  # diverging, centred on zero, so wet cells that are net sources show as sources
  scale_fill_gradientn(colours = c(rev(BLUES), REDS[-1]),
                      name = expression("Soil CH"[4]*" flux (nmol m"^-2*" ground s"^-1*")"),
                      limits = max(abs(soil_map$mean_flux_nmol))*c(-1, 1),
                      labels = function(b) sub("^-", "\u2212", format(b, trim = TRUE, drop0trailing = TRUE)),
                      guide = guide_colorbar(title.position = "top")) +
  map_scales + labs(x = "metres", y = "metres") +
  geom_path(data = STAND, aes(PX, PY), inherit.aes = FALSE,
            colour = "grey25", linewidth = 0.45) +
  coord_equal(xlim = c(-6, 206), ylim = c(-6, 206), expand = FALSE) + th_map

# ---- (b) tree map (source = red; per m2 WOODY SURFACE) ------------------------
tv <- tree_pts$flux_nmol_m2_s
# Colour on the arcsinh scale of the flux axes (code/lib/flux_transform.R): the
# predictions are right-skewed, and on a linear ramp most stems share the palest
# colours. Legend units stay nmol m-2 s-1; the breaks are near-evenly spaced on the
# transformed scale, anchored at zero.
asinh_col <- flux_asinh_trans()
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
.rot <- read.csv("data/tables/site.csv")$rotation_rad   # plot-axis rotation from grid north
# Placed bottom-right, inside the uncensused quadrat staircase (PX 140-200, PY 0-60),
# which carries no stems and no soil cells; "N" above the tip on a white halo.
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
# The stem series is mean flux at breast height per m2 of woody surface, comparable
# with the map in (b); per-m2-ground stem flux appears only in the budget panel (d).
CM <- read.csv("outputs/data/canonical_monthly.csv", stringsAsFactors = FALSE)
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
# steps: 0 -> soil -> stems 0-2 m -> whole surface -> net. The whole-surface bar is
# the decline case (the lower of the named pair); the whisker spans the grid
# (excluding the uptake form), and the net bar carries the same range. Whiskers are
# scenario ranges, not uncertainty intervals; the caption says so.
STEPS <- c("Soil\nuptake", "Stems,\n0–2 m", "Whole\nwoody\nsurface", "Net\nbudget", "Foliage\n(unknown)")
base <- soil_ann + tree_meas
HW <- 0.38                      # bar half-width
wf <- tibble(
  x    = 1:4,
  ymin = c(soil_ann, soil_ann, base, soil_ann + tree_dec),
  ymax = c(0, base, soil_ann + tree_dec, 0),
  fill = c("sink", "src", "src_lt", "net"),
  # the whole-surface bar is an extrapolation: drawn dashed
  lty  = c("solid", "solid", "dashed", "solid"))
# foliage: an unmeasured term of unknown MAGNITUDE and SIGN; the box size is
# arbitrary, chosen only to be readable
fol_h <- 40
fol   <- tibble(x = 5)
lab1 <- function(v) sprintf(if (abs(v) < 20) "+%.1f" else "+%.0f", v)
neg  <- function(v) sub("-", "\u2212", sprintf("%.0f", v))
pd <- ggplot(wf) +
  geom_rect(aes(xmin = x - HW, xmax = x + HW, ymin = ymin, ymax = ymax,
                fill = fill, linetype = lty), color = "grey30", linewidth = 0.35) +
  scale_linetype_identity() +
  # connectors: soil -> 0-2 m -> whole surface -> net
  annotate("segment", x = 1 + HW, xend = 2 - HW, y = soil_ann, yend = soil_ann,
           color = "grey55", linetype = "dashed", linewidth = 0.3) +
  annotate("segment", x = 2 + HW, xend = 3 - HW, y = base, yend = base,
           color = "grey55", linetype = "dashed", linewidth = 0.3) +
  annotate("segment", x = 3 + HW, xend = 4 - HW, y = soil_ann + tree_dec, yend = soil_ann + tree_dec,
           color = "grey55", linetype = "dashed", linewidth = 0.3) +
  # scenario range (200 combinations excluding the uptake form), on the whole-surface
  # bar and on the net bar; the net whisker is drawn black-under-white so it reads
  # on the dark bar and on the white background
  annotate("errorbar", x = 3, ymin = soil_ann + tree_lo, ymax = soil_ann + tree_hi,
           width = 0.14, colour = "grey25", linewidth = 0.45) +
  annotate("errorbar", x = 4, ymin = soil_ann + tree_lo, ymax = soil_ann + tree_hi,
           width = 0.14, colour = "black", linewidth = 1.1) +
  annotate("errorbar", x = 4, ymin = soil_ann + tree_lo, ymax = soil_ann + tree_hi,
           width = 0.14, colour = "white", linewidth = 0.45) +
  geom_rect(data = fol, aes(xmin = x - HW, xmax = x + HW, ymin = -fol_h, ymax = fol_h),
            fill = "grey88", color = "grey45", linetype = "dashed", linewidth = 0.4,
            inherit.aes = FALSE) +
  annotate("text", x = 5, y = 0, label = "?", size = 4.4, fontface = "bold", color = "grey35") +
  geom_hline(yintercept = 0, color = "grey60") +
  scale_fill_manual(values = c(sink = SINK, src = SRC, src_lt = SRC_LT, net = NETC), guide = "none") +
  scale_x_continuous(breaks = 1:5, labels = STEPS, limits = c(0.5, 5.5), expand = c(0, 0)) +
  scale_y_continuous(labels = function(b) sub("^-", "\u2212", format(b, trim = TRUE, drop0trailing = TRUE))) +
  # labels computed from the variables above -- never hardcode, they drift.
  # The whole-surface label is the TOTAL (it includes the 0-2 m band) with its
  # scenario range, above the whisker; the net label sits inside the bar, above it.
  annotate("text", x = 1, y = soil_ann, label = neg(soil_ann),
           vjust = -0.7, size = 3.4, colour = "white", fontface = "bold") +
  annotate("text", x = 2, y = base, label = sprintf("+%.1f", tree_meas), vjust = -0.8, size = 3.4) +
  annotate("text", x = 3, y = soil_ann + tree_hi,
           label = sprintf("%s\n(%.1f to %.0f)", lab1(tree_dec), tree_lo, tree_hi),
           vjust = -0.3, size = 3.2, lineheight = 0.9) +
  annotate("text", x = 4, y = (soil_ann + tree_hi) / 2,
           label = sprintf("%s\n(%s to\n%s)", neg(soil_ann + tree_dec), neg(soil_ann + tree_lo), neg(soil_ann + tree_hi)),
           size = 3.2, lineheight = 0.9, colour = "white", fontface = "bold") +
  labs(x = NULL, y = expression("CH"[4]*" (mg m"^-2*" yr"^-1*", per m"^2*" ground)")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 8, lineheight = 0.9),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

# ---- assemble (net-sink and basis notes are in the caption) --------------------
fig <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1.05, 1)) +
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
cat(sprintf("Net budget %.0f (measured only); %.0f (declines, drawn) to %.0f (holds); grid %.0f to %.0f mg CH4 m-2 yr-1 -- sink throughout\n",
            soil_ann + tree_meas, soil_ann + tree_dec, soil_ann + tree_hold, soil_ann + tree_grid_hi, soil_ann + tree_grid_lo))
