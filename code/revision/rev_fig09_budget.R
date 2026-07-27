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

# ---- inputs ------------------------------------------------------------------
# ALL FOUR PANELS NOW DRAW ON THE SAME SOURCES. Previously the map panel drew the
# interpolated soil surface while the seasonal and budget panels took soil from
# MONTHLY_FLUXES.csv -- two independent estimates of one quantity, 4.7% apart --
# and the budget panel stacked a 5.41 mg measured band (canonical) under an 85 mg
# whole-surface bar built on a 5.79 mg band (grid). Soil is the interpolated
# surface everywhere (decision 2026-07-27); the tree term is the canonical band.
source("code/revision/rev_geometry.R")
soil_map <- read.csv("outputs/tables/soil_surface_annual.csv")     # clipped to the stand
tree_pts <- read.csv("outputs/tables/tree_flux_predictions.csv")
tree_pts <- tree_pts[tree_pts$in_stand & tree_pts$located, ]       # mappable subset
soil_map$mean_flux_nmol <- soil_map$flux_nmol_m2_s

# ---- budget numbers: READ, never hardcoded ----------------------------------
# All values come from rev_budget_canonical.R, which computes them from the locked
# models and the inventory. Run that script first if any input has changed.
B <- read.csv("outputs/revision/canonical_budget.csv", stringsAsFactors = FALSE)
val <- function(q) { v <- B$value[B$quantity == q]
  if (!length(v)) stop("canonical_budget.csv is missing: ", q); v }
soil_ann  <- val("soil_mg_m2_yr")
tree_meas <- val("tree_measured_mg_m2_yr")

# whole-woody-surface scenarios, from the sensitivity grid
GR <- read.csv("outputs/revision/scaling_full_grid.csv", stringsAsFactors = FALSE)
CF <- GR[GR$flux == "constant", ]
tree_scen <- median(CF$total_mg[grepl("2.11", CF$WAI)])   # bottom-up WAI for this stand
tree_lo   <- min(CF$total_mg); tree_hi <- max(CF$total_mg)
tree_grid_lo <- min(GR$total_mg); tree_grid_hi <- max(GR$total_mg)
pct_extrapolated <- median(CF$pct_extrapolated[grepl("2.11", CF$WAI)])

off <- function(x) 100 * abs(x)/abs(soil_ann)
rf_tree_r2 <- val("tree_r2_oob"); rf_soil_r2 <- val("soil_r2_oob")
n_trees <- nrow(tree_pts)
cat(sprintf("read canonical: soil %.2f | tree measured %.2f | scenario %.1f (%.1f-%.1f) | grid %.1f-%.1f\n",
            soil_ann, tree_meas, tree_scen, tree_lo, tree_hi, tree_grid_lo, tree_grid_hi))

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
th_map <- theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        axis.text = element_text(size = 6, color = "grey45"),
        axis.title = element_text(size = 6.5, color = "grey35"),
        panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4),
        legend.position = "bottom", legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.9, "cm"), legend.title = element_text(size = 7),
        legend.text = element_text(size = 6))

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
                      name = expression("Soil CH"[4]*"  (nmol m"^-2*" s"^-1*", per m"^2*" ground)"),
                      limits = max(abs(soil_map$mean_flux_nmol))*c(-1, 1),
                      guide = guide_colorbar(title.position = "top")) +
  map_scales + labs(x = "metres", y = "metres") +
  geom_path(data = STAND, aes(PX, PY), inherit.aes = FALSE,
            colour = "grey25", linewidth = 0.45) +
  coord_equal(xlim = c(-6, 206), ylim = c(-6, 224), expand = FALSE) + th_map

# ---- (b) tree map (source = red; per m2 WOODY SURFACE; quantile-spread color) -
tv <- tree_pts$flux_nmol_m2_s
# LINEAR colour scale. This was previously rescaled to equal quantiles of the
# prediction distribution, which is pathological here: RF predictions are clumped
# (30% of stems sit within 1% of one value), so two of the five colour stops
# landed 1% apart and the ramp jumped from pale to mid-red at the median. The map
# read as uniformly hot when most stems sit at the low-middle of the range.
pb <- ggplot(tree_pts %>% arrange(flux_nmol_m2_s), aes(PX, PY, color = flux_nmol_m2_s, size = dbh_m)) +
  geom_point(alpha = 0.9) +
  scale_color_gradientn(colours = REDS, name = expression("Tree CH"[4]*"  (nmol m"^-2*" s"^-1*", per m"^2*" woody surface, 0-2 m band mean)"),
                        limits = c(0, max(tv)),
                        guide = guide_colorbar(title.position = "top")) +
  scale_size_continuous(range = c(0.15, 1.9), guide = "none") +
  map_scales + labs(x = "metres", y = "metres") +
  geom_path(data = STAND, aes(PX, PY), inherit.aes = FALSE,
            colour = "grey25", linewidth = 0.45) +
  coord_equal(xlim = c(-6, 206), ylim = c(-6, 224), expand = FALSE) + th_map

# ---- north arrow, rotated by the plot's own rotation ------------------------
# The axes are metres, so no scale bar is needed. The plot is rotated from north,
# so the needle is drawn at that angle rather than pointing up the page.
.rot <- local({ e <- new.env()
  load("data/processed/integrated/rf_workflow_input_data_with_2023.RData", envir = e)
  e$rf_workflow_data$rotation_angle })
add_ctx <- function(p) {
  cx <- 186; cy <- 214; L <- 11; w <- 4.2      # needle centre and size, in metres
  rr <- function(dx, dy) c(cx + dx*cos(-.rot) - dy*sin(-.rot),
                           cy + dx*sin(-.rot) + dy*cos(-.rot))
  tip <- rr(0, L); base <- rr(0, -L); lw <- rr(-w, -L*0.35); rw <- rr(w, -L*0.35)
  p +
    annotate("polygon", x = c(tip[1], lw[1], base[1]), y = c(tip[2], lw[2], base[2]),
             fill = "grey15", colour = "grey15", linewidth = 0.2) +
    annotate("polygon", x = c(tip[1], rw[1], base[1]), y = c(tip[2], rw[2], base[2]),
             fill = "white", colour = "grey15", linewidth = 0.2) +
    annotate("text", x = rr(0, L*1.5)[1], y = rr(0, L*1.5)[2], label = "N",
             size = 2.6, fontface = "bold", colour = "grey15")
}
pa <- add_ctx(pa); pb <- add_ctx(pb)

# ---- (c) MONTHLY series: soil (blue) + tree (red), points joined by lines --------
# Reads the canonical monthly series, which uses the height-integrated 0-2 m band for
# the tree term. MONTHLY_FLUXES.csv carries the superseded single-value tree scaling.
CM <- read.csv("outputs/revision/canonical_monthly.csv", stringsAsFactors = FALSE)
# The tree series is MEAN FLUX AT BREAST HEIGHT, per m2 of woody surface -- a
# measured quantity, directly comparable to the map in panel (b). It used to be
# the 0-2 m band's contribution per m2 of GROUND, which exists only because we
# summed stem area and divided by plot area; placed beside soil it invited a
# like-for-like reading when it is diluted by that denominator. Per-ground-area
# belongs in the budget panel, and that is where it now lives exclusively.
mon <- CM %>%
  transmute(month,
            `Soil\n(per sq.m ground)` = soil_nmol_m2_s,
            `Tree at 1.3 m\n(per sq.m bark)` = tree_bh_nmol_m2_s) %>%
  pivot_longer(-month, names_to = "src", values_to = "nmol") %>%
  mutate(src = factor(src, levels = c("Soil\n(per sq.m ground)",
                                      "Tree at 1.3 m\n(per sq.m bark)")))
pc <- ggplot(mon, aes(month, nmol, colour = src)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~src, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_colour_manual(values = setNames(c(SINK, SRC), levels(mon$src)), guide = "none") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(x = NULL, y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(size = 6.5),
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 7, color = "black"),
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
  # labels computed from the variables above -- never hardcode, they drift
  annotate("text", x = 1, y = soil_ann,             label = sprintf("%.0f", soil_ann),
           vjust = 1.5, size = 2.9) +
  annotate("text", x = 2, y = soil_ann + tree_meas, label = sprintf("+%.1f", tree_meas),
           vjust = -0.8, size = 2.9) +
  # total, not an increment: this bar INCLUDES the measured band drawn at x = 2
  annotate("text", x = 3, y = soil_ann + tree_scen, label = sprintf("%.0f total", tree_scen),
           vjust = -0.8, size = 2.9) +
  annotate("text", x = 4, y = soil_ann + tree_scen, label = sprintf("%.0f", soil_ann + tree_scen),
           vjust = 1.5, size = 2.9) +
  # bounding range on the scenario bar (WAI 1.50 - 3.07, constant flux)
  annotate("linerange", x = 3, ymin = soil_ann + tree_lo, ymax = soil_ann + tree_hi,
           colour = "grey25", linewidth = 0.5) +
  # range sits BELOW the scenario bar: at x = 3.46 it overprinted the net-budget bar
  annotate("text", x = 3, y = soil_ann + tree_lo,
           label = sprintf("WAI range %.0f-%.0f", tree_lo, tree_hi),
           size = 2.2, colour = "grey25", vjust = 1.9) +
  labs(x = NULL, y = expression("CH"[4]*" (mg m"^-2*" yr"^-1*", per m"^2*" ground)")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(size = 7.5), panel.grid.minor = element_blank())

# ---- assemble (net-sink / RF-R2 / basis notes go in the figure CAPTION) -------
fig <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & theme(plot.tag = element_text(face = "bold", size = 13))
ggsave(file.path(out, "fig_budget_maps.png"), fig, width = 10.5, height = 9, dpi = 300, bg = "white")
cat("Wrote fig_budget_maps.png\n")
cat(sprintf("Soil %.1f | tree measured %.2f (%.2f%% of soil) | tree scenario %.1f (%.1f%%, %.0f%% extrapolated)\n",
            soil_ann, tree_meas, off(tree_meas), tree_scen, off(tree_scen), pct_extrapolated))
cat(sprintf("Constant-flux WAI bounds %.1f-%.1f | full %d-combination range %.1f-%.1f\n",
            tree_lo, tree_hi, nrow(GR), tree_grid_lo, tree_grid_hi))
cat(sprintf("Net budget %.0f (measured only) to %.0f (scenario) mg CH4 m-2 yr-1 -- sink throughout\n",
            soil_ann+tree_meas, soil_ann+tree_scen))
