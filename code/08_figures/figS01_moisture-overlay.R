# ==============================================================================
# figS01_moisture-overlay.R -- Figure S1: study plot, sampling design and soil moisture
# ------------------------------------------------------------------------------
# Stage E; run from the repository root. Drawn in UTM 18N metres (EPSG:26918).
#   (a) sampling design on the soil-moisture surface: inventory stems (canonical
#       inventory, located in-stand stems, sized by basal area), monthly-survey trees
#       and their plots, soil collars, the stream, the censused plot outline and the
#       moisture survey's extent; inset: location in the northeastern US
#   (b) inventory stems by species group (share of basal area)
#
# Basemap: hillshade and 2 m contours from the 2016 Connecticut statewide lidar (USGS
#   3DEP 1 m, data/raw/inventory/spatial_data/lidar_3dep_1m_utm18n.tif).
# Moisture: the analysis surface. 04_drivers/04_moisture_surface.R fits a thin-plate
#   spline to the December 2020 survey points and saves the fit; it is predicted here on
#   a 1 m grid (values outside the survey's convex hull, dashed, are extrapolated) and
#   bounded to the survey range.
# Stream: the 21 GPS flags in walked order. Lidar flow routing shows the defined
#   channel is the middle reach (flags 4-12: contributing area >4,000 m2), draining
#   west from a low point at flags 8-12 (~198 m); the southern and northern flags sit
#   on near-flat wet ground (contributing area <1,100 m2). The outlet is traced
#   downslope along lidar flow directions from flag 10. Only the defined channel and
#   outlet (1-2 m wide) are set to 100% VWC, after the fit, so the measured bank
#   readings are not altered.
# Monthly-survey plots: each plot's ellipse is the 95% covariance ellipse of its plot
#   points and the monthly-survey trees nearest to it (trees more than 2 SD from the
#   trees' median position are left out as GPS outliers).
#
# Inputs: data/raw/inventory/spatial_data/{plots.csv, River.xlsx,
#   YM_trees_measured.csv.csv, lidar_3dep_1m_utm18n.tif};
#   outputs/models/moisture_surface_tps.rds; outputs/tables/inventory_stems.csv
# Output: outputs/figures/original/supplementary/figS1_moisture_overlay.png
# ==============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readxl)})

# ---- monthly-survey plots, trees and the stream flags (raw GPS files) -----------
SPATIAL <- "data/raw/inventory/spatial_data"
plots_data <- read.csv(file.path(SPATIAL, "plots.csv"))
trees_data <- read.csv(file.path(SPATIAL, "YM_trees_measured.csv.csv"))
river_data <- read_excel(file.path(SPATIAL, "River.xlsx"))

# each monthly-survey tree to its nearest plot point, after dropping GPS outliers
# (more than 2 SD from the trees' median position)
assign_trees_to_plots <- function(trees_data, plots_data) {
  tc <- trees_data
  d <- sqrt((tc$Longitude - median(tc$Longitude, na.rm = TRUE))^2 +
            (tc$Latitude - median(tc$Latitude, na.rm = TRUE))^2)
  keep <- d <= mean(d, na.rm = TRUE) + 2 * sd(d, na.rm = TRUE)
  if (any(!keep)) {
    cat("Removing", sum(!keep), "outlier tree(s):\n")
    print(tc[!keep, c("Label", "Longitude", "Latitude")])
  }
  tc <- tc[keep, ]
  tc$Nearest_Plot <- NA; tc$Nearest_Site <- NA; tc$Distance_to_Plot <- NA
  for (i in seq_len(nrow(tc))) {
    dist <- sqrt((tc$Longitude[i] - plots_data$Longitude)^2 + (tc$Latitude[i] - plots_data$Latitude)^2)
    j <- which.min(dist)
    tc$Nearest_Plot[i] <- paste(plots_data$Site[j], plots_data$Plot[j], sep = "_")
    tc$Nearest_Site[i] <- plots_data$Site[j]
    tc$Distance_to_Plot[i] <- round(dist[j] * 111000, 2)   # metres
  }
  tc
}
trees_with_plots <- assign_trees_to_plots(trees_data, plots_data)

# 95% covariance ellipse of each plot's points and its assigned trees (>= 3 points)
create_plot_tree_ellipses <- function(plots_data, trees_with_plots) {
  plots_data$Site_Plot <- paste(plots_data$Site, plots_data$Plot, sep = "_")
  out <- list()
  for (site_plot in unique(plots_data$Site_Plot)) {
    ps <- plots_data[plots_data$Site_Plot == site_plot, ]
    ts <- trees_with_plots[trees_with_plots$Nearest_Plot == site_plot, ]
    pts <- rbind(data.frame(Longitude = ps$Longitude, Latitude = ps$Latitude),
                 if (nrow(ts) > 0) data.frame(Longitude = ts$Longitude, Latitude = ts$Latitude))
    if (nrow(pts) < 3) next
    e <- eigen(cov(cbind(pts$Longitude, pts$Latitude)))
    a <- seq(0, 2*pi, length.out = 100); k <- sqrt(qchisq(0.95, df = 2))
    out[[site_plot]] <- data.frame(
      Longitude = mean(pts$Longitude) + k * sqrt(e$values[1]) * cos(a) * e$vectors[1,1] + k * sqrt(e$values[2]) * sin(a) * e$vectors[1,2],
      Latitude  = mean(pts$Latitude)  + k * sqrt(e$values[1]) * cos(a) * e$vectors[2,1] + k * sqrt(e$values[2]) * sin(a) * e$vectors[2,2],
      Site_Plot = site_plot, Site = ps$Site[1], Plot = ps$Plot[1], n_points = nrow(ps), n_trees = nrow(ts))
  }
  do.call(rbind, out)
}
plot_tree_ellipses <- create_plot_tree_ellipses(plots_data, trees_with_plots)
cat(sprintf("monthly-survey trees %d on %d plot ellipses; stream flags %d\n",
            nrow(trees_with_plots), length(unique(plot_tree_ellipses$Site_Plot)), nrow(river_data)))

# ---- the figure ----------------------------------------------------------------
suppressPackageStartupMessages({library(fields); library(sf); library(terra); library(ggnewscale)
                                library(ggspatial); library(patchwork); library(maps)})
source("code/lib/geometry.R")
UTM <- 26918
MS <- readRDS("outputs/models/moisture_surface_tps.rds")
tr_ms <- geo_transforms()
to_utm <- function(lon, lat) st_coordinates(st_transform(st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326), UTM))

# stems: the canonical inventory the upscaling uses (in-stand, located)
INVc <- canonical_inventory(); INVc <- INVc[INVc$located, ]
llc <- tr_ms$fwd(INVc$PX, INVc$PY)
TR  <- as.data.frame(to_utm(llc$lon, llc$lat)); TR$BA <- pi * (INVc$dbh_m / 2)^2; TR$species <- INVc$species
cat(sprintf("map stems: %d located of %d in the canonical inventory\n", nrow(INVc), nrow(canonical_inventory())))
MT  <- as.data.frame(to_utm(trees_with_plots$Longitude, trees_with_plots$Latitude))
PL  <- as.data.frame(to_utm(plots_data$Longitude, plots_data$Latitude))
EL  <- cbind(as.data.frame(to_utm(plot_tree_ellipses$Longitude, plot_tree_ellipses$Latitude)), g = plot_tree_ellipses$Site_Plot)
RV  <- as.data.frame(to_utm(river_data$Longitude, river_data$Latitude))
SR  <- stand_ring_lonlat(); SRu <- as.data.frame(to_utm(c(SR$lon, SR$lon[1]), c(SR$lat, SR$lat[1])))
sh  <- MS$points[chull(MS$points$PX, MS$points$PY), ]; sh <- rbind(sh, sh[1, ])
SHu <- as.data.frame(to_utm(tr_ms$fwd(sh$PX, sh$PY)$lon, tr_ms$fwd(sh$PX, sh$PY)$lat))

allx <- c(TR$X, RV$X, EL$X, SRu$X); ally <- c(TR$Y, RV$Y, EL$Y, SRu$Y)
E <- c(floor(min(allx)) - 25, ceiling(max(allx)) + 15, floor(min(ally)) - 15, ceiling(max(ally)) + 15)

dem <- rast("data/raw/inventory/spatial_data/lidar_3dep_1m_utm18n.tif"); crs(dem) <- paste0("EPSG:", UTM)
dem[dem < -100] <- NA
demc <- crop(dem, ext(E[1] - 50, E[2] + 50, E[3] - 50, E[4] + 50))
hs  <- shade(terrain(demc, "slope", unit = "radians"), terrain(demc, "aspect", unit = "radians"), 40, 315)
HS  <- as.data.frame(crop(hs, ext(E)), xy = TRUE); names(HS)[3] <- "hs"
demS <- focal(demc, w = 5, fun = "mean", na.rm = TRUE)
CN  <- st_as_sf(as.contour(crop(demS, ext(E)), levels = seq(180, 320, by = 2)))

# moisture on a 1 m grid
G <- expand.grid(X = seq(E[1] + 0.5, E[2] - 0.5, by = 1), Y = seq(E[3] + 0.5, E[4] - 0.5, by = 1))
ll <- st_coordinates(st_transform(st_as_sf(G, coords = c("X", "Y"), crs = UTM), 4326))
pl <- tr_ms$inv(ll[, 2], ll[, 1])
G$VWC <- pmin(pmax(as.numeric(predict(MS$fit, cbind(pl$PX, pl$PY))), MS$floor), max(MS$points$vwc))

# stream: defined channel = flags 4-12; outlet traced downslope from flag 10 along D8 directions
fd <- terrain(demS, v = "flowdir")
step <- list(`1` = c(1, 0), `2` = c(1, -1), `4` = c(0, -1), `8` = c(-1, -1), `16` = c(-1, 0), `32` = c(-1, 1), `64` = c(0, 1), `128` = c(1, 1))
p <- as.numeric(RV[10, c("X", "Y")]); OUT <- matrix(p, 1)
for (k in 1:600) {
  d <- extract(fd, matrix(p, 1))[[1]]
  if (is.na(d) || !(as.character(d) %in% names(step))) break
  p <- p + step[[as.character(d)]] * res(fd)[1]
  if (p[1] < E[1] || p[1] > E[2] || p[2] < E[3] || p[2] > E[4]) break
  OUT <- rbind(OUT, p)
}
OUT <- as.data.frame(OUT); names(OUT) <- c("X", "Y")
cat(sprintf("stream outlet traced %d m from flag 10\n", nrow(OUT) - 1))
# The defined channel is drawn as one smooth course through the main-stem flags
# (4-8, then the low point at 10) and on along the traced outlet; the flags clustered
# at the confluence (9, 11, 12) and the flags on the flat ground north and south are
# shown as points. Smoothing removes the 1 m D8 stair-steps and a few metres of GPS jitter.
P0 <- rbind(RV[c(4, 5, 6, 7, 8, 10), c("X", "Y")], OUT[-1, ])
tt <- c(0, cumsum(sqrt(diff(P0$X)^2 + diff(P0$Y)^2)))
dfk <- max(4, round(nrow(P0) / 4))
tq <- seq(0, max(tt), by = 0.5)
CH <- data.frame(X = predict(smooth.spline(tt, P0$X, df = dfk), tq)$y, Y = predict(smooth.spline(tt, P0$Y, df = dfk), tq)$y)
seg_d <- function(px, py, ax, ay, bx, by) { vx <- bx - ax; vy <- by - ay; t <- pmin(1, pmax(0, ((px - ax) * vx + (py - ay) * vy) / (vx^2 + vy^2 + 1e-12))); sqrt((px - ax - t * vx)^2 + (py - ay - t * vy)^2) }
chan_d <- function(P, L) { d <- rep(Inf, nrow(P)); for (j in seq_len(nrow(L) - 1)) d <- pmin(d, seg_d(P$X, P$Y, L$X[j], L$Y[j], L$X[j + 1], L$Y[j + 1])); d }
STREAM_HALF_WIDTH_M <- 0.75   # a 1-2 m brook
G$VWC[chan_d(G, CH) <= STREAM_HALF_WIDTH_M] <- 100
cat(sprintf("Moisture surface (04_moisture_surface.R fit) over the map: %d cells, VWC %.1f-%.1f%%\n", nrow(G), min(G$VWC), max(G$VWC)))

# Moisture is shown over the rectangle spanning every tree and collar drawn (inventory
# stems, monthly-survey trees, soil collars), aligned with the plot grid (plot-local
# PX/PY); beyond the dashed survey hull those values are extrapolated. Outside the
# rectangle only the terrain is drawn.
mt_p <- tr_ms$inv(trees_with_plots$Latitude, trees_with_plots$Longitude)
pl_p <- tr_ms$inv(plots_data$Latitude, plots_data$Longitude)
RX <- range(c(INVc$PX, mt_p$PX, pl_p$PX), na.rm = TRUE) + c(-5, 5); RY <- range(c(INVc$PY, mt_p$PY, pl_p$PY), na.rm = TRUE) + c(-5, 5)   # 5 m margin
in_rect <- pl$PX >= RX[1] & pl$PX <= RX[2] & pl$PY >= RY[1] & pl$PY <= RY[2]
G <- G[in_rect, ]

# shared layers
# hillshade as a fixed grey image, so it takes no fill scale
hsc <- crop(hs, ext(E)); hm <- as.matrix(hsc, wide = TRUE)
hm <- (hm - min(hm, na.rm = TRUE)) / diff(range(hm, na.rm = TRUE)); hm[is.na(hm)] <- 1
HSIMG <- matrix(grDevices::grey(0.25 + 0.75 * hm), nrow(hm)); he <- as.vector(ext(hsc))
terrain_layers <- list(annotation_raster(HSIMG, xmin = he[1], xmax = he[2], ymin = he[3], ymax = he[4]))
HSIMG_L <- matrix(grDevices::grey(0.55 + 0.45 * hm), nrow(hm))
terrain_light <- list(annotation_raster(HSIMG_L, xmin = he[1], xmax = he[2], ymin = he[3], ymax = he[4]))
map_frame <- list(
  coord_sf(crs = UTM, xlim = E[1:2], ylim = E[3:4], expand = FALSE, datum = NA),
  theme_void(),
  theme(legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9),
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(4, 4, 4, 4)))
CONT <- geom_sf(data = CN, colour = "grey30", linewidth = 0.15, alpha = 0.5)
PT_BREAKS <- c("Inventory stem", "Monthly-survey tree", "Soil collar", "Stream GPS flag")
LN_BREAKS <- c("Censused plot", "Moisture survey extent", "Monthly-survey plot", "Stream channel")

# (a) sampling design and soil moisture
pa <- ggplot() + terrain_layers +
  geom_raster(data = G, aes(X, Y, fill = VWC), alpha = 0.6) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Soil moisture\n(VWC, %)", limits = c(0, 100)) +
  CONT +
  geom_path(data = SRu, aes(X, Y), colour = "black", linewidth = 0.5) +
  geom_path(data = SHu, aes(X, Y), colour = "grey10", linewidth = 0.5, linetype = "22") +
  geom_path(data = CH, aes(X, Y), colour = "#2166ac", linewidth = 1.1, lineend = "round") +
  geom_point(data = TR, aes(X, Y, size = BA), colour = "grey15", alpha = 0.55, stroke = 0) +
  scale_size_area(max_size = 3, name = "Basal area (m²)", breaks = c(0.01, 0.05, 0.1, 0.2)) +
  geom_path(data = EL, aes(X, Y, group = g), colour = "#d95f02", linewidth = 0.7) +
  geom_point(data = MT, aes(X, Y), shape = 8, colour = "black", size = 2, stroke = 0.6) +
  geom_point(data = PL, aes(X, Y), shape = 21, fill = "white", colour = "black", size = 2.2, stroke = 0.7) +
  geom_point(data = RV, aes(X, Y), shape = 21, fill = "#6baed6", colour = "#08519c", size = 1.4, stroke = 0.4) +
  # legend-only layers, placed off the map: one key per point type and per line type
  geom_point(data = data.frame(X = E[1] - 1e4, Y = E[3] - 1e4, k = factor(PT_BREAKS, PT_BREAKS)), aes(X, Y, shape = k)) +
  geom_path(data = data.frame(X = E[1] - 1e4 + rep(0:1, 4), Y = E[3] - 1e4, k = factor(rep(LN_BREAKS, each = 2), LN_BREAKS)),
            aes(X, Y, group = k, linetype = k)) +
  scale_shape_manual(name = NULL, breaks = PT_BREAKS, values = c(16, 8, 21, 21),
    guide = guide_legend(order = 3, override.aes = list(shape = c(16, 8, 21, 21), colour = c("grey15", "black", "black", "#08519c"),
                                                        fill = c(NA, NA, "white", "#6baed6"), size = c(2, 2, 2.2, 1.8), alpha = 1, stroke = c(0, 0.6, 0.7, 0.4)))) +
  scale_linetype_manual(name = NULL, breaks = LN_BREAKS, values = c("solid", "22", "solid", "solid"),
    guide = guide_legend(order = 4, override.aes = list(colour = c("black", "grey10", "#d95f02", "#2166ac"), linewidth = c(0.5, 0.5, 0.7, 1.1)))) +
  guides(fill = guide_colourbar(order = 1), size = guide_legend(order = 2)) +
  annotation_scale(location = "bl", width_hint = 0.2, style = "ticks", line_col = "black", text_col = "black") +
  annotation_north_arrow(location = "tr", height = unit(0.8, "cm"), width = unit(0.6, "cm"), style = north_arrow_minimal()) +
  map_frame

# (b) species: eight groups, each its own colour, in Paul Tol's muted scheme plus two
# darks, assigned by tree: hemlock green, white pine pale cyan, oaks near-black, red maple
# wine, sugar maple blue, black birch indigo, yellow birch sand, mountain laurel rose.
# Checked by simulation (OKLab dE x100, Machado 2009): normal 17.8, deutan 9.0, protan
# 9.5 (tritan 5.8, sugar maple vs black birch; tritanopia is rare). The remaining species
# are hollow, since a ninth colour fails the normal-vision floor; they are listed with
# their shares in the caption. Legend labels carry each group's share of basal area,
# computed from the full canonical inventory.
SP_GROUPS <- c("Tsuga canadensis", "Pinus strobus", "Quercus spp.", "Acer rubrum", "Acer saccharum", "Betula lenta",
               "Betula alleghaniensis", "Kalmia latifolia", "Other")
SP_COLS <- c("#117733", "#88CCEE", "#222222", "#882255", "#0072B2", "#332288", "#DDCC77", "#CC6677", "white")
grp_of <- function(sp) factor(ifelse(grepl("^Quercus", sp), "Quercus spp.", ifelse(sp %in% SP_GROUPS, sp, "Other")), levels = SP_GROUPS)
TR$grp <- grp_of(TR$species)
INVall <- canonical_inventory(); BAall <- pi * (INVall$dbh_m / 2)^2
PCT <- tapply(BAall, grp_of(INVall$species), sum) / sum(BAall) * 100
fmt <- function(x) ifelse(x < 1, sprintf("%.1f%%", x), sprintf("%.0f%%", x))
SP_LABS <- lapply(SP_GROUPS, function(g) {
  pc <- fmt(PCT[[g]])
  if (g == "Quercus spp.") bquote(italic("Quercus")~"spp. ("*.(pc)*")")
  else if (g == "Other") bquote("Other species ("*.(pc)*")")
  else bquote(italic(.(g))~"("*.(pc)*")") })
cat("basal-area shares (%):\n"); print(round(PCT, 2))
TRb <- TR[order(TR$grp == "Other", TR$grp == "Kalmia latifolia", TR$BA, decreasing = c(TRUE, TRUE, FALSE), method = "radix"), ]
pb <- ggplot() + terrain_light + CONT +
  geom_path(data = SRu, aes(X, Y), colour = "black", linewidth = 0.5) +
  geom_path(data = CH, aes(X, Y), colour = "#0d366b", linewidth = 1.1, lineend = "round") +
  geom_point(data = TRb, aes(X, Y, size = BA, fill = grp), shape = 21, colour = "grey20", stroke = 0.12, alpha = 0.9) +
  scale_fill_manual(values = setNames(SP_COLS, SP_GROUPS), breaks = SP_GROUPS, name = "Species (share of basal area)", drop = FALSE,
    labels = SP_LABS,
    guide = guide_legend(override.aes = list(size = 3.2, alpha = 1, stroke = 0.4))) +
  scale_size(range = c(0.45, 3.4), guide = "none") +
  annotation_scale(location = "bl", width_hint = 0.2, style = "ticks", line_col = "black", text_col = "black") +
  map_frame

# location inset: northeastern US, Connecticut shaded, site marked
st <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
ne <- st[st$ID %in% c("connecticut", "massachusetts", "rhode island", "new york", "new jersey", "pennsylvania",
                      "vermont", "new hampshire", "maine"), ]
site <- st_as_sf(data.frame(lon = mean(river_data$Longitude), lat = mean(river_data$Latitude)), coords = c("lon", "lat"), crs = 4326)
inset <- ggplot() +
  geom_sf(data = ne, fill = "grey92", colour = "grey55", linewidth = 0.2) +
  geom_sf(data = ne[ne$ID == "connecticut", ], fill = "grey70", colour = "grey40", linewidth = 0.3) +
  geom_sf(data = site, shape = 21, fill = "red", colour = "black", size = 2.2, stroke = 0.5) +
  coord_sf(crs = 5070, expand = FALSE) + theme_void() +
  theme(panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.4), plot.background = element_rect(fill = "white", colour = NA))

# layout: each map with its legends in a right-hand column; the inset sits above (a)'s legends
legA <- cowplot::get_legend(pa); legB <- cowplot::get_legend(pb)
pa0 <- pa + theme(legend.position = "none"); pb0 <- pb + theme(legend.position = "none")
colA <- cowplot::plot_grid(inset, legA, ncol = 1, rel_heights = c(0.42, 1))
rowA <- cowplot::plot_grid(pa0, colA, nrow = 1, rel_widths = c(1, 0.28))
rowB <- cowplot::plot_grid(pb0, legB, nrow = 1, rel_widths = c(1, 0.28))
final_extended_plot <- cowplot::plot_grid(rowA, rowB, ncol = 1, labels = c("(a)", "(b)"), label_size = 14)

ggsave("outputs/figures/original/supplementary/figS1_moisture_overlay.png", final_extended_plot, width = 11, height = 15.5, dpi = 300, bg = "white")
