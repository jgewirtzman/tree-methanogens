#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_scaling_profiles.R
# ------------------------------------------------------------------------------
# What the upscaling assumptions ACTUALLY LOOK LIKE, height on the y axis.
#
# The grid says what each assumption does to a number. This says what each one IS,
# because "gaussian_75" and "linear_bounded_median" are labels a reader cannot
# picture, and the argument of the whole upscaling turns on these shapes agreeing
# inside the measured band and disagreeing above it.
#
# PANELS a-c ARE ONE TREE, NOT THE INVENTORY. An earlier version drew the stand
# profile -- the sum over ~8,000 stems of differing height -- and it was
# unreadable for a reason worth recording: trees drop out of successive height
# slabs as z passes their own height, so the curve acquires steps and bumps that
# belong to the DBH distribution rather than to the shape being illustrated.
# Worse, the 0-2 m band is the measured pi*DBH*2 while everything above is the
# WAI-scaled shape, so the stand curve carries a step at 2 m purely by
# construction. None of that tells a reader what "uniform_crown" means. These
# panels therefore use one canopy-height tree; only panel (d) sums the inventory,
# where the stand IS the point.
#
#   a  the six flux forms above 2 m AND the model's own measured profile below it.
#      The band is not one of the six: it is the RF's height response, a step
#      function breaking where the forest splits, rising to 1.9x the 2 m value at
#      the stem base. Showing it makes clear the forms are anchored to a measured
#      decline rather than to a guess.
#   b  the tree itself, bole and branch area stacked, for every bole x branch
#      combination. The two bole rows are near-identical, which is what "bole
#      shape does not matter" looks like.
#   c  the same tree with area filled by relative flux, so emission is visibly the
#      integral of width x fill.
#   d  the inventory sum, split at 2 m.
#
# Shapes are read from exports of rev_scaling_full_grid.R, never re-derived here.
#
# Output: outputs/revision/fig_scaling_profiles.png
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
outdir <- "outputs/revision"
FX <- read.csv(file.path(outdir, "scaling_flux_shapes.csv"),   stringsAsFactors = FALSE)
BP <- read.csv(file.path(outdir, "scaling_band_profile.csv"),  stringsAsFactors = FALSE)
KN <- read.csv(file.path(outdir, "scaling_shape_kernels.csv"), stringsAsFactors = FALSE)
AR <- read.csv(file.path(outdir, "scaling_area_profiles.csv"), stringsAsFactors = FALSE)
GR <- read.csv(file.path(outdir, "scaling_full_grid.csv"),     stringsAsFactors = FALSE)
source("code/revision/rev_geometry.R")

FO <- c("constant","rf_learned_capped","power","exponential","linear_floored",
        "linear_bounded_median")
BAND <- 2; CONV <- 86400*365.25*16e-6
H_TREE <- 25      # closed-canopy height at this site
D_TREE <- 0.556   # 95th-percentile DBH, the stem the height anchor is fitted to
REP_WAI <- grep("2.11", unique(AR$WAI), value = TRUE)[1]
REP_BRANCH <- "uniform_crown"; REP_BOLE <- "cone"; REP_FLUX <- "constant"

COLF <- setNames(c("#b2182b","#d6604d","#f4a582","#92c5de","#4393c3","#2166ac"), FO)
th <- theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.1, colour = "grey30"),
        panel.grid.minor = element_blank(), legend.title = element_blank(),
        legend.key.size = unit(0.32, "cm"), legend.text = element_text(size = 6.4))
bandrect <- annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = BAND,
                     fill = "grey85", alpha = 0.6)

# ---- a: the six forms above 2 m ----------------------------------------------
# rf_learned_capped is drawn as its area-weighted stand mean, and that is a
# summary of 8,006 DIFFERENT curves: the form fits a log-linear slope to each
# stem's own RF profile. It is not the same as constant -- median slope -0.607,
# so half the stems fall essentially to zero by 25 m -- but 2.3% have a slightly
# positive fitted slope and are held at 1.0 by the cap, which is what lifts the
# aggregate curve away from zero aloft. Every other form is one function of height
# applied to all stems.
FXa <- FX %>% filter(z >= BAND)
pa <- ggplot() +
  geom_vline(xintercept = 1, colour = "grey70", linewidth = 0.3, linetype = "dashed") +
  geom_path(data = FXa, aes(ratio, z, colour = flux), linewidth = 0.75) +
  scale_colour_manual(values = COLF, breaks = FO) +
  scale_y_continuous(limits = c(BAND, 26), expand = c(0, 0)) +
  labs(title = "a   the six extrapolation forms, above 2 m", y = "height (m)",
       x = "flux relative to that stem's 2 m value",
       subtitle = "rf_learned_capped is a per-stem slope shown as its stand mean; the rest are single functions") +
  th + theme(legend.position = c(0.70, 0.68),
             legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# ---- a2: the MEASURED band, per species and pooled ---------------------------
# The pooled curve is what the budget uses but it describes no individual species:
# the rise from 2 m down to the stem base runs 1.6x in Quercus rubra to 3.5x in
# Pinus strobus. Showing only the aggregate would present a stand-level average as
# if it were a biological profile. The SHAPE, though, is common to every species --
# a step up near the base and a shallow dip between 88 and 163 cm, which is the RF
# resolving three training heights and not a species-mixing artifact.
BPsp <- BP %>% filter(series != "all stems (area-weighted)")
BPag <- BP %>% filter(series == "all stems (area-weighted)")
pa2 <- ggplot() +
  geom_vline(xintercept = 1, colour = "grey70", linewidth = 0.3, linetype = "dashed") +
  geom_path(data = BPsp, aes(ratio, z, colour = series), linewidth = 0.6) +
  geom_path(data = BPag, aes(ratio, z), colour = "grey10", linewidth = 1.3) +
  annotate("text", x = max(BP$ratio)*0.97, y = 1.72, hjust = 1,
           label = "black = all stems, area-weighted", size = 2.3, colour = "grey10") +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(0, BAND), expand = c(0, 0)) +
  labs(title = "a2   the measured band, 0-2 m", y = "height (m)",
       x = "flux relative to that stem's 2 m value",
       subtitle = "the model's own height response, per species and pooled -- not one of the six forms") +
  th + theme(legend.position = "right", legend.text = element_text(size = 5.8),
             legend.key.size = unit(0.28, "cm"))

# ---- one canopy-height tree, built from the raw kernels ----------------------
NU <- length(unique(KN$u)); dz <- H_TREE/NU
CONIC <- pi * D_TREE * H_TREE / 2                    # W&W conic stem surface
tb <- KN %>% filter(part == "bole")   %>% select(u, bole = shape,   wb = w)
tr <- KN %>% filter(part == "branch") %>% select(u, branch = shape, wr = w)
TREE <- tb %>% inner_join(tr, by = "u", relationship = "many-to-many") %>%
  mutate(z = u*H_TREE,
         a_bole   = wb*(1/4.35)*CONIC/dz,     # W&W branch:stem 3.35 -> bole 1/4.35
         a_branch = wr*(3.35/4.35)*CONIC/dz) %>%
  pivot_longer(c(a_bole, a_branch), names_to = "part", values_to = "area_per_m") %>%
  mutate(part = factor(sub("^a_", "", part), levels = c("branch","bole")))

pb <- ggplot(TREE, aes(area_per_m, z, fill = part)) +
  bandrect +
  geom_area(orientation = "y", position = "stack", colour = NA) +
  scale_fill_manual(values = c(bole = "#8c6d46", branch = "#4d9221")) +
  facet_grid(bole ~ branch) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 1, 2)) +
  labs(title = sprintf("b   one %d m canopy tree: bole and branch area, stacked", H_TREE),
       subtitle = "rows = bole shape, columns = branch placement; the two rows are near-identical",
       x = expression("woody area (m"^2*" per m of height)"), y = "height (m)") +
  th + theme(legend.position = "bottom", strip.text = element_text(size = 6.4),
             axis.text = element_text(size = 5.8))

# ---- c: the same tree, filled by relative flux -------------------------------
# TWO scenarios, not one. Under constant flux the fill above 2 m is uniform by
# definition, so a single panel shows no gradient and the reader learns nothing
# from the colour. Pairing the reference case with a declining one makes the
# actual tension visible: woody area is concentrated high in the crown, exactly
# where every declining form says flux is lowest, so the two scenarios put their
# emission in completely different places on the same tree.
C_FORMS <- c("constant", "power")
fl_at <- function(z, form) {
  fx <- FX[FX$flux == form, ]
  ifelse(z <= BAND,
         approx(BPag$z, BPag$ratio, xout = pmin(z, BAND), rule = 2)$y,
         approx(fx$z, fx$ratio, xout = z, rule = 2)$y)
}
ONE <- TREE %>% filter(bole == REP_BOLE, branch == REP_BRANCH) %>%
  group_by(z) %>% summarise(area_per_m = sum(area_per_m), .groups = "drop")
ONE2 <- bind_rows(lapply(C_FORMS, function(f)
  ONE %>% mutate(form = factor(f, levels = C_FORMS), relflux = fl_at(z, f))))
pc <- ggplot(ONE2, aes(y = z)) +
  geom_ribbon(aes(xmin = 0, xmax = area_per_m, fill = relflux)) +
  scale_fill_gradientn(colours = c("#fff7f3","#fddbc7","#f4a582","#d6604d","#b2182b"),
                       name = "flux\n(rel. 2 m)") +
  geom_hline(yintercept = BAND, colour = "grey25", linewidth = 0.35, linetype = "dashed") +
  facet_wrap(~form, nrow = 1) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1, 1.5)) +
  labs(title = "c   one tree: width = area, fill = flux",
       subtitle = sprintf("%s bole, %s branch; emission is the integral of width x fill",
                          REP_BOLE, REP_BRANCH),
       x = expression("woody area (m"^2*" per m of height)"), y = "height (m)") +
  th + theme(legend.position = "right", legend.title = element_text(size = 6.2),
             strip.text = element_text(size = 7))

# ---- d: the inventory sum ----------------------------------------------------
P <- AR %>% filter(WAI == REP_WAI, bole == REP_BOLE, branch == REP_BRANCH) %>%
  select(z, area_m2_per_m, fluxarea_per_m) %>%
  left_join(FX %>% filter(flux == REP_FLUX) %>% select(z, ratio), "z") %>%
  mutate(ratio = ifelse(is.na(ratio), 1, ratio),
         integrand = fluxarea_per_m * ratio / STAND_AREA_M2 * CONV)
dzP <- diff(P$z)[1]; MEAS_BAND <- unique(GR$measured_mg)[1]
inb <- P$z <= BAND; rb <- sum(P$integrand[inb])*dzP
if (rb > 0) P$integrand[inb] <- P$integrand[inb]*(MEAS_BAND/rb)
tot <- sum(P$integrand)*dzP; band_tot <- sum(P$integrand[inb])*dzP
gv <- GR$total_mg[GR$WAI == REP_WAI & GR$bole == REP_BOLE &
                  GR$branch == REP_BRANCH & GR$flux == REP_FLUX]
pd <- ggplot(P, aes(y = z)) +
  geom_ribbon(aes(xmin = 0, xmax = integrand, fill = z <= BAND)) +
  scale_fill_manual(values = c(`TRUE` = "#2166ac", `FALSE` = "#f4a582"),
                    labels = c(`TRUE` = "measured (0-2 m)", `FALSE` = "extrapolated (>2 m)")) +
  geom_hline(yintercept = BAND, colour = "grey25", linewidth = 0.35, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  annotate("text", x = Inf, y = 24, hjust = 1.05, size = 2.9, colour = "grey20",
           label = sprintf("total %.0f mg CH4 m-2 yr-1\n%.0f%% above 2 m", tot,
                           100*(1 - band_tot/tot))) +
  labs(title = "d   summed over the whole inventory",
       subtitle = "steps are trees dropping out where height exceeds their own; the shaded area IS the total",
       x = expression("contribution (mg CH"[4]*" m"^-2*" yr"^-1*" per m of height)"),
       y = "height (m)") +
  th + theme(legend.position = c(0.64, 0.60),
             legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

fig <- (pa | pa2) / pb / (pc | pd) + plot_layout(heights = c(0.85, 1.05, 1)) +
  plot_annotation(
    title = "What the upscaling assumptions look like",
    subtitle = sprintf("(a)-(c) one canopy-height tree, so the shapes are legible; (d) the inventory sum from the %d-combination grid",
                       nrow(GR)),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.2, colour = "grey25")))
ggsave(file.path(outdir, "fig_scaling_profiles.png"), fig,
       width = 12, height = 9.5, dpi = 200, bg = "white")

cat(sprintf("one tree: H %d m, DBH %.3f m, conic stem surface %.1f m2\n", H_TREE, D_TREE, CONIC))
cat(sprintf("  area split bole %.2f : branch %.2f (W&W branch:stem 3.35)\n", 1/4.35, 3.35/4.35))
cat(sprintf("measured band, pooled: %.2fx the 2 m value at the stem base\n", BPag$ratio[1]))
cat("  per species at the base:\n")
print(as.data.frame(BPsp %>% group_by(series) %>% slice_min(z, n = 1) %>%
      transmute(species = series, base_ratio = round(ratio, 2)) %>% arrange(-base_ratio)),
      row.names = FALSE)
cat(sprintf("inventory sum, representative cell: %.1f mg (grid %.1f, residual %+.1f%% from binning)\n",
            tot, gv, 100*(tot-gv)/gv))
cat("\nwritten: outputs/revision/fig_scaling_profiles.png\n")
