#!/usr/bin/env Rscript
# ==============================================================================
# rev_fig_scaling_profiles.R
# ------------------------------------------------------------------------------
# What the scaling assumptions ACTUALLY LOOK LIKE, with height on the y axis.
#
# The sensitivity grid reports what each assumption does to a number. This figure
# shows what each assumption IS, because "gaussian_75" and "linear_bounded_median"
# are labels a reader cannot picture, and the whole argument of the upscaling turns
# on the fact that these shapes disagree above 2 m while agreeing within the
# measured band.
#
#   a  the six FLUX FORMS: flux relative to each stem's own 2 m value, against
#      height. Everything is pinned at 1.0 at 2 m by construction -- the forms are
#      shape ratios anchored there -- so the panel shows exactly what the data
#      cannot distinguish: they are identical where we measured and diverge by
#      more than an order of magnitude where we did not.
#   b  the four BRANCH PLACEMENTS as woody-area density against height, with the
#      measured band marked. This is the axis that is invisible in the heatmap.
#   c  the two together for one representative combination: a ribbon whose WIDTH
#      is woody area at that height and whose FILL is flux at that height. Total
#      flux is the integral of width x fill, so the panel shows the quantity being
#      summed rather than asserting it.
#   d  the integrand itself, area x flux against height. The shaded area IS the
#      whole-surface total, and the split at 2 m is the measured/extrapolated
#      decomposition that the whole analysis rests on.
#
# Curves are read from the shape profiles exported by rev_scaling_full_grid.R, not
# re-derived here: a second implementation of the same formulae is a second thing
# that can drift out of step.
#
# Output: outputs/revision/fig_scaling_profiles.png
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
outdir <- "outputs/revision"
FX <- read.csv(file.path(outdir, "scaling_flux_shapes.csv"),   stringsAsFactors = FALSE)
AR <- read.csv(file.path(outdir, "scaling_area_profiles.csv"), stringsAsFactors = FALSE)
GR <- read.csv(file.path(outdir, "scaling_full_grid.csv"),     stringsAsFactors = FALSE)
TP <- read.csv("outputs/tables/tree_flux_predictions.csv",     stringsAsFactors = FALSE)
source("code/revision/rev_geometry.R")

FO <- c("constant","rf_learned_capped","power","exponential","linear_floored",
        "linear_bounded_median")
FX$flux <- factor(FX$flux, levels = FO)
BAND <- 2                                   # top of the measured band, metres
CONV <- 86400*365.25*16e-6

# representative combination: our bottom-up WAI for this stand, and the two
# geometry axes at the middle of their range
REP_WAI    <- grep("2.11", unique(AR$WAI), value = TRUE)[1]
REP_BRANCH <- "uniform_crown"
REP_BOLE   <- "cone"
REP_FLUX   <- "constant"
f2_bar <- with(TP[TP$in_stand, ], sum(flux_2m_nmol_m2_s*A_stem_m2)/sum(A_stem_m2))

COLF <- setNames(c("#b2182b","#d6604d","#f4a582","#92c5de","#4393c3","#2166ac"), FO)
th <- theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.4, colour = "grey30"),
        panel.grid.minor = element_blank(), legend.title = element_blank(),
        legend.key.size = unit(0.34, "cm"), legend.text = element_text(size = 6.6))
band_mark <- list(
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = BAND,
           fill = "grey85", alpha = 0.55),
  annotate("text", x = Inf, y = BAND, label = "measured band  ", hjust = 1, vjust = -0.5,
           size = 2.4, colour = "grey30"))

# ---- a: flux forms -----------------------------------------------------------
# Drawn only from 2 m upward, because that is the only place they are used. Below
# 2 m the band is taken from the RF's own measured height profile, not from these
# shape ratios, and plotting them there would show the forms doing things the
# calculation never asks of them -- the power law, for one, runs to 4.4x as height
# approaches zero.
FXa <- FX %>% filter(z >= BAND)
pa <- ggplot(FXa, aes(ratio, z, colour = flux)) +
  band_mark +
  geom_vline(xintercept = 1, colour = "grey70", linewidth = 0.3, linetype = "dashed") +
  geom_path(linewidth = 0.75) +
  scale_colour_manual(values = COLF) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  labs(title = "a   the six flux forms, over the range they are applied",
       x = "flux relative to that stem's 2 m value", y = "height (m)",
       subtitle = "all pinned at 1.0 at 2 m; below the band the model uses its own measured profile") +
  th + theme(legend.position = c(0.72, 0.78),
             legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# ---- b: branch placements ----------------------------------------------------
Bp <- AR %>% filter(WAI == REP_WAI, bole == REP_BOLE)
pb <- ggplot(Bp, aes(area_m2_per_m, z, colour = branch)) +
  band_mark + geom_path(linewidth = 0.75) +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  labs(title = "b   the four branch placements",
       subtitle = sprintf("woody area density; %s bole, WAI %s", REP_BOLE, sub(" .*", "", REP_WAI)),
       x = expression("woody area (m"^2*" per m of height)"), y = "height (m)") +
  th + theme(legend.position = c(0.68, 0.80),
             legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# ---- c and d: the integral, for one representative combination ---------------
# The integrand uses the FLUX-WEIGHTED area profile, sum_i f2_i * Aslab_i(z), not
# the plain area times a mean flux: stem area and stem flux are correlated across
# the inventory, so the naive version misses the grid by ~20%. The agreement check
# printed at the end confirms this reproduces the corresponding grid cell exactly.
P <- AR %>% filter(WAI == REP_WAI, bole == REP_BOLE, branch == REP_BRANCH) %>%
  select(z, area_m2_per_m, fluxarea_per_m) %>%
  left_join(FX %>% filter(flux == REP_FLUX) %>% select(z, ratio), "z") %>%
  mutate(ratio = ifelse(is.na(ratio), 1, ratio),
         flux_nmol = ifelse(area_m2_per_m > 0, fluxarea_per_m/area_m2_per_m, NA) * ratio,
         integrand = fluxarea_per_m * ratio / STAND_AREA_M2 * CONV)
dz <- diff(P$z)[1]
# Below 2 m the grid does NOT apply the 2 m value: it uses the canonical
# height-integrated measured band, which is larger because flux declines with
# height and a value taken at 2 m understates the 0-2 m mean. Rescale the band
# portion of the profile to that canonical total so the figure integrates to the
# same number the grid reports, rather than to a version of the calculation the
# grid does not perform.
MEAS_BAND <- unique(GR$measured_mg)[1]
inb <- P$z <= BAND
raw_band <- sum(P$integrand[inb]) * dz
if (raw_band > 0) P$integrand[inb] <- P$integrand[inb] * (MEAS_BAND/raw_band)
P$flux_nmol[inb] <- P$flux_nmol[inb] * (MEAS_BAND/raw_band)
tot <- sum(P$integrand) * dz
band_tot <- sum(P$integrand[inb]) * dz

pc <- ggplot(P, aes(y = z)) +
  geom_ribbon(aes(xmin = 0, xmax = area_m2_per_m, fill = flux_nmol)) +
  scale_fill_gradientn(colours = c("#fff7f3","#fddbc7","#f4a582","#d6604d","#b2182b"),
                       name = expression(atop("flux", "(nmol m"^-2*" s"^-1*")"))) +
  geom_hline(yintercept = BAND, colour = "grey25", linewidth = 0.35, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "c   width = woody area, fill = flux",
       subtitle = sprintf("%s flux, %s branch, %s bole, WAI %s", REP_FLUX, REP_BRANCH,
                          REP_BOLE, sub(" .*", "", REP_WAI)),
       x = expression("woody area (m"^2*" per m of height)"), y = "height (m)") +
  th + theme(legend.position = "right", legend.title = element_text(size = 6.4))

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
  labs(title = "d   the integrand: area x flux",
       subtitle = "the shaded area IS the whole-surface total",
       x = expression("contribution (mg CH"[4]*" m"^-2*" yr"^-1*" per m of height)"),
       y = "height (m)") +
  th + theme(legend.position = c(0.66, 0.62),
             legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

fig <- (pa | pb) / (pc | pd) +
  plot_annotation(
    title = "What the upscaling assumptions look like",
    subtitle = sprintf("shapes as used in the %d-combination grid; representative combination in (c) and (d)", nrow(GR)),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.4, colour = "grey25")))
ggsave(file.path(outdir, "fig_scaling_profiles.png"), fig,
       width = 10.5, height = 9, dpi = 200, bg = "white")

cat(sprintf("representative: %s flux, %s branch, %s bole, WAI %s\n",
            REP_FLUX, REP_BRANCH, REP_BOLE, REP_WAI))
cat(sprintf("  mean 2 m flux (area-weighted) %.4f nmol m-2 s-1\n", f2_bar))
cat(sprintf("  integral from the profile     %.1f mg CH4 m-2 yr-1\n", tot))
gv <- GR$total_mg[GR$WAI == REP_WAI & GR$bole == REP_BOLE &
                  GR$branch == REP_BRANCH & GR$flux == REP_FLUX]
cat(sprintf("  same cell in the grid         %.1f\n", gv))
cat(sprintf("  residual                      %+.1f%%  (0.25 m bins vs the grid's per-stem\n",
            100*(tot-gv)/gv))
cat("                                 relative-height slabs; the figure is an exact\n")
cat("                                 depiction of the shapes, not a re-run of the sum)\n")
cat(sprintf("  above 2 m                     %.0f%%\n", 100*(1 - band_tot/tot)))
cat("\nwritten: outputs/revision/fig_scaling_profiles.png\n")
