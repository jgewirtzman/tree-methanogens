#!/usr/bin/env Rscript
# ==============================================================================
# figS27_scaling-profiles.R -- Figure S28: the scaling assumptions, height on the y axis
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) stem CH4 flux with height (absolute, nmol m-2 s-1): the modelled 0-2 m profile
#       (black, area-weighted over all stems; grey, per species; shaded, the measured
#       band) and the six forms of flux above 2 m
#   (b) bole and branch woody area with height for ONE canopy-height tree (H =
#       CANOPY_H_M, DBH = the 95th percentile over stems > 10 cm) under each branch
#       placement, with the named pair's bole shape. One tree, not the stand: summed
#       over stems of differing height the profile carries steps from the DBH
#       distribution, and a step at 2 m by construction.
#   (c) emission per m of height summed over the inventory, for the named pair
#       (scaling_headline.csv: decline and hold), measured band vs extrapolated; the
#       printed totals are the grid's (the drawn 1 m binned integral approximates them)
#   (d) leave-one-height-out cross-validation: RMSE of each form
# Shapes are read from the exports of 06_upscale/07_scaling_full_grid.R, never re-derived.
# Inputs (outputs/data/): scaling_flux_shapes.csv, scaling_band_profile.csv,
#   scaling_shape_kernels.csv, scaling_area_profiles.csv, scaling_full_grid.csv,
#   scaling_headline.csv, scaling_flux_anchors.csv (07_scaling_full_grid.R),
#   height_form_crossvalidation.csv (05_model/13_height_form_crossvalidation.R);
#   code/lib/geometry.R (inventory, canopy height, branch:stem area split)
# Output: outputs/figures/generated/fig_scaling_profiles.png
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
source("code/lib/outputs.R")
FX <- read.csv(out_path("scaling_flux_shapes.csv"),   stringsAsFactors = FALSE)
BP <- read.csv(out_path("scaling_band_profile.csv"),  stringsAsFactors = FALSE)
KN <- read.csv(out_path("scaling_shape_kernels.csv"), stringsAsFactors = FALSE)
AR <- read.csv(out_path("scaling_area_profiles.csv"), stringsAsFactors = FALSE)
GR <- read.csv(out_path("scaling_full_grid.csv"),     stringsAsFactors = FALSE)
source("code/lib/geometry.R")

FO <- c("constant","exp_band_slope","power","exponential","linear_floored",
        "linear_bounded_median")
stopifnot(setequal(FO, unique(FX$flux)), setequal(FO, unique(GR$flux)))
BAND <- 2; CONV <- 86400*365.25*16e-6
# flux purple, as stem flux in Figures 2 and 7
FLUX_PURPLE <- "#756BB1"; FLUX_PALE <- "#f2eff7"
H_TREE <- CANOPY_H_M   # closed-canopy height; defined once in geometry.R
# 95th-percentile DBH, the stem the height allometry is anchored to (geometry.R)
D_TREE <- as.numeric(local({ i <- canonical_inventory(); quantile(i$dbh_m[i$dbh_m > 0.10], .95, na.rm = TRUE) }))
# the representative WAI, selected by its label (meaning), not by its digits
REP_WAI <- grep("bottom-up, this stand", unique(AR$WAI), value = TRUE)[1]
stopifnot(length(REP_WAI) == 1, !is.na(REP_WAI))
# geometry and the two flux forms of the named pair, read from the producer
HL <- read.csv(out_path("scaling_headline.csv"), stringsAsFactors = FALSE)
stopifnot(identical(HL$case, c("decline", "hold")), length(unique(HL$branch)) == 1,
          length(unique(HL$bole)) == 1, all(HL$WAI == REP_WAI))
REP_BRANCH <- HL$branch[1]; REP_BOLE <- HL$bole[1]

COLF <- setNames(c("#b2182b","#d6604d","#f4a582","#92c5de","#4393c3","#2166ac"), FO)
th <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.title = element_blank(),
        legend.key.size = unit(0.45, "cm"), legend.text = element_text(size = 9.5))
FLUX_LAB <- c(constant = "Constant", exp_band_slope = "Per-stem exponential", power = "Power",
              exponential = "Exponential", linear_floored = "Linear, floored at zero",
              linear_bounded_median = "Linear into uptake")
BRANCH_LAB <- c(uniform_all = "Uniform,\nwhole stem", uniform_top50 = "Uniform,\nupper half",
                gaussian_50 = "Gaussian,\n50% height", gaussian_75 = "Gaussian,\n75% height")
bandrect <- annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = BAND,
                     fill = "grey85", alpha = 0.6)

# ---- a: band and extrapolations together, one continuous 0-26 m panel --------
# Absolute flux: the forms agree where flux was measured and diverge where it was not,
# and the values are of order 0.05-0.14 nmol m-2 s-1.
AN   <- read.csv(out_path("scaling_flux_anchors.csv"))
BPag <- BP %>% filter(series == "all stems (area-weighted)")
BPsp <- BP %>% filter(series != "all stems (area-weighted)")
Bb <- BPag %>% transmute(z, v = absolute)
Ff <- FX %>% filter(z >= BAND) %>% transmute(flux, z, v = ratio*AN$f2_stand_mean)
Ss <- BPsp %>% transmute(series, z, v = absolute)

pa <- ggplot() + bandrect +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_path(data = Ss, aes(v, z, group = series), colour = "grey72", linewidth = 0.35) +
  geom_path(data = Ff, aes(v, z, colour = flux), linewidth = 0.9) +
  geom_path(data = Bb, aes(v, z), colour = "grey10", linewidth = 1.3) +
  scale_colour_manual(values = COLF, breaks = FO, labels = FLUX_LAB) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  labs(y = "Height (m)", x = expression(CH[4]~"flux (nmol m"^-2*" s"^-1*")")) +
  th + theme(legend.position = "right")

# ---- one canopy-height tree, built from the raw kernels ----------------------
NU <- length(unique(KN$u)); dz <- H_TREE/NU
CONIC <- pi * D_TREE * H_TREE / 2                    # W&W conic stem surface
tb <- KN %>% filter(part == "bole")   %>% select(u, bole = shape,   wb = w)
tr <- KN %>% filter(part == "branch") %>% select(u, branch = shape, wr = w)
TREE <- tb %>% inner_join(tr, by = "u", relationship = "many-to-many") %>%
  mutate(z = u*H_TREE,
         a_bole   = wb*BOLE_FRACTION*CONIC/dz,     # W&W branch:stem area split (geometry.R)
         a_branch = wr*BRANCH_FRACTION*CONIC/dz) %>%
  pivot_longer(c(a_bole, a_branch), names_to = "part", values_to = "area_per_m") %>%
  mutate(part = factor(sub("^a_", "", part), levels = c("branch","bole")))

pb <- ggplot(TREE %>% filter(bole == REP_BOLE), aes(area_per_m, z, fill = part)) +
  bandrect +
  geom_area(orientation = "y", position = "stack", colour = NA) +
  scale_fill_manual(values = c(bole = "#8c6d46", branch = "#4d9221")) +
  facet_wrap(~branch, nrow = 1, labeller = as_labeller(BRANCH_LAB)) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 1, 2)) +
  labs(x = expression("Woody area (m"^2*" per m of height)"), y = "Height (m)") +
  th + theme(legend.position = "bottom", strip.text = element_text(size = 9.5),
             axis.text = element_text(size = 9))

# the named pair (decline = exponential, hold = constant)
C_FORMS <- HL$flux
CASE_LAB <- setNames(sprintf("%s: flux %s above 2 m", FLUX_LAB[HL$flux],
                             ifelse(HL$case == "decline", "declines", "holds")), HL$flux)
fx_of <- function(form) FX[FX$flux == form, ]

# ---- c: the inventory sum, for the named pair ------------------------------
# at the bottom-up WAI and the named pair's branch placement and bole shape
MEAS_BAND <- unique(GR$measured_mg)[1]
mkP <- function(form) {
  AR %>% filter(WAI == REP_WAI, bole == REP_BOLE, branch == REP_BRANCH) %>%
    select(z, fluxarea_per_m) %>%
    left_join(fx_of(form) %>% select(z, ratio), "z") %>%
    # the forms apply above 2 m only; the 0-2 m band is the measured profile for every form
    mutate(ratio = ifelse(is.na(ratio) | z <= BAND, 1, ratio),
           integrand = fluxarea_per_m * ratio / STAND_AREA_M2 * CONV,
           form = factor(form, levels = C_FORMS))
}
PD <- bind_rows(lapply(C_FORMS, mkP))
dzP <- diff(sort(unique(PD$z)))[1]
PD <- PD %>% group_by(form) %>%
  mutate(integrand = ifelse(z <= BAND,
           integrand * MEAS_BAND/sum(integrand[z <= BAND]*dzP), integrand)) %>% ungroup()
LAB <- PD %>% group_by(form) %>%
  summarise(tot = sum(integrand)*dzP,
            above = 100*(1 - sum(integrand[z <= BAND])*dzP/sum(integrand*dzP)), .groups="drop")
# The printed totals are the grid's (scaling_full_grid.csv), the numbers the text quotes;
# the binned integral drawn here is a 1 m approximation of the same sum.
LAB_BIN <- LAB
LAB <- LAB %>% mutate(form_chr = as.character(form)) %>%
  left_join(GR %>% filter(WAI == REP_WAI, branch == REP_BRANCH, bole == REP_BOLE) %>%
              transmute(form_chr = flux, g_tot = total_mg, g_above = pct_extrapolated), "form_chr") %>%
  mutate(tot = g_tot, above = g_above)
stopifnot(all(is.finite(LAB$tot)))
pd <- ggplot(PD, aes(y = z)) +
  geom_ribbon(aes(xmin = 0, xmax = integrand, fill = z <= BAND)) +
  # purple = measured, grey = extrapolated
  scale_fill_manual(values = c(`TRUE` = FLUX_PURPLE, `FALSE` = "grey72"),
                    labels = c(`TRUE` = "measured (0-2 m)", `FALSE` = "extrapolated (>2 m)")) +
  geom_hline(yintercept = BAND, colour = "grey25", linewidth = 0.3, linetype = "dashed") +
  geom_text(data = LAB, aes(x = Inf, y = 24, label = sprintf("%.1f mg in total\n%.0f%% above 2 m", tot, above)),
            hjust = 1.05, size = 3.4, colour = "grey20", inherit.aes = FALSE) +
  facet_wrap(~form, nrow = 1, labeller = as_labeller(CASE_LAB)) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = expression("mg CH"[4]*" m"^-2*" yr"^-1*" per m of height"), y = "Height (m)") +
  th + theme(legend.position = "bottom", strip.text = element_text(size = 9.5))

# ---- d: leave-one-height-out cross-validation of the forms ---------------------
CV <- read.csv(out_path("height_form_crossvalidation.csv"), stringsAsFactors = FALSE) %>%
  mutate(test = factor(test, levels = c("A_up", "C_int", "B_down"),
                       labels = c("Upward\n(0.5, 1.25 m predict 2 m)", "Interpolation\n(0.5, 2 m predict 1.25 m)",
                                  "Downward\n(1.25, 2 m predict 0.5 m)")),
         form = recode(form, constant = "Constant", linear = "Linear", power = "Power", exp_zero = "Exponential"))
CV$form <- factor(CV$form, levels = rev(c("Constant", "Linear", "Power", "Exponential")))
pdcv <- ggplot(CV, aes(rmse, form)) +
  geom_segment(aes(x = min(CV$rmse) * 0.7, xend = rmse, yend = form), colour = "grey75") +
  geom_point(size = 2.6, colour = FLUX_PURPLE) +
  scale_x_log10() + facet_wrap(~test, ncol = 1) +
  labs(x = expression("RMSE (nmol m"^-2*" s"^-1*", log scale)"), y = NULL) +
  th + theme(strip.text = element_text(size = 9.5))

# panels: (a) flux forms, (b) woody area for one canopy tree, (c) inventory sum, (d) cross-validation
fig <- (pa | pb) / (pd | pdcv) + plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(out_path("fig_scaling_profiles.png"), fig, width = 15, height = 14, dpi = 300, bg = "white")

cat(sprintf("one tree: H %g m, DBH %.3f m, conic stem surface %.1f m2\n", H_TREE, D_TREE, CONIC))
cat(sprintf("  area split bole %.2f : branch %.2f (W&W branch:stem %.2f)\n", BOLE_FRACTION, BRANCH_FRACTION, BRANCH_STEM_RATIO))
cat(sprintf("measured band, pooled: %.2fx the 2 m value at the stem base\n", BPag$ratio[1]))
cat("  per species at the base:\n")
print(as.data.frame(BPsp %>% group_by(series) %>% slice_min(z, n = 1) %>%
      transmute(species = series, base_ratio = round(ratio, 2)) %>% arrange(-base_ratio)),
      row.names = FALSE)
cat(sprintf("\n=== NAMED PAIR drawn in (c): WAI %s | branch %s | bole %s ===\n",
            REP_WAI, REP_BRANCH, REP_BOLE))
cat(sprintf("  %-8s %-12s grid total %.1f mg CH4 m-2 yr-1 (%.1f%% of soil, %.0f%% extrapolated)\n",
            HL$case, HL$flux, HL$total_mg, HL$pct_of_soil, HL$pct_extrapolated), sep = "")
cat("\n  panel (d) labels (grid) against this figure's 1 m binned integral:\n")
print(as.data.frame(LAB %>% transmute(form, grid_mg = round(tot,1), grid_pct_above = round(above),
                                      binned_mg = round(LAB_BIN$tot,1))), row.names = FALSE)
cat("\nwritten: outputs/figures/generated/fig_scaling_profiles.png\n")
