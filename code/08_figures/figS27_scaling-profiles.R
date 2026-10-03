source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# figS27_scaling-profiles.R
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
# Shapes are read from exports of 07_scaling_full_grid.R, never re-derived here.
#
# Output: outputs/figures/generated/fig_scaling_profiles.png
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
outdir <- "outputs"
FX <- read.csv(out_path("scaling_flux_shapes.csv"),   stringsAsFactors = FALSE)
BP <- read.csv(out_path("scaling_band_profile.csv"),  stringsAsFactors = FALSE)
KN <- read.csv(out_path("scaling_shape_kernels.csv"), stringsAsFactors = FALSE)
AR <- read.csv(out_path("scaling_area_profiles.csv"), stringsAsFactors = FALSE)
GR <- read.csv(out_path("scaling_full_grid.csv"),     stringsAsFactors = FALSE)
source("code/lib/geometry.R")

FO <- c("constant","exp_band_slope","power","exponential","linear_floored",
        "linear_bounded_median")
BAND <- 2; CONV <- 86400*365.25*16e-6
# FLUX PURPLE, the colour this project already uses for CH4 flux (rev_fig02,
# rev_fig02a, rev_fig07). The integral panels previously used red and blue, which
# encoded nothing -- the two flux forms are distinguished by their facet labels,
# not by hue, so an arbitrary two-colour scheme just invited a reader to look for
# a meaning that was not there.
FLUX_PURPLE <- "#756BB1"; FLUX_PALE <- "#f2eff7"
H_TREE <- CANOPY_H_M   # closed-canopy height; defined once in geometry.R
# 95th-percentile DBH, the stem the height anchor is fitted to. Was hardcoded 0.556;
# derived here so it cannot drift from the inventory the anchor is actually fitted on.
D_TREE <- as.numeric(local({ i <- canonical_inventory(); quantile(i$dbh_m[i$dbh_m > 0.10], .95, na.rm = TRUE) }))
# Select the representative WAI by MEANING, not by digits. This was
# grep("2.11", ...), which is the same digit-matching failure that
# 07_scaling_full_grid.R:339 documents as fixed on the PRODUCER side -- but this
# CONSUMER was not updated, so once the bottom-up estimate was corrected to 2.82 the
# grep returned NA, panel (d) filtered to zero rows and the subtitle printed "NA".
REP_WAI <- grep("bottom-up, this stand", unique(AR$WAI), value = TRUE)[1]
stopifnot(length(REP_WAI) == 1, !is.na(REP_WAI))
REP_BRANCH <- "gaussian_75"; REP_BOLE <- "cone"; REP_FLUX <- "exp_band_slope"

COLF <- setNames(c("#b2182b","#d6604d","#f4a582","#92c5de","#4393c3","#2166ac"), FO)
th <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.title = element_blank(),
        legend.key.size = unit(0.45, "cm"), legend.text = element_text(size = 9.5))
FLUX_LAB <- c(constant = "Constant", exp_band_slope = "Per-stem exponential (named)", power = "Power",
              exponential = "Exponential", linear_floored = "Linear, floored at zero",
              linear_bounded_median = "Linear into uptake")
BRANCH_LAB <- c(uniform_all = "Uniform,\nwhole stem", uniform_top50 = "Uniform,\nupper half",
                gaussian_50 = "Gaussian,\n50% height", gaussian_75 = "Gaussian,\n75% height")
bandrect <- annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = BAND,
                     fill = "grey85", alpha = 0.6)

# ---- a: ONE panel, band and extrapolations together, in three normalisations --
# Reverted to a single continuous 0-26 m panel: splitting the band onto its own
# axis broke the one thing the figure exists to show, which is that the forms
# agree where we measured and diverge where we did not.
#
# Three normalisations side by side, because the choice is not obvious and each
# answers a different question. Relative to 2 m matches how the forms are defined
# and anchored. Relative to 1.25 m centres on the middle measured height instead
# of the top one, so the band's own decline is not compressed against the anchor.
# Absolute is the only version in which a reader can see that all of this concerns
# fluxes of order 0.05-0.1 nmol m-2 s-1.
AN   <- read.csv(out_path("scaling_flux_anchors.csv"))
BPag <- BP %>% filter(series == "all stems (area-weighted)")
BPsp <- BP %>% filter(series != "all stems (area-weighted)")
FXa  <- FX %>% filter(z >= BAND)
SC125 <- AN$f2_stand_mean / AN$f125_stand_mean

mk <- function(kind) {
  if (kind == "relative to 2 m") {
    list(b  = BPag %>% transmute(z, v = ratio),
         f  = FXa  %>% transmute(flux, z, v = ratio),
         sp = BPsp %>% transmute(series, z, v = ratio))
  } else if (kind == "relative to 1.25 m") {
    list(b  = BPag %>% transmute(z, v = rel125),
         f  = FXa  %>% transmute(flux, z, v = ratio*SC125),
         sp = BPsp %>% transmute(series, z, v = rel125))
  } else {
    list(b  = BPag %>% transmute(z, v = absolute),
         f  = FXa  %>% transmute(flux, z, v = ratio*AN$f2_stand_mean),
         sp = BPsp %>% transmute(series, z, v = absolute))
  }
}
# relative-to-1.25 m was tested and dropped: it shifts the crossing point without
# adding information, and compresses the extrapolations. Relative to 2 m is how the
# forms are anchored; absolute is the only view showing these are fluxes of order
# 0.05-0.14 nmol m-2 s-1, which every ratio view hides.
KINDS <- c("absolute (nmol m-2 s-1)")   # 2026-10-02: one panel, absolute (the relative view is in the code history)
L  <- setNames(lapply(KINDS, mk), KINDS)
gr <- function(w) bind_rows(lapply(KINDS, function(k) L[[k]][[w]] %>% mutate(kind = k))) %>%
        mutate(kind = factor(kind, levels = KINDS))
Bb <- gr("b"); Ff <- gr("f"); Ss <- gr("sp")

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
         a_bole   = wb*(1/4.35)*CONIC/dz,     # W&W branch:stem 3.35 -> bole 1/4.35
         a_branch = wr*(3.35/4.35)*CONIC/dz) %>%
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

# ---- c: one tree, area filled by flux AND the resulting integral ------------
# Two views per flux form. The filled ribbon shows WHERE the surface area is and
# how hot it is; the integral shows what that actually contributes per metre of
# height. They are different questions and the second is the one that ends up in
# the budget: a crown can hold most of the area and still contribute little if the
# flux there is assumed low, which is exactly what separates the two forms.
C_FORMS <- c("constant", "exp_band_slope")
fx_of <- function(form) FX[FX$flux == form, ]
fl_at <- function(z, form) {
  fx <- fx_of(form)
  ifelse(z <= BAND,
         approx(BPag$z, BPag$ratio, xout = pmin(z, BAND), rule = 2)$y,
         approx(fx$z, fx$ratio, xout = z, rule = 2)$y)
}
ONE <- TREE %>% filter(bole == REP_BOLE, branch == REP_BRANCH) %>%
  group_by(z) %>% summarise(area_per_m = sum(area_per_m), .groups = "drop")
SEC_YR <- 86400*365.25; NMOL_MG <- 16e-6
ONE2 <- bind_rows(lapply(C_FORMS, function(f)
  ONE %>% mutate(form = factor(f, levels = C_FORMS),
                 relflux = fl_at(z, f),
                 absflux = relflux * AN$f2_stand_mean,
                 contrib = area_per_m * absflux * SEC_YR * NMOL_MG)))
pc1 <- ggplot(ONE2, aes(y = z)) +
  geom_ribbon(aes(xmin = 0, xmax = area_per_m, fill = relflux)) +
  # rates keep the red ramp: here the colour encodes MAGNITUDE, which a sequential
  # warm scale reads more directly than purple. Purple is reserved for the
  # integrals, where it identifies the quantity rather than scaling with it.
  scale_fill_gradientn(colours = c("#fff7f3","#fddbc7","#f4a582","#d6604d","#b2182b"),
                       name = "flux\n(rel. 2 m)") +
  geom_hline(yintercept = BAND, colour = "grey25", linewidth = 0.3, linetype = "dashed") +
  facet_wrap(~form, nrow = 1) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 1, 2)) +
  labs(title = "c   one tree: where the area is, and how hot",
       subtitle = sprintf("%s bole, %s branch", REP_BOLE, REP_BRANCH),
       x = expression("woody area (m"^2*" per m)"), y = "height (m)") +
  th + theme(legend.position = "bottom", legend.title = element_text(size = 6),
             strip.text = element_text(size = 7))
# lines, not filled ribbons: the fill in (c) already encodes area, so filling the
# integral too reads as a second area when it is a rate per metre of height
pc2 <- ggplot(ONE2, aes(y = z)) +
  geom_path(aes(x = contrib, colour = form), linewidth = 0.9) +
  scale_colour_manual(values = setNames(rep(FLUX_PURPLE, length(C_FORMS)), C_FORMS)) +
  geom_hline(yintercept = BAND, colour = "grey25", linewidth = 0.3, linetype = "dashed") +
  facet_wrap(~form, nrow = 1) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "     the integral: what each metre contributes",
       subtitle = "one tree, mg CH4 per year from each metre of height",
       x = expression("mg CH"[4]*" yr"^-1*" per m"), y = NULL) +
  th + theme(legend.position = "none", strip.text = element_text(size = 7))

# ---- d: the inventory sum, for the same two forms ---------------------------
# Which scenario is being drawn was previously implicit. It is now read from
# scaling_headline.csv, defined once in 07_scaling_full_grid.R: the conservative
# flux form, our own bottom-up WAI for this stand, a crown centred at 0.75H, and
# the bole shape that makes no difference. Both forms are shown, because the RANGE
# is the result and a single bar would read as an estimate.
HL <- read.csv(out_path("scaling_headline.csv"), stringsAsFactors = FALSE)
MEAS_BAND <- unique(GR$measured_mg)[1]
mkP <- function(form) {
  AR %>% filter(WAI == REP_WAI, bole == REP_BOLE, branch == REP_BRANCH) %>%
    select(z, fluxarea_per_m) %>%
    left_join(fx_of(form) %>% select(z, ratio), "z") %>%
    mutate(ratio = ifelse(is.na(ratio), 1, ratio),
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
# The printed totals are the GRID's (scaling_full_grid.csv), the numbers the text quotes.
# The binned integral drawn here is a 1 m approximation of the same sum and read
# 13.9 mg against the grid's 15.4 for the named cell; the shape is the figure's job,
# the number is the grid's.
GRD <- read.csv("outputs/data/scaling_full_grid.csv", stringsAsFactors = FALSE)
LAB_BIN <- LAB
LAB <- LAB %>% mutate(form_chr = as.character(form)) %>%
  left_join(GRD %>% filter(WAI == HL$WAI, branch == HL$branch, bole == HL$bole) %>%
              transmute(form_chr = flux, g_tot = total_mg, g_above = pct_extrapolated), "form_chr") %>%
  mutate(tot = g_tot, above = g_above)
stopifnot(all(is.finite(LAB$tot)))
pd <- ggplot(PD, aes(y = z)) +
  geom_ribbon(aes(xmin = 0, xmax = integrand, fill = z <= BAND)) +
  # purple = measured, grey = extrapolated. The colour now carries the one
  # distinction that matters in this panel -- what was observed against what was
  # assumed -- instead of an unrelated warm/cool pairing.
  scale_fill_manual(values = c(`TRUE` = FLUX_PURPLE, `FALSE` = "grey72"),
                    labels = c(`TRUE` = "measured (0-2 m)", `FALSE` = "extrapolated (>2 m)")) +
  geom_hline(yintercept = BAND, colour = "grey25", linewidth = 0.3, linetype = "dashed") +
  geom_text(data = LAB, aes(x = Inf, y = 24, label = sprintf("%.1f mg in total\n%.0f%% above 2 m", tot, above)),
            hjust = 1.05, size = 3.4, colour = "grey20", inherit.aes = FALSE) +
  facet_wrap(~form, nrow = 1, labeller = as_labeller(FLUX_LAB)) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = expression("mg CH"[4]*" m"^-2*" yr"^-1*" per m of height"), y = "Height (m)") +
  th + theme(legend.position = "bottom", strip.text = element_text(size = 9.5))

# ---- d: why a constant above 2 m: leave-one-height-out cross-validation ---------
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
cat(sprintf("  area split bole %.2f : branch %.2f (W&W branch:stem 3.35)\n", 1/4.35, 3.35/4.35))
cat(sprintf("measured band, pooled: %.2fx the 2 m value at the stem base\n", BPag$ratio[1]))
cat("  per species at the base:\n")
print(as.data.frame(BPsp %>% group_by(series) %>% slice_min(z, n = 1) %>%
      transmute(species = series, base_ratio = round(ratio, 2)) %>% arrange(-base_ratio)),
      row.names = FALSE)
cat("\n=== HEADLINE SCENARIO drawn in (c) and (d) ===\n")
cat(sprintf("  flux %s | WAI %s | branch %s | bole %s\n",
            HL$flux, HL$WAI, HL$branch, HL$bole))
cat(sprintf("  grid total %.1f mg CH4 m-2 yr-1 (%.1f%% of soil, %.0f%% extrapolated)\n",
            HL$total_mg, HL$pct_of_soil, HL$pct_extrapolated))
cat("\n  panel (d) labels (grid) against this figure's 1 m binned integral:\n")
print(as.data.frame(LAB %>% transmute(form, grid_mg = round(tot,1), grid_pct_above = round(above),
                                      binned_mg = round(LAB_BIN$tot,1))), row.names = FALSE)
cat("\nwritten: outputs/figures/generated/fig_scaling_profiles.png\n")
