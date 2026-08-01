#!/usr/bin/env Rscript
# ==============================================================================
# rev_surface_area_model.R
# ------------------------------------------------------------------------------
# An explicit, documented woody-surface-area model for the stand upscaling, and the
# vertical distribution scenarios that go with it.
#
# ------------------------------------------------------------------------------
# THE AREA MODEL
# ------------------------------------------------------------------------------
# Total woody surface area is NOT computed from our inventory -- it is set by the
# woody area index of Gauci et al. 2024 (Nature 631:796-800), the empirical value for
# global temperate forests:
#
#       A_total = WAI x A_plot,   WAI = 3.07 m2 woody surface per m2 GROUND
#
# We cannot measure WAI here: it requires terrestrial LiDAR / destructive branch
# harvest, neither of which was done at this site. WAI is therefore an imported
# stand-level constant, and the budget inherits whatever error it carries.
#
# What our inventory DOES give is the bole component, from DBH and a taper assumption:
#
#       conic taper:  d(z) = DBH * (H - z)/(H - 1.3)      (DBH measured at 1.3 m)
#       bole area  :  A_bole = integral_0^H pi*d(z) dz = pi*DBH*(H - 1.3)/2 + base
#
# Branch area is then the residual, A_branch = A_total - sum(A_bole). This is a
# DEFINITION, not a measurement: it absorbs all error in WAI, in the taper, and in the
# assumed canopy height H.
#
# ------------------------------------------------------------------------------
# VERTICAL DISTRIBUTION SCENARIOS  (area density a(z), normalised to integrate to 1)
# ------------------------------------------------------------------------------
#   uniform      a(z) = const on [0,H]
#                The simplest assumption and the one behind the manuscript's headline
#                scaling. Not physically motivated; it is the "no information" case.
#
#   conic        a(z) proportional to (H - z)
#                Bole only, main stem tapering linearly to a point at H. This is a real
#                geometric model, but it describes ONLY the bole -- it is the right
#                shape if branch area is negligible, which the reconciliation below
#                shows it is not.
#
#   bole+crown   a(z) = conic bole  +  branch area spread UNIFORMLY between the crown
#                base (z_cb = CROWN_BASE_FRAC * H) and the top.
#                This is the composite the area decomposition actually implies: a
#                tapering stem carrying most of the length, plus a much larger branch
#                area confined to the crown. Branch area is placed uniformly through
#                the crown rather than peaked, because nothing in our data or in the
#                imported WAI constrains its shape within the crown -- a peaked
#                (e.g. Gaussian) profile would be an invention.
#
# NOTE: an earlier version of this analysis used a Gaussian crown profile centred at
# 0.75H. That shape had no empirical basis and has been removed.
#
# ------------------------------------------------------------------------------
# WHEN THE DISTRIBUTION MATTERS
# ------------------------------------------------------------------------------
# If flux is assumed CONSTANT above 2 m, the whole-surface budget depends only on TOTAL
# area above 2 m; a(z) cancels. The vertical distribution enters only when flux is
# allowed to vary with height. This is reported explicitly below.
#
# Output: outputs/revision/surface_area_model.txt
#         outputs/revision/fig_surface_area_distribution.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(ranger)})
source("code/lib/rev_geometry.R")
source("code/lib/rev_species_levels.R")
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load("outputs/models/RF_MODELS.RData"); load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
INV <- canonical_inventory()

# CANONICAL stems and ground (2026-07-29). This read
# rf_workflow_data$PLACEHOLDER_INVENTORY over A_PLOT <- 200*200: fg19 only, so the
# entire bytag survey (961 stems) was missing -- 7,010 stems instead of 8,006 -- over
# the nominal 40,000 m2 square instead of the 37,200 m2 censused stand. Every
# per-ground figure below was 7.0% low and every WAI*area 7.5% high, and the stem-area
# index came out 0.0888 against the canonical 0.11629. The allometry and the DBH repair
# were also local copies; both now come from rev_geometry.R.
A_PLOT <- STAND_AREA_M2; WAI <- 3.07; CONV <- 86400*365.25*16e-6
CANOPY_H <- 20            # m, assumed; varied in the sensitivity grid
CROWN_BASE_FRAC <- 0.5    # crown base at half canopy height
rule <- function(s) cat("\n", strrep("=",78), "\n ", s, "\n", strrep("=",78), "\n", sep="")

# fix_dbh() and the local band definition removed: canonical_inventory() already
# returns repaired diameters plus band_m (0.75 m for Kalmia, a shrub; 2 m otherwise)
# and A_band_m2. Re-repairing repaired data with species == "..." rather than %in%
# turned the 41 NA-species stems into NA diameters, which is why this reported 7,998
# stems and an NA band area.
sum_piD <- sum(pi*INV$dbh)
A_band  <- sum(INV$A_band_m2)

# ==============================================================================
rule("AREA DECOMPOSITION")
# ==============================================================================
A_total <- WAI*A_PLOT
cat(sprintf("
  Imported total (WAI %.2f x %d m2)         %8.0f m2   index %.3f
  Measured band, 0-2 m (cylinder)           %8.0f m2   index %.3f   (%.1f%% of total)
", WAI, A_PLOT, A_total, WAI, A_band, A_band/A_PLOT, 100*A_band/A_total))

cat("\n  Bole / branch split under a conic taper, by assumed canopy height:\n\n")
cat(sprintf("    %6s %12s %12s %12s %10s\n","H (m)","bole (m2)","bole index","branch idx","branch %"))
for (H in c(15,20,25)) {
  Ab <- sum(pi*INV$dbh*(H-1.3)/2)
  cat(sprintf("    %6d %12.0f %12.3f %12.3f %9.0f%%\n",
              H, Ab, Ab/A_PLOT, (A_total-Ab)/A_PLOT, 100*(A_total-Ab)/A_total))
}
cat("
  Under a conic taper the bole accounts for only ~11-15% of the imported WAI, so the
  branch residual dominates. This is a consequence of the imported WAI, not an
  independent finding: any error in WAI lands almost entirely in the branch term.
")

# ==============================================================================
rule("VERTICAL AREA DENSITY a(z)")
# ==============================================================================
H <- CANOPY_H; zcb <- CROWN_BASE_FRAC*H
z  <- seq(0.01, H, length.out = 600)
A_bole_tot <- sum(pi*INV$dbh*(H-1.3)/2); A_branch_tot <- A_total - A_bole_tot

dens <- list(
  uniform    = rep(1, length(z)),
  conic      = pmax(H - z, 0),
  `bole+crown` = { b <- pmax(H - z, 0); b <- b/sum(b)*A_bole_tot
                   c_ <- as.numeric(z >= zcb); c_ <- c_/sum(c_)*A_branch_tot
                   b + c_ })
dens <- lapply(dens, function(v) v/sum(v))

cat(sprintf("
  Canopy height %d m, crown base %.0f m (%.0f%% of H).
  Fraction of total woody area lying ABOVE 2 m, and its area-weighted mean height:

", H, zcb, 100*CROWN_BASE_FRAC))
cat(sprintf("    %-14s %14s %20s\n","scenario","frac above 2m","mean height (m)"))
for (nm in names(dens))
  cat(sprintf("    %-14s %14.3f %20.1f\n", nm,
              sum(dens[[nm]][z>2]), sum(dens[[nm]]*z)))

# ==============================================================================
rule("DOES THE DISTRIBUTION MATTER? Depends entirely on the flux assumption.")
# ==============================================================================
pair <- tree_train_complete %>% filter(!is.na(measurement_height_cm), is.finite(stem_flux_corrected))
fb <- sapply(c(50,125,200), function(hh) mean(pair$stem_flux_corrected[pair$measurement_height_cm==hh]))
b_decl <- coef(lm(log(pmax(fb,1e-9)) ~ I(c(50,125,200)/100)))[2]

flux_forms <- list(
  `constant above 2 m`      = function(zz) rep(1, length(zz)),
  `log-linear decline`      = function(zz) exp(b_decl*(zz-2)),
  `decline to 50% asymptote`= function(zz) 0.5 + 0.5*exp(b_decl*(zz-2)))
cat(sprintf("\n    %-26s %12s %12s %12s\n","flux assumption above 2m", names(dens)[1], names(dens)[2], names(dens)[3]))
for (fn in names(flux_forms)) {
  vals <- sapply(names(dens), function(dn) {
    k <- z > 2; w <- dens[[dn]][k]/sum(dens[[dn]][k]); sum(w*flux_forms[[fn]](z[k])) })
  cat(sprintf("    %-26s %12.4f %12.4f %12.4f\n", fn, vals[1], vals[2], vals[3]))
}
cat("
  Read the first row: under a CONSTANT flux above 2 m every distribution gives exactly
  1.0000 -- the vertical distribution cancels completely. It only matters once flux is
  allowed to vary with height, and then it matters a great deal. This is the argument
  for making the constant-flux case the headline value: it is the one number in this
  analysis that is INSENSITIVE to the assumption we are least able to defend.
")

# ==============================================================================
rule("RF HEIGHT RESPONSE ACROSS SOIL MOISTURE (the transport mechanism)")
# ==============================================================================
d <- tree_train_complete; jm <- d[d$month==7,]
qs <- quantile(d$soil_moisture_at_tree, c(.1,.5,.9), na.rm=TRUE)
g <- expand.grid(species = levels(factor(d$species_clean)), height_cm = c(50,125,200),
                 soil_moisture_at_tree = qs, stringsAsFactors = FALSE)
g$dbh_m <- median(d$dbh_m,na.rm=TRUE)
g$soil_temp_C_mean <- mean(jm$soil_temp_C_mean,na.rm=TRUE)
g$air_temp_C_mean  <- mean(jm$air_temp_C_mean, na.rm=TRUE)
g$species <- factor(g$species, levels=levels(factor(d$species_clean)))
g$pred <- predict(TreeRF, g, num.threads=1)$predictions
mo <- g %>% group_by(soil_moisture_at_tree, height_cm) %>%
  summarise(f=mean(pred), .groups="drop") %>%
  pivot_wider(names_from=height_cm, values_from=f, names_prefix="h") %>%
  mutate(ratio_200_50 = h200/h50)
cat(sprintf("\n  Stand-mean RF prediction, by soil moisture percentile:\n\n"))
cat(sprintf("    %-12s %10s %10s %10s %14s\n","moisture","50 cm","125 cm","200 cm","ratio 200/50"))
for (i in seq_len(nrow(mo)))
  cat(sprintf("    %-12s %10.4f %10.4f %10.4f %14.3f\n",
      sprintf("%.3f (%s)", mo$soil_moisture_at_tree[i], c("p10","p50","p90")[i]),
      mo$h50[i], mo$h125[i], mo$h200[i], mo$ratio_200_50[i]))
cat("
  If the height decline is a transport phenomenon driven by soil water, the RF should
  show a STEEPER decline (smaller ratio) at high moisture. Whether it does is reported
  above; because moisture and species are confounded in the height campaign, this is
  suggestive of the mechanism rather than a clean test of it.
")

# --- figure -------------------------------------------------------------------
pd <- bind_rows(lapply(names(dens), function(nm)
  data.frame(z=z, a=dens[[nm]]/max(dens[[nm]]), scenario=nm)))
p <- ggplot(pd, aes(a, z, colour=scenario)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=2, fill="grey88") +
  geom_path(linewidth=.8) +
  annotate("text", x=.97, y=1, label="measured band (0-2 m)", hjust=1, size=2.6, colour="grey30") +
  scale_colour_manual(values=c(uniform="#999999", conic="#2166AC", `bole+crown`="#B2182B")) +
  labs(title="Assumed vertical distribution of woody surface area",
       subtitle=sprintf("canopy %d m, crown base %.0f m; WAI = %.2f imported from Gauci et al. 2024 (no LiDAR at this site)", H, zcb, WAI),
       x="relative area density a(z)", y="height (m)") +
  theme_bw(base_size=9)
ggsave(file.path(outdir,"fig_surface_area_distribution.png"), p, width=7, height=5, dpi=200, bg="white")
cat("\nWritten: outputs/revision/fig_surface_area_distribution.png\n")
