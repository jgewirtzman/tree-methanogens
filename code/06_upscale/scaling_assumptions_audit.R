source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# scaling_assumptions_audit.R
# ------------------------------------------------------------------------------
# What, exactly, is the stand upscaling assuming? Four questions, answered from the
# data where the data can answer them and flagged as unconstrained where it cannot.
#
#   Q1  What is the height-flux profile, overall and by species?
#   Q2  What allometry converts DBH to surface area?
#   Q3  How is woody surface area distributed vertically (bole vs branches)?
#   Q4  How does our own surface area relate to the Gauci et al. 2024 WAI of 3.07?
#
# Output: outputs/audit/scaling_assumptions_audit.txt (transcript of this run)
#         outputs/figures/generated/fig_height_profiles_by_species.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(ranger)})
source("code/lib/geometry.R")
source("code/lib/species_levels.R")
set.seed(42)
outdir <- "outputs"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")
load("data/processed/integrated/rf_workflow_input_data_with_2023.RData")
INV <- canonical_inventory()

# CANONICAL stems and ground (2026-07-29). This read
# rf_workflow_data$PLACEHOLDER_INVENTORY over A_PLOT <- 200*200: fg19 only, so the
# entire bytag survey (961 stems) was missing -- 7,010 stems instead of 8,006 -- over
# the nominal 40,000 m2 square instead of the 37,200 m2 censused stand. Every
# per-ground figure below was 7.0% low and every WAI*area 7.5% high, and the stem-area
# index came out 0.0888 against the canonical 0.11629. The allometry and the DBH repair
# were also local copies; both now come from geometry.R.
A_PLOT <- STAND_AREA_M2; WAI_GAUCI <- 3.07; CONV <- 86400*365.25*16e-6

rule <- function(s) cat("\n", strrep("=", 78), "\n ", s, "\n", strrep("=", 78), "\n", sep="")

# ==============================================================================
rule("Q1  HEIGHT-FLUX PROFILE: what the data actually contain")
# ==============================================================================
d  <- tree_train_complete
hh <- d[!is.na(d$measurement_height_cm), ]
cat(sprintf("
  Rows with a measured height : %d of %d training rows (%.0f%%)
  The other %d rows are assigned the 125 cm default, so the height effect in the RF
  is identified from the %d multi-height rows ONLY.
  All multi-height data come from month %s (a single campaign): heights %s cm.
  Trees measured at all three heights: %d
",
  nrow(hh), nrow(d), 100*nrow(hh)/nrow(d), nrow(d)-nrow(hh), nrow(hh),
  paste(sort(unique(hh$month)), collapse=","),
  paste(sort(unique(hh$measurement_height_cm)), collapse="/"),
  sum(table(hh$tree_id) == 3)))

# --- paired within-tree profile: the cleanest view, removes between-tree variance
pair <- hh %>% dplyr::select(tree_id, species_clean, h = measurement_height_cm,
                             f = stem_flux_corrected) %>%
  filter(is.finite(f)) %>% group_by(tree_id, species_clean, h) %>%
  summarise(f = mean(f), .groups="drop")

cat("\n  --- OVERALL profile (all trees pooled) ---\n")
ov <- pair %>% group_by(h) %>% summarise(mean=mean(f), median=median(f), n=n(), .groups="drop")
print(as.data.frame(ov), row.names=FALSE, digits=3)

wide <- pair %>% pivot_wider(names_from=h, values_from=f, names_prefix="h") %>%
  filter(is.finite(h50), is.finite(h125), is.finite(h200))
cat(sprintf("\n  Paired trees with all 3 heights: %d\n", nrow(wide)))
cat(sprintf("  Within-tree change 50->200 cm: mean %+.3f nmol m-2 s-1, median %+.3f\n",
            mean(wide$h200-wide$h50), median(wide$h200-wide$h50)))
cat(sprintf("  Trees DECLINING with height: %d/%d (%.0f%%)  [sign test p = %.4g]\n",
            sum(wide$h200 < wide$h50), nrow(wide),
            100*mean(wide$h200 < wide$h50),
            binom.test(sum(wide$h200 < wide$h50), nrow(wide))$p.value))
cat(sprintf("  Wilcoxon signed-rank (50 vs 200 cm): p = %.4g\n",
            wilcox.test(wide$h50, wide$h200, paired=TRUE)$p.value))

cat("\n  --- BY SPECIES (paired, within-tree) ---\n")
sp <- wide %>% group_by(species_clean) %>%
  summarise(n=n(), f50=mean(h50), f125=mean(h125), f200=mean(h200),
            d_50_200 = mean(h200-h50), pct_declining = 100*mean(h200<h50), .groups="drop") %>%
  arrange(d_50_200)
print(as.data.frame(sp), row.names=FALSE, digits=3)
cat("
  NOTE: n per species is 1-19. Only ~6 species have n >= 14. Species-specific height
  profiles are therefore weakly identified, and the RF is fitting a largely SHARED
  height response, not 13 independent ones. Reporting species-specific height curves as
  if they were separately estimated would overstate what these data support.
")

# --- what the RF actually does with height, by species
cat("\n  --- RF height response by species (other covariates at July means) ---\n")
jm <- d[d$month == 7, ]
grid <- expand.grid(species = levels(factor(d$species_clean)),
                    height_cm = c(50,125,200), stringsAsFactors = FALSE)
grid$dbh_m <- median(d$dbh_m, na.rm=TRUE)
grid$soil_moisture_at_tree <- mean(jm$soil_moisture_at_tree, na.rm=TRUE)
grid$soil_temp_C_mean <- mean(jm$soil_temp_C_mean, na.rm=TRUE)
grid$air_temp_C_mean  <- mean(jm$air_temp_C_mean,  na.rm=TRUE)
grid$species <- factor(grid$species, levels = levels(factor(d$species_clean)))
grid$pred <- predict(TreeRF, grid, num.threads=1)$predictions
rfw <- grid %>% dplyr::select(species, height_cm, pred) %>%
  pivot_wider(names_from=height_cm, values_from=pred, names_prefix="rf") %>%
  mutate(ratio_200_50 = rf200/rf50)
print(as.data.frame(rfw), row.names=FALSE, digits=3)
cat(sprintf("\n  RF ratio(200/50) range across species: %.2f to %.2f  (identical ratio for\n  every species would mean the RF learned no species x height interaction)\n",
            min(rfw$ratio_200_50), max(rfw$ratio_200_50)))

# ==============================================================================
rule("Q2  ALLOMETRY: DBH -> surface area")
# ==============================================================================
# fix_dbh() and the local band definition removed: canonical_inventory() already
# returns repaired diameters plus band_m (0.75 m for Kalmia, a shrub; 2 m otherwise)
# and A_band_m2. Re-repairing repaired data with species == "..." rather than %in%
# turned the 41 NA-species stems into NA diameters, which is why this reported 7,998
# stems and an NA band area.
INV$band_h <- INV$band_m
INV$A_stem <- INV$A_band_m2
sum_piD <- sum(pi*INV$dbh)

cat(sprintf("
  MEASURED BAND (0-2 m). The only allometry used is a cylinder:
        A_stem = pi * DBH * 2      (0.75 m for Kalmia latifolia, a shrub)
  Taper over 2 m is negligible. Bark rugosity is deliberately NOT corrected for: the
  chamber fluxes are themselves computed per unit 2D geometric chamber footprint, so
  applying them to a 2D geometric stem area is internally consistent. Rugosity is
  effectively baked into the measured flux, exactly as soil microtopography is baked
  into a soil flux reported per m2 of ground.

    n stems                    %d
    sum(pi*DBH)                %.0f m
    stem area 0-2 m            %.0f m2   -> index %.4f m2 m-2 ground
    basal area                 %.2f m2   -> %.1f m2 ha-1
", nrow(INV), sum_piD, sum(INV$A_stem), sum(INV$A_stem)/A_PLOT,
   sum(pi*(INV$dbh/2)^2), 1e4*sum(pi*(INV$dbh/2)^2)/A_PLOT))

cat("
  WHOLE WOODY SURFACE. There is NO height allometry in the current scaling, because
  there are no tree heights in this dataset to fit or validate one against (the
  inventory is tag/species/DBH/coordinates; the felled oak is the only tree with a
  measured height). Total woody area is taken from the literature WAI and allocated
  across stems in proportion to pi*DBH -- i.e. area assumed proportional to DBH^1.
  That allocation is itself an assumption: bole area scales ~DBH^1 x H, and branch area
  scales closer to DBH^2, so allocating on DBH^1 shifts area toward small stems
  relative to a branch-dominated reality.
")
alt <- data.frame(rule = c("pi*DBH (used)","DBH^1.5","DBH^2 (basal area)"),
  share_top10pct = c(
    { w<-pi*INV$dbh;      sum(sort(w,dec=TRUE)[1:round(.1*nrow(INV))])/sum(w) },
    { w<-INV$dbh^1.5;     sum(sort(w,dec=TRUE)[1:round(.1*nrow(INV))])/sum(w) },
    { w<-INV$dbh^2;       sum(sort(w,dec=TRUE)[1:round(.1*nrow(INV))])/sum(w) }))
alt$share_top10pct <- round(100*alt$share_top10pct,1)
cat("\n  Share of total woody area assigned to the largest 10% of stems:\n")
print(alt, row.names=FALSE)

# ==============================================================================
rule("Q3  VERTICAL DISTRIBUTION OF WOODY SURFACE AREA")
# ==============================================================================
cat("
  The scaling currently assumes woody surface area is distributed UNIFORMLY with
  height from 0 to the canopy. This is the single least defensible assumption in the
  chain, and the arithmetic below shows why.
")
for (H in c(15,20,25)) {
  bole <- sum_piD*H; bole_i <- bole/A_PLOT
  cat(sprintf("
    Assumed canopy height %2d m
      bole-only area (sum pi*DBH*H)   %7.0f m2   index %.2f
      Gauci WAI                                        %.2f
      => implied BRANCH area                           %.2f  (%.0f%% of woody surface)
", H, bole, bole_i, WAI_GAUCI, WAI_GAUCI-bole_i, 100*(WAI_GAUCI-bole_i)/WAI_GAUCI))
}
cat("
  So under the WAI we are importing, roughly two-thirds to three-quarters of the woody
  surface is BRANCH, not bole. Branch area is concentrated in the crown, i.e. high on
  the tree -- yet flux DECLINES with height in our measurements. A uniform-with-height
  assumption therefore places a large majority of the surface at heights where we have
  no data, and if anything BIASES THE BUDGET HIGH by giving crown surface the same flux
  weight as near-ground surface.

  Three vertical-distribution scenarios, as multipliers on the whole-surface budget
  (using the fitted log-linear height decline, canopy 20 m):
")
hhq <- c(50,125,200)/100
fbar <- sapply(hhq, function(z) mean(pair$f[abs(pair$h/100-z)<1e-6]))
b <- coef(lm(log(pmax(fbar,1e-9)) ~ hhq))[2]
Hc <- 20; zz <- seq(0.05, Hc, length.out=400)
prof <- exp(b*(zz-2))                       # flux relative to the 2 m value
wts <- list(
  `uniform with height (current)`       = rep(1, length(zz)),
  `bole-like (area ~ constant, tapering)` = pmax(1 - zz/Hc, 0),
  `crown-weighted (branch-dominated)`   = dnorm(zz, mean=0.75*Hc, sd=0.15*Hc))
cat(sprintf("    %-40s %10s %10s\n", "vertical distribution", "rel. flux", "vs uniform"))
u <- weighted.mean(prof, wts[[1]])
for (nm in names(wts)) {
  v <- weighted.mean(prof, wts[[nm]])
  cat(sprintf("    %-40s %10.4f %9.2fx\n", nm, v, v/u))
}
cat("
  The crown-weighted case is the physically plausible one if most woody area is branch,
  and it cuts the whole-surface budget substantially relative to uniform. We have no
  measurements to choose among these, and no measurement of woody area distribution
  with height at this site or in this dataset.
")

# ==============================================================================
rule("Q4  RELATION TO THE GAUCI et al. 2024 WAI")
# ==============================================================================
cat(sprintf("
  WAI = %.2f m2 woody area per m2 GROUND is a stand-level index taken from Gauci et al.
  2024 (Nature 631: 796-800) for temperate forests. It is imported wholesale: nothing
  in our data constrains it, because computing our own WAI would require tree heights
  and branch architecture we never measured.

  What our data DO give, and how it compares:
      our measured stem band (0-2 m)      %.4f m2 m-2   (%.1f%% of WAI)
      imported WAI                        %.4f m2 m-2
      ratio                               %.1fx

  So the whole-woody-surface scenario extrapolates our measurements over roughly %.0f
  times the area we actually sampled, into a vertical domain (crown branches) that
  differs from the sampled domain (lower bole) in bark type, diameter, tissue age,
  microclimate and likely microbial community. The extrapolation factor, not the flux
  estimate, dominates the uncertainty in the whole-surface number.
", WAI_GAUCI, sum(INV$A_stem)/A_PLOT, 100*(sum(INV$A_stem)/A_PLOT)/WAI_GAUCI,
   WAI_GAUCI, WAI_GAUCI/(sum(INV$A_stem)/A_PLOT), WAI_GAUCI/(sum(INV$A_stem)/A_PLOT)))

# ==============================================================================
# figure: species height profiles, observed (paired) and RF
# ==============================================================================
obs_long <- wide %>% pivot_longer(c(h50,h125,h200), names_to="h", values_to="f") %>%
  mutate(h = as.numeric(sub("h","",h))) %>% group_by(species_clean, h) %>%
  summarise(f=mean(f), n=n(), .groups="drop") %>% mutate(src="observed (paired)")
rf_long <- grid %>% transmute(species_clean=as.character(species), h=height_cm,
                              f=pred, n=NA_integer_, src="random forest")
both <- bind_rows(obs_long, rf_long) %>%
  group_by(species_clean) %>% filter(any(!is.na(n) & n >= 3)) %>% ungroup()

p <- ggplot(both, aes(h, f, colour=src)) +
  geom_hline(yintercept=0, colour="grey70", linewidth=.3) +
  geom_line(linewidth=.7) + geom_point(size=1.6) +
  facet_wrap(~species_clean, scales="free_y") +
  scale_colour_manual(values=c("observed (paired)"="#B2182B","random forest"="#2166AC"), name=NULL) +
  labs(title="Height-flux profile by species: measured vs random forest",
       subtitle=paste("paired within-tree means, July campaign only (n = 1-19 trees per species);",
                      "\nspecies-specific profiles are weakly identified at these sample sizes"),
       x="measurement height (cm)", y=expression(CH[4]~flux~(nmol~m^-2~s^-1))) +
  theme_bw(base_size=8) + theme(legend.position="top")
ggsave(out_path("fig_height_profiles_by_species.png"), p, width=9, height=6, dpi=200, bg="white")
cat("\nWritten: outputs/figures/generated/fig_height_profiles_by_species.png\n")
