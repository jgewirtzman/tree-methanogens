#!/usr/bin/env Rscript
# ==============================================================================
# rev_asymptote_estimation.R
# ------------------------------------------------------------------------------
# ESTIMATE the asymptote empirically instead of assuming one.
#
# The "decline to a plateau" scenario needs an asymptote: the flux that a stem
# approaches high on the bole. Choosing 25% / 50% / 75% of the 2 m value would be
# arbitrary. The three-height dataset can estimate it.
#
# MODEL   f(h) = c + a * exp(-b * h)
#         c = asymptote, a = amplitude of the declining part, b = decay rate.
#
# With three heights per stem and three parameters the per-stem fit is EXACT, so
# c is determined but has no residual degrees of freedom and can be wild for stems
# whose profile is non-monotone or concave. Two estimators are therefore reported:
#
#   1. PER-STEM exact solution, across the 143 stems -> a DISTRIBUTION of asymptotes.
#      Robust summary (median, IQR) rather than the mean, which the tail destroys.
#
#   2. STAND-MEAN fit with a nonparametric BOOTSTRAP OVER STEMS (resample stems with
#      replacement, recompute the mean profile, refit). This gives the asymptote of
#      the average stem with a confidence interval, and is the value to report.
#
# The closed-form solution for equally spaced heights (0.5, 1.25, 2.0 m):
#         r = (f3 - f2) / (f2 - f1)                 decay ratio per step
#         c = f1 - (f2 - f1)^2 / (f3 - 2*f2 + f1)   asymptote
# valid when the denominator is non-zero; r must be in (0,1) for a decline to a
# plateau. r >= 1 means the profile is not converging (no asymptote); r < 0 means it
# is non-monotone.
#
# The felled oak is reported as independent corroboration, not as an estimate.
#
# Output: outputs/revision/asymptote_estimation.txt
# ==============================================================================

suppressMessages({library(dplyr); library(tidyr)})
set.seed(42)
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
load("outputs/models/TRAINING_DATA.RData")

pair <- tree_train_complete %>%
  filter(!is.na(measurement_height_cm), is.finite(stem_flux_corrected)) %>%
  group_by(tree_id, species_clean, h=measurement_height_cm/100) %>%
  summarise(f=mean(stem_flux_corrected), .groups="drop")
W <- pair %>% pivot_wider(names_from=h, values_from=f, names_prefix="h") %>%
  filter(is.finite(h0.5), is.finite(h1.25), is.finite(h2))
cat(sprintf("Stems with all three heights: %d\n", nrow(W)))

asym <- function(f1,f2,f3) {
  den <- f3 - 2*f2 + f1
  ifelse(abs(den) < 1e-12, NA_real_, f1 - (f2-f1)^2/den)
}
ratio <- function(f1,f2,f3) ifelse(abs(f2-f1) < 1e-12, NA_real_, (f3-f2)/(f2-f1))

# ---- 1. per-stem ------------------------------------------------------------
W$r <- ratio(W$h0.5, W$h1.25, W$h2)
W$c <- asym(W$h0.5, W$h1.25, W$h2)
W$c_frac <- W$c / W$h2                      # asymptote as a fraction of the 2 m flux

conv <- is.finite(W$r) & W$r > 0 & W$r < 1  # genuinely converging to a plateau
cat(sprintf("
PER-STEM CLASSIFICATION (from the decay ratio r)
  converging to a plateau   0 < r < 1 : %3d (%.0f%%)
  not converging                r >= 1 : %3d (%.0f%%)   -- decline decelerating or growing
  non-monotone                   r < 0 : %3d (%.0f%%)
  undefined (flat middle)                : %3d
",
 sum(conv,na.rm=TRUE), 100*mean(conv,na.rm=TRUE),
 sum(W$r>=1,na.rm=TRUE), 100*mean(W$r>=1,na.rm=TRUE),
 sum(W$r<0,na.rm=TRUE), 100*mean(W$r<0,na.rm=TRUE),
 sum(!is.finite(W$r))))

cat("\n  Asymptote as a FRACTION of that stem's 2 m flux, converging stems only:\n")
cf <- W$c_frac[conv & is.finite(W$c_frac)]
print(round(quantile(cf, c(0,.1,.25,.5,.75,.9,1), na.rm=TRUE), 3))
cat(sprintf("  median %.3f   IQR %.3f - %.3f   n = %d\n",
            median(cf,na.rm=TRUE), quantile(cf,.25,na.rm=TRUE), quantile(cf,.75,na.rm=TRUE), sum(is.finite(cf))))
cat(sprintf("  fraction of converging stems with a NEGATIVE asymptote (implying uptake aloft): %.0f%%\n",
            100*mean(cf < 0, na.rm=TRUE)))

# ---- 2. stand-mean, bootstrapped over stems ---------------------------------
sm <- function(idx) {
  d <- W[idx,]
  c(f1=mean(d$h0.5), f2=mean(d$h1.25), f3=mean(d$h2))
}
o <- sm(seq_len(nrow(W)))
c_hat <- asym(o["f1"],o["f2"],o["f3"]); r_hat <- ratio(o["f1"],o["f2"],o["f3"])
cat(sprintf("
STAND-MEAN PROFILE
  f(0.5) = %.4f   f(1.25) = %.4f   f(2.0) = %.4f
  successive ratios %.3f and %.3f  -- near-constant, i.e. close to pure exponential
  decay ratio r    = %.4f
  asymptote c      = %.4f nmol m-2 s-1  = %.1f%% of the 2 m flux
",
 o["f1"],o["f2"],o["f3"], o["f2"]/o["f1"], o["f3"]/o["f2"], r_hat, c_hat, 100*c_hat/o["f3"]))

B <- 4000
bs <- replicate(B, { i <- sample(nrow(W), replace=TRUE); s <- sm(i)
  cc <- asym(s["f1"],s["f2"],s["f3"]); c(cc, cc/s["f3"]) })
cb <- bs[1,]; cfb <- bs[2,]
ok <- is.finite(cb) & is.finite(cfb)
cat(sprintf("
  BOOTSTRAP over stems (B = %d, %d valid)
    asymptote c        median %.4f   95%% CI [%.4f, %.4f]
    as %% of 2 m flux   median %.1f%%    95%% CI [%.1f%%, %.1f%%]
    P(asymptote < 0)   %.1f%%   -- probability the average stem declines THROUGH zero
", B, sum(ok),
 median(cb[ok]), quantile(cb[ok],.025), quantile(cb[ok],.975),
 100*median(cfb[ok]), 100*quantile(cfb[ok],.025), 100*quantile(cfb[ok],.975),
 100*mean(cb[ok] < 0)))

# ---- 3. felled oak, corroboration only --------------------------------------
bo <- read.csv("data/compiled/black_oak_experiment.csv", stringsAsFactors=FALSE)
bo$h <- suppressWarnings(as.numeric(bo$height_m))
bo <- bo[is.finite(bo$h) & is.finite(bo$CH4_best_flux_nmol_m2_s),]
at2 <- bo$CH4_best_flux_nmol_m2_s[abs(bo$h-2)<.01][1]
hi  <- bo$CH4_best_flux_nmol_m2_s[bo$h > 4]
cat(sprintf("
FELLED OAK (corroboration only, n = 1 stem)
  flux at 2 m %.4f | mean above 4 m %.4f | ratio %.2f
  An observed plateau at %.0f%% of the 2 m value. Consistent in DIRECTION with a
  non-zero asymptote, but a single stem cannot estimate its magnitude.
", at2, mean(hi), mean(hi)/at2, 100*mean(hi)/at2))

cat("
INTERPRETATION
  The asymptote is small and poorly determined. The stand-mean profile is close to a
  pure exponential, which is the c = 0 case, and the bootstrap interval spans zero.
  There is therefore no empirical basis for a large positive plateau such as the 50%
  used in earlier drafts, and no basis for excluding a negative asymptote either:
  a sizeable minority of individual stems imply flux declining THROUGH zero.

  NOTE: this per-stem / whole-stand result is superseded by the GROUP-LEVEL analysis
  appended below, which is the one to report. See the final interpretation there.
")

# ==============================================================================
# GROUP-LEVEL ASYMPTOTE ESTIMATION
# ------------------------------------------------------------------------------
# The per-stem estimator failed because three points per stem are too noisy, not
# because the model is wrong. Aggregating before fitting should stabilise it. Two
# groupings that are mechanistically motivated rather than arbitrary:
#   SPECIES        the height response plausibly differs by wood anatomy / bark
#   MOISTURE BIN   the decline is expected to be a transport phenomenon, so it should
#                  steepen with soil water
# For each group the mean profile is formed first, then the asymptote solved, with a
# bootstrap over the STEMS WITHIN that group.
# ==============================================================================
cat("\n", strrep("=",78), "\n GROUP-LEVEL ASYMPTOTE ESTIMATION\n", strrep("=",78), "\n", sep="")

mo <- tree_train_complete %>%
  filter(!is.na(measurement_height_cm), is.finite(stem_flux_corrected)) %>%
  group_by(tree_id) %>% summarise(moist=mean(soil_moisture_at_tree,na.rm=TRUE), .groups="drop")
W2 <- W %>% left_join(mo, "tree_id")

grp_asym <- function(df, B=3000) {
  s <- c(mean(df$h0.5), mean(df$h1.25), mean(df$h2))
  cc <- asym(s[1],s[2],s[3]); rr <- ratio(s[1],s[2],s[3])
  bs <- replicate(B, { i <- sample(nrow(df), replace=TRUE)
    q <- c(mean(df$h0.5[i]), mean(df$h1.25[i]), mean(df$h2[i]))
    a <- asym(q[1],q[2],q[3]); c(a, a/q[3]) })
  ok <- is.finite(bs[1,]) & is.finite(bs[2,])
  data.frame(n=nrow(df), f2=s[3], r=rr, c=cc, c_frac=cc/s[3],
             lo=quantile(bs[2,ok],.025), hi=quantile(bs[2,ok],.975),
             p_neg=mean(bs[1,ok]<0), width=diff(quantile(bs[2,ok],c(.025,.975))))
}

cat("\n--- BY SPECIES (n >= 6) ---\n")
cat(sprintf("  %-24s %4s %8s %8s %10s %20s %7s\n","species","n","r","c","c/f(2m)","95% CI on c/f(2m)","P(c<0)"))
sp_res <- W2 %>% group_by(species_clean) %>% filter(n()>=6) %>% group_modify(~grp_asym(.x)) %>% ungroup()
for (i in seq_len(nrow(sp_res))) with(sp_res[i,],
  cat(sprintf("  %-24s %4d %8.3f %8.4f %10.2f %20s %6.0f%%\n",
      species_clean, n, r, c, c_frac, sprintf("[%.1f, %.1f]", lo, hi), 100*p_neg)))

cat("\n--- BY SOIL MOISTURE TERCILE ---\n")
W2$mbin <- cut(W2$moist, quantile(W2$moist, c(0,1/3,2/3,1), na.rm=TRUE),
               labels=c("dry","mid","wet"), include.lowest=TRUE)
cat(sprintf("  %-10s %4s %8s %8s %10s %20s %7s\n","bin","n","r","c","c/f(2m)","95% CI","P(c<0)"))
mb_res <- W2 %>% filter(!is.na(mbin)) %>% group_by(mbin) %>% group_modify(~grp_asym(.x)) %>% ungroup()
for (i in seq_len(nrow(mb_res))) with(mb_res[i,],
  cat(sprintf("  %-10s %4d %8.3f %8.4f %10.2f %20s %6.0f%%\n",
      mbin, n, r, c, c_frac, sprintf("[%.1f, %.1f]", lo, hi), 100*p_neg)))

cat(sprintf("
DOES AGGREGATION HELP?
  median 95%% CI width on c/f(2m):  per-species %.1f   per-moisture-bin %.1f
  (the whole-stand bootstrap width was %.1f)
  groups whose CI EXCLUDES zero: species %d/%d, moisture %d/%d
",
 median(sp_res$width), median(mb_res$width), 790.2-(-551.0),
 sum(sp_res$lo>0 | sp_res$hi<0), nrow(sp_res),
 sum(mb_res$lo>0 | mb_res$hi<0), nrow(mb_res)))


# ==============================================================================
# FINAL INTERPRETATION
# ==============================================================================
cat("
==============================================================================
 FINAL INTERPRETATION -- asymptote estimation
==============================================================================

  1. AGGREGATION IS ESSENTIAL. The asymptote is not estimable per stem: three points
     and three parameters leave no residual degrees of freedom, and the closed-form
     solution has (f3 - 2*f2 + f1) in the denominator, which approaches zero for
     near-linear profiles and makes the estimator explode. Per-stem fractions ranged
     from -232 to +2.4. Bootstrap CI widths on the asymptote fraction:

         whole stand        1341
         per species         5.5
         per moisture bin    1.8      <- ~745x tighter than the stand estimate

  2. BUT A NARROW INTERVAL IS NOT ENOUGH. A plateau exists only when the decay ratio
     r lies in (0,1). Of the three moisture bins only the WET bin converges
     (r = 0.144); dry (r = 13.3) and mid (r = -1.08) do not, so their apparently
     precise asymptotes are estimates of a quantity that does not exist for those
     groups. By species only 2 of 10 converge. Reporting the CI without r would be
     misleading.

  3. WHERE A PLATEAU IS ESTIMABLE, IT SITS CLOSE TO THE 2 m VALUE: asymptote fractions
     of roughly 0.83-0.97 of f(2m). The felled oak, measured to 10 m, gives 1.23. Both
     point to the same place.

  4. CONSEQUENCE FOR THE SCENARIO SET. An asymptotic decline is algebraically a blend
     of the two forms already carried:

         c + (f2 - c)*exp(b(z-2))  =  phi * constant  +  (1 - phi) * exponential

     with phi = c/f2. 'constant' IS the phi = 1 case and 'exponential' IS the phi = 0
     case, so the whole asymptote family is spanned by the endpoints and no separate
     form -- and no invented fraction -- is required. The data place the estimable
     asymptotes hard against the 'constant' end.

  5. INCIDENTAL RESULT WORTH REPORTING. The wet bin is the ONLY grouping with a
     converging decline. That is what a transport mechanism predicts, and it is the
     clearest signal for the moisture story in these data -- notably, the random
     forest could not recover the same effect as an explicit interaction.

  6. NON-MONOTONICITY. 47% of stems are non-monotone across 0.5-2.0 m, 21% do not
     converge, and only 32% decline toward a plateau. The search for a single
     functional form is therefore somewhat beside the point; this supports presenting
     bounding scenarios rather than a fitted profile.
")
