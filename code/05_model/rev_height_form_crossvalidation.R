source("code/lib/outputs.R")
#!/usr/bin/env Rscript
# ==============================================================================
# rev_height_form_crossvalidation.R
# ------------------------------------------------------------------------------
# Which height-flux form extrapolates best? Tested on the 143 trees measured at all
# three heights, NOT on the single climbed tree.
#
# WHY NOT THE FELLED OAK. An earlier version of this analysis ranked forms by RMSE
# against the climbed black oak's 6/8/10 m points. That was not a quantitatively
# meaningful test: a single chamber measurement on one tree at one time is not an
# estimate of "the flux at 6 m", and the oak carries a 10x hotspot at 4 m showing how
# much within-tree variation a single point can hide. The oak is retained as
# QUALITATIVE corroboration of the direction of the result, and nothing more.
#
# THE TEST. Each tree contributes three heights (0.5, 1.25, 2.0 m). For each tree and
# each candidate form we fit the form to TWO heights and predict the THIRD:
#
#   Test A  UPWARD extrapolation   fit 0.5 + 1.25  ->  predict 2.0   (the relevant one:
#           this is the direction the whole-surface scaling extrapolates)
#   Test B  DOWNWARD extrapolation fit 1.25 + 2.0  ->  predict 0.5   (relevant to the
#           below-0.5 m treatment in the band integral)
#   Test C  INTERPOLATION          fit 0.5 + 2.0   ->  predict 1.25  (a sanity check;
#           all reasonable forms should do well here)
#
# With two points and two free parameters each form fits its two points EXACTLY, so the
# held-out prediction is driven purely by functional form, with no fitting noise. The
# comparison replicates across 143 trees, so it has real statistical content.
#
# FORMS (2 parameters unless noted)
#   constant   f = c                       (1 param: mean of the two fitted points)
#   linear     f = a + b*h
#   exp_zero   f = a*exp(b*h)              requires positive flux at fitted points
#   power      f = a*h^b                   requires positive flux at fitted points
#
# Trees with non-positive flux at a fitting height cannot support the log-based forms;
# the number affected is reported rather than silently dropped, and a common-subset
# comparison is given so the ranking is not an artefact of differing sample sets.
#
# Output: outputs/audit/height_form_crossvalidation.txt
#         outputs/figures/generated/fig_height_form_cv.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr)})
set.seed(42)
outdir <- "outputs"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
load("outputs/models/TRAINING_DATA.RData")

pair <- tree_train_complete %>%
  filter(!is.na(measurement_height_cm), is.finite(stem_flux_corrected)) %>%
  group_by(tree_id, species_clean, h = measurement_height_cm/100) %>%
  summarise(f = mean(stem_flux_corrected), .groups="drop")
W <- pair %>% pivot_wider(names_from=h, values_from=f, names_prefix="h") %>%
  filter(is.finite(h0.5), is.finite(h1.25), is.finite(h2))
cat(sprintf("Trees with all three heights: %d\n", nrow(W)))
cat(sprintf("Trees with all three fluxes positive: %d (%.0f%%)\n",
            sum(W$h0.5>0 & W$h1.25>0 & W$h2>0),
            100*mean(W$h0.5>0 & W$h1.25>0 & W$h2>0)))

# --- predictors: fit form to (h1,f1),(h2,f2); return prediction at h3 --------
pred_form <- function(form, h1,f1,h2,f2,h3) {
  switch(form,
    constant = (f1+f2)/2,
    linear   = { b <- (f2-f1)/(h2-h1); f1 + b*(h3-h1) },
    exp_zero = if (f1>0 && f2>0) { b <- (log(f2)-log(f1))/(h2-h1); exp(log(f1)+b*(h3-h1)) } else NA_real_,
    power    = if (f1>0 && f2>0) { b <- (log(f2)-log(f1))/(log(h2)-log(h1)); exp(log(f1)+b*(log(h3)-log(h1))) } else NA_real_)
}
FORMS <- c("constant","linear","exp_zero","power")
TESTS <- list(
  A_up   = list(fit=c("h0.5","h1.25"), hf=c(0.5,1.25), out="h2",    ho=2.0,
                lab="UPWARD  fit 0.5+1.25 -> predict 2.0 m"),
  B_down = list(fit=c("h1.25","h2"),   hf=c(1.25,2.0), out="h0.5",  ho=0.5,
                lab="DOWNWARD fit 1.25+2.0 -> predict 0.5 m"),
  C_int  = list(fit=c("h0.5","h2"),    hf=c(0.5,2.0),  out="h1.25", ho=1.25,
                lab="INTERPOLATION fit 0.5+2.0 -> predict 1.25 m"))

allres <- list()
for (tn in names(TESTS)) {
  T <- TESTS[[tn]]
  f1 <- W[[T$fit[1]]]; f2 <- W[[T$fit[2]]]; obs <- W[[T$out]]
  P <- sapply(FORMS, function(fm)
    mapply(function(a,b) pred_form(fm, T$hf[1],a,T$hf[2],b,T$ho), f1, f2))
  # common subset: rows where EVERY form produced a prediction
  ok <- apply(is.finite(P),1,all) & is.finite(obs)
  cat(sprintf("\n%s\n  %s\n  evaluable on %d of %d trees (log-based forms need positive flux at both fitted heights)\n",
              strrep("-",78), T$lab, sum(ok), nrow(W)))
  cat(sprintf("  %-10s %10s %10s %10s %10s %9s\n","form","RMSE","MAE","medAE","bias","rho"))
  rr <- lapply(FORMS, function(fm) {
    p <- P[ok,fm]; o <- obs[ok]
    data.frame(test=tn, form=fm, rmse=sqrt(mean((p-o)^2)), mae=mean(abs(p-o)),
               medae=median(abs(p-o)), bias=mean(p-o),
               rho=suppressWarnings(cor(p,o,method="spearman")), n=sum(ok))
  }) %>% bind_rows() %>% arrange(rmse)
  for (i in seq_len(nrow(rr)))
    cat(sprintf("  %-10s %10.4f %10.4f %10.4f %+10.4f %9.3f%s\n", rr$form[i], rr$rmse[i],
                rr$mae[i], rr$medae[i], rr$bias[i], rr$rho[i], if(i==1) "   <- best" else ""))
  allres[[tn]] <- rr
}
res <- bind_rows(allres)
write.csv(res, out_path("height_form_crossvalidation.csv"), row.names=FALSE)

cat(sprintf("\n%s\n SUMMARY\n%s\n", strrep("=",78), strrep("=",78)))
cat("
  The UPWARD test (A) is the one that matters: it is the same direction and the same
  kind of extrapolation the whole-woody-surface scaling performs. Its ranking, computed
  across 143 independent trees, is the defensible basis for choosing a form -- or for
  concluding that no form is preferred.
")
best_up <- res %>% filter(test=="A_up") %>% arrange(rmse)
cat(sprintf("\n  Upward-extrapolation ranking: %s\n",
            paste(sprintf("%s (%.4f)", best_up$form, best_up$rmse), collapse="  <  ")))
cat(sprintf("  Bias (predicted - observed) at 2 m: %s\n",
            paste(sprintf("%s %+.4f", best_up$form, best_up$bias), collapse="   ")))

# paired comparison of the two leading forms, so the ranking is not just a point estimate
o <- W[[TESTS$A_up$out]]
Pu <- sapply(FORMS, function(fm) mapply(function(a,b)
  pred_form(fm,0.5,a,1.25,b,2.0), W$h0.5, W$h1.25))
oku <- apply(is.finite(Pu),1,all) & is.finite(o)
e <- abs(Pu[oku,] - o[oku])
cat("\n  Paired tests on absolute error, upward extrapolation (Wilcoxon signed-rank):\n")
for (i in 1:(length(FORMS)-1)) for (j in (i+1):length(FORMS)) {
  p <- suppressWarnings(wilcox.test(e[,FORMS[i]], e[,FORMS[j]], paired=TRUE)$p.value)
  cat(sprintf("    %-9s vs %-9s  median |err| %.4f vs %.4f   p = %.4g%s\n",
      FORMS[i], FORMS[j], median(e[,FORMS[i]]), median(e[,FORMS[j]]), p,
      if (p<0.05) "  *" else ""))
}

# --- figure ------------------------------------------------------------------
pd <- bind_rows(lapply(names(TESTS), function(tn)
  res %>% filter(test==tn) %>% mutate(lab=TESTS[[tn]]$lab)))
p <- ggplot(pd, aes(reorder(form,-rmse), rmse, fill=form)) +
  geom_col() + coord_flip() + facet_wrap(~lab, scales="free_x", ncol=1) +
  scale_fill_brewer(palette="Set2", guide="none") +
  labs(title="Height-form comparison by leave-one-height-out on 143 three-height trees",
       subtitle="each form fitted to two heights and used to predict the third; lower is better",
       x=NULL, y="RMSE (nmol m-2 s-1)") +
  theme_bw(base_size=9)
ggsave(out_path("fig_height_form_cv.png"), p, width=8, height=6, dpi=200, bg="white")
cat("\nWritten: outputs/figures/generated/fig_height_form_cv.png\n")
