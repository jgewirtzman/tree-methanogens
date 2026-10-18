# ==============================================================================
# run_mode.R -- every run-size setting, for the full and the lite run
# ------------------------------------------------------------------------------
#   Rscript code/run_all.R          full: the canonical run (the values below, RUN_FULL)
#   Rscript code/run_all.R --lite   lite: every script, in minutes, in a separate copy
#                                   of the repository (run_all.R builds it)
#
# The mode is the environment variable PIPELINE_RUN_MODE ("full" when unset). run_all.R
# --lite sets it only for the copy it builds, whose root holds the marker file
# .lite_run. Lite settings anywhere else stop the script, so a lite run can never
# write the canonical outputs. A script reads RUN$<name>; nothing else in the code
# sets the size of a forest, a cross-validation or a resampling loop.
#
# Lite numbers are not results. They show that every script runs and that the pieces
# fit together, nothing more.
# ==============================================================================

RUN_FULL <- list(
  # the stem and soil forests (05_model/02_rf_models.R, through rf_spec.R). Scripts that
  # refit the fitted model (04, 05, 07, 15) take num.trees from the forest, so follow it.
  rf_num_trees            = 10000,
  # forests refitted inside audits and the Fig S26 importance splits
  # (05_model/09, 10, 12; 08_figures/figS24)
  refit_num_trees         = 800,
  grouped_cv_repeats      = 100,          # 05_model/04_rf_grouped_cv.R (and 15, through 04)
  selection_cv_repeats    = 30,           # 05_model/07_rf_predictor_selection.R
  audit_seeds             = c(42, 1, 2),  # 05_model/05_audit_training_population.R
  pooling_repeats         = 4,            # 05_model/10_rf_species_pooling.R
  family_repeats          = 3,            # 05_model/12_model_family_comparison.R
  tuning_num_trees        = 1000,         # 05_model/16_rf_hyperparameter_tuning.R
  tuning_repeats          = 10,           #   outer repeats
  tuning_node_sizes       = c(1, 3, 5, 10, 20),   # min.node.size grid (must hold the fitted 5)
  daily_driver_stems      = 300,          # 06_upscale/11_sens_daily_drivers.R: sampled stems
  daily_driver_day_step   = 1,            #   every day (lite: every 30th)
  gene_cv_repeats         = 30,           # 09_tables_stats/19_gene-rf.R
  importance_splits       = 50,           # 08_figures/figS24 (Fig S26c, d): 70/30 splits
  importance_permutations = 20            #   permutations per predictor and split
)

RUN_LITE <- modifyList(RUN_FULL, list(
  rf_num_trees = 200, refit_num_trees = 100,
  grouped_cv_repeats = 5, selection_cv_repeats = 3, audit_seeds = 42,
  pooling_repeats = 1, family_repeats = 1,
  tuning_num_trees = 100, tuning_repeats = 1, tuning_node_sizes = c(5, 20),
  daily_driver_stems = 30, daily_driver_day_step = 30,
  gene_cv_repeats = 3, importance_splits = 5, importance_permutations = 3))
stopifnot(identical(names(RUN_LITE), names(RUN_FULL)))

RUN_MODE <- Sys.getenv("PIPELINE_RUN_MODE", "full")
if (!RUN_MODE %in% c("full", "lite"))
  stop("PIPELINE_RUN_MODE must be 'full' or 'lite', not '", RUN_MODE, "'", call. = FALSE)
if (RUN_MODE == "lite") local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "code/pipeline.csv")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, ".lite_run")))
    stop("PIPELINE_RUN_MODE=lite outside a lite copy (no .lite_run at ", d, "); ",
         "run lite with  Rscript code/run_all.R --lite", call. = FALSE)
})
RUN <- if (RUN_MODE == "lite") RUN_LITE else RUN_FULL
