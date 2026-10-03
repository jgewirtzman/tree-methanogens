# ==============================================================================
# rf_spec.R -- forest size and threads for the stem and soil random forests
# ------------------------------------------------------------------------------
# RF_NUM_TREES. With 800 trees the stem model's OOB R2 ranged 0.228-0.252 across six
# seeds and row orders on identical data (2026-10-01), so the printed value depended on
# the seed. At 10,000 trees the range is 0.240-0.246 (soil: 0.453-0.455), and numbers
# reported to two decimals no longer move between runs. Used by the two models and the
# two helper forests in 05_model/02_rf_models.R; scripts that refit "the model" read
# num.trees from the fitted object and so inherit it.
#
# RF_THREADS. ranger draws one seed per tree from `seed`, so a fit is identical at any
# thread count (checked: 1, 4 and 8 threads give identical predictions). Threads only
# change the run time.
# ==============================================================================
RF_NUM_TREES <- 10000
RF_THREADS   <- max(1L, min(8L, parallel::detectCores() - 2L))
