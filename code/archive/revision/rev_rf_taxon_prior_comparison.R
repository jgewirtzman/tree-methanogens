#!/usr/bin/env Rscript
# ==============================================================================
# rev_rf_taxon_prior_comparison.R
# ------------------------------------------------------------------------------
# How should the model inform species it never measured?
#
# The pipeline uses `taxon_prior_asinh`: the MEDIAN residual of an environment-only
# RF, taken at the finest taxonomic level with n >= 5 (species -> genus -> family).
# It is the 4th most important tree predictor, so the choice matters. Three problems:
#   * self-inclusion   - a row's own residual is inside its group's median
#   * no shrinkage     - Prunus (n=11) is trusted as much as Tsuga (n=209)
#   * coarse fallback  - Linnaean ranks treated as equal-sized steps
#
# COMPARED HERE (all computed INSIDE each training fold, never on held-out data):
#   none        no taxon prior at all                       (baseline)
#   current     median residual, species -> genus -> family (as implemented)
#   loo         same, leave-one-out                         (removes self-inclusion)
#   hier        partial pooling: (1|genus/species) BLUPs    (shrinkage by n)
#   traits      per-species functional traits as features   (mechanistic)
#   hier+traits both
#
# EVALUATION: grouped 5-fold CV, holding out WHOLE TREES. Upscaling predicts flux for
# 7,015 inventory stems the model has never seen, so a metric that lets a tree's own
# siblings sit in the training set (OOB, or random CV) is measuring the wrong thing --
# it reads ~0.33 where grouped CV reads ~0.15.
#
# Output: outputs/revision/taxon_prior_comparison.csv
# ==============================================================================

suppressMessages({library(ranger); library(dplyr); library(lme4)})
set.seed(42)

load("outputs/models/RF_MODELS.RData")
load("outputs/models/TRAINING_DATA.RData")

X0 <- as.data.frame(X_tree)
d  <- tree_train_complete
y  <- asinh(d$stem_flux_corrected)
keep <- is.finite(y)
X0 <- X0[keep, , drop = FALSE]; d <- d[keep, ]; y <- y[keep]
species <- d$species; tree_id <- d$tree_id

# taxonomy (species -> genus -> family)
genus  <- sub(" .*", "", species)
fam_map <- c(Acer="Sapindaceae", Betula="Betulaceae", Carya="Juglandaceae",
             Fagus="Fagaceae", Fraxinus="Oleaceae", Kalmia="Ericaceae",
             Pinus="Pinaceae", Prunus="Rosaceae", Quercus="Fagaceae",
             Sassafras="Lauraceae", Tsuga="Pinaceae")
family <- unname(fam_map[genus])

# --- per-species traits -----------------------------------------------------
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
TRAIT_COLS <- c("wood_density_gcm3","bark_density_gcm3","bark_wood_ratio","porosity_num",
                "wood_heartwood_pH","wood_sapwood_pH","try_rooting_depth",
                "try_plant_longevity","gymnosperm")
sp_traits <- NULL
if (file.exists(TRAITS)) {
  tr <- read.csv(TRAITS, stringsAsFactors = FALSE)
  have <- intersect(TRAIT_COLS, names(tr))
  sp_traits <- tr %>% group_by(spcode) %>%
    summarise(across(all_of(have), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(across(-spcode, ~ifelse(is.nan(.x), NA_real_, .x)))
  cat("species traits loaded:", nrow(sp_traits), "species x", length(have), "traits\n")
  cat("  ", paste(have, collapse = ", "), "\n\n")
}
# map species binomial -> 4-letter code used in the trait table
code_of <- function(sp) {
  p <- strsplit(sp, " ")
  vapply(p, function(z) if (length(z) >= 2)
    toupper(paste0(substr(z[1],1,2), substr(z[2],1,2))) else NA_character_, character(1))
}
spcode <- code_of(species)

# --- prior builders: fitted on TRAIN rows only, applied to any rows ----------
env_cols <- intersect(c("air_temp_C_mean","soil_temp_C_mean","soil_moisture_at_tree",
                        "dbh_within_z","height_cm"), names(X0))

resid_of <- function(tr) {                       # environment-only residuals
  Xe <- X0[tr, env_cols, drop = FALSE]
  rf <- ranger(x = Xe, y = y[tr], num.trees = 200, min.node.size = 5,
               num.threads = 1, seed = 42)
  y[tr] - rf$predictions
}

prior_median <- function(tr, loo = FALSE) {
  r <- resid_of(tr)
  key <- list(species = species[tr], genus = genus[tr], family = family[tr])
  agg <- function(k) { s <- split(r, k); s <- s[vapply(s, length, 1L) >= 5]
                       vapply(s, median, numeric(1), na.rm = TRUE) }
  ms <- agg(key$species); mg <- agg(key$genus); mf <- agg(key$family)
  out <- rep(0, length(species))
  hit <- function(v, tab) ifelse(!is.na(match(v, names(tab))), tab[match(v, names(tab))], NA)
  out <- ifelse(!is.na(hit(species, ms)), hit(species, ms),
         ifelse(!is.na(hit(genus, mg)),  hit(genus, mg),
         ifelse(!is.na(hit(family, mf)), hit(family, mf), 0)))
  out <- as.numeric(ifelse(is.na(out), 0, out))
  if (loo) {   # remove each TRAIN row's own contribution from its species median
    for (i in tr) {
      same <- tr[species[tr] == species[i]]
      if (length(same) >= 6) out[i] <- median(r[match(setdiff(same, i), tr)], na.rm = TRUE)
    }
    out[is.na(out)] <- 0
  }
  out
}

prior_hier <- function(tr) {                     # partial pooling, shrinks by n
  r <- resid_of(tr)
  df <- data.frame(r = r, sp = species[tr], gn = genus[tr])
  m <- tryCatch(suppressMessages(lmer(r ~ 1 + (1 | gn/sp), data = df)),
                error = function(e) NULL)
  if (is.null(m)) return(rep(0, length(species)))
  re <- ranef(m)
  b_sp <- setNames(re$`sp:gn`[, 1], rownames(re$`sp:gn`))
  b_gn <- setNames(re$gn[, 1], rownames(re$gn))
  key <- paste(species, genus, sep = ":")
  out <- ifelse(!is.na(match(key, names(b_sp))), b_sp[match(key, names(b_sp))],
         ifelse(!is.na(match(genus, names(b_gn))), b_gn[match(genus, names(b_gn))], 0))
  as.numeric(ifelse(is.na(out), 0, out))
}

trait_block <- function() {
  if (is.null(sp_traits)) return(NULL)
  m <- sp_traits[match(spcode, sp_traits$spcode), setdiff(names(sp_traits), "spcode"), drop = FALSE]
  as.data.frame(lapply(m, function(v) { v[is.na(v)] <- median(v, na.rm = TRUE); as.numeric(v) }))
}
TB <- trait_block()

# --- grouped CV -------------------------------------------------------------
ut <- unique(tree_id); set.seed(1)
gmap <- setNames(sample(rep(1:5, length.out = length(ut))), ut)
fold <- unname(gmap[tree_id])
r2 <- function(a, b) 1 - sum((a - b)^2) / sum((a - mean(a))^2)

CONFIGS <- c("none","current","loo","hier","traits","hier+traits")
res <- lapply(CONFIGS, function(cfg) {
  pr <- rep(NA_real_, length(y))
  for (k in 1:5) {
    tr <- which(fold != k); te <- which(fold == k)
    Xb <- X0[, setdiff(names(X0), "taxon_prior_asinh"), drop = FALSE]
    if (cfg %in% c("current","loo","hier","hier+traits")) {
      p <- switch(cfg,
        current       = prior_median(tr, loo = FALSE),
        loo           = prior_median(tr, loo = TRUE),
        prior_hier(tr))
      Xb$taxon_prior_asinh <- p
    }
    if (cfg %in% c("traits","hier+traits") && !is.null(TB)) Xb <- cbind(Xb, TB)
    m <- ranger(x = Xb[tr, , drop = FALSE], y = y[tr], num.trees = 800, min.node.size = 5,
                mtry = max(1, floor(sqrt(ncol(Xb)))), num.threads = 1, seed = 42)
    pr[te] <- predict(m, Xb[te, , drop = FALSE])$predictions
  }
  data.frame(prior = cfg, r2 = r2(y, pr),
             spearman = cor(y, pr, method = "spearman"),
             mean_ratio = mean(sinh(pr)) / mean(sinh(y)))
}) %>% bind_rows()

cat("\n", strrep("=", 66), "\n", sep = "")
cat("  TAXON PRIOR COMPARISON - grouped 5-fold CV (whole trees held out)\n")
cat(strrep("=", 66), "\n\n", sep = "")
print(res %>% arrange(desc(r2)), row.names = FALSE, digits = 3)

dir.create("outputs/revision", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "outputs/revision/taxon_prior_comparison.csv", row.names = FALSE)
cat("\nWritten: outputs/revision/taxon_prior_comparison.csv\n")
