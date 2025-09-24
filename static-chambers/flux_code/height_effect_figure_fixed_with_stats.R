# =========================
# Full analysis & reporting
# =========================

# ---- Libraries ----
library(ggplot2)
library(dplyr)
library(lme4)
library(gridExtra)
library(grid)
library(purrr)    # for map2 used in breaks
library(tidyr)
library(readr)
library(glue)

# =========================
#    Data prep & styling
# =========================

# Species name mapping
species_mapping <- c(
  "ACRU" = "Acer rubrum",
  "ACSA" = "Acer saccharum", 
  "BEAL" = "Betula alleghaniensis",
  "BELE" = "Betula lenta",
  "BEPA" = "Betula papyrifera",
  "FAGR" = "Fagus grandifolia",
  "FRAM" = "Fraxinus americana",
  "PIST" = "Pinus strobus",
  "QURU" = "Quercus rubra",
  "TSCA" = "Tsuga canadensis",
  "CAOV" = "Carya ovata",
  "KALA" = "Kalmia latifolia",
  "PRSE" = "Prunus serotina",
  "QUAL" = "Quercus alba",
  "QUVE" = "Quercus velutina",
  "SAAL" = "Sassafras albidum"
)

# Define publication colors and shapes - solid and hollow points
colors_pub <- c(
  "Significant Negative" = "#D73027",      # Deep red
  "Significant Positive" = "#4575B4",      # Deep blue  
  "Non-significant Negative" = "#F46D43",  # Light red-orange
  "Non-significant Positive" = "#74ADD1"   # Light blue
)

shapes_pub <- c(
  "Significant Negative" = 16,    # Solid circle
  "Significant Positive" = 16,    # Solid circle
  "Non-significant Negative" = 1, # Hollow circle
  "Non-significant Positive" = 1  # Hollow circle
)

# Color for boxplots - single subtle color
boxplot_color <- "#756BB1"  # Muted purple

# =========================
#       Top panel data
# =========================

# Expecting a data.frame `final_dataset` in your environment
plot_data_top <- final_dataset %>%
  mutate(height_numeric = as.numeric(measurement_height)) %>%
  filter(!is.na(CH4_best.flux) & !is.na(height_numeric)) %>%
  # Fix TSLA -> TSCA typo, remove Kalmia, convert height to meters, and map to Latin names
  mutate(
    species = ifelse(species == "TSLA", "TSCA", species),
    height_m = height_numeric / 100  # Convert cm to meters
  ) %>%
  filter(species != "KALA") %>%  # Remove Kalmia latifolia
  mutate(
    tree_unique = paste(plot, tree_id, sep = "_"),
    species_latin = species_mapping[species]
  )

# Calculate tree counts per species
tree_counts <- plot_data_top %>%
  group_by(species, species_latin) %>%
  summarise(n_trees = n_distinct(tree_unique), .groups = 'drop') %>%
  mutate(species_label = paste0(species_latin, " (n=", n_trees, ")"))

# Add labels to plot data
plot_data_top <- plot_data_top %>%
  left_join(tree_counts %>% select(species, species_label), by = "species") %>%
  mutate(species_label = factor(species_label, levels = sort(unique(species_label))))

# Calculate actual data breaks for each species
species_breaks <- plot_data_top %>%
  group_by(species, species_label) %>%
  summarise(
    data_min = min(CH4_best.flux, na.rm = TRUE),
    data_max = max(CH4_best.flux, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    data_mid = (data_min + data_max) / 2,
    breaks_list = map2(data_min, data_max, ~ c(.x, (.x + .y)/2, .y))
  )

# Create a named list of breaks for each species
breaks_lookup <- setNames(species_breaks$breaks_list, species_breaks$species_label)

# Top panel: Boxplots with breaks calculated from original data (unchanged visual)
p_top <- ggplot(plot_data_top, aes(x = factor(height_m), y = CH4_best.flux)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.7) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, fill = boxplot_color, color = "gray30") +
  geom_jitter(alpha = 0.4, width = 0.2, size = 0.8, color = "gray20") +
  facet_wrap(~ species_label, scales = "free", ncol = 3) +
  coord_flip() +
  scale_y_continuous(
    breaks = function(x) {
      current_panel <- unique(plot_data_top$species_label[
        plot_data_top$CH4_best.flux >= min(x, na.rm = TRUE) & 
          plot_data_top$CH4_best.flux <= max(x, na.rm = TRUE)
      ])
      if(length(current_panel) == 1 && current_panel %in% names(breaks_lookup)) {
        return(breaks_lookup[[current_panel]])
      } else {
        max_val <- max(x, na.rm = TRUE)
        min_val <- min(x, na.rm = TRUE)
        if(is.na(max_val) || is.na(min_val)) return(c(0))
        if(max_val == min_val) return(c(min_val))
        mid_val <- (min_val + max_val) / 2
        return(c(min_val, mid_val, max_val))
      }
    },
    labels = function(x) {
      ifelse(abs(x) < 0.001, "0", 
             ifelse(abs(x) < 0.01, sprintf("%.3f", x),
                    ifelse(abs(x) < 0.1, sprintf("%.2f", x), 
                           sprintf("%.1f", x))))
    }
  ) +
  labs(
    x = "Height (m)",
    y = expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")")
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8, face = "italic"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

# Prepare data for bottom panel (coefficients) - unchanged structure
analysis_data <- plot_data_top %>%
  mutate(tree_unique = paste(plot, tree_id, sep = "_"))

# =========================================
#     STATS (with provenance tracking)
# =========================================

# Quiet try wrapper
trySuppressWarnings <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}

# Safely get a 95% CI for the "height_numeric" slope with provenance
get_ci_for_height <- function(model) {
  # Try PROFILE first (more accurate), fallback to WALD
  ci_method_used <- NA_character_
  ci_try <- trySuppressWarnings(try(
    confint(model, parm = "height_numeric", method = "profile", oldNames = FALSE),
    silent = TRUE
  ))
  if (!inherits(ci_try, "try-error") && !any(is.na(ci_try))) {
    ci_method_used <- "profile"
  } else {
    ci_try <- confint(model, parm = "height_numeric", method = "Wald", oldNames = FALSE)
    ci_method_used <- "wald_fallback"
  }
  list(
    low = as.numeric(ci_try[1]),
    high = as.numeric(ci_try[2]),
    ci_method = ci_method_used
  )
}

# Core function (fits model, extracts estimate/SE/t, CI + provenance)
test_species_height_effect <- function(species_name, data) {
  species_data <- data %>% filter(species == species_name)
  
  # Repeated heights check at tree level
  n_obs <- nrow(species_data)
  n_trees <- length(unique(species_data$tree_unique))
  n_trees_multi_height <- species_data %>%
    group_by(tree_unique) %>%
    summarise(n_heights = n_distinct(height_numeric), .groups = "drop") %>%
    filter(n_heights > 1) %>%
    nrow()
  
  if (n_trees_multi_height < 3) {
    return(tibble(
      species = species_name,
      n_obs = n_obs,
      n_trees = n_trees,
      n_trees_multi_height = n_trees_multi_height,
      height_coef = NA_real_,
      height_se = NA_real_,
      height_t = NA_real_,
      height_p = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      model_type = "insufficient_data",
      ci_method = NA_character_,
      p_method = NA_character_,
      fit_status = "skipped_insufficient_repeated_heights"
    ))
  }
  
  # Try mixed-effects model first
  mixed_try <- try({
    model <- lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = species_data)
    ct <- summary(model)$coefficients
    height_row <- which(rownames(ct) == "height_numeric")
    
    ci <- get_ci_for_height(model)
    
    tibble(
      species = species_name,
      n_obs = n_obs,
      n_trees = n_trees,
      n_trees_multi_height = n_trees_multi_height,
      height_coef = unname(ct[height_row, "Estimate"]),
      height_se   = unname(ct[height_row, "Std. Error"]),
      height_t    = unname(ct[height_row, "t value"]),
      # Large-sample z approx p-value for reference only
      height_p    = 2 * (1 - pnorm(abs(unname(ct[height_row, "t value"])))),
      ci_low      = ci$low,
      ci_high     = ci$high,
      model_type  = "mixed_effects",
      ci_method   = ci$ci_method,
      p_method    = "z_approx_from_lmer_t",
      fit_status  = "ok"
    )
  }, silent = TRUE)
  
  if (!inherits(mixed_try, "try-error")) return(mixed_try)
  
  # Fallback to simple linear model if mixed model fails
  lm_try <- try({
    model <- lm(CH4_best.flux ~ height_numeric, data = species_data)
    ct <- summary(model)$coefficients
    height_row <- which(rownames(ct) == "height_numeric")
    ci_mat <- confint(model, "height_numeric")
    
    tibble(
      species = species_name,
      n_obs = n_obs,
      n_trees = n_trees,
      n_trees_multi_height = n_trees_multi_height,
      height_coef = unname(ct[height_row, "Estimate"]),
      height_se   = unname(ct[height_row, "Std. Error"]),
      height_t    = unname(ct[height_row, "t value"]),
      height_p    = unname(ct[height_row, "Pr(>|t|)"]),
      ci_low      = as.numeric(ci_mat[1]),
      ci_high     = as.numeric(ci_mat[2]),
      model_type  = "linear_model",
      ci_method   = "t_interval_lm",
      p_method    = "t_test_lm",
      fit_status  = "mixed_failed_lm_ok"
    )
  }, silent = TRUE)
  
  if (!inherits(lm_try, "try-error")) return(lm_try)
  
  # If all fails
  tibble(
    species = species_name,
    n_obs = n_obs,
    n_trees = n_trees,
    n_trees_multi_height = n_trees_multi_height,
    height_coef = NA_real_,
    height_se = NA_real_,
    height_t = NA_real_,
    height_p = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    model_type = "failed",
    ci_method = NA_character_,
    p_method = NA_character_,
    fit_status = "all_failed"
  )
}

# Run per-species fits
species_list <- sort(unique(plot_data_top$species))
species_results <- map_dfr(species_list, test_species_height_effect, analysis_data)

# =========================
# Bottom panel data (CIs)
# =========================

plot_data_bottom <- species_results %>%
  filter(!is.na(height_coef)) %>%
  left_join(tree_counts %>% select(species, species_label), by = "species") %>%
  mutate(
    species_label = factor(species_label, levels = sort(unique(species_label))),
    # CI-based significance (exclude 0)
    significant = (!is.na(ci_low) & !is.na(ci_high)) & ((ci_low > 0) | (ci_high < 0)),
    negative_trend = ifelse(is.na(height_coef), FALSE, height_coef < 0),
    color_group = case_when(
      significant & negative_trend ~ "Significant Negative",
      significant & !negative_trend ~ "Significant Positive", 
      !significant & negative_trend ~ "Non-significant Negative",
      TRUE ~ "Non-significant Positive"
    )
  ) %>%
  arrange(height_coef)

# =========================
# Bottom panel (plot)
# =========================

p_bottom <- ggplot(plot_data_bottom, aes(x = reorder(species_label, height_coef),
                                         y = height_coef, color = color_group, shape = color_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3, alpha = 0.8) +
  scale_color_manual(values = colors_pub, name = "") +
  scale_shape_manual(values = shapes_pub, name = "") +
  labs(x = "", y = expression("Height effect")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8, face = "italic"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 50)
  )

# =========================
# Optional: individual top panels
# =========================

create_species_plot <- function(species_name, species_data, breaks_data) {
  current_breaks <- breaks_data$breaks_list[[which(breaks_data$species == species_name)]]
  current_label <- breaks_data$species_label[breaks_data$species == species_name]
  
  ggplot(species_data, aes(x = factor(height_m), y = CH4_best.flux)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.7) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, fill = boxplot_color, color = "gray30") +
    geom_jitter(alpha = 0.4, width = 0.2, size = 0.8, color = "gray20") +
    coord_flip() +
    scale_y_continuous(
      breaks = current_breaks,
      labels = function(x) {
        ifelse(abs(x) < 0.001, "0", 
               ifelse(abs(x) < 0.01, sprintf("%.3f", x),
                      ifelse(abs(x) < 0.1, sprintf("%.2f", x), 
                             sprintf("%.1f", x))))
      }
    ) +
    labs(x = "", y = "") +  # Remove individual axis labels
    ggtitle(current_label) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8, face = "italic"),
      plot.title = element_text(size = 8, face = "italic", margin = margin(b = 2)),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)  # Tighter margins
    )
}

plot_list <- list()
for(sp in species_list) {
  sp_data <- plot_data_top %>% filter(species == sp)
  if(nrow(sp_data) > 0) {
    plot_list[[sp]] <- create_species_plot(sp, sp_data, species_breaks)
  }
}

p_top_final <- grid.arrange(
  grobs = plot_list, 
  ncol = 3,
  left = textGrob("Height (m)", rot = 90, vjust = 1),
  bottom = textGrob(expression("CH"[4]*" flux (nmol m"^-2*" s"^-1*")"), vjust = 0),
  padding = unit(0.1, "line")
)

combined_plot <- grid.arrange(p_top_final, p_bottom, ncol = 1, heights = c(2.5, 1))
print(combined_plot)

# =========================
#     Understanding & QA
# =========================

# 1) Clean, readable per-species table with provenance
results_tidy <- species_results %>%
  left_join(tree_counts %>% select(species, species_label), by = "species") %>%
  transmute(
    species, species_label,
    model_type, fit_status, ci_method, p_method,
    n_obs, n_trees, n_trees_multi_height,
    estimate = height_coef,
    se = height_se,
    t_value = height_t,
    p_value = height_p,
    ci_low, ci_high,
    ci_width = ci_high - ci_low,
    ci_excludes_zero = if_else(!is.na(ci_low) & !is.na(ci_high) & (ci_low > 0 | ci_high < 0), TRUE, FALSE),
    direction = case_when(
      is.na(estimate) ~ NA_character_,
      estimate > 0 ~ "positive",
      estimate < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  arrange(desc(ci_excludes_zero), estimate)

cat("\n=== Per-species results (first 50 rows) ===\n")
print(results_tidy, n = 50)

# 2) Quick summaries
cat("\n=== Coverage ===\n")
coverage <- species_results %>%
  summarise(
    total_species = n_distinct(species),
    analyzed = sum(!is.na(height_coef)),
    skipped_insufficient = sum(model_type == "insufficient_data"),
    failed = sum(model_type == "failed")
  )
print(coverage)

cat("\n=== Model types used ===\n")
print(table(species_results$model_type, useNA = "ifany"))

cat("\n=== CI methods used ===\n")
print(table(results_tidy$ci_method, useNA = "ifany"))

cat("\n=== Fit statuses ===\n")
print(table(species_results$fit_status, useNA = "ifany"))

cat("\n=== CI-based significance counts by direction ===\n")
sig_counts <- results_tidy %>%
  filter(!is.na(estimate)) %>%
  count(ci_excludes_zero, direction)
print(sig_counts)

cat("\n=== Species excluded (insufficient repeated heights) ===\n")
excluded <- species_results %>% filter(model_type == "insufficient_data") %>%
  left_join(tree_counts %>% select(species, species_label), by = "species") %>%
  select(species, species_label, n_obs, n_trees, n_trees_multi_height, fit_status)
print(excluded, n = 50)

# 3) Ranked lists
cat("\n=== Top 5 most positive (by estimate) ===\n")
results_tidy %>%
  filter(!is.na(estimate)) %>%
  arrange(desc(estimate)) %>%
  head(5) %>%
  select(species, species_label, estimate, ci_low, ci_high, model_type, ci_method, fit_status) %>%
  print()

cat("\n=== Top 5 most negative (by estimate) ===\n")
results_tidy %>%
  filter(!is.na(estimate)) %>%
  arrange(estimate) %>%
  head(5) %>%
  select(species, species_label, estimate, ci_low, ci_high, model_type, ci_method, fit_status) %>%
  print()

# 4) Plain-language summaries
cat("\n=== Plain-language summaries (CI-based) ===\n")
summaries <- results_tidy %>%
  mutate(
    summary = case_when(
      is.na(estimate) ~ glue("{species_label}: not analyzed ({fit_status})."),
      ci_excludes_zero & estimate > 0 ~
        glue("{species_label}: CH4 flux increases with height; slope = {signif(estimate,3)} (95% CI {signif(ci_low,3)}, {signif(ci_high,3)}) [{model_type}, {ci_method}]."),
      ci_excludes_zero & estimate < 0 ~
        glue("{species_label}: CH4 flux decreases with height; slope = {signif(estimate,3)} (95% CI {signif(ci_low,3)}, {signif(ci_high,3)}) [{model_type}, {ci_method}]."),
      TRUE ~
        glue("{species_label}: no clear height effect; slope = {signif(estimate,3)} (95% CI {signif(ci_low,3)}, {signif(ci_high,3)}) [{model_type}, {ci_method}].")
    )
  ) %>%
  pull(summary)
cat(paste0("- ", summaries), sep = "\n")

# 5) CI width & precision checks
cat("\n=== CI width summary (narrower = more precise) ===\n")
results_tidy %>%
  filter(!is.na(ci_width)) %>%
  summarise(
    min_width = min(ci_width, na.rm = TRUE),
    q25 = quantile(ci_width, 0.25, na.rm = TRUE),
    median_width = median(ci_width, na.rm = TRUE),
    q75 = quantile(ci_width, 0.75, na.rm = TRUE),
    max_width = max(ci_width, na.rm = TRUE)
  ) %>% print()

cat("\n=== Relationship between precision and data volume ===\n")
prec_check <- results_tidy %>%
  filter(!is.na(ci_width)) %>%
  select(species_label, model_type, ci_method, n_obs, n_trees, n_trees_multi_height, ci_width) %>%
  arrange(ci_width)
print(prec_check, n = 50)

# 6) Export CSV
write_csv(results_tidy, "species_height_effects_ci_based_with_provenance.csv")
cat('\nSaved: "species_height_effects_ci_based_with_provenance.csv"\n')

# 7) Borderline cases (CI near zero)
cat("\n=== Borderline cases (CI near zero) ===\n")
near_zero <- results_tidy %>%
  filter(!is.na(ci_low), !is.na(ci_high)) %>%
  mutate(dist_to_zero = pmin(abs(ci_low), abs(ci_high))) %>%
  arrange(dist_to_zero) %>%
  head(10) %>%
  select(species, species_label, estimate, ci_low, ci_high, dist_to_zero, model_type, ci_method)
print(near_zero)

# 8) Per-species diagnostic helper (set code, then run)
inspect_species <- NA_character_  # e.g., "ACRU" to inspect
if (!is.na(inspect_species)) {
  cat(glue("\n=== Diagnostics for {inspect_species} ===\n"))
  this_dat <- analysis_data %>% filter(species == inspect_species)
  cat(glue("Obs: {nrow(this_dat)}, Trees: {n_distinct(this_dat$tree_unique)}\n"))
  has_repeats <- this_dat %>%
    group_by(tree_unique) %>%
    summarise(nh = n_distinct(height_numeric), .groups = "drop") %>%
    filter(nh > 1) %>% nrow() >= 3
  
  if (has_repeats) {
    m_try <- try(lmer(CH4_best.flux ~ height_numeric + (1|tree_unique), data = this_dat), silent = TRUE)
    if (!inherits(m_try, "try-error")) {
      print(summary(m_try))
      cat("\nProfile/Wald CI for slope:\n")
      print(try(confint(m_try, parm = "height_numeric", method = "profile", oldNames = FALSE), silent = TRUE))
      print(confint(m_try, parm = "height_numeric", method = "Wald", oldNames = FALSE))
    } else {
      m_lm <- lm(CH4_best.flux ~ height_numeric, data = this_dat)
      print(summary(m_lm))
      cat("\n95% CI for slope (lm):\n")
      print(confint(m_lm, "height_numeric"))
    }
  } else {
    cat("Insufficient repeated heights for mixed model diagnostics.\n")
  }
}
