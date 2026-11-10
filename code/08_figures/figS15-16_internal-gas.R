# ==============================================================================
# figS15-16_internal-gas.R -- Figures S15 and S16: internal stem gases
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root. SI house style: code/lib/fig_style.R.
# Figure S15  internal CH4 by species (quasirandom points coloured by internal O2;
#             mean +/- SE), species ordered by mean CH4, pseudo-log axis
# Figure S16  (a) CO2 against CH4, (b) O2 against CH4, (c) heartwood mcrA against CH4,
#             (d) stem flux at 50, 125 and 200 cm against CH4; linear fits with R2 and
#             p; points coloured broadleaf / conifer (GROUP_COLS)
# Inputs: data/tables/ via code/lib/tables_read.R: survey_2021_trees() (internal_gas
#   2021 whole-stem samples, flux_stem 2021 heights, ddpcr + samples, trees, species)
# Outputs: outputs/figures/original/supplementary/figS7_internal_gas_beeswarm.png  (S15)
#          outputs/figures/original/supplementary/figS8_internal_gas_profiles.png  (S16)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(ggbeeswarm); library(patchwork); library(broom)
})
source("code/lib/tables_read.R")
source("code/lib/fig_style.R")
tree_all <- survey_2021_trees()   # 2021 survey per tree; species from trees.csv

# Axis labels: ppm below 1%, % above ("0", "10 ppm", "1,000 ppm", "1%", "10%")
format_mixed_labels <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (v == 0) return("0")
    if (v < 10000) return(paste(format(v, big.mark = ",", scientific = FALSE, trim = TRUE), "ppm"))
    paste0(format(v / 10000, scientific = FALSE, trim = TRUE, drop0trailing = TRUE), "%")
  }, character(1))
}
pct_labels <- function(x) paste0(format(x / 10000, trim = TRUE, drop0trailing = TRUE), "%")

# statistics annotation, one line: "R^2 = 0.22, p < 0.001" (plotmath)
stat_label <- function(r2, p) {
  r2s <- formatC(signif(r2, 2), digits = 2, format = "fg", flag = "#")
  ps  <- ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  sprintf("R^2 == '%s' * ',' ~ '%s'", r2s, ps)
}

# ==============================================================================
# Figure S15: internal CH4 by species
# ==============================================================================
ridge_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(O2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_name(species_id)) %>%
  filter(!is.na(Species_Latin)) %>%
  mutate(Species_Latin = reorder(Species_Latin, CH4_concentration, mean, na.rm = TRUE))

species_stats <- ridge_data %>%
  group_by(Species_Latin) %>%
  summarise(
    mean_ch4 = mean(CH4_concentration, na.rm = TRUE),
    se_ch4 = sd(CH4_concentration, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

p_beeswarm <- ggplot(ridge_data, aes(x = Species_Latin, y = CH4_concentration)) +
  geom_quasirandom(aes(color = O2_concentration),
                   size = 1.8, alpha = 0.75, width = 0.3) +
  geom_pointrange(data = species_stats,
                  aes(x = Species_Latin, y = mean_ch4,
                      ymin = mean_ch4 - se_ch4, ymax = mean_ch4 + se_ch4),
                  color = "black", size = 0.35, linewidth = 0.8) +
  scale_color_gradientn(
    colors = rev(c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")),
    values = scales::rescale(c(0, 100000, 160000, 190000, 200000, 220000)),
    name = expression(Internal ~ O[2] ~ "(%)"),
    labels = function(x) format(x / 10000, trim = TRUE, drop0trailing = TRUE),
    guide = guide_colorbar(barheight = unit(10, "lines"), barwidth = unit(0.8, "lines"),
                           ticks.colour = "black")
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(0, 10, 100, 1000, 10000, 100000),
    labels = format_mixed_labels
  ) +
  coord_flip() +
  labs(y = expression(CH[4] ~ concentration), x = NULL) +
  theme_si() +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major.y = element_line(colour = "grey94"))
si_save("outputs/figures/original/supplementary/figS7_internal_gas_beeswarm.png",
        p_beeswarm, height = 6)

# ==============================================================================
# Figure S16: internal gases, mcrA and flux
# ==============================================================================
# Prepare data for CH4 conc vs CO2 conc
conc_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(CO2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_name(species_id)) %>%
  filter(!is.na(Species_Latin))

# Prepare data for CH4 conc vs O2 conc
o2_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(O2_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_name(species_id)) %>%
  filter(!is.na(Species_Latin))

# Prepare data for heartwood mcrA vs CH4 conc
mcra_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(ddpcr_mcra_probe_Inner_loose) & !is.na(species_id)) %>%
  mutate(
    Species_Latin = species_name(species_id),
    mcra_copies = ddpcr_mcra_probe_Inner_loose + 1
  ) %>%
  filter(!is.na(Species_Latin))

# Prepare data for CH4 flux vs CH4 conc
flux_conc_data <- tree_all %>%
  filter(!is.na(CH4_concentration) & !is.na(species_id)) %>%
  mutate(Species_Latin = species_name(species_id)) %>%
  filter(!is.na(Species_Latin)) %>%
  pivot_longer(
    cols = c(CH4_best.flux_50cm, CH4_best.flux_125cm, CH4_best.flux_200cm),
    names_to = "height",
    values_to = "CH4_flux",
    values_drop_na = TRUE
  ) %>%
  mutate(
    height_clean = case_when(
      height == "CH4_best.flux_50cm" ~ "50 cm",
      height == "CH4_best.flux_125cm" ~ "125 cm",
      height == "CH4_best.flux_200cm" ~ "200 cm",
      TRUE ~ height
    ),
    height_clean = factor(height_clean, levels = c("50 cm", "125 cm", "200 cm"))
  )

# broadleaf / conifer grouping for the point colours
grp <- function(d) mutate(d, Group = factor(species_group(as.character(Species_Latin)),
                                            levels = names(GROUP_COLS)))
conc_data <- grp(conc_data); o2_data <- grp(o2_data)
mcra_data <- grp(mcra_data); flux_conc_data <- grp(flux_conc_data)

# Calculate correlations - CO2
pearson1 <- cor.test(log10(conc_data$CH4_concentration + 1),
                     log10(conc_data$CO2_concentration + 1), method = "pearson")
model1 <- lm(log10(CH4_concentration + 1) ~ log10(CO2_concentration + 1), data = conc_data)
stats1 <- glance(model1)
label1 <- stat_label(stats1$r.squared, pearson1$p.value)

# Calculate correlations - O2
pearson2 <- cor.test(log10(o2_data$CH4_concentration + 1),
                     log10(o2_data$O2_concentration + 1), method = "pearson")
model2 <- lm(log10(CH4_concentration + 1) ~ log10(O2_concentration + 1), data = o2_data)
stats2 <- glance(model2)
label2 <- stat_label(stats2$r.squared, pearson2$p.value)

# Calculate correlations - mcrA
pearson_mcra <- cor.test(log10(mcra_data$mcra_copies),
                         log10(mcra_data$CH4_concentration + 1), method = "pearson")
model_mcra <- lm(log10(CH4_concentration + 1) ~ log10(mcra_copies), data = mcra_data)
stats_mcra <- glance(model_mcra)
label_mcra <- stat_label(stats_mcra$r.squared, pearson_mcra$p.value)

# Calculate stats for each height
flux_stats <- flux_conc_data %>%
  group_by(height_clean) %>%
  do({
    data <- .
    model <- lm(log10(abs(CH4_flux) + 0.01) ~ log10(CH4_concentration + 1), data = data)
    pearson_test <- cor.test(log10(data$CH4_concentration + 1),
                             log10(abs(data$CH4_flux) + 0.01), method = "pearson")
    stats <- glance(model)
    data.frame(
      height_clean = unique(data$height_clean),
      label = stat_label(stats$r.squared, pearson_test$p.value)
    )
  }) %>%
  ungroup()
cat("S16 annotations:", label1, "|", label2, "|", label_mcra, "|",
    paste(flux_stats$label, collapse = " | "), "\n")

# ---- shared layers ----
ch4_breaks <- c(0, 10, 100, 1000, 10000, 100000)
pt_layer  <- function() geom_point(aes(color = Group), size = 1.5, alpha = 0.7)
fit_layer <- function() geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                                    color = "black", linewidth = 0.8)
grp_scale <- function() scale_color_manual(values = GROUP_COLS, name = NULL, drop = FALSE)
stat_text <- function(lab, x = -Inf, hjust = -0.05)
  annotate("text", x = x, y = Inf, label = lab, hjust = hjust, vjust = 1.4,
           size = 3, parse = TRUE)
y_ch4 <- scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                            breaks = ch4_breaks, labels = format_mixed_labels,
                            expand = expansion(mult = c(0.05, 0.14)))

p1 <- ggplot(conc_data, aes(x = CO2_concentration, y = CH4_concentration)) +
  pt_layer() + fit_layer() + stat_text(label1) + grp_scale() +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                     breaks = c(1000, 10000, 100000),
                     minor_breaks = c(2000, 5000, 20000, 50000, 200000),
                     labels = format_mixed_labels) +
  y_ch4 +
  labs(x = expression(CO[2] ~ concentration), y = expression(CH[4] ~ concentration)) +
  theme_si()

p2 <- ggplot(o2_data, aes(x = O2_concentration, y = CH4_concentration)) +
  pt_layer() + fit_layer() + stat_text(label2) + grp_scale() +
  scale_x_continuous(breaks = seq(0, 220000, by = 50000), labels = pct_labels) +
  y_ch4 +
  labs(x = expression(O[2] ~ concentration), y = NULL) +
  theme_si()

p3 <- ggplot(mcra_data, aes(x = mcra_copies, y = CH4_concentration)) +
  pt_layer() + fit_layer() + stat_text(label_mcra, x = 0) + grp_scale() +   # x = 0 is -Inf on the log axis
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  y_ch4 +
  labs(x = expression(atop(Heartwood ~ italic(mcrA), (copies ~ g^{-1}))), y = NULL) +
  theme_si()

p4 <- ggplot(flux_conc_data, aes(x = CH4_concentration, y = CH4_flux)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  pt_layer() + fit_layer() +
  geom_text(data = flux_stats, aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.4, size = 3, inherit.aes = FALSE, parse = TRUE) +
  facet_wrap(~ height_clean, ncol = 3) +
  grp_scale() +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 1),
                     breaks = c(10, 1000, 100000), minor_breaks = c(100, 10000),
                     labels = format_mixed_labels) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.01),
                     breaks = c(-0.1, 0, 0.1, 0.5, 1, 2, 5),
                     labels = c("−0.1", "0", "0.1", "0.5", "1", "2", "5"),
                     expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = expression(CH[4] ~ concentration), y = lab_flux) +
  theme_si() + theme(panel.spacing.x = unit(1, "lines"))

# si_tags() in fig_style.R returns NULL (plot_annotation & theme is not a patchwork),
# so its tag settings are applied here directly: bold "(a)" at SI_BASE + 1
final_plot <- ((p1 | p2 | p3) / p4) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = SI_BASE + 1, face = "bold"),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0))

si_save("outputs/figures/original/supplementary/figS8_internal_gas_profiles.png",
        final_plot, height = 6.6)
