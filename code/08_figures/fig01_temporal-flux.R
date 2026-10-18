# ==============================================================================
# fig01_temporal-flux.R -- Figure 1: soil and stem CH4 flux through the monthly survey
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) soil VWC at the collars, density by landscape position
#   (b) soil and stem CH4 flux by sampling campaign: deployments (jittered) and the
#       campaign mean +/- SE, arcsinh axis (code/lib/flux_transform.R)
#   (c) soil temperature and VWC through time (loess +/- SE)
# Inputs (data/tables/, via code/lib/tables_read.R)
#   flux_soil.csv, flux_stem.csv  monthly_2020_2021 deployments. A deployment is left out
#                                 only for procedural failure: soil chambers that failed
#                                 the closure screen (flux_soil.qc_pass). No r2 or range
#                                 filter. Stems: every deployment, untagged and dead
#                                 stems included.
#   collars.csv, trees.csv        landscape position of each collar and stem
#   soil_readings.csv             collar readings, one per collar and date
#                                 (code/lib/soil_collars.R collar_soil_readings())
# Output
#   outputs/figures/generated/fig1_final.png; summary statistics to the console
# ==============================================================================

library(ggplot2)
library(dplyr)
library(lubridate)
library(cowplot)
library(scales)
source("code/lib/tables_read.R")
source("code/lib/flux_transform.R")
source("code/lib/soil_collars.R")

collars <- read_table("collars"); trees <- trees_table()
soil_dataset <- read_table("flux_soil") %>% filter(campaign == "monthly_2020_2021", qc_pass) %>%
  transmute(UniqueID = measurement_id, Date = as.Date(date), Plot.Tag = collar_id,
            Plot.letter = collars$landscape_position[match(collar_id, collars$collar_id)],
            CH4_best.flux = CH4_flux_nmol_m2_s)
tree_dataset <- read_table("flux_stem") %>% filter(campaign == "monthly_2020_2021") %>%
  transmute(UniqueID = measurement_id, Date = as.Date(date), Plot.Tag = tree_id,
            Plot.letter = trees$landscape_position[match(tree_id, trees$tree_id)],
            CH4_best.flux = CH4_flux_nmol_m2_s)
stopifnot(!anyNA(soil_dataset$Plot.letter), !anyNA(tree_dataset$Plot.letter),
          !anyNA(soil_dataset$Date), !anyNA(tree_dataset$Date))
moisture_data <- collar_soil_readings()   # one row per collar and date

POSITIONS <- c(U = "Upland", I = "Intermediate", W = "Wetland")   # WD and WS pooled as W

# ===== PREPARE DATA FOR CH4 FLUX PLOT (MIDDLE PANEL) =====
# Sampling campaigns: a date more than 7 days after the previous one starts a new
# campaign, dated by its first day (soil and stems grouped separately).
with_campaign <- function(d, type) d %>%
  filter(!is.na(CH4_best.flux)) %>%
  mutate(Plot_Type = ifelse(Plot.letter %in% c("WD", "WS"), "W", Plot.letter),
         Data_Type = type, Plot.Tag = as.character(Plot.Tag)) %>%
  arrange(Date) %>%
  mutate(Date_group = cumsum(c(TRUE, diff(as.numeric(Date)) > 7))) %>%
  group_by(Date_group) %>%
  mutate(Date_interval = min(Date)) %>%
  ungroup() %>%
  rename(CH4_flux = CH4_best.flux)
soil_raw_data <- with_campaign(soil_dataset, "Soil")
tree_raw_data <- with_campaign(tree_dataset, "Tree")

combined_facet_raw <- bind_rows(soil_raw_data, tree_raw_data) %>%
  mutate(Plot_Type = factor(POSITIONS[Plot_Type], levels = POSITIONS))
# campaign mean and SE by data type and landscape position
combined_facet_data <- combined_facet_raw %>%
  group_by(Data_Type, Plot_Type, Date_interval) %>%
  summarise(mean_flux = mean(CH4_flux), se_flux = sd(CH4_flux) / sqrt(n()),
            n_measurements = n(), .groups = "drop")

# ===== PREPARE DATA FOR VWC DENSITY (TOP) AND SOIL TEMP/MOISTURE (BOTTOM) PANELS =====
collar_position <- function(letter) factor(POSITIONS[ifelse(letter %in% c("WS", "WD"), "W", letter)],
                                           levels = POSITIONS)
vwc_data <- moisture_data %>%
  mutate(Date = date, Plot_Type = collar_position(plot_letter), Soil_Moisture = soil_vwc_pct) %>%
  filter(!is.na(Soil_Moisture))
soil_temp_moisture_data <- moisture_data %>%
  mutate(Date = date, Plot_Type = collar_position(plot_letter),
         Soil_Temperature = soil_temp_C, Soil_Moisture = soil_vwc_pct) %>%
  filter(!is.na(Date), !is.na(Soil_Temperature), !is.na(Soil_Moisture))

# ===== CREATE VWC DENSITY PLOT (TOP) =====
vwc_plot <- ggplot(vwc_data, aes(x = Soil_Moisture, fill = Plot_Type)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ Plot_Type, strip.position = "top") +
  scale_fill_brewer(type = "seq", palette = "Blues") +
  scale_x_continuous(position = "bottom") +
  geom_text(data = data.frame(Plot_Type = factor("Upland",
              levels = c("Upland", "Intermediate", "Wetland")),
              x = 10, y = -Inf),
            aes(x = x, y = y), label = "VWC (%)", size = 3,
            hjust = 0.5, vjust = -0.5) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x.bottom = element_line(color = "black", linewidth = 0.3),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 16, unit = "pt"),
    strip.placement = "outside",
    panel.spacing = unit(0.1, "lines")
  ) +
  guides(fill = guide_legend(title = NULL))

# ===== CREATE CH4 FLUX PLOT (MIDDLE) =====
ch4_flux_plot <- ggplot() +
  geom_jitter(data = combined_facet_raw,
              aes(x = Date_interval, y = CH4_flux, color = Data_Type),
              alpha = 0.3, size = 1, width = 8, height = 0) +
  geom_errorbar(data = combined_facet_data,
                aes(x = Date_interval, y = mean_flux, color = Data_Type,
                    ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 15, alpha = 1, position = position_dodge(width = 10)) +
  geom_point(data = combined_facet_data,
             aes(x = Date_interval, y = mean_flux, color = Data_Type, fill = Data_Type),
             size = 3, alpha = 0.6, pch = 21, position = position_dodge(width = 10)) +
  geom_line(data = combined_facet_data,
            aes(x = Date_interval, y = mean_flux, color = Data_Type, group = Data_Type),
            linewidth = .8, alpha = 0.3) +
  facet_wrap(~ Plot_Type, nrow = 1) +
  scale_color_manual(values = c("Soil" = "#8B4513", "Tree" = "#228B22"), labels = c(Soil = "Soil", Tree = "Stem")) +
  scale_fill_manual(values = c("Soil" = "#8B4513", "Tree" = "#228B22"), labels = c(Soil = "Soil", Tree = "Stem")) +
  labs(
    y = bquote(bold("CH"[4] ~ "flux (nmol m"^-2 ~ "s"^-1 * ")"))
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 1),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 16, unit = "pt"),
    panel.spacing = unit(1, "lines")
  ) +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_y_continuous(
    trans = flux_asinh_trans(),
    breaks = c(-10, -1, -0.1, 0, 0.1, 1, 10, 100),
    labels = c("−10", "−1", "−0.1", "0", "0.1", "1", "10", "100")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.9)

# ===== CREATE SOIL TEMP/MOISTURE PLOT (BOTTOM) =====
# loess +/- SE; VWC drawn on the temperature axis divided by 2.5 (secondary axis)

temp_moisture_plot <- ggplot(soil_temp_moisture_data, aes(x = Date)) +
  # Temperature with SE band
  geom_smooth(aes(y = Soil_Temperature, fill = "Temperature", color = "Temperature"),
              method = "loess", se = TRUE, alpha = 0.2, linewidth = 0.75) +
  # Moisture with SE band (scaled to temperature axis)
  geom_smooth(aes(y = Soil_Moisture / 2.5, fill = "Moisture", color = "Moisture"),
              method = "loess", se = TRUE, alpha = 0.2, linewidth = 0.75) +
  facet_wrap(~ Plot_Type) +
  scale_y_continuous(
    name = expression("Soil T (" * degree * "C)"),
    sec.axis = sec_axis(~ . * 2.5, name = "VWC (%)")
  ) +
  scale_color_manual(values = c("Temperature" = alpha("#D73027", 0.5),
                                "Moisture" = alpha("#4575B4", 0.5))) +
  scale_fill_manual(values = c("Temperature" = "#D73027", "Moisture" = "#4575B4")) +
  labs(x = "Date") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 16, unit = "pt"),
    panel.spacing = unit(1, "lines"),
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt"),
    legend.box.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  xlab(NULL) +
  guides(color = guide_legend(title = NULL), fill = "none") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months")

# ===== COMBINE ALL PLOTS =====
final_plot <- plot_grid(
  vwc_plot,
  ch4_flux_plot,
  temp_moisture_plot,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(.75, 1.5, 1),
  labels = c("(a)", "(b)", "(c)"),
  label_size = 11,
  label_fontface = "bold"
)

# Save the plot
ggsave("outputs/figures/generated/fig1_final.png",
       final_plot, width = 8, height = 6, dpi = 300, bg = "white")   # near print size (17 cm) so text prints at ~7-9 pt

# ===== CALCULATE STATISTICS FOR FIGURE CAPTION =====

# Overall statistics
cat("\n===== OVERALL STATISTICS =====\n")
cat("Total tree flux measurements:", nrow(tree_raw_data), "\n")
cat("Total soil flux measurements:", nrow(soil_raw_data), "\n")
cat("Date range:", format(min(combined_facet_data$Date_interval, na.rm = TRUE), "%B %Y"),
    "to", format(max(combined_facet_data$Date_interval, na.rm = TRUE), "%B %Y"), "\n")

# CH4 flux statistics by landscape position
cat("\n===== CH4 FLUX BY LANDSCAPE POSITION (nmol m-2 s-1) =====\n")
flux_stats <- combined_facet_raw %>%
  group_by(Plot_Type, Data_Type) %>%
  summarise(
    mean = mean(CH4_flux, na.rm = TRUE),
    median = median(CH4_flux, na.rm = TRUE),
    min = min(CH4_flux, na.rm = TRUE),
    max = max(CH4_flux, na.rm = TRUE),
    q25 = quantile(CH4_flux, 0.25, na.rm = TRUE),
    q75 = quantile(CH4_flux, 0.75, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )
print(flux_stats)

# VWC statistics by landscape position
cat("\n===== VOLUMETRIC WATER CONTENT BY LANDSCAPE POSITION (%) =====\n")
vwc_stats <- vwc_data %>%
  group_by(Plot_Type) %>%
  summarise(
    mean = mean(Soil_Moisture, na.rm = TRUE),
    sd = sd(Soil_Moisture, na.rm = TRUE),
    median = median(Soil_Moisture, na.rm = TRUE),
    min = min(Soil_Moisture, na.rm = TRUE),
    max = max(Soil_Moisture, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )
print(vwc_stats)

# Number of unique trees and soil plots
cat("\n===== SAMPLING UNITS =====\n")
cat("Unique trees sampled:", length(unique(tree_raw_data$Plot.Tag)), "\n")
cat("Unique soil plots sampled:", length(unique(soil_raw_data$Plot.Tag)), "\n")

# Number of measurement dates
cat("\n===== TEMPORAL COVERAGE =====\n")
cat("Number of measurement campaigns (grouped by 7-day intervals):",
    length(unique(combined_facet_data$Date_interval)), "\n")
cat("Tree measurement dates:", length(unique(tree_raw_data$Date_interval)), "\n")
cat("Soil measurement dates:", length(unique(soil_raw_data$Date_interval)), "\n")

# Temperature ranges
cat("\n===== SOIL TEMPERATURE RANGES =====\n")
temp_stats <- soil_temp_moisture_data %>%
  group_by(Plot_Type) %>%
  summarise(
    temp_min = min(Soil_Temperature, na.rm = TRUE),
    temp_max = max(Soil_Temperature, na.rm = TRUE),
    temp_mean = mean(Soil_Temperature, na.rm = TRUE),
    .groups = 'drop'
  )
print(temp_stats)


# ===== ADDITIONAL STATISTICS FOR RESULTS TEXT =====

cat("\n===== FLUX STATISTICS IN DIFFERENT UNITS =====\n")
# nmol m-2 s-1 -> ug CH4 m-2 h-1: x 16.04 g mol-1 x 3600 s h-1 / 1000
conversion_factor <- 57.744
flux_converted <- flux_stats %>%
  mutate(across(c(mean, median, min, max, q25, q75), ~ .x * conversion_factor, .names = "{.col}_ug"))
for (pos in POSITIONS) for (type in c("Tree", "Soil")) {
  r <- flux_converted %>% filter(Plot_Type == pos, Data_Type == type)
  cat(sprintf("\n%s %s (ug CH4 m-2 hr-1):\n", pos, if (type == "Tree") "tree" else "soil"))
  cat("  Mean:", round(r$mean_ug, 1), "\n")
  cat("  Median:", round(r$median_ug, 1), "\n")
  cat("  Range:", round(r$min_ug, 1), "to", round(r$max_ug, 1), "\n")
  cat("  IQR:", round(r$q25_ug, 1), "to", round(r$q75_ug, 1), "\n")
}

# Check proportion of positive vs negative fluxes
cat("\n===== FLUX DIRECTIONALITY =====\n")
direction_stats <- combined_facet_raw %>%
  group_by(Plot_Type, Data_Type) %>%
  summarise(
    n_positive = sum(CH4_flux > 0, na.rm = TRUE),
    n_negative = sum(CH4_flux < 0, na.rm = TRUE),
    n_zero = sum(CH4_flux == 0, na.rm = TRUE),
    pct_positive = n_positive / n() * 100,
    .groups = 'drop'
  )
print(direction_stats)

# Statistical significance tests
cat("\n===== STATISTICAL COMPARISONS =====\n")
cat("Testing if upland tree fluxes are significantly > 0:\n")
upland_trees_only <- combined_facet_raw %>%
  filter(Plot_Type == "Upland", Data_Type == "Tree")
t_test_result <- t.test(upland_trees_only$CH4_flux, mu = 0, alternative = "greater")
cat("  t-statistic:", round(t_test_result$statistic, 3), "\n")
cat("  p-value:", format.pval(t_test_result$p.value), "\n")
