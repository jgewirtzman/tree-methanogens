# ==============================================================================
# Soil and Tree Flux Time Series (Figure 1)
# ==============================================================================
# Purpose: Combined soil and tree flux faceted time series plots.
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - semirigid flux datasets (from data/processed/flux/)
#   - soilmoisture_total.csv (from data/raw/field_data/)
#
# Outputs:
#   - outputs/figures/main/fig1_temporal_flux_timeseries.png
# ==============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(cowplot)
library(scales)

# Read the datasets
soil_dataset <- read.csv('data/processed/flux/semirigid_tree_final_complete_dataset_soil.csv')
tree_dataset <- read.csv('data/processed/flux/semirigid_tree_final_complete_dataset.csv')
moisture_data <- read.csv('data/raw/field_data/ipad_data/Cleaned data/soilmoisture_total.csv')
names(moisture_data)[grep("Date", names(moisture_data))[1]] <- "Date"  # fix BOM-mangled column name

# Convert Date columns to proper date format
soil_dataset$Date <- as.Date(soil_dataset$Date)
tree_dataset$Date <- as.Date(tree_dataset$Date)

# ===== PREPARE DATA FOR CH4 FLUX PLOT (MIDDLE PANEL) =====

# Prepare soil data
soil_plot_data <- soil_dataset %>%
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2 >= 0.7 | CH4_HM.r2 >= 0.7),
         CH4_best.flux >= -100 & CH4_best.flux <= 200) %>%
  mutate(
    Plot_Type = case_when(
      Plot.letter %in% c("WD", "WS") ~ "W",
      Plot.letter == "I" ~ "I",
      Plot.letter == "U" ~ "U",
      TRUE ~ Plot.letter
    )
  ) %>%
  arrange(Date) %>%
  mutate(
    Date_numeric = as.numeric(Date),
    Date_group = cumsum(c(TRUE, diff(Date_numeric) > 7))
  ) %>%
  group_by(Date_group) %>%
  mutate(Date_interval = min(Date)) %>%
  ungroup() %>%
  group_by(Plot_Type, Date_interval) %>%
  summarise(
    mean_flux = mean(CH4_best.flux, na.rm = TRUE),
    se_flux = sd(CH4_best.flux, na.rm = TRUE) / sqrt(n()),
    n_measurements = n(),
    .groups = 'drop'
  ) %>%
  filter(!is.na(mean_flux), !is.na(Date_interval)) %>%
  mutate(Data_Type = "Soil")

# Prepare tree data
tree_plot_data <- tree_dataset %>%
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2.x >= 0.7 | CH4_HM.r2.x >= 0.7),
         CH4_best.flux.x >= -100 & CH4_best.flux.x <= 200) %>%
  mutate(
    Plot_Type = case_when(
      Plot.Letter %in% c("WD", "WS") ~ "W",
      Plot.Letter == "I" ~ "I",
      Plot.Letter == "U" ~ "U",
      TRUE ~ Plot.Letter
    )
  ) %>%
  arrange(Date) %>%
  mutate(
    Date_numeric = as.numeric(Date),
    Date_group = cumsum(c(TRUE, diff(Date_numeric) > 7))
  ) %>%
  group_by(Date_group) %>%
  mutate(Date_interval = min(Date)) %>%
  ungroup() %>%
  group_by(Plot_Type, Date_interval) %>%
  summarise(
    mean_flux = mean(CH4_best.flux.x, na.rm = TRUE),
    se_flux = sd(CH4_best.flux.x, na.rm = TRUE) / sqrt(n()),
    n_measurements = n(),
    .groups = 'drop'
  ) %>%
  filter(!is.na(mean_flux), !is.na(Date_interval)) %>%
  mutate(Data_Type = "Tree")

# Combine datasets for CH4 flux
combined_facet_data <- bind_rows(soil_plot_data, tree_plot_data) %>%
  mutate(Plot_Type = factor(Plot_Type, levels = c("U", "I", "W"),
                            labels = c("Upland", "Intermediate", "Wetland")))

# Prepare raw data for jittering
soil_raw_data <- soil_dataset %>%
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2 >= 0.7 | CH4_HM.r2 >= 0.7),
         CH4_best.flux >= -100 & CH4_best.flux <= 200) %>%
  mutate(
    Plot_Type = case_when(
      Plot.letter %in% c("WD", "WS") ~ "W",
      Plot.letter == "I" ~ "I",
      Plot.letter == "U" ~ "U",
      TRUE ~ Plot.letter
    ),
    Data_Type = "Soil",
    Plot.Tag = as.character(Plot.Tag)
  ) %>%
  arrange(Date) %>%
  mutate(
    Date_numeric = as.numeric(Date),
    Date_group = cumsum(c(TRUE, diff(Date_numeric) > 7))
  ) %>%
  group_by(Date_group) %>%
  mutate(Date_interval = min(Date)) %>%
  ungroup() %>%
  filter(!is.na(CH4_best.flux), !is.na(Date_interval)) %>%
  rename(CH4_flux = CH4_best.flux)

tree_raw_data <- tree_dataset %>%
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2.x >= 0.7 | CH4_HM.r2.x >= 0.7),
         CH4_best.flux.x >= -100 & CH4_best.flux.x <= 200) %>%
  mutate(
    Plot_Type = case_when(
      Plot.Letter %in% c("WD", "WS") ~ "W",
      Plot.Letter == "I" ~ "I",
      Plot.Letter == "U" ~ "U",
      TRUE ~ Plot.Letter
    ),
    Data_Type = "Tree",
    Plot.Tag = as.character(Plot.Tag)
  ) %>%
  arrange(Date) %>%
  mutate(
    Date_numeric = as.numeric(Date),
    Date_group = cumsum(c(TRUE, diff(Date_numeric) > 7))
  ) %>%
  group_by(Date_group) %>%
  mutate(Date_interval = min(Date)) %>%
  ungroup() %>%
  filter(!is.na(CH4_best.flux.x), !is.na(Date_interval)) %>%
  rename(CH4_flux = CH4_best.flux.x)

combined_facet_raw <- bind_rows(soil_raw_data, tree_raw_data) %>%
  mutate(Plot_Type = factor(Plot_Type, levels = c("U", "I", "W"),
                            labels = c("Upland", "Intermediate", "Wetland")))

# ===== PREPARE DATA FOR VWC DENSITY PLOT (TOP PANEL) =====
vwc_data <- moisture_data %>%
  mutate(
    Date = mdy(Date),
    Plot_Type = case_when(
      Plot.letter == "U" ~ "Upland",
      Plot.letter == "I" ~ "Intermediate",
      Plot.letter %in% c("WS", "WD") ~ "Wetland"
    ),
    Plot_Type = factor(Plot_Type, levels = c("Upland", "Intermediate", "Wetland")),
    Soil_Moisture = as.numeric(VWC)
  ) %>%
  filter(!is.na(Soil_Moisture))

# ===== PREPARE DATA FOR SOIL TEMP/MOISTURE PLOT (BOTTOM PANEL) =====
soil_temp_moisture_data <- moisture_data %>%
  mutate(
    Date = mdy(Date),
    Plot_Type = case_when(
      Plot.letter == "U" ~ "Upland",
      Plot.letter == "I" ~ "Intermediate",
      Plot.letter %in% c("WS", "WD") ~ "Wetland"
    ),
    Plot_Type = factor(Plot_Type, levels = c("Upland", "Intermediate", "Wetland")),
    Soil_Temperature = as.numeric(Soil.temp),
    Soil_Moisture = as.numeric(VWC)
  ) %>%
  filter(!is.na(Date), !is.na(Soil_Temperature), !is.na(Soil_Moisture))

# ===== CREATE VWC DENSITY PLOT (TOP) =====
vwc_plot <- ggplot(vwc_data, aes(x = Soil_Moisture, fill = Plot_Type)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ Plot_Type, strip.position = "top") +
  scale_fill_brewer(type = "seq", palette = "Blues") +
  scale_x_continuous(position = "bottom") +
  labs(
    x = "Volumetric Water Content (%)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_blank(),
    axis.title.x = element_text(size = 8.8),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x.bottom = element_line(color = "black", linewidth = 0.3),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5, unit = "pt"),
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
  scale_color_manual(values = c("Soil" = "#8B4513", "Tree" = "#228B22")) +
  scale_fill_manual(values = c("Soil" = "#8B4513", "Tree" = "#228B22")) +
  labs(
    y = bquote("CH"[4] ~ "Flux (nmol m"^-2 ~ "s"^-1 ~ ")")
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 1),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"),
    panel.spacing = unit(1, "lines")
  ) +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.005),
    breaks = c(-10, -1, -0.1, 0, 0.1, 1, 10, 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.9)

# ===== CREATE SOIL TEMP/MOISTURE PLOT (BOTTOM) =====
# Uses geom_smooth(se = TRUE) for identical curves to original + SE bands

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
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"),
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

# Display the final combined plot
print(final_plot)

# Save the plot
ggsave("outputs/figures/main/fig1_temporal_flux_timeseries.png",
       final_plot, width = 10, height = 7, dpi = 300)

# ===== CALCULATE STATISTICS FOR FIGURE CAPTION =====

# Overall statistics
cat("\n===== OVERALL STATISTICS =====\n")
cat("Total tree flux measurements:", nrow(tree_dataset %>%
                                            filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2.x >= 0.7 | CH4_HM.r2.x >= 0.7),
                                                   CH4_best.flux.x >= -100 & CH4_best.flux.x <= 200)), "\n")
cat("Total soil flux measurements:", nrow(soil_dataset %>%
                                            filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2 >= 0.7 | CH4_HM.r2 >= 0.7),
                                                   CH4_best.flux >= -100 & CH4_best.flux <= 200)), "\n")
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

# Convert nmol m-2 s-1 to ug CH4 m-2 hr-1
# Factor = 57.744
conversion_factor <- 57.744

flux_converted <- flux_stats %>%
  mutate(
    mean_ug = mean * conversion_factor,
    median_ug = median * conversion_factor,
    min_ug = min * conversion_factor,
    max_ug = max * conversion_factor,
    q25_ug = q25 * conversion_factor,
    q75_ug = q75 * conversion_factor
  )

cat("\nUpland tree emissions (ug CH4 m-2 hr-1):\n")
upland_tree <- flux_converted %>% filter(Plot_Type == "Upland", Data_Type == "Tree")
cat("  Mean:", round(upland_tree$mean_ug, 1), "\n")
cat("  Range:", round(upland_tree$min_ug, 1), "to", round(upland_tree$max_ug, 1), "\n")
cat("  IQR:", round(upland_tree$q25_ug, 1), "to", round(upland_tree$q75_ug, 1), "\n")

cat("\nUpland soil consumption (ug CH4 m-2 hr-1):\n")
upland_soil <- flux_converted %>% filter(Plot_Type == "Upland", Data_Type == "Soil")
cat("  Mean:", round(upland_soil$mean_ug, 1), "\n")
cat("  Range:", round(upland_soil$min_ug, 1), "to", round(upland_soil$max_ug, 1), "\n")

cat("\nIntermediate tree emissions (ug CH4 m-2 hr-1):\n")
int_tree <- flux_converted %>% filter(Plot_Type == "Intermediate", Data_Type == "Tree")
cat("  Mean:", round(int_tree$mean_ug, 1), "\n")
cat("  Range:", round(int_tree$min_ug, 1), "to", round(int_tree$max_ug, 1), "\n")

cat("\nIntermediate soil consumption (ug CH4 m-2 hr-1):\n")
int_soil <- flux_converted %>% filter(Plot_Type == "Intermediate", Data_Type == "Soil")
cat("  Mean:", round(int_soil$mean_ug, 1), "\n")
cat("  Range:", round(int_soil$min_ug, 1), "to", round(int_soil$max_ug, 1), "\n")

cat("\nTransitional wetland tree emissions (ug CH4 m-2 hr-1):\n")
wet_tree <- flux_converted %>% filter(Plot_Type == "Wetland", Data_Type == "Tree")
cat("  Mean:", round(wet_tree$mean_ug, 1), "\n")
cat("  Range:", round(wet_tree$min_ug, 1), "to", round(wet_tree$max_ug, 1), "\n")

cat("\nTransitional wetland soil fluxes (ug CH4 m-2 hr-1):\n")
wet_soil <- flux_converted %>% filter(Plot_Type == "Wetland", Data_Type == "Soil")
cat("  Mean:", round(wet_soil$mean_ug, 1), "\n")
cat("  Median:", round(wet_soil$median_ug, 1), "\n")
cat("  Range:", round(wet_soil$min_ug, 1), "to", round(wet_soil$max_ug, 1), "\n")

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
