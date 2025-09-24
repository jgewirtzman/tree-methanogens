# Load required libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(scales)

# Read the datasets
soil_dataset <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/semirigid_tree_final_complete_dataset_soil.csv')
tree_dataset <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/untitled folder/semirigid_tree_final_complete_dataset.csv')

# Convert Date columns to proper date format
soil_dataset$Date <- as.Date(soil_dataset$Date)
tree_dataset$Date <- as.Date(tree_dataset$Date)

# ===== SOIL CH4 FLUX PLOT =====

# Prepare the soil data (same processing as your original tree code)
soil_plot_data <- soil_dataset %>%
  # Filter: keep CH4 flux if CO2 flux OR CH4 flux has R2 > 0.7, and filter CH4 flux outliers
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2 >= 0.7 | CH4_HM.r2 >= 0.7),
         CH4_best.flux >= -100 & CH4_best.flux <= 200) %>%
  # Group WD and WS as W (Wetland)
  mutate(
    Plot_Type = case_when(
      Plot.letter %in% c("WD", "WS") ~ "W",
      Plot.letter == "I" ~ "I", 
      Plot.letter == "U" ~ "U",
      TRUE ~ Plot.letter
    )
  ) %>%
  # Create measurement intervals by grouping dates within 7 days of each other
  arrange(Date) %>%
  mutate(
    Date_numeric = as.numeric(Date),
    Date_group = cumsum(c(TRUE, diff(Date_numeric) > 7)) # New group when gap > 7 days
  ) %>%
  group_by(Date_group) %>%
  mutate(Date_interval = min(Date)) %>% # Use earliest date in each group as representative
  ungroup() %>%
  # Calculate mean flux per plot type per measurement interval
  group_by(Plot_Type, Date_interval) %>%
  summarise(
    mean_flux = mean(CH4_best.flux, na.rm = TRUE),
    se_flux = sd(CH4_best.flux, na.rm = TRUE) / sqrt(n()),
    n_measurements = n(),
    .groups = 'drop'
  ) %>%
  # Remove any rows with missing data
  filter(!is.na(mean_flux), !is.na(Date_interval))

# Prepare raw soil data for jittered points
soil_raw_data <- soil_dataset %>%
  # Filter: keep CH4 flux if CO2 flux OR CH4 flux has R2 > 0.7, and filter CH4 flux outliers
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2 >= 0.7 | CH4_HM.r2 >= 0.7),
         CH4_best.flux >= -100 & CH4_best.flux <= 200) %>%
  # Group WD and WS as W (Wetland)
  mutate(
    Plot_Type = case_when(
      Plot.letter %in% c("WD", "WS") ~ "W",
      Plot.letter == "I" ~ "I", 
      Plot.letter == "U" ~ "U",
      TRUE ~ Plot.letter
    )
  ) %>%
  # Create measurement intervals by grouping dates within 7 days of each other
  arrange(Date) %>%
  mutate(
    Date_numeric = as.numeric(Date),
    Date_group = cumsum(c(TRUE, diff(Date_numeric) > 7)) # New group when gap > 7 days
  ) %>%
  group_by(Date_group) %>%
  mutate(Date_interval = min(Date)) %>% # Use earliest date in each group as representative
  ungroup() %>%
  filter(!is.na(CH4_best.flux), !is.na(Date_interval))

# Create the soil CH4 flux plot
soil_ch4_plot <- ggplot() +
  # Add jittered raw points
  geom_jitter(data = soil_raw_data, 
              aes(x = Date_interval, y = CH4_best.flux, color = Plot_Type),
              alpha = 0.4, size = 1, width = 1) +
  # Add mean points with error bars
  geom_errorbar(data = soil_plot_data, 
                aes(x = Date_interval, y = mean_flux, color = Plot_Type,
                    ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 2, alpha = 0.8) +
  geom_point(data = soil_plot_data, 
             aes(x = Date_interval, y = mean_flux, color = Plot_Type),
             size = 3, alpha = 0.9) +
  geom_smooth(data = soil_plot_data, 
              aes(x = Date_interval, y = mean_flux, color = Plot_Type),
              method = "loess", se = FALSE, size = 1.2) +
  scale_color_manual(
    values = c("I" = "grey", "U" = "forestgreen", "W" = "blue"),
    labels = c("I" = "Intermediate", "U" = "Upland", "W" = "Wetland")
  ) +
  labs(
    x = "Date",
    y = expression(paste("Mean CH"[4], " Flux (Soil)")),
    color = "Plot Type",
    title = "Soil CH4 Flux Over Time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.02))

# Display the soil plot
print(soil_ch4_plot)

# Print summary statistics for soil
cat("Summary of soil data by plot type:\n")
soil_plot_data %>%
  group_by(Plot_Type) %>%
  summarise(
    n_intervals = n(),
    mean_flux_overall = mean(mean_flux, na.rm = TRUE),
    min_flux = min(mean_flux, na.rm = TRUE),
    max_flux = max(mean_flux, na.rm = TRUE)
  ) %>%
  print()

# ===== COMBINED TREES + SOIL PLOT =====

# Prepare tree data (using CH4_best.flux.x)
tree_plot_data <- tree_dataset %>%
  # Filter: keep CH4 flux if CO2 flux OR CH4 flux has R2 > 0.7, and filter CH4 flux outliers
  filter((CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7 | CH4_LM.r2.x >= 0.7 | CH4_HM.r2.x >= 0.7),
         CH4_best.flux.x >= -100 & CH4_best.flux.x <= 200) %>%
  # Group WD and WS as W (Wetland)
  mutate(
    Plot_Type = case_when(
      Plot.Letter %in% c("WD", "WS") ~ "W",
      Plot.Letter == "I" ~ "I", 
      Plot.Letter == "U" ~ "U",
      TRUE ~ Plot.Letter
    ),
    Data_Type = "Trees"
  ) %>%
  # Create measurement intervals by grouping dates within 7 days of each other
  arrange(Date) %>%
  mutate(
    Date_numeric = as.numeric(Date),
    Date_group = cumsum(c(TRUE, diff(Date_numeric) > 7)) # New group when gap > 7 days
  ) %>%
  group_by(Date_group) %>%
  mutate(Date_interval = min(Date)) %>% # Use earliest date in each group as representative
  ungroup() %>%
  # Calculate mean flux per plot type per measurement interval
  group_by(Plot_Type, Date_interval, Data_Type) %>%
  summarise(
    mean_flux = mean(CH4_best.flux.x, na.rm = TRUE),
    se_flux = sd(CH4_best.flux.x, na.rm = TRUE) / sqrt(n()),
    n_measurements = n(),
    .groups = 'drop'
  ) %>%
  # Remove any rows with missing data
  filter(!is.na(mean_flux), !is.na(Date_interval))

# Add Data_Type to soil data and combine
soil_plot_data_labeled <- soil_plot_data %>%
  mutate(Data_Type = "Soil")

# Combine both datasets
combined_plot_data <- bind_rows(tree_plot_data, soil_plot_data_labeled)

# Prepare raw data for combined plot
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
    Data_Type = "Trees"
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

soil_raw_data_labeled <- soil_raw_data %>%
  mutate(Data_Type = "Soil",
         Plot.Tag = as.character(Plot.Tag)) %>%  # Convert to character to match tree data
  rename(CH4_flux = CH4_best.flux) %>%
  filter(CH4_flux >= -100 & CH4_flux <= 200)

# Also convert Plot.Tag in tree data to character for consistency
tree_raw_data <- tree_raw_data %>%
  mutate(Plot.Tag = as.character(Plot.Tag))

combined_raw_data <- bind_rows(tree_raw_data, soil_raw_data_labeled)

# Create the combined plot
combined_ch4_plot <- ggplot() +
  # Add jittered raw points
  geom_jitter(data = combined_raw_data, 
              aes(x = Date_interval, y = CH4_flux, color = Plot_Type, shape = Data_Type),
              alpha = 0.4, size = 1, width = 1) +
  # Add mean points with error bars
  geom_errorbar(data = combined_plot_data, 
                aes(x = Date_interval, y = mean_flux, color = Plot_Type,
                    ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 2, alpha = 0.8) +
  geom_point(data = combined_plot_data, 
             aes(x = Date_interval, y = mean_flux, color = Plot_Type, shape = Data_Type),
             size = 3, alpha = 0.9) +
  geom_smooth(data = combined_plot_data, 
              aes(x = Date_interval, y = mean_flux, color = Plot_Type, linetype = Data_Type),
              method = "loess", se = FALSE, size = 1.2) +
  scale_color_manual(
    values = c("I" = "grey", "U" = "forestgreen", "W" = "blue"),
    labels = c("I" = "Intermediate", "U" = "Upland", "W" = "Wetland")
  ) +
  scale_shape_manual(
    values = c("Trees" = 16, "Soil" = 17),
    labels = c("Trees" = "Trees", "Soil" = "Soil")
  ) +
  scale_linetype_manual(
    values = c("Trees" = "solid", "Soil" = "dashed"),
    labels = c("Trees" = "Trees", "Soil" = "Soil")
  ) +
  labs(
    x = "Date",
    y = expression(paste("Mean CH"[4], " Flux")),
    color = "Plot Type",
    shape = "Data Type",
    linetype = "Data Type",
    title = "Combined Trees and Soil CH4 Flux Over Time"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.02))

# Display the combined plot
print(combined_ch4_plot)

# Print summary statistics for combined data
cat("\nSummary of combined data by plot type and data type:\n")
combined_plot_data %>%
  group_by(Plot_Type, Data_Type) %>%
  summarise(
    n_intervals = n(),
    mean_flux_overall = mean(mean_flux, na.rm = TRUE),
    min_flux = min(mean_flux, na.rm = TRUE),
    max_flux = max(mean_flux, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  print()

# ===== FACETED PLOT BY PLOT TYPE (like your reference image) =====

# Create a comprehensive combined dataset for faceting
combined_facet_data <- bind_rows(
  # Soil data
  soil_plot_data %>% mutate(Data_Type = "Soil"),
  # Tree data  
  tree_plot_data %>% mutate(Data_Type = "Tree")
) %>%
  # Create ordered factor for plot types to control facet order
  mutate(Plot_Type = factor(Plot_Type, levels = c("U", "I", "W"), 
                            labels = c("Upland", "Intermediate", "Wetland")))

# Prepare combined raw data with jittering
combined_facet_raw <- combined_raw_data %>%
  mutate(Plot_Type = factor(Plot_Type, levels = c("U", "I", "W"), 
                            labels = c("Upland", "Intermediate", "Wetland")))

# Create faceted plot
faceted_ch4_plot <- ggplot() +
  # Add jittered raw points first (underneath means)
  geom_jitter(data = combined_facet_raw,
              aes(x = Date_interval, y = CH4_flux, color = Data_Type),
              alpha = 0.3, size = 0.8, width = 8, height = 0) +
  # Add error bars
  geom_errorbar(data = combined_facet_data,
                aes(x = Date_interval, y = mean_flux, color = Data_Type,
                    ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 5, alpha = 0.7, position = position_dodge(width = 10)) +
  # Add mean points
  geom_point(data = combined_facet_data,
             aes(x = Date_interval, y = mean_flux, color = Data_Type),
             size = 2, alpha = 0.8, position = position_dodge(width = 10)) +
  # Facet by plot type with correct order
  facet_wrap(~ Plot_Type, scales = "free_y", nrow = 1) +
  # Color scheme
  scale_color_manual(
    values = c("Soil" = "#8B4513", "Tree" = "#228B22"),  # Earthy brown for soil, forest green for trees
    labels = c("Soil" = "Soil", "Tree" = "Tree")
  ) +
  # Labels and theme
  labs(
    x = "Date",
    y = expression(paste("CH"[4], " Flux (μg C/m²/hr)")),
    color = "Type",
    title = "CH4 Flux by Plot Type: Soil vs Trees"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  # Use pseudo-log transformation to expand values near zero
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.01, base = 10))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7)

# Display the faceted plot
print(faceted_ch4_plot)

# Optional: Save plots
# ggsave("soil_ch4_flux_plot.png", soil_ch4_plot, width = 12, height = 8, dpi = 300)
# ggsave("combined_ch4_flux_plot.png", combined_ch4_plot, width = 12, height = 8, dpi = 300)
# ggsave("faceted_ch4_flux_plot.png", faceted_ch4_plot, width = 14, height = 8, dpi = 300)



# Clean the data to have consistent "Tree" values
combined_facet_data_clean <- combined_facet_data %>%
  mutate(Data_Type = ifelse(Data_Type %in% c("Tree", "Trees"), "Tree", Data_Type))
combined_facet_raw_clean <- combined_facet_raw %>%
  mutate(Data_Type = ifelse(Data_Type %in% c("Tree", "Trees"), "Tree", Data_Type))

# Option 1: Remove scales = "free_y" to use same y-axis
faceted_ch4_plot_same_axis <- ggplot() +
  # Add jittered raw points first (underneath means)
  geom_jitter(data = combined_facet_raw_clean,
              aes(x = Date_interval, y = CH4_flux, color = Data_Type),
              alpha = 0.3, size = 1.5, width = 8, height = 0) +
  # Add error bars
  geom_errorbar(data = combined_facet_data_clean,
                aes(x = Date_interval, y = mean_flux, color = Data_Type,
                    ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 15, alpha = 1, position = position_dodge(width = 10)) +
  # Add mean points
  geom_point(data = combined_facet_data_clean,
             aes(x = Date_interval, y = mean_flux, color = Data_Type),
             size = 3, alpha = 0.8, position = position_dodge(width = 10)) +
  # Facet by plot type with SAME y-axis (removed scales = "free_y")
  facet_wrap(~ Plot_Type, nrow = 1) +
  # Simplified color scheme - only the categories that exist
  scale_color_manual(
    values = c("Soil" = "#8B4513", "Tree" = "#228B22")
  ) +
  # Labels and theme
  labs(
    x = "Date",
    y = expression(paste("CH"[4], " Flux (nmol m"^-2, " s"^-1, ")")),
    color = "Type") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  # Use pseudo-log transformation with custom breaks
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.005),
    breaks = c(-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.9)

print(faceted_ch4_plot_same_axis)


# Clean the data to have consistent "Tree" values
combined_facet_data_clean <- combined_facet_data %>%
  mutate(Data_Type = ifelse(Data_Type %in% c("Tree", "Trees"), "Tree", Data_Type))
combined_facet_raw_clean <- combined_facet_raw %>%
  mutate(Data_Type = ifelse(Data_Type %in% c("Tree", "Trees"), "Tree", Data_Type))

# Linear y-axis with free scales between facets
faceted_ch4_plot_free_linear <- ggplot() +
  # Add jittered raw points first (underneath means)
  geom_jitter(data = combined_facet_raw_clean,
              aes(x = Date_interval, y = CH4_flux, color = Data_Type),
              alpha = 0.3, size = 1.5, width = 8, height = 0) +
  # Add error bars
  geom_errorbar(data = combined_facet_data_clean,
                aes(x = Date_interval, y = mean_flux, color = Data_Type,
                    ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 15, alpha = 1, position = position_dodge(width = 10)) +
  # Add mean points
  geom_point(data = combined_facet_data_clean,
             aes(x = Date_interval, y = mean_flux, color = Data_Type),
             size = 3, alpha = 0.8, position = position_dodge(width = 10)) +
  # Facet by plot type with FREE y-axis scales
  facet_wrap(~ Plot_Type, nrow = 1, scales = "free_y") +
  # Simplified color scheme - only the categories that exist
  scale_color_manual(
    values = c("Soil" = "#8B4513", "Tree" = "#228B22")
  ) +
  # Labels and theme
  labs(
    x = "Date",
    y = expression(paste("CH"[4], " Flux (nmol m"^-2, " s"^-1, ")")),
    color = "Type") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  # Linear y-axis (no transformation)
  scale_y_continuous() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.9)

print(faceted_ch4_plot_free_linear)

# Clean the data to have consistent "Tree" values
combined_facet_data_clean <- combined_facet_data %>%
  mutate(Data_Type = ifelse(Data_Type %in% c("Tree", "Trees"), "Tree", Data_Type))
combined_facet_raw_clean <- combined_facet_raw %>%
  mutate(Data_Type = ifelse(Data_Type %in% c("Tree", "Trees"), "Tree", Data_Type))

# Option 1: Remove scales = "free_y" to use same y-axis
faceted_ch4_plot_same_axis <- ggplot() +
  # Add jittered raw points first (underneath means)
  geom_jitter(data = combined_facet_raw_clean,
              aes(x = Date_interval, y = CH4_flux, color = Data_Type),
              alpha = 0.3, size = 1.5, width = 8, height = 0) +
  # Add error bars
  geom_errorbar(data = combined_facet_data_clean,
                aes(x = Date_interval, y = mean_flux, color = Data_Type,
                    ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 15, alpha = 1, position = position_dodge(width = 10)) +
  # Add mean points
  geom_point(data = combined_facet_data_clean,
             aes(x = Date_interval, y = mean_flux, color = Data_Type),
             size = 3, alpha = 0.8, position = position_dodge(width = 10)) +
  # Facet by plot type with SAME y-axis (removed scales = "free_y")
  facet_wrap(~ Plot_Type, nrow = 1) +
  # Simplified color scheme - only the categories that exist
  scale_color_manual(
    values = c("Soil" = "#8B4513", "Tree" = "#228B22")
  ) +
  # Labels and theme
  labs(
    x = "Date",
    y = expression(paste("CH"[4], " Flux (nmol m"^-2, " s"^-1, ")")),
    color = "Type") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  # Add lines connecting the mean points
  geom_line(data = combined_facet_data_clean,
            aes(x = Date_interval, y = mean_flux, color = Data_Type, group = Data_Type),
            size = .8, alpha = 0.3) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  # Use pseudo-log transformation with custom breaks
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.005),
    breaks = c(-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.9)

print(faceted_ch4_plot_same_axis)
