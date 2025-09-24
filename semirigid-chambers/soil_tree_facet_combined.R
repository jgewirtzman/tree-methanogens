# Load required libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(cowplot)
library(scales)

# Read the datasets (using your existing data loading code)
soil_dataset <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/semirigid_tree_final_complete_dataset_soil.csv')
tree_dataset <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/flux_code/untitled folder/semirigid_tree_final_complete_dataset.csv')
moisture_data <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/ipad_data/Cleaned data/soilmoisture_total.csv')

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
  #geom_point(aes(y = 0), position = position_jitter(height = 0.001), alpha = 0.3, size = 1) +
  facet_wrap(~ Plot_Type, strip.position = "bottom") +
  scale_fill_brewer(type = "seq", palette = "Blues") +
  scale_x_continuous(position = "top") +
  labs(
    x = "Volumetric Water Content (%)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    panel.grid = element_blank(),
    panel.grid.y = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank()
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
            size = .8, alpha = 0.3) +
  facet_wrap(~ Plot_Type, nrow = 1) +
  scale_color_manual(values = c("Soil" = "#8B4513", "Tree" = "#228B22")) +
  scale_fill_manual(values = c("Soil" = "#8B4513", "Tree" = "#228B22")) +
  labs(
    y = bquote("CH"[4] ~ "Flux (nmol m"^-2 ~ "s"^-1 ~ ")")  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    panel.border = element_rect(color = "grey30", fill = NA, size = 1),  # Changed from grey80 and size 0.5
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.005),
    breaks = c(-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.9)


# ===== CREATE SOIL TEMP/MOISTURE PLOT (BOTTOM) =====
temp_moisture_plot <- ggplot(soil_temp_moisture_data, aes(x = Date)) +
  #geom_jitter(aes(y = Soil_Temperature, color = "Temperature"), alpha = 0.5, width=3, size=1) +
  geom_smooth(aes(y = Soil_Temperature, color = "Temperature"), method = "loess", se = FALSE, alpha = 0.7, size=0.75) +
  #geom_jitter(aes(y = Soil_Moisture / 2.5, color = "Moisture"), alpha = 0.5, width=3, size=1) +
  geom_smooth(aes(y = Soil_Moisture / 2.5, color = "Moisture"), method = "loess", se = FALSE, alpha = 0.7, size=0.75) +
  facet_wrap(~ Plot_Type) +
  scale_y_continuous(
    name = "Soil T (°C)",
    sec.axis = sec_axis(~ . * 2.5, name = "VWC (%)")
  ) +
  scale_color_manual(values = c("Temperature" = alpha("#D73027", 0.5),
                                "Moisture" = alpha("#4575B4", 0.5))) +
  labs(
    x = "Date"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +  xlab(NULL) +
  guides(color = guide_legend(title = NULL)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months")

# ===== COMBINE ALL PLOTS (no legend gathering) =====
final_plot <- plot_grid(
  vwc_plot,
  ch4_flux_plot, 
  temp_moisture_plot,
  ncol = 1,
  align = "v",  # Try "hv" instead of "v"
  axis = "lr", # Try "tblr" instead of "lr"
  rel_heights = c(.75, 1.5, 1)
)

# Display the final combined plot
print(final_plot)

# Optional: Save the plot
# ggsave("combined_three_panel_plot.png", final_plot, width = 16, height = 12, dpi = 300)



# Load required libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(cowplot)
library(scales)

# [Keep all your existing data preparation code here - I'm showing just the plot modifications]

# ===== CREATE VWC DENSITY PLOT (TOP) - MODIFIED =====
vwc_plot <- ggplot(vwc_data, aes(x = Soil_Moisture, fill = Plot_Type)) +
  geom_density(alpha = 0.2) +
  #geom_point(aes(y = 0), position = position_jitter(height = 0.001), alpha = 0.3, size = 1) +
  facet_wrap(~ Plot_Type, strip.position = "bottom") +
  scale_fill_brewer(type = "seq", palette = "Blues") +
  scale_x_continuous(position = "top") +
  labs(
    x = "Volumetric Water Content (%)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    # NEW: Reduce margins and spacing
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"),  # Reduce bottom margin
    strip.placement = "outside",
    panel.spacing = unit(0.1, "lines")  # Reduce spacing between facet panels
  ) +
  guides(fill = guide_legend(title = NULL))

# ===== CREATE CH4 FLUX PLOT (MIDDLE) - MODIFIED =====
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
            size = .8, alpha = 0.3) +
  facet_wrap(~ Plot_Type, nrow = 1) +
  scale_color_manual(values = c("Soil" = "#8B4513", "Tree" = "#228B22")) +
  scale_fill_manual(values = c("Soil" = "#8B4513", "Tree" = "#228B22")) +
  labs(
    y = bquote("CH"[4] ~ "Flux (nmol m"^-2 ~ "s"^-1 ~ ")")  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    panel.border = element_rect(color = "grey30", fill = NA, size = 1),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    # NEW: Reduce margins and spacing
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"),  # Remove top margin
    panel.spacing = unit(0.1, "lines")  # Reduce spacing between facet panels
  ) +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.005),
    breaks = c(-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20", alpha = 0.9)

# ===== CREATE SOIL TEMP/MOISTURE PLOT (BOTTOM) - MODIFIED =====
temp_moisture_plot <- ggplot(soil_temp_moisture_data, aes(x = Date)) +
  #geom_jitter(aes(y = Soil_Temperature, color = "Temperature"), alpha = 0.5, width=3, size=1) +
  geom_smooth(aes(y = Soil_Temperature, color = "Temperature"), method = "loess", se = FALSE, alpha = 0.7, size=0.75) +
  #geom_jitter(aes(y = Soil_Moisture / 2.5, color = "Moisture"), alpha = 0.5, width=3, size=1) +
  geom_smooth(aes(y = Soil_Moisture / 2.5, color = "Moisture"), method = "loess", se = FALSE, alpha = 0.7, size=0.75) +
  facet_wrap(~ Plot_Type) +
  scale_y_continuous(
    name = "Soil T (°C)",
    sec.axis = sec_axis(~ . * 2.5, name = "VWC (%)")
  ) +
  scale_color_manual(values = c("Temperature" = alpha("#D73027", 0.5),
                                "Moisture" = alpha("#4575B4", 0.5))) +
  labs(
    x = "Date"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    # NEW: Reduce top margin and legend spacing
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"),
    panel.spacing = unit(0.1, "lines"),
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt"),  # Reduce legend top margin
    legend.box.margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt")  # Reduce legend box margin
  ) +  
  xlab(NULL) +
  guides(color = guide_legend(title = NULL)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months")

# ===== COMBINE ALL PLOTS WITH TIGHTER SPACING =====
final_plot <- plot_grid(
  vwc_plot,
  ch4_flux_plot, 
  temp_moisture_plot,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(.75, 1.5, 1)
)

# Display the final combined plot
print(final_plot)

# Optional: Save the plot
# ggsave("combined_three_panel_plot.png", final_plot, width = 16, height = 12, dpi = 300)