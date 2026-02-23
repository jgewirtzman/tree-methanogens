# ==============================================================================
# Simple Tree Flux Diagnostic Plots
# ==============================================================================
# Purpose: Simple diagnostic plots of tree flux patterns.
# ==============================================================================

library(ggplot2)
ggplot(final_dataset, aes(x=Date, y=CH4_best.flux.x, color=`Plot Letter`))+
  geom_point()+geom_smooth()+ylim(0,1)

unique(final_dataset$Date)
head(final_dataset)
str(final_dataset)

















library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)

# Prepare the data
plot_data <- final_dataset %>%
  # Filter out points where CO2 flux r2 < 0.7
  filter(CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7) %>%
  # Group WD and WS as W (Wetland)
  mutate(
    Plot_Type = case_when(
      `Plot Letter` %in% c("WD", "WS") ~ "W",
      `Plot Letter` == "I" ~ "I", 
      `Plot Letter` == "U" ~ "U",
      TRUE ~ `Plot Letter`
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
    mean_flux = mean(CH4_best.flux.x, na.rm = TRUE),
    se_flux = sd(CH4_best.flux.x, na.rm = TRUE) / sqrt(n()),
    n_measurements = n(),
    .groups = 'drop'
  ) %>%
  # Remove any rows with missing data
  filter(!is.na(mean_flux), !is.na(Date_interval))

# Prepare raw data for jittered points
raw_data <- final_dataset %>%
  # Filter out points where CO2 flux r2 < 0.7
  filter(CO2_LM.r2 >= 0.7 | CO2_HM.r2 >= 0.7) %>%
  # Group WD and WS as W (Wetland)
  mutate(
    Plot_Type = case_when(
      `Plot Letter` %in% c("WD", "WS") ~ "W",
      `Plot Letter` == "I" ~ "I", 
      `Plot Letter` == "U" ~ "U",
      TRUE ~ `Plot Letter`
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
  filter(!is.na(CH4_best.flux.x), !is.na(Date_interval))

# Create the plot
ggplot() +
  # Add jittered raw points
  geom_jitter(data = raw_data, 
              aes(x = Date_interval, y = CH4_best.flux.x, color = Plot_Type),
              alpha = 0.4, size = 1, width = 1) +
  # Add mean points with error bars
  geom_errorbar(data = plot_data, 
                aes(x = Date_interval, y = mean_flux, color = Plot_Type,
                    ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 2, alpha = 0.8) +
  geom_point(data = plot_data, 
             aes(x = Date_interval, y = mean_flux, color = Plot_Type),
             size = 3, alpha = 0.9) +
  geom_smooth(data = plot_data, 
              aes(x = Date_interval, y = mean_flux, color = Plot_Type),
              method = "loess", se = FALSE, size = 1.2) +
  scale_color_manual(
    values = c("I" = "grey", "U" = "forestgreen", "W" = "blue"),
    labels = c("I" = "Intermediate", "U" = "Upland", "W" = "Wetland")
  ) +
  labs(
    x = "Date",
    y = expression(paste("Mean CH"[4], " Flux")),
    color = "Plot Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11)
  ) +
  scale_x_datetime(date_labels = "%b %Y", date_breaks = "2 months") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.02)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Optional: Print summary statistics
cat("Summary of data by plot type:\n")
plot_data %>%
  group_by(Plot_Type) %>%
  summarise(
    n_intervals = n(),
    mean_flux_overall = mean(mean_flux, na.rm = TRUE),
    min_flux = min(mean_flux, na.rm = TRUE),
    max_flux = max(mean_flux, na.rm = TRUE)
  ) %>%
  print()
