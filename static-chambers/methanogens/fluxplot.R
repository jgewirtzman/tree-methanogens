ph_gas<-read.csv("/Users/jongewirtzman/Downloads/hf005-01-trace-gas - hf005-01-trace-gas.csv")

library(ggplot2)
ggplot(ph_gas, aes(x=date, y=co2_flux, color=treatment))+
  geom_line()

data<-ph_gas


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(readr)
library(tidyr)

# Convert date to proper Date format
data$date <- as.Date(data$date)

# Function to calculate mean and standard error by date and treatment
calculate_stats <- function(data, value_column) {
  data %>%
    filter(!is.na(!!sym(value_column))) %>%
    group_by(date, treatment) %>%
    summarize(
      mean = mean(!!sym(value_column), na.rm = TRUE),
      se = sd(!!sym(value_column), na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
}

# Calculate statistics for each gas
co2_stats <- calculate_stats(data, "co2_flux")
ch4_stats <- calculate_stats(data, "ch4_flux")
n2o_stats <- calculate_stats(data, "n2o_flux")

# Create a function to generate plots
create_flux_plot <- function(stats_data, y_label, title, color_palette = c("#1B9E77", "#D95F02", "#7570B3")) {
  # Convert treatment to factor for better plotting
  stats_data$treatment <- as.factor(stats_data$treatment)
  
  ggplot(stats_data, aes(x = date, y = mean, color = treatment, group = treatment)) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 20) +
    scale_color_manual(values = color_palette, 
                       name = "Treatment",
                       labels = paste("Treatment", 1:3)) +
    labs(title = title,
         x = "Date",
         y = y_label) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# Create the plots
co2_plot <- create_flux_plot(
  co2_stats, 
  y_label = expression(CO[2]~Flux~(mg~C~m^{-2}~hr^{-1})), 
  title = expression(CO[2]~Flux~by~Date~and~Treatment)
)

ch4_plot <- create_flux_plot(
  ch4_stats, 
  y_label = expression(CH[4]~Flux~(mg~C~m^{-2}~hr^{-1})), 
  title = expression(CH[4]~Flux~by~Date~and~Treatment)
)

n2o_plot <- create_flux_plot(
  n2o_stats, 
  y_label = expression(N2O~Flux~(mg~N~m^{-2}~hr^{-1})), 
  title = expression(N2O~Flux~by~Date~and~Treatment)
)

# Print plots
print(co2_plot)
print(ch4_plot)
print(n2o_plot)

# If you want to save the plots:
ggsave("co2_flux_plot.png", co2_plot, width = 10, height = 6)
ggsave("ch4_flux_plot.png", ch4_plot, width = 10, height = 6)
ggsave("n2o_flux_plot.png", n2o_plot, width = 10, height = 6)

# To create a combined plot with multiple panels
library(patchwork)

# Adjust individual plots for combined view
co2_plot_combined <- co2_plot + theme(legend.position = "none") + ggtitle("CO2 Flux")
ch4_plot_combined <- ch4_plot + theme(legend.position = "none") + ggtitle("CH4 Flux")
n2o_plot_combined <- n2o_plot + theme(legend.position = "none") + ggtitle("N2O Flux")

# Common legend
legend <- cowplot::get_legend(co2_plot + theme(legend.position = "bottom"))

# Combine plots
combined_plot <- (co2_plot_combined / ch4_plot_combined / n2o_plot_combined) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save the combined plot
ggsave("combined_gas_flux_plot.png", combined_plot, width = 10, height = 12)
