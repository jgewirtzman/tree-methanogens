# Load required libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(cowplot)

# Read the data
data <- read.csv('/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens/to-organize/yale-forest-ch4/ipad_data/Cleaned data/soilmoisture_total.csv')

# Simple data prep
soil_data <- data %>%
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

# Single plot with secondary y-axis
ggplot(soil_data, aes(x = Date)) +
  geom_jitter(aes(y = Soil_Temperature, color = "Temperature"), alpha = 0.5, width=3) +
  geom_smooth(aes(y = Soil_Temperature, color = "Temperature"), method = "loess", se = FALSE, alpha = 0.7) +
  geom_jitter(aes(y = Soil_Moisture / 2.5, color = "Moisture"), alpha = 0.5, width=3) +
  geom_smooth(aes(y = Soil_Moisture / 2.5, color = "Moisture"), method = "loess", se = FALSE, alpha = 0.7) +
  facet_wrap(~ Plot_Type) +
  scale_y_continuous(
    name = "Temperature (°C)",
    sec.axis = sec_axis(~ . * 2.5, name = "VWC (%)")
  ) +
  scale_color_manual(values = c("Temperature" = alpha("#D73027", 0.7), "Moisture" = alpha("#4575B4", 0.7))) +
  labs(
    x = "Date",
    color = "Measurement"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
