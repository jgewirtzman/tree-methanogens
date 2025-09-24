# Load required libraries
library(ggplot2)
library(dplyr)
library(lubridate)
library(ggridges)

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
    Soil_Moisture = as.numeric(VWC)
  ) %>%
  filter(!is.na(Soil_Moisture))

# Distribution plot faceted by plot type with raw points
ggplot(soil_data, aes(x = Soil_Moisture, fill = Plot_Type)) +
  geom_density(alpha = 0.7) +
  geom_point(aes(y = 0), position = position_jitter(height = 0.001), 
             alpha = 0.3, size = 1) +
  facet_wrap(~ Plot_Type) +
  scale_fill_brewer(type = "seq", palette = "Blues") +
  labs(
    #title = "Distribution of Soil Moisture (VWC) by Plot Type",
    x = "Volumetric Water Content (%)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
