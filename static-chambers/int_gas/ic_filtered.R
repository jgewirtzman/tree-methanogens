library(tidyverse)

gas<-read.csv("/Users/jongewirtzman/Google Drive/Research/YMF Tree Microbiomes & Methane/Tree Methane Lab/Data/2021 Data/int_conc_prelim.csv")

filtered_data <- gas %>%
  filter(grepl("^[A-Za-z]{4}$", Species.ID))

plot(log(filtered_data$CH4.ppm)~filtered_data$O2.ppm)
abline(lm(log(filtered_data$CH4.ppm)~filtered_data$O2.ppm))

       