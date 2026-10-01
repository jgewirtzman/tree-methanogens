# ==============================================================================
# mcrA Linear Models
# ==============================================================================
# Purpose: Linear models relating mcrA gene abundance to CH4 flux.
#
# Pipeline stage: 3 — Analysis
#
# Inputs:
#   - tree_methane_all_data.csv (from data/raw/external/)
# ==============================================================================

library(tidyverse)
library(sjPlot)

all_data<-read.csv("../../../data/raw/external/tree_methane_all_data.csv")

model<-lm(CH4_int~mcra_probe_strict, data=all_data)
summary(model)

tab_model(model, standardize = "refit")


tab_model(model_mixed, standardize = "refit")


sd(scale(all_data$mcra_probe_strict), na.rm=T)


ggplot(aes(x=log(CH4_int), y=log(ch4_125)), data=all_data)+
  geom_point()+
  facet_wrap(~species)+
  geom_smooth(method="lm")

ggplot(aes(x=log(CH4_int), y=log(ch4_125)), data=all_data)+
  geom_point()+
  geom_smooth(method="lm")

##

ggplot(aes(x=log(mcra_probe_strict), y=log(ch4_125)), data=all_data)+
  geom_point()+
  facet_wrap(~species)+
  geom_smooth(method="lm")

ggplot(aes(x=log(mcra_probe_strict), y=log(ch4_125)), data=all_data)+
  geom_point()+
  geom_smooth(method="lm")

##
ggplot(aes(x=log(mcra_probe_strict), y=log(CH4_int)), data=all_data)+
  geom_point()+
  facet_wrap(~species*material)+
  geom_smooth(method="lm")

ggplot(aes(x=log(mcra_probe_strict), y=log(CH4_int)), data=all_data)+
  geom_point()+
  geom_smooth(method="lm")
