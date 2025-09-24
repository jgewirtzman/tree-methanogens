# Install necessary packages (if not installed already)
install.packages("lavaan")
install.packages("semPlot")

# Load the packages
library(lavaan)
library(semPlot)
library(dplyr)


all_data<-read.csv("/Users/jongewirtzman/Downloads/tree_methane_all_data.csv")

# Filter for only core_type "Inner" and "Outer"
filtered_data <- all_data %>%
  filter(core_type %in% c("Inner", "Outer") & species != "BEAL")


filtered_data$mcrA<-(filtered_data$mcra_probe_loose+1)
filtered_data$CH4<-(filtered_data$CH4_int)
filtered_data$Flux<-(filtered_data$ch4_125)


# Define the path model
path_model <- '
 # mcrA affects CH4, modified by species (interaction)
  CH4 ~ mcrA + pmoa_loose + species

  # CH4 affects methane flux, modified by species (interaction)
  Flux ~ CH4 + species
'

# Fit the path model using your data
fit <- sem(path_model, data = filtered_data)

# Summarize the model fit with fit statistics and standardized coefficients
summary(fit, standardized = TRUE, fit.measures = TRUE)

# Visualize the path model with standardized coefficients
semPaths(fit, "standardized", layout = "tree", whatLabels = "std", edge.label.cex = 1.2)

# Extract important fit statistics
fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea"))

# Calculate and view indirect effects
indirect_effects <- parameterEstimates(fit, standardized = TRUE)
indirect_effects[indirect_effects$op == ":=",]

