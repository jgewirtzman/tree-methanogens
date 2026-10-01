# ==============================================================================
# 04_process_internal_gas_2024.R -- calibrate the October 2024 heartwood /
# sapwood GC run
# ------------------------------------------------------------------------------
# stem_gas_isotopes_picarro_run.csv (GC run 2024-10-10/11, paired heartwood and
# sapwood samples) arrived with concentration columns computed elsewhere, and
# they include negative CH4 and CO2. Same fault as the 2021 survey before
# 03_process_internal_gas.R was recalibrated: a fit that does not respect the
# low end. This mirrors tree-gas-traits code/00b_process_ymf_2024.R exactly:
#   - standards: the certified tanks SB1, SB3, SB4, SB5 only (the a/b/c
#     dilutions in this run have no recorded concentrations)
#   - CH4: weighted linear; CO2: weighted quadratic; weights 1/(conc + c0)^2
#   - O2 and N2O: not calibrated (no O2 standards in this run)
# The input file is left untouched. Calibrated columns go to a new file.
#
# Run from the repo root:  Rscript code/01_import/04_process_internal_gas_2024.R
# Reads:  data/processed/internal_gas/stem_gas_isotopes_picarro_run.csv
#         data/raw/internal_gas/Internal Concentration.xlsx (certified SB tanks)
# Writes: data/processed/internal_gas/stem_gas_2024_calibrated.csv
#         outputs/audit/internal_gas_2024_calibration.txt
# ==============================================================================
source("code/lib/outputs.R")
suppressMessages({ library(readr); library(dplyr); library(readxl) })

IN  <- "data/processed/internal_gas/stem_gas_isotopes_picarro_run.csv"
OUT <- "data/processed/internal_gas/stem_gas_2024_calibrated.csv"
d   <- read_csv(IN, show_col_types = FALSE)

cert <- read_excel("data/raw/internal_gas/Internal Concentration.xlsx",
                   sheet = "Standard Concentrations") %>%
  transmute(Sample, CH4 = suppressWarnings(as.numeric(`[CH4] (ppm)`)),
            CO2 = suppressWarnings(as.numeric(`[CO2] (ppm)`))) %>%
  filter(Sample %in% c("SB1", "SB3", "SB4", "SB5"))

stds <- d %>% filter(Sample.ID %in% cert$Sample, !is.na(CH4_Area)) %>%
  left_join(cert, by = c("Sample.ID" = "Sample"))
stopifnot(setequal(unique(stds$Sample.ID), cert$Sample))

ch4_fit <- lm(CH4 ~ CH4_Area, data = stds, weights = 1 / (CH4 + 0.3)^2)
co2_fit <- lm(CO2 ~ CO2_Area + I(CO2_Area^2), data = stds, weights = 1 / (CO2 + 50)^2)

out <- d %>% mutate(
  CH4_calibrated_ppm = ifelse(is.na(CH4_Area), NA_real_,
                              as.numeric(predict(ch4_fit, newdata = data.frame(CH4_Area = CH4_Area)))),
  CO2_calibrated_ppm = ifelse(is.na(CO2_Area), NA_real_,
                              as.numeric(predict(co2_fit, newdata = data.frame(CO2_Area = CO2_Area)))))
write_csv(out, OUT)

check <- stds %>%
  mutate(CH4_pred = predict(ch4_fit, .), CO2_pred = predict(co2_fit, .)) %>%
  group_by(Sample.ID) %>%
  summarise(CH4_cert = first(CH4), CH4_err = round(median(CH4_pred / CH4 - 1), 3),
            CO2_cert = first(CO2), CO2_err = round(median(CO2_pred / CO2 - 1), 3), .groups = "drop")
s <- out %>% filter(Sample.Type == "Sample")
sink(out_path("internal_gas_2024_calibration.txt"))
cat("2024 heartwood/sapwood GC calibration (back-predicted certified tanks, relative error)\n\n")
print(as.data.frame(check), row.names = FALSE)
cat(sprintf("\nbiological samples: %d\n", nrow(s)))
cat(sprintf("  CH4: lab column min %.2f, %d negative | recalibrated min %.2f, %d negative\n",
            min(s$CH4_concentration, na.rm = TRUE), sum(s$CH4_concentration < 0, na.rm = TRUE),
            min(s$CH4_calibrated_ppm, na.rm = TRUE), sum(s$CH4_calibrated_ppm < 0, na.rm = TRUE)))
cat(sprintf("  CO2: lab column min %.0f, %d negative | recalibrated min %.0f, %d negative\n",
            min(s$CO2_concentration, na.rm = TRUE), sum(s$CO2_concentration < 0, na.rm = TRUE),
            min(s$CO2_calibrated_ppm, na.rm = TRUE), sum(s$CO2_calibrated_ppm < 0, na.rm = TRUE)))
cat(sprintf("  Spearman, lab vs recalibrated CH4: %.3f\n",
            cor(s$CH4_concentration, s$CH4_calibrated_ppm, method = "spearman", use = "complete.obs")))
sink()
cat(readLines(out_path("internal_gas_2024_calibration.txt")), sep = "\n")
