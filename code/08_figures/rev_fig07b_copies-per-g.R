source("code/lib/outputs.R")
# ==============================================================================
# REVISION — Fig 7b fix: black-oak mcrA in copies per gram (was copies/uL)
# ==============================================================================
# Referee 2 (Fig 7b): "Why was mcrA copies per uL used instead of mcrA copies per
# gram? The former is variable depending on both elution volume and sample quantity."
#
# We now have the black-oak extraction masses ("Sample Mass Added to Tube (mg)")
# from Black Oak Tree Project Data.xlsx (extracted to data/processed/molecular/black_oak/),
# which join to black_oak_mcrA.csv by Sample ID. Convert:
#     copies/g = Conc(copies/uL) * elution_uL / (mass_mg/1000)
# Wood was freeze-dried before weighing, so black-oak wood copies/g is DRY-weight.
#
# NEW file; edits nothing existing. Run: Rscript code/revision/R_fig7b_copies_per_gram.R
# Outputs: outputs/revision/fig7b_copies_per_gram.png, fig7b_copies_per_gram.csv
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out_dir <- "outputs/revision"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ELUTION_UL <- 75   # CANONICAL pipeline value (02_harmonize_all_data.R:317: Conc*75/mass*1000).
                   # Earlier draft used 37.5 (=75/2) -> 2x low; corrected to match manuscript.

mcra <- read_csv("data/raw/ddpcr/black_oak_mcrA.csv", show_col_types = FALSE)
conc_col <- grep("Conc", names(mcra), value = TRUE)[1]
mass <- read_csv("data/processed/molecular/black_oak/bo_extraction_mass.csv", show_col_types = FALSE)

d <- mcra %>%
  transmute(sample_id = `Sample ID`, height_cm = `Height (cm)`, component = Component,
            conc_uL = .data[[conc_col]]) %>%
  filter(component %in% c("Heartwood", "Sapwood")) %>%
  left_join(mass %>% transmute(sample_id = `Sample ID`,
                               mass_mg = `Sample Mass Added to Tube (mg)`), by = "sample_id") %>%
  mutate(copies_per_g = conc_uL * ELUTION_UL / (mass_mg/1000))

matched <- sum(!is.na(d$mass_mg)); total <- nrow(d)
write.csv(d, out_path("fig7b_copies_per_gram.csv"), row.names = FALSE)

# ---- Figure: copies/uL vs copies/g profiles (show the fix) --------------------
p_uL <- ggplot(d, aes(conc_uL, height_cm/100, color = component)) +
  geom_path() + geom_point(size = 2) +
  scale_color_manual(values = c(Heartwood = "#8c510a", Sapwood = "#2166ac")) +
  labs(title = "Fig 7b (original units): mcrA copies/uL", x = "mcrA (copies/uL)",
       y = "Height (m)", color = NULL) + theme_bw(base_size = 10)
p_g <- ggplot(d %>% filter(!is.na(copies_per_g)), aes(copies_per_g, height_cm/100, color = component)) +
  geom_path() + geom_point(size = 2) +
  scale_color_manual(values = c(Heartwood = "#8c510a", Sapwood = "#2166ac")) +
  labs(title = "Fig 7b (fixed): mcrA copies/g dry wood", x = "mcrA (copies / g dry)",
       y = "Height (m)", color = NULL) + theme_bw(base_size = 10)
if (requireNamespace("gridExtra", quietly = TRUE))
  ggsave(out_path("fig7b_copies_per_gram.png"),
         gridExtra::arrangeGrob(p_uL, p_g, nrow = 1), width = 11, height = 5, dpi = 150)

cat(sprintf("Matched %d/%d black-oak mcrA samples to extraction mass.\n", matched, total))
cat("Wood copies/g here is DRY-weight (cores freeze-dried before weighing).\n")
cat(sprintf("Elution assumed %.1f uL (CONFIRM for black-oak extraction).\n", ELUTION_UL))
cat("Wrote fig7b_copies_per_gram.png + .csv\n")
