source("code/lib/outputs.R")
# ==============================================================================
# REVISION — Fig 7c with the CH4 flux UNIT added (Referee 2 #6)
# ==============================================================================
# R2: Fig 7c (felled-oak stem CH4 flux vs height) is "missing the CH4 flux unit."
#
# DIAGNOSIS (code/08_figures/09_felled_oak_profiles.R:221):
#   Panel 3 draws  ylab(expression(CH[4]~"Flux"))  -- no unit. The flux is the
#   goFlux CH4_best.flux (same variable/pipeline as Fig 2; range here -0.005 to
#   0.40, median 0.039), i.e. nmol m-2 s-1.
#
# FIX (the whole fix): y-axis label -> CH4 flux (nmol m-2 s-1).
#   >>> ONE-LINE PATCH to apply in 09_felled_oak_profiles.R:221 (for the merge):
#       ylab(expression(CH[4]~"flux (nmol m"^-2~"s"^-1*")"))
#
# This NEW file regenerates panel c standalone with the corrected label (and the
# same internal-CH4 point coloring), so the fix can be seen before merging. It
# reproduces only the CH4 calibration path from the source script; it does not
# edit the pipeline or overwrite any tracked output.
# Run: Rscript code/revision/R_fig7c_flux_unit.R
# Output: outputs/revision/fig7c_flux_unit_fixed.png
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(viridis) })
out <- "outputs/revision"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- CH4 calibration (mirrors 09_felled_oak_profiles.R, CH4 path only) --------
GHG_standards <- read.csv("data/raw/field_data/black_oak/ghg_standards.csv", fileEncoding = "UTF-8-BOM")
GC_data       <- read.csv("data/raw/field_data/black_oak/quve_gc_data.csv",  fileEncoding = "UTF-8-BOM")
colnames(GC_data)       <- make.names(colnames(GC_data))
colnames(GHG_standards) <- make.names(colnames(GHG_standards))

GHG_data <- GHG_standards %>% left_join(GC_data, by = c("Sample" = "Sample.Type"))
GHG_data$CH4.Area <- as.numeric(GHG_data$CH4.Area)
GHG_data$CH4_ppm  <- as.numeric(gsub(",", "", GHG_data$X.CH4...ppm.))
GHG_data <- GHG_data %>% filter(!is.na(CH4.Area), !is.na(CH4_ppm))

low_range_standards <- c("N2", "SB1", "SB2", "SB3")
low_range_data <- GHG_data %>% filter(Sample %in% low_range_standards)
CH4_low_curve  <- lm(CH4_ppm ~ CH4.Area, data = low_range_data)
CH4_high_curve <- lm(CH4_ppm ~ CH4.Area, data = GHG_data)
max_CH4_low_area <- 1000

GC_data <- GC_data %>% mutate(
  CH4_concentration = if_else(CH4.Area <= max_CH4_low_area,
    predict(CH4_low_curve,  newdata = data.frame(CH4.Area = CH4.Area)),
    predict(CH4_high_curve, newdata = data.frame(CH4.Area = CH4.Area))))

int_gas <- GC_data %>% filter(!is.na(Lab.ID), Tree.Tissue == "Trunk Gas")
int_gas$Tree.Height <- as.numeric(int_gas$Tree.Height)

# ---- flux data + internal-CH4 color (mirrors source) -------------------------
flux_data <- read.csv("data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv")
flux_data$Height_m <- suppressWarnings(as.numeric(flux_data$Height_m))
flux_df <- flux_data %>% filter(!is.na(Height_m)) %>%
  select(Height_m, CH4_best.flux) %>% rename(height = Height_m, flux = CH4_best.flux) %>%
  left_join(int_gas %>% group_by(Tree.Height) %>%
              summarize(ch4 = mean(CH4_concentration, na.rm = TRUE), .groups = "drop"),
            by = c("height" = "Tree.Height"))

height_breaks <- sort(unique(int_gas$Tree.Height))
shared_theme <- theme_classic(base_size = 12) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 11),
        legend.position = "bottom", legend.title = element_text(size = 10),
        legend.text = element_text(size = 9), legend.box.margin = margin(-5, 0, 0, 0))

# ---- Panel c, corrected label ------------------------------------------------
p3 <- ggplot(flux_df, aes(y = flux, x = height, fill = ch4)) +
  geom_smooth(se = FALSE, color = "black") +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.6, alpha = 0.85) +
  shared_theme +
  scale_fill_viridis_c(option = "E", direction = 1, name = expression(CH[4]~"(ppm)")) +
  coord_flip() + xlab("Height (m)") +
  ylab(expression(CH[4]~"flux (nmol m"^-2~"s"^-1*")")) +        # <-- UNIT ADDED
  scale_x_continuous(breaks = height_breaks, minor_breaks = NULL) +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top"))

ggsave(out_path("fig7c_flux_unit_fixed.png"), p3, width = 4.2, height = 5.1,
       dpi = 300, bg = "white")
cat("Wrote", out_path("fig7c_flux_unit_fixed.png"), "\n")
cat("Fix: y-axis label 'CH4 Flux' -> 'CH4 flux (nmol m-2 s-1)'.\n")
cat("Apply in code/08_figures/09_felled_oak_profiles.R:221 (with approval).\n")
