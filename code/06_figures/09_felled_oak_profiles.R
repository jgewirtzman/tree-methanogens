# ==============================================================================
# Felled Oak Vertical Profiles (Figure 7)
# ==============================================================================
# Purpose: Creates a 3-panel figure showing vertical profiles of internal CH4
#   concentration, mcrA gene abundance, and CH4 flux along the trunk of a
#   felled black oak (QUVE).
#
# Pipeline stage: 4 — Publication Figures
#
# Inputs:
#   - quve_gc_data.csv (from data/raw/field_data/black_oak/)
#   - ghg_standards.csv (from data/raw/field_data/black_oak/)
#   - o2_standards.csv (from data/raw/field_data/black_oak/)
#   - black_oak_mcrA.csv (from data/raw/ddpcr/)
#   - ymf_black_oak_flux_compiled.csv (from data/raw/field_data/black_oak/)
#
# Outputs:
#   - fig7_felled_oak_profiles.png
#
# Required packages: tidyverse, patchwork, viridis, broom
# ==============================================================================

# Load required libraries
library(tidyverse)
library(patchwork)
library(broom)

# ==============================================================================
# STEP 1: Process GC data from raw (inline from GC_processing_quve.R)
# ==============================================================================

# Load raw data (fileEncoding handles BOM in CSV files)
O2_standards <- read.csv("data/raw/field_data/black_oak/o2_standards.csv", fileEncoding = "UTF-8-BOM")
GHG_standards <- read.csv("data/raw/field_data/black_oak/ghg_standards.csv", fileEncoding = "UTF-8-BOM")
GC_data <- read.csv("data/raw/field_data/black_oak/quve_gc_data.csv", fileEncoding = "UTF-8-BOM")

# Clean column names
colnames(GC_data) <- make.names(colnames(GC_data))
colnames(O2_standards) <- make.names(colnames(O2_standards))
colnames(GHG_standards) <- make.names(colnames(GHG_standards))

# Set N2 O2 area to NA
GC_data$O2.Area[which(GC_data$Sample.Type == "N2")] <- NA

# --- O2 Standard Curves ---
O2_data <- O2_standards %>%
  left_join(GC_data, by = c("Sample" = "Sample.Type"))

O2_data$O2.Area <- as.numeric(O2_data$O2.Area)
O2_data$O2_ppm <- as.numeric(gsub(",", "", O2_data$X.O2...ppm.))

# Define low/high range O2 standards
low_range_o2_standards <- c("Outdoor Air 2", "Outdoor Air 3", "Outdoor Air 4",
                             "Oxygen Standard 1", "Oxygen Standard 2",
                             "Oxygen Standard 3", "Oxygen Standard 4")
high_range_o2_standards <- c("Outdoor Air 5", "Outdoor Air 6", "Outdoor Air 7",
                              "Oxygen Standard 5")

low_range_o2_data <- O2_data %>% filter(Sample %in% low_range_o2_standards)
all_high_range_o2_data <- bind_rows(
  low_range_o2_data,
  O2_data %>% filter(Sample %in% high_range_o2_standards)
)

O2_low_curve <- lm(O2_ppm ~ O2.Area, data = low_range_o2_data)
O2_high_curve <- lm(O2_ppm ~ O2.Area, data = all_high_range_o2_data)
max_O2_low_area <- max(low_range_o2_data$O2.Area, na.rm = TRUE)

# --- GHG Standard Curves ---
GHG_data <- GHG_standards %>%
  left_join(GC_data, by = c("Sample" = "Sample.Type"))

GHG_data$N2O.Area <- as.numeric(GHG_data$N2O.Area)
GHG_data$N2O_ppm <- as.numeric(gsub(",", "", GHG_data$X.N2O...ppm.))
GHG_data$CH4.Area <- as.numeric(GHG_data$CH4.Area)
GHG_data$CH4_ppm <- as.numeric(gsub(",", "", GHG_data$X.CH4...ppm.))
GHG_data$CO2.Area <- as.numeric(GHG_data$CO2.Area)
GHG_data$CO2_ppm <- as.numeric(gsub(",", "", GHG_data$X.CO2...ppm.))

GHG_data <- GHG_data %>%
  filter(!is.na(N2O.Area) & !is.na(N2O_ppm)) %>%
  filter(!is.na(CH4.Area) & !is.na(CH4_ppm)) %>%
  filter(!is.na(CO2.Area) & !is.na(CO2_ppm))

low_range_standards <- c("N2", "SB1", "SB2", "SB3")

low_range_data <- GHG_data %>% filter(Sample %in% low_range_standards)

N2O_low_curve <- lm(N2O_ppm ~ N2O.Area, data = low_range_data)
N2O_high_curve <- lm(N2O_ppm ~ N2O.Area, data = GHG_data)
CH4_low_curve <- lm(CH4_ppm ~ CH4.Area, data = low_range_data)
CH4_high_curve <- lm(CH4_ppm ~ CH4.Area, data = GHG_data)
CO2_low_curve <- lm(CO2_ppm ~ CO2.Area, data = low_range_data)
CO2_high_curve <- lm(CO2_ppm ~ CO2.Area, data = GHG_data)

# Thresholds for low/high range selection
max_N2O_low_area <- 1500
max_CH4_low_area <- 1000
max_CO2_low_area <- 2000

# Apply calibration to all GC data
GC_data <- GC_data %>%
  mutate(
    O2_concentration = if_else(
      O2.Area <= max_O2_low_area,
      predict(O2_low_curve, newdata = data.frame(O2.Area = O2.Area)),
      predict(O2_high_curve, newdata = data.frame(O2.Area = O2.Area))
    ),
    N2O_concentration = if_else(
      N2O.Area <= max_N2O_low_area,
      predict(N2O_low_curve, newdata = data.frame(N2O.Area = N2O.Area)),
      predict(N2O_high_curve, newdata = data.frame(N2O.Area = N2O.Area))
    ),
    CH4_concentration = if_else(
      CH4.Area <= max_CH4_low_area,
      predict(CH4_low_curve, newdata = data.frame(CH4.Area = CH4.Area)),
      predict(CH4_high_curve, newdata = data.frame(CH4.Area = CH4.Area))
    ),
    CO2_concentration = if_else(
      CO2.Area <= max_CO2_low_area,
      predict(CO2_low_curve, newdata = data.frame(CO2.Area = CO2.Area)),
      predict(CO2_high_curve, newdata = data.frame(CO2.Area = CO2.Area))
    )
  )

# ==============================================================================
# STEP 2: Prepare figure data
# ==============================================================================

# Filter for trunk gas samples
quve_gas <- GC_data %>% filter(!is.na(Lab.ID))
int_gas <- quve_gas %>% filter(Tree.Tissue == "Trunk Gas")
int_gas$Tree.Height <- as.numeric(int_gas$Tree.Height)

# Load mcrA ddPCR data (readr handles encoding/BOM properly)
mcra <- readr::read_csv("data/raw/ddpcr/black_oak_mcrA.csv", show_col_types = FALSE)
# Rename columns by position/pattern to avoid encoding issues with µ
conc_col <- grep("Conc", colnames(mcra), value = TRUE)[1]
height_col <- grep("Height", colnames(mcra), value = TRUE)[1]
mcra <- mcra %>% rename(Height_cm = !!height_col, Conc_copies_uL = !!conc_col)
mcra_clean <- mcra %>% filter(Component %in% c("Heartwood", "Sapwood"))

# Load flux data from goFlux results
flux_data <- read.csv("data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv")

# Clean Height_m: remove non-numeric entries ("Seam", "4 (restarted)") and convert
flux_data$Height_m <- suppressWarnings(as.numeric(flux_data$Height_m))

# Build flux dataframe: extract CH4 flux by height
flux_df <- flux_data %>%
  filter(!is.na(Height_m)) %>%
  select(Height_m, CH4_best.flux) %>%
  rename(height = Height_m, flux = CH4_best.flux)

# Merge CH4 concentration with flux by height
flux_df <- flux_df %>%
  left_join(
    int_gas %>%
      group_by(Tree.Height) %>%
      summarize(ch4 = mean(CH4_concentration, na.rm = TRUE), .groups = "drop"),
    by = c("height" = "Tree.Height")
  )

# Define common height breaks for y-axis
height_breaks <- sort(unique(int_gas$Tree.Height))

# ==============================================================================
# STEP 3: Create 3-panel figure
# ==============================================================================

# Shared theme for publication quality
shared_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0),
    plot.margin = margin(5, 8, 5, 5)
  )

# Panel 1: Height vs CH4 concentration, colored by O2%
p1a <- ggplot(int_gas, aes(x = Tree.Height, y = CH4_concentration)) +
  geom_smooth(se = FALSE, color = "black") +
  geom_point(aes(fill = O2_concentration / 1e6 * 100), size = 3, shape = 21,
             color = "black", stroke = 0.6, alpha = 0.85) +
  shared_theme +
  scale_fill_distiller(palette = "RdYlBu", direction = 1,
                       name = expression(O[2]~"(%)")) +
  coord_flip() +
  ylab(expression(CH[4]~"(ppm)")) +
  xlab("Height (m)") +
  scale_x_continuous(breaks = height_breaks, minor_breaks = NULL) +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top"))

# Panel 2: Height vs mcrA copies/µL (heartwood vs sapwood)
p2 <- ggplot(mcra_clean, aes(y = Conc_copies_uL, x = Height_cm / 100)) +
  geom_smooth(se = FALSE, aes(color = Component)) +
  geom_jitter(size = 3, shape = 21, aes(fill = Component),
              color = "black", stroke = 0.6, alpha = 0.85) +
  shared_theme +
  scale_color_manual(values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4")) +
  scale_fill_manual(values = c("Heartwood" = "#a6611a", "Sapwood" = "#1f78b4")) +
  coord_flip() +
  ylab(expression(italic(mcrA)~"(copies/"*mu*"L)")) +
  xlab("") +
  scale_x_continuous(breaks = height_breaks, minor_breaks = NULL) +
  guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))

# Panel 3: Height vs CH4 flux, colored by internal CH4 concentration
p3 <- ggplot(flux_df, aes(y = flux, x = height, fill = ch4)) +
  geom_smooth(se = FALSE, color = "black") +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.6, alpha = 0.85) +
  shared_theme +
  scale_fill_viridis_c(option = "E", direction = -1,
                       name = expression(CH[4]~"(ppm)")) +
  coord_flip() +
  xlab("") +
  ylab(expression(CH[4]~"Flux")) +
  scale_x_continuous(breaks = height_breaks, minor_breaks = NULL) +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top"))

# Combine — each panel keeps its own legend
fig7 <- p1a + p2 + p3 + plot_layout(nrow = 1)

print(fig7)

ggsave("../../outputs/figures/main/fig7_felled_oak_profiles.png",
       fig7, width = 12, height = 5.1, dpi = 300)
