# ==============================================================================
# REVISION — assemble final Fig 7 with corrected panels b and c.
# Reuses the EXACT panel objects by sourcing the original figure script (panel a,
# unchanged) and the two revision fix scripts (panel b = copies/g, panel c = flux
# unit). No existing script is edited. NEW file.
# Output: revision/outputs/fig7_assembled.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

# Original script builds p1a (panel a, unchanged) + p2/p3 (old) and writes old fig7.
source("code/06_figures/09_felled_oak_profiles.R")
# Corrected panels (each defines its panel object; harmlessly re-saves its own PNG).
source("revision/analysis/R_fig7b_copies_per_gram.R")   # defines p_g  (copies/g)
source("revision/analysis/R_fig7c_flux_unit.R")         # redefines p3 (flux unit fixed)

# strip development titles carried in from the standalone fix scripts
p_g <- p_g + labs(title = NULL)
p3  <- p3  + labs(title = NULL)

fig7 <- p1a + p_g + p3 + plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))
ggsave("revision/outputs/fig7_assembled.png", fig7, width = 12, height = 5.1, dpi = 300, bg = "white")
cat("Wrote revision/outputs/fig7_assembled.png\n")
