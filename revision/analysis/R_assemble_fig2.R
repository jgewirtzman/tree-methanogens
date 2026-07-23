# ==============================================================================
# REVISION — assemble final Fig 2 with corrected panel a (consistent axis).
# Sources the original figure script (unchanged panels b = p_middle, c = p_bottom)
# and the revision fix (panel a = object `p`, shared signed pseudo-log axis).
# No existing script is edited. NEW file.
# Output: revision/outputs/fig2_assembled.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

# original builds p_top_gg (old panel a) + p_middle + p_bottom, writes old fig2
source("code/01_flux_processing/static/04_height_effect_analysis.R")
# corrected panel a -> object `p`
source("revision/analysis/R_fig2a_consistent_axis.R")

fig2 <- (p / p_middle / p_bottom) +
  plot_layout(heights = c(5, 1.2, 1.2)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                  theme = theme(plot.tag = element_text(size = 11, face = "bold")))
ggsave("revision/outputs/fig2_assembled.png", fig2, width = 7, height = 7.5,
       units = "in", dpi = 300, bg = "white")
cat("Wrote revision/outputs/fig2_assembled.png\n")
