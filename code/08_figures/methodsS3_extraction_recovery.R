#!/usr/bin/env Rscript
# ==============================================================================
# methodsS3_extraction_recovery.R -- where target copies go during DNA extraction
# (schematic for SI Methods S3, §8 Recovery)
# ------------------------------------------------------------------------------
# Follows 100 target copies spiked onto ground wood through the extraction, using
# the method study's measurements (Arnold et al. 2024, Methods Ecol. Evol.):
#   - lysis-buffer control: gBlocks added to lysis buffer alone recover 30.4%
#   - field wood: recovery averaged 81.7% of that control (ground wood spiked
#     immediately before extraction)
#   - lysate processed: 250 of the 800 µL lysate (method-study Dryad sheets:
#     "Super Attempted" = 250 µL)
# Derived: cleanup keeps 30.4 / 31.25 = 97%; wood releases 81.7%. The eluate then
# holds 81.7 x 0.3125 x 0.973 = 24.8 of 100 copies. Freeze-drying and grinding
# losses (measured separately on spiked dowels) come before this and are not drawn.
# Writes outputs/figures/generated/methodsS3_extraction_recovery.png
# ==============================================================================
suppressMessages(library(ggplot2))
source("code/lib/outputs.R")

CONTROL   <- 0.304          # lysis-buffer control recovery
REL_WOOD  <- 0.817          # field wood relative to the control
TRANSFER  <- 250 / 800      # fraction of lysate processed
CLEANUP   <- CONTROL / TRANSFER
released    <- 100 * REL_WOOD
transferred <- released * TRANSFER
eluate      <- transferred * CLEANUP
lost <- c(wood = 100 - released, lysate = released - transferred, cleanup = transferred - eluate)

band <- function(x0, x1, top0, bot0, top1, bot1, n = 60) {
  s <- (1 - cos(seq(0, pi, length.out = n))) / 2
  x <- x0 + (x1 - x0) * seq(0, 1, length.out = n)
  data.frame(x = c(x, rev(x)), y = c(top0 + (top1 - top0) * s, rev(bot0 + (bot1 - bot0) * s)))
}
W <- 0.05; GAP <- 8
X <- c(0, 1, 2, 3)
# main nodes stack from the top (y = 100) down; loss nodes sit below the next main node
nodes <- data.frame(
  x    = c(X[1], X[2], X[3], X[4], X[2], X[3], X[4]),
  top  = c(100, 100, 100, 100, 100 - released - GAP, 100 - transferred - GAP, 100 - eluate - GAP),
  h    = c(100, released, transferred, eluate, lost["wood"], lost["lysate"], max(lost["cleanup"], 0.8)),
  kind = c(rep("kept", 4), rep("lost", 3)))
nodes$bot <- nodes$top - nodes$h
flows <- rbind(
  cbind(band(X[1] + W, X[2], 100, 100 - released, 100, 100 - released), id = 1, kind = "kept"),
  cbind(band(X[1] + W, X[2], 100 - released, 0, nodes$top[5], nodes$bot[5]), id = 2, kind = "lost"),
  cbind(band(X[2] + W, X[3], 100, 100 - transferred, 100, 100 - transferred), id = 3, kind = "kept"),
  cbind(band(X[2] + W, X[3], 100 - transferred, 100 - released, nodes$top[6], nodes$bot[6]), id = 4, kind = "lost"),
  cbind(band(X[3] + W, X[4], 100, 100 - eluate, 100, 100 - eluate), id = 5, kind = "kept"),
  cbind(band(X[3] + W, X[4], 100 - eluate, 100 - transferred, nodes$top[7], nodes$top[7] - max(lost["cleanup"], 0.8)), id = 6, kind = "lost"))

KEPT <- "#3a9e8c"; LOST <- "#b4b2a9"
f1 <- function(v) formatC(v, format = "f", digits = 1)
heads <- data.frame(x = X, y = 108,
  lab = c("Spiked onto ground wood\n100 copies",
          sprintf("Lysate, 800 µL\n%s released", f1(released)),
          sprintf("250 µL to cleanup\n%s", f1(transferred)),
          sprintf("Eluate, 75 µL\n%s measured", f1(eluate))))
losses <- data.frame(
  x = c(X[2], X[3], X[4]) + W + 0.04,
  y = c(nodes$top[5] - lost["wood"] / 2, nodes$top[6] - lost["lysate"] / 2, nodes$top[7] - 2),
  lab = c(sprintf("Kept by the wood\n%s", f1(lost["wood"])),
          sprintf("Lysate not processed\n%s (550 of 800 µL)", f1(lost["lysate"])),
          sprintf("Cleanup loss\n%s", f1(lost["cleanup"]))))

p <- ggplot() +
  geom_polygon(data = flows, aes(x, y, group = id, fill = kind), alpha = 0.45, colour = NA) +
  geom_rect(data = nodes, aes(xmin = x, xmax = x + W, ymin = bot, ymax = top, fill = kind), colour = NA) +
  geom_text(data = heads, aes(x + W / 2, y, label = lab), size = 3.1, lineheight = 0.95, vjust = 0) +
  geom_text(data = losses, aes(x, y, label = lab), size = 2.9, hjust = 0, lineheight = 0.95, colour = "grey30") +
  scale_fill_manual(values = c(kept = KEPT, lost = LOST),
                    labels = c(kept = "Copies carried forward", lost = "Copies lost"), name = NULL) +
  scale_x_continuous(limits = c(-0.25, 3.75)) +
  scale_y_continuous(limits = c(min(nodes$bot) - 4, 122)) +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom", plot.margin = margin(6, 6, 6, 6))
ggsave(out_path("methodsS3_extraction_recovery.png"), p, width = 7.2, height = 3.9, dpi = 300, bg = "white")
cat(sprintf("released %.1f | transferred %.1f | eluate %.1f | cleanup %.3f\n", released, transferred, eluate, CLEANUP))
