#!/usr/bin/env Rscript
# ==============================================================================
# methodsS3_extraction_recovery.R -- where target copies are lost, core to eluate
# (schematic for SI Methods S3, §8 Recovery)
# ------------------------------------------------------------------------------
# Follows 100 target copies spiked onto ground wood through the extraction, using
# the method study's measurements (Arnold et al. 2024, Methods Ecol. Evol.):
#   - lysis-buffer control: gBlocks added to lysis buffer alone recover 30.4%
#   - field wood: recovery averaged 81.7% of that control (ground wood spiked
#     immediately before extraction)
#   - lysate processed: 250 of the 800 µL lysate (method-study Dryad sheets:
#     "Super Attempted" = 250 µL)
#   - freeze-drying + cryo-grinding (spiked sterile dowels, relative to the control,
#     averaged over gBlocks and cells): 43.0% spiked before extraction -> 19.2% before
#     freeze-drying, so about 45% survive
# Derived: cleanup keeps 30.4 / 31.25 = 97%; wood releases 81.7%. Of 100 copies in a
# core, 44.7 reach the ground powder and 0.447 x 0.817 x 0.3125 x 0.973 x 100 = 11.1
# reach the eluate. Dowels ground less finely than cores (paper), so the first step
# may overstate field losses.
# Writes outputs/figures/generated/methodsS3_extraction_recovery.png
# ==============================================================================
suppressMessages(library(ggplot2))
source("code/lib/outputs.R")

CONTROL   <- 0.304          # lysis-buffer control recovery
REL_WOOD  <- 0.817          # field wood relative to the control
TRANSFER  <- 250 / 800      # fraction of lysate processed
CLEANUP   <- CONTROL / TRANSFER
# freeze-drying + cryo-grinding (sterile dowels; recovery relative to the control,
# averaged over gBlocks and cultured cells): 43.0% spiked before extraction -> 19.2%
# spiked before freeze-drying
GRIND     <- 19.2 / 43.0
ground      <- 100 * GRIND
released    <- ground * REL_WOOD
transferred <- released * TRANSFER
eluate      <- transferred * CLEANUP
lost <- c(grind = 100 - ground, wood = ground - released, lysate = released - transferred,
          cleanup = transferred - eluate)

band <- function(x0, x1, top0, bot0, top1, bot1, n = 60) {
  s <- (1 - cos(seq(0, pi, length.out = n))) / 2
  x <- x0 + (x1 - x0) * seq(0, 1, length.out = n)
  data.frame(x = c(x, rev(x)), y = c(top0 + (top1 - top0) * s, rev(bot0 + (bot1 - bot0) * s)))
}
W <- 0.05; GAP <- 6
X <- 0:4
main <- c(100, ground, released, transferred, eluate)
nodes <- data.frame(x = X, top = 100, h = main, kind = "kept")
lossn <- data.frame(x = X[2:5], top = 100 - main[2:5] - GAP, h = pmax(lost, 0.8), kind = "lost")
nodes <- rbind(nodes, lossn); nodes$bot <- nodes$top - nodes$h
flows <- do.call(rbind, lapply(1:4, function(i) {
  keep <- cbind(band(X[i] + W, X[i + 1], 100, 100 - main[i + 1], 100, 100 - main[i + 1]), id = 2 * i - 1, kind = "kept")
  lo   <- cbind(band(X[i] + W, X[i + 1], 100 - main[i + 1], 100 - main[i], lossn$top[i], lossn$top[i] - lossn$h[i]),
                id = 2 * i, kind = "lost")
  rbind(keep, lo)
}))

KEPT <- "#3a9e8c"; LOST <- "#b4b2a9"
f1 <- function(v) formatC(v, format = "f", digits = 1)
heads <- data.frame(x = X, y = 108,
  lab = c("Wood core\n100 copies",
          sprintf("Ground powder\n%s", f1(ground)),
          sprintf("Lysate, 800 µL\n%s released", f1(released)),
          sprintf("250 µL to cleanup\n%s", f1(transferred)),
          sprintf("Eluate, 75 µL\n%s measured", f1(eluate))))
losses <- data.frame(
  x = X[2:5] + W + 0.04,
  y = lossn$top - pmax(lossn$h / 2, 3),
  lab = c(sprintf("Freeze-drying and grinding\n%s (dowel experiment)", f1(lost["grind"])),
          sprintf("Kept by the wood\n%s", f1(lost["wood"])),
          sprintf("Lysate not processed\n%s (550 of 800 µL)", f1(lost["lysate"])),
          sprintf("Cleanup loss\n%s", f1(lost["cleanup"]))))

p <- ggplot() +
  geom_polygon(data = flows, aes(x, y, group = id, fill = kind), alpha = 0.45, colour = NA) +
  geom_rect(data = nodes, aes(xmin = x, xmax = x + W, ymin = bot, ymax = top, fill = kind), colour = NA) +
  geom_text(data = heads, aes(x + W / 2, y, label = lab), size = 3.1, lineheight = 0.95, vjust = 0) +
  geom_text(data = losses, aes(x, y, label = lab), size = 2.9, hjust = 0, lineheight = 0.95, colour = "grey30") +
  scale_fill_manual(values = c(kept = KEPT, lost = LOST),
                    labels = c(kept = "Copies carried forward", lost = "Copies lost"), name = NULL) +
  scale_x_continuous(limits = c(-0.3, 4.85)) +
  scale_y_continuous(limits = c(min(nodes$bot) - 4, 122)) +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom", plot.margin = margin(6, 6, 6, 6))
ggsave(out_path("methodsS3_extraction_recovery.png"), p, width = 8.2, height = 4.2, dpi = 300, bg = "white")
cat(sprintf("ground %.1f | released %.1f | transferred %.1f | eluate %.1f | cleanup %.3f\n", ground, released, transferred, eluate, CLEANUP))
