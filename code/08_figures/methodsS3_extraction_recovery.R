#!/usr/bin/env Rscript
# ==============================================================================
# methodsS3_extraction_recovery.R -- where target copies are lost, core to eluate
# (schematic for SI Methods S3, §8 Recovery)
# ------------------------------------------------------------------------------
# Follows target copies from the wood core to the eluate, using the method study's
# measurements (Arnold et al. 2024, Methods Ecol. Evol.; bioRxiv 10.1101/2023.11.14.567064):
#   - lysis-buffer control: gBlocks added to lysis buffer alone recover 30.4%
#   - field wood: recovery averaged 81.7% of that control (ground wood spiked
#     immediately before extraction)
#   - lysate processed: 250 of the 800 µL lysate (method-study Dryad sheets:
#     "Super Attempted" = 250 µL; 127-238 µL actually taken)
# Derived: cleanup keeps 30.4 / 31.25 = 97%. Of 100 copies in the ground powder,
# 81.7 are released, 25.5 processed and 24.8 reach the eluate.
# Freeze-drying and grinding lose copies before this (about half in spiked dowels),
# but the dowels ground less finely than cores, so that step is drawn as an unfilled
# dashed outline without a number.
# Writes outputs/figures/generated/methodsS3_extraction_recovery.png
# ==============================================================================
suppressMessages(library(ggplot2))
source("code/lib/outputs.R")

CONTROL   <- 0.304          # lysis-buffer control recovery
REL_WOOD  <- 0.817          # field wood relative to the control
TRANSFER  <- 250 / 800      # fraction of lysate processed
CLEANUP   <- CONTROL / TRANSFER
# Freeze-drying + cryo-grinding lose copies before extraction (spiked dowels: about
# half), but how well that carries over to field cores is uncertain, so that step is
# drawn as an unfilled dashed outline with no number. Quantities start at the powder.
released    <- 100 * REL_WOOD
transferred <- released * TRANSFER
eluate      <- transferred * CLEANUP
lost <- c(wood = 100 - released, lysate = released - transferred, cleanup = transferred - eluate)

band <- function(x0, x1, top0, bot0, top1, bot1, n = 60) {
  s <- (1 - cos(seq(0, pi, length.out = n))) / 2
  x <- x0 + (x1 - x0) * seq(0, 1, length.out = n)
  data.frame(x = c(x, rev(x)), y = c(top0 + (top1 - top0) * s, rev(bot0 + (bot1 - bot0) * s)))
}
W <- 0.05; GAP <- 6; UNK <- 45      # UNK: drawn height of the unquantified loss (not data)
X <- 0:4
main <- c(100, released, transferred, eluate)          # at X[2:5]
nodes <- data.frame(x = X[2:5], top = 100, h = main, kind = "kept")
lossn <- data.frame(x = X[3:5], top = 100 - main[2:4] - GAP, h = pmax(lost, 0.8), kind = "lost")
nodes <- rbind(nodes, lossn); nodes$bot <- nodes$top - nodes$h
flows <- do.call(rbind, lapply(1:3, function(i) {
  xa <- X[i + 1] + W; xb <- X[i + 2]
  keep <- cbind(band(xa, xb, 100, 100 - main[i + 1], 100, 100 - main[i + 1]), id = 2 * i - 1, kind = "kept")
  lo   <- cbind(band(xa, xb, 100 - main[i + 1], 100 - main[i], lossn$top[i], lossn$top[i] - lossn$h[i]),
                id = 2 * i, kind = "lost")
  rbind(keep, lo)
}))
# upstream, unquantified: core (dashed) -> powder (kept, outlined) and -> loss (outlined)
core   <- data.frame(xmin = X[1], xmax = X[1] + W, ymin = 100 - 100 - UNK, ymax = 100)
up_keep <- band(X[1] + W, X[2], 100, 0, 100, 0)
up_lost <- band(X[1] + W, X[2], 0, -UNK, -GAP, -GAP - UNK)
up_node <- data.frame(xmin = X[2], xmax = X[2] + W, ymin = -GAP - UNK, ymax = -GAP)

KEPT <- "#3a9e8c"; LOST <- "#b4b2a9"
f1 <- function(v) formatC(v, format = "f", digits = 1)
heads <- data.frame(x = X, y = 108,
  lab = c("Wood core",
          "Ground powder\n100 copies",
          sprintf("Lysate, 800 µL\n%s released", f1(released)),
          sprintf("250 µL to cleanup\n%s", f1(transferred)),
          sprintf("Eluate, 75 µL\n%s measured", f1(eluate))))
losses <- data.frame(
  x = c(X[2] + W + 0.04, X[3:5] + W + 0.04),
  y = c(-GAP - UNK / 2, lossn$top - pmax(lossn$h / 2, 3)),
  lab = c("Freeze-drying and grinding\nsize uncertain (about half\nin spiked dowels)",
          sprintf("Kept by the wood\n%s", f1(lost["wood"])),
          sprintf("Lysate not processed\n%s (550 of 800 µL;\n127–238 µL taken in practice)", f1(lost["lysate"])),
          sprintf("Cleanup loss\n%s", f1(lost["cleanup"]))))

p <- ggplot() +
  geom_polygon(data = up_keep, aes(x, y), fill = NA, colour = KEPT, linetype = "dashed", linewidth = 0.35) +
  geom_polygon(data = up_lost, aes(x, y), fill = NA, colour = "grey45", linetype = "dashed", linewidth = 0.35) +
  geom_rect(data = core, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, colour = KEPT, linetype = "dashed", linewidth = 0.35) +
  geom_rect(data = up_node, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, colour = "grey45", linetype = "dashed", linewidth = 0.35) +
  geom_polygon(data = flows, aes(x, y, group = id, fill = kind), alpha = 0.45, colour = NA) +
  geom_rect(data = nodes, aes(xmin = x, xmax = x + W, ymin = bot, ymax = top, fill = kind), colour = NA) +
  geom_text(data = heads, aes(x + W / 2, y, label = lab), size = 3.1, lineheight = 0.95, vjust = 0) +
  geom_text(data = losses, aes(x, y, label = lab), size = 2.8, hjust = 0, lineheight = 0.95, colour = "grey30") +
  scale_fill_manual(values = c(kept = KEPT, lost = LOST),
                    labels = c(kept = "Copies carried forward", lost = "Copies lost"), name = NULL) +
  scale_x_continuous(limits = c(-0.3, 5.0)) +
  scale_y_continuous(limits = c(-GAP - UNK - 6, 122)) +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom", plot.margin = margin(6, 6, 6, 6))
ggsave(out_path("methodsS3_extraction_recovery.png"), p, width = 8.2, height = 4.2, dpi = 300, bg = "white")
cat(sprintf("released %.1f | transferred %.1f | eluate %.1f | cleanup %.3f\n", released, transferred, eluate, CLEANUP))
