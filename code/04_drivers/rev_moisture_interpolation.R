#!/usr/bin/env Rscript
# ==============================================================================
# rev_moisture_interpolation.R
# ------------------------------------------------------------------------------
# Which interpolator should build the soil-moisture SHAPE?
#
# The pipeline uses akima::interp() at its default linear = TRUE, which fits flat
# triangles between survey points. That is why the Figure 9 soil panel shows
# planar facets meeting at visible creases, worst near the hull edge where the
# triangles get long and thin. No comparison was ever recorded: the file
# code/06_upscale/03_interpolation_methods.R is titled a comparison of Akima, TPS
# and GAM, and contains TPS and GAM implementations, but they were written to
# extend the surface beyond the convex hull rather than chosen on accuracy.
#
# WHAT THE SURFACE HAS TO DO, and it is not what a naive CV rewards. Since
# rev_predict_soil_surface.R now uses a fixed pattern times a monthly level, the
# interpolator sets only the RELATIVE spatial pattern; the mean is imposed by the
# climatology. So the surface must get the large-scale wet/dry structure right
# and be smooth enough not to invent features -- matching any individual survey
# point is not the objective, and an interpolator that nails every point by
# wiggling between them is worse, not better.
#
# Compared here by leave-one-out over the survey points, plus roughness, because
# the two together separate "captures the gradient" from "chases noise".
#
# A CAVEAT ON THE INPUT worth stating in the methods: river points are inserted
# with an assumed VWC of 100%, which is not a measurement. Every method inherits
# that, and it anchors the wet end of the gradient.
#
# Output: outputs/data/moisture_interpolation.csv / .txt
#         outputs/figures/generated/fig_moisture_interpolation.png
# ==============================================================================
suppressPackageStartupMessages({library(dplyr); library(akima); library(fields)
                                library(mgcv); library(ggplot2); library(patchwork)})
set.seed(42)

SRC <- "data/raw/inventory/spatial_data/soil_moisture_20201216.csv"
if (!file.exists(SRC)) stop("missing ", SRC)
D <- read.csv(SRC, stringsAsFactors = FALSE)
xc <- intersect(c("Longitude","lon","x","X"), names(D))[1]
yc <- intersect(c("Latitude","lat","y","Y"), names(D))[1]
D <- D[is.finite(D[[xc]]) & is.finite(D[[yc]]) & is.finite(D$VWC), ]
D <- data.frame(x = D[[xc]], y = D[[yc]], v = D$VWC)
D <- D[!duplicated(round(cbind(D$x, D$y), 8)), ]
cat(sprintf("survey points: %d | VWC %.1f - %.1f (mean %.1f)\n",
            nrow(D), min(D$v), max(D$v), mean(D$v)))
# work in metres so distances and roughness are interpretable
mlon <- 111320*cos(mean(D$y)*pi/180); mlat <- 111132
D$mx <- (D$x - min(D$x))*mlon; D$my <- (D$y - min(D$y))*mlat
cat(sprintf("extent %.0f x %.0f m | median nearest-neighbour spacing %.1f m\n\n",
            diff(range(D$mx)), diff(range(D$my)),
            median(sapply(seq_len(nrow(D)), function(i)
              min(sqrt((D$mx[-i]-D$mx[i])^2 + (D$my[-i]-D$my[i])^2))))))

# ---- interpolators, each a function(train, newxy) -> values ------------------
METH <- list(
  `akima linear (current)` = function(tr, nw)
    akima::interpp(tr$mx, tr$my, tr$v, nw$mx, nw$my, linear = TRUE, extrap = FALSE)$z,
  `akima spline` = function(tr, nw)
    akima::interpp(tr$mx, tr$my, tr$v, nw$mx, nw$my, linear = FALSE, extrap = TRUE)$z,
  `nearest neighbour` = function(tr, nw)
    tr$v[apply(fields::rdist(cbind(nw$mx, nw$my), cbind(tr$mx, tr$my)), 1, which.min)],
  `IDW power 1` = function(tr, nw) idw(tr, nw, 1),
  `IDW power 2` = function(tr, nw) idw(tr, nw, 2),
  `thin-plate spline` = function(tr, nw) {
    m <- suppressWarnings(fields::Tps(cbind(tr$mx, tr$my), tr$v))
    as.numeric(predict(m, cbind(nw$mx, nw$my)))
  },
  `GAM s(x,y)` = function(tr, nw) {
    m <- mgcv::gam(v ~ s(mx, my, k = min(60, nrow(tr) - 5)), data = tr, method = "REML")
    as.numeric(predict(m, newdata = nw))
  })
idw <- function(tr, nw, p) {
  d <- fields::rdist(cbind(nw$mx, nw$my), cbind(tr$mx, tr$my)); d[d < 1e-9] <- 1e-9
  w <- 1/d^p; as.numeric((w %*% tr$v)/rowSums(w))
}

# ---- leave-one-out ----------------------------------------------------------
cat("running leave-one-out over survey points...\n")
res <- lapply(names(METH), function(nm) {
  p <- rep(NA_real_, nrow(D))
  for (i in seq_len(nrow(D)))
    p[i] <- tryCatch(METH[[nm]](D[-i, ], D[i, , drop = FALSE]), error = function(e) NA_real_)
  ok <- is.finite(p)
  data.frame(method = nm, n_pred = sum(ok),
             rmse = sqrt(mean((D$v[ok]-p[ok])^2)),
             mae  = mean(abs(D$v[ok]-p[ok])),
             r    = suppressWarnings(cor(D$v[ok], p[ok])),
             bias = mean(p[ok]) - mean(D$v[ok]))
})
R <- bind_rows(res)

# ---- roughness of the fitted surface on a regular grid ----------------------
gr <- expand.grid(mx = seq(min(D$mx), max(D$mx), length.out = 80),
                  my = seq(min(D$my), max(D$my), length.out = 80))
surf <- list()
for (nm in names(METH)) {
  z <- tryCatch(METH[[nm]](D, gr), error = function(e) rep(NA_real_, nrow(gr)))
  surf[[nm]] <- z
  M <- matrix(z, 80, 80)
  # mean absolute second difference: how much the surface bends per cell
  r1 <- mean(abs(diff(M, differences = 2)), na.rm = TRUE)
  r2 <- mean(abs(diff(t(M), differences = 2)), na.rm = TRUE)
  R$roughness[R$method == nm] <- mean(c(r1, r2), na.rm = TRUE)
  R$pct_na[R$method == nm] <- 100*mean(!is.finite(z))
}
R <- R %>% arrange(rmse)
cat("\n=== LEAVE-ONE-OUT OVER SURVEY POINTS ===\n")
cat("roughness = mean |2nd difference| on an 80x80 grid; lower is smoother\n")
cat("pct_na = share of grid the method cannot fill (hull limitation)\n\n")
print(as.data.frame(R %>% mutate(across(where(is.numeric), ~round(.x, 3)))), row.names = FALSE)

# ---- figure ------------------------------------------------------------------
pl <- lapply(names(METH), function(nm) {
  g <- gr; g$z <- surf[[nm]]
  ggplot(g, aes(mx, my, fill = z)) + geom_raster() +
    geom_point(data = D, aes(mx, my), inherit.aes = FALSE, size = 0.25, colour = "grey20", alpha = 0.6) +
    scale_fill_viridis_c(option = "mako", na.value = "white", name = "VWC") +
    coord_equal() + labs(title = nm, x = NULL, y = NULL) +
    theme_bw(base_size = 8) + theme(legend.position = "none",
      axis.text = element_blank(), axis.ticks = element_blank())
})
ggsave("outputs/figures/generated/fig_moisture_interpolation.png",
       wrap_plots(pl, ncol = 4), width = 13, height = 6.5, dpi = 180, bg = "white")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
write.csv(R, "outputs/data/moisture_interpolation.csv", row.names = FALSE)
con <- file("outputs/audit/moisture_interpolation.txt", "w")
cat(sprintf("MOISTURE INTERPOLATION COMPARISON\nbuilt %s\n%d survey points\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), nrow(D)), file = con)
utils::write.table(round(as.data.frame(R %>% select(-method)), 4), con, sep = "\t", quote = FALSE)
close(con)
cat("\nwritten: outputs/moisture_interpolation.{csv,txt} and fig_moisture_interpolation.png\n")
