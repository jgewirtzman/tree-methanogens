# ==============================================================================
# flux_transform.R -- the inverse hyperbolic sine used for CH4 flux
# ------------------------------------------------------------------------------
# asinh(x / c) is linear for |x| << c and logarithmic for |x| >> c, and keeps the
# sign, so uptake (negative flux) and zero stay on the axis. The cofactor c is
# 0.1 nmol m-2 s-1 for flux throughout (Figs 1-3, 7, 9, S19, S22, S27; the gene
# models of Methods S4); other cofactors are passed explicitly.
#
#   flux_asinh(x)        asinh(x / c)                 (natural-log units)
#   flux_asinh10(x)      asinh(x / c) / log(10)       (log10 units, for gene models)
#   flux_asinh_trans()   the ggplot2 scale transform of flux_asinh (base10 = TRUE:
#                        of flux_asinh10)
# ==============================================================================
FLUX_COFACTOR <- 0.1

flux_asinh   <- function(x, cofactor = FLUX_COFACTOR) asinh(x / cofactor)
flux_asinh10 <- function(x, cofactor = FLUX_COFACTOR) asinh(x / cofactor) / log(10)

flux_asinh_trans <- function(cofactor = FLUX_COFACTOR, base10 = FALSE) {
  if (base10)
    scales::trans_new("arcsinh", function(x) asinh(x / cofactor) / log(10),
                      function(x) cofactor * sinh(x * log(10)))
  else
    scales::trans_new("arcsinh", function(x) asinh(x / cofactor), function(x) cofactor * sinh(x))
}
