# ==============================================================================
# outputs.R -- where outputs go: one definition, routed by the file's extension
# ------------------------------------------------------------------------------
#   out_path("canonical_budget.csv")  -> outputs/data/canonical_budget.csv
#   out_path("rf_grouped_cv.txt")     -> outputs/audit/rf_grouped_cv.txt
#   out_path("fig9_budget.png")       -> outputs/figures/generated/fig9_budget.png
#
# Pass `kind` to override the routing (e.g. out_path("x.png", "figures/qc")).
# The directory is created on demand.
# ==============================================================================

OUT_ROOT <- "outputs"

# data   : canonical CSV/RDS products -- the numbers everything else cites
# audit  : human-readable plain-text transcripts
# figures/generated : raw generator output, before the assembler numbers it
OUT_KIND <- c(csv = "data",  rds = "data",  rda = "data",
              txt = "audit", md  = "audit",
              png = "figures/generated", pdf = "figures/generated")

out_path <- function(name, kind = NULL) {
  if (is.null(kind)) {
    ext  <- tolower(sub(".*\\.([A-Za-z0-9]+)$", "\\1", name))
    kind <- OUT_KIND[[ext]]
    if (is.null(kind))
      stop("out_path(): no route for extension '", ext, "' in '", name,
           "'. Add it to OUT_KIND or pass kind= explicitly.", call. = FALSE)
  }
  d <- file.path(OUT_ROOT, kind)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, name)
}
