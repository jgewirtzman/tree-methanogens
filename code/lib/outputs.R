# ==============================================================================
# Where outputs go. One definition, routed by artifact kind.
#
# WHY THIS EXISTS. outputs/revision/ was a flat dump of 196 files -- 70 .png,
# 68 .csv, 51 .txt -- three unrelated kinds in one directory, mirroring the
# code/revision/ mistake on the output side. Worse, 39 scripts wrote through a
# local `out <- "outputs/revision"` variable and 15 of those emitted more than
# one kind through it, so the destination could not be fixed by retargeting a
# variable. Routing on the artifact's own extension handles both the literal and
# the variable forms uniformly.
#
#   out_path("canonical_budget.csv")  -> outputs/data/canonical_budget.csv
#   out_path("rf_grouped_cv.txt")     -> outputs/audit/rf_grouped_cv.txt
#   out_path("fig9_budget.png")       -> outputs/figures/generated/fig9_budget.png
#
# Pass `kind` to override the routing (e.g. out_path("x.png", "figures/qc")).
# The directory is created on demand, so callers no longer need dir.create().
# ==============================================================================

OUT_ROOT <- "outputs"

# data   : canonical CSV/RDS products -- the numbers everything else cites
# audit  : referee-facing plain-text transcripts
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
