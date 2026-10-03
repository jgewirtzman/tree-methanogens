# make_figures.R -- stage E only: every figure generator, then the assembler.
# The list and order live in code/pipeline.csv; this is a shortcut for
#   Rscript code/run_all.R --only E
system2("Rscript", c("code/run_all.R", "--only", "E"))
