# methanogen_families.R -- the 16S families counted as methanogens: the one list,
# sourced by every script that classifies methanogens (figures 5, 6, S5, S21; the
# clade census, Tables S2 and S5, the manuscript statistics, the FAPROTAX dump).
# 09_tables_stats/03_clade-census.R reports which of them occur; all nine do, and no
# other methanogen family is in either 16S run. Anaerobic methanotrophs (ANME, e.g.
# Methanoperedenaceae) are not methanogens; Bathyarchaeia are not included
# (methanogenic potential debated).
METHANOGEN_FAMILIES <- c("Methanobacteriaceae", "Methanomassiliicoccaceae", "Methanoregulaceae",
                         "Methanocellaceae", "Methanosaetaceae", "Methanomicrobiaceae",
                         "Methanosarcinaceae", "Methanomethyliaceae", "Methanocorpusculaceae")
