# ==============================================================================
# tables_read.R -- read the canonical tables (data/tables/) for analysis scripts
# ------------------------------------------------------------------------------
# Readers and reshapes only: no values are filled, imputed or re-derived here. Every
# value comes from one row of one table in data/tables/ (definitions in
# code/tables/dictionary*.csv). Sourced by code/07_molecular/, 08_figures/,
# 09_tables_stats/ and zenodo/ (run from the repository root).
#
#   read_table(name)        data/tables/<name>.csv as a data.frame (character
#                           columns kept as character; "" read as NA)
#   species_name(code)      species code -> binomial (species.csv)
#   SPECIES_NAMES           named vector code -> binomial, all codes in species.csv
#   SPECIES_ABBREV          named vector code -> abbreviated binomial ("A. rubrum")
#   tree_species(tree_id)   tree -> species code (trees.csv)
#   tree_lookup(spelling, campaign, source_file = NULL)
#                           tree_id for a recorded spelling (tree_tags.csv);
#                           reports spellings that match no tree (never drops silently)
#   flux_stem_campaign(campaign)
#                           flux_stem.csv rows of one campaign + species_code and
#                           species (from trees.csv)
#   asv_table(run)          16S ASV counts + taxonomy of one run (as the raw OTU table)
#   asv_reads(run, samples_only)
#                           16S counts without organelle ASVs; tot = the denominator of
#                           every 16S relative abundance
#   survey_2021_trees()     one row per tree measured in the 2021 survey (flux,
#                           DNA, gas, cores or soil holes), wide, in the column names
#                           the analysis scripts use (see the function)
# ==============================================================================

TABLES_DIR <- "data/tables"

read_table <- function(name) {
  f <- file.path(TABLES_DIR, paste0(name, ".csv"))
  if (!file.exists(f)) stop("canonical table not found: ", f, call. = FALSE)
  utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"), check.names = FALSE)
}

.SPECIES <- read_table("species")
SPECIES_NAMES <- stats::setNames(.SPECIES$species, .SPECIES$species_code)
species_name <- function(code) unname(SPECIES_NAMES[as.character(code)])
# abbreviated binomial for labels: "A. rubrum"; named by code like SPECIES_NAMES
SPECIES_ABBREV <- stats::setNames(sub("^(\\S)\\S*\\s+", "\\1. ", .SPECIES$species), .SPECIES$species_code)

.TREES <- NULL
trees_table <- function() {
  if (is.null(.TREES)) .TREES <<- read_table("trees")
  .TREES
}
tree_species <- function(tree_id) {
  tr <- trees_table()
  tr$species_code[match(as.character(tree_id), tr$tree_id)]
}

# tree_id for spellings as recorded in a campaign's files (tree_tags.csv). A spelling
# must resolve to one tree for the campaign (and source file, if given); spellings
# with no tree, or more than one, are reported and returned as NA.
tree_lookup <- function(spelling, campaign, source_file = NULL, what = "") {
  tt <- read_table("tree_tags")
  tt <- tt[tt$campaign %in% campaign, ]
  if (!is.null(source_file)) tt <- tt[tt$source_file %in% source_file, ]
  tt <- unique(tt[, c("tag_recorded", "tree_id")])
  amb <- unique(tt$tag_recorded[duplicated(tt$tag_recorded)])
  tt <- tt[!tt$tag_recorded %in% amb, ]
  sp <- as.character(spelling)
  out <- tt$tree_id[match(sp, tt$tag_recorded)]
  miss <- unique(sp[is.na(out) & !is.na(sp)])
  if (length(miss))
    cat(sprintf("tree_lookup%s: %d spelling(s) with no single tree in tree_tags.csv (%s): %s\n",
                if (nzchar(what)) paste0(" [", what, "]") else "", length(miss),
                paste(campaign, collapse = "/"), paste(utils::head(miss, 30), collapse = ", ")))
  out
}

flux_stem_campaign <- function(campaign) {
  fs <- read_table("flux_stem")
  fs <- fs[fs$campaign %in% campaign, ]
  fs$species_code <- tree_species(fs$tree_id)
  fs$species <- species_name(fs$species_code)
  fs
}

# ------------------------------------------------------------------------------
# survey_2021_trees(): the 2021 survey per tree (one row per tree_id), built from
#   trees.csv          tree_id, species_id (= species_code), species, dbh (cm; dbh_source),
#                      plot (landscape_position: Ridge, Lowland, Lab), dead
#   flux_stem.csv      2021_multiheight: CH4_best.flux_<h>cm, CO2_best.flux_<h>cm,
#                      Temp_Air_<h>cm (chamber temperature) for h = 50, 75, 125, 200,
#                      RootCrown; mean of the deployments at that height
#   soil_readings.csv  the tree's own 2021 survey holes: ORP_mean, OrganicDepth_mean,
#                      MineralDepth_mean (means over holes)
#   attach_soil()      (soil = TRUE) VWC_mean (%) and SoilTemp_mean (degC) at the tree at
#                      its 2021 visit (code/lib/soil_fill.R, unit "tree"); soil_match and
#                      soil_reading_date say how they were found
#   internal_gas.csv   2021_survey whole-stem gas: <gas>_concentration (ppm)
#   wood_cores.csv     <position>_moisture_fresh_percent, _moisture_dry_percent,
#                      _density_final (mean over the tree's cores; filled sections as
#                      in the table, flagged there)
#   ddpcr.csv + samples.csv
#                      ddpcr_<target>_<Compartment>_<loose|strict>, copies per g
#                      (qc_pass rows; mean over a tree's samples and runs of one
#                      compartment). target: mcra (EvaGreen mcrA), mcra_probe, pmoa,
#                      mmox, 16s_arc, 16s_bact. Compartment: Inner, Outer, Organic,
#                      Mineral. loose = copies_per_g; strict = the strict-threshold
#                      concentration through the same conversion (ddpcr_constants.R)
# ------------------------------------------------------------------------------
survey_2021_trees <- function(soil = FALSE) {
  source("code/lib/ddpcr_constants.R", local = TRUE)
  tr <- trees_table()
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  wide <- function(id, key, val) {         # mean of val per id x key, keys as columns
    k <- unique(key); ids <- unique(id)
    m <- sapply(k, function(kk) sapply(ids, function(i) mean_na(val[id == i & key == kk])))
    m <- matrix(m, nrow = length(ids), dimnames = list(NULL, k))
    data.frame(tree_id = ids, m, check.names = FALSE, stringsAsFactors = FALSE)
  }

  # flux (2021 multi-height campaign)
  fx <- read_table("flux_stem"); fx <- fx[fx$campaign == "2021_multiheight", ]
  h <- ifelse(!is.na(fx$height_label) & fx$height_label == "root crown", "RootCrown",
              as.character(fx$height_cm))
  stopifnot(!anyNA(h))
  flux <- Reduce(function(a, b) merge(a, b, by = "tree_id", all = TRUE), list(
    wide(fx$tree_id, paste0("CO2_best.flux_", h, "cm"), fx$CO2_flux_umol_m2_s),
    wide(fx$tree_id, paste0("CH4_best.flux_", h, "cm"), fx$CH4_flux_nmol_m2_s),
    wide(fx$tree_id, paste0("Temp_Air_", h, "cm"), fx$chamber_temp_C)))

  # soil holes at the tree, 2021 survey (tree readings dated 2021; every soil reading is
  # dated). The 2023 field-sheet readings at trees are not survey holes.
  sr <- read_table("soil_readings")
  sr <- sr[!is.na(sr$tree_id) & (is.na(sr$date) | substr(sr$date, 1, 4) == "2021"), ]
  soilh <- data.frame(tree_id = unique(sr$tree_id))
  soilh$ORP_mean          <- sapply(soilh$tree_id, function(i) mean_na(sr$orp_mV[sr$tree_id == i]))
  soilh$OrganicDepth_mean <- sapply(soilh$tree_id, function(i) mean_na(sr$organic_depth_cm[sr$tree_id == i]))
  soilh$MineralDepth_mean <- sapply(soilh$tree_id, function(i) mean_na(sr$mineral_depth_cm[sr$tree_id == i]))

  # internal gas (2021 whole-stem samples)
  g <- read_table("internal_gas"); g <- g[g$campaign == "2021_survey" & !is.na(g$tree_id), ]
  stopifnot(!anyDuplicated(g$tree_id))
  gas <- data.frame(tree_id = g$tree_id, CO2_concentration = g$CO2_ppm, CH4_concentration = g$CH4_ppm,
                    N2O_concentration = g$N2O_ppm, O2_concentration = g$O2_ppm)

  # wood cores
  wc <- read_table("wood_cores")
  wood <- Reduce(function(a, b) merge(a, b, by = "tree_id", all = TRUE), list(
    wide(wc$tree_id, paste0(wc$position, "_moisture_fresh_percent"), wc$moisture_fresh_pct),
    wide(wc$tree_id, paste0(wc$position, "_moisture_dry_percent"), wc$moisture_dry_pct),
    wide(wc$tree_id, paste0(wc$position, "_density_final"), wc$density_g_cm3)))

  # ddPCR, copies per g
  dd <- read_table("ddpcr"); smp <- read_table("samples")
  smp <- smp[smp$campaign == "2021_survey", ]
  dd <- dd[dd$qc_pass & dd$sample_id %in% smp$sample_id, ]
  s <- smp[match(dd$sample_id, smp$sample_id), ]
  comp <- c(inner = "Inner", outer = "Outer", organic = "Organic", mineral = "Mineral")[s$position]
  target <- ifelse(dd$gene == "mcrA" & dd$assay == "probe", "mcra_probe",
            c(mcrA = "mcra", pmoA = "pmoa", mmoX = "mmox", `16S_archaea` = "16s_arc",
              `16S_bacteria` = "16s_bact")[dd$gene])
  strict <- ddpcr_copies_per_g(dd$copies_per_uL_rxn_strict, s$extraction_mass_mg, comp)
  ok <- !is.na(s$tree_id)
  genes <- merge(
    wide(s$tree_id[ok], paste0("ddpcr_", target, "_", comp, "_loose")[ok], dd$copies_per_g[ok]),
    wide(s$tree_id[ok], paste0("ddpcr_", target, "_", comp, "_strict")[ok], strict[ok]),
    by = "tree_id", all = TRUE)

  ids <- sort(unique(c(flux$tree_id, gas$tree_id, wood$tree_id, genes$tree_id, soilh$tree_id)))
  base <- tr[match(ids, tr$tree_id), ]
  out <- data.frame(tree_id = ids, species_id = base$species_code, species = species_name(base$species_code),
                    dbh = base$dbh_cm, dbh_source = base$dbh_source, plot = base$landscape_position,
                    dead = base$dead, stringsAsFactors = FALSE)
  for (d in list(gas, flux, soilh, wood, genes)) out <- merge(out, d, by = "tree_id", all.x = TRUE)

  if (soil) {   # soil temperature and moisture at the tree: code/lib/soil_fill.R only
    if (!exists("attach_soil", mode = "function")) source("code/lib/soil_fill.R")
    # time of the tree's 2021 survey visit: the date its soil holes were read, else its first
    # 2021 deployment (attach_soil takes the tree's own readings within 3 days of it)
    first_start <- tapply(fx$start_time, fx$tree_id, min)
    hole_date <- tapply(sr$date, sr$tree_id, function(d) if (all(is.na(d))) NA_character_ else min(d, na.rm = TRUE))
    out$time <- unname(ifelse(is.na(hole_date[out$tree_id]), first_start[out$tree_id], hole_date[out$tree_id]))
    out <- attach_soil(out, id_col = "tree_id", time_col = "time", unit = "tree")
    out$VWC_mean <- out$soil_vwc_pct; out$SoilTemp_mean <- out$soil_temp_C
    out$time <- NULL
  }
  out
}

# ------------------------------------------------------------------------------
# dna_samples(campaign): one row per DNA sample (samples.csv) of the campaign(s), with
#   species_code (trees.csv), library_id (samples_16S.csv, the sequenced library),
#   facility qPCR: X16S_per_ul, ITS_per_ul (copies per uL extract; qpcr.csv), ITS_undetected,
#   ddPCR per target (mcra_probe, mcra_eva, pmoa, mmox, arc, bact), mean over the sample's
#   runs that passed QC (ddpcr.csv):
#     <target>_g          copies per g (loose threshold; ddpcr.copies_per_g; mass_basis
#                         dry wood / fresh soil)
#     <target>_g_strict   copies per g at the strict threshold (same conversion)
#     <target>_ul         copies per uL of DNA extract (loose; Conc x 25/2.5)
#   material: Wood (2021 inner/outer; felled-oak wood) or Soil (2021 organic/mineral)
#   as samples.material capitalised; core_type: Inner/Outer/Organic/Mineral (2021) or
#   the felled-oak position as recorded.
# ------------------------------------------------------------------------------
dna_samples <- function(campaign = "2021_survey") {
  source("code/lib/ddpcr_constants.R", local = TRUE)
  smp <- read_table("samples"); smp <- smp[smp$campaign %in% campaign, ]
  cap <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  out <- data.frame(sample_id = smp$sample_id, tree_id = smp$tree_id, campaign = smp$campaign,
                    species_code = tree_species(smp$tree_id), material = cap(smp$material),
                    core_type = cap(smp$position), height_cm = smp$height_cm,
                    extraction_mass_mg = smp$extraction_mass_mg, stringsAsFactors = FALSE)
  s16 <- read_table("samples_16S"); s16 <- s16[s16$library_type == "sample", ]
  out$library_id <- s16$library_id[match(out$sample_id, s16$sample_id)]
  q <- read_table("qpcr")
  q16 <- q[q$target == "16S", ]; qits <- q[q$target == "ITS", ]
  out$X16S_per_ul <- q16$copies_per_uL_extract[match(out$sample_id, q16$sample_id)]
  out$ITS_per_ul <- qits$copies_per_uL_extract[match(out$sample_id, qits$sample_id)]
  out$ITS_undetected <- qits$undetected[match(out$sample_id, qits$sample_id)]
  dd <- read_table("ddpcr"); dd <- dd[dd$qc_pass & dd$sample_id %in% out$sample_id, ]
  tg <- ifelse(dd$gene == "mcrA", ifelse(dd$assay == "probe", "mcra_probe", "mcra_eva"),
               c(pmoA = "pmoa", mmoX = "mmox", `16S_archaea` = "arc", `16S_bacteria` = "bact")[dd$gene])
  mass <- smp$extraction_mass_mg[match(dd$sample_id, smp$sample_id)]
  mat <- ifelse(smp$material[match(dd$sample_id, smp$sample_id)] == "soil", "Soil", "Wood")
  dd$g_strict <- ddpcr_copies_per_g(dd$copies_per_uL_rxn_strict, mass, mat)
  dd$ul <- dd$copies_per_uL_rxn * DDPCR_REACTION_UL / DDPCR_TEMPLATE_UL
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  for (t in sort(unique(tg))) {
    k <- tg == t
    for (v in c("copies_per_g", "g_strict", "ul")) {
      col <- paste0(t, c(copies_per_g = "_g", g_strict = "_g_strict", ul = "_ul")[[v]])
      out[[col]] <- as.vector(tapply(dd[[v]][k], dd$sample_id[k], mean_na)[out$sample_id])
    }
  }
  out
}

# ------------------------------------------------------------------------------
# asv_table(run): the 16S ASV table of one sequencing run as a data.frame laid out like
#   the raw OTU_table.txt (row names = feature_id; one integer column per library,
#   named by samples_16S.library_id; then Kingdom ... Species, capitalised as in the raw
#   file). run = "main" (asv_counts_16S.csv) or "felled_oak" (..._felled_oak.csv).
#   The two runs were denoised separately: their feature_ids are not the same sequences.
# ------------------------------------------------------------------------------
asv_table <- function(run = c("main", "felled_oak")) {
  run <- match.arg(run)
  sfx <- if (run == "main") "" else "_felled_oak"
  cnt <- utils::read.csv(file.path(TABLES_DIR, paste0("asv_counts_16S", sfx, ".csv")),
                         check.names = FALSE, stringsAsFactors = FALSE)
  tax <- read_table(paste0("asv_taxonomy_16S", sfx))
  stopifnot(identical(cnt$feature_id, tax$feature_id))
  rownames(cnt) <- cnt$feature_id; cnt$feature_id <- NULL
  tax <- tax[, c("kingdom", "phylum", "class", "order", "family", "genus", "species")]
  names(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  cbind(cnt, tax)
}

# asv_reads(run, samples_only): the 16S counts every relative abundance is computed from.
#   Organelle ASVs (family Mitochondria, order Chloroplast; is_organelle() in
#   code/lib/build_phyloseq.R, the one rule) are removed, so a library's total is its non-organelle
#   reads: the denominator of every 16S relative abundance. samples_only = TRUE keeps the
#   libraries of real samples only (samples_16S.library_type "sample": drops the Zymo mock
#   and the "other" libraries). Returns list(cnt = ASV x library count matrix,
#   tax = Kingdom ... Species of those ASVs, tot = colSums(cnt), the non-organelle reads).
asv_reads <- function(run = c("main", "felled_oak"), samples_only = FALSE) {
  run <- match.arg(run)
  if (!exists("is_organelle", mode = "function")) source("code/lib/build_phyloseq.R")
  a <- asv_table(run)
  taxc <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax <- a[, taxc]; cnt <- as.matrix(a[, setdiff(names(a), taxc), drop = FALSE])
  if (samples_only) {
    L <- read_table("samples_16S"); L <- L[L$run == run, ]
    stopifnot(setequal(L$library_id, colnames(cnt)))
    cnt <- cnt[, colnames(cnt) %in% L$library_id[L$library_type == "sample"], drop = FALSE]
  }
  org <- is_organelle(tax)
  cnt <- cnt[!org, , drop = FALSE]; tax <- tax[!org, , drop = FALSE]
  list(cnt = cnt, tax = tax, tot = colSums(cnt, na.rm = TRUE))
}

# ------------------------------------------------------------------------------
# survey_16s_meta(): per sequenced 2021 DNA sample, row names = 16S library_id (the
#   column names of asv_table("main") and of the PICRUSt2 outputs), in the layout the
#   16S/PICRUSt scripts use: seq_id (= tree_id), core_type, material, species.x (code),
#   X16S_per_ul, ITS_per_ul, and <target>_loose = ddPCR copies per 20 uL well at the loose
#   threshold (QX200 "Copies/20uLWell" = Conc x 20; mean over the sample's QC-passing runs).
# ------------------------------------------------------------------------------
survey_16s_meta <- function() {
  d <- dna_samples("2021_survey"); d <- d[!is.na(d$library_id), ]
  out <- data.frame(seq_id = d$tree_id, core_type = d$core_type, material = d$material,
                    species.x = d$species_code, X16S_per_ul = d$X16S_per_ul, ITS_per_ul = d$ITS_per_ul,
                    row.names = d$library_id, stringsAsFactors = FALSE)
  for (t in c("mcra_probe", "mcra_eva", "pmoa", "mmox", "arc", "bact"))
    out[[paste0(t, "_loose")]] <- d[[paste0(t, "_ul")]] * 2   # per uL extract x 2 = per 20 uL well
  out
}

# ------------------------------------------------------------------------------
# isotope_whole_tree(): the whole-tree d13C set (isotopes.csv in_whole_tree_set; the rule
#   is code/lib/isotope_samples.R, applied when the table is built), one row per tree, in
#   the layout of isotope_whole_tree_samples(): SampleName, tree_id, species (trees.csv
#   species_code), d13CH4, ch4_ppm, d13CO2, co2_ppm.
# isotope_atmosphere(): atmosphere vials for reference lines ("Amb..." runs, CH4 from the
#   floor ISO_CH4_FLOOR to 5 ppm; higher = contaminated), as isotope_atmosphere_samples().
# ------------------------------------------------------------------------------
isotope_whole_tree <- function() {
  iso <- read_table("isotopes"); iso <- iso[iso$in_whole_tree_set, ]
  stopifnot(!anyDuplicated(iso$tree_id), !anyNA(iso$tree_id))
  data.frame(SampleName = iso$sample_name, tree_id = iso$tree_id, species = tree_species(iso$tree_id),
             d13CH4 = iso$d13CH4_permil, ch4_ppm = iso$CH4_ppm, d13CO2 = iso$d13CO2_permil,
             co2_ppm = iso$CO2_ppm, stringsAsFactors = FALSE)
}
isotope_atmosphere <- function(floor_ppm = 1.5, cap_ppm = 5) {
  iso <- read_table("isotopes")
  iso <- iso[iso$sample_type == "atmosphere or room air" & grepl("^Amb", iso$sample_name) &
             !is.na(iso$d13CH4_permil) & iso$CH4_ppm >= floor_ppm & iso$CH4_ppm <= cap_ppm, ]
  data.frame(SampleName = iso$sample_name, d13CH4 = iso$d13CH4_permil, ch4_ppm = iso$CH4_ppm,
             stringsAsFactors = FALSE)
}
