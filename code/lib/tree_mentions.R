# ==============================================================================
# tree_mentions.R -- every spelling of every tree in the files the pipeline reads
# ------------------------------------------------------------------------------
# tree_mentions() returns one row per (source_file, source_column, tag_recorded,
# species_recorded, dbh_recorded_cm) with what that record says about the tree:
#   kind      census_tag (numeric plot census tag), name (2021/2024 field name),
#             untag (monthly untagged stem), untag2023, ash2023, felled
#   key       identity before tag-collision resolution: the tag, the NCBI-style
#             name (code/lib/tree_ids.R), or the label
#   measured  TRUE where the source is a measurement of the tree (flux, internal
#             gas, wood core, DNA, isotope); FALSE for tree lists and soil means
# plus species_code, dbh_cm, date, landscape_position, GPS latitude/longitude and
# dead evidence where the source carries them. Spellings that name no single tree
# (generic "Untagged" rows of the raw monthly sheet, 16S MISC libraries) are
# returned with key = NA and a note.
# ==============================================================================
suppressPackageStartupMessages({library(dplyr); library(readxl)})
source("code/lib/tree_ids.R"); source("code/lib/monthly_cols.R")

.num <- function(x) suppressWarnings(as.numeric(x))
.rd  <- function(f, ...) utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, ...)
.m   <- function(campaign, file, column, tag, ...) {
  d <- data.frame(campaign = campaign, source_file = file, source_column = column,
                  tag_recorded = as.character(tag), stringsAsFactors = FALSE)
  extra <- list(...)
  for (k in names(extra)) d[[k]] <- rep_len(extra[[k]], nrow(d))
  d
}
.first_line <- function(x) trimws(gsub("_x000D_|\r", "", sub("\n.*", "", x)))
.is_tag <- function(x) grepl("^[0-9]+(\\.0)?$", trimws(x))

MON <- "monthly_2020_2021"; F21 <- "2021_multiheight"; S21 <- "2021_survey"
C23 <- "2023_cross_species"; OAK <- "2022_felled_oak"; G24 <- "2024_heartwood_sapwood"
FELLED_OAK_ID <- "QUVE_felled_2022"

tree_mentions <- function() {
  M <- list()

  # ---- monthly survey 2020-21 ------------------------------------------------------
  ym <- .rd("data/raw/inventory/spatial_data/YM_trees_measured.csv.csv", fileEncoding = "UTF-8-BOM")
  ym$Label <- trimws(ym$Label)
  ym$tag <- ifelse(.is_tag(ym$Sample), ym$Sample, NA)          # Sample is the tag (Label 3355 is 3359's row)
  ut <- .rd("data/processed/flux/untagged_rescue/untagged_monthly_trees.csv")
  ut$site_sp <- paste(ut$site, ut$species)
  ym$untag <- ifelse(is.na(ym$tag), ut$tree_id[match(paste(ym$Site, ym$Species), ut$site_sp)], NA)
  stopifnot(all(!is.na(ym$tag) | !is.na(ym$untag)))
  ym$key <- ifelse(!is.na(ym$tag), ym$tag, ym$untag)
  ym_att <- ym %>% transmute(key, species_code = fix_species_code(Species), dbh_cm = .num(D_stem),
                             gps_latitude = .num(Latitude), gps_longitude = .num(Longitude))
  y1 <- ym[ym$Label != "", ]
  M$ym_label <- .m(MON, "data/raw/inventory/spatial_data/YM_trees_measured.csv.csv", "Label", y1$Label,
                   species_recorded = y1$Species, dbh_recorded_cm = y1$D_stem, key = y1$key, measured = FALSE)
  M$ym_sample <- .m(MON, "data/raw/inventory/spatial_data/YM_trees_measured.csv.csv", "Sample", ym$Sample,
                    species_recorded = ym$Species, dbh_recorded_cm = ym$D_stem, key = ym$key, measured = FALSE)

  mon <- monthly_plain(.rd("data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv"))
  names(mon)[names(mon) == "Plot Tag"] <- "Plot.Tag"; names(mon)[names(mon) == "Plot Letter"] <- "Plot.Letter"
  mon_att <- mon %>% group_by(key = as.character(Plot.Tag)) %>%
    summarise(landscape_position = names(sort(table(Plot.Letter), decreasing = TRUE))[1],
              date = min(as.Date(substr(Date, 1, 10))), dead_flag = any(dead %in% TRUE), .groups = "drop")
  M$mon <- .m(MON, "data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv", "Plot.Tag",
              mon_att$key, key = mon_att$key, measured = TRUE)
  m0 <- .rd("data/processed/flux/semirigid_tree_final_complete_dataset.csv")
  M$mon0 <- .m(MON, "data/processed/flux/semirigid_tree_final_complete_dataset.csv", "Plot Tag",
               unique(as.character(m0$`Plot Tag`)), measured = TRUE)

  tf <- suppressMessages(read_excel("data/raw/field_data/ipad_data/treeflux_total.xlsx", sheet = "Sheet1"))
  tf$PT <- trimws(as.character(tf$`Plot Tag`))
  tf_tags <- unique(tf$PT[!is.na(tf$PT)])
  k <- ifelse(.is_tag(tf_tags), sub("\\.0$", "", tf_tags), NA)
  M$tf <- .m(MON, "data/raw/field_data/ipad_data/treeflux_total.xlsx", "Plot Tag", tf_tags, key = k, measured = TRUE,
             note = ifelse(is.na(k), "untagged stems are identified per deployment by the rescue (code/02_flux/rescue/)", NA))
  # dead notes on tagged trees: visits with "dead|snag" over all visits
  tfd <- tf %>% filter(.is_tag(PT)) %>% mutate(key = sub("\\.0$", "", PT)) %>% group_by(key) %>%
    summarise(n_dead = sum(grepl("dead|snag", Notes, ignore.case = TRUE)), n = n(), .groups = "drop") %>%
    filter(n_dead > 0)

  M$ut <- .m(MON, "data/processed/flux/untagged_rescue/untagged_monthly_trees.csv", "tree_id", ut$tree_id,
             key = ut$tree_id, measured = TRUE)
  uf <- .rd("data/raw/flux_windows/untagged/untagged_fluxes.csv") %>% distinct(Sample, site, species, Dstem)
  M$uf <- .m(MON, "data/raw/flux_windows/untagged/untagged_fluxes.csv", "Sample", uf$Sample,
             species_recorded = uf$species, dbh_recorded_cm = uf$Dstem,
             key = ut$tree_id[match(paste(uf$site, uf$species), ut$site_sp)], measured = TRUE)
  ua <- .rd("data/processed/flux/untagged_rescue/untagged_auxfile.csv") %>% distinct(Sample, site, species, Dstem)
  M$ua <- .m(MON, "data/processed/flux/untagged_rescue/untagged_auxfile.csv", "Sample", ua$Sample,
             species_recorded = ua$species, dbh_recorded_cm = ua$Dstem,
             key = ut$tree_id[match(paste(ua$site, ua$species), ut$site_sp)], measured = TRUE)

  mon_rec <- mon_att %>% left_join(ym_att, by = "key") %>%
    left_join(tfd %>% transmute(key, notes_dead = sprintf("monthly field notes: dead (%d of %d visits)", n_dead, n)), by = "key") %>%
    mutate(campaign = MON, kind = ifelse(.is_tag(key), "census_tag", "untag"),
           dead_evidence = paste_na(ifelse(grepl("_dead$", key), "monthly label: _dead", NA), notes_dead),
           dbh_date = as.Date(NA), dbh_source = "monthly survey tree list (YM_trees_measured)",
           coord_source = "GPS (monthly survey tree list)") %>% select(-notes_dead, -dead_flag)

  # ---- 2021 multi-height and molecular survey -------------------------------------------
  h <- .rd("data/processed/flux/tree_flux_2021_multiheight.csv")
  h_att <- h %>% group_by(tag_recorded = tree_id) %>%
    summarise(species_code = fix_species_code(names(sort(table(species), decreasing = TRUE))[1]),
              landscape_position = names(sort(table(plot), decreasing = TRUE))[1],
              date = min(as.Date(substr(UniqueID, 1, 8), "%Y%m%d")), .groups = "drop")
  M$h <- .m(F21, "data/processed/flux/tree_flux_2021_multiheight.csv", "tree_id", h_att$tag_recorded, measured = TRUE)
  ax <- .rd("data/processed/flux/goflux_auxfile.csv")
  M$ax <- .m(F21, "data/processed/flux/goflux_auxfile.csv", "tree_id", unique(ax$tree_id), measured = TRUE)
  fdw <- .rd("data/processed/flux/flux_data_with_dimensions_weather_and_final_airtemp.csv")
  M$fdw <- .m(F21, "data/processed/flux/flux_data_with_dimensions_weather_and_final_airtemp.csv", "tree_id",
              unique(fdw$tree_id), measured = TRUE)
  ds <- .rd("data/raw/field_data/static_chamber_field/Gewirtzman_datasheets_fixing_dates.csv", fileEncoding = "UTF-8-BOM")
  M$ds <- .m(F21, "data/raw/field_data/static_chamber_field/Gewirtzman_datasheets_fixing_dates.csv", "Tree_ID",
             unique(ds$Tree_ID), measured = TRUE)
  dbh <- .rd("data/raw/dbh/Trees_DBH.csv", fileEncoding = "UTF-8-BOM")
  M$dbh <- .m(S21, "data/raw/dbh/Trees_DBH.csv", "Tree ID", unique(dbh$`Tree ID`), measured = FALSE)
  cs <- .rd("data/processed/tree_data/tree_dbh_consensus_comprehensive.csv")
  M$cs <- .m(S21, "data/processed/tree_data/tree_dbh_consensus_comprehensive.csv", "consensus_tree_id",
             unique(cs$consensus_tree_id), measured = FALSE)
  sm <- .rd("data/raw/field_data/static_chamber_field/soil_tree_level_means.csv")
  M$sm <- .m(S21, "data/raw/field_data/static_chamber_field/soil_tree_level_means.csv", "Tree_ID",
             unique(sm$Tree_ID), measured = FALSE)
  g <- .rd("data/processed/internal_gas/sample_data_only.csv")
  M$g <- .m(S21, "data/processed/internal_gas/sample_data_only.csv", "Tree.ID", unique(g$Tree.ID), measured = TRUE)
  c1 <- .rd("data/raw/tree_cores/Tree_Core_Sectioning_Data.csv")
  M$c1 <- .m(S21, "data/raw/tree_cores/Tree_Core_Sectioning_Data.csv", "seq_id", unique(c1$seq_id), measured = TRUE)
  c2 <- .rd("data/raw/tree_cores/Tree Core Sectioning Data.csv")
  M$c2 <- .m(S21, "data/raw/tree_cores/Tree Core Sectioning Data.csv", "Tree Core ID", unique(c2$`Tree Core ID`), measured = TRUE)
  cf <- .rd("data/processed/tree_cores/tree_core_filled_complete.csv")
  M$cf <- .m(S21, "data/processed/tree_cores/tree_core_filled_complete.csv", "tree_id", unique(cf$tree_id), measured = TRUE)
  dm <- .rd("data/raw/ddpcr/ddPCR_meta_all_data.csv")
  M$dm <- .m(S21, "data/raw/ddpcr/ddPCR_meta_all_data.csv", "seq_id", unique(dm$seq_id), measured = TRUE)
  for (f in c("data/raw/ddpcr/ddPCR_meta_all_data.csv", "data/processed/molecular/processed_ddpcr_data.csv",
              "data/raw/ddpcr_raw/ddPCR_tree_full_data.csv")) {
    s <- unique(.first_line(.rd(f)$Inner.Core.Sample.ID))
    M[[f]] <- .m(S21, f, "Inner.Core.Sample.ID (first line)", s, key = ncbi_tree_id(sub("^[A-Z]{4} ?_", "", s)), measured = TRUE)
  }
  s16 <- read.delim("data/raw/16s/tree_16s_mapping_dada2_corrected.txt", stringsAsFactors = FALSE)
  s16 <- unique(s16$seq_id)
  M$s16 <- .m(S21, "data/raw/16s/tree_16s_mapping_dada2_corrected.txt", "seq_id", s16,
              key = ifelse(s16 == "MISC", NA, ncbi_tree_id(s16)), measured = TRUE,
              note = ifelse(s16 == "MISC", "unassigned libraries (blanks, mock, muck and water samples)",
                     ifelse(s16 == "RO7/8", "soil sampled between RO7 and RO8", NA)))
  M$s16$key[M$s16$tag_recorded == "RO7/8"] <- NA

  pic <- unique(unlist(lapply(setdiff(list.files("data/raw/internal_gas/picarro", "_results.csv$", full.names = TRUE),
                                      list.files("data/raw/internal_gas/picarro", "merged", full.names = TRUE)),
                              function(f) utils::read.csv(f, check.names = FALSE)$SampleName)))
  pic <- pic[!grepl("^Amb|^Atm|room|^V[0-9]|FLUSH|RINS", pic, ignore.case = TRUE) &
               !pic %in% c("SB1","SB3a","SB3b","SB3","SB4a","SB4","SB5a","SB5","S3a","S3b","S3c","SA1")]  # SA1: a standard here
  hs <- grepl("[HS]$", pic)                                # heartwood/sapwood pair (2024)
  M$pic <- .m(ifelse(hs, G24, S21), "data/raw/internal_gas/picarro/*_results.csv", "SampleName", pic,
              key = ncbi_tree_id(sub("_Redo$", "", ifelse(hs, sub("[HS]$", "", pic), pic))), measured = TRUE)

  # ---- 2024 heartwood/sapwood gas -----------------------------------------------------------
  g4 <- .rd("data/processed/internal_gas/stem_gas_2024_calibrated.csv") %>% filter(Sample.Type == "Sample", !is.na(Tree.No))
  M$g4 <- .m(G24, "data/processed/internal_gas/stem_gas_2024_calibrated.csv", "Tree.No", unique(g4$Tree.No), measured = TRUE)

  # ---- 2023 cross-species survey ------------------------------------------------------------
  y0 <- .rd("data/raw/field_data/static_chamber_field/Compiled YMF Data 2023 - Sheet1.csv")
  # a row without a flux start time was not measured (one dead birch "too damaged to measure"; tag 131)
  y0$timed <- !(is.na(y0$`Flux Start (System) Time`) | trimws(y0$`Flux Start (System) Time`) %in% c("", "NA"))
  y_rec <- y0 %>% transmute(tag_recorded = trimws(`Tree Tag`), species_recorded = `Species Code`, timed,
                            dbh_recorded_cm = .num(`DBH (cm)`), date = as.Date(Date, "%m/%d/%y"),
                            dead_evidence = paste_na(
                              ifelse(grepl("^dead$", trimws(`Bark Missing`), ignore.case = TRUE) |
                                     grepl("^dead$", trimws(`Wounding/Holes`), ignore.case = TRUE), "2023 condition score: dead", NA),
                              ifelse(grepl("dead|snag", Notes, ignore.case = TRUE), "2023 notes: dead", NA)))
  M$y0 <- .m(C23, "data/raw/field_data/static_chamber_field/Compiled YMF Data 2023 - Sheet1.csv", "Tree Tag",
             y_rec$tag_recorded, species_recorded = y_rec$species_recorded, dbh_recorded_cm = y_rec$dbh_recorded_cm,
             measured = y_rec$timed, note = ifelse(y_rec$timed, NA, "2023 row without a flux start time: not measured"))
  for (f in c("data/processed/flux/tree_flux_2023_cross_species.csv", "data/processed/flux/ymf2023_goflux_auxfile.csv")) {
    y <- .rd(f)
    M[[f]] <- .m(C23, f, "Tree Tag", trimws(y$`Tree Tag`), species_recorded = y$`Species Code`,
                 dbh_recorded_cm = .num(y$`DBH (cm)`), measured = TRUE)
  }

  # ---- felled black oak ---------------------------------------------------------------------
  M$oak <- .m(OAK, "data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv", "(whole file: one tree)", "felled black oak", key = FELLED_OAK_ID,
              kind = "felled", measured = TRUE)

  A <- bind_rows(M)
  for (v in c("key", "kind", "note", "species_recorded")) if (!v %in% names(A)) A[[v]] <- NA
  A$dbh_recorded_cm <- .num(A$dbh_recorded_cm)
  A$species_recorded <- as.character(A$species_recorded)

  # ---- classify and key every spelling ------------------------------------------------------
  t0 <- A$tag_recorded
  is_mon <- A$campaign == MON; is_23 <- A$campaign == C23
  ulab <- function(z) gsub("[^A-Za-z0-9]", "", z)
  A$kind <- dplyr::case_when(
    !is.na(A$kind) ~ A$kind,
    .is_tag(t0) ~ "census_tag",
    is_mon ~ "untag",
    is_23 & grepl("^Ash[0-9]+$", t0) ~ "ash2023",
    is_23 & grepl("^untagged", t0, ignore.case = TRUE) ~ "untag2023",
    TRUE ~ "name")
  need <- is.na(A$key) & is.na(A$note)
  A$key[need & A$kind == "census_tag"] <- sub("\\.0$", "", trimws(t0[need & A$kind == "census_tag"]))
  A$key[need & A$kind == "untag"] <- ut$tree_id[match(ulab(t0[need & A$kind == "untag"]), ulab(ut$tree_id))]
  A$key[need & A$kind == "ash2023"] <- t0[need & A$kind == "ash2023"]
  A$key[need & A$kind == "name"] <- ncbi_tree_id(t0[need & A$kind == "name"])
  A$key[A$tag_recorded %in% c("Untagged", "untagged") & A$campaign == S21] <- NA
  A$note[A$tag_recorded %in% c("Untagged", "untagged") & A$campaign == S21] <- "untagged trees in the 2021 DBH list"
  A$key[A$kind == "untag2023"] <- NA                       # resolved per record below
  A$fix_species <- fix_species_code(A$species_recorded)

  list(mentions = A, monthly = mon_rec, flux21 = h_att, y2023 = y_rec[y_rec$timed, ], gas24 = g4,
       consensus = cs, dbh21 = dbh, datasheets = ds, gas21 = g, ddpcr = dm[, c("seq_id", "species")], cores = c1[, c("seq_id", "Tree Core Species Code")], untagged = ut)
}

paste_na <- function(...) {
  x <- do.call(cbind, lapply(list(...), as.character))
  apply(x, 1, function(r) { r <- r[!is.na(r) & r != ""]; if (length(r)) paste(r, collapse = "; ") else NA_character_ })
}
