#!/usr/bin/env Rscript
# ==============================================================================
# 04_trees.R -- inventory_stems.csv, trees.csv and tree_tags.csv (data/tables/)
# ------------------------------------------------------------------------------
# inventory_stems  one row per census stem (outputs/tables/inventory_stems.csv,
#                  from 08_inventory_build.R), with stem_id = <source>-<tag>-<stem>
#                  (stem = the census stem label, else 1, 2, ... in file order) and
#                  the physical tree it belongs to (tree_id).
# trees            one row per physical tree measured for flux, internal gas, wood
#                  cores or DNA in any campaign.
# tree_tags        every spelling of every measured tree in the files the pipeline
#                  reads -> tree_id. species_recorded / dbh_recorded_cm are filled
#                  only where the spelling alone names more than one tree in that
#                  column (e.g. 2023 tag 1985 was written on two trees).
#
# TREE IDs
#   census tag      the plot's stem tag (the 2019 plot census and the 2018 by-tag
#                   survey use one tag series; the monthly and 2023 surveys use the same
#                   tags). Where one tag belongs to more than one physical tree, each
#                   gets <tag>_a, <tag>_b, ... (census trees first: 2019 census, then
#                   by-tag survey, in file order; then measured trees not in the census).
#   2021 survey     the NCBI spelling of the field name (code/lib/tree_ids.R).
#   monthly untagged  the survey label (UNTAG_<position>_<species>[_dead]).
#   2023 untagged   UNTAG_2023_<species>[_<n>] (n by date, then DBH); the forked
#                   one UNTAG_2023_BELE_forked. 2023 ash: as recorded (Ash1-Ash30).
#   felled oak      QUVE_felled_2022.
#
# ONE TAG, SEVERAL PLACES. Census stems sharing a tag are one tree when they have
# the same species and lie within LINK_M of each other (multi-stem trees; the
# 2018 multi-stem sheet puts every stem at the tree's point). A measured tree
# matches a census tree with the same tag when the species agree and the DBHs are
# within a factor of 1.5, or the species disagree but the DBHs agree to
# max(2.5 cm, 10%). Proposals and evidence are written to
# outputs/audit/tag_collisions.csv; the decisions applied are those in
# code/lib/tag_collision_decisions.csv (hand-editable; every record of every
# collision must be listed there, or the script stops).
#
# DEAD. One definition, applied here; the per-deployment flag in flux_stem.csv
# (dead_stem) is read from this table. Dated marks: the monthly label (_dead; every
# visit of an untagged snag), monthly field notes (dead|snag, per visit), a 2021 name
# containing "dead" (per deployment), the 2023 condition score "dead" (bark missing or
# wounding) or 2023 notes (dead|snag). A tree is dead from the first date any source
# marks it dead, but only if it is marked dead on more than one date, or is not
# measured again after that date; measurements before that date are of a live stem.
# A tree marked dead once and measured again later without the mark stays live.
#   dead          TRUE / FALSE by that rule
#   dead_from     first date marked dead (dead trees)
#   dead_evidence every mark with its date; dead_basis: how the rule decided.
#
# DBH. One rule for the 2021 survey: the 2021 flux datasheets write a stem diameter
# for each chamber height (50, 125, 200 cm; decreasing up the stem). DBH is the one
# written for the 125 cm chamber (nearest breast height; the median of the three);
# trees without a 125 cm row (Kalmia, 75 cm and root crown) take the median of their
# rows; trees with no datasheet diameter take the 2021 tree DBH list (Trees_DBH.csv,
# via the consensus file). Monthly survey: its tree list; 2023: the field sheet;
# census-only: the census stem. Measured DBH only: a tree with no recorded DBH anywhere
# has dbh_cm NA and dbh_source "none: no DBH measurement recorded" (nothing is imputed;
# the trees are listed in the build report).
#
# SPECIES. species_source records where species_code comes from: a field record
# (campaign), the plot census, the 2021 records (majority of the sources listed)
# or the tree-name prefix. Trees whose 2021 species disagrees with their name prefix
# (e.g. RM105, recorded QURU) are listed in tree_species_conflicts.csv.
#
# Writes data/tables/{inventory_stems,trees,tree_tags}.csv
#        outputs/audit/{tag_collisions,tree_species_conflicts,tree_dead_status,
#                       tree_tags_unresolved}.csv, outputs/audit/trees_build.txt
# ==============================================================================
source("code/lib/outputs.R")
suppressPackageStartupMessages(library(tidyverse))
source("code/lib/geometry.R"); source("code/lib/tree_mentions.R")
OUT <- "data/tables"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DECISIONS <- "code/lib/tag_collision_decisions.csv"
LINK_M <- 3                       # census stems of one tree: same species, within 3 m

first_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else x[NA_integer_][1] }
paste_u <- function(x) { x <- unique(na.omit(unlist(strsplit(x, "; ")))); if (length(x)) paste(x, collapse = "; ") else NA_character_ }

R <- tree_mentions()
A <- R$mentions %>% filter(!is.na(tag_recorded), trimws(tag_recorded) != "")
tr <- geo_transforms()

# ==============================================================================
# 1. INVENTORY STEMS
# ==============================================================================
I <- read.csv("outputs/tables/inventory_stems.csv", stringsAsFactors = FALSE)
I$row <- seq_len(nrow(I))
stem_k <- function(stem) {
  k <- as.character(stem); used <- k[!is.na(k)]; i <- 0L
  for (j in which(is.na(k))) { repeat { i <- i + 1L; if (!as.character(i) %in% used) break }; k[j] <- as.character(i) }
  k
}
I <- I %>% group_by(source, tag) %>% mutate(stem_k = stem_k(stem)) %>% ungroup() %>%
  mutate(stem_id = paste(source, ifelse(is.na(tag), "NA", tag), stem_k, sep = "-"),
         species_code = fix_species_code(species_code))
stopifnot(!anyDuplicated(I$stem_id))
ll <- tr$fwd(ifelse(I$located, I$PX, NA), ifelse(I$located, I$PY, NA))
I$latitude <- ll$lat; I$longitude <- ll$lon

# ==============================================================================
# 2. MEASURED TREES CARRYING A CENSUS TAG
# ==============================================================================
mon <- R$monthly
# 2021 flux datasheets: the diameter written for the 125 cm chamber, else the median
DS_DBH <- R$datasheets %>% mutate(d = .num(DBH), h = trimws(`Measurement height`)) %>% filter(is.finite(d)) %>%
  group_by(Tree_ID) %>%
  summarise(d = if (any(h == "125")) d[h == "125"][1] else median(d),
            dbh_source = if (any(h == "125")) "2021 flux datasheet: stem diameter at the 125 cm chamber"
                         else "2021 flux datasheet: median of the diameters at its chamber heights", .groups = "drop")
dbh21 <- DS_DBH
MR <- bind_rows(
  mon %>% filter(kind == "census_tag") %>%
    transmute(campaign, tag = key, species_code, dbh_cm, date, landscape_position, gps_latitude, gps_longitude,
              dead_evidence, dbh_date, dbh_source, coord_source),
  R$flux21 %>% filter(.is_tag(tag_recorded)) %>%
    transmute(campaign = F21, tag = tag_recorded, species_code, dbh_cm = dbh21$d[match(tag_recorded, dbh21$Tree_ID)],
              date, landscape_position, dbh_source = dbh21$dbh_source[match(tag_recorded, dbh21$Tree_ID)]),
  R$y2023 %>% filter(.is_tag(tag_recorded)) %>%
    transmute(campaign = C23, tag = tag_recorded, species_code = fix_species_code(species_recorded),
              dbh_cm = dbh_recorded_cm, date, dead_evidence, dbh_date = date, dbh_source = "2023 field sheet") %>%
    group_by(campaign, tag, species_code, dbh_cm) %>%
    summarise(date = min(date), dbh_date = min(dbh_date), dbh_source = first(dbh_source),
              dead_evidence = paste_u(dead_evidence), .groups = "drop")) %>%
  mutate(tag = as.integer(tag), record_key = paste(campaign, tag, species_code, dbh_cm, sep = "|"))
stopifnot(!anyDuplicated(MR$record_key))

# ---- census clusters ---------------------------------------------------------------
cluster_stems <- function(s) {                 # s: census stems sharing one tag
  n <- nrow(s); cl <- seq_len(n)
  if (n > 1) {
    find <- function(i) { while (cl[i] != i) i <- cl[i]; i }
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      same_sp <- !is.na(s$species_code[i]) && identical(s$species_code[i], s$species_code[j])
      d <- sqrt((s$PX[i] - s$PX[j])^2 + (s$PY[i] - s$PY[j])^2)
      if (same_sp && s$located[i] && s$located[j] && is.finite(d) && d <= LINK_M) cl[find(j)] <- find(i)
    }
    cl <- sapply(seq_len(n), find)
    for (i in which(!s$located)) {               # unlocated: join the only same-species tree
      o <- unique(cl[s$located & s$species_code %in% s$species_code[i]])
      if (length(o) == 1) cl[i] <- o
    }
  }
  match(cl, unique(cl))
}
compat <- function(sp1, d1, sp2, d2) {          # measured record vs census stem / another record
  same <- !is.na(sp1) & !is.na(sp2) & sp1 == sp2
  r <- pmax(d1, d2) / pmin(d1, d2); dd <- abs(d1 - d2)
  ok_d <- is.finite(r)
  ifelse(same & (!ok_d | r <= 1.5), ifelse(!ok_d, "species agree (no DBH)", ifelse(r <= 1.25, "agree", "species agree, DBH differs")),
  ifelse(!same & ok_d & dd <= pmax(2.5, 0.1 * pmax(d1, d2)), "DBH agrees, species differs", NA))
}

tags <- sort(unique(c(I$tag[!is.na(I$tag)], MR$tag)))
multi <- tags[tags %in% I$tag[duplicated(paste(I$tag))] | tags %in% MR$tag]   # anything beyond one plain stem
assign_rows <- list()
for (t in multi) {
  s <- I %>% filter(tag %in% t) %>% arrange(source != "fg19", row)
  m <- MR %>% filter(tag == t) %>% arrange(match(campaign, c(MON, F21, C23)))
  cl <- if (nrow(s)) cluster_stems(s) else integer(0)
  # each census tree: species and DBH of its largest stem
  ct <- if (nrow(s)) s %>% mutate(cl = cl) %>% group_by(cl) %>% slice_max(dbh_cm, n = 1, with_ties = FALSE) %>% ungroup() else NULL
  # best census tree for each measured record, best pairs first; one record per campaign per tree
  m$cl <- NA_integer_; m$how <- NA_character_
  if (!is.null(ct) && nrow(m)) {
    P <- expand.grid(i = seq_len(nrow(m)), j = seq_len(nrow(ct)))
    P$h <- compat(m$species_code[P$i], m$dbh_cm[P$i], ct$species_code[P$j], ct$dbh_cm[P$j])
    P <- P[!is.na(P$h), ]
    P$sc <- match(P$h, c("agree", "species agree (no DBH)", "species agree, DBH differs", "DBH agrees, species differs")) * 1000 +
      coalesce(abs(m$dbh_cm[P$i] - ct$dbh_cm[P$j]), 0)
    for (k in order(P$sc)) {
      i <- P$i[k]; c0 <- ct$cl[P$j[k]]
      if (!is.na(m$cl[i]) || any(m$cl %in% c0 & m$campaign == m$campaign[i])) next
      m$cl[i] <- c0; m$how[i] <- P$h[k]
    }
  }
  nxt <- if (length(cl)) max(cl) else 0L
  for (i in which(is.na(m$cl))) {              # measured trees not in the census: same tree across campaigns?
    prev <- which(!is.na(m$cl) & m$cl > (if (length(cl)) max(cl) else 0L) & m$campaign != m$campaign[i])
    h <- if (length(prev)) compat(m$species_code[i], m$dbh_cm[i], m$species_code[prev], m$dbh_cm[prev]) else NA
    if (length(prev) && any(!is.na(h))) { b <- prev[which(!is.na(h))[1]]; m$cl[i] <- m$cl[b]; m$how[i] <- paste(h[which(!is.na(h))[1]], "(other campaign)") }
    else { nxt <- nxt + 1L; m$cl[i] <- nxt; m$how[i] <- "not in the census" }
  }
  n_tree <- length(unique(c(cl, m$cl)))
  ids <- if (n_tree == 1) as.character(t) else paste0(t, "_", letters[seq_len(n_tree)])
  key_of <- sort(unique(c(cl, m$cl)))
  # distance of each stem to its tree's largest stem, and to the nearest stem of another tree
  dist_own <- dist_other <- rep(NA_real_, nrow(s))
  for (i in seq_len(nrow(s))) {
    d <- sqrt((s$PX[i] - s$PX)^2 + (s$PY[i] - s$PY)^2)
    b <- ct$stem_id[ct$cl == cl[i]]; dist_own[i] <- d[s$stem_id == b]
    o <- d[cl != cl[i]]; if (length(o) && any(is.finite(o))) dist_other[i] <- min(o, na.rm = TRUE)
  }
  trivial <- n_tree == 1 && all(m$how %in% c("agree", "species agree (no DBH)", "not in the census", "agree (other campaign)")) &&
             (nrow(s) <= 1 || (all(s$source == "bytag") && all(dist_own %in% 0)))
  assign_rows[[length(assign_rows) + 1]] <- bind_rows(
    if (nrow(s)) tibble(tag = t, record_key = s$stem_id, record = paste("census", s$source), stem_id = s$stem_id,
                        species_code = s$species_code, dbh_cm = s$dbh_cm, PX = s$PX, PY = s$PY,
                        dist_to_tree_m = round(dist_own, 2), dist_to_other_tree_m = round(dist_other, 2),
                        proposed_tree_id = ids[match(cl, key_of)],
                        how = ifelse(dist_own %in% 0 & s$source == "bytag" & grepl("multi-stem", s$position),
                                     "multi-stem sheet", "census stem")),
    if (nrow(m)) tibble(tag = t, record_key = m$record_key, record = m$campaign, species_code = m$species_code,
                        dbh_cm = m$dbh_cm, proposed_tree_id = ids[match(m$cl, key_of)], how = m$how)) %>%
    mutate(n_trees = n_tree, collision = !trivial)
}
TA <- bind_rows(assign_rows)

# unsure: a judgement worth checking by hand
TA$unsure <- FALSE
for (i in seq_len(nrow(TA))) {
  g <- TA[TA$tag == TA$tag[i], ]; cen <- grepl("^census", g$record)
  u <- grepl("DBH agrees, species differs|species agree, DBH differs", TA$how[i])
  if (grepl("^census", TA$record[i]))
    u <- u || n_distinct(g$record[cen & g$proposed_tree_id == TA$proposed_tree_id[i]]) > 1 || (is.finite(TA$dist_to_tree_m[i]) && TA$dist_to_tree_m[i] > 1) ||
      (is.finite(TA$dist_to_other_tree_m[i]) && TA$dist_to_other_tree_m[i] <= 5 &&
         any(g$species_code[cen & g$proposed_tree_id != TA$proposed_tree_id[i]] %in% TA$species_code[i]))
  if (TA$how[i] == "not in the census" && any(cen))
    u <- u || any(abs(g$dbh_cm[cen] - TA$dbh_cm[i]) <= 0.3 * TA$dbh_cm[i], na.rm = TRUE) ||
      any(g$species_code[cen] %in% TA$species_code[i])
  TA$unsure[i] <- u
}

# ---- decisions ------------------------------------------------------------------------
COL <- TA %>% filter(collision)
COLS <- c("tag", "record_key", "record", "stem_id", "species_code", "dbh_cm", "PX", "PY", "dist_to_tree_m",
          "dist_to_other_tree_m", "how", "n_trees", "unsure", "proposed_tree_id")
write_csv(COL[, COLS], out_path("tag_collisions.csv", "audit"), na = "")
if (!file.exists(DECISIONS)) stop("missing ", DECISIONS, " (proposals: outputs/audit/tag_collisions.csv)")
DEC <- read.csv(DECISIONS, stringsAsFactors = FALSE, colClasses = "character")
miss <- setdiff(COL$record_key, DEC$record_key)
if (length(miss)) {
  write_csv(COL %>% filter(record_key %in% miss), out_path("tag_collisions_undecided.csv", "audit"))
  stop(length(miss), " collision records have no decision in ", DECISIONS,
       " (listed in outputs/audit/tag_collisions_undecided.csv)")
}
TA$tree_id <- TA$proposed_tree_id
hit <- match(TA$record_key, DEC$record_key)
TA$tree_id[!is.na(hit)] <- DEC$tree_id[hit[!is.na(hit)]]
TA$decision_status <- DEC$status[hit]
stale <- setdiff(DEC$record_key, TA$record_key)
write_csv(TA %>% filter(collision) %>% select(all_of(COLS), tree_id, decision_status),
          out_path("tag_collisions.csv", "audit"), na = "")

# tree_id of every census stem
I$tree_id <- ifelse(is.na(I$tag), paste0("NOTAG_", I$source, "_", cumsum(is.na(I$tag))), as.character(I$tag))
h <- match(I$stem_id, TA$record_key); I$tree_id[!is.na(h)] <- TA$tree_id[h[!is.na(h)]]
MR$tree_id <- as.character(MR$tag)
h <- match(MR$record_key, TA$record_key); MR$tree_id[!is.na(h)] <- TA$tree_id[h[!is.na(h)]]

# ==============================================================================
# 3. TREES
# ==============================================================================
# ---- census-tag trees ------------------------------------------------------------------
cen <- I %>% filter(tree_id %in% MR$tree_id) %>% group_by(tree_id) %>%
  arrange(desc(located), desc(dbh_cm), .by_group = TRUE) %>%
  summarise(stem_id = first(stem_id), c_lat = first(latitude), c_lon = first(longitude),
            c_src = first(source), c_sp = first(species_code), c_dbh = first(dbh_cm), .groups = "drop")
T1 <- MR %>% arrange(match(campaign, c(MON, F21, C23))) %>% group_by(tree_id) %>%
  summarise(sp_src = campaign[!is.na(species_code)][1], species_code = first_non_na(species_code),
            dbh_date = dbh_date[!is.na(dbh_cm)][1], dbh_source = dbh_source[!is.na(dbh_cm)][1],
            dbh_cm = first_non_na(dbh_cm),
            landscape_position = first_non_na(landscape_position),
            gps_latitude = first_non_na(gps_latitude), gps_longitude = first_non_na(gps_longitude),
            dead_evidence = paste_u(dead_evidence), campaigns_m = paste(unique(campaign), collapse = ";"),
            sp_all = paste(unique(na.omit(species_code)), collapse = "/"), .groups = "drop") %>%
  left_join(cen, by = "tree_id") %>%
  mutate(species_source = ifelse(!is.na(species_code), paste("field record:", sp_src),
                                 ifelse(c_src == "fg19", "plot census 2019", "plot by-tag survey 2018")),
         species_code = coalesce(species_code, c_sp),
         dbh_source = ifelse(is.na(dbh_cm) & !is.na(c_dbh), ifelse(c_src == "fg19", "plot census 2019", "plot by-tag survey 2018"), dbh_source),
         dbh_cm = coalesce(dbh_cm, c_dbh),
         latitude = coalesce(c_lat, gps_latitude), longitude = coalesce(c_lon, gps_longitude),
         coord_source = ifelse(!is.na(c_lat), ifelse(c_src == "fg19", "plot census 2019 (PX, PY)", "plot by-tag survey 2018 (PX, PY)"),
                        ifelse(!is.na(gps_latitude), "GPS (monthly survey tree list)", NA)),
         in_inventory = stem_id, census_species = c_sp)

# ---- named trees: 2021 survey and 2024 resample ------------------------------------------
NM <- A %>% filter(kind == "name", !is.na(key))
named <- NM %>% group_by(key) %>% filter(any(measured)) %>% ungroup()
votes <- bind_rows(
  R$flux21 %>% transmute(key = ncbi_tree_id(tag_recorded), sp = species_code, src = "2021 flux"),
  R$gas21 %>% transmute(key = ncbi_tree_id(Tree.ID), sp = fix_species_code(Species.ID), src = "2021 gas"),
  R$ddpcr %>% transmute(key = ncbi_tree_id(seq_id), sp = fix_species_code(species), src = "ddPCR"),
  R$cores %>% transmute(key = ncbi_tree_id(seq_id), sp = fix_species_code(`Tree Core Species Code`), src = "cores"),
  R$dbh21 %>% transmute(key = ncbi_tree_id(`Tree ID`), sp = fix_species_code(Species), src = "2021 DBH list"),
  R$datasheets %>% transmute(key = ncbi_tree_id(Tree_ID), sp = fix_species_code(Tree_sp), src = "2021 datasheets")) %>%
  filter(!is.na(key), !is.na(sp)) %>% distinct()
sp21 <- votes %>% group_by(key, sp) %>% summarise(n = n(), srcs = paste(sort(unique(src)), collapse = ", "), .groups = "drop") %>%
  arrange(key, desc(n)) %>% group_by(key) %>%
  summarise(species_code = first(sp), species_source = paste0("2021 records: ", first(srcs)),
            sp_all = paste(unique(sp), collapse = "/"), .groups = "drop")
plot21 <- bind_rows(R$flux21 %>% transmute(key = ncbi_tree_id(tag_recorded), p = landscape_position, o = 1),
                    R$datasheets %>% transmute(key = ncbi_tree_id(Tree_ID), p = Plot, o = 2),
                    R$consensus %>% transmute(key = ncbi_tree_id(consensus_tree_id), p = Plot, o = 3),
                    R$dbh21 %>% transmute(key = ncbi_tree_id(`Tree ID`), p = Plot, o = 4)) %>%
  filter(!is.na(p), p != "") %>% arrange(key, o) %>% group_by(key) %>% summarise(landscape_position = first(p), .groups = "drop")
dbhc <- bind_rows(
  DS_DBH %>% transmute(key = ncbi_tree_id(Tree_ID), dbh_cm = d, dbh_source),
  R$consensus %>% filter(source != "datasheets") %>%
    transmute(key = ncbi_tree_id(consensus_tree_id), dbh_cm = .num(DBH), dbh_source = paste0("2021 tree DBH list (", source, ")"))) %>%
  filter(is.finite(dbh_cm)) %>% distinct(key, .keep_all = TRUE)
camp21 <- named %>% filter(measured) %>% group_by(key) %>% summarise(campaigns = paste(sort(unique(campaign)), collapse = ";"), .groups = "drop")
T2 <- tibble(tree_id = unique(named$key)) %>%
  left_join(sp21, by = c("tree_id" = "key")) %>% left_join(plot21, by = c("tree_id" = "key")) %>%
  left_join(dbhc, by = c("tree_id" = "key")) %>% left_join(camp21, by = c("tree_id" = "key")) %>%
  mutate(species_source = ifelse(is.na(species_code), "tree-name prefix", species_source),
         species_code = coalesce(species_code, prefix_species(tree_id)),
         dead_evidence = ifelse(tree_id %in% named$key[grepl("dead", named$tag_recorded, ignore.case = TRUE)], "2021 name: dead", NA))

# ---- monthly untagged, 2023 untagged and ash, felled oak -------------------------------------
T3 <- mon %>% filter(kind == "untag") %>%
  transmute(tree_id = key, species_code, species_source = "monthly survey (tree list; untagged: field notes)", dbh_cm, dbh_date, dbh_source, landscape_position,
            latitude = gps_latitude, longitude = gps_longitude, coord_source, dead_evidence, campaigns = MON)
y <- R$y2023 %>% mutate(species_code = fix_species_code(species_recorded))
yu <- y %>% filter(grepl("^untagged", tag_recorded, ignore.case = TRUE)) %>% distinct(tag_recorded, species_code, dbh_recorded_cm, .keep_all = TRUE) %>%
  arrange(species_code, date, dbh_recorded_cm) %>% group_by(species_code) %>%
  mutate(n_sp = sum(!grepl("forked", tag_recorded)),
         tree_id = case_when(grepl("forked", tag_recorded) ~ paste0("UNTAG_2023_", species_code, "_forked"),
                             n_sp > 1 ~ paste0("UNTAG_2023_", species_code, "_", cumsum(!grepl("forked", tag_recorded))),
                             TRUE ~ paste0("UNTAG_2023_", species_code))) %>% ungroup()
ya <- y %>% filter(grepl("^Ash[0-9]+$", tag_recorded)) %>% mutate(tree_id = tag_recorded)
T4 <- bind_rows(yu, ya) %>% group_by(tree_id) %>%
  summarise(species_code = first(species_code), dbh_cm = first(dbh_recorded_cm), dbh_date = first(date),
            dead_evidence = paste_u(dead_evidence), .groups = "drop") %>%
  mutate(dbh_source = "2023 field sheet", species_source = "field record: 2023_cross_species", campaigns = C23)
oak <- read.csv("data/raw/field_data/black_oak/ymf_black_oak_flux_compiled.csv", check.names = FALSE)
T5 <- tibble(tree_id = FELLED_OAK_ID, species_code = "QUVE", species_source = "felled-oak sheet",
             dbh_cm = .num(oak$Stem_Diam_mm[oak$Height_m == "1.25"])[1], dbh_date = as.Date("2022-10-04"),
             dbh_source = "felled-oak sheet, stem diameter at 1.25 m (column labelled mm; values are cm)",
             campaigns = OAK)

# ==============================================================================
# 4. EVERY SPELLING -> tree_id
# ==============================================================================
A$tree_id <- NA_character_
ct <- A$kind == "census_tag" & !is.na(A$key)
for (i in which(ct)) {
  cmp <- if (A$campaign[i] == S21) F21 else A$campaign[i]
  c <- MR %>% filter(tag == as.integer(A$key[i]), campaign == cmp)
  if (!nrow(c)) next
  if (n_distinct(c$tree_id) == 1) { A$tree_id[i] <- c$tree_id[1]; next }
  sp <- fix_species_code(A$species_recorded[i]); d <- A$dbh_recorded_cm[i]
  sc <- ifelse(c$species_code %in% sp, 0, 1000) + ifelse(is.finite(d), abs(c$dbh_cm - d), 500)
  if (!is.na(sp) || is.finite(d)) A$tree_id[i] <- c$tree_id[which.min(sc)]
}
A$tree_id[A$kind == "untag"] <- A$key[A$kind == "untag"]
A$tree_id[A$kind %in% c("ash2023", "felled")] <- A$key[A$kind %in% c("ash2023", "felled")]
A$tree_id[A$kind == "name" & A$key %in% T2$tree_id] <- A$key[A$kind == "name" & A$key %in% T2$tree_id]
for (i in which(A$kind == "untag2023")) {
  c <- yu %>% filter(species_code %in% fix_species_code(A$species_recorded[i]))
  if (grepl("forked", A$tag_recorded[i])) c <- c %>% filter(grepl("forked", tag_recorded))
  if (nrow(c)) { b <- which.min(abs(c$dbh_recorded_cm - A$dbh_recorded_cm[i])); if (abs(c$dbh_recorded_cm[b] - A$dbh_recorded_cm[i]) < 1) A$tree_id[i] <- c$tree_id[b] }
}

# campaigns of census-tag trees from every measured spelling
camp_ct <- A %>% filter(measured, !is.na(tree_id)) %>% group_by(tree_id) %>%
  summarise(campaigns_a = paste(sort(unique(campaign)), collapse = ";"), .groups = "drop")

# ==============================================================================
# 5. ASSEMBLE
# ==============================================================================
TREES <- bind_rows(
  T1 %>% left_join(camp_ct, by = "tree_id") %>% mutate(campaigns = coalesce(campaigns_a, campaigns_m)),
  T2 %>% mutate(coord_source = NA_character_),
  T3, T4, T5) %>%
  mutate(campaigns = sapply(strsplit(campaigns, ";"), function(z)
           paste(intersect(c(MON, F21, S21, OAK, C23, G24), z), collapse = ";")),
         dead = !is.na(dead_evidence)) %>%
  transmute(tree_id, species_code, species_source, dbh_cm = round(dbh_cm, 1), dbh_date = as.Date(dbh_date), dbh_source,
            landscape_position, latitude = round(latitude, 7), longitude = round(longitude, 7), coord_source,
            dead, dead_evidence, in_inventory, campaigns) %>% arrange(tree_id)
if (anyDuplicated(TREES$tree_id)) stop("duplicate tree_id: ", paste(unique(TREES$tree_id[duplicated(TREES$tree_id)]), collapse = ", "))

# ---- tree_tags ------------------------------------------------------------------------------------
TT <- A %>% filter(!is.na(tree_id), tree_id %in% TREES$tree_id) %>%
  distinct(campaign, source_file, source_column, tag_recorded, species_recorded, dbh_recorded_cm, tree_id) %>%
  group_by(campaign, source_file, source_column, tag_recorded) %>%
  mutate(ambiguous = n_distinct(tree_id) > 1) %>% ungroup() %>%
  mutate(species_recorded = ifelse(ambiguous, species_recorded, NA), dbh_recorded_cm = ifelse(ambiguous, dbh_recorded_cm, NA)) %>%
  distinct(campaign, source_file, source_column, tag_recorded, species_recorded, dbh_recorded_cm, tree_id) %>%
  arrange(source_file, source_column, tag_recorded)
k <- paste(TT$campaign, TT$source_file, TT$source_column, TT$tag_recorded, TT$species_recorded, TT$dbh_recorded_cm)
if (anyDuplicated(k)) { print(TT[k %in% k[duplicated(k)], ]); stop("tree_tags: a spelling maps to two trees with the same qualifiers") }
UNRES <- A %>% filter(is.na(tree_id) | !tree_id %in% TREES$tree_id) %>%
  distinct(campaign, source_file, source_column, tag_recorded, kind, key, measured, note) %>%
  mutate(reason = case_when(!is.na(note) ~ note,
                            kind == "name" & !measured ~ "tree list or soil means only: no flux, gas, core or DNA record",
                            kind == "name" ~ "name not measured", kind == "census_tag" ~ "tag with no measured record in this campaign",
                            TRUE ~ "unresolved"))

# ==============================================================================
# 5b. DEAD: dated marks and the rule (header)
# ==============================================================================
tt_tree <- function(file, column, tag, sp = NA, dbh = NA) {   # spelling -> tree_id (NA if none)
  t <- TT[TT$source_file == file & TT$source_column == column, ]
  tag <- trimws(as.character(tag)); sp <- rep_len(sp, length(tag)); dbh <- rep_len(dbh, length(tag))
  vapply(seq_along(tag), function(i) {
    c <- t[trimws(t$tag_recorded) %in% tag[i], ]
    if (nrow(c) > 1) c <- c[c$species_recorded %in% sp[i] & abs(.num(c$dbh_recorded_cm) - .num(dbh[i])) < 0.05, ]
    if (nrow(c) == 1) c$tree_id else NA_character_ }, character(1))
}
F_TF <- "data/raw/field_data/ipad_data/treeflux_total.xlsx"
tf <- suppressMessages(readxl::read_excel(F_TF, sheet = "Sheet1"))
tf <- tibble(tag = trimws(as.character(tf$`Plot Tag`)), date = as.Date(tf$Date), notes = as.character(tf$Notes)) %>%
  filter(.is_tag(tag), !is.na(date)) %>% mutate(tree_id = tt_tree(F_TF, "Plot Tag", tag))
F_MU <- "data/processed/flux/semirigid_tree_final_complete_dataset_with_untagged.csv"
mu <- read.csv(F_MU, stringsAsFactors = FALSE) %>% filter(grepl("^UNTAG_", Plot.Tag)) %>%
  transmute(tree_id = tt_tree(F_MU, "Plot.Tag", Plot.Tag), date = as.Date(substr(Date, 1, 10)))
F_H <- "data/processed/flux/tree_flux_2021_multiheight.csv"
hh <- read.csv(F_H, stringsAsFactors = FALSE) %>%
  transmute(name = tree_id, tree_id = tt_tree(F_H, "tree_id", tree_id), date = as.Date(substr(start.time, 1, 10)))
F_Y <- "data/raw/field_data/static_chamber_field/Compiled YMF Data 2023 - Sheet1.csv"
yy <- read.csv(F_Y, check.names = FALSE, stringsAsFactors = FALSE)
yy <- yy %>% filter(!(is.na(`Flux Start (System) Time`) | trimws(`Flux Start (System) Time`) %in% c("", "NA"))) %>%
  transmute(tree_id = tt_tree(F_Y, "Tree Tag", `Tree Tag`, `Species Code`, `DBH (cm)`), date = as.Date(Date, "%m/%d/%y"),
            cond = grepl("^dead$", trimws(`Bark Missing`), ignore.case = TRUE) | grepl("^dead$", trimws(`Wounding/Holes`), ignore.case = TRUE),
            note = grepl("dead|snag", Notes, ignore.case = TRUE))
VISITS <- bind_rows(tf %>% select(tree_id, date), mu, hh %>% select(tree_id, date), yy %>% select(tree_id, date)) %>%
  filter(!is.na(tree_id)) %>% distinct()
MARKS <- bind_rows(
  tf %>% filter(grepl("dead|snag", notes, ignore.case = TRUE)) %>% transmute(tree_id, date, what = "monthly field notes"),
  mu %>% filter(grepl("_dead$", tree_id)) %>% transmute(tree_id, date, what = "monthly label _dead"),
  hh %>% filter(grepl("dead", name, ignore.case = TRUE)) %>% transmute(tree_id, date, what = "2021 name"),
  yy %>% filter(cond) %>% transmute(tree_id, date, what = "2023 condition score"),
  yy %>% filter(note) %>% transmute(tree_id, date, what = "2023 notes")) %>%
  filter(!is.na(tree_id)) %>% distinct()
DEAD <- MARKS %>% group_by(tree_id) %>%
  summarise(first = min(date), n_dates = n_distinct(date), .groups = "drop")
DEAD$dead_evidence <- sapply(DEAD$tree_id, function(t) {
  m <- MARKS[MARKS$tree_id == t, ]
  paste(sapply(split(m$date, m$what), function(d) paste(format(sort(unique(d))), collapse = ", ")) %>%
          { paste0(names(.), ": dead ", .) }, collapse = "; ") })
DEAD$n_after <- sapply(seq_len(nrow(DEAD)), function(i) sum(VISITS$tree_id == DEAD$tree_id[i] & VISITS$date > DEAD$first[i]))
DEAD <- DEAD %>% mutate(
  dead = n_dates > 1 | n_after == 0,
  dead_from = as.Date(ifelse(dead, first, NA)),
  dead_basis = case_when(n_dates > 1 ~ sprintf("marked dead on %d dates; dead from %s", n_dates, first),
                         n_after == 0 ~ sprintf("marked dead once (%s); not measured after", first),
                         TRUE ~ sprintf("marked dead once (%s); measured again on %d later date(s) without the mark: live", first, n_after)))
# every tree the field records call dead has a dated mark
old_dead <- TREES$tree_id[!is.na(TREES$dead_evidence)]
if (length(setdiff(old_dead, DEAD$tree_id))) stop("dead evidence without a dated mark: ", paste(setdiff(old_dead, DEAD$tree_id), collapse = ", "))
TREES <- TREES %>% select(-dead, -dead_evidence) %>% left_join(DEAD %>% select(tree_id, dead, dead_from, dead_evidence, dead_basis), by = "tree_id") %>%
  mutate(dead = coalesce(dead, FALSE)) %>%
  relocate(dead, dead_from, dead_evidence, dead_basis, .after = coord_source)

# ==============================================================================
# 5c. DBH NOT RECORDED: left missing (no imputation), and said so in dbh_source
# ==============================================================================
NO_DBH <- "none: no DBH measurement recorded"
stopifnot(all(is.na(TREES$dbh_source[is.na(TREES$dbh_cm)])), all(!is.na(TREES$dbh_source[!is.na(TREES$dbh_cm)])))
TREES$dbh_source[is.na(TREES$dbh_cm)] <- NO_DBH
FLUX_CAMPAIGNS <- c(MON, F21, C23, OAK)
NODBH <- TREES %>% filter(is.na(dbh_cm)) %>%
  transmute(tree_id, species_code, campaigns,
            flux_measured = sapply(strsplit(campaigns, ";"), function(z) any(z %in% FLUX_CAMPAIGNS)))

# ==============================================================================
# 6. INVENTORY TABLE
# ==============================================================================
INV <- I %>% transmute(stem_id, tree_id, source, tag, stem, species_code, quadrat, PX, PY,
                       latitude = round(latitude, 7), longitude = round(longitude, 7), position, located,
                       dbh_cm_raw, dbh_cm, dbh_shifts, in_notch, in_square, in_stand)
stopifnot(all(TREES$in_inventory[!is.na(TREES$in_inventory)] %in% INV$stem_id),
          all(INV$tree_id[match(na.omit(TREES$in_inventory), INV$stem_id)] == TREES$tree_id[!is.na(TREES$in_inventory)]))

write_csv(INV, file.path(OUT, "inventory_stems.csv"), na = "")
write_csv(TREES, file.path(OUT, "trees.csv"), na = "")
write_csv(TT, file.path(OUT, "tree_tags.csv"), na = "")
write_csv(UNRES, out_path("tree_tags_unresolved.csv", "audit"), na = "")

# ==============================================================================
# 7. AUDIT
# ==============================================================================
# species: field record vs census vs other campaigns
SPC <- bind_rows(
  T1 %>% filter(grepl("/", sp_all) | (!is.na(census_species) & census_species != species_code)) %>%
    transmute(tree_id, species_code, field_species = sp_all, census_species, source = "census tag"),
  T2 %>% filter(grepl("/", sp_all)) %>% transmute(tree_id, species_code, field_species = sp_all, source = "2021 survey"))

# 2021 species vs the tree-name prefix (e.g. RM105: recorded QURU, prefix RM = ACRU)
PFX <- TREES %>% filter(tree_id %in% T2$tree_id) %>% mutate(prefix = prefix_species(tree_id)) %>%
  filter(!is.na(prefix), prefix != species_code, !tree_id %in% SPC$tree_id) %>%
  transmute(tree_id, species_code, field_species = species_code, census_species = NA_character_,
            source = paste0("2021 survey: name prefix says ", prefix))
SPC <- bind_rows(SPC, PFX)
write_csv(SPC, out_path("tree_species_conflicts.csv", "audit"), na = "")

# dead: every tree with a mark, its evidence and the outcome
DS <- TREES %>% filter(tree_id %in% DEAD$tree_id) %>%
  transmute(tree_id, species_code, campaigns, dead, dead_from, dead_basis, dead_evidence)
write_csv(DS, out_path("tree_dead_status.csv", "audit"), na = "")

sink(out_path("trees_build.txt", "audit")); op <- options(width = 160)
cat(sprintf("inventory_stems %d | trees %d | tree_tags %d spellings\n\n", nrow(INV), nrow(TREES), nrow(TT)))
cat("Trees per campaign:\n")
cc <- table(unlist(strsplit(TREES$campaigns, ";"))); print(cc)
cat(sprintf("\nTag collisions: %d tags, %d records; %d tags split into >1 tree, %d tags merged into one tree; %d records unsure\n",
            n_distinct(COL$tag), nrow(COL), n_distinct(COL$tag[COL$n_trees > 1]), n_distinct(COL$tag[COL$n_trees == 1]),
            sum(COL$unsure)))
cat(sprintf("  decisions file: %d rows (%s); overrides of the proposal: %d; stale rows: %d\n", nrow(DEC),
            paste(names(table(DEC$status)), table(DEC$status), collapse = ", "),
            sum(TA$tree_id != TA$proposed_tree_id), length(stale)))
cat(sprintf("\nTrees marked dead in some record: %d; dead by the rule: %d\n", nrow(DS), sum(DS$dead)))
print(as.data.frame(DS %>% select(-campaigns)), row.names = FALSE)
cat(sprintf("\nNo DBH measurement recorded (dbh_cm NA; not imputed): %d trees, %d of them flux-measured\n",
            nrow(NODBH), sum(NODBH$flux_measured)))
if (nrow(NODBH)) print(as.data.frame(NODBH), row.names = FALSE)
cat("\nDBH sources:\n"); print(as.data.frame(TREES %>% count(dbh_source)), row.names = FALSE)
cat("\nSpecies sources:\n"); print(as.data.frame(TREES %>% count(species_source)), row.names = FALSE)
cat(sprintf("\nSpecies disagreements between records: %d trees (outputs/audit/tree_species_conflicts.csv)\n", nrow(SPC)))
cat(sprintf("\nSpellings not mapped to a measured tree: %d (outputs/audit/tree_tags_unresolved.csv)\n", nrow(UNRES)))
print(as.data.frame(UNRES %>% count(reason)), row.names = FALSE)
cat("\nTrees without species, DBH or coordinates:\n")
cat(sprintf("  species %d | DBH %d | coordinates %d\n", sum(is.na(TREES$species_code)), sum(is.na(TREES$dbh_cm)),
            sum(is.na(TREES$latitude))))
options(op); sink()
cat(readLines(out_path("trees_build.txt", "audit")), sep = "\n")
