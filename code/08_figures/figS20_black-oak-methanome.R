# ==============================================================================
# figS20_black-oak-methanome.R -- Figure S21: methane-cycling taxa across tissues of
# the felled black oak (QUVE)
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
# Heatmap of methanogen, Known methanotroph and Putative methanotroph genera (Table S2
# classification, include_conditional = TRUE as Fig. 5) in each tissue of the felled oak.
#
# Value per taxon (family | genus) and tissue: relative abundance in percent,
#   100 x (reads of the taxon summed over the tissue's libraries) /
#         (non-organelle 16S reads summed over the same libraries)
# i.e. pooled within tissue, not a mean of per-library shares. Non-organelle reads are
# asv_reads() (organelles by is_organelle(), the denominator of every 16S relative
# abundance). Colour is on a log10 scale; a taxon with no reads in a tissue is grey.
#
# Inputs (data/tables/ via code/lib/tables_read.R): asv_reads("felled_oak"),
#   dna_samples("2022_felled_oak") (library -> tissue)
# Output: outputs/figures/generated/black_oak_methanome_revised.png
# ==============================================================================
library(tidyverse)
source("code/lib/tables_read.R")
source("code/lib/methanogen_families.R")
source("code/lib/load_methanotroph_definitions.R")

# tissue of each library, grouped for display: branch bark with Bark, branch interior with
# Branch, root epidermis and interior with Coarse Root, heart rot as Rot
MEDIA <- c(Heartwood = "Heartwood", Sapwood = "Sapwood", Litter = "Leaf Litter", Bark = "Bark",
           `Branch bark` = "Bark", Foliage = "Foliage", Mineral = "Mineral Soil", Organic = "Organic Soil",
           Branch = "Branch", `Branch interior` = "Branch", `Coarse root` = "Coarse Root",
           `Coarse root epidermis` = "Coarse Root", `Coarse root interior` = "Coarse Root",
           `Fine root` = "Fine Root", `Heart rot` = "Rot")
media_order <- c("Heartwood", "Sapwood", "Rot", "Bark", "Branch", "Foliage", "Leaf Litter",
                 "Coarse Root", "Fine Root", "Organic Soil", "Mineral Soil")   # interior -> exterior -> ground

A  <- asv_reads("felled_oak", samples_only = TRUE)            # non-organelle counts; A$tot per library
bo <- dna_samples("2022_felled_oak")
media <- unname(MEDIA[bo$core_type[match(colnames(A$cnt), bo$library_id)]])
stopifnot(!anyNA(media), setequal(media, media_order))

st <- classify_methanotrophs(data.frame(Family = A$tax$Family, Genus = A$tax$Genus, stringsAsFactors = FALSE),
                             load_methanotroph_defs(), include_conditional = TRUE)
grp <- case_when(A$tax$Family %in% METHANOGEN_FAMILIES ~ "Methanogen",
                 st %in% "Known" ~ "Methanotroph (Known)",
                 st %in% "Putative" ~ "Methanotroph (Putative)")
k <- !is.na(grp) & !is.na(A$tax$Family)

tissue_tot <- tapply(A$tot, media, sum)                        # pooled denominator per tissue
taxon_reads <- rowsum(A$cnt[k, , drop = FALSE],                # reads per taxon (family | genus) x library
                      paste(grp[k], A$tax$Family[k], coalesce(A$tax$Genus[k], ""), sep = " | "))
tis_reads <- t(rowsum(t(taxon_reads), media))                 # taxon x tissue
pct <- sweep(tis_reads, 2, tissue_tot[colnames(tis_reads)], "/") * 100

hm <- as.data.frame(as.table(pct), stringsAsFactors = FALSE) %>%
  setNames(c("key", "Media", "pct")) %>%
  separate(key, c("Group", "Family", "Genus"), sep = " \\| ", fill = "right") %>%
  mutate(Media = factor(Media, levels = media_order),
         Row_Label = paste(Family, coalesce(Genus, ""), sep = " | "),
         fill_pct = ifelse(pct > 0, pct, NA_real_))
nz <- range(hm$fill_pct, na.rm = TRUE)
cat(sprintf("Fig S21: %d taxa x %d tissues; nonzero %.2g%% to %.2g%%; %d of %d cells zero\n",
            n_distinct(hm$Row_Label), length(media_order), nz[1], nz[2], sum(hm$pct == 0), nrow(hm)))

BREAKS <- 10^(floor(log10(nz[1])):ceiling(log10(nz[2])))
pct_lab <- function(x) paste0(format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE), "%")
taxon_labels <- function(v) do.call(expression, lapply(v, function(x) {
  p <- strsplit(x, " | ", fixed = TRUE)[[1]]; gen <- if (length(p) > 1) trimws(p[2]) else ""
  if (gen == "") bquote(.(paste0(p[1], " (genus unresolved)"))) else bquote(.(paste0(p[1], ": ")) * italic(.(gen)))
}))

groups <- c("Methanogen", "Methanotroph (Known)", "Methanotroph (Putative)")
groups <- groups[groups %in% hm$Group]
hm$Panel <- factor(hm$Group, levels = groups,
                   labels = paste0("(", letters[seq_along(groups)], ") ", groups))
# one plot faceted by group: columns align across panels and rows keep equal height
p <- ggplot(hm, aes(Media, Row_Label, fill = fill_pct)) +
  geom_tile(colour = "white") +
  geom_point(data = filter(hm, pct == 0), aes(shape = "Not detected"), alpha = 0) +   # legend key for grey
  facet_grid(Panel ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = taxon_labels) +
  scale_fill_viridis_c(option = "C", trans = "log10", limits = range(BREAKS), breaks = BREAKS,
                       labels = pct_lab, na.value = "grey85", oob = scales::squish,
                       name = "Relative abundance\n(% of non-organelle\n16S reads)") +
  scale_shape_manual(values = 15, name = NULL) +
  guides(fill = guide_colourbar(order = 1, barheight = unit(40, "mm")),
         shape = guide_legend(order = 2, override.aes = list(alpha = 1, size = 6, colour = "grey85"))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
        panel.spacing.y = unit(4, "mm"))

ggsave("outputs/figures/generated/black_oak_methanome_revised.png", p, width = 11, height = 5, dpi = 300, bg = "white")
