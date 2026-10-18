# ==============================================================================
# fig05_methane-cycling.R -- Figure 5: methanogens and methanotrophs in the 16S data
# ------------------------------------------------------------------------------
# Stage E. Run from the repository root.
#   (a) methanogen relative abundance by species and compartment (mean +/- SE)
#   (b) methanogen family composition, families shaded by methanogenic pathway
#   (c) methanotroph relative abundance, Known and Putative stacked (error bar on
#       the total)
#   (d) methanotroph family composition, Known and Putative separately
# Methanogens: the families of code/lib/methanogen_families.R. Methanotrophs: the
# Known / Putative classification of code/lib/load_methanotroph_definitions.R
# (include_conditional = TRUE).
# Inputs: the 16S phyloseq object (code/lib/build_16s_phyloseq(): data/tables/
#   asv_counts_16S, asv_taxonomy_16S, samples_16S, samples, trees; organelles removed),
#   rarefied with rarefy_16s(); species names from data/tables/species.csv
# Output: outputs/figures/generated/fig5_final.png
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(patchwork)
library(viridis)
source("code/lib/build_phyloseq.R")
source("code/lib/tables_read.R")
source("code/lib/methanogen_families.R")
source("code/lib/load_methanotroph_definitions.R")

# ==============================================================================
# STEP 1: 16S relative abundance (rarefied), survey wood and soil samples
# ==============================================================================

no_mito <- build_16s_phyloseq()$ps
taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
ps.rare <- rarefy_16s(no_mito)
ps.ra <- transform_sample_counts(ps.rare, function(x) x / sum(x) * 100)

# Filter samples
samp_df <- data.frame(sample_data(ps.ra), stringsAsFactors = FALSE)
samp_df$compartment <- case_when(
  samp_df$core_type == "Inner" ~ "Heartwood",
  samp_df$core_type == "Outer" ~ "Sapwood",
  samp_df$core_type == "Mineral" ~ "Mineral Soil",
  samp_df$core_type == "Organic" ~ "Organic Soil",
  TRUE ~ NA_character_
)
samp_df <- samp_df %>% filter(!is.na(compartment), !is.na(species.x), species.x != "")
sample_data(ps.ra) <- sample_data(samp_df)
ps.filt <- prune_samples(!is.na(sample_data(ps.ra)$compartment), ps.ra)

# Extract tables
otu_df <- as.data.frame(otu_table(ps.filt))
tax_df <- as.data.frame(tax_table(ps.filt))
samp_meta <- data.frame(sample_data(ps.filt))

species_mapping <- SPECIES_ABBREV   # "A. rubrum", from data/tables/species.csv

# ==============================================================================
# STEP 2: Methanogen abundance, total and per family
# ==============================================================================

methanogen_families <- METHANOGEN_FAMILIES

mg_asvs <- rownames(tax_df)[tax_df$Family %in% methanogen_families]
samp_meta$methanogen_pct <- colSums(otu_df[mg_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]

for (fam in methanogen_families) {
  fam_asvs <- rownames(tax_df)[which(tax_df$Family == fam)]
  samp_meta[[paste0("mg_", fam)]] <- colSums(otu_df[fam_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]
}

# ==============================================================================
# STEP 3: Methanotroph abundance (Known, Putative), total and per family x status
# ==============================================================================

mt_defs <- load_methanotroph_defs()

# Classify each ASV as Known, Putative, or NA
tax_df$mt_status <- classify_methanotrophs(tax_df, mt_defs, include_conditional = TRUE)

known_asvs   <- rownames(tax_df)[tax_df$mt_status == "Known" & !is.na(tax_df$mt_status)]
putative_asvs <- rownames(tax_df)[tax_df$mt_status == "Putative" & !is.na(tax_df$mt_status)]
all_mt_asvs  <- c(known_asvs, putative_asvs)

cat("Known methanotroph ASVs:", length(known_asvs), "\n")
cat("Putative methanotroph ASVs:", length(putative_asvs), "\n")

# Per-sample abundance: total, known, putative
samp_meta$methanotroph_pct <- colSums(otu_df[all_mt_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]
samp_meta$mt_known_pct     <- colSums(otu_df[known_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]
samp_meta$mt_putative_pct  <- colSums(otu_df[putative_asvs, , drop = FALSE], na.rm = TRUE)[rownames(samp_meta)]

# per family x status: columns mt_<status>_<family>
tax_df$mt_family <- NA_character_
tax_df$mt_family[!is.na(tax_df$mt_status)] <- tax_df$Family[!is.na(tax_df$mt_status)]
mt_fam_status_unique <- unique(tax_df[!is.na(tax_df$mt_status), c("mt_family", "mt_status")])
for (i in seq_len(nrow(mt_fam_status_unique))) {
  fam <- mt_fam_status_unique$mt_family[i]
  stat <- mt_fam_status_unique$mt_status[i]
  fam_stat_asvs <- rownames(tax_df)[which(tax_df$mt_family == fam & tax_df$mt_status == stat)]
  samp_meta[[paste0("mt_", stat, "_", fam)]] <- colSums(otu_df[fam_stat_asvs, , drop = FALSE],
                                                        na.rm = TRUE)[rownames(samp_meta)]
}

# ==============================================================================
# STEP 4: Prepare plot data
# ==============================================================================

plot_data <- samp_meta
plot_data$species_label <- species_mapping[plot_data$species.x]
plot_data <- plot_data %>% filter(!is.na(species_label))
plot_data$compartment <- factor(plot_data$compartment,
                                 levels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))
sp_order <- sort(unique(plot_data$species_label))
plot_data$species_label <- factor(plot_data$species_label, levels = sp_order)

# ==============================================================================
# Panel (a) — Methanogen relative abundance
# ==============================================================================

summary_mg_a <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_pct = mean(methanogen_pct, na.rm = TRUE),
            se_pct = sd(methanogen_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# solid fill — compartment is already given by the facet strip (no ambiguous colour mapping)
p_mg_a <- ggplot(summary_mg_a, aes(x = species_label, y = mean_pct)) +
  geom_col(width = 0.7, fill = "#B2182B", color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(mean_pct - se_pct, 0), ymax = mean_pct + se_pct),
                width = 0.2, linewidth = 0.3) +
  facet_wrap(~ compartment, nrow = 1) +
  labs(x = NULL, y = "Methanogen relative\nabundance (%)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12.8),
        legend.position = "none",
        strip.text = element_text(size = 11.2, face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(5, 10, 0, 10))

# ==============================================================================
# Panel (b) — Methanogen family composition
# ==============================================================================

mg_family_cols <- paste0("mg_", methanogen_families)
mg_family_cols <- intersect(mg_family_cols, colnames(plot_data))

mg_comp <- plot_data %>%
  select(species_label, compartment, all_of(mg_family_cols)) %>%
  pivot_longer(cols = all_of(mg_family_cols), names_to = "Family", values_to = "Abundance") %>%
  mutate(Family = sub("^mg_", "", Family))

summary_mg_b <- mg_comp %>%
  group_by(species_label, compartment, Family) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(species_label, compartment) %>%
  mutate(total = sum(mean_abund),
         proportion = if_else(total > 0, mean_abund / total * 100, 0)) %>%
  ungroup()

# Methanogen family colours, shaded by methanogenic pathway (H2, methyl, acetate)
mg_pathway <- c(Methanobacteriaceae = "H", Methanocellaceae = "H", Methanocorpusculaceae = "H",
                Methanomicrobiaceae = "H", Methanoregulaceae = "H",
                Methanomassiliicoccaceae = "M", Methanomethyliaceae = "M",
                Methanosaetaceae = "A", Methanosarcinaceae = "A")
stopifnot(setequal(names(mg_pathway), methanogen_families))
mg_families_sorted <- names(mg_pathway)
mg_colors <- c(Methanobacteriaceae = "#F4C430", Methanocellaceae = "#E8973A", Methanocorpusculaceae = "#D4733A",
               Methanomicrobiaceae = "#B5562C", Methanoregulaceae = "#8C3B1E",
               Methanomassiliicoccaceae = "#3FA796", Methanomethyliaceae = "#1F6E63",
               Methanosaetaceae = "#8E6BBF", Methanosarcinaceae = "#5B3F8C")
mg_labels <- setNames(paste0(mg_families_sorted, c(H = " (H\u2082)", M = " (methyl)", A = " (acetate)")[mg_pathway[mg_families_sorted]]),
                      mg_families_sorted)
summary_mg_b$Family <- factor(summary_mg_b$Family, levels = mg_families_sorted)

p_mg_b <- ggplot(summary_mg_b, aes(x = species_label, y = proportion, fill = Family)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = mg_colors, labels = mg_labels, name = "Methanogen\nfamily") +
  labs(x = NULL, y = "Proportion of\nmethanogens (%)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12.8),
        legend.position = "right",
        legend.justification = c(0, 1),
        legend.key.size = unit(0.38, "cm"),
        legend.text = element_text(size = 10.4),
        legend.title = element_text(size = 11.2, face = "bold"),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(0, 10, 0, 10))

# ==============================================================================
# Panel (c) — Methanotroph relative abundance (Known + Putative stacked)
# ==============================================================================

summary_mt_known <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_pct = mean(mt_known_pct, na.rm = TRUE),
            se_pct = sd(mt_known_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(Status = "Known")

summary_mt_putative <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_pct = mean(mt_putative_pct, na.rm = TRUE),
            se_pct = sd(mt_putative_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(Status = "Putative")

summary_mt_a <- bind_rows(summary_mt_known, summary_mt_putative) %>%
  mutate(Status = factor(Status, levels = c("Known", "Putative")))

# Error bars on total (Known + Putative combined)
summary_mt_total <- plot_data %>%
  group_by(species_label, compartment) %>%
  summarize(mean_total = mean(methanotroph_pct, na.rm = TRUE),
            se_total = sd(methanotroph_pct, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# colour encodes CLASSIFICATION only (Known vs Putative); compartment is given by the facet strip
panel_c_fills <- c("Known" = "#2166AC", "Putative" = "grey70")
summary_mt_a$fill_group <- factor(summary_mt_a$Status, levels = c("Known", "Putative"))

p_mt_a <- ggplot(summary_mt_a, aes(x = species_label, y = mean_pct, fill = fill_group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_errorbar(data = summary_mt_total,
                aes(x = species_label, y = mean_total,
                    ymin = pmax(mean_total - se_total, 0),
                    ymax = mean_total + se_total),
                inherit.aes = FALSE, width = 0.2, linewidth = 0.3) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = panel_c_fills, name = "Classification") +
  labs(x = NULL, y = "Methanotroph relative\nabundance (%)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12.8),
        legend.position = "right",
        legend.justification = c(0, 0.5),
        legend.text = element_text(size = 10.4),
        legend.title = element_text(size = 11.2, face = "bold"),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(0, 10, 0, 10))

# ==============================================================================
# Panel (d) — Methanotroph family composition (Known families first, then Putative)
# ==============================================================================

mt_fam_stat_cols <- grep("^mt_(Known|Putative)_", colnames(plot_data), value = TRUE)

mt_comp <- plot_data %>%
  select(species_label, compartment, all_of(mt_fam_stat_cols)) %>%
  pivot_longer(cols = all_of(mt_fam_stat_cols),
               names_to = "Family_Status", values_to = "Abundance") %>%
  mutate(
    Status = sub("^mt_(Known|Putative)_.*", "\\1", Family_Status),
    Family = sub("^mt_(Known|Putative)_", "", Family_Status)
  )

summary_mt_b <- mt_comp %>%
  group_by(species_label, compartment, Family, Status) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(species_label, compartment) %>%
  mutate(total = sum(mean_abund),
         proportion = if_else(total > 0, mean_abund / total * 100, 0)) %>%
  ungroup()

# Drop families with zero total abundance (not present in plotted samples)
fam_totals <- summary_mt_b %>% group_by(Family, Status) %>%
  summarize(total_abund = sum(mean_abund), .groups = "drop") %>%
  filter(total_abund > 0)
summary_mt_b <- summary_mt_b %>%
  semi_join(fam_totals, by = c("Family", "Status"))

# Distinct color for each Family x Status combination
all_mt_fams <- sort(unique(summary_mt_b$Family))
summary_mt_b$Family_status <- paste0(summary_mt_b$Family, " (", summary_mt_b$Status, ")")

# Order levels: for each family, Known first then Putative
status_levels <- c()
for (fam in all_mt_fams) {
  if (any(summary_mt_b$Family_status == paste0(fam, " (Known)")))
    status_levels <- c(status_levels, paste0(fam, " (Known)"))
  if (any(summary_mt_b$Family_status == paste0(fam, " (Putative)")))
    status_levels <- c(status_levels, paste0(fam, " (Putative)"))
}
summary_mt_b$Family_status <- factor(summary_mt_b$Family_status, levels = status_levels)

mt_colors <- setNames(viridis(length(status_levels), option = "D", direction = -1, alpha = 0.7), status_levels)

p_mt_b <- ggplot(summary_mt_b, aes(x = species_label, y = proportion, fill = Family_status)) +
  geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = mt_colors, name = "Methanotroph\nfamily") +
  labs(x = NULL, y = "Proportion of\nmethanotrophs (%)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.6, face = "italic"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12.8),
        legend.position = "right",
        legend.justification = c(0, 1),
        legend.key.size = unit(0.38, "cm"),
        legend.text = element_text(size = 10.4),
        legend.title = element_text(size = 11.2, face = "bold"),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        plot.margin = margin(0, 10, 5, 10))

# ==============================================================================
# Combine all panels vertically with shared x-axis
# ==============================================================================

# Vertical stack: a (methanogen abundance) / b (methanogen composition) /
#                 c (methanotroph abundance) / d (methanotroph composition)
# Compartment facet headers only on panel (a), species labels only on panel (d)
# Add bold tags to each panel individually (ggplot2 4.x compatible)
tag_theme <- theme(plot.tag = element_text(size = 12, face = "bold"))
p_mg_a <- p_mg_a + tag_theme
p_mg_b <- p_mg_b + tag_theme
p_mt_a <- p_mt_a + tag_theme
p_mt_b <- p_mt_b + tag_theme

fig_combined <- p_mg_a / p_mg_b / p_mt_a / p_mt_b +
  plot_layout(heights = c(1, 1.2, 1, 1.2)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave("outputs/figures/generated/fig5_final.png", fig_combined, width = 10.5, height = 9.5, dpi = 300, bg = "white")   # near print size

cat("Wrote outputs/figures/generated/fig5_final.png\n")
