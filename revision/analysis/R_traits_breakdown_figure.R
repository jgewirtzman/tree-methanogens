# ==============================================================================
# REVISION — "How it breaks" figure: trait effects collapse without the 2 conifers
# ==============================================================================
# Panel for the SI walk-through (heatmap -> ordination -> THIS -> "mostly
# unexplained"). Dumbbell: each key species-level relationship, Spearman rho with
# ALL species vs ANGIOSPERM-ONLY. Most collapse toward 0 once the 2 conifers are
# removed; only moisture->methanogen survives.
#
# NEW file. Run: Rscript revision/analysis/R_traits_breakdown_figure.R
# Output: revision/outputs/traits_breakdown.png
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",
  BELE="Betula lenta",BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",
  FRAM="Fraxinus americana",KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",
  QUAL="Quercus alba",QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",
  TSCA="Tsuga canadensis")

source("revision/analysis/_prep_species_data.R")
micro <- tree_level_complete %>% group_by(species_id) %>%
  summarise(meth = median(log_meth), mcra = median(log_mcra), balance = median(log_ratio), .groups="drop")
flux <- flux_by_species %>% transmute(species_id, flux = median_flux)
tr <- read_csv(TRAITS, show_col_types = FALSE) %>% group_by(spcode) %>% slice(1) %>% ungroup() %>%
  transmute(species_id = spcode, wood_density = wood_density_gcm3, bark_density = bark_density_gcm3, gymnosperm)
niche <- read_csv("revision/outputs/tree_species_moisture_niche.csv", show_col_types = FALSE) %>%
  mutate(species_id = names(sp_map)[match(species, sp_map)]) %>% transmute(species_id, vwc = vwc_mean)
sp <- micro %>% full_join(flux, by="species_id") %>% left_join(tr, by="species_id") %>%
  left_join(niche, by="species_id") %>% mutate(clade = if_else(gymnosperm==1,"Gymnosperm","Angiosperm"))

rho <- function(d, y, x) { d <- d %>% filter(is.finite(.data[[y]]), is.finite(.data[[x]]))
  if (nrow(d) < 5) return(NA_real_)
  suppressWarnings(cor(d[[x]], d[[y]], method = "spearman")) }

rel <- tribble(~label, ~y, ~x,
  "Moisture -> methanogen (mcrA)", "mcra", "vwc",
  "Wood density -> methanotroph",  "meth", "wood_density",
  "Bark density -> methanotroph",  "meth", "bark_density",
  "Bark density -> balance",       "balance", "bark_density",
  "Wood density -> balance",       "balance", "wood_density",
  "Bark density -> stem flux",     "flux", "bark_density",
  "Wood density -> stem flux",     "flux", "wood_density")

d <- rel %>% rowwise() %>%
  mutate(all = rho(sp, y, x), ang = rho(sp %>% filter(clade=="Angiosperm"), y, x)) %>%
  ungroup() %>%
  # "retained" = direction kept and >=60% of the effect size survives; else collapses/flips
  mutate(retained = sign(ang) == sign(all) & abs(ang) >= 0.30 & abs(ang) >= 0.6*abs(all),
         tag = ifelse(retained, "retained", "collapses / flips"),
         label = fct_reorder(label, all))

long <- d %>% select(label, all, ang) %>%
  pivot_longer(c(all, ang), names_to = "set", values_to = "rho") %>%
  mutate(set = recode(set, all = "All species (incl. 2 conifers)", ang = "Angiosperms only"))

p <- ggplot(d) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_segment(aes(x = all, xend = ang, y = label, yend = label), color = "grey55", linewidth = 0.8) +
  geom_point(data = long, aes(rho, label, shape = set, fill = set), size = 3.2) +
  geom_text(data = d, aes(x = pmax(all, ang) + 0.05, y = label, label = tag,
            color = retained), hjust = 0, size = 2.7, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "#1a9850", `FALSE` = "#b2182b")) +
  scale_shape_manual(values = c(`All species (incl. 2 conifers)` = 21, `Angiosperms only` = 24), name = NULL) +
  scale_fill_manual(values = c(`All species (incl. 2 conifers)` = "grey30", `Angiosperms only` = "#f1a340"), name = NULL) +
  scale_x_continuous(limits = c(-0.65, 0.95), breaks = seq(-0.6, 0.8, 0.2),
                     labels = function(z) sprintf("%.1f", z)) +
  labs(x = "Spearman rho (species level)", y = NULL,
       title = "Trait-community associations mostly collapse without the 2 conifers",
       subtitle = "Grey = all species; orange = angiosperms only. Density effects on the balance and flux collapse or flip;\nonly the moisture->methanogen and (weaker) density->methanotroph directions persist. Among angiosperms (n<=13)\nNONE is statistically significant -> most among-species variation is not yet explained by measured traits.") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.3),
        plot.margin = margin(8, 30, 8, 8))
ggsave(file.path(out, "traits_breakdown.png"), p, width = 9, height = 5, dpi = 300, bg = "white")
cat("Wrote", file.path(out, "traits_breakdown.png"), "\n")
print(as.data.frame(d %>% transmute(label, rho_all = round(all,2), rho_angio = round(ang,2), tag)), row.names = FALSE)
