# ==============================================================================
# REVISION — Trait x methane-cycling HEATMAP with TWO significance markers.
# Supersedes traits_heatmap.png + the (misleading) traits_breakdown.png by folding
# the clade-robustness INTO one heatmap:
#   *  univariate Spearman p < 0.05 (all n<=10 species)
#   †  survives control for the gymnosperm/angiosperm split
#      (rank-partial: rank(response) ~ rank(trait) + gymnosperm, uses all n, 1 df)
# The partial-clade test is the defensible small-n robustness check (keeps all
# species, costs 1 df) for "does the trait predict BEYOND the deepest clade split".
# Cells that are †-only from a near-zero univariate rho are suppression artifacts
# (few df) and are NOT interpreted; the honest reading is: the density->methanotroph
# suppression is largely the conifer contrast (loses significance under †), whereas
# PLANT LONGEVITY -> balance/methanotroph is robust (univariate + angiosperm-only +
# clade-controlled). n<=10, EXPLORATORY / hypothesis-generating.
# NEW file; original traits_heatmap.R untouched.
# Output: revision/outputs/traits_heatmap_robust.png
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse) })
out <- "revision/outputs"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
TRAITS <- "/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv"
sp_map <- c(ACRU="Acer rubrum",ACSA="Acer saccharum",BEAL="Betula alleghaniensis",BELE="Betula lenta",
  BEPA="Betula papyrifera",CAOV="Carya ovata",FAGR="Fagus grandifolia",FRAM="Fraxinus americana",
  KALA="Kalmia latifolia",PIST="Pinus strobus",PRSE="Prunus serotina",QUAL="Quercus alba",
  QURU="Quercus rubra",QUVE="Quercus velutina",SAAL="Sassafras albidum",TSCA="Tsuga canadensis")

source("revision/analysis/_prep_species_data.R")
resp <- analysis_ratio %>% transmute(species_id, `Stem CH4 flux`=median_flux, `Balance (mcrA:MOB)`=median_log_ratio) %>%
  left_join(analysis_mcra %>% transmute(species_id, `mcrA (methanogen)`=log10(value+1)), by="species_id") %>%
  left_join(analysis_meth %>% transmute(species_id, `Methanotroph`=log10(value+1)), by="species_id")
keep <- c("wood_density_gcm3","bark_density_gcm3","bark_wood_ratio","porosity_num","gymnosperm",
          "wood_sapwood_pH","wood_heartwood_pH","try_wood_CN_ratio","try_stem_C","try_antifungal",
          "try_CWD_stem_decomp_rate_k","try_rooting_depth","try_fine_root_diameter",
          "try_fine_root_tissue_density","try_plant_height","try_plant_longevity")
trait_lab <- c(wood_density_gcm3="Wood density", bark_density_gcm3="Bark density", bark_wood_ratio="Bark:wood ratio",
  porosity_num="Wood porosity", gymnosperm="Gymnosperm", wood_sapwood_pH="Sapwood pH", wood_heartwood_pH="Heartwood pH",
  try_wood_CN_ratio="Wood C:N", try_stem_C="Stem C", try_antifungal="Antifungal chem.", try_CWD_stem_decomp_rate_k="CWD decomp rate",
  try_rooting_depth="Rooting depth", try_fine_root_diameter="Fine-root diameter", try_fine_root_tissue_density="Fine-root tissue density",
  try_plant_height="Plant height", try_plant_longevity="Plant longevity", vwc_realized="Realized soil moisture (VWC)")
cat_of <- c(`Wood density`="Structure",`Bark density`="Structure",`Bark:wood ratio`="Structure",`Wood porosity`="Structure",
  `Gymnosperm`="Structure",`Sapwood pH`="Chemistry",`Heartwood pH`="Chemistry",`Wood C:N`="Chemistry",`Stem C`="Chemistry",
  `Antifungal chem.`="Chemistry",`CWD decomp rate`="Chemistry",`Rooting depth`="Roots",`Fine-root diameter`="Roots",
  `Fine-root tissue density`="Roots",`Plant height`="Whole-plant",`Plant longevity`="Whole-plant",`Realized soil moisture (VWC)`="Moisture")

traits <- read_csv(TRAITS, show_col_types=FALSE) %>% filter(spcode %in% resp$species_id) %>%
  group_by(species_id=spcode) %>% slice(1) %>% ungroup() %>% select(species_id, any_of(keep))
niche <- read_csv("revision/outputs/tree_species_moisture_niche.csv", show_col_types=FALSE) %>%
  mutate(species_id = names(sp_map)[match(species, sp_map)]) %>% transmute(species_id, vwc_realized=vwc_mean)
dat <- resp %>% left_join(traits, by="species_id") %>% left_join(niche, by="species_id")
pred_cols <- c(keep, "vwc_realized"); resp_cols <- c("mcrA (methanogen)","Methanotroph","Balance (mcrA:MOB)","Stem CH4 flux")

grid <- expand_grid(trait=pred_cols, response=resp_cols) %>% rowwise() %>% mutate(
    x=list(dat[[trait]]), y=list(dat[[response]]), g=list(dat$gymnosperm)) %>%
  mutate(ok=list(is.finite(unlist(x)) & is.finite(unlist(y)))) %>%
  mutate(n=sum(unlist(ok)),
    rho = if(n>=6) suppressWarnings(cor(unlist(x)[unlist(ok)], unlist(y)[unlist(ok)], method="spearman")) else NA_real_,
    p_uni = if(n>=6) suppressWarnings(cor.test(unlist(x)[unlist(ok)], unlist(y)[unlist(ok)], method="spearman")$p.value) else NA_real_,
    p_clade = if(n>=6 && trait!="gymnosperm") tryCatch({
        xi<-unlist(x)[unlist(ok)]; yi<-unlist(y)[unlist(ok)]; gi<-unlist(g)[unlist(ok)]
        summary(lm(rank(yi) ~ rank(xi) + gi))$coef["rank(xi)",4] }, error=function(e) NA_real_) else NA_real_) %>%
  ungroup() %>% mutate(trait_label=trait_lab[trait], category=cat_of[trait_label],
    mark = ifelse(!is.na(p_uni)&p_uni<0.05,"*",""),
    # clade-robust = survives gymnosperm control AND univariate rho is meaningful (|rho|>=0.3;
    # excludes near-zero-rho suppression artifacts of the 1-df partial at n<=10)
    robust = !is.na(p_clade) & p_clade<0.05 & abs(rho)>=0.3,
    label = ifelse(is.na(rho),"", sprintf("%.2f%s", rho, mark)))
write.csv(grid %>% select(category,trait_label,response,n,rho,p_uni,p_clade),
          file.path(out,"traits_heatmap_robust_matrix.csv"), row.names=FALSE)

row_order <- grid %>% filter(response=="Stem CH4 flux") %>%
  arrange(factor(category, levels=c("Structure","Chemistry","Roots","Moisture","Whole-plant")), rho) %>% pull(trait_label)
grid <- grid %>% mutate(trait_label=factor(trait_label, levels=row_order),
  response=factor(response, levels=resp_cols), category=factor(category, levels=c("Structure","Chemistry","Roots","Moisture","Whole-plant")))

p <- ggplot(grid, aes(response, trait_label, fill=rho)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_tile(data=grid %>% filter(robust), fill=NA, color="black", linewidth=1.1) +
  geom_text(aes(label=label), size=2.6) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0, limits=c(-1,1), name="Spearman\nrho") +
  facet_grid(category ~ ., scales="free_y", space="free_y", switch="y") +
  scale_x_discrete(position="top") +
  labs(x=NULL, y=NULL, title="Plant traits structure the stem methane-cycling community (n<=10 species)",
       caption=paste("Cell = Spearman rho (species level); columns follow the microbial chain: methanogens -> methanotrophs -> balance -> net flux.",
         "\n*  univariate p<0.05.   Black outline = survives control for the gymnosperm/angiosperm split (rank-partial, all n, 1 df; |rho|>=0.3).",
         "\nEXPLORATORY (n<=10, nothing survives FDR): density->methanotroph is largely the conifer contrast (no outline); plant longevity is clade-robust.")) +
  theme_minimal(base_size=9) +
  theme(axis.text.x.top=element_text(angle=15, hjust=0, size=8.5), axis.text.y=element_text(size=8),
        strip.text.y.left=element_text(angle=0, face="bold", size=8), strip.placement="outside",
        panel.grid=element_blank(), plot.title=element_text(size=12, margin=margin(b=22)),
        plot.caption=element_text(size=7.3, hjust=0, margin=margin(t=10)), plot.margin=margin(t=10,r=45,b=6,l=6))
ggsave(file.path(out,"traits_heatmap_robust.png"), p, width=9.2, height=7.6, dpi=300, bg="white")
cat("clade-robust (dagger) cells:\n")
print(grid %>% filter(!is.na(p_clade) & p_clade<0.05) %>% select(trait_label,response,rho,p_uni,p_clade), row.names=FALSE)
cat("Wrote traits_heatmap_robust.png\n")
