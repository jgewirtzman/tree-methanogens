source("code/lib/outputs.R")
# ==============================================================================
# figS05_16s-ddpcr-association.R
# ------------------------------------------------------------------------------
# Does the 16S abundance of each classified group track the functional gene
# measured by ddPCR in the SAME sample? A check on the 16S classification, not a
# test of function: relative abundance and a gene ratio can co-vary with habitat as
# well as through gene carriage, which is why non-methanotroph control groups are
# shown alongside -- they set the size of association that co-occurrence alone
# produces.
#
#   response   gene copies per bacterial 16S copy (both by ddPCR, same extract), so
#              that a relative quantity is compared with a relative abundance
#   predictor  16S relative abundance of the group (% of reads, unrarefied)
#   test       Spearman, within compartment (wood and soil never pooled)
#   multiple   Benjamini-Hochberg across every test drawn in the figure
#   cells with fewer than 5 samples where the group is present are not tested
#
# Groups follow the capacity rule set (OPEN_DECISIONS, 2026-09-30): a taxon counts if it
# carries methane monooxygenase genes; tiers Known / Putative / listed-not-counted, classified by code/lib/load_methanotroph_definitions.R.
# Join of 16S and ddPCR samples as in code/09_tables_stats/05_table_ddpcr-16s-concordance.R.
#
# NEW file. Outputs: figSI_16s_ddpcr_association.png, 16s_ddpcr_association.csv
# ==============================================================================
suppressMessages({library(dplyr); library(ggplot2)})
options(warn = -1)
# rho and dagger glyphs need a UTF-8 character locale; under the C locale (non-interactive
# shells) they are silently dropped from the PNG.
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
num <- function(x) suppressWarnings(as.numeric(x))

# ---- 16S: relative abundance of each group, per sample -----------------------
otu  <- read.delim("data/raw/16s/OTU_table.txt", header = TRUE, row.names = 1, check.names = FALSE)
taxc <- intersect(c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), names(otu))
tax  <- otu[, taxc]
cnt  <- as.matrix(sapply(otu[, setdiff(names(otu), taxc), drop = FALSE], num))
rownames(cnt) <- rownames(otu); tot <- colSums(cnt, na.rm = TRUE)

source("code/lib/load_methanotroph_definitions.R")
st  <- classify_methanotrophs(tax, load_methanotroph_defs())
unr <- function(g) is.na(g) | g == "" | tolower(g) %in% c("none", "unclassified")
fam <- function(f) !is.na(tax$Family) & tax$Family == f
gen <- function(g) !is.na(tax$Genus)  & tax$Genus  == g
MG  <- c("Methanobacteriaceae","Methanomassiliicoccaceae","Methanoregulaceae","Methanocellaceae",
         "Methanosaetaceae","Methanomicrobiaceae","Methanosarcinaceae","Methanomethyliaceae",
         "Methanocorpusculaceae")

# group, tier, the ASVs it covers, and which gene(s) it is tested against
G <- list(
  list("Known, aerobic (pMMO / sMMO)",          "Known",                     st %in% "Known" & !fam("Methylomirabilaceae"), c("pmoA","mmoX")),
  list("NC10 (Methylomirabilaceae; anaerobic)", "Known",                     st %in% "Known" & fam("Methylomirabilaceae"),  c("pmoA","mmoX")),
  list("Beijerinckiaceae, putative",           "Putative",                  st %in% "Putative" & fam("Beijerinckiaceae"), c("pmoA","mmoX")),
  list("Methylacidiphilaceae, genus unresolved","Putative",                 fam("Methylacidiphilaceae") & unr(tax$Genus), c("pmoA","mmoX")),
  list("Lichenibacterium (1174-901-12)",       "Listed, not counted",       gen("1174-901-12"),                         c("pmoA","mmoX")),
  list("Roseiarcus",                           "Listed, not counted",       gen("Roseiarcus"),                          c("pmoA","mmoX")),
  list("Methylobacterium",                     "Listed, not counted",       gen("Methylobacterium-Methylorubrum"),      c("pmoA","mmoX")),
  list("Methanobacteriaceae",                  "Methanogen",                fam("Methanobacteriaceae"),                 "mcrA"),
  list("Methanomassiliicoccaceae",             "Methanogen",                fam("Methanomassiliicoccaceae"),            "mcrA"),
  list("Other methanogen families",            "Methanogen",                tax$Family %in% setdiff(MG, c("Methanobacteriaceae","Methanomassiliicoccaceae")), "mcrA"),
  list("Ammonia-oxidising archaea (Nitrososphaeraceae)", "Control",         fam("Nitrososphaeraceae"),                  c("pmoA","mmoX","mcrA")))
ab <- sapply(G, function(g) 100 * colSums(cnt[g[[3]], , drop = FALSE], na.rm = TRUE) / tot)
colnames(ab) <- sapply(G, `[[`, 1)
s16 <- data.frame(key = sub("[.]16S[.]S[0-9]*$", "", colnames(cnt)), ab, check.names = FALSE)

# ---- ddPCR: gene copies per bacterial 16S copy --------------------------------
o <- read.csv("data/raw/external/tree-microbiome/tree_data_methanogen_group.csv", check.names = FALSE)
o <- o[, names(o) != ""]
o$pmoA <- num(o$pmoa_loose); o$mmoX <- num(o$mmox_loose); o$mcrA <- num(o$mcra_probe_loose)
o$b16  <- num(o$X16S_per_ul); o$key <- paste0(o$seq_id, o$core_type)
m <- merge(o[, c("key","core_type","pmoA","mmoX","mcrA","b16")], s16, by = "key")
COMPS <- c(Inner = "Heartwood", Outer = "Sapwood", Mineral = "Mineral soil", Organic = "Organic soil")
m$comp <- COMPS[m$core_type]

# ---- tests --------------------------------------------------------------------
R <- do.call(rbind, lapply(G, function(g) do.call(rbind, lapply(g[[4]], function(gene)
  do.call(rbind, lapply(COMPS, function(k) {
    s  <- m[m$comp == k & is.finite(m[[gene]]) & is.finite(m$b16) & m$b16 > 0, ]
    x  <- s[[g[[1]]]]; y <- s[[gene]] / s$b16; np <- sum(x > 0)
    out <- data.frame(group = g[[1]], tier = g[[2]], gene = gene, compartment = k,
                      n = nrow(s), n_present = np, rho = NA_real_, p = NA_real_)
    if (np >= 5) { ct <- cor.test(x, y, method = "spearman", exact = FALSE)
                   out$rho <- unname(ct$estimate); out$p <- ct$p.value }
    out }))))))
R$q <- p.adjust(R$p, "BH")
R$sig <- with(R, ifelse(is.na(q), "", ifelse(q < .001, "***", ifelse(q < .01, "**",
               ifelse(q < .05, "*", ifelse(p < .05, "†", ""))))))
write.csv(R, out_path("16s_ddpcr_association.csv"), row.names = FALSE)

# ---- figure -------------------------------------------------------------------
TIERS <- c("Known","Putative","Listed, not counted","Methanogen","Control")
R$tier  <- factor(R$tier, levels = TIERS)
R$group <- factor(R$group, levels = rev(unique(sapply(G, `[[`, 1))))
R$gene  <- factor(R$gene, levels = c("pmoA","mmoX","mcrA"),
                  labels = c("italic(pmoA)","italic(mmoX)","italic(mcrA)"))
R$compartment <- factor(R$compartment, levels = COMPS)
R$tested <- !is.na(R$rho)
R$label  <- ifelse(R$tested, sprintf("%+.2f%s\n[%d]", R$rho, R$sig, R$n_present),
                   sprintf("n+ %d", R$n_present))
R$ink    <- ifelse(R$tested & abs(R$rho) > 0.33, "inverse", "primary")

p <- ggplot(R, aes(compartment, group)) +
  geom_tile(data = R[R$tested, ], aes(fill = rho), colour = "white", linewidth = 1) +
  geom_tile(data = R[!R$tested, ], fill = NA, colour = "#b9b8b3", linewidth = .35,
            linetype = "22", width = .9, height = .85) +
  geom_text(aes(label = label, colour = ink), size = 2.9, lineheight = .9) +
  scale_colour_manual(values = c(primary = "#1a1a19", inverse = "#ffffff"), guide = "none") +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0,
                       limits = c(-.5, .5), oob = scales::squish,
                       name = "Spearman ρ", breaks = c(-.5, -.25, 0, .25, .5)) +
  facet_grid(tier ~ gene, scales = "free", space = "free", labeller = labeller(gene = label_parsed)) +
  labs(x = NULL, y = NULL) +   # key (BH across all tests; dagger; dashed cells) is in the SI caption
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
        strip.text.x = element_text(size = 11),
        axis.text.x  = element_text(angle = 30, hjust = 1),
        legend.position = "bottom", legend.key.width = unit(1.4, "cm"),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(out_path("figSI_16s_ddpcr_association.png"), p, width = 9.5, height = 6.4, dpi = 300, bg = "white",
       device = ragg::agg_png)   # ragg: the default png device drops the rho and dagger glyphs
cat("Wrote figSI_16s_ddpcr_association.png and 16s_ddpcr_association.csv\n")
cat(sprintf("%d tests; %d significant after BH (q<.05); %d nominal only\n",
            sum(R$tested), sum(R$q < .05, na.rm = TRUE), sum(R$p < .05 & R$q >= .05, na.rm = TRUE)))
