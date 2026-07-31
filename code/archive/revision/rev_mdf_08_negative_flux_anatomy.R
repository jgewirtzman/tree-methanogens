#!/usr/bin/env Rscript
# ==============================================================================
# rev_mdf_08_negative_flux_anatomy.R      STEP 8
# ------------------------------------------------------------------------------
# ARE ANY OF THE NEGATIVE STEM FLUXES MEANINGFUL UPTAKE?
#
# Three detection treatments are carried side by side:
#   RAW        every measurement, no criterion
#   p < 0.05   fitted slope distinguishable from zero (statistically correct test)
#   MDF W95    above instrument resolution (conservative screen)
#
# For each: the full flux distribution, and then the negatives specifically -- how
# many, how large, and whether they concentrate in any species, height, campaign or
# season. A handful of tiny negatives scattered at random is noise. A coherent set of
# substantial negatives in a particular species or at a particular height would be a
# finding.
#
# BENCHMARKS for "meaningful"
#   our soil uptake, median            about -0.82 nmol m-2 s-1
#   our positive stem fluxes, median   reported below
#   instrument resolution (MDF W95)    reported below
# An uptake value that is orders of magnitude below the soil sink and barely above
# the detection limit is not a meaningful sink, whatever its p-value.
#
# Output: outputs/revision/negative_flux_anatomy.txt
#         outputs/revision/fig_negative_anatomy.png
# ==============================================================================

suppressMessages({library(dplyr); library(ggplot2); library(tidyr)})
outdir <- "outputs/revision"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rule <- function(s) cat("\n",strrep("=",78),"\n ",s,"\n",strrep("=",78),"\n",sep="")

ALL <- read.csv("outputs/revision/mdf_allan_all_campaigns.csv", stringsAsFactors=FALSE)
G <- bind_rows(
  read.csv("data/processed/flux/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="height_2021"),
  read.csv("data/processed/flux/CH4_best_flux_lgr_results_soil.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="semirigid_soil"),
  read.csv("deprecated/static-chambers/2023/CH4_best_flux_lgr_results.csv", stringsAsFactors=FALSE) %>%
    mutate(campaign="static_2023")) %>%
  filter(is.finite(best.flux), is.finite(flux.term), is.finite(nb.obs), nb.obs>2)

SIG <- ALL %>% mutate(campaign=ifelse(campaign=="untagged","semirigid_tree",campaign)) %>%
  group_by(campaign) %>% summarise(sigma=median(allan_sd_CH4,na.rm=TRUE), .groups="drop")
SIG <- bind_rows(SIG, data.frame(campaign="static_2023", sigma=1.298/0.86))
G <- G %>% left_join(SIG,"campaign") %>% filter(is.finite(sigma))
G$MDF95 <- G$sigma*qnorm(.975)/G$nb.obs*G$flux.term
G$is_soil <- G$campaign=="semirigid_soil"

# --- attach species / height where available ----------------------------------
ax <- read.csv("data/processed/flux/goflux_auxfile.csv", stringsAsFactors=FALSE)
G$height_cm <- ax$measurement_height[match(G$UniqueID, ax$UniqueID)]
G$species   <- ax$species[match(G$UniqueID, ax$UniqueID)]
a23 <- tryCatch(read.csv("data/processed/flux/ymf2023_goflux_auxfile.csv",
                         stringsAsFactors=FALSE, check.names=TRUE), error=function(e) NULL)
if (!is.null(a23) && "UniqueID" %in% names(a23)) {
  sp23 <- grep("Species.Code|Species_Code|species", names(a23), value=TRUE)[1]
  if (!is.na(sp23)) {
    i <- match(G$UniqueID, a23$UniqueID)
    G$species[is.na(G$species)] <- a23[[sp23]][i][is.na(G$species)]
  }
}
G$date <- as.Date(sub("^(\\d{8}).*$","\\1", G$UniqueID), format="%Y%m%d")
G$month <- as.integer(format(G$date, "%m"))

TR <- G[!G$is_soil,]
TREAT <- list(raw = rep(TRUE, nrow(TR)),
              `p<0.05` = TR$LM.p.val < 0.05,
              `MDF W95` = abs(TR$best.flux) > TR$MDF95)

# ==============================================================================
rule("1. FLUX DISTRIBUTION UNDER EACH TREATMENT (stems)")
# ==============================================================================
qs <- c(0,.01,.05,.25,.5,.75,.95,.99,1)
for (nm in names(TREAT)) {
  m <- TREAT[[nm]]; m[is.na(m)] <- FALSE; f <- TR$best.flux[m]
  cat(sprintf("\n  %-9s n=%4d   mean %7.4f   median %7.4f\n", nm, length(f), mean(f), median(f)))
  cat("    "); print(round(quantile(f, qs), 4))
}

# ==============================================================================
rule("2. THE NEGATIVES: how many, and how large?")
# ==============================================================================
cat(sprintf("\n  benchmarks:  soil uptake median %.3f | stem MDF W95 median %.4f | stem positive median %.4f\n",
    median(G$best.flux[G$is_soil]), median(TR$MDF95), median(TR$best.flux[TR$best.flux>0])))
cat(sprintf("\n  %-9s %8s %8s %10s %10s %10s %12s\n",
    "treatment","n neg","% neg","median neg","most neg","|neg|/MDF","vs soil sink"))
for (nm in names(TREAT)) {
  m <- TREAT[[nm]]; m[is.na(m)] <- FALSE
  n <- TR[m & TR$best.flux < 0, ]
  if (!nrow(n)) { cat(sprintf("  %-9s %8d\n", nm, 0)); next }
  cat(sprintf("  %-9s %8d %7.1f%% %10.4f %10.4f %10.2f %11.1f%%\n", nm, nrow(n),
      100*nrow(n)/sum(m), median(n$best.flux), min(n$best.flux),
      median(abs(n$best.flux)/n$MDF95),
      100*abs(median(n$best.flux))/abs(median(G$best.flux[G$is_soil]))))
}
cat("
  The last column is the median negative stem flux as a percentage of the median soil
  uptake. A stem 'sink' two orders of magnitude weaker than the soil sink, sitting a
  hair above its own detection limit, is not a meaningful sink.
")

# ==============================================================================
rule("3. DO THE NEGATIVES CONCENTRATE ANYWHERE?")
# ==============================================================================
for (nm in c("p<0.05","MDF W95")) {
  m <- TREAT[[nm]]; m[is.na(m)] <- FALSE
  d <- TR[m, ]; d$neg <- d$best.flux < 0
  cat(sprintf("\n  ---- %s (n=%d, %d negative) ----\n", nm, nrow(d), sum(d$neg)))
  cat("\n  BY CAMPAIGN:\n")
  print(d %>% group_by(campaign) %>% summarise(n=n(), n_neg=sum(neg),
    pct=round(100*mean(neg),1), med_neg=round(median(best.flux[neg]),4), .groups="drop") %>%
    as.data.frame(), row.names=FALSE)
  if (any(!is.na(d$height_cm))) {
    cat("\n  BY MEASUREMENT HEIGHT (height campaign only):\n")
    print(d %>% filter(!is.na(height_cm)) %>% group_by(height_cm) %>%
      summarise(n=n(), n_neg=sum(neg), pct=round(100*mean(neg),1), .groups="drop") %>%
      as.data.frame(), row.names=FALSE)
  }
  if (any(!is.na(d$species))) {
    ss <- d %>% filter(!is.na(species)) %>% group_by(species) %>%
      summarise(n=n(), n_neg=sum(neg), pct=round(100*mean(neg),1), .groups="drop") %>%
      filter(n>=8) %>% arrange(desc(pct))
    cat("\n  BY SPECIES (n>=8):\n"); print(as.data.frame(ss), row.names=FALSE)
  }
  cat("\n  BY MONTH:\n")
  print(d %>% filter(!is.na(month)) %>% group_by(month) %>%
    summarise(n=n(), n_neg=sum(neg), pct=round(100*mean(neg),1), .groups="drop") %>%
    as.data.frame(), row.names=FALSE)
}

# ==============================================================================
rule("4. THE ACTUAL NEGATIVE VALUES AT MDF W95")
# ==============================================================================
m <- TREAT[["MDF W95"]]; m[is.na(m)] <- FALSE
n <- TR[m & TR$best.flux<0, ] %>% arrange(best.flux)
cat(sprintf("\n  All %d significantly-negative stem measurements:\n\n", nrow(n)))
cat(sprintf("  %-30s %10s %9s %8s %7s %6s %s\n","UniqueID","flux","MDF95","ratio","r2","ht","species"))
for (i in seq_len(nrow(n))) with(n[i,],
  cat(sprintf("  %-30s %10.4f %9.4f %8.2f %7.3f %6s %s\n", substr(UniqueID,1,30), best.flux,
      MDF95, abs(best.flux)/MDF95, ifelse(is.na(LM.r2),NA,LM.r2),
      ifelse(is.na(height_cm),"-",height_cm), ifelse(is.na(species),"-",species))))

p <- bind_rows(lapply(names(TREAT), function(nm) {
  mm <- TREAT[[nm]]; mm[is.na(mm)] <- FALSE
  data.frame(treatment=nm, flux=TR$best.flux[mm]) })) %>%
  mutate(treatment=factor(treatment, levels=names(TREAT)))
gp <- ggplot(p, aes(flux, treatment, fill=treatment)) +
  geom_violin(scale="width", linewidth=.2, alpha=.75) +
  geom_vline(xintercept=0, linetype="22", colour="grey30") +
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=.005, base=10)) +
  labs(title="Stem CH4 flux distribution under three detection treatments",
       subtitle="pseudo-log axis; dashed line = zero",
       x=expression(CH[4]~flux~(nmol~m^-2~s^-1)), y=NULL) +
  theme_bw(base_size=9) + theme(legend.position="none")
ggsave(file.path(outdir,"fig_negative_anatomy.png"), gp, width=8, height=4, dpi=200, bg="white")
cat("\n  Written: outputs/revision/fig_negative_anatomy.png\n")
