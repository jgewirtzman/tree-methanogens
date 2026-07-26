suppressMessages({library(ranger); library(dplyr); library(lme4)}); set.seed(42)
load("outputs/models/TRAINING_DATA.RData")
X0 <- as.data.frame(X_tree); d <- tree_train_complete
y <- asinh(d$stem_flux_corrected); k <- is.finite(y)
X0<-X0[k,]; d<-d[k,]; y<-y[k]; sp<-d$species; gn<-sub(" .*","",sp)
env <- intersect(c("air_temp_C_mean","soil_temp_C_mean","soil_moisture_at_tree","dbh_within_z","height_cm"), names(X0))
r2 <- function(a,b) 1-sum((a-b)^2)/sum((a-mean(a))^2)
TR <- read.csv("/Users/jongewirtzman/My Drive/Research/YMF Tree Microbiomes & Methane/tree-gas-traits/data/processed/ymf_with_traits.csv", stringsAsFactors=FALSE)
tc <- intersect(c("wood_density_gcm3","bark_density_gcm3","porosity_num","wood_heartwood_pH","try_rooting_depth","try_plant_longevity","gymnosperm"), names(TR))
spt <- TR %>% group_by(spcode) %>% summarise(across(all_of(tc), ~mean(.x,na.rm=TRUE)), .groups="drop")
ep <- strsplit(sp, " ")
code <- toupper(vapply(ep, function(z) if(length(z)>=2) paste0(substr(z[1],1,2), substr(z[2],1,2)) else NA_character_, character(1)))
TB <- spt[match(code, spt$spcode), tc]
TB <- as.data.frame(lapply(TB, function(v){v[is.na(v)] <- median(v, na.rm=TRUE); as.numeric(v)}))
cat("=== LEAVE-ONE-SPECIES-OUT: predicting a species never measured ===\n\n")
big <- names(which(table(sp) >= 30))
out <- bind_rows(lapply(c("none","current","hier","traits"), function(cfg){
  pr <- rep(NA_real_, length(y))
  for (s in big) {
    te <- which(sp==s); tr <- which(sp!=s)
    Xb <- X0[, setdiff(names(X0),"taxon_prior_asinh"), drop=FALSE]
    if (cfg %in% c("current","hier")) {
      rf0 <- ranger(x=X0[tr,env,drop=FALSE], y=y[tr], num.trees=200, min.node.size=5, num.threads=1, seed=42)
      r <- y[tr]-rf0$predictions
      if (cfg=="current") {
        z <- split(r, gn[tr]); z <- z[sapply(z,length)>=5]; mg <- sapply(z, median, na.rm=TRUE)
        p <- unname(ifelse(!is.na(match(gn,names(mg))), mg[match(gn,names(mg))], 0))
      } else {
        m <- tryCatch(suppressMessages(lmer(r ~ 1 + (1|gn), data=data.frame(r=r, gn=gn[tr]))), error=function(e) NULL)
        b <- if (is.null(m)) c() else setNames(ranef(m)$gn[,1], rownames(ranef(m)$gn))
        p <- if (length(b)) unname(ifelse(!is.na(match(gn,names(b))), b[match(gn,names(b))], 0)) else rep(0,length(gn))
      }
      Xb$taxon_prior_asinh <- as.numeric(ifelse(is.na(p),0,p))
    }
    if (cfg=="traits") Xb <- cbind(Xb, TB)
    m <- ranger(x=Xb[tr,,drop=FALSE], y=y[tr], num.trees=800, min.node.size=5,
                mtry=max(1,floor(sqrt(ncol(Xb)))), num.threads=1, seed=42)
    pr[te] <- predict(m, Xb[te,,drop=FALSE])$predictions
  }
  i <- !is.na(pr)
  data.frame(prior=cfg, n=sum(i), r2=r2(y[i],pr[i]), rho=cor(y[i],pr[i],method="spearman"),
             ratio=mean(sinh(pr[i]))/mean(sinh(y[i])))
}))
print(out %>% arrange(desc(r2)), row.names=FALSE, digits=3)
cat("\n  species held out:", paste(big, collapse=", "), "\n")
