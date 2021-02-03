# this script will reconstruc the ancestral state based on our model averaged corHMM models

## functions 
WtQ <- function(Q, Weights){
  Weights <- Weights[!is.na(Q)]/sum(Weights[!is.na(Q)])
  AvgQ <- sum(Q[!is.na(Q)] * Weights)
  return(AvgQ)
}


# model average rates from corHMM
getModelAvgRate <- function(file){
  load(file)
  rate.mat <- obj$res_ARD.ARD.2$index.mat
  AICc <- unlist(lapply(obj, function(x) x$AICc))
  AICwt <- exp(-0.5 * AICc - min(AICc))/sum(exp(-0.5 * AICc - min(AICc)))
  # Solutions <- lapply(obj, function(x) x$solution)
  Solutions <- lapply(obj, function(x) c(na.omit(c(x$solution))))
  Solutions[[1]] <-  c(Solutions[[1]][1], NA, Solutions[[1]][2], NA,
                       NA, Solutions[[1]][1], NA, Solutions[[1]][2])
  Solutions[[2]] <-  c(Solutions[[2]][1], NA, Solutions[[2]][2], NA,
                       NA, Solutions[[2]][1], NA, Solutions[[2]][2])
  Rates <- do.call(rbind, Solutions)
  p.wt <- apply(Rates, 2, function(x) WtQ(x, AICwt))
  rate.mat[!is.na(rate.mat)] <- p.wt 
  return(rate.mat)
}


getTipRecon <- function(file){
  load(file)
  phy <- obj$res_ER$phy
  data <- obj$res_ER$data
  root.p <- obj$res_ARD$root.p
  index.mat <- obj$res_ARD.ARD.2$index.mat
  p <- getModelAvgRate(file)[sapply(1:max(index.mat, na.rm = TRUE), function(x) match(x, index.mat))]
  res <- corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = index.mat, node.states = "marginal", p = p, root.p = root.p, get.tip.states = TRUE)
  return(res)
}

## imports
require(corHMM)

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
files <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))

for(i in 1:length(files)){
  file <- files[i]
  file.name <- gsub("res_corhmm", "recon_corhmm", file)
  file.name <- gsub("corRes", "corRecon", file.name)
  res <- getTipRecon(file)
  save(res, file = file.name)
}








