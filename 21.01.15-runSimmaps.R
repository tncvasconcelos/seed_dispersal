# this script will generate the simmaps based on model averaged corHMM models generated for each dataset.

## functions
WtQ <- function(Q, Weights){
  Weights <- Weights[!is.na(Q)]/sum(Weights[!is.na(Q)])
  AvgQ <- sum(Q[!is.na(Q)] * Weights)
  return(AvgQ)
}

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

testMkToHMM <- function(file){
  load(file)
  phy <- obj[[1]]$phy
  dat <- obj[[1]]$data
  p.mk <- as.vector(na.omit(c(obj[[2]]$solution)))
  R1 <- getStateMat4Dat(dat, "ARD")$rate.mat
  rate.mat <- getFullMat(list(R1, R1), R1)
  p <- c(p.mk[1], p.mk[2], p.mk[1], p.mk[2], 1, 1)
  MK_LogLik <- obj[[2]]$loglik
  MK2HMM_LogLik <- corHMM(phy, dat, 2, rate.mat = rate.mat, p = p, node.states = "none")
  return(round(MK_LogLik, 5) == round(MK2HMM_LogLik$loglik, 5))
}

getSimmaps <- function(file, save.file, nMap){
  load(file)
  Q <- getModelAvgRate(file)
  phy <- obj[[1]]$phy
  dat <- obj[[1]]$data
  simmaps <- makeSimmap(phy, dat, Q, 2, nSim = nMap)
  save(simmaps, file = save.file)
  return(simmaps)
}

## imports
require(corHMM)

## run

# file organization
# wd <- "~/2021_SeedDispersal/"
wd <- getwd()
setwd(wd)
Rsaves <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))
labels <- unlist(lapply(strsplit(dir("res_corhmm/"), "-"), function(x) x[1]))
save.files <- paste0(wd, "/simmaps/", labels, "-Simmaps.R")

# test that the we can make a HMM structure with a MK model and get the same likelihood (this will allow us to model average MK and HMMs)
testMkToHMM(Rsaves[1])
testMkToHMM(Rsaves[2])
testMkToHMM(Rsaves[3])
testMkToHMM(Rsaves[4])

# model averaged corhmm runs
getModelAvgRate(Rsaves[1])
getModelAvgRate(Rsaves[2])
getModelAvgRate(Rsaves[3])
getModelAvgRate(Rsaves[4])

# run the simmaps
maps <- vector("list", length(Rsaves))
for(i in 1:length(Rsaves)){
  maps[[i]] <- getSimmaps(file = Rsaves[i], save.file = save.files[i], 1000)
}




