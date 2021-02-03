# this script will generate the simmaps based on model averaged corHMM models generated for each dataset.
# then this script will run the OU models on the simmas generated in runSimmaps


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

getSimmaps <- function(file, dat, save.file, nMap){
  load(file)
  Q <- getModelAvgRate(file)
  phy <- obj[[1]]$phy
  phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% dat$sp])
  simmaps <- makeSimmap(phy, dat[,c(1,2)], Q, 2, nSim = nMap)
  save(simmaps, file = save.file)
  return(simmaps)
}

getCSVs <- function(wd){
  CSVs <- dir(paste0(wd, "/trait_data/"))
  CSVs <- CSVs[grep("niche", CSVs)]
  CSVs <- paste0(wd, "/trait_data/", CSVs)
  return(CSVs)
}

getData <- function(csv){
  dat <- read.csv(csv)
  dat.temp.se <- data.frame(sp = dat$species, reg = dat$Dispersal_mode, temp = dat$temp, se_temp = dat$se_temp)
  dat.temp.se <- dat.temp.se[which(apply(dat.temp.se, 1, function(x) !any(is.na(x)))),]
  dat.prec.se <- data.frame(sp = dat$species, reg = dat$Dispersal_mode, prec = dat$prec, se_prec = dat$se_prec)
  dat.prec.se <- dat.prec.se[which(apply(dat.prec.se, 1, function(x) !any(is.na(x)))),]
  dat.temp <- data.frame(sp = dat$species, reg = dat$Dispersal_mode, temp = dat$temp)
  dat.temp <- dat.temp[which(apply(dat.temp, 1, function(x) !any(is.na(x)))),]
  dat.prec <- data.frame(sp = dat$species, reg = dat$Dispersal_mode, prec = dat$prec)
  dat.prec <- dat.prec[which(apply(dat.prec, 1, function(x) !any(is.na(x)))),]
  return(list(dat.temp.se = dat.temp.se,
              dat.prec.se = dat.prec.se,
              dat.temp = dat.temp,
              dat.prec = dat.prec))
}

organizeDat <- function(dat, simmap){
  mapping <- unlist(lapply(simmap$maps, function(x) names(x[length(x)])))
  nTip <- length(simmap$tip.label)
  TipStates <- mapping[match(match(dat$sp, simmap$tip.label), simmap$edge[,2])]
  dat$reg <- TipStates
  return(dat)
}

singleRun <- function(dat, simmap, model, mserr){
  data <- organizeDat(dat, simmap)
  obj <- OUwie(simmap, data, model, simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr)
  return(obj)
}

## imports
require(OUwie)
require(corHMM)
require(parallel)

## run

# file organization
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
Rsaves <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))
labels <- unlist(lapply(strsplit(dir("res_corhmm/"), "-"), function(x) x[1]))
CSVs <- getCSVs(wd)

# # test that the we can make a HMM structure with a MK model and get the same likelihood (this will allow us to model average MK and HMMs)
# testMkToHMM(Rsaves[1])
# testMkToHMM(Rsaves[2])
# testMkToHMM(Rsaves[3])
# testMkToHMM(Rsaves[4])
# 
# # model averaged corhmm runs
# getModelAvgRate(Rsaves[1])
# getModelAvgRate(Rsaves[2])
# getModelAvgRate(Rsaves[3])
# getModelAvgRate(Rsaves[4])

# input params 
ncores <- 44
nmap <- 100
iter <- 2
# i = j = 1
# k = 6
# run the simmaps
for(i in 1:length(CSVs)){
  csv <- CSVs[i]
  data <- getData(csv)
  # for j in each trait dataset
  for(j in 1:length(data)){
    file.name <- paste0(wd, "/simmaps/", labels[i], "-", names(data)[j], "-", format(Sys.time(), "%y_%m_%d"), "-simmap-", iter, ".Rsave")
    simmaps <- getSimmaps(file = Rsaves[i], dat = data[[j]], save.file = file.name, nMap = nmap)
    if(dim(data[[j]])[2] == 4){
      mserr = "known"
    }else{
      mserr = "none"
    }
    models <- c("BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
    for(k in 1:length(models)){
      obj <- mclapply(simmaps, function(x) singleRun(data[[j]], x, models[k], mserr), mc.cores = ncores)
      # save the modeling results of a dataset
      file.name <- paste0(wd, "/res_ouwie/", labels[i], "/", labels[i], "-", names(data)[j], "-", format(Sys.time(), "%y_%m_%d"), "-OURes-", models[k], "-", iter, ".Rsave")
      save(obj, file = file.name)
      obj <- NULL
    }
  }
}



data <- organizeDat(data[[j]], simmaps[[1]])
tmp <- data[,c(1,2,4)]
obj <- obj[[1]]

testA <- OUwie.boot(phy=simmaps[[1]], data=data, model=obj$model, nboot=2, alpha=obj$solution[1,], sigma.sq=obj$solution[2,], theta=obj$solution[3,], theta0=obj$theta[1,1], simmap.tree=FALSE, scaleHeight=TRUE, mserr="known", algorithm="three.point")

testB <- OUwie.boot(phy=simmaps[[1]], data=data, model=obj$model, nboot=2, alpha=obj$solution[1,], sigma.sq=obj$solution[2,], theta=obj$solution[3,], theta0=obj$theta[1,1], simmap.tree=TRUE, scaleHeight=TRUE, mserr="known", algorithm="three.point")

sim.data<-OUwie.sim(simmaps[[1]],tmp,simmap.tree=TRUE,scaleHeight=TRUE,alpha=obj$solution[1,], sigma.sq=obj$solution[2,], theta=obj$solution[3,], theta0=obj$theta[1,1], mserr="known")



