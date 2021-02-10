# this script will generate the simmaps based on model averaged corHMM models generated for each dataset.
# then this script will run the OU models on the simmas generated in runSimmaps


## functions
getSimmaps <- function(file, dat, save.file, nMap){
  load(file)
  AICc <- unlist(lapply(obj, function(x) x$AICc))
  AICc <- (AICc + min(abs(AICc)))
  dAICc <- AICc - min(AICc)
  AICwt <- exp(-0.5 * dAICc)/sum(exp(-0.5 * dAICc))
  MapSample <- sample(x = names(AICwt), size = nMap, replace = TRUE, prob = AICwt)
  IndexMaps <- sapply(names(AICwt), function(x) length(which(MapSample %in% x)))
  phy <- obj[[1]]$phy
  phy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% dat$sp])
  simmaps.list <- vector("list", length(IndexMaps))
  for(i in 1:length(IndexMaps)){
    if(IndexMaps[i] == 0) next
    Q <- obj[[i]]$solution
    if(dim(Q)[1] == 4){
      rate.cat = 2
    }else{
      rate.cat = 1
    }
    simmaps.list[[i]] <- makeSimmap(phy, dat[,c(1,2)], Q, rate.cat, nSim = IndexMaps[i])
  }
  names(simmaps.list) <- names(IndexMaps)
  simmaps <- do.call(c, simmaps.list)
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
  dat.temp.se <- data.frame(sp = dat$species, reg = dat$Dispersal_mode, temp = dat$temp, se_temp = dat$within_sp_var_temp)
  dat.temp.se <- dat.temp.se[which(apply(dat.temp.se, 1, function(x) !any(is.na(x)))),]
  dat.prec.se <- data.frame(sp = dat$species, reg = dat$Dispersal_mode, prec = dat$prec, se_prec = dat$within_sp_var_prec)
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
  obj <- try(OUwie(simmap, data, model, simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr, ub = 10))
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

# input params 
ncores <- 40
nmap <- 100
iter <- 2
# i = j = 1
# k = 1
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



