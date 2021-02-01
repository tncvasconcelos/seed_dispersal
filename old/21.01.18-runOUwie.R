# this script will run the OU models on the simmas generated in runSimmaps

## functions
# gets all the CSV files we want to analyze 
getCSVs <- function(wd){
  CSVs <- dir(paste0(wd, "/trait_data/"))
  CSVs <- CSVs[grep("niche", CSVs)]
  CSVs <- paste0(wd, "/trait_data/", CSVs)
  return(CSVs)
}

# organizes the data from the csvs into a list of all the datasets we want to analyze
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

# this function will match data to simmaps and simmaps to data. it will also change the regimes of the data to match the tips from the simmaps 
organizeDat <- function(dat, simmaps){
  simmap <- simmaps[[1]]
  Phy2Drop <- simmap$tip.label[!simmap$tip.label %in% dat$sp]
  simmaps.pruned <- lapply(simmaps, function(x) drop.tip.simmap(x, Phy2Drop))
  simmap <- simmaps.pruned[[1]]
  dat <- dat[dat$sp %in% simmap$tip.label,]
  mapping <- unlist(lapply(simmap$maps, function(x) names(x[length(x)])))
  nTip <- length(simmap$tip.label)
  TipStates <- mapping[match(match(dat$sp, simmap$tip.label), simmap$edge[,2])]
  dat$reg <- TipStates
  return(list(dat = dat,
              simmaps = simmaps.pruned))
}

## imports
require(OUwie)
require(corHMM)
require(phytools)
require(parallel)


simmap <- simmaps[[1]]
mapping <- unlist(lapply(simmap$maps, function(x) names(x[length(x)])))
nTip <- length(simmap$tip.label)
TipStates <- mapping[match(match(simmap$tip.label, simmap$tip.label), simmap$edge[,2])]
dat <- data.frame(sp = simmap$tip.label, 
                  reg = TipStates,
                  dat = rnorm(length(TipStates), 10, 2))

OUwie(simmap, dat, "BMS", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = "none")

## run
wd <- "~/2021_SeedDispersal"
setwd(wd)
# wd <- getwd()

CSVs <- getCSVs(wd)
Simmaps <- paste0(wd, "/simmaps/", dir("simmaps/"))
labels <- unique(unlist(lapply(strsplit(sort(dir("trait_data/")), "_"), function(x) x[1])))
ncores <- 1
i = j = 1
# for i in each clade
for(i in 1:length(CSVs)){
  csv <- CSVs[i]
  load(Simmaps[i])
  simmaps <- simmaps[1:4] # temporary for testing purposes
  data <- getData(csv)
  # for j in each trait dataset
  for(j in 1:length(data)){
    data_j <- organizeDat(dat = data[[j]], simmap = simmaps)
    if(dim(data_j$dat)[2] == 4){
      mserr = "known"
    }else{
      mserr = "none"
    }
    # BM1
    OUwie(data_j[[2]][[1]], data_j[[1]], "BMS", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr)
    
    BM1 <- mclapply(data_j[[2]], function(x) OUwie(x, data_j[[1]], "BM1", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr), mc.cores = ncores)
    # BMS
    BMS <- mclapply(data_j[[2]], function(x) OUwie(x, data_j[[1]], "BMS", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr), mc.cores = ncores)
    # OU1
    OU1 <- mclapply(data_j[[2]], function(x) OUwie(x, data_j[[1]], "OU1", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr), mc.cores = ncores)
    # OUM
    OUM <- mclapply(data_j[[2]], function(x) OUwie(x, data_j[[1]], "OUM", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr), mc.cores = ncores)
    # OUMA
    OUMA <- mclapply(data_j[[2]], function(x) OUwie(x, data_j[[1]], "OUMA", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr), mc.cores = ncores)
    # OUMV
    OUMV <- mclapply(data_j[[2]], function(x) OUwie(x, data_j[[1]], "OUMV", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr), mc.cores = ncores)
    # OUMVA
    OUMVA <- mclapply(data_j[[2]], function(x) OUwie(x, data_j[[1]], "OUMVA", simmap.tree = TRUE, algorithm = "three.point", scaleHeight = TRUE, mserr = mserr), mc.cores = ncores)
    obj_j <- list(BM1, BMS, OU1, OUM, OUMA, OUMV, OUMVA)
    names(obj_j) <- paste0(c("BM1_", "BMS_", "OU1_", "OUM_", "OUMA_", "OUMV_", "OUMVA_"), colnames(data[[j]])[3])
    # save the modeling results of a dataset
    file.name <- paste0(labels[i], "-", colnames(data[[j]])[3], "-", format(Sys.time(), "%y_%m_%d"), "-OURes.Rsave")
    save(obj_j, file = file.name)
  }
}


# csv -> trait data -> mod with simap -> run OUwie



