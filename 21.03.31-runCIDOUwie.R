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
  dat.temp.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, temp = dat$mean_temp, se_temp = dat$within_sp_var_temp)
  dat.temp.se <- dat.temp.se[which(apply(dat.temp.se, 1, function(x) !any(is.na(x)))),]
  dat.temp <- data.frame(sp = dat$species, reg = dat$Fruit_type, temp = dat$mean_temp)
  dat.temp <- dat.temp[which(apply(dat.temp, 1, function(x) !any(is.na(x)))),]
  
  dat.prec.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, prec = dat$mean_prec, se_prec = dat$within_sp_var_prec)
  dat.prec.se <- dat.prec.se[which(apply(dat.prec.se, 1, function(x) !any(is.na(x)))),]
  dat.prec <- data.frame(sp = dat$species, reg = dat$Fruit_type, prec = dat$mean_prec)
  dat.prec <- dat.prec[which(apply(dat.prec, 1, function(x) !any(is.na(x)))),]
  
  dat.pet.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, pet = dat$mean_pet, se_pet = dat$within_sp_var_pet)
  dat.pet.se <- dat.pet.se[which(apply(dat.pet.se, 1, function(x) !any(is.na(x)))),]
  dat.pet <- data.frame(sp = dat$species, reg = dat$Fruit_type, pet = dat$mean_pet)
  dat.pet <- dat.pet[which(apply(dat.pet, 1, function(x) !any(is.na(x)))),]
  
  dat.arid.se <- data.frame(sp = dat$species, reg = dat$Fruit_type, arid = dat$mean_aridity, se_temp = dat$within_sp_var_aridity)
  dat.arid.se <- dat.arid.se[which(apply(dat.arid.se, 1, function(x) !any(is.na(x)))),]
  dat.arid <- data.frame(sp = dat$species, reg = dat$Fruit_type, arid = dat$mean_aridity)
  dat.arid <- dat.arid[which(apply(dat.arid, 1, function(x) !any(is.na(x)))),]
  return(list(dat.temp = dat.temp,
              dat.temp.se = dat.temp.se,
              dat.prec = dat.prec,
              dat.prec.se = dat.prec.se,
              dat.pet = dat.pet,
              dat.pet.se = dat.pet.se,
              dat.arid = dat.arid,
              dat.arid.se = dat.arid.se
  ))
}

organizeDat <- function(dat, simmap){
  mapping <- unlist(lapply(simmap$maps, function(x) names(x[length(x)])))
  nTip <- length(simmap$tip.label)
  TipStates <- mapping[match(match(dat$sp, simmap$tip.label), simmap$edge[,2])]
  dat$reg <- TipStates
  return(dat)
}

singleRunCID <- function(dat, simmap, model, mserr){
  for(i in 1:length(simmap$maps)){
    names(simmap$maps[[i]])[names(simmap$maps[[i]]) == "1"] <- "1"
    names(simmap$maps[[i]])[names(simmap$maps[[i]]) == "2"] <- "1"
    names(simmap$maps[[i]])[names(simmap$maps[[i]]) == "3"] <- "2"
    names(simmap$maps[[i]])[names(simmap$maps[[i]]) == "4"] <- "2"
  }
  simmap$mapped.edge[,1] <- simmap$mapped.edge[,1] + simmap$mapped.edge[,2]
  simmap$mapped.edge[,2] <- simmap$mapped.edge[,3] + simmap$mapped.edge[,4]
  simmap$mapped.edge <- simmap$mapped.edge[,-c(3,4)]
  data <- organizeDat(dat, simmap)
  obj <- try(OUwie(simmap, data, model, simmap.tree = TRUE, algorithm = "three.point", scaleHeight = FALSE, mserr = mserr))
  return(obj)
}

singleRunCD <- function(dat, simmap, model, mserr){
  for(i in 1:length(simmap$maps)){
    names(simmap$maps[[i]])[names(simmap$maps[[i]]) == "1"] <- "1"
    names(simmap$maps[[i]])[names(simmap$maps[[i]]) == "2"] <- "2"
    names(simmap$maps[[i]])[names(simmap$maps[[i]]) == "3"] <- "1"
    names(simmap$maps[[i]])[names(simmap$maps[[i]]) == "4"] <- "2"
  }
  simmap$mapped.edge[,1] <- simmap$mapped.edge[,1] + simmap$mapped.edge[,3]
  simmap$mapped.edge[,2] <- simmap$mapped.edge[,2] + simmap$mapped.edge[,4]
  simmap$mapped.edge <- simmap$mapped.edge[,-c(3,4)]
  data <- organizeDat(dat, simmap)
  obj <- try(OUwie(simmap, data, model, simmap.tree = TRUE, algorithm = "three.point", scaleHeight = FALSE, mserr = mserr))
  return(obj)
}

## imports
require(OUwie)
require(corHMM)
require(parallel)

wd <- getwd()
setwd(wd)
Rsaves <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))
labels <- unlist(lapply(strsplit(dir("res_corhmm/"), "-"), function(x) x[1]))
CSVs <- getCSVs(wd)[-5]
map_files <- paste0(wd, "/simmaps/", dir("simmaps/"))

ncores <- 42

# i in clade
for(i in 3){
  csv <- CSVs[i]
  name.clade <- gsub(".*trait_data/", "", csv)
  name.clade <- gsub("_niche.csv", "", name.clade)
  data <- getData(csv)
  completed.files <- dir(paste0("res_CID/", name.clade))
  # for j in each dataset (temp, prec, se or no, etc.)
  for(j in 1:length(data)){
    dat <- data[[j]]
    mserr <- ifelse(dim(dat)[2] == 3, "none", "known")
    name.dat <- names(data)[j]
    name.clade.dat <- paste0(name.clade, "-", name.dat, "-")
    if(length(grep(name.clade.dat, completed.files)) == 140){
      # we've already finished this dataset completely
      cat(name.clade.dat, "already complete...\n")
      next
    }
    simmaps_file <- map_files[grep(name.clade.dat, map_files)]
    # for each iteration of 100 simmaps (1000 total iterations)
    for(k in 1:length(simmaps_file)){
      load(simmaps_file[k])
      simmapnames <- names(simmaps)
      simmapnames <- gsub("res_", "", simmapnames)
      simmapnames <- gsub("\\.[0-9]*", "", simmapnames)
      simmaps <- simmaps[nchar(gsub("[0-9]", "", simmapnames)) >= 4]
      models <- c("BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
      # for l in each type of model being run (only hidden state models allowed)
      for(l in 1:length(models)){
        CIDChk <- length(grep(paste0(name.clade.dat, ".*", "-CID-", models[l], "-", k, ".Rsave"), completed.files))
        CDChk <- length(grep(paste0(name.clade.dat, ".*", "-CD-", models[l], "-", k, ".Rsave"), completed.files))
        if(CIDChk == 0){
          obj.CID <- mclapply(simmaps, function(x) singleRunCID(dat, x, models[l], mserr), mc.cores = ncores)
          file.name.CID <- paste0("/space_2/jamesboyko/2021_SeedDispersal/res_CID/", name.clade, "/", name.clade.dat, format(Sys.time(), "%y_%m_%d"), "-CID-", models[l], "-", k, ".Rsave")
          save(obj.CID, file = file.name.CID)
        }
        if(CDChk == 0){
          obj.CD <- mclapply(simmaps, function(x) singleRunCD(dat, x, models[l], mserr), mc.cores = ncores)
          file.name.CD <- paste0("/space_2/jamesboyko/2021_SeedDispersal/res_CID/", name.clade, "/", name.clade.dat, format(Sys.time(), "%y_%m_%d"), "-CD-", models[l], "-", k, ".Rsave")
          save(obj.CD, file = file.name.CD)
        }
      }
    }
  }
}


# test code
# data <- getData(csv)
# dat <- data$dat.arid.se
# # for(iter in 1:10){
# #   file.name <- paste0("/space_2/jamesboyko/2021_SeedDispersal/simmaps_pruned/Rosaceae-dat.arid.se-", format(Sys.time(), "%y_%m_%d"), "-simmap-", iter, ".Rsave")
# #   getSimmaps(Rsaves[4], dat = dat, save.file = file.name, nMap = 100)
# # }
simmaps_file <- dir("simmaps_pruned/", full.names = TRUE)
name.clade.dat <- "Rosaceae-dat.arid.se"

simmaps_file <- map_files[grep("Solanaceae-dat.arid.se", map_files)]
name.clade.dat <- "Solanaceae-dat.arid.se"



for(i in 1:length(simmaps_file)){
  load(simmaps_file[i])
  simmapnames <- names(simmaps)
  simmaps <- simmaps[nchar(gsub("[0-9]", "", simmapnames)) >= 9]
  models <- c("BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
  for(k in 1:length(models)){
    obj.CID <- mclapply(simmaps, function(x) singleRunCID(dat, x, models[k], mserr), mc.cores = ncores)
    obj.CD <- mclapply(simmaps, function(x) singleRunCD(dat, x, models[k], mserr), mc.cores = ncores)
    # save the modeling results of a dataset
    paste0("/space_2/jamesboyko/2021_SeedDispersal/res_CID/", name.clade.dat, "-CID.Rsave")
    save(obj.CID, file = paste0("/space_2/jamesboyko/2021_SeedDispersal/res_CID/", name.clade.dat, "-CID-", format(Sys.time(), "%y_%m_%d"), "-OURes-", models[k], "-", i, ".Rsave"))
    save(obj.CD, file = paste0("/space_2/jamesboyko/2021_SeedDispersal/res_CID/", name.clade.dat, "-CD-", format(Sys.time(), "%y_%m_%d"), "-OURes-", models[k], "-", i, ".Rsave"))
    obj <- NULL
  }
}


### evaluation

files <- dir()
models <- c("BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
iters <- paste0("-", 1:10, ".Rsave")
sola.files <- files[grep("Solanaceae", files)]
rosa.files <- files[grep("Rosaceae", files)]


data.frame_i <- c()
for(iter in 1:10){
  sola.files_i <- sola.files[grep(iters[iter], sola.files)]
  sola.files_i.CD <- sola.files_i[grep("CD", sola.files_i)]
  sola.files_i.CID <- sola.files_i[grep("CID", sola.files_i)]
  for(i in 1:length(sola.files_i.CD)){
    load(sola.files_i.CD[i])
    load(sola.files_i.CID[i])
    AIC.CD <- simplify2array(lapply(obj.CD, "[[", "AICc"))
    AIC.CID <- simplify2array(lapply(obj.CID, "[[", "AICc"))
    data.frame_CD <- data.frame(OUFit = models[i], Model = "CD", AICc = AIC.CD, dAICc = AIC.CD - AIC.CID)
    data.frame_CID <- data.frame(OUFit = models[i], Model = "CID", AICc = AIC.CID, dAICc = AIC.CID - AIC.CD)
    data.frame_i <- rbind(data.frame_i, data.frame_CD, data.frame_CID)
    cat("\r", paste0(round(i/length(sola.files_i.CD) * 100), "%"),"of iteration", iter, "complete...     ")
  }
  cat("\n")
}
sola.CID.fits <- data.frame_i
save(sola.CID.fits, file = "Solanaceae.CID.Fits.Rsave")



data.frame_i <- c()
for(iter in 1:10){
  rosa.files_i <- rosa.files[grep(iters[iter], rosa.files)]
  rosa.files_i.CD <- rosa.files_i[grep("CD", rosa.files_i)]
  rosa.files_i.CID <- rosa.files_i[grep("CID", rosa.files_i)]
  for(i in 1:length(rosa.files_i.CD)){
    load(rosa.files_i.CD[i])
    load(rosa.files_i.CID[i])
    index <- (unlist(lapply(obj.CD, class)) != "try-error") | (unlist(lapply(obj.CID, class)) != "try-error")
    obj.CD <- obj.CD[index]
    obj.CID <- obj.CID[index]
    AIC.CD <- simplify2array(lapply(obj.CD, "[[", "AICc"))
    AIC.CID <- simplify2array(lapply(obj.CID, "[[", "AICc"))
    data.frame_CD <- data.frame(OUFit = models[i], Model = "CD", AICc = AIC.CD, dAICc = AIC.CD - AIC.CID)
    data.frame_CID <- data.frame(OUFit = models[i], Model = "CID", AICc = AIC.CID, dAICc = AIC.CID - AIC.CD)
    data.frame_i <- rbind(data.frame_i, data.frame_CD, data.frame_CID)
    cat("\r", paste0(round(i/length(rosa.files_i.CD) * 100), "%"),"of iteration", iter, "complete...     ")
  }
  cat("\n")
}
rosa.CID.fits <- data.frame_i
save(rosa.CID.fits, file = "Rosaceae.CID.Fits.Rsave")



