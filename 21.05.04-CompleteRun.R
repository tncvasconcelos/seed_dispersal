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

# organizeDat <- function(dat, simmap){
#   mapping <- unlist(lapply(simmap$maps, function(x) names(x[length(x)])))
#   nTip <- length(simmap$tip.label)
#   TipStates <- mapping[match(match(rownames(dat), simmap$tip.label), simmap$edge[,2])]
#   dat$reg <- TipStates
#   return(dat)
# }

singleRun <- function(dat, simmap, model, mserr){
  data <- organizeDat(dat, simmap)
  obj <- try(OUwie(simmap, data, model, simmap.tree = TRUE, algorithm = "three.point", scaleHeight = FALSE, mserr = mserr))
  return(obj)
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

ComeTogether <- function(name.clade, name.dat, wd){
  cat("Starting", name.clade, "for", name.dat, "...\n")
  folder.path.a <- paste0(wd, "/res_CID/", name.clade)
  Rsaves.a <- dir(folder.path.a, full.names = TRUE)
  Rsaves.a <- Rsaves.a[grep(paste0("-dat.", name.dat, "-"), Rsaves.a)]
  folder.path.b <- paste0(wd, "/res_ouwie/", name.clade)
  Rsaves.b <- dir(folder.path.b, full.names = TRUE)
  Rsaves.b <- Rsaves.b[grep(paste0("-dat.", name.dat, "-"), Rsaves.b)]
  map.files <- dir(paste0(wd, "/simmaps/"))
  map.files <- map.files[grep(paste0(name.clade, "-dat.", name.dat, "-"), map.files)]
  
  Rsaves.CID <- Rsaves.a[grep("-CID-", Rsaves.a)]
  Rsaves.CD <- Rsaves.a[grep("-CD-", Rsaves.a)]
  Rsaves.HYB <- Rsaves.b
  
  for(iter in 1:10){
    Rsaves.CID_i <- Rsaves.CID[grep(paste0("-", iter, ".Rsave"), Rsaves.CID)]
    Rsaves.CD_i <- Rsaves.CD[grep(paste0("-", iter, ".Rsave"), Rsaves.CD)]
    Rsaves.HYB_i <- Rsaves.HYB[grep(paste0("-", iter, ".Rsave"), Rsaves.HYB)]
    models <- c("-BM1-", "-BMS-", "-OU1-", "-OUM-", "-OUMA-", "-OUMV-", "-OUMVA-")
    obj_HYB_iter <- obj_CID_iter <- obj_CD_iter <- vector("list", length(models))
    names(obj_HYB_iter) <- names(obj_CID_iter) <- names(obj_CD_iter) <- gsub("-", "", models)
    count <- 1
    cat("Loading models for iteration", iter, "of", "10.\n")
    for(model in models){
      load(Rsaves.HYB_i[grep(model, Rsaves.HYB_i)])  #obj
      obj_HYB_iter[[count]] <- obj
      load(Rsaves.CID_i[grep(model, Rsaves.CID_i)])  #obj.CID
      obj_CID_iter[[count]] <- obj.CID
      load(Rsaves.CD_i[grep(model, Rsaves.CD_i)])  #obj.CD
      obj_CD_iter[[count]] <- obj.CD
      count <- count + 1
    }
    for(i in 1:100){
      name.map <- names(obj_HYB_iter[[1]])[i]
      HYB_map_i <- lapply(obj_HYB_iter, function(x) x[[i]])
      if(name.map %in% names(obj_CID_iter[[1]])){
        CID_map_i <- lapply(obj_CID_iter, function(x) x[[grep(paste0(name.map, "$"), names(obj_CID_iter[[1]]))]])
        CD_map_i <- lapply(obj_CD_iter, function(x) x[[grep(paste0(name.map, "$"), names(obj_CID_iter[[1]]))]])
        out <- list(HYB = HYB_map_i,
                    CID = CID_map_i,
                    CD = CD_map_i)
      }else{
        out <- list(HYB = HYB_map_i)
      }
      file.name <- paste0(name.clade, "-", name.dat, "-", name.map, "-", iter, ".Rsave")
      save(out, file = paste0(wd, "/res_bymap/", name.clade, "/", file.name))
    }
  }
}

### functions to examine the parameter and modeling results

getAIC <- function(obj){
  AIC <- unlist(lapply(obj, function(x) lapply(x, "[[", "AIC")))
  dAIC <- AIC - min(AIC)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  return(data.frame(AIC = AIC, dAIC=dAIC, AICwt=AICwt))
}

getSummaryTables <- function(name.clade, name.dat, wd){
  folder.path <- paste0(wd, "/res_bymap/", name.clade)
  Rsaves <- dir(folder.path, full.names = TRUE)
  Rsaves <- Rsaves[grep(paste0("-", name.dat, "-"), Rsaves)]
  AICTables <- vector("list", length(Rsaves))
  for(i in 1:length(Rsaves)){
    cat("\r", round(i/length(Rsaves) * 100), "%:", Rsaves[i], "    ")
    load(Rsaves[i])
    if(any(unlist(lapply(out, function(x) lapply(x, class))) == "try-error") | 
       any(unlist(lapply(out, function(x) lapply(x, is.null))))){
      AICTables[[i]] <- NA
    }else{
      AICTables[[i]] <- getAIC(out)
      nMaxRegime <- max(unlist(lapply(out, function(x) lapply(x, "[[", "param.count"))))
      ParNames <- paste0(rep(c("alpha_", "sigma.sq_", "theta_"), nMaxRegime/3), rep(1:(nMaxRegime/3), each = 3))
      ParTable_i <- matrix(data = NA, nrow = dim(AICTables[[i]])[1], ncol = length(ParNames), 
                           dimnames = list(rownames(AICTables[[i]]), ParNames))
      count <- 1
      for(j in 1:length(out)){
        for(k in 1:length(out[[j]])){
          ParTable_i[count, ] <- getParVector(out[[j]][[k]]$solution, ParNames)
          count <- count + 1
        }
      }
      AICTables[[i]] <- cbind(AICTables[[i]], ParTable_i)
    }
  }
  cat("\nSaving Tables. \n")
  file.name <- paste0(wd, "/res_tables/", name.clade, "-", name.dat, ".Rsave")
  save(AICTables, file = file.name)
}

getParVector <- function(solution, par.vec){
  out.vector <- rep(NA, length(par.vec))
  names(out.vector) <- par.vec
  solution.par <- paste0(rownames(solution), "_", rep(colnames(solution), each = 3))
  solution.tmp <- as.vector(solution)
  names(solution.tmp) <- solution.par
  out.vector[match(names(solution.tmp), names(out.vector))] <- solution.tmp
  return(out.vector)
}

getSummaryTablesNames <- function(name.clade, name.dat, wd){
  folder.path <- paste0(wd, "/res_bymap/", name.clade)
  Rsaves <- dir(folder.path, full.names = TRUE)
  Rsaves <- Rsaves[grep(paste0("-", name.dat, "-"), Rsaves)]
  AICTableNames <- vector("character", length(Rsaves))
  for(i in 1:length(Rsaves)){
    AICTableNames[i] <- lapply(strsplit(Rsaves[i], "/"), function(x) x[length(x)])[[1]]
  }
  cat("\nSaving Table Names. \n")
  file.name <- paste0(wd, "/res_tables/names-", name.clade, "-", name.dat, ".Rsave")
  save(AICTableNames, file = file.name)
}










## imports
require(OUwie)
require(corHMM)
require(parallel)

## run
# file organization
# wd <- "~/2021_SeedDispersal"
wd <- getwd()
setwd(wd)
Rsaves <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))
labels <- unlist(lapply(strsplit(dir("res_corhmm/"), "-"), function(x) x[1]))
CSVs <- getCSVs(wd)
CSVs <- CSVs[-4]

# for now skipping non-measurement error
# bolu4 is running Mela - i = 3 (90 cores)
# bolu2 is running Eric - i = 2 (42 cores)
# bolu1 is running Apoc - i = 1 (84 cores)
# bolu2 is running Sola - i = 5 (42 cores)
# bolu1 is running Rosa - i = 4 (84 cores)

ncores <- 60
nmap <- 100
i <- 3

for(iter in 8:10){
  csv <- CSVs[i]
  data <- getData(csv)
  # for j in each trait dataset
  for(j in 1:length(data)){
    if(dim(data[[j]])[2] == 4){
      mserr = "known"
    }else{
      # mserr = "none"
      next # change this to do msser 
    }
    print(names(data)[j])
    file.name <- paste0(wd, "/simmaps/", labels[i], "-", names(data)[j], "-", format(Sys.time(), "%y_%m_%d"), "-simmap-", iter, ".Rsave")
    simmaps <- getSimmaps(file = Rsaves[i], dat = data[[j]], save.file = file.name, nMap = nmap)
    simmapnames <- names(simmaps)
    simmapnames <- gsub("res_", "", simmapnames)
    simmapnames <- gsub("\\.[0-9]*", "", simmapnames)
    simmaps.pruned <- simmaps[nchar(gsub("[0-9]", "", simmapnames)) >= 4]
    models <- c("BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
    for(k in 1:length(models)){
      obj <- mclapply(simmaps, function(x) singleRun(data[[j]], x, models[k], mserr), mc.cores = ncores)
      # save the modeling results of a dataset
      file.name <- paste0(wd, "/res_ouwie/", labels[i], "/", labels[i], "-", names(data)[j], "-", format(Sys.time(), "%y_%m_%d"), "-OURes-", models[k], "-", iter, ".Rsave")
      save(obj, file = file.name)
      obj.CID <- mclapply(simmaps.pruned, function(x) singleRunCID(data[[j]], x, models[k], mserr), mc.cores = ncores)
      file.name.CID <- paste0(wd, "/res_CID/", labels[i], "/", labels[i], "-", names(data)[j], "-", format(Sys.time(), "%y_%m_%d"), "-CID-", models[k], "-", iter, ".Rsave")
      save(obj.CID, file = file.name.CID)
      obj.CD <- mclapply(simmaps.pruned, function(x) singleRunCD(data[[j]], x, models[k], mserr), mc.cores = ncores)
      file.name.CD <- paste0(wd, "/res_CID/", labels[i], "/", labels[i], "-", names(data)[j], "-", format(Sys.time(), "%y_%m_%d"), "-CD-", models[k], "-", iter, ".Rsave")
      save(obj.CD, file = file.name.CD)
    }
  }
}

# Ericaceae completed on bolu2 and being summarized
# Apocynaceae completed on bolu1 and being summarized
# Rosaceae completed on bolu1 and being summarized
# Solanaceae completed on bolu2 and being summarized
# Melastomataceae completed on bolu4 and summarized

require(parallel)
dat.types <- c("temp.se", "prec.se", "pet.se", "arid.se")
wd <- getwd()
# bring everything into a by map result
mclapply(dat.types, function(x) ComeTogether(name.clade = "Melastomataceae", name.dat = x, wd = wd), mc.cores = 4)

# summarize the per map resutls
dat.types <- c("temp.se", "prec.se", "pet.se", "arid.se")
wd <- getwd()
mclapply(dat.types, function(x) getSummaryTables(name.clade = "Melastomataceae", name.dat = x, wd = wd), mc.cores = 4)
mclapply(dat.types, function(x) getSummaryTablesNames(name.clade = "Melastomataceae", name.dat = x, wd = wd), mc.cores = 4)








# evaluation script for the res_tables for each dataset and set of simmaps. mainly responisble for tip averaging.
# generates tip averages found in tip average tables folder
# this is done locally
organizeAICTable <- function(AICTable){
  if(class(AICTable) == "logical"){
    return(NA)
  }
  # get the alpha and sigma parameters
  alpha_index <- grep("alpha", colnames(AICTable))
  sigma.sq_index <- grep("sigma.sq", colnames(AICTable))
  theta_index <- grep("theta", colnames(AICTable))
  
  ## organize HYB
  HYBTable <- AICTable[grep("HYB", rownames(AICTable)),]
  # organize alpha
  HYBTable[1,alpha_index] <- 1e-10 
  HYBTable[2,alpha_index] <- 1e-10 
  HYBTable[3,alpha_index] <- rep(HYBTable[3,alpha_index[1]], length(alpha_index))
  # organize sigma.sq
  HYBTable[1,sigma.sq_index] <- rep(HYBTable[1,sigma.sq_index[1]], length(sigma.sq_index))
  HYBTable[3,sigma.sq_index] <- rep(HYBTable[3,sigma.sq_index[1]], length(sigma.sq_index))
  # organize theta
  HYBTable[1,theta_index] <- rep(HYBTable[1,theta_index[1]], length(theta_index))
  HYBTable[3,theta_index] <- rep(HYBTable[3,theta_index[1]], length(theta_index))
  
  if(dim(AICTable)[1] == 21){
    ## organize CID
    CIDTable <- AICTable[grep("CID", rownames(AICTable)),]
    # organize alpha
    CIDTable[1,alpha_index] <- 1e-10 
    CIDTable[2,alpha_index] <- 1e-10 
    CIDTable[3,alpha_index] <- rep(CIDTable[3,alpha_index[1]], length(alpha_index))
    # organize sigma.sq
    CIDTable[1,sigma.sq_index] <- rep(CIDTable[1,sigma.sq_index[1]], length(sigma.sq_index))
    CIDTable[3,sigma.sq_index] <- rep(CIDTable[3,sigma.sq_index[1]], length(sigma.sq_index))
    # organize theta
    CIDTable[1,theta_index] <- rep(CIDTable[1,theta_index[1]], length(theta_index))
    CIDTable[3,theta_index] <- rep(CIDTable[3,theta_index[1]], length(theta_index))
    CIDTable.cp <- CIDTable
    CIDTable[,alpha_index] <- CIDTable.cp[,alpha_index[c(1,1,2,2)]] # state 1 reprersentst state A (i.e. 1 and 2)
    CIDTable[,sigma.sq_index] <- CIDTable.cp[,sigma.sq_index[c(1,1,2,2)]] # state 2 represetnts state B (i.e. 3 and 4)
    CIDTable[,theta_index] <- CIDTable.cp[,theta_index[c(1,1,2,2)]] # so we take value 1 and put it in 1,2 and 2 goes 3,4
    ## organize CD
    CDTable <- AICTable[grep("CD", rownames(AICTable)),]
    # organize alpha
    CDTable[1,alpha_index] <- 1e-10 
    CDTable[2,alpha_index] <- 1e-10 
    CDTable[3,alpha_index] <- rep(CDTable[3,alpha_index[1]], length(alpha_index))
    # organize sigma.sq
    CDTable[1,sigma.sq_index] <- rep(CDTable[1,sigma.sq_index[1]], length(sigma.sq_index))
    CDTable[3,sigma.sq_index] <- rep(CDTable[3,sigma.sq_index[1]], length(sigma.sq_index))
    # organize theta
    CDTable[1,theta_index] <- rep(CDTable[1,theta_index[1]], length(theta_index))
    CDTable[3,theta_index] <- rep(CDTable[3,theta_index[1]], length(theta_index))
    CDTable.cp <- CDTable
    CDTable[,alpha_index] <- CDTable.cp[,alpha_index[c(1,2,1,2)]] # state 1 is sstate 1A and 1B
    CDTable[,sigma.sq_index] <- CDTable.cp[,sigma.sq_index[c(1,2,1,2)]] # state 2 is 2A and 2B
    CDTable[,theta_index] <- CDTable.cp[,theta_index[c(1,2,1,2)]] # 1 goes into 1 and 3, 2 goes into 2 and 4
    out.table <- rbind(HYBTable, CIDTable, CDTable)
  }else{
    out.table <- HYBTable
  }
  
  return(out.table)
}

getAvgParTablePerMap <- function(RsaveResult, RsaveMapNames, corObject){
  load(RsaveResult)
  load(RsaveMapNames)
  organizedAICTables <- lapply(AICTables, organizeAICTable)
  TipReconTables <- vector("list", length = length(organizedAICTables))
  for(i in 1:length(organizedAICTables)){
    if(class(organizedAICTables[[i]]) == "logical"){
      next
    }
    if(any(organizedAICTables[[i]]$AIC < -100000)){
      next
    }
    map_name_i <- strsplit(AICTableNames[i], "-")[[1]][length(strsplit(AICTableNames[i], "-")[[1]])-1]
    CharNum <- nchar(gsub("[0-9]", "", map_name_i))
    MatchedNames <- unlist(sapply(names(corObject), function(x) grep(paste0(x, "+"), map_name_i)))
    Simmap_name_i <- names(MatchedNames[ifelse(CharNum > 8, 2, 1)])
    corObject_i <- corObject[[match(Simmap_name_i, names(corObject))]]
    AvgParPerMap <- colSums(organizedAICTables[[i]][,4:dim(organizedAICTables[[i]])[2]] * organizedAICTables[[i]]$AICwt)
    AvgParPerMap <- matrix(AvgParPerMap, 3, dim(corObject_i$tip.states)[2], dimnames = list(c("alpha", "sigma.sq", "theta")))
    TipAverageTable_i <- matrix(NA, dim(corObject_i$tip.states)[1], 3, dimnames = list(rownames(corObject_i$tip.states), c("alpha", "sigma.sq", "theta")))
    for(j in 1:dim(corObject_i$tip.states)[1]){
      tip_recon_j <- corObject_i$tip.states[j,]
      sp_j_pars <- colSums(t(AvgParPerMap) * tip_recon_j)
      TipAverageTable_i[j,] <- sp_j_pars
    }
    TipReconTables[[i]] <- TipAverageTable_i
  }
  TipReconTables <- TipReconTables[unlist(lapply(TipReconTables, function(x) !is.null(x)))]
  # combining all the different map results into 1 andremoving the tails of the species distributions 
  AvgTipReconTable <- matrix(NA, dim(corObject_i$tip.states)[1], 3, dimnames = list(rownames(corObject_i$tip.states), c("alpha", "sigma.sq", "theta")))
  for(i in 1:dim(AvgTipReconTable)[1]){
    sp_i_table <- do.call(rbind, lapply(TipReconTables, function(x) x[i,]))
    no.to.remove <- round(0.05 * dim(sp_i_table)[1])
    alpha_i <- mean(sort(sp_i_table[,1])[no.to.remove:(dim(sp_i_table)[1]-no.to.remove)])
    sigma.sq_i <- mean(sort(sp_i_table[,2])[no.to.remove:(dim(sp_i_table)[1]-no.to.remove)])
    theta_i <- mean(sort(sp_i_table[,3])[no.to.remove:(dim(sp_i_table)[1]-no.to.remove)])
    AvgTipReconTable[i,] <- c(alpha_i, sigma.sq_i, theta_i)
  }
  return(AvgTipReconTable)
}

getAvgTipReconForRsave <- function(RsaveResult, Rsaves_names, cor_folder){
  RsaveFileName <- unlist(strsplit(RsaveResult, "/"))[length(unlist(strsplit(RsaveResult, "/")))]
  name.clade.dat <- gsub(".Rsave", "", RsaveFileName)
  RsaveMapNames <- Rsaves_names[grep(paste0(name.clade.dat, ".Rsave"), Rsaves_names)]
  name.clade <- gsub("-.*", "", name.clade.dat)
  cor_file <- cor_folder[grep(name.clade, cor_folder)]
  load(cor_file)
  AvgTipReconTbl <- getAvgParTablePerMap(RsaveResult, RsaveMapNames, obj)
  return(AvgTipReconTbl)
}

# actually running things
Rsaves <- dir("~/2021_SeedDispersal/res_tables/", full.names = TRUE)
Rsaves_names <- Rsaves[grep("names-", Rsaves)]
Rsaves_res <- Rsaves[-grep("names-", Rsaves)]
cor_folder <- dir("~/2021_SeedDispersal/res_corhmm/", full.names = TRUE)

for(i in 1:length(Rsaves_res)){
  print(Rsaves_res[i])
  tmp <- getAvgTipReconForRsave(Rsaves_res[i], Rsaves_names, cor_folder)
  FileName <- strsplit(Rsaves_res[i], "/")[[1]][length(strsplit(Rsaves_res[i], "/")[[1]])]
  FileName <- gsub(".Rsave", "-TipAvgTable.csv", FileName)
  FileName <- paste0("~/2021_SeedDispersal/tip_avg_tables/", FileName)
  write.csv(tmp, file = FileName)
}


