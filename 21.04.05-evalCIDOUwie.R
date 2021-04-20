## functiosn to bring all modeling results together per one simmap
name.clade <- "Apocynaceae"
name.dat <- "arid"

ComeTogether <- function(name.clade, name.dat){
  cat("Starting", name.clade, "for", name.dat, "...\n")
  folder.path.a <- paste0("/space_2/jamesboyko/2021_SeedDispersal/res_CID/", name.clade)
  Rsaves.a <- dir(folder.path.a, full.names = TRUE)
  Rsaves.a <- Rsaves.a[grep(paste0("-dat.", name.dat, "-"), Rsaves.a)]
  folder.path.b <- paste0("/space_2/jamesboyko/2021_SeedDispersal/res_ouwie/", name.clade)
  Rsaves.b <- dir(folder.path.b, full.names = TRUE)
  Rsaves.b <- Rsaves.b[grep(paste0("-dat.", name.dat, "-"), Rsaves.b)]
  map.files <- dir("/space_2/jamesboyko/2021_SeedDispersal/simmaps/")
  map.files <- map.files[grep(paste0(name.clade, "-dat.", name.dat, "-"), map.files)]
  
  Rsaves.CID <- Rsaves.a[grep("-CID-", Rsaves.a)]
  Rsaves.CD <- Rsaves.a[grep("-CD-", Rsaves.a)]
  Rsaves.HYB <- Rsaves.b
  
  for(iter in 1:10){
    # the iteration for the map file related to HYB models are correct, but the iteration of the CD and CID files are based (incorrectly) on the order in whcih a simmap appears in the map.files. i.e. something named iter=10 in the simmap file, could actually correspond to iter=7 in the Rsaves HYB
    CID_Iter <- paste0(grep(paste0(iter, ".Rsave"), map.files), ".Rsave")
    Rsaves.CID_i <- Rsaves.CID[grep(CID_Iter, Rsaves.CID)]
    Rsaves.CD_i <- Rsaves.CD[grep(CID_Iter, Rsaves.CD)]
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
      save(out, file = paste0("/space_2/jamesboyko/2021_SeedDispersal/res_bymap/", name.clade, "/", file.name))
    }
  }
}


clades <- c("Melastomataceae")
dat.types <- c("temp", "prec", "pet", "arid", "temp.se", "prec.se", "pet.se", "arid.se")

for(clade in clades){
  for(dat.type in dat.types){
    ComeTogether(name.clade = clade, name.dat = dat.type)
  }
}

### functions to examine the parameter and modeling results

getAIC <- function(obj){
  AIC <- unlist(lapply(obj, function(x) lapply(x, "[[", "AIC")))
  dAIC <- AIC - min(AIC)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  return(data.frame(AIC = AIC, dAIC=dAIC, AICwt=AICwt))
}

getSummaryTables <- function(name.clade, name.dat){
  folder.path <- paste0("/space_2/jamesboyko/2021_SeedDispersal/res_bymap/", name.clade)
  Rsaves <- dir(folder.path, full.names = TRUE)
  Rsaves <- Rsaves[grep(paste0("-", name.dat, "-"), Rsaves)]
  AICTables <- vector("list", length(Rsaves))
  for(i in 1:length(Rsaves)){
    cat("\r", round(i/length(Rsaves) * 100), "%:", Rsaves[i], "    ")
    load(Rsaves[i])
    if(any(unlist(lapply(out, function(x) lapply(x, class))) == "try-error")){
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
  file.name <- paste0("/space_2/jamesboyko/2021_SeedDispersal/res_tables/", name.clade, "-", name.dat, ".Rsave")
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



clades <- c("Apocynaceae", "Ericaceae", "Melastomataceae", "Rosaceae", "Solanaceae")
dat.types <- c("temp", "prec", "pet", "arid", "temp.se", "prec.se", "pet.se", "arid.se")

for(clade in clades){
  for(dat.type in dat.types){
    getSummaryTables(name.clade = clade, name.dat = dat.type)
  }
}







load("~/2021_SeedDispersal/res_CID/Rosaceae.CID.Fits.Rsave")
load("~/2021_SeedDispersal/res_CID/Solanaceae.CID.Fits.Rsave")

load("~/2021_SeedDispersal/res_tables/Rosaceae-arid-TRUE.Rsave")
p1 <- ggplot(rosa.CID.fits, aes(x=OUFit, y=dAICc, fill=Model)) + 
      geom_boxplot()

p2 <- ggplot(rosa.CID.fits, aes(x=OUFit, y=AICc, fill=Model)) + 
  geom_boxplot()

grid.arrange(p1, p2, nrow=1)

load("~/2021_SeedDispersal/res_tables/Solanaceae-arid-TRUE.Rsave")
p1 <- ggplot(sola.CID.fits, aes(x=OUFit, y=dAICc, fill=Model)) + 
  geom_boxplot()

p2 <- ggplot(sola.CID.fits, aes(x=OUFit, y=AICc, fill=Model)) + 
  geom_boxplot()

minAICc <- unlist(lapply(ResTables, function(x) min(x[,1])))
OtherFits <- data.frame(OUFit = "???", Model = "O.G.", AICc = minAICc, dAICc = 0)
sola.CID.fits <- rbind(sola.CID.fits, OtherFits)
p3 <- ggplot(sola.CID.fits, aes(x=OUFit, y=AICc, fill=Model)) + 
  geom_boxplot()


grid.arrange(p1, p2, p3, nrow=1)

CD.rosa <- rosa.CID.fits[rosa.CID.fits$Model == "CD", ]
CID.rosa <- rosa.CID.fits[rosa.CID.fits$Model == "CID", ]

AICcTable <- data.frame(CD = CD.rosa$AICc, CID = CID.rosa$AICc)
colMeans(AICcTable)