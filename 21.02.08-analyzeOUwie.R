# this script will load, summarize, and analyze the results of 1000 OUwie runs on our simmaps 

## functions
# a function which will return a list of tables that have the model parameters and AIC wt 
getSummaryTable <- function(folder, dat.type, se){
  Rsaves <- dir(folder)
  if(se == TRUE){
    Rsaves <- Rsaves[grep("se-", Rsaves)]
  }else{
    Rsaves <- Rsaves[-grep("se-", Rsaves)]
  }
  files <- paste0(folder, "/", Rsaves[grep(dat.type, Rsaves)])
  models <- c("-BM1-", "-BMS-", "-OU1-", "-OUM-", "-OUMA-", "-OUMV-", "-OUMVA-")
  ModelList <- vector("list", length(models))
  names(ModelList) <- models
  count <- 1
  cat("Loading Models...\n")

  for(i in models){
    ToLoad <- files[grep(i, files)]
    iters <- unlist(lapply(strsplit(ToLoad, "-"), function(x) x[length(x)]))
    iters <- sort(iters)
    obj_j <- list()
    for(j in 1:length(iters)){
      load(ToLoad[grep(iters[j], ToLoad)])
      obj_j <- c(obj_j, obj)
    }
    ModelList[[count]] <- obj_j
    count <- count + 1
  }
  
  nMap <- length(ModelList[[1]])
  corModels <- names(ModelList[[1]])
  ErrMat <- matrix(unlist(lapply(ModelList, function(x) lapply(x, function(y) class(y) == "try-error"))), nMap, 7)
  LikMat <- matrix(unlist(lapply(ModelList, function(x) lapply(x, function(y) try(abs(y$loglik), silent = TRUE) > 1e5))), nMap, length(models))
  ErrMat <- ErrMat | LikMat
  ErrIndex <- which(!apply(ErrMat, 1, any))
  cat(nMap - length(ErrIndex), "models resulted in try-errors.\n")
  ResTables <- list()
  count <- 1
  # get the relavent params and stats for each map i and model in the set
  for(i in ErrIndex){
    obj_i <- lapply(ModelList, function(x) x[[i]])
    dat.check <- do.call(cbind, lapply(obj_i, function(x) x$phy$maps))
    if(!all(apply(dat.check, 1, function(x) length(unique(x)) == 1))){
      cat("\nNot all the data matches, check", i)
    }
    lapply(ModelList, function(x) x[[i]])
    nRegime <- length(obj_i[[7]]$tot.states)
    ResTable <- matrix(0, length(obj_i), (nRegime * 3) + 2)
    rownames(ResTable) <- unlist(lapply(strsplit(names(obj_i), "_"), function(x) x[1]))
    colnames(ResTable) <- c("AICc", "AICcWt", 
                            paste0("Alpha_", 1:nRegime), paste0("Sigma_", 1:nRegime), paste0("Optim_", 1:nRegime))
    AICcs <- unlist(lapply(obj_i, function(x) x$AICc))
    AICcs <- AICcs + abs(min(AICcs))
    dAICcs <- AICcs - min(AICcs)
    AICcWt <- exp(-0.5 * dAICcs)/sum(exp(-0.5 * dAICcs))
    Solutions <- lapply(obj_i, function(x) x$solution)
    BM1 <- c(rep(Solutions[[1]][1], nRegime), c(rep(Solutions[[1]][2], nRegime)), c(rep(Solutions[[1]][3], nRegime)))
    BMS <- c(t(Solutions[[2]]))
    OU1 <- c(rep(Solutions[[3]][1], nRegime), c(rep(Solutions[[3]][2], nRegime)), c(rep(Solutions[[3]][3], nRegime)))
    OUM <- c(t(Solutions[[4]]))
    OUMA <- c(t(Solutions[[5]]))
    OUMV <- c(t(Solutions[[6]]))
    OUMVA <- c(t(Solutions[[7]]))
    ResTable[,1] <- unlist(lapply(obj_i, function(x) x$AICc))
    ResTable[,2] <- AICcWt
    ResTable[,3:((nRegime * 3) + 2)] <- rbind(BM1, BMS, OU1, OUM, OUMA, OUMV, OUMVA)
    ResTables[[count]] <- ResTable
    count <- count + 1
  }
  names(ResTables) <- corModels[ErrIndex]
  file.name <- folder
  clade <- lapply(strsplit(folder, "/"), function(x) x[length(x)])[[1]]
  replacer <- paste0(clade, "-", dat.type, "-", se, ".Rsave")
  file.name <- gsub("res_ouwie", "res_tables", file.name)
  file.name <- gsub(clade, replacer, file.name)
  cat("\nSaving ResTables...", file.name)
  save(ResTables, file = file.name)
  return(ResTables)
}

# model averages a single simmap given the simmap table containing pars and AICc 
summarizeTable <- function(table){
  obj <- vector("numeric", dim(table)[2]-2)
  names(obj) <- colnames(table)[3:dim(table)[2]]
  for(i in 3:dim(table)[2]){
    Param <- table[,i]
    isParam <- !is.na(table[,i])
    AICcs <- table[,1][isParam]
    AICcs <- AICcs + abs(min(AICcs))
    dAICcs <- AICcs - min(AICcs)
    AICcWt <- exp(-0.5 * dAICcs)/sum(exp(-0.5 * dAICcs))
    obj[i-2] <- sum(Param[isParam] * AICcWt)
  }
  return(obj)
}

# produces a table where the states, hidden states, and value are given per user request
getFigureTable <- function(SumTable, param){
  ParamTable <- do.call(rbind, lapply(SumTable, summarizeTable))
  SubsetTable <- ParamTable[,grep(param, colnames(ParamTable))]
  n <- dim(SubsetTable)[1]
  FigureTable <- data.frame(
    ObsState = rep(rep(c("Abiotic", "Biotic"), each = n), 2),
    HidState = rep(c("A", "B"), each = n*2),
    Val = as.vector(SubsetTable))
  return(FigureTable)
}

# quick tip averging function
averageTipStatePars <- function(pars, tip.states, biota){
  Alpha <- colSums(t(tip.states) * pars[grep("Alpha", names(pars))])
  Sigma <- colSums(t(tip.states) * pars[grep("Sigma", names(pars))])
  Optim <- colSums(t(tip.states) * pars[grep("Optim", names(pars))])
  out <- data.frame("ObsState" = biota, Alpha = Alpha, Sigma = Sigma, Optim = Optim, row.names = rownames(tip.states))
  return(out)
}


# get model averaged tip rates
getTipRates <- function(cor_file, AvgParams){
  load(cor_file)
  corNames <- names(obj)
  AICcs <- unlist(lapply(obj, function(x) x$AICc))
  AICcs <- AICcs + abs(min(AICcs))
  dAICcs <- AICcs - min(AICcs)
  AICcWt <- exp(-0.5 * dAICcs)/sum(exp(-0.5 * dAICcs))
  Index <- vector("numeric", length = length(AvgParams))
  for(i in 1:length(corNames)){
    ind_i <- grepl(corNames[i], names(AvgParams)) & (nchar(names(AvgParams)) < nchar(corNames[i])+3)
    Index[ind_i] <- i
  }
  tmp <- vector("list", length(AvgParams))
  for(i in 1:length(AvgParams)){
    tip.states_i <- obj[[Index[i]]]$tip.states
    biota_i <- obj[[Index[i]]]$data[,2]
    pars_i <- AvgParams[[i]]
    tmp[[i]] <- averageTipStatePars(pars_i, tip.states_i, biota_i)
  }

  Alpha <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[,2])))
  Sigma <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[,3])))
  Optim <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[,4])))
  out <- data.frame(ObsSt = tmp[[1]][,1], Alpha = Alpha, Sigma = Sigma, Optim = Optim, row.names = rownames(tmp[[1]]))
  return(out)
}


## imports
require(corHMM)
require(OUwie)

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
Folders <- paste0(wd, "/res_ouwie/", dir("res_ouwie/"))
cor_folder <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))
dat.types <- c("temp", "prec", "pet", "arid")
params <- c("Alpha", "Sigma", "Optim")

# make plots for a given clade and standard error
se <- FALSE
for(se in c(TRUE, FALSE)){
  for(i in 1:length(dat.types)){
    dat.type <- dat.types[i]
    # file.name <- paste0(wd, "/figures/", dat.type, "-SE.", se, ".pdf")
    file.name <- paste0(wd, "/tables/", dat.type, "-SE.", se, ".csv")
    print(file.name)
    res <- vector("list", length(Folders))
    for(k in 1:length(Folders)){
      ou_folder <- Folders[k]
      Clade <- strsplit(ou_folder, "/")[[1]][length(strsplit(ou_folder, "/")[[1]])]
      cat("Starting", Clade, "...\n")
      cor_file <- cor_folder[grep(Clade, cor_folder)]
      # extract the results of the OU models
      SumTable <- getSummaryTable(ou_folder, dat.type, se)
      # model average the OU models
      AvgParams <- lapply(SumTable, summarizeTable)
      # do tip averaging to get it in terms of biotic and abiotic
      tab <- getTipRates(cor_file, AvgParams)
      AvgAICcWt <- colMeans(do.call(rbind, lapply(SumTable, function(x) x[,2])))
      res[[k]] <- cbind(Clade, tab)
    }
    res_table <- do.call(rbind, res)  
    write.csv(res_table, file.name)
  }
}



## imports
require(corHMM)
require(OUwie)

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
clades <- unlist(lapply(strsplit(dir("res_corhmm/"), "-"), function(x) x[1]))
cor_folder <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))
table_folder <- paste0(wd, "/res_tables/", dir("res_tables/"))
dat.types <- c("temp", "prec", "pet", "arid")
params <- c("Alpha", "Sigma", "Optim")

for(se in c(TRUE, FALSE)){
  for(i in 1:length(dat.types)){
    dat.type <- dat.types[i]
    # file.name <- paste0(wd, "/figures/", dat.type, "-SE.", se, ".pdf")
    file.name <- paste0(wd, "/tables/", dat.type, "-SE.", se, ".csv")
    res <- vector("list", length(clades))
    for(k in 1:length(clades)){
      clade_i <- clades[k]
      cat("\nStarting", clade_i, "...\n")
      print(file.name)
      cor_file <- cor_folder[grep(clade_i, cor_folder)]
      load(cor_file)
      # extract the results of the OU models
      sumfile <- table_folder[grep(clade_i, table_folder)]
      sumfile <- sumfile[grep(dat.type, sumfile)]
      sumfile <- sumfile[grep(se, sumfile)]
      load(sumfile)
      SumTable <- ResTables
      # model average the OU models
      AvgParams <- lapply(SumTable, summarizeTable)
      # do tip averaging to get it in terms of biotic and abiotic
      tab <- getTipRates(cor_file, AvgParams)
      AvgAICcWt <- colMeans(do.call(rbind, lapply(SumTable, function(x) x[,2])))
      res[[k]] <- cbind(clade_i, tab)
    }
    res_table <- do.call(rbind, res)  
    write.csv(res_table, file.name)
  }
}


out <- vector("numeric", length(table_folder))
name <- unlist(lapply(strsplit(table_folder, "/"), function(x) x[length(x)]))
names(out) <- name
for(i in 1:length(table_folder)){
  load(table_folder[i])
  out[i] <- length(ResTables)
}

