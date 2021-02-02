# this script will load, summarize, and analyze the results of 1000 OUwie runs on our simmaps 

## functions
# a function which will return a list of tables that have the model parameters and AIC wt 
getSummaryTable <- function(file, round = FALSE){
  load(file)
  obj <- obj_j
  nMap <- length(obj[[1]])
  nRegime <- length(obj[[7]][[1]]$tot.states)
  ResTables <- vector("list", nMap)
  # get the relavent params and stats for each map i and model in the set
  for(i in 1:nMap){
    obj_i <- lapply(obj, function(x) x[[i]])
    ResTable <- matrix(0, length(obj_i), (nRegime * 3) + 2)
    rownames(ResTable) <- unlist(lapply(strsplit(names(obj_i), "_"), function(x) x[1]))
    colnames(ResTable) <- c("AICc", "AICcWt", 
                            paste0("Alpha_", 1:nRegime), paste0("Sigma_", 1:nRegime), paste0("Optim_", 1:nRegime))
    AICcs <- unlist(lapply(obj_i, function(x) x$AICc))
    AICcWt <- exp(-0.5 * AICcs - min(AICcs))/sum(exp(-0.5 * AICcs - min(AICcs)))
    Solutions <- lapply(obj_i, function(x) x$solution)
    BM1 <- c(rep(Solutions[[1]][1], nRegime), c(rep(Solutions[[1]][2], nRegime)), c(rep(Solutions[[1]][3], nRegime)))
    BMS <- c(t(Solutions[[2]]))
    OU1 <- c(rep(Solutions[[3]][1], nRegime), c(rep(Solutions[[3]][2], nRegime)), c(rep(Solutions[[3]][3], nRegime)))
    OUM <- c(t(Solutions[[4]]))
    OUMA <- c(t(Solutions[[5]]))
    OUMV <- c(t(Solutions[[6]]))
    OUMVA <- c(t(Solutions[[7]]))
    ResTable[,1] <- AICcs
    ResTable[,2] <- AICcWt
    ResTable[,3:((nRegime * 3) + 2)] <- rbind(BM1, BMS, OU1, OUM, OUMA, OUMV, OUMVA)
    if(round == TRUE){
      ResTables[[i]] <- round(ResTable, 2)
    }else{
      ResTables[[i]] <- ResTable
    }
  }
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
    AICcWt <- exp(-0.5 * AICcs - min(AICcs))/sum(exp(-0.5 * AICcs - min(AICcs)))
    obj[i-2] <- sum(Param[isParam] * AICcWt)
  }
  return(obj)
}


## imports
require(corHMM)
require(OUwie)

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
Rsaves <- paste0(wd, "/res_ouwie/", dir("res_ouwie/"))
labels <- unlist(lapply(strsplit(dir("res_ouwie/"), "-"), function(x) paste0(x[1], x[2])))

file <- Rsaves[1]


SummaryTables <- getSummaryTable(file, FALSE)
AvgParams <- colMeans(do.call(rbind, lapply(SummaryTables, summarizeTable)))
AvgAICcWt <- colMeans(do.call(rbind, lapply(SummaryTables, function(x) x[,2])))

round(AvgParams, 3)
round(AvgAICcWt, 3)




