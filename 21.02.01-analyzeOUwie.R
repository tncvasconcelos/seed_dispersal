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
  ModelList <- vector("list", length(files))
  names(ModelList) <- models
  count <- 1
  cat("Loading Models...\n")
  for(i in models){
    ToLoad <- files[grep(i, files)]
    load(ToLoad)
    cat("\r", ToLoad, "loaded...                  ")
    ModelList[[count]] <- obj
    count <- count + 1
  }
  nMap <- length(ModelList[[1]])
  ErrMat <- matrix(unlist(lapply(ModelList, function(x) lapply(x, function(y) class(y) == "try-error"))), nMap, 7)
  ErrIndex <- which(!apply(ErrMat, 1, any))
  nRegime <- length(ModelList[[7]][[ErrIndex[1]]]$tot.states)
  ResTables <- vector("list", length(ErrIndex))
  count <- 1
  # get the relavent params and stats for each map i and model in the set
  for(i in ErrIndex){
    obj_i <- lapply(ModelList, function(x) x[[i]])
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
    ResTable[,1] <- AICcs
    ResTable[,2] <- AICcWt
    ResTable[,3:((nRegime * 3) + 2)] <- rbind(BM1, BMS, OU1, OUM, OUMA, OUMV, OUMVA)
    ResTables[[count]] <- ResTable
    count <- count + 1
  }
  cat("\nDone.\n")
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

# get model averaged tip rates
getTipRates <- function(cor_file, AvgParams){
  pars <- matrix(AvgParams, 4, 3, dimnames = list(c("1A", "2A", "1B", "2B"), c("Alpha", "Sigma", "Optim")))
  out <- matrix(NA, dim(res$tip.states)[1], 4, dimnames = list(rownames(res$tip.states), c("ObsSt", "Alpha", "Sigma", "Optim")))
  out <- as.data.frame(out)
  for(i in 1:dim(res$tip.states)[1]){
    sp_i <- rownames(res$tip.states)[i]
    out[i,1] <- res$data[i, 2]
    out[i,2] <- sum(res$tip.states[i,] * pars[,1])
    out[i,3] <- sum(res$tip.states[i,] * pars[,2])
    out[i,4] <- sum(res$tip.states[i,] * pars[,3])
  }
  return(out)
}


## imports
require(corHMM)
require(OUwie)
require(viridis)
require(ggplot2)
require(gridExtra)

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
Folders <- paste0(wd, "/res_ouwie/", dir("res_ouwie/"))
cor_folder <- paste0(wd, "/recon_corhmm/", dir("recon_corhmm/"))
dat.types <- c("temp", "prec")
params <- c("Alpha", "Sigma", "Optim")

# param <- params[3]
ou_folder <- Folders[3]
clade <- strsplit(ou_folder, "/")[[1]][length(strsplit(ou_folder, "/")[[1]])]
cor_file <- cor_folder[grep(clade, cor_folder)]
se <- FALSE

# make plots for a given clade and standard error
count <- 1
nPlots <- length(dat.types) * length(params)
plots <- vector("list", nPlots)

for(j in 1:length(dat.types)){
  dat.type <- dat.types[j]
  SumTable <- getSummaryTable(ou_folder, dat.type, se)
  AvgParams <- colMeans(do.call(rbind, lapply(SumTable, summarizeTable)))
  AvgAICcWt <- colMeans(do.call(rbind, lapply(SumTable, function(x) x[,2])))
  tab <- getTipRates(cor_file, AvgParams)
  for(i in 1:3){
    tab_i<- tab[,c(1, i+1)]
    colnames(tab_i)[2] <- "Val"
    param <- colnames(tab)[i+1]
    cols <- viridis(2)
    plots[[count]] <- 
      ggplot(tab_i, aes(x=ObsSt, y=Val)) + 
      labs(x = dat.type, y = param) +
      theme(axis.text = element_text(size = 11), legend.justification=c(0,0), legend.position=c(0,0.85)) +
      scale_fill_manual(values=cols) + 
      theme(text = element_text(size = 20)) + 
      ggtitle((clade)) + 
      geom_boxplot()
    count <- count + 1
  }
}

file.name <- paste0(wd, "/figures/", clade, "-SE=", se, ".pdf")
pdf(file = file.name, width = 12, height = 10)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], 
             plots[[4]], plots[[5]], plots[[6]],
             nrow=2, ncol=3)
dev.off()





