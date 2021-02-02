# this script will load, summarize, and analyze the results of 1000 OUwie runs on our simmaps 

## functions
# a function which will return a list of tables that have the model parameters and AIC wt 
getSummaryTable <- function(file){
  load(file)
  obj <- obj_j
  nMap <- length(obj[[1]])
  ErrMat <- matrix(unlist(lapply(obj, function(x) lapply(x, function(y) class(y) == "try-error"))), nMap, 7)
  ErrIndex <- which(!apply(ErrMat, 1, any))
  nRegime <- length(obj[[7]][[ErrIndex[1]]]$tot.states)
  ResTables <- vector("list", length(ErrIndex))
  count <- 1
  # get the relavent params and stats for each map i and model in the set
  for(i in ErrIndex){
    obj_i <- lapply(obj, function(x) x[[i]])
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
getFigureTable <- function(clade, dat.type, se, param, Rsaves){
  files <- Rsaves[grep(clade, Rsaves)]
  if(se == TRUE){
    files <- files[grep("se-", files)]
  }else{
    files <- files[-grep("se-", files)]
  }
  file <- files[grep(dat.type, files)]
  SummaryTables <- getSummaryTable(file)
  ParamTable <- do.call(rbind, lapply(SummaryTables, summarizeTable))
  SubsetTable <- ParamTable[,grep(param, colnames(ParamTable))]
  n <- dim(SubsetTable)[1]
  FigureTable <- data.frame(
    ObsState = rep(rep(c("Abiotic", "Biotic"), each = n), 2),
    HidState = rep(c("A", "B"), each = n*2),
    Val = as.vector(SubsetTable))
  return(FigureTable)
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
Rsaves <- paste0(wd, "/res_ouwie/", dir("res_ouwie/"))
labels <- unlist(lapply(strsplit(dir("res_ouwie/"), "-"), function(x) paste0(x[1], x[2])))
models <- c("BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
params <- c(paste0("Alpha_", 1:4), paste0("Sigma_", 1:4), paste0("Optim_", 1:4))

# gets params and AIC tables
AICTable <- matrix(0, length(labels), length(models))
ParamTable <- matrix(0, length(labels), length(params))
rownames(AICTable) <- rownames(ParamTable) <- labels
colnames(AICTable) <- models
colnames(ParamTable) <- params
for(i in 1:length(Rsaves)){
  if(i != 13){
    file <- Rsaves[i]
    SummaryTables <- getSummaryTable(file)
    AvgParams <- colMeans(do.call(rbind, lapply(SummaryTables, summarizeTable)))
    AvgAICcWt <- colMeans(do.call(rbind, lapply(SummaryTables, function(x) x[,2])))
    
    ParamTable[i,] <- round(AvgParams, 2)
    AICTable[i,] <- round(AvgAICcWt, 2)
  }
}

# nice figures maybe
# states 1 and 3 are abiotic. states 2 and 4 are biotic
clades <- gsub("dat", "", unique(unlist(lapply(strsplit(labels, "\\."), function(x) x[1]))))
params <- c("Alpha", "Sigma", "Optim")
dat.types <- c("temp", "prec")

# sepcify
se <- FALSE
clade <- clades[1]
param <- params[1]
dat.type <- dat.types[1]

plots <- vector("list", 3)
for(i in 1:3){
  clade <- clades[i]
  tab <- getFigureTable(clade = clade, dat.type = dat.type, se = se, param = param, Rsaves = Rsaves)
  cols <- viridis(2)
  plots[[i]] <- ggplot(tab, aes(x=ObsState, y=Val, fill=HidState)) + 
    labs(x = param, y = dat.type) +
    theme(axis.text = element_text(size = 11), legend.justification=c(0,0), legend.position=c(0,0.85)) +
    scale_fill_manual(values=cols) + 
    theme(text = element_text(size = 20)) + 
    ggtitle(clade) + 
    geom_boxplot()
}

grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow=1, ncol=3)






