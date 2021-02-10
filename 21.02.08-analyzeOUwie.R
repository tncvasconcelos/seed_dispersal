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
    ModelList[[count]] <- obj
    count <- count + 1
  }
  nMap <- length(ModelList[[1]])
  corModels <- names(ModelList[[1]])
  ErrMat <- matrix(unlist(lapply(ModelList, function(x) lapply(x, function(y) class(y) == "try-error"))), nMap, 7)
  ErrIndex <- which(!apply(ErrMat, 1, any))
  cat(nMap - length(ErrIndex), "models resulted in try-errors.\n")
  ResTables <- list()
  count <- 1
  # get the relavent params and stats for each map i and model in the set
  for(i in ErrIndex){
    obj_i <- lapply(ModelList, function(x) x[[i]])
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
require(viridis)
require(ggplot2)
require(gridExtra)

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
Folders <- paste0(wd, "/res_ouwie/", dir("res_ouwie/"))
cor_folder <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))
dat.types <- c("temp", "prec")
params <- c("Alpha", "Sigma", "Optim")

# make plots for a given clade and standard error
se <- TRUE
dat.type <- dat.types[2]
file.name <- paste0(wd, "/figures/", dat.type, "-SE.", se, ".pdf")

res <- vector("list", length(Folders))
for(k in 1:length(Folders)){
  ou_folder <- Folders[k]
  Clade <- strsplit(ou_folder, "/")[[1]][length(strsplit(ou_folder, "/")[[1]])]
  cat("\nStarting", Clade, "...\n")
  cor_file <- cor_folder[grep(clade, cor_folder)]
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
res_table$Optim <- exp(res_table$Optim)
if(dat.type == "temp"){
  res_table$Optim <- res_table$Optim - 273
  ylabC <- "Temp. Optima (\u00B0C)"
  main <- paste0("Data: Temperature & SE: ", se)
}else{
  ylabC <- "Precip. Optima (mm)"
  main <- paste0("Data: Precipitation & SE: ", se)
}
cols <- viridis(2)

pdf(file = file.name, width = 12, height = 10)
grid.arrange(
  ggplot(res_table, aes(x=Clade, y=Alpha, fill = ObsSt)) + 
    labs(x = "", y = "Alpha") +
    scale_fill_manual(values=cols) + 
    scale_colour_manual(values=cols) + 
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) + 
    ggtitle(main) +
    geom_boxplot(),
  
  ggplot(res_table, aes(x=Clade, y=Sigma, fill = ObsSt)) + 
    labs(x = "", y = "Sigma") +
    scale_fill_manual(values=cols) + 
    theme(text = element_text(size = 20)) + 
    geom_boxplot(),
  
  ggplot(res_table, aes(x=Clade, y=Optim, fill = ObsSt)) + 
    labs(x = "Clades", y = ylabC) +
    scale_fill_manual(values=cols) + 
    theme(text = element_text(size = 20)) + 
    geom_boxplot(),
  nrow=3, ncol=1
)
dev.off()






grid.arrange(plots[[1]], plots[[2]], plots[[3]], 
             plots[[4]], plots[[5]], plots[[6]],
             nrow=2, ncol=3)

# plot the phylogeny figure
load(cor_file)
Tmax <- max(branching.times(res$phy))
phy <- ladderize(res$phy)
cols <- viridis(2)[ifelse(res$data[,2] == "Abiotic", 1, 2)]
Xadd <- (0.1 * Tmax)
plot(phy, show.tip.label = FALSE, type = "fan")
# tiplabels(pch = 16, col = cols, cex = 0.5, offset = 4)
offset = 5
offset <- offset * tab$Optim/max(tab$Optim)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip <- 1:lastPP$Ntip
XX <- lastPP$xx[tip]
YY <- lastPP$yy[tip]
tmp <- rect2polar(XX, YY)
tmp.init <- polar2rect(tmp$r + 0, tmp$angle)
tmp.final <- polar2rect(tmp$r + offset, tmp$angle)
XX.init <- tmp.init$x
YY.init <- tmp.init$y
XX.final <- tmp.final$x
YY.final <- tmp.final$y
segments(x0 = XX.init, y0 = YY.init, x1 = XX.final, YY.final, col = cols)

# points(x = XX.init, y = YY.init, pch = 16, col = cols)
# points(x = XX.final, y = YY.final, pch = 16, col = cols)

