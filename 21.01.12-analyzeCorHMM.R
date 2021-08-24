# analyze the corHMM results

## imports
require(MASS)

## functions
## functions
# saves the results table and returns it
getResultsTable <- function(file){
  load(file)
  # the results table
  ResTable <- matrix(0, length(obj), 5)
  rownames(ResTable) <- c("ER", "ARD", "ER/ER.1", "ER/ER.2", "ARD/ARD.1", "ARD/ARD.2", "ER/ARD.1", "ER/ARD.2")
  colnames(ResTable) <- c("k.rate", "AICc", "AICcWt", "MeanRate", "ASR")
  count <- 1
  for(i in 1:length(obj)){
    CorRes_i <- obj[[i]]
    Model <- CorRes_i$solution
    AICc <- round(CorRes_i$AICc, 2)
    AICwt <- 0
    MeanRate <- round(mean(CorRes_i$solution, na.rm = TRUE),2)
    Model[is.na(Model)] <- 0
    ASR <-  round(CorRes_i$states[1,][which.max(CorRes_i$states[1,])], 2)
    ASR <-  paste(names(ASR), paste(ASR*100, "%", sep = ""))
    k.rate <- max(CorRes_i$index.mat, na.rm = TRUE)
    diag(Model) <- -rowSums(Model)
    Eq <- c(Null(Model))/sum(Null(Model))
    EntUncond <- sum(Eq*-log2(Eq))
    EntCond <- mean(rowSums(CorRes_i$states * -log2(CorRes_i$states)))
    MutInfo <- round(EntUncond - EntCond, 2)
    #PropInfo <- round(MutInfo/EntUncond*100, 2)
    ResTable[count,] <- c(k.rate, AICc, AICwt, MeanRate, ASR)
    count <- count + 1
  }
  AICcs <- as.numeric(ResTable[,2])
  ResTable[,3] <- round(exp(-0.5 * AICcs - min(AICcs))/sum(exp(-0.5 * AICcs - min(AICcs))),2 )
  ResTable[ResTable[,3] < 0.001,3] <- "<0.001"
  ResTable[ResTable[,4] < 0.001,4] <- "<0.001"
  table.name <- gsub(".Rsave", ".csv", file)
  table.name <- gsub("res_corhmm", "table_corhmm", table.name)
  table.name <- gsub("corRes", "ResTable", table.name)
  write.csv(ResTable, file = table.name)
  return(ResTable)
}

getModelRateMats <- function(file){
  load(file)
  rate.mats <- lapply(obj, function(x) x$solution)
  MaxRate <- max(unlist(lapply(rate.mats, function(x) dim(x)[1])))
  res <- rep("~", MaxRate)
  for(i in 1:length(obj)){
    rate.mat_i <- round(rate.mats[[i]], 3)
    rate.mat_i[!is.na(rate.mat_i) & rate.mat_i == 0] <- "<0.001"
    if(dim(rate.mat_i)[1] < MaxRate){
      rate.mat_i <- rbind(rownames(rate.mat_i), rate.mat_i)
      rate.mat_i <- cbind(rate.mat_i, matrix("!", MaxRate - dim(rate.mat_i)[1] + 2, MaxRate - dim(rate.mat_i)[2]))
      rate.mat_i <- rbind(rate.mat_i, " ")
      rate.mat_i <- rbind(rate.mat_i, "~")
      res <- rbind(res, rate.mat_i)
    }else{
      rate.mat_i <- rbind(rownames(rate.mat_i), rate.mat_i)
      rate.mat_i <- rbind(rate.mat_i, " ")
      rate.mat_i <- rbind(rate.mat_i, "~")
      res <- rbind(res, rate.mat_i)
    }
  }
  table.name <- gsub(".Rsave", ".csv", file)
  table.name <- gsub("res_corhmm", "table_corhmm", table.name)
  table.name <- gsub("corRes", "QMat", table.name)
  write.csv(res, file = table.name)
  return(res)
}

getWeightedASR <- function(file){
  load(file)
  # the results table
  AICcs <- unlist(lapply(obj, function(x) x$AICc))
  AICwt <- exp(-0.5 * AICcs - min(AICcs))/sum(exp(-0.5 * AICcs - min(AICcs)))
  res <- matrix(0, dim(obj[[1]]$states)[1], dim(obj[[1]]$states)[2])
  for(i in 1:length(obj)){
    States <- colnames(obj[[i]]$solution)
    if(length(grep("R2", States)) == 0){
      ASR_i <- obj[[i]]$states[,grep("R1", States)]
    }else{
      ASR_i <- obj[[i]]$states[,grep("R1", States)] + obj[[i]]$states[,grep("R2", States)]
    }
    res <- res + (ASR_i * AICwt[i])
  }
  table.name <- gsub(".Rsave", ".csv", file)
  table.name <- gsub("res_corhmm", "table_corhmm", table.name)
  table.name <- gsub("corRes", "ASR", table.name)
  write.csv(res, file = table.name)
  return(res)
}

# set wd
wd <- "~/2021_SeedDispersal/"
setwd(wd)

Rsaves <- paste0(wd, "res_corhmm/", dir("res_corhmm/"))
labels <- unlist(lapply(strsplit(dir("res_corhmm/"), "-"), function(x) x[1]))

for(i in Rsaves[2]){
  getResultsTable(i)
  getModelRateMats(i)
  getWeightedASR(i)
}
