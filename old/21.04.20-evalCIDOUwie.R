# evaluation script for the res_tables for each dataset and set of simmaps. mainly responisble for tip averaging.
# generates tip averages found in tip average tables folder

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




