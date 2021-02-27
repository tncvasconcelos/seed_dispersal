# this script will run corHMM analyses on the 4 + 1 datasets
# set wd
wd <- "~/2021_SeedDispersal/"
setwd(wd)
wd <- paste0(getwd(), "/")
setwd(wd)

## imports
require(corHMM)

## functions

clean_dat <- function(csv, tre, remove=TRUE){
  # read files
  dat <- read.csv(csv)[,c(1,2)]
  phy <- read.tree(tre)
  missing_data <- c("no_info","doubtful")
  to_remove <- "remove" # these ones are to be removed anyway because they are outgroups
  # if unknown data is to be removed, remove it 
  # else, replace with unknown (?)
  if(remove == TRUE){
    phy <- drop.tip(phy, dat[,1][which(dat[,2] %in% c(missing_data, to_remove))]) # removes both missing data and outgroups
    dat <- dat[dat[,1] %in% phy$tip.label,]
    dat <- dat[order(match(dat[,1],phy$tip.label)),]
  }else{
    phy <- drop.tip(phy, dat[,1][which(dat[,2] %in% to_remove)]) # removes only outgroups
    dat[dat[,2] %in% missing_data,2] <- "?"
    dat <- dat[order(match(dat[,1],phy$tip.label)),]
  }
  if(length(phy$tip.label) != length(dat[,1])){
    return(cat("James messed up, the data and phylogeny don't match"))
  }
  return(list(dat = dat, phy = phy))
}

runCorHMM <- function(data_list, name, nStarts=0, nCores=1){
  dat <- data_list$dat
  phy <- data_list$phy
  ER <- getStateMat4Dat(dat, model = "ER")$rate.mat
  ARD <- getStateMat4Dat(dat, model = "ARD")$rate.mat
  # ER 
  res_ER <- corHMM(phy = phy, data = dat, rate.cat = 1, model = "ER", node.states = "marginal", nstarts = nStarts, n.cores = nCores)
  
  # ARD
  res_ARD <- corHMM(phy = phy, data = dat, rate.cat = 1, model = "ARD", node.states = "marginal", nstarts = nStarts, n.cores = nCores)
  
  # ER/ER - the final number in the matrix name (ER.ER ** 1 **) refers to the number of independent parameters being estimated for the hidden part of the model (the parameter process)
  mat_ER.ER.1 <- getFullMat(list(ER, ER), ER)
  res_ER.ER.1 <- corHMM(phy = phy, data = dat, rate.cat = 2, rate.mat = mat_ER.ER.1, node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  
  mat_ER.ER.2 <- getFullMat(list(ER, ER), ARD)
  res_ER.ER.2 <- corHMM(phy = phy, data = dat, rate.cat = 2, rate.mat = mat_ER.ER.2, node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  
  # ARD/ARD
  mat_ARD.ARD.1 <- getFullMat(list(ARD, ARD), ER)
  res_ARD.ARD.1 <- corHMM(phy = phy, data = dat, rate.cat = 2, rate.mat = mat_ARD.ARD.1, node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  
  mat_ARD.ARD.2 <- getFullMat(list(ARD, ARD), ARD)
  res_ARD.ARD.2 <- corHMM(phy = phy, data = dat, rate.cat = 2, rate.mat = mat_ARD.ARD.2, node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  
  # ER/ARD
  mat_ER.ARD.1 <- getFullMat(list(ER, ARD), ER)
  res_ER.ARD.1 <- corHMM(phy = phy, data = dat, rate.cat = 2, rate.mat = mat_ER.ARD.1, node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  
  mat_ER.ARD.2 <- getFullMat(list(ER, ARD), ARD)
  res_ER.ARD.2 <- corHMM(phy = phy, data = dat, rate.cat = 2, rate.mat = mat_ER.ARD.2, node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  obj <- list(res_ER = res_ER,
              res_ARD = res_ARD,
              res_ER.ER.1 = res_ER.ER.1,
              res_ER.ER.2 = res_ER.ER.2,
              res_ARD.ARD.1 = res_ARD.ARD.1,
              res_ARD.ARD.2 = res_ARD.ARD.2,
              res_ER.ARD.1 = res_ER.ARD.1,
              res_ER.ARD.2 = res_ER.ARD.2)
  save(obj, file = paste0(name, "-", format(Sys.time(), "%y_%m_%d"), "-corRes.Rsave"))
  return(obj)
}

## import and organize the data
labels <- unique(unlist(lapply(strsplit(sort(dir("trait_data/")), "_"), function(x) x[1])))
csv <- paste0(wd, "trait_data/", sort(dir("trait_data/")))
csv <- csv[grep("_trait_data", csv)]
tre <- paste0(wd, "trees/", sort(dir("trees/")))
data <- mapply(clean_dat, csv, tre, TRUE)
colnames(data) <- labels
getStateMat4Dat(data[,2]$dat)

# collect and save the modeling results
Results <- vector("list", length(labels))
names(Results) <- labels
for(i in 1:length(labels)){
  Results[[i]] <- runCorHMM(data[,i], name = labels[i], nStarts = 49, nCores = 50)
}

Apoc <- runCorHMM(data[,1], name = labels[1], nStarts = 49, nCores = 10)

