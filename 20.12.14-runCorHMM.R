# this script will run corHMM analyses on the 4 datasets
# set wd
wd <- "~/2021_SeedDispersal/"
# wd <- paste0(getwd(), "/")

## imports
require(corHMM)

## functions
clean_dat <- function(csv, tre, remove=TRUE){
  # read files
  dat <- read.csv(csv)[,c(1,2)]
  phy <- read.tree(tre)
  missing_data <- c("remove","no_info","doubtful")
  # if unknown data is to be removed, remove it 
  # else, replace with unknown (?)
  if(remove == TRUE){
    phy <- drop.tip(phy, dat[,1][which(dat[,2] %in% missing_data)])
    dat <- dat[dat[,1] %in% phy$tip.label,]
    dat <- dat[order(match(dat[,1],phy$tip.label)),]
  }else{
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
  # ER 
  res_ER <- corHMM(phy = phy, data = dat, rate.cat = 1, model = "ER", node.states = "marginal", nstarts = nStarts, n.cores = nCores)
  # ARD
  res_ARD <- corHMM(phy = phy, data = dat, rate.cat = 1, model = "ARD", node.states = "marginal", nstarts = nStarts, n.cores = nCores)
  # ER/ER
  res_ER.ER <- corHMM(phy = phy, data = dat, rate.cat = 2, model = "ER", node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  # ARD/ARD
  res_ARD.ARD <- corHMM(phy = phy, data = dat, rate.cat = 2, model = "ARD", node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  # ER/ARD
  R1 <- getStateMat4Dat(dat, model = "ER")$rate.mat
  R2 <- getStateMat4Dat(dat, model = "ARD")$rate.mat
  mat_ER.ARD <- getFullMat(list(R1, R2))
  res_ER.ARD <- corHMM(phy = phy, data = dat, rate.cat = 2, rate.mat = mat_ER.ARD, node.states = "marginal", get.tip.states = TRUE, nstarts = nStarts, n.cores = nCores)
  obj <- list(res_ER = res_ER,
              res_ARD = res_ARD,
              res_ER.ER = res_ER.ER,
              res_ARD.ARD = res_ARD.ARD,
              res_ER.ARD = res_ER.ARD)
  save(obj, file = paste0(name, "-", format(Sys.time(), "%y_%m_%d"), "-corRes.Rsave"))
  return(obj)
}

## import and organize the data
labels <- unlist(lapply(strsplit(sort(dir("trait_data/")), "_"), function(x) x[1]))
csv <- paste0(wd, "trait_data/", sort(dir("trait_data/")))
tre <- paste0(wd, "trees/", sort(dir("trees/")))
data <- mapply(clean_dat, csv, tre, TRUE)
colnames(data) <- labels

# collect and save the modeling results
Results <- vector("list", length(labels))
names(Results) <- labels
for(i in 1:length(labels)){
  Results[[i]] <- runCorHMM(data[,i], name = labels[i], nStarts = 9, nCores = 10)
}



