# this script is a cleanup script which will check if a run completed and if not, re run it 
# load("Ericaceae-dat.arid-21_02_14-OURes-OUMA-3.Rsave")
# load("Ericaceae-dat.arid-21_02_14-simmap-3.Rsave")

## functions

# rerun an entire clade both optimization and try errors
reRunClade <- function(folder, map_files){
  Rsaves <- dir(folder)
  files <- paste0(folder, "/", Rsaves)
  for(i in 1:length(files)){
    save.obj = FALSE
    load(files[i])
    ou.file.name <- strsplit(files[i], "/")[[1]][length(strsplit(files[i], "/")[[1]])]
    cat("\nStarting...", ou.file.name, ":", i, "of", length(files), ".")
    # re run the try errors
    ErrTry <- which(unlist(lapply(obj, function(x) class(x) == "try-error")))
    if(length(ErrTry) > 0){
      GoodObj <- obj[[(1:100)[-ErrTry][1]]] # an object that contains things we need
      map.file <- matchMapFile(ou.file.name, map_files)
      load(map.file) # load the correct simmap file
      map_i <- simmaps
      cat("\n", length(ErrTry), "try-errors found.")
      dat <- cbind(rownames(GoodObj$data), GoodObj$data)
      model <- GoodObj$model
      save.obj = TRUE
      ncores <- ifelse(length(ErrTry) > 10, 10, length(ErrTry)) 
      out.try <- mclapply(ErrTry, function(x) try(reRunObjTry(dat = dat, model = model, index = x, map = map_i), silent = TRUE), mc.cores = ncores)
      for(j in 1:length(ErrTry)){
        obj[[ErrTry[j]]] <- out.try[[j]]
      }
      map_i <- NULL
    }
    
    # re run the likelihood errors
    ErrLik <- which(unlist(lapply(obj, function(x) try(abs(x$loglik), silent = TRUE) > 1e5)))
    ErrTry <- which(unlist(lapply(obj, function(x) class(x) == "try-error")))
    ErrLik <- ErrLik[!ErrLik %in% ErrTry]
    cat("\n", length(ErrLik), "optim-errors found &", length(ErrTry), "try-errors found.\n")
    if(length(ErrLik) > 0){
      vals <- rowMeans(obj[[(1:100)[!1:100 %in% ErrLik][1]]]$solution)
      save.obj = TRUE
      ncores <- ifelse(length(ErrLik) > 10, 10, length(ErrLik))
      out.lik <- mclapply(ErrLik, function(x) try(reRunObjLik(obj, x, vals), silent = TRUE), mc.cores = ncores)
      for(j in 1:length(ErrLik)){
        obj[[ErrLik[j]]] <- out.lik[[j]]
      }
    }
    if(save.obj){
      save(obj, file = files[i])
    }
  }
}

# rerun try error objects
reRunObjTry <- function(dat, model, index, map){
  mapping <- unlist(lapply(map[[index]]$maps, function(x) names(x[length(x)])))
  nTip <- length(map[[index]]$tip.label)
  TipStates <- mapping[match(match(dat[,1], map[[index]]$tip.label), map[[index]]$edge[,2])]
  dat[,2] <- TipStates
  res <- OUwie(phy=map[[index]],
               data = dat, 
               model = model,
               simmap.tree = TRUE, 
               mserr = ifelse(dim(dat)[2] == 4, "known", "none"),
               algorithm = "three.point",
               ub = 20, 
               quiet = TRUE)
  return(res)
}

# rerun optimizaion failures 
reRunObjLik <- function(obj, index, vals){
  if(any(is.na(vals))){
    vals <- vals[2]
  }else{
    vals <- vals[c(1,2)]
  }
  res <- OUwie(phy=obj[[index]]$phy, 
               data = cbind(rownames(obj[[index]]$data), obj[[index]]$data), 
               model = obj[[index]]$model,
               simmap.tree = TRUE, 
               mserr = ifelse(dim(obj[[index]]$data)[2] == 3, "known", "none"),
               algorithm = "three.point",
               ub = 20, 
               starting.vals = vals,
               quiet = TRUE)
  return(res)
}

# match a simmap file to the out file being rerun 
matchMapFile<- function(ou.file.name, map_files){
  clade <- strsplit(ou.file.name, "-")[[1]][1]
  dat.type <- strsplit(ou.file.name, "-")[[1]][2]
  iter <- strsplit(ou.file.name, "-")[[1]][length(strsplit(ou.file.name, "-")[[1]])]
  se <- ifelse(length(grep("se", ou.file.name))==1, TRUE, FALSE)
  if(se){
    map_file <- map_files[grep(".se-", map_files)]
  }else{
    map_file <- map_files[-grep(".se-", map_files)]
  }
  map_file <- map_file[grep(clade, map_file)]
  map_file <- map_file[grep(dat.type, map_file)]
  map_file <- map_file[grep(iter, map_file)]
  return(map_file)
}


## imports
require(corHMM)
require(OUwie)
require(parallel)

## run
# wd <- "~/2021_SeedDispersal"
wd <- getwd()
setwd(wd)
Folders <- paste0(wd, "/res_ouwie/", dir("res_ouwie/"))
cor_folder <- paste0(wd, "/res_corhmm/", dir("res_corhmm/"))
map_files <- paste0(wd, "/simmaps/", dir("simmaps/"))


reRunClade(Folders[4], map_files)


