makeTable <- function(files.tables, dat.type, clade, se){
  ToLoad <- files.tables[grep(dat.type, files.tables)]
  ToLoad <- ToLoad[grep(se, ToLoad)]
  ToLoad <- ToLoad[grep(clade, ToLoad)]
  load(ToLoad)
  n.fit.maps <- length(ResTables)
  models.maps <- gsub("res_", "", names(ResTables))
  string.number.loc <- str_locate_all(models.maps, "[1-9]")
  for(i in 1:length(models.maps)){
    if(length(string.number.loc[[i]]) == 0){
      models.maps[i] <- models.maps[i]
    }else{
      if(length(grep("\\.", models.maps[i])) == 0){
        models.maps[i] <- str_sub(models.maps[i], 1L, string.number.loc[[i]][1,1]-1)
      }else{
        models.maps[i] <- str_sub(models.maps[i], 1L, string.number.loc[[i]][1,1])
      }
    }
  }
  names.maps <- rep(models.maps, each = 7)
  AICwt <- unlist(lapply(ResTables, function(x) x[,2]))
  unlist(lapply(ResTables, function(x) x[,2]))
  names.ou <- gsub("-", "", unlist(lapply(ResTables, function(x) rownames(x))))
  obj <- data.frame(dat.type=dat.type, clade=clade, se=se, 
                    model.cor=names.maps, model.ou=names.ou, AICcwt=AICwt)
  return(obj)
}


require(corHMM)
require(OUwie)
require(ggplot2)
require(stringr)

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
files.tables <- paste0(wd, "/res_tables/", dir("res_tables/"))
dat.types <- c("temp", "prec", "pet", "arid")
clades <- unique(unlist(lapply(strsplit(dir("res_tables/"), "-"), function(x)x[1])))

file.name <- paste0("~/2021_SeedDispersal/figures/AICcWtsByCor/", "AICSubsetPlots.pdf")
pdf(file = file.name, width = 8, height = 8)
for(se in c(TRUE, FALSE)){
  for(clade in clades){
    for(dat.type in dat.types){
      obj <- makeTable(files.tables, dat.type, clade, se)
      # file.name <- paste(c(obj$clade[1], obj$dat.type[1], obj$se[1], "pdf"),collapse = ".")
      # file.name <- paste0("~/2021_SeedDispersal/figures/AICcWtsByCor/", file.name)
      # pdf(file = file.name, width = 12, height = 10)
      print(ggplot(obj, aes(x=model.cor, y=AICcwt, fill=model.ou)) + 
              ggtitle(paste(c(obj$clade[1], obj$dat.type[1], obj$se[1]),collapse = "-")) +
              geom_boxplot())
      # dev.off()
    }
  }
}
dev.off()

file.name <- paste0("~/2021_SeedDispersal/figures/AICcWts/", "AICPlots.pdf")
pdf(file = file.name, width = 8, height = 8)
for(se in c(TRUE, FALSE)){
  for(clade in clades){
    for(dat.type in dat.types){
      obj <- makeTable(files.tables, dat.type, clade, se)
      print(ggplot(obj, aes(x=model.ou, y=AICcwt)) + 
              ggtitle(paste(c(obj$clade[1], obj$dat.type[1], obj$se[1]),collapse = "-")) +
              geom_boxplot())
    }
  }
}
dev.off()



# run some predictive tests about the models
big.obj <- c()
for(se in c(TRUE, FALSE)){
  for(clade in clades){
    for(dat.type in dat.types){
      obj <- makeTable(files.tables, dat.type, clade, se)
      big.obj <- rbind(big.obj, obj)
    }
  }
}

TopModels <- big.obj[as.logical(round(big.obj$AICcwt)),c(4,5)]
TopModels[,1] <- factor(TopModels[,1], levels = c("ER", "ARD", "ER.ER.", "ER.ARD.", "ARD.ARD."))
TopModels[,2] <- factor(TopModels[,2], levels=c("OU1", "OUM", "OUMV", "OUMA", "OUMVA"))

c.table <- xtabs(~model.cor + model.ou, TopModels)
X <- chisq.test(c.table)
round(X$expected)
c.table
round(X$residuals)

