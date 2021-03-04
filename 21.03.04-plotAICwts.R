makeTable <- function(ResTables, dat.type, clade, se){
  ToLoad <- ResTables[grep(dat.type, ResTables)]
  ToLoad <- ToLoad[grep(se, ToLoad)]
  ToLoad <- ToLoad[grep(clade, ToLoad)]
  load(ToLoad)
  n.fit.maps <- length(ResTables)
  models.maps <- gsub("res_", "", names(ResTables))
  names.maps <- unlist(lapply(strsplit(models.maps, "\\."), function(x) paste(x[1:(length(x)-1)], collapse=".")))
  names.maps <- rep(names.maps, each = 7)
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

## run
wd <- "~/2021_SeedDispersal"
# wd <- getwd()
setwd(wd)
files.tables <- paste0(wd, "/res_tables/", dir("res_tables/"))
dat.types <- c("temp", "prec", "pet", "arid")
clades <- unique(unlist(lapply(strsplit(dir("res_tables/"), "-"), function(x)x[1])))



obj <- dim(obj)[1]
colnames(obj)
boxplot(AICcwt~model.ou:model.cor, data = obj)

file.name
obj <- makeTable(files.tables, dat.types[1], clades[1], FALSE)
ggplot(obj, aes(x=model.cor, y=AICcwt, fill=model.ou)) + 
  ggtitle(paste(c(obj$clade[1], obj$dat.type[1], obj$se[1]),collapse = "-")) +
  geom_boxplot()


ggplot2



