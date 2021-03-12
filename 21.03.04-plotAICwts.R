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
library(gridExtra)

## run
# wd <- "~/Desktop/climate_niche_seed_dispersal/seed_dispersal"
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


#-------------------------
# testing different plot arrangement
#-------------------------

pal <- "BuPu"

pdf(paste0(wd, "/figures/AICcWts/aicw.pdf"), width=12, height=10)

se <- TRUE
# Apocynaceae
clade <- clades[1]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)
  
plot_apo1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="AI") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_apo2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAT") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_apo3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAP") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

# Ericaceae
clade <- clades[2]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_eri1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="AI") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_eri2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAT") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_eri3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAP") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

# Melastomataceae
clade <- clades[3]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_melas1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="AI") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_melas2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAT") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_melas3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAP") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

# Rosaceae
clade <- clades[4]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_rosa1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="AI") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_rosa2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAT") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_rosa3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAP") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

# Solanaceae
clade <- clades[5]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_sol1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="AI") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_sol2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAT") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)

plot_sol3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=1.01, size=5, hjust=0, label="MAP") + 
  xlab("") +
  scale_fill_brewer(palette=pal)

grid.arrange(plot_apo1, plot_apo2, plot_apo3, 
             plot_eri1, plot_eri2, plot_eri3,
             plot_melas1, plot_melas2, plot_melas3,
             plot_rosa1, plot_rosa2, plot_rosa3,
             plot_sol1, plot_sol2, plot_sol3, ncol=3, nrow = 5)

dev.off()





### testing model number 
wd <- "~/2021_SeedDispersal"
setwd(wd)
files.tables <- paste0(wd, "/res_tables/", dir("res_tables/"))
dat.types <- c("temp", "prec", "pet", "arid")
clades <- unique(unlist(lapply(strsplit(dir("res_tables/"), "-"), function(x)x[1])))


se <- TRUE
no.of.complete <- matrix(NA, length(clades), length(dat.types), dimnames = list(clades, dat.types))
for(clade in clades){
  for(dat.type in dat.types){
    ToLoad <- files.tables[grep(dat.type, files.tables)]
    ToLoad <- ToLoad[grep(se, ToLoad)]
    ToLoad <- ToLoad[grep(clade, ToLoad)]
    load(ToLoad)
    no.of.complete[rownames(no.of.complete) %in% clade, colnames(no.of.complete) %in% dat.type] <- length(ResTables)
  }
}

