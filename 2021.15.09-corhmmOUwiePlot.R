# wd <- "~/Desktop/climate_niche_seed_dispersal/seed_dispersal"
setwd(wd)

library(ggrepel)
require(ggplot2)
library(gridExtra)
library(tidyverse)

#############################################
# Loading corHMM results
#############################################
corhmm.dir <- paste0(wd, "/table_corhmm")

apocynaceae <- read.csv(paste0(corhmm.dir,"/Apocynaceae-21_08_24-ResTable.csv"))
ericaceae <- read.csv(paste0(corhmm.dir,"/Ericaceae-21_01_11-ResTable.csv"))
melastomataceae <- read.csv(paste0(corhmm.dir,"/Melastomataceae-21_01_11-ResTable.csv"))
rosaceae <- read.csv(paste0(corhmm.dir,"/Rosaceae-21_01_11-ResTable.csv"))
solanaceae <- read.csv(paste0(corhmm.dir,"/Solanaceae-21_01_11-ResTable.csv"))

col = "#59121C"

apocynaceae$X <- factor(apocynaceae$X, levels = apocynaceae$X)
apo <- ggplot(data=apocynaceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw(base_size = 8) + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = apocynaceae,
                  fontface ="plain", color = "black", size = 3
  ) +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 7)) 


ericaceae$X <- factor(ericaceae$X, levels = ericaceae$X)
eri <- ggplot(data=ericaceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw(base_size = 8) + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = ericaceae,
                  fontface ="plain", color = "black", size = 3
  ) +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 7)) 

melastomataceae$X <- factor(melastomataceae$X, levels = melastomataceae$X)
melas <- ggplot(data=melastomataceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw(base_size = 8) + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = melastomataceae,
                  fontface ="plain", color = "black", size = 3
  ) +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 7)) 

rosaceae$X <- factor(rosaceae$X, levels = rosaceae$X)
rosa <- ggplot(data=rosaceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw(base_size = 8) + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = rosaceae,
                  fontface ="plain", color = "black", size = 3
  ) +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 7)) 

solanaceae$X <- factor(solanaceae$X, levels = solanaceae$X)
sol <- ggplot(data=solanaceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw(base_size = 8) + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = solanaceae,
                  fontface ="plain", color = "black", size = 3
  ) +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 7)) 

#############################################
# Loading OUwie results
#############################################

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
clades <- subset(clades, clades!="names")


makeTable <- function(files.tables, dat.type, clade, se=NULL) {
  ToLoad <- files.tables[grep(dat.type, files.tables)]
  ToLoad <- ToLoad[grep(clade, ToLoad)]
  load(ToLoad)
  OUwie_models <- unique(gsub(paste(c("HYB.","CD.","CID."), collapse="|"), "", rownames(one_sim_map)))
  all_aicwt <- as.data.frame(matrix(ncol=2, nrow=0))
  colnames(all_aicwt) <- c("model.ou","AICcwt")
  for(model_index in seq_along(OUwie_models)) {
    one_model <- OUwie_models[model_index]
    for(i in 1:length(AICTables)){
      one_sim_map <- AICTables[[i]]
      if(class(one_sim_map) == "data.frame") {
        tmp_aicwts <- one_sim_map[grep(paste0(one_model,"$"), rownames(one_sim_map)),]$AICwt
        tmp_result <- data.frame(model.ou=rep(one_model, length(tmp_aicwts)), AICcwt=tmp_aicwts)
        all_aicwt <- rbind(all_aicwt, tmp_result)
      }
    }
  }
  #all_aicwt <- subset(all_aicwt, !is.na(all_aicwt))
  return(all_aicwt)
}

pal <- "BuPu"

pal_out <- hcl.colors(7, palette = "BuPu", alpha = 0.75)

se <- TRUE
# Apocynaceae
clade <- clades[1]

files.tables <- subset(files.tables, !grepl("names", files.tables))

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$outlier.high] <- "black"
obj$outlier.color[obj$outlier.low] <- "black"

obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_apo1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  #geom_errorbar() +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="humidity") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.05) #+
  #scale_fill_brewer(palette=pal)


#------
dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA

obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_apo2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  #geom_boxplot(lwd=0.3, outlier.shape = NA) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="temperature") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_apo3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="precipitation") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)


# Ericaceae
clade <- clades[2]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_eri1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="humidity") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)


dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_eri2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="temperature") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)


dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_eri3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="precipitation") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)


# Melastomataceae
clade <- clades[3]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_melas1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="humidity") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_melas2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="temperature") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)


dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_melas3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="precipitation") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)


# Rosaceae
clade <- clades[4]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_rosa1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="humidity") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)


dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_rosa2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="temperature") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_rosa3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="precipitation") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)



# Solanaceae
clade <- clades[5]

dat.type <- dat.types[4]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_sol1 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="humidity") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[1]
obj <- makeTable(files.tables, dat.type, clade, se)
obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                                outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_sol2 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="temperature") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) +
  scale_fill_brewer(palette=pal)

dat.type <- dat.types[2]
obj <- makeTable(files.tables, dat.type, clade, se)

obj <- obj %>% group_by(model.ou) %>% 
  mutate(outlier.high = AICcwt > quantile(AICcwt, .75) + 1.50*IQR(AICcwt),
         outlier.low = AICcwt < quantile(AICcwt, .25) - 1.50*IQR(AICcwt))
obj <- obj %>% mutate(outlier.color = case_when(outlier.high ~ "black",
                                            outlier.low ~ "black"))
obj$outlier.color <- NA
obj$outlier.color[obj$model.ou=="OUMVA"] <- pal_out[1]
obj$outlier.color[obj$model.ou=="OUMV"] <- pal_out[2]
obj$outlier.color[obj$model.ou=="OUMA"] <- pal_out[3]
obj$outlier.color[obj$model.ou=="OUM"] <- pal_out[4]
obj$outlier.color[obj$model.ou=="OU1"] <- pal_out[5]
obj$outlier.color[obj$model.ou=="BMS"] <- pal_out[6]
obj$outlier.color[obj$model.ou=="BM1"] <- pal_out[7]

plot_sol3 <- ggplot(obj, aes(x=model.ou, y=AICcwt, fill=model.ou)) + 
  geom_violin(lwd=0, fill="white") + 
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.5, y=0.9, size=3, hjust=0, label="precipitation") + 
  xlab("") +
  theme(axis.text.y = element_text(colour = 'black', size = 6),
        axis.text.x = element_text(colour = 'black', size = 6)) +
  #geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  geom_jitter(color = obj$outlier.color, width = .2, na.rm = T, size=0.05, alpha=0.2) + 
  scale_fill_brewer(palette=pal)


#####################################################
# Making grid plot
#####################################################
pdf(paste0(wd, "/figures/corhmmOuwie2.pdf"), width=7, height=9)

grid.arrange(
  apo, plot_apo1, plot_apo2, plot_apo3,
  eri, plot_eri1, plot_eri2, plot_eri3,
  melas, plot_melas1, plot_melas2, plot_melas3,
  rosa, plot_rosa1, plot_rosa2, plot_rosa3,
  sol, plot_sol1, plot_sol2, plot_sol3,
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1,1,1),
                        c(2,3,4),
                        c(5,5,5),
                        c(6,7,8),
                        c(9,9,9),
                        c(10,11,12),
                        c(13,13,13),
                        c(14,15,16),
                        c(17,17,17),
                        c(18,19,20))
)

dev.off()


pdf(paste0(wd, "/figures/corhmmOuwie_Apocynaceae.pdf"), width=8, height=1.8)
grid.arrange( apo, plot_apo1, plot_apo2, plot_apo3,
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1,1,1),
                        c(2,3,4))
)
dev.off()

pdf(paste0(wd, "/figures/corhmmOuwie_Ericaceae.pdf"), width=8, height=1.8)
grid.arrange( eri, plot_eri1, plot_eri2, plot_eri3,
              widths = c(1, 1, 1),
              layout_matrix = rbind(c(1,1,1),
                                    c(2,3,4))
)
dev.off()

pdf(paste0(wd, "/figures/corhmmOuwie_Melastomataceae.pdf"), width=8, height=1.8)
grid.arrange( melas, plot_melas1, plot_melas2, plot_melas3,
              widths = c(1, 1, 1),
              layout_matrix = rbind(c(1,1,1),
                                    c(2,3,4))
)
dev.off()

pdf(paste0(wd, "/figures/corhmmOuwie_Rosaceae.pdf"), width=8, height=1.8)
grid.arrange( rosa, plot_rosa1, plot_rosa2, plot_rosa3,
              widths = c(1, 1, 1),
              layout_matrix = rbind(c(1,1,1),
                                    c(2,3,4))
)
dev.off()

pdf(paste0(wd, "/figures/corhmmOuwie_Solanaceae.pdf"), width=8, height=1.8)
grid.arrange( sol, plot_sol1, plot_sol2, plot_sol3,
              widths = c(1, 1, 1),
              layout_matrix = rbind(c(1,1,1),
                                    c(2,3,4))
)
dev.off()
