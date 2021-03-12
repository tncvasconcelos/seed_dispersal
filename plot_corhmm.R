
# wd <- "~/Desktop/climate_niche_seed_dispersal/seed_dispersal"
setwd(wd)

library(ggrepel)
require(ggplot2)
library(gridExtra)

corhmm.dir <- paste0(wd, "/table_corhmm")

apocynaceae <- read.csv(paste0(corhmm.dir,"/Apocynaceae-21_02_25-ResTable.csv"))
ericaceae <- read.csv(paste0(corhmm.dir,"/Ericaceae-21_01_11-ResTable.csv"))
melastomataceae <- read.csv(paste0(corhmm.dir,"/Melastomataceae-21_01_11-ResTable.csv"))
rosaceae <- read.csv(paste0(corhmm.dir,"/Rosaceae-21_01_11-ResTable.csv"))
solanaceae <- read.csv(paste0(corhmm.dir,"/Solanaceae-21_01_11-ResTable.csv"))

col = "#20A387"
 
pdf(paste0(wd, "/figures/corhmm_all.pdf"), width=8, height=5)

apocynaceae$X <- factor(apocynaceae$X, levels = apocynaceae$X)
apo <- ggplot(data=apocynaceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw() + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = apocynaceae,
                  fontface ="plain", color = "black", size = 4
  )

ericaceae$X <- factor(ericaceae$X, levels = ericaceae$X)
eri <- ggplot(data=ericaceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw() + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = ericaceae,
                  fontface ="plain", color = "black", size = 4
  )

melastomataceae$X <- factor(melastomataceae$X, levels = melastomataceae$X)
melas <- ggplot(data=melastomataceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw() + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = melastomataceae,
                  fontface ="plain", color = "black", size = 4
  )

rosaceae$X <- factor(rosaceae$X, levels = rosaceae$X)
rosa <- ggplot(data=rosaceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw() + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = rosaceae,
                  fontface ="plain", color = "black", size = 4
  )

solanaceae$X <- factor(solanaceae$X, levels = solanaceae$X)
sol <- ggplot(data=solanaceae, aes(x=X, y=AICc, group=1)) +
  geom_line(color=col)+
  theme_bw() + 
  xlab("") +
  geom_point() +
  geom_text_repel(aes(label = AICcWt), data = solanaceae,
                  fontface ="plain", color = "black", size = 4
  )

grid.arrange(apo, eri, melas, rosa, sol, ncol=1, nrow = 5)

dev.off()

