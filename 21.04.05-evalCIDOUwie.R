require(ggplot2)
require(gridExtra)

load("~/2021_SeedDispersal/res_CID/Rosaceae.CID.Fits.Rsave")
load("~/2021_SeedDispersal/res_CID/Solanaceae.CID.Fits.Rsave")

load("~/2021_SeedDispersal/res_tables/Rosaceae-arid-TRUE.Rsave")
p1 <- ggplot(rosa.CID.fits, aes(x=OUFit, y=dAICc, fill=Model)) + 
      geom_boxplot()

p2 <- ggplot(rosa.CID.fits, aes(x=OUFit, y=AICc, fill=Model)) + 
  geom_boxplot()

grid.arrange(p1, p2, nrow=1)

load("~/2021_SeedDispersal/res_tables/Solanaceae-arid-TRUE.Rsave")
p1 <- ggplot(sola.CID.fits, aes(x=OUFit, y=dAICc, fill=Model)) + 
  geom_boxplot()

p2 <- ggplot(sola.CID.fits, aes(x=OUFit, y=AICc, fill=Model)) + 
  geom_boxplot()

minAICc <- unlist(lapply(ResTables, function(x) min(x[,1])))
OtherFits <- data.frame(OUFit = "???", Model = "O.G.", AICc = minAICc, dAICc = 0)
sola.CID.fits <- rbind(sola.CID.fits, OtherFits)
p3 <- ggplot(sola.CID.fits, aes(x=OUFit, y=AICc, fill=Model)) + 
  geom_boxplot()


grid.arrange(p1, p2, p3, nrow=1)

CD.rosa <- rosa.CID.fits[rosa.CID.fits$Model == "CD", ]
CID.rosa <- rosa.CID.fits[rosa.CID.fits$Model == "CID", ]

AICcTable <- data.frame(CD = CD.rosa$AICc, CID = CID.rosa$AICc)
colMeans(AICcTable)