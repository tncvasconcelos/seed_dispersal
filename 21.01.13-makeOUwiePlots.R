# this script takes OUwie output tables and makes figures 

require(viridis)
require(ggplot2)
require(gridExtra)

wd <- "~/2021_SeedDispersal/"
files <- paste0("~/2021_SeedDispersal/tables/", dir("~/2021_SeedDispersal/tables/"))
cols <- viridis(2)

for(i in 1:length(files)){
  file_i <- files[i]
  data.name <- strsplit(file_i, "/")[[1]][length(strsplit(file_i, "/")[[1]])]
  dat.type <- strsplit(data.name, "-")[[1]][1]
  se <- strsplit(strsplit(data.name, "-")[[1]][2], "\\.")[[1]][2]
  res_table <- read.csv(file_i)
  
  if(dat.type == "temp"){
    res_table$Optim <- exp(res_table$Optim) - 273
    ylabC <- "Temp. Optima (\u00B0C)"
    main <- paste0("Data: Temperature & SE: ", se)
  }
  if(dat.type == "prec"){
    res_table$Optim <- exp(res_table$Optim)
    ylabC <- "Precip. Optima (mm)"
    main <- paste0("Data: Precipitation & SE: ", se)
  }
  if(dat.type == "arid"){
    res_table$Optim <- exp(res_table$Optim)
    ylabC <- "Aridity"
    main <- paste0("Data: Global Aridity Index & SE: ", se)
  }
  if(dat.type == "pet"){
    res_table$Optim <- exp(res_table$Optim)
    ylabC <- "Potential Evapo-Transpiration"
    main <- paste0("Data: Potential Evapo-Transpiration & SE: ", se)
  }
  
  file.name <- paste0(wd, "figures/", dat.type, "-SE.", se, ".pdf")
  
  pdf(file = file.name, width = 12, height = 10)
  grid.arrange(
    ggplot(res_table, aes(x=Clade, y=Alpha, fill = ObsSt)) + 
      labs(x = "", y = "Alpha") +
      scale_fill_manual(values=cols) + 
      scale_colour_manual(values=cols) + 
      theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) + 
      ggtitle(main) +
      geom_boxplot(),
    
    ggplot(res_table, aes(x=Clade, y=Sigma, fill = ObsSt)) + 
      labs(x = "", y = "Sigma") +
      scale_fill_manual(values=cols) + 
      theme(text = element_text(size = 20)) + 
      geom_boxplot(),
    
    ggplot(res_table, aes(x=Clade, y=Optim, fill = ObsSt)) + 
      labs(x = "Clades", y = ylabC) +
      scale_fill_manual(values=cols) + 
      theme(text = element_text(size = 20)) + 
      geom_boxplot(),
    nrow=3, ncol=1
  )
  dev.off()
}



grid.arrange(plots[[1]], plots[[2]], plots[[3]], 
             plots[[4]], plots[[5]], plots[[6]],
             nrow=2, ncol=3)

# plot the phylogeny figure
load(cor_file)
Tmax <- max(branching.times(res$phy))
phy <- ladderize(res$phy)
cols <- viridis(2)[ifelse(res$data[,2] == "Abiotic", 1, 2)]
Xadd <- (0.1 * Tmax)
plot(phy, show.tip.label = FALSE, type = "fan")
# tiplabels(pch = 16, col = cols, cex = 0.5, offset = 4)
offset = 5
offset <- offset * tab$Optim/max(tab$Optim)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip <- 1:lastPP$Ntip
XX <- lastPP$xx[tip]
YY <- lastPP$yy[tip]
tmp <- rect2polar(XX, YY)
tmp.init <- polar2rect(tmp$r + 0, tmp$angle)
tmp.final <- polar2rect(tmp$r + offset, tmp$angle)
XX.init <- tmp.init$x
YY.init <- tmp.init$y
XX.final <- tmp.final$x
YY.final <- tmp.final$y
segments(x0 = XX.init, y0 = YY.init, x1 = XX.final, YY.final, col = cols)

# points(x = XX.init, y = YY.init, pch = 16, col = cols)
# points(x = XX.final, y = YY.final, pch = 16, col = cols)

