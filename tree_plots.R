# Checking distribution of seed dispersal mode on trees
rm(list=ls())
library(ape)
setwd("~/Desktop/climate_niche_seed_dispersal/data")
setwd("~/2021_SeedDispersal/")

tree.dir <- paste0(getwd(), "/trees")
trait.dir <- paste0(getwd(), "/trait_data")

# Trees
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
trees <- lapply(trees, ladderize)
labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))
names(trees) <- labels

# Trait data
trait_files <- list.files(trait.dir, ".csv")
traits <- lapply(paste0(trait.dir, "/", trait_files), read.csv)
labels <- unlist(lapply(strsplit(trait_files, "_"), "[[", 1))
names(traits) <- labels

# 
# plots

for(group_index in 1:length(labels)) {
  group <- names(traits)[group_index]
  group_traits <- traits[[group]][,1:2] 
  group_tree <- trees[[grep(group, names(trees))]]

  to_drop <- c("remove","no_info","doubtful")
  if(any(to_drop %in% group_traits[,2])) {
    group_tree <- drop.tip(group_tree, group_traits[,1][which(group_traits[,2] %in% to_drop)])
  }
  group_traits <- group_traits[group_traits[,1] %in% group_tree$tip.label,]
  group_traits <- group_traits[order(match(group_traits[,1],group_tree$tip.label)),]
  
  mode <- group_traits[,2]
  names(mode) <- group_traits[,1] 
  colors_states <- c("midnightblue", "goldenrod")
 
  #tip.cols <- colors_states[as.factor(mode)]
 
  pdf(paste0(group, "_dispersal_mode_overview.pdf"), width= 4, height= 12)
  plot(group_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.08)
  par(fg="transparent")
  tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.3,lwd=0.2, frame = "n")
  par(fg="black")
  #tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)

  legend("bottomleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
  title(main=paste0(group))
  axisPhylo()
  dev.off()
}


