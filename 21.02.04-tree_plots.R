# Checking distribution of seed dispersal mode on trees
# rm(list=ls())
library(ape)
setwd("~/Desktop/climate_niche_seed_dispersal/seed_dispersal")
#setwd("~/2021_SeedDispersal/")

tree.dir <- paste0(getwd(), "/trees")
trait.dir <- paste0(getwd(), "/trait_data")
climate_data.dir <-  paste0(getwd(), "/climate_data")

# Trees
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
trees <- lapply(trees, ladderize)
labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))
names(trees) <- labels

# Trait data
trait_files <- list.files(trait.dir, ".csv")
traits <- lapply(paste0(trait.dir, "/", trait_files[grep("trait_data", trait_files)]), read.csv)
niches <- lapply(paste0(trait.dir, "/", trait_files[grep("niche", trait_files)]), read.csv)
labels <- unlist(lapply(strsplit(trait_files, "_"), "[[", 1))
names(traits) <- unique(labels)
names(niches) <- unique(labels)

# Climate data
summstats_files <- list.files(climate_data.dir, "summstats.csv")
summstats <- lapply(paste0(climate_data.dir, "/", summstats_files), read.csv)
names(summstats) <- unlist(lapply(strsplit(summstats_files, "_"), "[[", 1))

# Tree plots
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

  pdf(paste0(getwd(), "/figures/", group, "_dispersal_mode_overview.pdf"), width= 4, height= 12)
  plot(group_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.08)
  par(fg="transparent")
  tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.3,lwd=0.2, frame = "n")
  par(fg="black")
  #tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)

  legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
  title(main=paste0(group))
  axisPhylo()
  dev.off()
}


# Tree plots with climatic data
for(group_index in 1:length(niches)) {
  # Getting individual data for each group
  group <- names(niches)[group_index]
  group_tree <- trees[[grep(group, names(trees))]]
  group_niches <- niches[[grep(group, names(niches))]]

  # Matching datasets
  tree_pruned <- keep.tip(group_tree, group_niches$species)
  # Small corrections so that tips appear in the right order
  is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
  ordered_tips <- tree_pruned$edge[is_tip, 2]
  right_order <- tree_pruned$tip.label[ordered_tips]

  # Organizing so climatic data is in the same order as tips of the tree
  cleaned_table <- group_niches[match(right_order, group_niches$species),]
  # The vector has to be in the reverse order to match the tips when ploting
  se_temp <- rev(as.numeric(cleaned_table$se_temp))
  se_temp[which(is.na(se_temp))] <- 0
  temp <- rev(as.numeric(cleaned_table$temp))
  se_prec <- rev(as.numeric(cleaned_table$se_prec))
  se_prec[which(is.na(se_prec))] <- 0
  prec <- rev(as.numeric(cleaned_table$prec))
  
  # Setting dispersal mode for plot
  mode <- rev(cleaned_table$Dispersal_mode)
  names(mode) <- rev(cleaned_table$species)
  colors_states <- c("midnightblue", "goldenrod")
  mode_cols <- mode
  mode_cols[mode_cols=="Biotic"] <- "goldenrod"
  mode_cols[mode_cols=="Abiotic"] <- "midnightblue"
  mode_cols <- mode_cols[match(names(mode_cols),tree_pruned$tip.label)]

  if(group %in% c("Melastomataceae")) { # still no idea why Melas doesn't plot right at first
    mode_cols <- rev(mode_cols)
  }

  # Plotting
  pdf(paste0(getwd(), "/figures/", group, "_tree_and_niche.pdf"), height=6, width=4)
  layout.matrix <- matrix(c(1,2,3), nrow = 1, ncol = 3)
  layout(mat = layout.matrix,widths = c(2,1,1))
  par(mar=c(5,0.5,3,0.5))

  # Tree
  plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
  legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
  title(main=paste0(group))
  axisPhylo()

  # Temperature
  x <- 1:length(temp)
  plot(temp, x, xlim=range(c(temp-se_temp, temp+se_temp)),
     pch=19, yaxt = "n", xlab="Temperature (C*10)", ylab="", frame.plot=T, cex=0.3, col=mode_cols)
  suppressWarnings(arrows(temp-se_temp, x, temp+se_temp, x, length=0.05, angle=90, code=3, lwd=0.1, col=mode_cols))

  # Precipitation
  x <- 1:length(prec)
  plot(prec, x, xlim=range(c(prec-se_prec, prec+se_prec)),
     pch=19, yaxt = "n", xlab="Precipitation (mm)", ylab="", frame.plot=T, cex=0.3,col=mode_cols)
  suppressWarnings(arrows(prec-se_prec, x, prec+se_prec, x, length=0.05, angle=90, code=3, lwd=0.1, col=mode_cols))
  dev.off()
}



