# Checking distribution of seed dispersal mode on trees
# rm(list=ls())
library(ape)
library(phytools)
setwd("~/Desktop/climate_niche_seed_dispersal/seed_dispersal")
# setwd("~/2021_SeedDispersal/")

tree.dir <- paste0(getwd(), "/trees")
trait.dir <- paste0(getwd(), "/trait_data")
climate_data.dir <-  paste0(getwd(), "/climate_data")
corhmm.dir <- paste0(getwd(), "/table_corhmm")

# Trees
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
trees <- lapply(trees, ladderize)
labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))

# Trait data
trait_files <- list.files(trait.dir, ".csv")
traits <- lapply(paste0(trait.dir, "/", trait_files[grep("trait_data", trait_files)]), read.csv)
niches <- lapply(paste0(trait.dir, "/", trait_files[grep("niche", trait_files)][-5]), read.csv) # [-5] is a workaround to keep only the Rosaceae database with no outliers

# Climate data
summstats_files <- list.files(climate_data.dir, "summstats.csv")
summstats <- lapply(paste0(climate_data.dir, "/", summstats_files), read.csv)

# ASR results
asr_files <- list.files(corhmm.dir)[grep("ASR", list.files(corhmm.dir))]
asr_files <- subset(asr_files, asr_files!= "Apocynaceae-21_02_25-ASR.csv")
asr <- lapply(paste0(corhmm.dir, "/", asr_files), function(x) read.csv(x)[,-1])

# labeling 
names(asr) <- names(summstats) <- names(traits) <- names(trees) <- names(niches) <- labels

# Tree plots
for(group_index in 1:length(trees)) {
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
for(group_index in 1:length(trees)) {
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
  temp <- rev(as.numeric(cleaned_table$mean_temp))
  se_prec <- rev(as.numeric(cleaned_table$se_prec))
  se_prec[which(is.na(se_prec))] <- 0
  prec <- rev(as.numeric(cleaned_table$mean_prec))
  
  # Setting dispersal mode for plot
  mode <- rev(cleaned_table$Fruit_type)
  names(mode) <- rev(cleaned_table$species)
  colors_states <- c("midnightblue", "goldenrod")
  mode_cols <- mode
  mode_cols[mode_cols=="Fleshy"] <- "goldenrod"
  mode_cols[mode_cols=="Dry"] <- "midnightblue"
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



# Tree plots final
for(group_index in 1:length(trees)) {
  group <- names(traits)[group_index]
  group_traits <- traits[[group]][,1:2]
  group_tree <- trees[[grep(group, names(trees))]]
  group_niches <- niches[[grep(group, names(niches))]]

  to_drop <- c("remove","no_info","doubtful") # Species marked as any of these are dropped from the plots.
  if(any(to_drop %in% group_traits[,2])) {
    group_tree <- drop.tip(group_tree, group_traits[,1][which(group_traits[,2] %in% to_drop)])
  }
  group_traits <- group_traits[group_traits[,1] %in% group_tree$tip.label,]
  group_traits <- group_traits[order(match(group_traits[,1],group_tree$tip.label)),]
  full_mode <- group_traits[,2]
  names(full_mode) <- group_traits[,1]

  no_climatic_data <- group_tree$tip.label[which(!group_tree$tip.label %in% group_niches$species)]
  rows_to_add <- (nrow(group_niches) + 1) : (nrow(group_niches) + length(no_climatic_data))
  group_niches[rows_to_add,"species"] <- no_climatic_data
  group_niches[rows_to_add,"Fruit_type"] <- "no_data"

  # Matching datasets
  tree_pruned <- group_tree
  # Small corrections so that tips appear in the right order
  is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
  ordered_tips <- tree_pruned$edge[is_tip, 2]
  right_order <- tree_pruned$tip.label[ordered_tips]

  # Organizing so climatic data is in the same order as tips of the tree
  cleaned_table <- group_niches[match(right_order, group_niches$species),]
  # The vector has to be in the reverse order to match the tips when ploting
  #se_aridity <- rev(as.numeric(cleaned_table$within_sp_var_aridity))
  #se_aridity[which(is.na(se_aridity))] <- 0
  aridity <- as.numeric(cleaned_table$mean_aridity)
  names(aridity) <- cleaned_table$species
  aridity <- (exp(aridity)) * 0.0001
  aridity[which(is.na(aridity))] <- min(aridity[!is.na(aridity)])
  #
  #se_temp <- rev(as.numeric(cleaned_table$within_sp_var_temp))
  #se_temp[which(is.na(se_temp))] <- 0
  temp <- as.numeric(cleaned_table$mean_temp)
  names(temp) <- cleaned_table$species
  temp <- (exp(temp)) - 273.15
  temp[which(is.na(temp))] <- min(temp[!is.na(temp)])
  #
  #mean(exp(cleaned_table$mean_temp[cleaned_table$Fruit_type=="Dry"]) - 273.15)
  #se_prec <- rev(as.numeric(cleaned_table$within_sp_var_prec))
  #se_prec[which(is.na(se_prec))] <- 0
  #se_prec <- exp(se_prec)
  prec <- as.numeric(cleaned_table$mean_prec)
  names(prec) <- cleaned_table$species
  prec <- (exp(prec))
  prec[which(is.na(prec))] <- min(prec[!is.na(prec)])
  
  # Setting dispersal mode for plot
  if(group=="Ericaceae"){
    mode <- rev(cleaned_table$Fruit_type)
    names(mode) <- rev(cleaned_table$species)
  } else {
    mode <- cleaned_table$Fruit_type
    names(mode) <- cleaned_table$species
  }

  #pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
  #colors_states <- pal[c(15,5)]
  pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
  colors_states <- pal[c(23,5)]
  mode_cols <- mode
  mode_cols[mode_cols=="Fleshy"] <- colors_states[2]
  mode_cols[mode_cols=="Dry"] <- colors_states[1]
  mode_cols[mode_cols=="no_data"] <- "#FFFFFFAA"
  mode_cols <- mode_cols[match(names(mode_cols),tree_pruned$tip.label)]
  
  if(group %in% c("Melastomataceae")) { # still no idea why Melas doesn't plot right at first
    mode_cols <- rev(mode_cols)
  }

#------------------
#------------------

  pdf(paste0(getwd(), "/figures/", group, "_dispersal_mode_final.pdf"), width= 3, height= 12)
  layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
  layout(mat = layout.matrix,widths = c(1.5,0.5,0.5,0.5))
  par(mar=c(5,0.5,3,0.5))
  
  plot(group_tree, show.tip.label=F, edge.width=0.05, adj=1, cex=0.08)
  par(fg="transparent")
  tiplabels(pie=to.matrix(full_mode, sort(unique(full_mode))),piecol=colors_states,cex= 0.5,lwd=0.2, frame = "n")
  nodelabels(node = (Ntip(group_tree)+1):(Ntip(group_tree)+(Ntip(group_tree)-1)) , pie=as.matrix(asr[[group_index]]),cex=1.25, piecol=colors_states)
  
  par(fg="black")
  #tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)

  title(main=paste0(group))
  axisPhylo()
  
  # Aridity index
  aridity <- aridity[order(match(names(aridity), names(mode_cols)))]
  x <- 1:length(aridity)
  plot(aridity, x, lwd = 0.2, xlim=range(c(min(aridity), max(aridity))),
       pch=19, yaxt = "n", xlab="AI", ylab="", frame.plot=T, cex=0.5, col=mode_cols, las=2)
  segments(min(aridity), 1:length(aridity), aridity[1:length(aridity)], 1:length(aridity), col= mode_cols,lwd = 0.2)
  
  # Temperature
  temp <- temp[order(match(names(temp), names(mode_cols)))]
  x <- 1:length(temp)
  plot(temp, x, lwd = 0.2, xlim=range(c(min(temp), max(temp))),
       pch=19, yaxt = "n", xlab="MAT", ylab="", frame.plot=T, cex=0.5, col=mode_cols, las=2)
  segments(min(temp), 1:length(temp), temp[1:length(temp)], 1:length(temp), col= mode_cols,lwd = 0.2)
  
  # Precipitation
  prec <- prec[order(match(names(prec), names(mode_cols)))]
  x <- 1:length(prec)
  plot(prec, x, lwd = 0.2, xlim=range(c(min(prec), max(prec))),
       pch=19, yaxt = "n", xlab="MAP", ylab="", frame.plot=T, cex=0.5,col=mode_cols, las=2)
  segments(min(prec), 1:length(prec), prec[1:length(prec)], 1:length(prec), col= mode_cols,lwd = 0.2)
  
  dev.off()
}






