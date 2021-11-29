setwd("~/Desktop/climate_niche_seed_dispersal/seed_dispersal")
# rm(list=ls())
library(ape)
library(ggplot2)

#-----------------------------
#-----------------------------
#-----------------------------
# Old results (hybrid models only)
#results.dir <- paste0(getwd(),"/tables")
#results.files <- list.files(results.dir)[grep("TRUE",list.files(results.dir))]
#labels <- sub("-SE.TRUE.csv","", results.files)
#results <- lapply(paste0(results.dir, "/", results.files), read.csv)
#names(results) <- labels

tree.dir <- paste0(getwd(), "/trees")
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)

trees <- lapply(trees, ladderize)
tree_labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))
names(trees) <- tree_labels

get.node.age <- function (phy) {
  root.node <- length(phy$tip.label)+1
  seq.nodes <- phy$edge
  dists <- phy$edge.length
  res <- numeric(max(phy$edge))
  for (i in seq_len(nrow(seq.nodes))) {
    res[seq.nodes[i, 2]] <- res[seq.nodes[i,1]] + dists[i]
  }
  ages <- abs(round(res,2)-round(max(res),2))
  return(ages)
}

#result_list <- list()
#for(i in 1:length(results)){
#  results[[i]]$ObsSt[which(results[[i]]$ObsSt=="Abiotic")] <- "Dry"
#  results[[i]]$ObsSt[which(results[[i]]$ObsSt=="Biotic")] <- "Fleshy"
#  var1 <- results[[i]]
#  clade_names <- unique(results[[i]]$clade_i)
#  final <- matrix(nrow=0,ncol=7)
#  for(j in 1:length(clade_names)){
#    clade1 <- var1[var1$clade_i %in% clade_names[j],]
#    mean_dry_alpha <- round(mean(clade1[clade1$ObsSt=="Dry","Alpha"]),2)
#    se_dry_alpha <- round(sd(clade1[clade1$ObsSt=="Dry","Alpha"]) / sqrt(length(clade1[clade1$ObsSt#=="Dry","Alpha"])),2)
#    mean_fleshy_alpha <- round(mean(clade1[clade1$ObsSt=="Fleshy","Alpha"]),2)
#    se_fleshy_alpha <- round(sd(clade1[clade1$ObsSt=="Fleshy","Alpha"]) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Alpha"])),2)
    #
    # getting halflife in proportion to tree height
#    tmp_tree <- trees[which(names(trees)==clade_names[j])][[1]]
#    tree_height <- get.node.age(tmp_tree)[Ntip(tmp_tree)+1]
#    halflifes_dry <- log(2)/ clade1[clade1$ObsSt=="Dry","Alpha"]
#    mean_halflife_dry <- round(mean(halflifes_dry / tree_height ), 3)
#    se_halflife_dry <- round(sd(halflifes_dry / tree_height) / sqrt(length(halflifes_dry)), 2)
#    halflifes_fleshy <- log(2)/ clade1[clade1$ObsSt=="Fleshy","Alpha"]
#    mean_halflife_fleshy <- round(mean(halflifes_fleshy / tree_height), 2)
#    se_halflife_fleshy <- round(sd(halflifes_fleshy / tree_height) / sqrt(length(halflifes_fleshy)), 2)
    #
    #mean_halflife_dry <- round(mean(log(2)/ clade1[clade1$ObsSt=="Dry","Alpha"]), 2)
    #se_halflife_dry <- round(sd(log(2)/ clade1[clade1$ObsSt=="Dry","Alpha"])  / sqrt(length(clade1[clade1$ObsSt=="Dry","Alpha"])), 2)
    #mean_halflife_fleshy <- round(mean(log(2)/ clade1[clade1$ObsSt=="Fleshy","Alpha"]), 2)
    #se_halflife_fleshy <- round(sd(log(2)/ clade1[clade1$ObsSt=="Fleshy","Alpha"]) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Alpha"])), 2)
    #
#    mean_dry_sigma <- round(mean(clade1[clade1$ObsSt=="Dry","Sigma"]),2)
#    se_dry_sigma <- round((sd(clade1[clade1$ObsSt=="Dry","Sigma"])  / sqrt(length(clade1[clade1$ObsSt#=="Dry","Sigma"]))),2)
#    mean_fleshy_sigma <- round(mean(clade1[clade1$ObsSt=="Fleshy","Sigma"]),2)
#    se_fleshy_sigma <- round((sd(clade1[clade1$ObsSt=="Fleshy","Sigma"]) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Sigma"]))),2)
    #
#    mean_dry_st_var <- round(mean(clade1[clade1$ObsSt=="Dry","Sigma"] / (2*clade1[clade1$ObsSt=="Dry","Alpha"])),2) 
#    se_dry_st_var <- round(sd(clade1[clade1$ObsSt=="Dry","Sigma"] / (2*clade1[clade1$ObsSt=="Dry","Alpha"])) / sqrt(length(clade1[clade1$ObsSt=="Dry","Sigma"])),2) 
#    mean_fleshy_st_var <- round(mean(clade1[clade1$ObsSt=="Fleshy","Sigma"] / (2*clade1[clade1$ObsSt=="Fleshy","Alpha"])),2) 
#    se_fleshy_st_var <- round(sd(clade1[clade1$ObsSt=="Fleshy","Sigma"] / (2*clade1[clade1$ObsSt=="Fleshy","Alpha"])) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Sigma"])),2)
    #
#    transformed_theta <- exp(clade1[,"Optim"])
#    if(names(results)[i]=="temp") {
#    transformed_theta <- transformed_theta - 273.15
#    }
#    if(names(results)[i]=="arid") {
#    transformed_theta <- transformed_theta * 0.0001
#    }
#    mean_dry_theta_t <- round(mean(transformed_theta[clade1$ObsSt=="Dry"]),2)
#    se_dry_theta_t <- round(sd(transformed_theta[clade1$ObsSt=="Dry"]),2)
#    mean_fleshy_theta_t <- round(mean(transformed_theta[clade1$ObsSt=="Fleshy"]),2)
#    se_fleshy_theta_t <- round(sd(transformed_theta[clade1$ObsSt=="Fleshy"]),2)
    #
#    total_fleshy <- c(clade_names[j],"Fleshy", 
#      #paste0(mean_fleshy_theta," (", se_fleshy_theta, ")"),
#      paste0(mean_fleshy_theta_t," (", se_fleshy_theta_t, ")"),
#      paste0(mean_fleshy_sigma," (", se_fleshy_sigma, ")"),
#      paste0(mean_fleshy_alpha," (", se_fleshy_alpha, ")"),
#      paste0(mean_halflife_fleshy," (", se_halflife_fleshy, ")"),
#      paste0(mean_fleshy_st_var," (", se_fleshy_st_var, ")"))
#      
#    total_dry <- c(clade_names[j],"Dry", 
#      #paste0(mean_dry_theta," (", se_dry_theta, ")"),
#      paste0(mean_dry_theta_t," (", se_dry_theta_t, ")"),
#      paste0(mean_dry_sigma," (", se_dry_sigma, ")"),
#      paste0(mean_dry_alpha," (", se_dry_alpha, ")"),
#      paste0(mean_halflife_dry," (", se_halflife_dry, ")"),
#      paste0(mean_dry_st_var," (", se_dry_st_var, ")"))
#      
#  final <- rbind(final, rbind(total_fleshy, total_dry))
#  }
#  final <- as.data.frame(final)
#  colnames(final) <- c("clade","fruit_type","theta_t","sigma2","alpha","half-life","stationary_var")
#  result_list[[i]] <- final
#  names(result_list)[i] <- labels[i]
#}

#write.csv(result_list, file=paste0(getwd(), "/tables/results_summary.csv"))


#------------------------------------
#------------------------------------
#------------------------------------
# Table 1
# New results (with CID and CD models)
results.dir <- paste0(getwd(),"/tip_avg_tables")
results.files <- list.files(results.dir)[grep("se-TipAvgTable",list.files(results.dir))]
results <- lapply(paste0(results.dir, "/", results.files), read.csv)
names(results) <- results.files

#results.files <- list.files(results.dir)[grep("HMM",list.files(results.dir))]


# Trait information
trait.dir <- paste0(getwd(),"/trait_data")
trait.files <- list.files(trait.dir)[grep("niche.csv",list.files(trait.dir))]
traits <- lapply(paste0(trait.dir, "/", trait.files), read.csv)
names(traits) <- trait.files

clade_names <- unique(unlist(lapply(strsplit(results.files, "-"), "[[", 1)))
climate_datasets <- unique(unlist(lapply(strsplit(results.files, "-"), "[[", 2)))

# Trees
tree.dir <- paste0(getwd(), "/trees")
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
trees <- lapply(trees, ladderize)
tree_labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))
names(trees) <- tree_labels

# substitute Rosaceae for Rosaceae pruned
# results.dir <- paste0(getwd(),"/tip_avg_tables")
# results.files_pruned <- list.files(results.dir)[grep("se-Pruned",list.files(results.dir))]
# results_pruned <- lapply(paste0(results.dir, "/", results.files_pruned), read.csv)
# names(results_pruned) <- results.files_pruned
# results[grep("Rosaceae", names(results))] <- results_pruned

# Prune all tips with no climate data?
for(i in 1:length(tree_labels)) {
  tmp <- results[grep(tree_labels[i], names(results))]
  species_to_keep <- traits[grep(tree_labels[i], names(traits))][[1]]$species
  for(j in 1:length(tmp)){
    tmp[[j]] <- tmp[[j]][which(tmp[[j]]$X %in% species_to_keep),]  
  }
  results[grep(tree_labels[i], names(results))] <- tmp
}
 
get.node.age <- function (phy) {
  root.node <- length(phy$tip.label)+1
  seq.nodes <- phy$edge
  dists <- phy$edge.length
  res <- numeric(max(phy$edge))
  for (i in seq_len(nrow(seq.nodes))) {
    res[seq.nodes[i, 2]] <- res[seq.nodes[i,1]] + dists[i]
  }
  ages <- abs(round(res,2)-round(max(res),2))
  return(ages)
}

result_list_new <- list()
for(i in 1:length(climate_datasets)){
  tmp_climate_set <- results[grep(climate_datasets[i], names(results))]
  final <- matrix(nrow=0,ncol=7)
  
  for(j in 1:length(clade_names)){
    tmp_group_set <- tmp_climate_set[grep(clade_names[j], names(tmp_climate_set))][[1]]
    tmp_trait_set <- traits[grep(clade_names[j], names(traits))][[1]][,c("species","Fruit_type")]
    tmp_trait_group_set <- merge(tmp_group_set, tmp_trait_set, by.x="X", by.y="species")
    tmp_tree <- trees[which(names(trees)==clade_names[j])][[1]]

    mean_dry_alpha <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]),2)
    se_dry_alpha <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])),2)
    mean_fleshy_alpha <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]),2)
    se_fleshy_alpha <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])),2)
    #
    # getting halflife in proportion to tree height
    tree_height <- get.node.age(tmp_tree)[Ntip(tmp_tree)+1]
    halflifes_dry <- log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]
    mean_halflife_dry <- round(mean(halflifes_dry / tree_height ), 3)
    se_halflife_dry <- round(sd(halflifes_dry / tree_height) / sqrt(length(halflifes_dry)), 2)
    halflifes_fleshy <- log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]
    mean_halflife_fleshy <- round(mean(halflifes_fleshy / tree_height), 2)
    se_halflife_fleshy <- round(sd(halflifes_fleshy / tree_height) / sqrt(length(halflifes_fleshy)), 2)
    #
    #mean_halflife_dry <- round(mean(log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]), 2)
    #se_halflife_dry <- round(sd(log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])  / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])), 2)
    #mean_halflife_fleshy <- round(mean(log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]), 2)
    #se_halflife_fleshy <- round(sd(log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])), 2)
    #
    mean_dry_sigma <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"]),2)
    se_dry_sigma <- round((sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"])  / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"]))),2)
    mean_fleshy_sigma <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]),2)
    se_fleshy_sigma <- round((sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]))),2)
    #
    mean_dry_st_var <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])),2) 
    se_dry_st_var <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"])),2) 
    mean_fleshy_st_var <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])),2) 
    se_fleshy_st_var <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"])),2)
    #
    transformed_theta <- exp(tmp_trait_group_set[,"theta"])
    
    if(climate_datasets[i]=="temp.se") {
      transformed_theta <- transformed_theta - 273.15
    }
    if(climate_datasets[i]=="arid.se") {
      transformed_theta <- transformed_theta * 0.0001
    }
    
    mean_dry_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
    se_dry_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
    mean_fleshy_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
    se_fleshy_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
    #
    total_fleshy <- c(clade_names[j],"Fleshy", 
                      #paste0(mean_fleshy_theta," (", se_fleshy_theta, ")"),
                      paste0(mean_fleshy_theta_t," (", se_fleshy_theta_t, ")"),
                      paste0(mean_fleshy_sigma," (", se_fleshy_sigma, ")"),
                      paste0(mean_fleshy_alpha," (", se_fleshy_alpha, ")"),
                      paste0(mean_halflife_fleshy," (", se_halflife_fleshy, ")"),
                      paste0(mean_fleshy_st_var," (", se_fleshy_st_var, ")"))
    
    total_dry <- c(clade_names[j],"Dry", 
                   #paste0(mean_dry_theta," (", se_dry_theta, ")"),
                   paste0(mean_dry_theta_t," (", se_dry_theta_t, ")"),
                   paste0(mean_dry_sigma," (", se_dry_sigma, ")"),
                   paste0(mean_dry_alpha," (", se_dry_alpha, ")"),
                   paste0(mean_halflife_dry," (", se_halflife_dry, ")"),
                   paste0(mean_dry_st_var," (", se_dry_st_var, ")"))
    
    final <- rbind(final, rbind(total_fleshy, total_dry))
  }
  final <- as.data.frame(final)
  colnames(final) <- c("clade","fruit_type","theta_t","sigma2","alpha","half-life","stationary_var")
  result_list_new[[i]] <- final
  names(result_list_new)[i] <- climate_datasets[i]
}

#result_list_new <- list()
#for(i in 1:length(climate_datasets)){
#  tmp_climate_set <- results[grep(climate_datasets[i], names(results))]
#  result_list_new2 <- list()
#  for(j in 1:length(clade_names)){
#    tmp_group_set <- tmp_climate_set[grep(clade_names[j], names(tmp_climate_set))][[1]]
#    tmp_trait_set <- traits[grep(clade_names[j], names(traits))][[1]][,c("species","Fruit_type")]
#    tmp_trait_group_set <- merge(tmp_group_set, tmp_trait_set, by.x="X", by.y="species")
#    tmp_tree <- trees[which(names(trees)==clade_names[j])][[1]]
#    
#    dry_st_var <-tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])
#    fleshy_st_var <-tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])
#    
#    transformed_theta <- exp(tmp_trait_group_set[,"theta"])
#    
#    if(climate_datasets[i]=="temp.se") {
#      transformed_theta <- transformed_theta - 273.15
#    }
#    if(climate_datasets[i]=="arid.se") {
#      transformed_theta <- transformed_theta * 0.0001
#    }
#    
#    dry_theta <- transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]
#    fleshy_theta <- transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]
#    
#    table_results <- rbind(cbind(fleshy_st_var, fleshy_theta, "Fleshy"),
#    cbind(dry_st_var, dry_theta, "Dry"))
#    
#    result_list_new2[[j]] <- as.data.frame(table_results)
#    names(result_list_new2)[j] <- clade_names[j]
#  }
#
#  result_list_new[[i]] <- result_list_new2
#  names(result_list_new)[i] <- climate_datasets[i]
#}


#arid <- do.call(rbind, result_list_new$prec.se)
#plot(arid$fleshy_st_var, arid$fleshy_theta)
#length(arid$fleshy_st_var)

#------------------------------------
#------------------------------------
#------------------------------------

# Other summaries - Mean and Optima comparison (Figure 4)
trait.dir <- paste0(getwd(),"/trait_data")
trait_files <- list.files(trait.dir)[grep("niche",list.files(trait.dir))]
results <- lapply(paste0(trait.dir, "/", trait_files), read.csv)
labels <- sub("_niche.csv","", trait_files)
names(results) <- labels

# workaround to keep only the Rosaceae dataset without outliers
results$Rosaceae <- NULL
names(results)[grep("Rosaceae",names(results))] <- "Rosaceae"
labels <- names(results)

#pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
#colors_states <- pal[c(15,5)]

#for(clade_index in labels) {
#  #pdf(paste0(trait.dir, "/actual_mean/", clade_index,"_actual_mean.pdf"))
#  clade1 <- results[[clade_index]]
#  for(var_index in c("mean_aridity","mean_prec","mean_temp")){
#    actual_unit <- clade1[,var_index]
#    if(var_index=="mean_aridity") {
#      actual_unit <- exp(actual_unit) * 0.0001
#    }
#    if(var_index=="mean_prec") {
#      actual_unit <- exp(actual_unit) 
#    }
#    if(var_index=="mean_temp") {
#      actual_unit <- exp(actual_unit) - 273.15
#    }
#    dry <- actual_unit[which(clade1$Fruit_type=="Dry")]
#    fleshy <- actual_unit[which(clade1$Fruit_type=="Fleshy")]
#    summ_boxplot <- boxplot(dry, fleshy, col=colors_states, names=c("Dry","Fleshy"))
#    title(var_index)
#    #sink(paste0(trait.dir, "/actual_mean/", clade_index ,"_", var_index,"_actual_mean_boxplot.txt"))
#    print(summ_boxplot)
#    #sink()
#  }
#  #dev.off()
#}

all_clades <- matrix(nrow=0, ncol=6) 
for(var_index in c("mean_aridity","mean_prec","mean_temp")){
  for(clade_index in labels) {
  clade1 <- results[[clade_index]]
    actual_unit <- clade1[,var_index]
    if(var_index=="mean_aridity") {
      actual_unit <- exp(actual_unit) * 0.0001
    }
    if(var_index=="mean_prec") {
      actual_unit <- exp(actual_unit) 
    }
    if(var_index=="mean_temp") {
      actual_unit <- exp(actual_unit) - 273.15
    }
    dry <- actual_unit[which(clade1$Fruit_type=="Dry")]
    fleshy <- actual_unit[which(clade1$Fruit_type=="Fleshy")]
    dry <- dry[!is.na(dry)]
    fleshy <- fleshy[!is.na(fleshy)]
    all_clades <- rbind(all_clades, 
                  c(clade_index,var_index,"fleshy",length(fleshy), 
                    round(mean(fleshy),2), round(sd(fleshy) / sqrt(length(fleshy)),2)),
    c(clade_index,var_index,"dry",length(dry), 
      round(mean(dry),2), round(sd(dry) / sqrt(length(dry)),2)))
  }
}
colnames(all_clades) <- c("clade","variable","fruit_type","n","mean_trait", "sd")
#write.csv(all_clades, paste0(getwd(), "/actual_mean/all_clades_actual_mean_se.csv"))

# Making similar table for optimas
results.dir <- paste0(getwd(),"/tip_avg_tables")
results.files <- list.files(results.dir)[grep("se-TipAvgTable",list.files(results.dir))]
results <- lapply(paste0(results.dir, "/", results.files), read.csv)
names(results) <- results.files

# Trait information
trait.dir <- paste0(getwd(),"/trait_data")
trait.files <- list.files(trait.dir)[grep("niche.csv",list.files(trait.dir))]
traits <- lapply(paste0(trait.dir, "/", trait.files), read.csv)
names(traits) <- trait.files

clade_names <- unique(unlist(lapply(strsplit(results.files, "-"), "[[", 1)))
climate_datasets <- unique(unlist(lapply(strsplit(results.files, "-"), "[[", 2)))

# Trees
tree.dir <- paste0(getwd(), "/trees")
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
trees <- lapply(trees, ladderize)
tree_labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))
names(trees) <- tree_labels

# Prune all tips with no climate data?
for(i in 1:length(tree_labels)) {
  tmp <- results[grep(tree_labels[i], names(results))]
  species_to_keep <- traits[grep(tree_labels[i], names(traits))][[1]]$species
  for(j in 1:length(tmp)){
    tmp[[j]] <- tmp[[j]][which(tmp[[j]]$X %in% species_to_keep),]  
  }
  results[grep(tree_labels[i], names(results))] <- tmp
}

result_list_new <- list()
for(i in 1:length(climate_datasets)){
  tmp_climate_set <- results[grep(climate_datasets[i], names(results))]
  final <- matrix(nrow=0,ncol=5)
  
  for(j in 1:length(clade_names)){
    tmp_group_set <- tmp_climate_set[grep(clade_names[j], names(tmp_climate_set))][[1]]
    tmp_trait_set <- traits[grep(clade_names[j], names(traits))][[1]][,c("species","Fruit_type")]
    tmp_trait_group_set <- merge(tmp_group_set, tmp_trait_set, by.x="X", by.y="species")
    tmp_tree <- trees[which(names(trees)==clade_names[j])][[1]]

    transformed_theta <- exp(tmp_trait_group_set[,"theta"])
    if(climate_datasets[i]=="temp.se") {
      transformed_theta <- transformed_theta - 273.15
    }
    if(climate_datasets[i]=="arid.se") {
      transformed_theta <- transformed_theta * 0.0001
    }
    mean_dry_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
    se_dry_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
    mean_fleshy_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
    se_fleshy_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
    #
    total_fleshy <- c(climate_datasets[i], clade_names[j],"fleshy", 
                      #paste0(mean_fleshy_theta," (", se_fleshy_theta, ")"),
                      mean_fleshy_theta_t, se_fleshy_theta_t)
    
    total_dry <- c(climate_datasets[i], clade_names[j],"dry", 
                   #paste0(mean_dry_theta," (", se_dry_theta, ")"),
                   mean_dry_theta_t, se_dry_theta_t)
                  
    final <- rbind(final, rbind(total_fleshy, total_dry))
  }
  final <- as.data.frame(final)
  colnames(final) <- c("dataset","clade","fruit_type","theta_t","theta_sd")
  result_list_new[[i]] <- final
  names(result_list_new)[i] <- climate_datasets[i]
}

# Table 1 - results
# write.csv(result_list_new, file=paste0(getwd(), "/tables/results_summary_new.csv"))

# Standartise both tables
full_table1 <- do.call(rbind , result_list_new)
full_table1 <- full_table1[-which(full_table1$dataset=="pet.se"),]

full_table2 <- as.data.frame(all_clades)

full_table <- cbind(full_table1, full_table2)

full_table <- full_table[,c("dataset","clade","fruit_type","theta_t","theta_sd","mean_trait","sd")]

climate_datasets <- unique(full_table$dataset)
comp_table <- list()
for(i in 1:length(climate_datasets)) {
  tmp <- full_table[full_table$dataset==climate_datasets[i],]
  tmp1 <- tmp[,c("dataset","clade","fruit_type","theta_t","theta_sd")]
  tmp2 <-  tmp[,c("dataset","clade","fruit_type","mean_trait","sd")]
  tmp1$class <- "optimum"
  tmp2$class <- "mean"

  colnames(tmp1) <- colnames(tmp2)
  
  tmp1$mean_trait <- as.numeric(tmp1$mean_trait)
  tmp1$sd <- as.numeric(tmp1$sd)
  
  tmp2$mean_trait <- as.numeric(tmp2$mean_trait)
  tmp2$sd <- as.numeric(tmp2$sd)
  
  comp_table[[i]]  <- rbind(tmp1, tmp2)
}


#------------------------
#-------Figure 4---------
## Making optima vs mean tables
## Plot optima vs mean
#comparison.dir <- paste0(getwd(), "/mean_optimum_comparisons")
#comparison.files <- list.files(comparison.dir, ".csv")
#comparison.tables <- lapply(paste0(comparison.dir,"/" ,comparison.files), read.csv)
#names(comparison.tables) <- c("ai","map","mat")

comparison.tables <- comp_table
names(comparison.tables) <- c("ai","map","mat")

# AI
comparison.tables$ai$clade_plus_fruit <- paste0(comparison.tables$ai$clade, " (",comparison.tables$ai$fruit_type,")")


p <- ggplot(comparison.tables$ai, aes(factor(clade_plus_fruit),
                                      x=mean_trait,xmin=mean_trait-2*sd,xmax=mean_trait+2*sd,color=factor(fruit_type), shape = class))+ #color=factor(fruit_type):factor(class)
  geom_pointrange()+
  theme_bw() #+ 
  #facet_wrap(~fruit_type, ncol=1, nrow=2) 

pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
colors_states <- pal[c(23,5)]

#pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
#colors_states <- pal[c(15,5)]
p + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top") + 
  scale_shape_manual(values = c(0, 16))


pdf(paste0(getwd(), "/figures/ai_comp.pdf"),  height=3, width=4)
p + scale_color_manual(values=colors_states) +
  theme(legend.position  = "top")
dev.off()

# MAP
comparison.tables$map$clade_plus_fruit <- paste0(comparison.tables$map$clade, " (",comparison.tables$map$fruit_type,")")

p <- ggplot(comparison.tables$map, aes(factor(clade_plus_fruit),
                                      x=mean_trait,xmin=mean_trait-2*sd,xmax=mean_trait+2*sd,color=factor(class)))+
  geom_pointrange()+
  theme_bw() #+ 
#facet_wrap(~fruit_type, ncol=1, nrow=2) 

pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
colors_states <- pal[c(15,5)]
p + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top")


pdf(paste0(getwd(), "/figures/map_comp.pdf"), height=3, width=4)
p + scale_color_manual(values=colors_states) +
  theme(legend.position  = "top")
dev.off()

# MAT
comparison.tables$mat$clade_plus_fruit <- paste0(comparison.tables$mat$clade, " (",comparison.tables$mat$fruit_type,")")

p <- ggplot(comparison.tables$mat, aes(factor(clade_plus_fruit),
                                      x=mean_trait,xmin=mean_trait-2*sd,xmax=mean_trait+2*sd,color=factor(class)))+
  geom_pointrange()+
  theme_bw() #+ 
#facet_wrap(~fruit_type, ncol=1, nrow=2) 


pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
colors_states <- pal[c(15,5)]
p + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top")


pdf(paste0(getwd(), "/figures/mat_comp.pdf"), height=3, width=4)
p + scale_color_manual(values=colors_states) +
  theme(legend.position  = "top")
dev.off()

#---------------------------------------
#---------------------------------------
#---------------------------------------

## Theta vs. Stationary_var (Figure 5)

result_list <- list()
for(i in 1:length(climate_datasets)){
  tmp_climate_set <- results[grep(climate_datasets[i], names(results))]
  final <- matrix(nrow=0,ncol=4)
  
  for(j in 1:length(clade_names)){
    tmp_group_set <- tmp_climate_set[grep(clade_names[j], names(tmp_climate_set))][[1]]
    tmp_trait_set <- traits[grep(clade_names[j], names(traits))][[1]][,c("species","Fruit_type")]
    tmp_trait_group_set <- merge(tmp_group_set, tmp_trait_set, by.x="X", by.y="species")
    tmp_tree <- trees[which(names(trees)==clade_names[j])][[1]]
    
    mean_dry_alpha <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]),2)
    se_dry_alpha <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])),2)
    mean_fleshy_alpha <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]),2)
    se_fleshy_alpha <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])),2)
    #
    # getting halflife in proportion to tree height
    tree_height <- get.node.age(tmp_tree)[Ntip(tmp_tree)+1]
    halflifes_dry <- log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]
    mean_halflife_dry <- round(mean(halflifes_dry / tree_height ), 3)
    se_halflife_dry <- round(sd(halflifes_dry / tree_height) / sqrt(length(halflifes_dry)), 2)
    halflifes_fleshy <- log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]
    mean_halflife_fleshy <- round(mean(halflifes_fleshy / tree_height), 2)
    se_halflife_fleshy <- round(sd(halflifes_fleshy / tree_height) / sqrt(length(halflifes_fleshy)), 2)
    #
    #mean_halflife_dry <- round(mean(log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]), 2)
    #se_halflife_dry <- round(sd(log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])  / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])), 2)
    #mean_halflife_fleshy <- round(mean(log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]), 2)
    #se_halflife_fleshy <- round(sd(log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])), 2)
    #
    mean_dry_sigma <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"]),2)
    se_dry_sigma <- round((sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"])  / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"]))),2)
    mean_fleshy_sigma <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]),2)
    se_fleshy_sigma <- round((sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]))),2)
    #
    #
    mean_dry_st_var <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])),2) 
    se_dry_st_var <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"])),2) 
    mean_fleshy_st_var <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])),2) 
    se_fleshy_st_var <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"])),2)
    #
    
    transformed_theta <- exp(tmp_trait_group_set[,"theta"])
    
    if(climate_datasets[i]=="temp.se") {
      transformed_theta <- transformed_theta - 273.15
    }
    if(climate_datasets[i]=="arid.se") {
      transformed_theta <- transformed_theta * 0.0001
    }
    mean_dry_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
    se_dry_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
    mean_fleshy_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
    se_fleshy_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
    #
    total_fleshy <- c(clade_names[j],"Fleshy",mean_fleshy_theta_t, mean_fleshy_st_var)
    total_dry <- c(clade_names[j],"Dry", mean_dry_theta_t, mean_dry_st_var)
    
    final <- rbind(final, rbind(total_fleshy, total_dry))
  }
  final <- as.data.frame(final)
  colnames(final) <- c("clade","fruit_type","theta_t","stationary_var")
  result_list[[i]] <- final
  names(result_list)[i] <- climate_datasets[i]
}

#plot(result_list$arid.se$theta_t ~ result_list$arid.se$stationary_var)
#model<- lm(as.numeric(result_list$arid.se$theta_t) ~ as.numeric(result_list$arid.se$stationary_var))
#summary(model)
#abline(model, col="blue")


#write.csv(result_list, file=paste0(wd, "/tables/results_summary.csv"))

result_list$arid.se$theta_t <- as.numeric(result_list$arid.se$theta_t)
result_list$arid.se$stationary_var <- as.numeric(result_list$arid.se$stationary_var)

result_list$prec.se$theta_t <- as.numeric(result_list$prec.se$theta_t)
result_list$prec.se$stationary_var <- as.numeric(result_list$prec.se$stationary_var)

result_list$temp.se$theta_t <- as.numeric(result_list$temp.se$theta_t)
result_list$temp.se$stationary_var <- as.numeric(result_list$temp.se$stationary_var)


library(ggpmisc)
library(gridExtra)

pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
colors_states <- pal[c(15,5)]

my.formula <- y ~ x
pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
colors_states <- pal[c(23,5)]
reg_ai <- ggplot(data=result_list$arid.se, aes(y=stationary_var, x=theta_t, color=fruit_type)) +
  #geom_smooth(method = "lm", se=F, linetype="dashed" ,formula = my.formula, col=hcl.colors(50, palette = "Inferno", alpha = 1)[10]) +
  geom_smooth(method = stats::lm,  se=F,  fullrange=TRUE, linetype="dashed" ,formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label..)), 
               parse = TRUE) +         
  geom_point(aes(color=fruit_type),size=3) +
  theme_bw() + scale_color_manual(values=colors_states) + ylim(-0.5, 3)
pdf(paste0(getwd(), "/figures/theta_stavar_reg_ai.pdf"), height=2.5, width=4)
reg_ai 
dev.off()

my.formula <- y ~ x
pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
colors_states <- pal[c(23,5)]
reg_mat <- ggplot(data=result_list$temp.se, aes(y=stationary_var , x=theta_t)) +
  geom_smooth(method = stats::lm, se=F, linetype="dashed" ,formula = my.formula, col=hcl.colors(50, palette = "Inferno", alpha = 1)[10]) +
  #stat_poly_eq(formula = my.formula, 
  #             aes(label = paste(..rr.label..)), 
  #             parse = TRUE) +         
  geom_point(aes(color=fruit_type),size=3) +
  theme_bw() + scale_color_manual(values=colors_states) + ylim(-0.5, 3)
pdf(paste0(getwd(), "/figures/theta_stavar_reg_mat.pdf"), height=2.5, width=4)
reg_mat 
dev.off()

my.formula <- y ~ x
pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
colors_states <- pal[c(23,5)]
reg_map <- ggplot(data=result_list$prec.se, aes(y=stationary_var, x= theta_t, color=fruit_type)) +
  #geom_smooth(method = "lm", se=F, linetype="dashed" ,formula = my.formula, col=hcl.colors(50, palette = "Inferno", alpha = 1)[10]) +
  geom_smooth(method = stats::lm,  se=F,  fullrange=TRUE, linetype="dashed" ,formula = my.formula) +
  
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label..)), 
               parse = TRUE) +         
  geom_point(aes(color=fruit_type),size=3) +
  theme_bw() + scale_color_manual(values=colors_states)+ ylim(-0.5, 3) + xlim(0, 2500)
pdf(paste0(getwd(), "/figures/theta_stavar_reg_map.pdf"), height=2.5, width=4.2)
reg_map 
dev.off()

#grid.arrange(reg_ai, reg_mat, reg_map, ncol=1, nrow = 3)



##################################
# 2021-10-29 Revision: making a nicer figure for results
##################################
library(ape)
{
result_list_new <- list()
for(i in 1:length(climate_datasets)){
  tmp_climate_set <- results[grep(climate_datasets[i], names(results))]
  final <- matrix(nrow=0,ncol=8)
  
  for(j in 1:length(clade_names)){
    tmp_group_set <- tmp_climate_set[grep(clade_names[j], names(tmp_climate_set))][[1]]
    tmp_trait_set <- traits[grep(clade_names[j], names(traits))][[1]][,c("species","Fruit_type")]
    tmp_trait_group_set <- merge(tmp_group_set, tmp_trait_set, by.x="X", by.y="species")
    tmp_tree <- trees[which(names(trees)==clade_names[j])][[1]]
    
    mean_dry_alpha <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]),2)
    se_dry_alpha <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])),2)
    mean_fleshy_alpha <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]),2)
    se_fleshy_alpha <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])),2)
    #
    tree_height <- get.node.age(tmp_tree)[Ntip(tmp_tree)+1]
    halflifes_dry <- log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]
    mean_halflife_dry <- round(mean(halflifes_dry / tree_height ), 3)
    se_halflife_dry <- round(sd(halflifes_dry / tree_height) / sqrt(length(halflifes_dry)), 2)
    halflifes_fleshy <- log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]
    mean_halflife_fleshy <- round(mean(halflifes_fleshy / tree_height), 2)
    se_halflife_fleshy <- round(sd(halflifes_fleshy / tree_height) / sqrt(length(halflifes_fleshy)), 2)
    #
    mean_dry_sigma <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"]),2)
    se_dry_sigma <- round((sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"])  / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"]))),2)
    mean_fleshy_sigma <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]),2)
    se_fleshy_sigma <- round((sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]))),2)
    #
    mean_dry_st_var <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])),2) 
    se_dry_st_var <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"])),2) 
    mean_fleshy_st_var <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])),2) 
    se_fleshy_st_var <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"])),2)
    #
    transformed_theta <- exp(tmp_trait_group_set[,"theta"])
    
    if(climate_datasets[i]=="temp.se") {
      transformed_theta <- transformed_theta - 273.15
    }
    if(climate_datasets[i]=="arid.se") {
      transformed_theta <- transformed_theta * 0.0001
    }
    
    mean_dry_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
    se_dry_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
    mean_fleshy_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
    se_fleshy_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
    #
    total_fleshy <- c(clade_names[j],"Fleshy", 
                      #paste0(mean_fleshy_theta," (", se_fleshy_theta, ")"),
                      mean_fleshy_theta_t, se_fleshy_theta_t, 
                      mean_halflife_fleshy, se_halflife_fleshy, 
                      mean_fleshy_st_var, se_fleshy_st_var)
    
    total_dry <- c(clade_names[j],"Dry", 
                   #paste0(mean_dry_theta," (", se_dry_theta, ")"),
                   mean_dry_theta_t, se_dry_theta_t,
                   mean_halflife_dry, se_halflife_dry,
                   mean_dry_st_var, se_dry_st_var)
    
    final <- rbind(final, rbind(total_fleshy, total_dry))
  }
  final <- as.data.frame(final)
  colnames(final) <- c("clade","fruit_type","theta_t","se_theta_t","half-life","se_half-life","stationary_var","se_stationary_var")
  result_list_new[[i]] <- final
  names(result_list_new)[i] <- climate_datasets[i]
}
level_order <- rev(c('Apocynaceae', 'Ericaceae', 'Melastomataceae',"Rosaceae","Solanaceae"))
}
#########################################
# Stationary variance
#########################################

#sds <- sd(as.numeric(result_list_new$arid$stationary_var))
#means <- mean(as.numeric(result_list_new$arid$stationary_var))


result_list_new$arid$stationary_var <- as.numeric(result_list_new$arid$stationary_var)
result_list_new$arid$se_stationary_var <- as.numeric(result_list_new$arid$se_stationary_var)
plot4 <- ggplot(result_list_new$arid, aes(x = factor(clade, level = level_order) , y = stationary_var, group = fruit_type, color = fruit_type)) +
  geom_point(size=2) +
  theme_bw() + 
  geom_errorbar(aes(ymin = stationary_var - se_stationary_var, ymax = stationary_var + se_stationary_var), width = 0.75) +
  #geom_line(linetype = 1, size=0.25)+
  coord_flip(ylim = c(-0.25, 3)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 



#pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
#colors_states <- pal[c(15,5)]
pdf("plot4.pdf",width=3, height=2)
plot4 + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top") + 
  scale_shape_manual(values = c(0, 16))
dev.off()

result_list_new$temp$stationary_var <- as.numeric(result_list_new$temp$stationary_var)
result_list_new$temp$se_stationary_var <- as.numeric(result_list_new$temp$se_stationary_var)
plot5 <- ggplot(result_list_new$temp, aes(x = factor(clade, level = level_order), y = stationary_var, group = fruit_type, color = fruit_type)) +
  geom_point(size=2) +
  theme_bw() +
  geom_errorbar(aes(ymin = stationary_var - se_stationary_var, ymax = stationary_var + se_stationary_var), width =  0.75) +
  #geom_line(linetype = 1, size=0.25)+
  coord_flip(ylim = c(-0.25, 3)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
#pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
#colors_states <- pal[c(15,5)]
pdf("plot5.pdf",width=3, height=2)
plot5 + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top") + 
  scale_shape_manual(values = c(0, 16))
dev.off()

result_list_new$prec$stationary_var <- as.numeric(result_list_new$prec$stationary_var)
result_list_new$prec$se_stationary_var <- as.numeric(result_list_new$prec$se_stationary_var)
plot6 <- ggplot(result_list_new$prec, aes(x = factor(clade, level = level_order), y = stationary_var, group = fruit_type, color = fruit_type)) +
  geom_point(size=2) +
  theme_bw() +
  geom_errorbar(aes(ymin = stationary_var - se_stationary_var, ymax = stationary_var + se_stationary_var), width = 0.75) +
  #geom_line(linetype = 1, size=0.25)+
  coord_flip(ylim = c(-0.25, 3)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
#pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
#colors_states <- pal[c(15,5)]
pdf("plot6.pdf",width=3, height=2)
plot6 + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top") + 
  scale_shape_manual(values = c(0, 16))
dev.off()

#########################################
# Half-life
#########################################
result_list_new$arid$`half-life` <- as.numeric(result_list_new$arid$`half-life`)
result_list_new$arid$`se_half-life` <- as.numeric(result_list_new$arid$`se_half-life`)
plot7 <- ggplot(result_list_new$arid, aes(x = factor(clade, level = level_order), y = `half-life`, group = fruit_type, color = fruit_type)) +
  geom_point(size=2) +
  theme_bw() +
  geom_errorbar(aes(ymin = `half-life` - `se_half-life`, ymax = `half-life` + `se_half-life`), width =  0.75) +
  #geom_line(linetype = 1, size=0.25)+
  coord_flip(ylim = c(0, 0.25)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
#pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
#colors_states <- pal[c(15,5)]
pdf("plot7.pdf", width=3, height=2)
plot7 + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top") + 
  scale_shape_manual(values = c(0, 16))
dev.off()

result_list_new$temp$`half-life` <- as.numeric(result_list_new$temp$`half-life`)
result_list_new$temp$`se_half-life` <- as.numeric(result_list_new$temp$`se_half-life`)
plot8 <- ggplot(result_list_new$temp, aes(x = factor(clade, level = level_order) , y = `half-life`, group = fruit_type, color = fruit_type)) +
  geom_point(size=2) +
  theme_bw() +
  geom_errorbar(aes(ymin = `half-life` - `se_half-life`, ymax = `half-life` + `se_half-life`), width =  0.75) +
  #geom_line(linetype = 1, size=0.25)+
  coord_flip(ylim = c(0, 0.25)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
#pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
#colors_states <- pal[c(15,5)]
pdf("plot8.pdf",width=3, height=2)
plot8 + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top") + 
  scale_shape_manual(values = c(0, 16))
dev.off()

result_list_new$prec$`half-life` <- as.numeric(result_list_new$prec$`half-life`)
result_list_new$prec$`se_half-life` <- as.numeric(result_list_new$prec$`se_half-life`)
plot9 <- ggplot(result_list_new$prec, aes(x = factor(clade, level = level_order), y = `half-life`, group = fruit_type, color = fruit_type)) +
  geom_point(size=2) +
  theme_bw() +
  geom_errorbar(aes(ymin = `half-life` - `se_half-life`, ymax = `half-life` + `se_half-life`), width =  0.75) +
  #geom_line(linetype = 1, size=0.25)+
  coord_flip(ylim = c(0, 0.25)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
#pal <- hcl.colors(30, palette = "Viridis", alpha = 1)
#colors_states <- pal[c(15,5)]
pdf("plot9.pdf",width=3, height=2)
plot9 + scale_color_manual(values=colors_states)+
  theme(legend.position  = "top") + 
  scale_shape_manual(values = c(0, 16))
dev.off()


####################
####################
{
  result_list_new <- list()
  for(i in 1:length(climate_datasets)){
    tmp_climate_set <- results[grep(climate_datasets[i], names(results))]
    final <- matrix(nrow=0,ncol=8)
    
    for(j in 1:length(clade_names)){
      tmp_group_set <- tmp_climate_set[grep(clade_names[j], names(tmp_climate_set))][[1]]
      tmp_trait_set <- traits[grep(clade_names[j], names(traits))][[1]][,c("species","Fruit_type")]
      tmp_trait_group_set <- merge(tmp_group_set, tmp_trait_set, by.x="X", by.y="species")
      tmp_tree <- trees[which(names(trees)==clade_names[j])][[1]]
      
      mean_dry_alpha <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]),2)
      se_dry_alpha <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])),2)
      mean_fleshy_alpha <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]),2)
      se_fleshy_alpha <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])),2)
      #
      tree_height <- get.node.age(tmp_tree)[Ntip(tmp_tree)+1]
      halflifes_dry <- log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"]
      mean_halflife_dry <- round(mean(halflifes_dry / tree_height ), 3)
      se_halflife_dry <- round(sd(halflifes_dry / tree_height) / sqrt(length(halflifes_dry)), 2)
      halflifes_fleshy <- log(2)/ tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"]
      mean_halflife_fleshy <- round(mean(halflifes_fleshy / tree_height), 2)
      se_halflife_fleshy <- round(sd(halflifes_fleshy / tree_height) / sqrt(length(halflifes_fleshy)), 2)
      #
      mean_dry_sigma <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"]),2)
      se_dry_sigma <- round((sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"])  / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"]))),2)
      mean_fleshy_sigma <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]),2)
      se_fleshy_sigma <- round((sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"]))),2)
      #
      mean_dry_st_var <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])),2) 
      se_dry_st_var <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","alpha"])) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Dry","sigma.sq"])),2) 
      mean_fleshy_st_var <- round(mean(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])),2) 
      se_fleshy_st_var <- round(sd(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"] / (2*tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","alpha"])) / sqrt(length(tmp_trait_group_set[tmp_trait_group_set$Fruit_type=="Fleshy","sigma.sq"])),2)
      #
      transformed_theta <- exp(tmp_trait_group_set[,"theta"])
      
      if(climate_datasets[i]=="temp.se") {
        transformed_theta <- transformed_theta - 273.15
      }
      if(climate_datasets[i]=="arid.se") {
        transformed_theta <- transformed_theta * 0.0001
      }
      
      mean_dry_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
      se_dry_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Dry"]),2)
      mean_fleshy_theta_t <- round(mean(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
      se_fleshy_theta_t <- round(sd(transformed_theta[tmp_trait_group_set$Fruit_type=="Fleshy"]),2)
      #
      total_fleshy <- c(clade_names[j],"Fleshy", 
                        #paste0(mean_fleshy_theta," (", se_fleshy_theta, ")"),
                        mean_fleshy_theta_t, se_fleshy_theta_t, 
                        mean_halflife_fleshy, se_halflife_fleshy, 
                        mean_fleshy_st_var, se_fleshy_st_var)
      
      total_dry <- c(clade_names[j],"Dry", 
                     #paste0(mean_dry_theta," (", se_dry_theta, ")"),
                     mean_dry_theta_t, se_dry_theta_t,
                     mean_halflife_dry, se_halflife_dry,
                     mean_dry_st_var, se_dry_st_var)
      
      final <- rbind(final, rbind(total_fleshy, total_dry))
    }
    final <- as.data.frame(final)
    colnames(final) <- c("clade","fruit_type","theta_t","se_theta_t","half-life","se_half-life","stationary_var","se_stationary_var")
    result_list_new[[i]] <- final
    names(result_list_new)[i] <- climate_datasets[i]
  }
  level_order <- rev(c('Apocynaceae', 'Ericaceae', 'Melastomataceae',"Rosaceae","Solanaceae"))
}
full_table_figure <- list()
vars <- c("arid","prec","temp")
for(i in 1:3) {
  one_var <- vars[i]
  one_table <- result_list_new[[grep(one_var, names(result_list_new))]]
  means_one_var <- subset(all_clades, grepl(one_var, all_clades[,2])) 
  means_one_var <- as.data.frame(means_one_var)
  means_one_var$fruit_type <- gsub("fleshy","Fleshy", means_one_var$fruit_type)
  means_one_var$fruit_type <- gsub("dry","Dry", means_one_var$fruit_type)
  means_one_var <- means_one_var[,c("n","mean_trait","sd")]
  full_table_figure[[i]] <- cbind(one_table, means_one_var)
}
names(full_table_figure) <- vars

result_list_new <- full_table_figure
pal <- hcl.colors(30, palette = "Inferno", alpha = 1)
colors_states <- pal[c(23,5)]
pal2 <- hcl.colors(30, palette = "Inferno", alpha = 0.5)
colors_states2 <- pal[c(25,7)]


#########################################
# Figure for theta vs. mean
# Hard-coded (of course)
#########################################
# aridity
#########################################
result_list_new$arid$theta_t <- as.numeric(result_list_new$arid$theta_t)
result_list_new$arid$se_theta_t <- as.numeric(result_list_new$arid$se_theta_t)
result_list_new$arid$mean_trait <- as.numeric(result_list_new$arid$mean_trait)
result_list_new$arid$sd <- as.numeric(result_list_new$arid$sd)

tmp1 <- result_list_new$arid[,c("clade","mean_trait","sd")]
tmp1$class <- paste0(result_list_new$arid$fruit_type, "_mean") 
tmp2 <- result_list_new$arid[,c("clade","theta_t","se_theta_t")]
tmp2$class <- paste0(result_list_new$arid$fruit_type, "_optma")
colnames(tmp2) <- colnames(tmp1) <- c("clade","mean", "sd", "class")
table_for_plot <- rbind(tmp1, tmp2)

plot1 <- ggplot(table_for_plot, aes(x = factor(clade, level = level_order), y = mean, group = class, color = class, linetype = class, shape=class)) +
  geom_point(size=c(2)) +
  theme_bw() +
  #geom_line(size=c(0.25))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width =  0.75, linetype=1) +
  coord_flip() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 

pdf("plot1.pdf",width=3, height=2)
plot1 + scale_color_manual(values=c(colors_states2[1],colors_states[1], colors_states2[2],colors_states[2])) +
  scale_shape_manual(values = c(1, 16, 1, 16)) +
  scale_linetype_manual("", values=c(2,1,2,1)) +
  theme(legend.position  = "top") 

dev.off()

#########################################
# temperature
#########################################
result_list_new$temp$theta_t <- as.numeric(result_list_new$temp$theta_t)
result_list_new$temp$se_theta_t <- as.numeric(result_list_new$temp$se_theta_t)
result_list_new$temp$mean_trait <- as.numeric(result_list_new$temp$mean_trait)
result_list_new$temp$sd <- as.numeric(result_list_new$temp$sd)

tmp1 <- result_list_new$temp[,c("clade","mean_trait","sd")]
tmp1$class <- paste0(result_list_new$temp$fruit_type, "_mean") 
tmp2 <- result_list_new$temp[,c("clade","theta_t","se_theta_t")]
tmp2$class <- paste0(result_list_new$temp$fruit_type, "_optma")
colnames(tmp2) <- colnames(tmp1) <- c("clade","mean", "sd", "class")
table_for_plot <- rbind(tmp1, tmp2)

plot2 <- ggplot(table_for_plot, aes(x = factor(clade, level = level_order) , y = mean, group = class, color = class, linetype = class, shape=class)) +
  geom_point(size=c(2)) +
  theme_bw() +
  #geom_line(size=c(0.25))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.75, linetype=1) +
  coord_flip() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 

pdf("plot2.pdf",width=3, height=2)
plot2 + scale_color_manual(values=c(colors_states2[1],colors_states[1], colors_states2[2],colors_states[2])) +
  scale_shape_manual(values = c(1, 16, 1, 16)) +
  scale_linetype_manual("", values=c(2,1,2,1)) +
  theme(legend.position  = "top") 

dev.off()

#########################################
# precipitation
#########################################
result_list_new$prec$theta_t <- as.numeric(result_list_new$prec$theta_t)
result_list_new$prec$se_theta_t <- as.numeric(result_list_new$prec$se_theta_t)
result_list_new$prec$mean_trait <- as.numeric(result_list_new$prec$mean_trait)
result_list_new$prec$sd <- as.numeric(result_list_new$prec$sd)

tmp1 <- result_list_new$prec[,c("clade","mean_trait","sd")]
tmp1$class <- paste0(result_list_new$prec$fruit_type, "_mean") 
tmp2 <- result_list_new$prec[,c("clade","theta_t","se_theta_t")]
tmp2$class <- paste0(result_list_new$prec$fruit_type, "_optma")
colnames(tmp2) <- colnames(tmp1) <- c("clade","mean", "sd", "class")
table_for_plot <- rbind(tmp1, tmp2)

plot3 <- ggplot(table_for_plot, aes(x = factor(clade, level = level_order) , y = mean, group = class, color = class, linetype = class, shape=class)) +
  geom_point(size=c(2)) +
  theme_bw() +
  #geom_line(size=c(0.25))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.75, linetype=1) +
  coord_flip() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 

pdf("plot3.pdf",width=3, height=2)
plot3 + scale_color_manual(values=c(colors_states2[1],colors_states[1], colors_states2[2],colors_states[2])) +
  scale_shape_manual(values = c(1, 16, 1, 16)) +
  scale_linetype_manual("", values=c(2,1,2,1)) +
  theme(legend.position  = "top") 

dev.off()
