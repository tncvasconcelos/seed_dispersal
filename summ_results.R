setwd("~/Desktop/climate_niche_seed_dispersal/seed_dispersal")

results.dir <- paste0(getwd(),"/tables")
results.files <- list.files(results.dir)[grep("TRUE",list.files(results.dir))]
labels <- sub("-SE.TRUE.csv","", results.files)
results <- lapply(paste0(results.dir, "/", results.files), read.csv)
names(results) <- labels

result_list <- list()
for(i in 1:length(results)){
  results[[i]]$ObsSt[which(results[[i]]$ObsSt=="Abiotic")] <- "Dry"
  results[[i]]$ObsSt[which(results[[i]]$ObsSt=="Biotic")] <- "Fleshy"
  var1 <- results[[i]]
  clade_names <- unique(results[[i]]$Clade)
  final <- matrix(nrow=0,ncol=7)
  for(j in 1:length(clade_names)){
    clade1 <- var1[var1$Clade %in% clade_names[j],]
    mean_dry_alpha <- round(mean(clade1[clade1$ObsSt=="Dry","Alpha"]),2)
    se_dry_alpha <- round(sd(clade1[clade1$ObsSt=="Dry","Alpha"]) / sqrt(length(clade1[clade1$ObsSt=="Dry","Alpha"])),2)
    mean_fleshy_alpha <- round(mean(clade1[clade1$ObsSt=="Fleshy","Alpha"]),2)
    se_fleshy_alpha <- round(sd(clade1[clade1$ObsSt=="Fleshy","Alpha"]) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Alpha"])),2)
    #
    mean_halflife_dry <- round(mean(log(2)/ clade1[clade1$ObsSt=="Dry","Alpha"]), 2)
    se_halflife_dry <- round(sd(log(2)/ clade1[clade1$ObsSt=="Dry","Alpha"])  / sqrt(length(clade1[clade1$ObsSt=="Dry","Alpha"])), 2)
    mean_halflife_fleshy <- round(mean(log(2)/ clade1[clade1$ObsSt=="Fleshy","Alpha"]), 2)
    se_halflife_fleshy <- round(sd(log(2)/ clade1[clade1$ObsSt=="Fleshy","Alpha"]) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Alpha"])), 2)
    #
    mean_dry_sigma <- round(mean(clade1[clade1$ObsSt=="Dry","Sigma"]),2)
    se_dry_sigma <- round((sd(clade1[clade1$ObsSt=="Dry","Sigma"])  / sqrt(length(clade1[clade1$ObsSt=="Dry","Sigma"]))),2)
    mean_fleshy_sigma <- round(mean(clade1[clade1$ObsSt=="Fleshy","Sigma"]),2)
    se_fleshy_sigma <- round((sd(clade1[clade1$ObsSt=="Fleshy","Sigma"]) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Sigma"]))),2)
    #
    mean_dry_theta <- round(mean(clade1[clade1$ObsSt=="Dry","Optim"]),2)
    se_dry_theta <- round(sd(clade1[clade1$ObsSt=="Dry","Optim"]) / sqrt(length(clade1[clade1$ObsSt=="Dry","Optim"])) ,2)
    mean_fleshy_theta <- round(mean(clade1[clade1$ObsSt=="Fleshy","Optim"]),2)
    se_fleshy_theta <- round(sd(clade1[clade1$ObsSt=="Fleshy","Optim"]) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Optim"])) ,2)
    #
    mean_dry_theta <- round(mean(clade1[clade1$ObsSt=="Dry","Optim"]),2)
    se_dry_theta <- round(sd(clade1[clade1$ObsSt=="Dry","Optim"]) / sqrt(length(clade1[clade1$ObsSt=="Dry","Optim"])) ,2)
    mean_fleshy_theta <- round(mean(clade1[clade1$ObsSt=="Fleshy","Optim"]),2)
    se_fleshy_theta <- round(sd(clade1[clade1$ObsSt=="Fleshy","Optim"]) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Optim"])) ,2)
    #
    mean_dry_st_var <- round(mean(clade1[clade1$ObsSt=="Dry","Sigma"] / (2*clade1[clade1$ObsSt=="Dry","Alpha"])),2) 
    se_dry_st_var <- round(sd(clade1[clade1$ObsSt=="Dry","Sigma"] / (2*clade1[clade1$ObsSt=="Dry","Alpha"])) / sqrt(length(clade1[clade1$ObsSt=="Dry","Sigma"])),2) 
    mean_fleshy_st_var <- round(mean(clade1[clade1$ObsSt=="Fleshy","Sigma"] / (2*clade1[clade1$ObsSt=="Fleshy","Alpha"])),2) 
    se_fleshy_st_var <- round(sd(clade1[clade1$ObsSt=="Fleshy","Sigma"] / (2*clade1[clade1$ObsSt=="Fleshy","Alpha"])) / sqrt(length(clade1[clade1$ObsSt=="Fleshy","Sigma"])),2)
    #
    transformed_theta <- exp(clade1[,"Optim"])
    if(names(results)[i]=="temp") {
    transformed_theta <- transformed_theta - 273.15
    }
    if(names(results)[i]=="arid") {
    transformed_theta <- transformed_theta * 0.0001
    }
    mean_dry_theta_t <- round(mean(transformed_theta[clade1$ObsSt=="Dry"]),2)
    se_dry_theta_t <- round(sd(transformed_theta[clade1$ObsSt=="Dry"]),2)
    mean_fleshy_theta_t <- round(mean(transformed_theta[clade1$ObsSt=="Fleshy"]),2)
    se_fleshy_theta_t <- round(sd(transformed_theta[clade1$ObsSt=="Fleshy"]),2)
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
  result_list[[i]] <- final
  names(result_list)[i] <- labels[i]
}

write.csv(result_list, file=paste0(wd, "/tables/results_summary.csv"))




