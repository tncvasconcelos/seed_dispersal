# Getting and cleaning distribution points
# setwd("~/Desktop/climate_niche_seed_dispersal/seed_dispersal/")
# rm(list=ls())

library(ape)
library(taxize)
library(rgbif)
tree.dir <- "~/Desktop/climate_niche_seed_dispersal/seed_dispersal/trees"

# Trees
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
trees <- lapply(trees, ladderize)
labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))
names(trees) <- labels

#length(unlist(lapply(trees, "[[","tip.label")))


#-------------------------------
# Download GBIF data
#-------------------------------

# fill in gbif.org credentials
user <- "" # your gbif.org username
pwd <- "" # your gbif.org password
email <- "" # your email

for(i in 1:length(trees)){
  taxized_names <- resolveGBIF(trees[[i]]$tip.label) # it should take about 2-3min
  rgbif::occ_download(rgbif::pred_in("scientificName", taxized_names),
                      pred_in("basisOfRecord", 'PRESERVED_SPECIMEN'),
                      pred("hasCoordinate", TRUE),
                      format = "SIMPLE_CSV", user=user,pwd=pwd,email=email)
}


#-------------------------------
# Cleaning GBIF data
#-------------------------------
library(maptools)
library(sp)
library(rgeos)
library(rworldmap)
library(data.table)

# ...download GBIF tables from your download page on GBIF

points.dir <- "~/Desktop/climate_niche_seed_dispersal/GBIF_points"
point_files <- list.files(points.dir, ".csv")
point_files <- point_files[!grepl("cleaned", point_files) & !grepl("summ_stats", point_files) ]
points <- lapply(paste0(points.dir, "/", point_files), fread,quote="")
names(points) <- unlist(lapply(strsplit(point_files, "_"), "[[", 1))
#setwd(points.dir)

#setwd("~/Desktop/rangers_all_files/rangers")
#source("/Users/thaisvasconcelos/Desktop/rangers_all_files/rangers/R/GetTraitDistribution.R")
#source("/Users/thaisvasconcelos/Desktop/rangers_all_files/rangers/R/GetRanges.R")
#template.map <- readRDS("R/template.map.Rdata")

all_cleaned_points <- list()
for(i in 1:length(points)) {
  cleaned_points <- points[[i]][-which(points[[i]]$species==""),]
  cleaned_points <- RemoveNoDecimal(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveCentroids(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveDuplicates(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveOutliers(cleaned_points, species="species", lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveSeaPoints(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveWrongCountries(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveHerbariaLocalities(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveCropsAndInvasive(cleaned_points, species="species")
  cleaned_points <- RemoveZeros(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  species_to_remove <- UnusualDistributions(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  if(length(species_to_remove)>0) {
    cleaned_points <- cleaned_points[-which(cleaned_points$species%in%species_to_remove),]
  }
  all_cleaned_points[[i]] <- cleaned_points
  write.csv(cleaned_points, file=paste0(names(points)[i], "_cleaned.csv"))
}
names(all_cleaned_points) <- names(points)


# Ploting
for(family_index in 1:length(all_cleaned_points)) {
  points_cleaned <- all_cleaned_points[[family_index]]
  spp <- unique(points_cleaned$species)
  pdf(paste0(points.dir, "/", names(all_cleaned_points)[family_index], "_points.pdf"))
  for(species_index in 1:length(spp)){
    tmp_subset <- as.data.frame(points_cleaned[points_cleaned$species==spp[species_index],])
    coord <- tmp_subset[,c("decimalLongitude","decimalLatitude")]
    coordinates(coord) <- ~ decimalLongitude + decimalLatitude
    plot(wrld_simpl)
    plot(coord, col="red", add=T)
    title(spp[species_index])
    print(species_index)
  }
  dev.off()
}

###
# to recover cleaned points back:
points.dir <- "~/Desktop/climate_niche_seed_dispersal/GBIF_points"
cleaned_point_files <- list.files(points.dir, ".csv")
cleaned_point_files <- cleaned_point_files[grepl("cleaned", cleaned_point_files)]
all_cleaned_points <- lapply(paste0(points.dir, "/", cleaned_point_files), fread)
names(all_cleaned_points) <- unlist(lapply(strsplit(cleaned_point_files, "_"), "[[", 1))

########
# Getting climate data
########
climate_data.dir <- "~/Desktop/climate_niche_seed_dispersal/seed_dispersal/climate_data"
for(family_index in 1:length(all_cleaned_points)) {
  # 1. Getting value for every point
    #allpoints <- ClimateFromPoints(all_cleaned_points[[family_index]], lon="decimalLongitude", lat="decimalLatitude", res=2.5)
    #write.csv(allpoints, file=paste0(names(all_cleaned_points)[family_index], "_allpoints.csv"))
  # 2. Getting summary statistics of climatic variables for each species
  # Thinning occurence data first
  thinned_points <- Thinning(all_cleaned_points[[family_index]], species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 1)
  summstats <- GetClimateSummStats_seed_dispersal(thinned_points, lat = "lat", lon="lon", res=2.5)
  write.csv(summstats, file=paste0(climate_data.dir, "/", names(all_cleaned_points)[family_index], "_summstats.csv"))
}


#-------------------------------
#-------------------------------
# Getting organized table for OUwie
climate_data.dir <- "~/Desktop/climate_niche_seed_dispersal/seed_dispersal/climate_data"
summstats_files <- list.files(climate_data.dir, "summstats.csv")
summstats <- lapply(paste0(climate_data.dir, "/", summstats_files), read.csv)
names(summstats) <- unlist(lapply(strsplit(summstats_files, "_"), "[[", 1))

# Trait data
trait.dir <- "~/Desktop/climate_niche_seed_dispersal/seed_dispersal/trait_data"
trait_files <- list.files(trait.dir, ".csv")
traits <- lapply(paste0(trait.dir, "/", trait_files[grep("trait_data",trait_files)]), read.csv)
labels <- unique(unlist(lapply(strsplit(trait_files, "_"), "[[", 1)))
names(traits) <- labels

for(group_index in 1:length(labels)) {
  group <- names(traits)[group_index]
  group_traits <- traits[[group]][,1:2]
  group_summstats <- summstats[[grep(group, names(summstats))]]

  # Matching datasets
  group_summstats$species <- sub(" ","_", group_summstats$species)
  merged_table <- merge(group_summstats, group_traits, by.x="species", by.y="Species")

  #cleaned_table <- merged_table[,c("species","mean_bio1","se_bio1","mean_bio10","mean_bio11","mean_bio12","se_bio12",
  #                                 "mean_bio13","mean_bio14","Dispersal_mode")]
  cleaned_table <- merged_table[,c("species","mean_bio1","se_bio1","mean_bio12","se_bio12","Dispersal_mode")]
  
  # Creating vectors for temperature and preciptation mean and breadth
  #cleaned_table$temp_breadth <- as.numeric(cleaned_table$mean_bio10) - as.numeric(cleaned_table$mean_bio11)
  #cleaned_table$prec_breadth <- as.numeric(cleaned_table$mean_bio13) - as.numeric(cleaned_table$mean_bio14)
  
  # reorganizing
  #cleaned_table <- cleaned_table[,c("species","mean_bio1","se_bio1","mean_bio12","se_bio12","temp_breadth","prec_breadth","Dispersal_mode")]
  #colnames(cleaned_table) <- c("species","temp","se_temp","prec","se_prec","temp_breadth","prec_breadth","Dispersal_mode")
  cleaned_table <- cleaned_table[,c("species","mean_bio1","se_bio1","mean_bio12","se_bio12","Dispersal_mode")]
  colnames(cleaned_table) <- c("species","temp","se_temp","prec","se_prec","Dispersal_mode")
  
  for(i in sequence(ncol(cleaned_table))) {
    cleaned_table[,i][which(cleaned_table[,i] == Inf)] <- NA
  }
  write.csv(cleaned_table, file=paste0(trait.dir,"/",group, "_niche.csv"))
}


## Maping dispersal modes

#setwd("~/Desktop/rangers_all_files/rangers")
#source("/Users/thaisvasconcelos/Desktop/rangers_all_files/rangers/R/GetTraitDistribution.R")
#source("/Users/thaisvasconcelos/Desktop/rangers_all_files/rangers/R/GetRanges.R")
#template.map <- readRDS("R/template.map.Rdata")
pal <- hcl.colors(30, palette = "SunsetDark", alpha = 1)

for(i in 1:length(cleaned_points)) {
  group <- names(traits)[i]
  traits_group <- traits[[i]][,1:2]
  points <- cleaned_points[[i]]
  ranges <- GetRanges(points)

  traits_group[,1] <- sub("_"," ", traits_group[,1])
  traits_subset <- traits_group[traits_group$Species %in% names(ranges),]
  ranges_subset <- ranges[names(ranges) %in% traits_subset$Species]

  fruit_distribution <- GetTraitDistribution(ranges_subset, traits_subset, type="binary")

  pdf(paste0(getwd(), "/", group, ".pdf"), height=8, width=6)
  par(mfrow=c(2,1))
  plot(template.map, col="gray", legend=F)
  plot(fruit_distribution$Biotic, add=T, col=rev(pal))
  title("Biotic")
  plot(template.map, col="gray", legend=F)
  plot(fruit_distribution$Abiotic, add=T, col=rev(pal))
  title("Abiotic")
  dev.off()
}
getwd()


