# Getting and cleaning distribution points
# setwd("~/Desktop/climate_niche_seed_dispersal/seed_dispersal/")
# rm(list=ls())

source(list.files(pattern="functionsClimate.R"))
library(ape)
library(taxize)
library(rgbif)
library(maptools)
library(sp)
library(rgeos)
library(rworldmap)
library(data.table)
data("wrld_simpl")

#-------------------------------
# Download GBIF data
#-------------------------------

# Trees
tree.dir <- "./trees"
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
trees <- lapply(trees, ladderize)
labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))
names(trees) <- labels

# Number of tips in full dataset
# length(unlist(lapply(trees, "[[","tip.label")))
# [1] 5564 

# Fill in GBIF credentials for download
user <- "" # username
pwd <- "" # password
email <- "" # email

{; for(family_index in 1:length(trees)) {
  taxized_names <- resolveGBIF(trees[[family_index]]$tip.label) # it should take about 2-3min for each tree
  rgbif::occ_download(rgbif::pred_in("scientificName", taxized_names),
                      pred_in("basisOfRecord", 'PRESERVED_SPECIMEN'),
                      pred("hasCoordinate", TRUE),
                      format = "SIMPLE_CSV", user=user,pwd=pwd,email=email)
}
beepr::beep("fanfare"); } 

#-------------------------------
# Cleaning GBIF data
#-------------------------------
# ...download GBIF tables from your download page on GBIF
points.dir <- "./GBIF_points"
point_files <- list.files(points.dir, ".csv")
point_files <- point_files[!grepl("cleaned", point_files)]
points <- lapply(paste0(points.dir, "/", point_files), fread,quote="")
names(points) <- unlist(lapply(strsplit(point_files, "_"), "[[", 1))

{; all_cleaned_points <- list()
for(family_index in 1:length(points)) {
  cleaned_points <- points[[family_index]][-which(points[[family_index]]$species==""),]
  cleaned_points <- RemoveNoDecimal(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveCentroids(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveDuplicates(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveOutliers(cleaned_points, species="species", lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveSeaPoints(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveWrongCountries(cleaned_points, lon="decimalLongitude", lat="decimalLatitude", wrld_simpl=wrld_simpl)
  #cleaned_points <- RemoveHerbariaLocalities(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveZeros(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  species_to_remove <- UnusualDistributions(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  if(length(species_to_remove) > 0) {
    cleaned_points <- cleaned_points[-which(cleaned_points$species%in%species_to_remove),]
  }
  all_cleaned_points[[family_index]] <- cleaned_points
  write.csv(cleaned_points, file=paste0(points.dir, "/", names(points)[family_index], "_cleaned.csv"))
}
names(all_cleaned_points) <- names(points)
beepr::beep("fanfare"); } 

# Plotting to inspect distributions
{; for(family_index in 1:length(all_cleaned_points)) {
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
beepr::beep("fanfare"); } 

#-------------------------------
# Getting climate data
#-------------------------------
# Load cleaned points back:
points.dir <- "./GBIF_points"
cleaned_point_files <- list.files(points.dir, ".csv")
cleaned_point_files <- cleaned_point_files[grepl("cleaned", cleaned_point_files)]
all_cleaned_points <- lapply(paste0(points.dir, "/", cleaned_point_files), fread)
names(all_cleaned_points) <- unlist(lapply(strsplit(cleaned_point_files, "_"), "[[", 1))

# Directory to save preliminary datasets:
climate_data.dir <- "./climate_data"
# Directory where climate layers are:
climate_layers.dir <- "./climate_layers"

{; for(family_index in 1:length(all_cleaned_points)) {
    # 1. Thinning occurence data first
    thinned_points <- Thinning(all_cleaned_points[[family_index]], species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 1)
    # 2. Getting summary statistics of climatic variables for each species
    allpoints <- ClimateFromPoint_custom(thinned_points, species="species",lon="lon", lat="lat", layerdir = climate_layers.dir)
    write.csv(allpoints, file=paste0(climate_data.dir, "/", names(all_cleaned_points)[family_index], "_allpoints.csv"))
    summstats <- GetClimateSummStats_custom(allpoints, type="raw")
    write.csv(summstats, file=paste0(climate_data.dir, "/", names(all_cleaned_points)[family_index], "_summstats_raw.csv"))
    summstats <- GetClimateSummStats_custom(allpoints, type="transformed")
    write.csv(summstats, file=paste0(climate_data.dir, "/", names(all_cleaned_points)[family_index], "_summstats.csv"))
  }
beepr::beep("fanfare"); } 

#-------------------------------
# Getting organized table for OUwie
#-------------------------------
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

{; for(family_index in 1:length(labels)) {
  group <- names(traits)[family_index]
  group_traits <- traits[[group]][,1:2]
  group_summstats <- summstats[[grep(group, names(summstats))]]

  # Matching datasets
  group_summstats$species <- sub(" ","_", group_summstats$species)
  merged_table <- merge(group_summstats, group_traits, by.x="species", by.y="Species")
  cleaned_table <- merged_table[,c("species","mean_temp","se_temp","within_sp_var_temp","mean_prec","se_prec","within_sp_var_prec",
                                   "mean_pet","se_pet","within_sp_var_pet","mean_aridity","se_aridity", "within_sp_var_aridity","Fruit_type")]

  if(any(is.na(cleaned_table$mean_prec))) { 
    cleaned_table <- cleaned_table[-which(is.na(cleaned_table$mean_prec)),]
  }
  if(any(cleaned_table$mean_prec==0)) {
    cleaned_table <- cleaned_table[-which(cleaned_table$mean_prec == 0),]
  }
  write.csv(cleaned_table, file=paste0(trait.dir,"/",group, "_niche.csv"))
}
beepr::beep("fanfare"); } 
