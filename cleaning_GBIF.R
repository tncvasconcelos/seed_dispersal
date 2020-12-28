# Getting and cleaning distribution points
# setwd("~/Desktop/climate_niche_seed_dispersal/seed_dispersal/")
# rm(list=ls())

library(ape)
library(taxize)
library(rgbif)
tree.dir <- paste0(getwd(), "/trees")
trait.dir <- paste0(getwd(), "/trait_data")

# Trees
tree_files <- list.files(tree.dir, ".tre")
trees <- lapply(paste0(tree.dir, "/", tree_files), read.tree)
trees <- lapply(trees, ladderize)
labels <- unlist(lapply(strsplit(tree_files, "_"), "[[", 1))
names(trees) <- labels

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
  all_cleaned_points[[i]] <- cleaned_points
  write.csv(cleaned_points, file=paste0(names(points)[i], "_cleaned.csv"))
}
names(all_cleaned_points) <- names(points)

### 
# to recover cleaned points back:
points.dir <- "~/Desktop/climate_niche_seed_dispersal/GBIF_points"
cleaned_point_files <- list.files(points.dir, ".csv")
cleaned_point_files <- cleaned_point_files[grepl("cleaned", cleaned_point_files)]
cleaned_points <- lapply(paste0(points.dir, "/", cleaned_point_files), fread)
names(cleaned_points) <- unlist(lapply(strsplit(cleaned_point_files, "_"), "[[", 1))

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

########
# Getting climate data
########
for(family_index in 1:length(all_cleaned_points)) {
  # 1. Getting value for every point 
  allpoints <- DataFromPoints(all_cleaned_points[[family_index]], lon="decimalLongitude", lat="decimalLatitude", res=2.5)
  write.csv(allpoints, file=paste0(names(all_cleaned_points)[family_index], "_allpoints.csv"))
  # 2. Getting summary statistics of climatic variables for each species
  # Thinning occurence data first
  thinned_points <- Thinning(all_cleaned_points[[family_index]], species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 1)
  summstats <- GetClimateSummStats(thinned_points, lat = "lat", lon="lon", res=2.5)
  write.csv(summstats, file=paste0(names(all_cleaned_points)[family_index], "_summstats.csv"))
}


