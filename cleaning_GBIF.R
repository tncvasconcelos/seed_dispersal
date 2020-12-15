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

# Trait data
trait_files <- list.files(trait.dir, ".csv")
traits <- lapply(paste0(trait.dir, "/", trait_files), read.csv)
labels <- unlist(lapply(strsplit(trait_files, "_"), "[[", 1))
names(traits) <- labels


for(group_index in 1:length(labels)) {
  group <- names(traits)[group_index]
  group_traits <- traits[[group]][,1:2] 
  group_tree <- trees[[grep(group, names(trees))]]
  
  to_drop <- c("remove","no_info","doubtful")
  if(any(to_drop %in% group_traits[,2])) {
    group_tree <- drop.tip(group_tree, group_traits[,1][which(group_traits[,2] %in% to_drop)])
    trees[[group_index]] <- group_tree # return "cleaned" tree
  }
}

#-------------------------------
# Download GBIF data
#-------------------------------

# fill in your gbif.org credentials 
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
wrld_map <- getMap(resolution="low") # leaving both maps in the code for now, should probably drop one of them later
data(wrld_simpl) 

# Download GBIF tables from your download page on GBIF

points.dir <- paste0(getwd(), "/GBIF_points")
point_files <- list.files(points.dir, ".csv")
points <- lapply(paste0(points.dir, "/", point_files), fread,quote="")
names(points) <- unlist(lapply(strsplit(point_files, "_"), "[[", 1))


# Full filtering 
family_index =4
for(family_index in 2:length(points)) {
points_cleaned <- points[[family_index]]
#------------------------------
# Removing collections with no species identification
#------------------------------
points_cleaned <- points_cleaned[-which(points_cleaned$species==""),]
#------------------------------
# Cleaning points in the sea 
#------------------------------
buffer = 2 # buffer in degrees
coords <- points_cleaned[,c("decimalLongitude","decimalLatitude")]
coordinates(coords) <- ~ decimalLongitude + decimalLatitude
proj4string(coords) <- proj4string(wrld_map)
country_plus_buffer <- buffer(wrld_map, buffer) # adding buffer around polygons
answer <- which(is.na(over(coords, country_plus_buffer)))
points_cleaned <- points_cleaned[-answer,]

#------------------------------
# Cleaning points by country 
#------------------------------
# That will remove points that are not located in the countries where they were
# said to be collected
buffer = 5 # buffer in degrees
countries <- as.character(wrld_simpl[2]$ISO2)
dubiousGBIF_ids <- c()
for(country_index in 1:length(countries)) {
  tmp_country <- countries[country_index]
  if(length(which(points_cleaned$countryCode %in% tmp_country)) > 0) {
    tmp_subset <- points_cleaned[points_cleaned$countryCode==tmp_country,]
    coords <- tmp_subset[,c("decimalLongitude","decimalLatitude")]
    coordinates(coords) <- ~ decimalLongitude + decimalLatitude
    proj4string(coords) <- proj4string(wrld_simpl)
    country_plus_buffer <- buffer(wrld_simpl[country_index,], buffer) # adding buffer of 1 degree around country
    answer <- which(is.na(over(coords, country_plus_buffer)))
    dubiousGBIF_ids <- c(dubiousGBIF_ids, tmp_subset$gbifID[answer])
  }
}
points_cleaned <- points_cleaned[-which(points_cleaned$gbifID %in% dubiousGBIF_ids),]

#------------------------------
# Removing centroids
#------------------------------
buffer = 75000 # here the buffer is in meters
# Note: probably should do something about small countries where the buffer of 50km may be too large
coords <- points_cleaned[,c("decimalLongitude","decimalLatitude")]
coordinates(coords) <- ~ decimalLongitude + decimalLatitude
proj4string(coords) <- proj4string(wrld_map)
centroids <- gCentroid(wrld_map, byid=TRUE)
centroids_plus_buffer <- buffer(centroids, buffer) # adding buffer around polygons
answer <- which(!is.na(over(coords, centroids_plus_buffer)))
points_cleaned <- points_cleaned[-answer,]

#------------------------------
# Removing duplicates
#------------------------------
spp <- unique(points_cleaned$species)
all_points <- list()
for(species_index in 1:length(spp)){
  tmp_subset <- as.data.frame(points_cleaned[points_cleaned$species==spp[species_index],])
  all_points[[species_index]] <- tmp_subset[-which(duplicated(tmp_subset[,c("decimalLongitude","decimalLatitude")])),]
  print(species_index)
}
points_cleaned <- do.call(rbind, all_points)

#------------------------------
# Removing points with no decimal cases
#------------------------------
length_decimal_lat <- nchar(sub("^[^.]*", "", points_cleaned$decimalLatitude))
length_decimal_lon <- nchar(sub("^[^.]*", "", points_cleaned$decimalLongitude))
points_cleaned <- points_cleaned[which(length_decimal_lat>=1 & length_decimal_lon>=1),]

#------------------------------
# Removing invasive and crop species
#------------------------------
invasive <- readRDS(paste0(points.dir, "/invasive_taxized_gbif.Rdata"))
points_cleaned <- points_cleaned[-which(points_cleaned$species %in% invasive), ]
crops <- readRDS(paste0(points.dir, "/crops_taxized_gbif.Rdata"))
points_cleaned <- points_cleaned[-which(points_cleaned$species %in% crops), ]

#------------------------------
# Removing outlier points
#------------------------------
spp <- unique(points_cleaned$species)
all_points <- list()
for(species_index in 1:length(spp)){
  sp0 <- points_cleaned[points_cleaned$species==spp[species_index],]
  out_lat <- boxplot.stats(sp0$decimalLatitude)$out
  out_lon <- boxplot.stats(sp0$decimalLongitude)$out
  sp0 <- sp0[!sp0$decimalLatitude %in% out_lat, ]
  sp0 <- sp0[!sp0$decimalLongitude %in% out_lon, ]
  all_points[[species_index]] <- sp0
  print(species_index)
}
points_cleaned <- do.call(rbind, all_points)
write.csv(points_cleaned, paste0(points.dir, "/", names(points)[family_index], ".csv"))

# Ploting
spp <- unique(points_cleaned$species)
pdf(paste0(points.dir, "/", names(points)[family_index], "_points.pdf"))
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


#----------------------------------- to be continued
#------------------------------
# Thinning to smooth sampling bias
#------------------------------
spp <- unique(points_cleaned$species)
all_points <- list()
for(species_index in 1:length(spp)){
  sp0 <- points_cleaned[points_cleaned$species==spp[species_index],]
  coords <- sp0[,c("decimalLatitude","decimalLongitude")]
  coordinates(coords) <- ~ decimalLatitude + decimalLongitude
  crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  r0 <- raster(coords)
  res(r0) <- 0.5 # cell resolution
  r0 <- extend(r0, extent(r0) + 1) # expand the extent of the RasterLayer a little
  thinned_points <- as.data.frame(dismo::gridSample(coords, r0, n = 1)) # n = maximum number of points per cell
  #species <- points_for_range$species[sequence(nrow(thinned_points))]
  #thinned_points <- cbind(species, thinned_points)
}


# https://datahub.io/core/country-list#data
# https://datahub.io/JohnSnowLabs/country-and-continent-codes-list

#country_code <- read.csv(paste0(points.dir, "/country-and-continent-codes-list.csv"))
#country_code[which(is.na(country_code$Continent_Code)),] <- "NA"
#country_code[which(is.na(country_code$Two_Letter_Country_Code)),] <- "NA"

