## functions for getting climatic data

#' Taxize GBIF
#' @param name A character vector with species names.
#' @export
#' @importFrom taxize gnr_datasources gnr_resolve
#' @importFrom pbapply pbsapply
resolveGBIF <- function(name) {
  gnr_resolve_x <- function(x) {
    sources <- taxize::gnr_datasources()
    tmp.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "GBIF Backbone Taxonomy"], best_match_only=TRUE)$matched_name)
    if(is.null(tmp.name)) {
      tmp.name <- paste0("UNMATCHED_",x)
    }
    return(tmp.name)
  }
  new.names <- pbapply::pbsapply(name, gnr_resolve_x)
  return(as.character(new.names))
}

#' Removes points in the sea
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number of degrees around continental areas where points are still kept after cleaning
#' @importFrom rworldmap getMap
#' @importFrom sp coordinates proj4string over
#' @importFrom raster buffer
#' @export
RemoveSeaPoints <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=0) {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  wrld_map <- rworldmap::getMap(resolution="low") # leaving both maps in the code for now, should probably drop one of them later
  coords <- tmp_points[,c("x","y")]
  sp::coordinates(coords) <- ~ x + y
  sp::proj4string(coords) <- sp::proj4string(wrld_map)
  country_plus_buffer <- raster::buffer(wrld_map, buffer) # adding buffer around polygons
  answer <- which(is.na(sp::over(coords, country_plus_buffer)))
  if(length(answer) > 0) {
    points <- points[-answer,]
    npoints_end <- nrow(points)
    print(paste0(npoints_start - npoints_end, " points removed."))
    return(points)
  } else {
    print("no points removed")
    return(points) }
}

#' Removes points that have 0 for both latitude and longitude
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @export
RemoveZeros <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
  if(!inherits(points, "data.frame")) {
    stop("Argument points is not a data.frame.")
  }
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  if(any(tmp_points$x==0 & tmp_points$y==0)) {
    points <- points[-which(tmp_points$x==0 & tmp_points$y==0),]
  }
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}

#' Removes outliers
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param species The name of the column in the data.frame with the names of species
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @export
#' @importFrom grDevices boxplot.stats
RemoveOutliers <- function(points, species="species", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  colnames(tmp_points)[colnames(tmp_points)==species] <- "species"
  spp <- unique(tmp_points$species)
  all_points <- list()
  for(species_index in 1:length(spp)){
    sp0 <- tmp_points[tmp_points$species==spp[species_index],]
    out_lat <- grDevices::boxplot.stats(sp0$y)$out
    out_lon <- grDevices::boxplot.stats(sp0$x)$out
    sp0 <- sp0[!sp0$y %in% out_lat, ]
    sp0 <- sp0[!sp0$x %in% out_lon, ]
    all_points[[species_index]] <- sp0
  }
  points <- do.call(rbind, all_points)
  colnames(points)[colnames(points)=="x"] <- lon
  colnames(points)[colnames(points)=="y"] <- lat
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}

#' Removes points that are located in the wrong country according to their GBIF labels
#' That will remove points that are not located in the countries where their labels say they were
#' collected
#' @param points A data.frame of distribution points with at least five columns where one column represents species names and other two decimal coordinates.
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number of degrees around each country where points are still consider part of that country
#' @details For now, the input data.frame must have a column named countryCode and one named gbifID, as the .csv files downloaded directly from GBIF.
#' @importFrom sp coordinates proj4string over
#' @importFrom raster buffer
#' @importFrom stats na.omit
#' @export
RemoveWrongCountries <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=5) {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  wrld_simpl <- readRDS("R/wrld_simpl.Rdata")
  countries <- as.character(wrld_simpl[2]$ISO2)
  dubiousGBIF_ids <- c()
  for(country_index in 1:length(countries)) {
    tmp_country <- countries[country_index]
    if(length(which(tmp_points$countryCode %in% tmp_country)) > 0) {
      tmp_subset <- tmp_points[tmp_points$countryCode==tmp_country,]
      coords <- stats::na.omit(tmp_subset[,c("x","y")])
      sp::coordinates(coords) <- ~ x + y
      sp::proj4string(coords) <- sp::proj4string(wrld_simpl)
      country_plus_buffer <- raster::buffer(wrld_simpl[country_index,], buffer) # adding buffer of 1 degree around country
      answer <- which(is.na(sp::over(coords, country_plus_buffer)))
      dubiousGBIF_ids <- c(dubiousGBIF_ids, tmp_subset$gbifID[answer])
    }
  }
  points_cleaned <- points[-which(points$gbifID %in% dubiousGBIF_ids),]
  npoints_end <- nrow(points_cleaned)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points_cleaned)
}

#' Removes points that are located in country centroids
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number in meters around each country centroid for points to be removed
#' @importFrom sp coordinates proj4string over
#' @importFrom raster buffer
#' @importFrom rgeos gCentroid
#' @export
RemoveCentroids <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=75000) {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  wrld_map <- rworldmap::getMap(resolution="low")
  # here the buffer is in meters
  # Note: probably should do something about small countries where the buffer of 75km may be too broad
  coords <- tmp_points[,c("x","y")]
  sp::coordinates(coords) <- ~ x + y
  sp::proj4string(coords) <- sp::proj4string(wrld_map)
  centroids <- rgeos::gCentroid(wrld_map, byid=TRUE)
  centroids_plus_buffer <- raster::buffer(centroids, buffer) # adding buffer around polygons
  answer <- which(!is.na(sp::over(coords, centroids_plus_buffer)))
  if(length(answer) > 0) {
    points <- points[-answer,]
    npoints_end <- nrow(points)
    print(paste0(npoints_start - npoints_end, " points removed."))
    return(points)
  } else {
    print("no points removed")
    return(points) }
}


#' Removes duplicated species, latitudes and longitudes keeping only one of the duplicated values per species
#' @param points A data.frame of distribution points with at least five columns where one column represents species names and other two decimal coordinates.
#' @param species The name of the column in the data.frame with the names of species
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @export
RemoveDuplicates <- function(points, species="species", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  colnames(tmp_points)[colnames(tmp_points)==species] <- "species"
  spp <- unique(tmp_points$species)
  all_points <- list()
  for(species_index in 1:length(spp)){
    tmp_subset <- as.data.frame(tmp_points[tmp_points$species==spp[species_index],])
    all_points[[species_index]] <- tmp_subset[-which(duplicated(tmp_subset[,c("x","y")])),]
  }
  points <- do.call(rbind, all_points)
  colnames(points)[colnames(points)=="x"] <- lon
  colnames(points)[colnames(points)=="y"] <- lat
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}


#' Removes points with coordinates without decimal cases (probably innacurate)
#' @param points A data.frame of distribution points with at least five columns where one column represents species names and other two decimal coordinates.
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @export
RemoveNoDecimal <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  length_decimal_lat <- nchar(sub("^[^.]*", "", tmp_points$y))
  length_decimal_lon <- nchar(sub("^[^.]*", "", tmp_points$x))
  points <- points[which(length_decimal_lat>=1 & length_decimal_lon>=1),]
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}


#' Removes points around large herbaria.
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and other two decimal coordinates.
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @export
RemoveHerbariaLocalities <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  herbaria_localities <- readRDS("R/allHerbaria_ADM1_badCoords.Rdata") # Jeremy provided this list
  # herbaria_localities <- read.csv("~/Desktop/MiSSEgradient/MiSSEGradient/Sampling/data/allHerbaria_ADM1_badCoords.txt", header=FALSE)
  for(herb.index in 1:length(herbaria_localities)){
    points <- tmp_points[which(!round(tmp_points$y,2)==herbaria_localities[herb.index,1] | !round(tmp_points$x,2)==herbaria_localities[herb.index,2]),]
  }
  colnames(points)[colnames(points)=="x"] <- lon
  colnames(points)[colnames(points)=="y"] <- lat
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}


#' Flag species that have weird distributions in the dataset
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and other two decimal coordinates.
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @export
UnusualDistributions <- function (points, species="species", lat="decimalLatitude", lon="decimalLongitude", buffer.polygon=c(5,10)) {
  tmp_points = as.data.frame(points)
  tmp_points = tmp_points[,c(which(colnames(tmp_points)==species),which(colnames(tmp_points)==lon),which(colnames(tmp_points)==lat))]
  colnames(tmp_points) <- c("species","lon","lat")
  spp <- unique(tmp_points$species)
  
  # WWFload is taken from speciesgeocodeR:
  WWFload <- function(x = NULL) {
    if (missing(x)) {
      x <- getwd()
    }
    download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                  destfile = file.path(x, "wwf_ecoregions.zip"), quiet=TRUE)
    unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
    file.remove(file.path(x, "wwf_ecoregions.zip"))
    wwf <- maptools::readShapeSpatial(file.path(x, "WWF_ecoregions", "official",
                                                "wwf_terr_ecos.shp"))
    return(wwf)
  }
  wwf <- WWFload(tempdir())
  result <- matrix(nrow=length(spp), ncol=2)
  for(species_index in 1:length(spp)) {
    one_species <- tmp_points[tmp_points$species==spp[species_index],]
    locations.spatial <- sp::SpatialPointsDataFrame(coords=one_species[,c("lon", "lat")], data=one_species)
    mappedregions <- sp::over(locations.spatial, wwf)
    
    print(species_index)
    if(length(which(table(mappedregions$REALM)!=0)) > 2) {
      result[species_index,1] <- spp[species_index]
      result[species_index,2] <- "flagged"
    } else {
      result[species_index,1] <- spp[species_index]
      result[species_index,2] <- "ok"
      
    }
  }
  
  flagged <- result[,1][which(result[,2]=="flagged")]
  answer <- c()
  for(flagged_index in 1:length(flagged)) {
    print(flagged_index)
    subset <- points[points$species==flagged[flagged_index],]
    subset1 <- subset[,c("decimalLongitude","decimalLatitude")]
    distances <- as.matrix(dist(subset1))
    mean_distances <- apply( distances, 1, function(x) mean( x[order(x)][2:4] ) )
    outliers <- subset[which(mean_distances > mean(distances)),]
    
    if(nrow(outliers)>0){
      subset <- subset[-which(mean_distances > mean(distances)),]
    }
    max.lat <- ceiling(max(subset[,"decimalLatitude"])) + buffer.polygon[1]
    min.lat <- floor(min(subset[,"decimalLatitude"])) - buffer.polygon[1]
    max.lon <- ceiling(max(subset[,"decimalLongitude"])) + buffer.polygon[2]
    min.lon <- floor(min(subset[,"decimalLongitude"])) - buffer.polygon[2]
    if(any(c(max.lat > min.lat + 75, max.lon > min.lon + 125))){
      answer <- c(answer,flagged[flagged_index])
    }
  }
  return(answer)
}


#' Extract climate information from species points.
#' @description Function to extract climate information from species points.
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param res A number (2.5, 5 or 10) indicating resolution of climatic layers from where the climate data is being extracted.
#' @importFrom raster extract getData addLayer
#' @importFrom sp coordinates
#' @export
ClimateFromPoints <- function (points, species="", lat = "lat", lon="lon", res=5) {
  tmp_points = points
  colnames(tmp_points)[which(colnames(tmp_points) == lon)] <- "lon"
  colnames(tmp_points)[which(colnames(tmp_points) == lat)] <- "lat"
  colnames(tmp_points)[which(colnames(tmp_points) == species)] <- "species"
  tmp_points <- tmp_points[,c("species","lat","lon")]
  
  bio <- raster::getData("worldclim", var="bio", res=res)
  alt <- raster::getData("worldclim", var="alt", res=res)
  bio <- raster::addLayer(bio, alt)
  
  vars <- c(names(bio))
  final_matrix <- matrix(nrow=nrow(tmp_points), ncol=length(vars))
  
  cat("Extracting climatic information of", nrow(tmp_points), "points",  "\n")
  sp::coordinates(tmp_points) <- ~ lon + lat
  for(var_index in 1:length(vars)) {
    layer <- bio[[which(names(bio)==vars[var_index])]]
    cat("\r",vars[var_index])
    cat("","\n")
    values <- raster::extract(layer, tmp_points)
    final_matrix[,var_index] <- values
  }
  colnames(final_matrix) <- vars
  result <- cbind(tmp_points, final_matrix)
  
  try(unlink("wc2-5", recursive = TRUE))
  try(unlink("wc5", recursive = TRUE))
  try(unlink("wc10", recursive = TRUE))
  
  return(as.data.frame(result))
}


#' Thinning distribution data to smooth sampling bias
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param n A number indicating how many points to keep in each cell after thinning
#' @importFrom raster extend extent
#' @importFrom dismo gridSample
#' @importFrom sp coordinates
#' @export
Thinning <- function(points, species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 1) {
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  spp <- unique(tmp_points[,species])
  results <- list()
  for(species_index in 1:length(spp)) {
    coords <- tmp_points[tmp_points[,species]==spp[species_index],c("y","x")]
    coords <- coords[!duplicated(coords[,"x"]) & !duplicated(coords[,"y"]),]
    if(nrow(coords) > 1) {
      sp::coordinates(coords) <- ~ y + x
      raster::crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      r0 <- raster::raster(coords)
      raster::res(r0) <- 1 # cell resolution
      r0 <- raster::extend(r0, raster::extent(r0) + 5) # expand the extent of the RasterLayer a little
      res <- cbind(spp[species_index], as.data.frame(dismo::gridSample(coords, r0, n))) # n = maximum number of points per cell
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    } else {
      res <- cbind(spp[species_index],coords)
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    }
  }
  results <- do.call(rbind, results)
  return(results)
}


##
GetClimateSummStats_seed_dispersal <- function (points) {
  tmp_points <- points[,-which(colnames(points) %in% c("lon","lat"))]
  vars <- c("bio1","bio12")
  allclimatevars <- list()
  spp <- unique(tmp_points$species)
  for(var_index in 1:length(vars)) {
    cat("\r",vars[var_index])
    cat("","\n")
    n_i <- c()
    sigma2_wi <- c()
    summ_stats <- matrix(nrow=length(spp), ncol=5)
    for(species_index in 1:length(spp)){
      sp1 <- tmp_points[tmp_points$species==spp[species_index],]
      cat("\r","Now doing species", species_index)
      cat("","\n")
      values <- sp1[,vars[var_index]]
      values <- values[!is.na(values)]
      if(vars[var_index] %in% c("bio1")) {
        values <-  (values / 10) + 273.15
      }
      
      values <- log(values)
      n_i[species_index] <- length(values) # sample size
      sigma2_wi[species_index] <- ifelse(is.na(var(values)),0,var(values))  # sample variance

      if(any(values== -Inf)){
        values <- values[-which(values== -Inf)]
      }
      n0 <- length(values)
      mean0 <- round(mean(values), 6)
      sd0 <- round(stats::sd(values), 6)
      se0 <- round(sd0/ sqrt(n0), 6)
      tmp_summ_stats <- c(n0, mean0, sd0, se0)
      summ_stats[species_index,] <- c(spp[species_index], tmp_summ_stats)
      colnames(summ_stats) <- c("species",paste0("n_",vars[var_index]), paste0("mean_",vars[var_index]),
                                paste0("sd_",vars[var_index]), paste0("se_",vars[var_index]))
    }
    sigma2_w <- sum(sigma2_wi*(n_i - 1)) / sum(n_i - 1)
    within_sp_var <-  round(sigma2_w/n_i, 6)
    summ_stats <- cbind(summ_stats, within_sp_var)
    colnames(summ_stats)[6] <- paste0("within_sp_var_",vars[var_index])
    allclimatevars[[var_index]] <- summ_stats
  }
  return(allclimatevars)
}


GetClimateSummStats_seed_dispersal <- function (points, type=c("raw","transformed")) {
  tmp_points <- points[,-which(colnames(points) %in% c("lon","lat"))]
  vars <- c("bio1","bio12")
  allclimatevars <- list()
  spp <- unique(tmp_points$species)
  for(var_index in 1:length(vars)) {
    cat("\r",vars[var_index])
    cat("","\n")
    n_i <- c()
    sigma2_wi <- c()
    summ_stats <- matrix(nrow=length(spp), ncol=5)
    for(species_index in 1:length(spp)){
      sp1 <- tmp_points[tmp_points$species==spp[species_index],]
      cat("\r","Now doing species", species_index)
      cat("","\n")
      values <- sp1[,vars[var_index]]
      values <- values[!is.na(values)]
      if(type=="raw") {
       if(vars[var_index] %in% c("bio1")) {
          values <-  (values / 10) 
        }
      n_i[species_index] <- length(values) # sample size
      sigma2_wi[species_index] <- ifelse(is.na(var(values)),0,var(values))  # sample variance
      
       if(any(values== -Inf)){
          values <- values[-which(values== -Inf)]
        }
      }
      if(type=="transformed") {
        if(vars[var_index] %in% c("bio1")) {
          values <-  (values / 10) + 273.15 # transform to Kelvin
        }
        
        values <- log(values) # log
        n_i[species_index] <- length(values) # sample size
        sigma2_wi[species_index] <- ifelse(is.na(var(values)),0,var(values))  # sample variance
        
        if(any(values== -Inf)){
          values <- values[-which(values== -Inf)]
        }
      }
      n0 <- length(values)
      mean0 <- round(mean(values), 6)
      sd0 <- round(stats::sd(values), 6)
      se0 <- round(sd0/ sqrt(n0), 6)
      tmp_summ_stats <- c(n0, mean0, sd0, se0)
      summ_stats[species_index,] <- c(spp[species_index], tmp_summ_stats)
      colnames(summ_stats) <- c("species",paste0("n_",vars[var_index]), paste0("mean_",vars[var_index]),
                                paste0("sd_",vars[var_index]), paste0("se_",vars[var_index]))
    }
    sigma2_w <- sum(sigma2_wi*(n_i - 1)) / sum(n_i - 1)
    within_sp_var <-  round(sigma2_w/n_i, 6)
    summ_stats <- cbind(summ_stats, within_sp_var)
    colnames(summ_stats)[6] <- paste0("within_sp_var_",vars[var_index])
    allclimatevars[[var_index]] <- summ_stats
  }
  return(allclimatevars)
}

