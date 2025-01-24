####---- Libraries ----####
library(devtools)
library(BioGeoBEARS)
library(Rcpp)
library(ape)
library(terra)
library(sf)

####---- Functions ----####

# gets pairwise minimum geographic distances (km) between regional polygons
getGeoDistances <- function(x, mergeSthAmericaAntarctica=FALSE, removeAntarctica=TRUE, removeNthAmerica=TRUE){
  
  if(mergeSthAmericaAntarctica==TRUE){
    
    SthAmerica<-  x$SthAmerica
    crs(SthAmerica) <- "+proj=longlat +datum=WGS84"
    Antarctica<-  x$Antarctica
    crs(Antarctica) <- "+proj=longlat +datum=WGS84"
    x$SthAmerica <- combineGeoms(SthAmerica, Antarctica)
    x <- x[which(names(x) != "Antarctica")]
  }
  if(removeAntarctica==TRUE){
    x <- x[which(names(x) != "Antarctica")]
  }
  if(removeNthAmerica==TRUE){
    x <- x[which(names(x) != "TemperateNorthAmerica")]
  }
  # get all the names
  polygons_names_all <- names(x)
  
  # remove missing
  suppressWarnings({x <- x[lapply(x,is.na)==0]})
  suppressWarnings({x <- x[lapply(x,is.null)==0]})
  polygons <- x
  
  # get subset names
  polygons_names <- names(polygons)
  
  #turn into a vect, set crs, and turn into sf
  polygons <- vect(polygons)
  polygons <- crop(polygons, ext(-179.99, 179.99, -90, 90))
  crs(polygons) <- "+proj=longlat +datum=WGS84"
  polygons <- st_as_sf(polygons)
  polygons <- st_make_valid(polygons)
  # get distances and name the matrix
  polygons_distances <- st_distance(polygons)
  
  colnames(polygons_distances) <- polygons_names
  rownames(polygons_distances) <- polygons_names
  
  # Get the combined set of elements from both matrices
  all_elements <- polygons_names_all
  
  # Create a new empty matrix with all elements
  new_polygons_distances <- matrix(NA, nrow = length(all_elements), ncol = length(all_elements), dimnames = list(all_elements, all_elements))
  
  # Copy existing distances from dist_matrix1 to new_dist_matrix1
  new_polygons_distances[rownames(polygons_distances), colnames(polygons_distances)] <- polygons_distances
  
  # return matrix in kms
  new_polygons_distances/1000
}

# rescales geographic distances between 0 and 10
scaleBGBDistanceMatrix <- function(x, upper=10){
  require(scales)
  x <- round(x, 2)
  x[which(x == 0)] <- 1
  x[which(x != 0)] <- scales::rescale(x[which(x != 0)], to=c(1, upper))
  x[which(is.na(x))] <- upper
  diag(x) <- 0
  x <- round(x, 2)
  colnames(x) <- NULL
  rownames(x) <- NULL
  return(x)
  
}

# Are polygons contiguous or do they have distance in betwee?
isContiguous <- function(x,  mergeSthAmericaAntarctica=TRUE, removeAntarctica=TRUE, removeNthAmerica=TRUE){
  
  contiguous <- vector("list", nrow(x))
  
  for(i in 1:nrow(x)){
    
    if(all(is.na(x[i,]))){contiguous[[i]] <- NA; next}
    
    connecting_names_1 <- colnames(x)[i]
    connecting_names_2 <- unique(colnames(x[,which(x[i,] ==0), drop=F]))
    x_tmp <- x[,which(x[i,] ==0)]
    
    while(!length(connecting_names_1) == length(connecting_names_2)){
      connecting_names_1 <- connecting_names_2
      x_tmp <- x[, names(which(apply(x_tmp, 1, function(row) any(row == 0)))), drop=F]
      connecting_names_2 <- colnames(x_tmp)
    }
    
    contiguous[[i]] <- connecting_names_2 
  }
  return(contiguous)
}

# Calculate connectivity between regions and assign values as high (contiguous), medium (seperated by moderate over-water distances), or low (seperated by vast over-water distances)
connectivityMatrix <- function(x, high, medium, low){
  
  new_x <- x
  new_x[] <- 0.01
  contiguous <- isContiguous(x)
  for(i in 1:length(contiguous)){
    new_x[which(x[i,] <= 3000), i] <- medium
    new_x[contiguous[[i]], i] <- high
    new_x[which(x[i,] > 3000), i] <- low
    
    print(new_x)
  }
  diag(new_x) <- 0
  
  return(new_x)
}


####---- Model Name ----####

# name model
# +w connectivity weights
# +n environmental weights
# +x geographic distance weights
# +j jump dispersal speciation
# s state space modified
model <- "timestrat_DEC_s_+w_+n_+x_+j"

####---- Biome Table ----####

# read in the biome data
koppen_pam <- readxl::read_xlsx("Dataset_S3.xlsx", sheet="Biome_Occupancy")

# change colnames to single alphabertical values
koppen_names <- colnames(koppen_pam)
colnames(koppen_pam) <- c("species", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")
# set rownames
rownames(koppen_pam) <- koppen_pam$species
koppen_pam <- koppen_pam[, 2:ncol(koppen_pam)]

# turn proportion of range values into presence/absence values (1/0) based on 10% threshold
koppen_mat <- as.matrix(koppen_pam)
koppen_mat[which(koppen_mat <= 0.1)] <- 0
koppen_mat[which(koppen_mat > 0.1)] <- 1

# make the BGB text file
bgb_header_txt <- c(nrow(koppen_mat), length(which(colSums(koppen_mat)>0)), paste("(",paste(colnames(koppen_mat), collapse=" "), ")", sep=""))
bgb_column_txt <- matrix(apply(koppen_mat, MARGIN=1, FUN=function(x){paste(x[1:length(x)], collapse="")}), ncol=1)
rownames( bgb_column_txt) <- rownames(koppen_mat)
write.table( bgb_column_txt, file="Outputs/koppen_regional_geog.txt", col.names=c(paste(bgb_header_txt, collapse="\t")), row.names=T, sep="\t",quote = FALSE, append=T)

####---- Ranges List ----####

# get a list of allowed states in each timeslice based on the paleobiome reconstructions
max_range_size = max(rowSums(koppen_mat))
areas = colnames(koppen_mat)

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = cladoRcpp::rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# How many states/ranges, by default: 232
length(states_list_0based)

# Make the list of ranges
states = NULL
for (i in 1:length(states_list_0based)){    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) ){
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  states = c(states, tmprange)
}

# How many states/ranges, by default: 232
length(states)

states_000 <- states_list_0based
states_020 <- states_list_0based[-c(which(sapply(c("e"), grepl, states)))]
states_040 <- states_list_0based[-c(which(sapply(c("e"), grepl, states)))]
states_060 <- states_list_0based[-c(which(sapply(paste0(c("e", "d"), collapse="|"), grepl, states)))]
states_080 <- states_list_0based[-c(which(sapply(paste0(c("e", "d"), collapse="|"), grepl, states)))]
states_100 <- states_list_0based[-c(which(sapply(paste0(c("e", "d"), collapse="|"), grepl, states)))]
states_120 <- states_list_0based[-c(which(sapply(paste0(c("a","c", "e", "d"), collapse="|"), grepl, states)))]

ranges_list  <- list(states_000, states_020,
                     states_040, states_060,
                     states_080, states_100)


####---- Time Periods ----####

# 20 Ma intervals from 120 - 0 (we don';'t include 0 explicitly in this object however)
timeperiods <- matrix(seq(from=20, to=120, by=20), ncol=1)
write.table(timeperiods,"Outputs/timeperiods_table_120Ma.txt", row.names=F,col.names=F, sep="\t", quote = FALSE)

####----Environmental Distances ----####

# koppen data from Beck et al 
koppen <- rast("External_datasets/Beck_KG_V1_present_0p5.tif")

# climate data from CHELSA
chelsa <- rast("External_datasets/chelsa_stack.tif")

# get distribution of species present in the phylogeny
richness <- rast("External_datasets/Grevilleoideae_richness_phylospecies.tif")
richness <- resample(richness , koppen)

# combine
climate_df <- as.data.frame(c(koppen, chelsa, richness), xy=TRUE)
climate_df <- climate_df[which(climate_df$sum > 0), ]

# Remove North America (which has Proteales outgroups - Platanus and Nelumbo)
climate_df <- climate_df[-which(climate_df$y > 20 & climate_df$x < 0), ]

# Combine biomes to new modification scheme
climate_df$biome <- NA
climate_df$biome[which(climate_df[,3] %in% c(0))] <- NA
climate_df$biome[which(climate_df[,3] %in% c(1, 2, 3, 11))] <- "Tropical"
climate_df$biome[which(climate_df[,3] %in% c(8, 9, 6))] <- "Mediterranean"
climate_df$biome[which(climate_df[,3] %in% c(4,5))] <- "Arid"
climate_df$biome[which(climate_df[,3] %in% c(6,7))] <- "SemiArid"
climate_df$biome[which(climate_df[,3] %in% c(14, 15, 29, 16, 12, 19, 10,18,  27, 25, 26))] <- "Subtropical"

# define regions based on Lat/Lons
climate_df$biome[which(climate_df$x < 0 &  climate_df$y < 20  & climate_df$biome %in% c("Subtropical", "Mediterranean","Tropical"))] <- "SthAmerica"
climate_df$biome[which(climate_df$x < 75 & climate_df$x > 0 & climate_df$y < -30 )] <- "Cape"
climate_df$biome[which(climate_df$x < 75 & climate_df$x > 0 & climate_df$y > -30 )] <- "Madagascar"
climate_df$biome[which(climate_df$x > 160   & climate_df$y < -30)]  <- "NewZealand"
climate_df$biome[which(climate_df$x > 162   & climate_df$y > -30)]  <- "NewCaledonia"
climate_df$biome[which(climate_df$x > 0 & climate_df$x < 130 & climate_df$y > 0)] <- "TropicalAsia"

# Log-transform preciptation variables
climate_df[, c("CHELSA_bio10_12", "CHELSA_bio10_13", "CHELSA_bio10_14", 
               "CHELSA_bio10_16", "CHELSA_bio10_17", "CHELSA_bio10_18", 
               "CHELSA_bio10_19")] <- log(climate_df[, c("CHELSA_bio10_12", "CHELSA_bio10_13", "CHELSA_bio10_14", 
                                                         "CHELSA_bio10_16", "CHELSA_bio10_17", "CHELSA_bio10_18", 
                                                         "CHELSA_bio10_19")]+0.00001)
# remove NAs
climate_df <- na.omit(climate_df)

# scale values
scale_climate_df <- scale(climate_df[,4:26])
pca <- prcomp(scale_climate_df)

# PCA centroids
pc_axes <- as.data.frame(pca$x)
pc_axes$biome <- climate_df$biome
pca_centroids <- aggregate(pc_axes[,1:23], list(biome = pc_axes$biome), mean)

# get the environmental distances!
enviro_dist <- as.matrix(dist(pca_centroids))
colnames(enviro_dist) <- pca_centroids$biome
rownames(enviro_dist) <- pca_centroids$biome
enviro_dist <- enviro_dist[c("Tropical", "Subtropical","Mediterranean",   "SemiArid",  "Arid", "SthAmerica",
                             "Madagascar", "Cape", "TropicalAsia" ,"NewCaledonia", "NewZealand"), 
                           c("Tropical", "Subtropical","Mediterranean",   "SemiArid",  "Arid", "SthAmerica",
                             "Madagascar", "Cape", "TropicalAsia" ,"NewCaledonia", "NewZealand")]

# rescale between 0 and 10 (for)
enviro_dist <- round(scales::rescale(as.matrix(enviro_dist), to=c(0,10)), 2)
colnames(enviro_dist) <- NULL
rownames(enviro_dist) <- NULL

# Non-time-stratified dists
write.table(enviro_dist,"Outputs/environmental_distances_table.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"),sep="\t", row.names=F, quote = FALSE)
write.table("","Outputs/environmental_distances_table.txt",append=T, col.names=F,sep="\t", row.names=F, quote = FALSE)
write.table("END","Outputs/environmental_distances_table.txt",append=T, col.names=F,sep="\t", row.names=F, quote = FALSE)

# Time-stratified 120 Ma (here, we just assume the distances remain the same at each timestep)
file.remove("Outputs/environmental_distances_table_timestratified_120Ma.txt")
write.table(enviro_dist,"Outputs/environmental_distances_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), sep="\t",row.names=F,  quote = FALSE)
for(i in 1:6){
  write.table(" ","Outputs/environmental_distances_table_timestratified_120Ma.txt", col.names=F, append=T,row.names=F, sep="\t", quote = FALSE)
  write.table(enviro_dist,"Outputs/environmental_distances_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"),append=T, row.names=F, sep="\t", quote = FALSE)
}
write.table("","Outputs/environmental_distances_table_timestratified_120Ma.txt", col.names=F,append=T,sep="\t", row.names=F, quote = FALSE)
write.table("END","Outputs/environmental_distances_table_timestratified_120Ma.txt", col.names=F,append=T,sep="\t", row.names=F, quote = FALSE)

####---- Geographic Distances ----####

# load in the paleobiome rasters and polygons of each region
# this takes a few minutes to load
source("Software/PALEOBIOMES.R")

# merge Polar Sth America / Antarctica as they were connected
geo_dists_120MA <- getGeoDistances(RegionalPolygons[[1]], mergeSthAmericaAntarctica=TRUE)
geo_dists_100MA <- getGeoDistances(RegionalPolygons[[2]], mergeSthAmericaAntarctica=TRUE)
geo_dists_080MA <- getGeoDistances(RegionalPolygons[[3]], mergeSthAmericaAntarctica=TRUE)
geo_dists_060MA <- getGeoDistances(RegionalPolygons[[4]], mergeSthAmericaAntarctica=TRUE)
geo_dists_040MA <- getGeoDistances(RegionalPolygons[[5]], mergeSthAmericaAntarctica=TRUE)
geo_dists_020MA <- getGeoDistances(RegionalPolygons[[6]], mergeSthAmericaAntarctica=FALSE)
geo_dists_000MA <- getGeoDistances(RegionalPolygons[[7]], mergeSthAmericaAntarctica=FALSE)


# Non Time-Stratified Geographic Distances
file.remove("Outputs/geographic_distances_table.txt")
geo_dists_000MA <- scaleBGBDistanceMatrix(geo_dists_000MA)
write.table(geo_dists_000MA,"Outputs/geographic_distances_table.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=F,row.names=F, sep="\t", quote = FALSE)
write.table("","Outputs/geographic_distances_table.txt", col.names=F,append=T,sep="\t", row.names=F, quote = FALSE)
write.table("END","Outputs/geographic_distances_table.txt", col.names=F,append=T,sep="\t", row.names=F, quote = FALSE)

# Time Stratified Geographic Distances
file.remove("Outputs/geographic_distances_table_timestratified_120Ma.txt")
write.table(geo_dists_000MA,"Outputs/geographic_distances_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=F,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/geographic_distances_table_timestratified_120Ma.txt", append=T, col.names=F,row.names=F, sep="\t", quote = FALSE)
geo_dists_020MA <- scaleBGBDistanceMatrix(geo_dists_020MA)
write.table(geo_dists_020MA,"Outputs/geographic_distances_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/geographic_distances_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
geo_dists_040MA <- scaleBGBDistanceMatrix(geo_dists_040MA)
write.table(geo_dists_040MA,"Outputs/geographic_distances_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/geographic_distances_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
geo_dists_060MA <- scaleBGBDistanceMatrix(geo_dists_060MA)
write.table(geo_dists_060MA,"Outputs/geographic_distances_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/geographic_distances_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
geo_dists_080MA <- scaleBGBDistanceMatrix(geo_dists_080MA)
write.table(geo_dists_080MA,"Outputs/geographic_distances_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/geographic_distances_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
geo_dists_100MA <- scaleBGBDistanceMatrix(geo_dists_100MA)
write.table(geo_dists_100MA,"Outputs/geographic_distances_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/geographic_distances_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)

write.table("","Outputs/geographic_distances_table_timestratified_120Ma.txt", col.names=F,append=T,sep="\t", row.names=F, quote = FALSE)
write.table("END","Outputs/geographic_distances_table_timestratified_120Ma.txt", col.names=F,append=T,sep="\t", row.names=F, quote = FALSE)

####---- Connectivity Matrices ----####

# merge Polar Sth America / Antarctica as they were connected
disp_probs_120MA <- connectivityMatrix(getGeoDistances(RegionalPolygons[[1]], mergeSthAmericaAntarctica=TRUE), high=1, medium=0.75, low=0.25)
disp_probs_100MA <- connectivityMatrix(getGeoDistances(RegionalPolygons[[2]], mergeSthAmericaAntarctica=TRUE), high=1, medium=0.75, low=0.25)
disp_probs_080MA <- connectivityMatrix(getGeoDistances(RegionalPolygons[[3]], mergeSthAmericaAntarctica=TRUE), high=1, medium=0.75, low=0.25)
disp_probs_060MA <- connectivityMatrix(getGeoDistances(RegionalPolygons[[4]], mergeSthAmericaAntarctica=TRUE), high=1, medium=0.75, low=0.25)
disp_probs_040MA <- connectivityMatrix(getGeoDistances(RegionalPolygons[[5]], mergeSthAmericaAntarctica=TRUE), high=1, medium=0.75, low=0.25)
disp_probs_020MA <- connectivityMatrix(getGeoDistances(RegionalPolygons[[6]], mergeSthAmericaAntarctica=FALSE), high=1, medium=0.75, low=0.25)
disp_probs_000MA <- connectivityMatrix(getGeoDistances(RegionalPolygons[[7]], mergeSthAmericaAntarctica=FALSE), high=1, medium=0.75, low=0.25)

file.remove("Outputs/connectivity_table_timestratified_120Ma.txt")
write.table(disp_probs_000MA,"Outputs/connectivity_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=F,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/connectivity_table_timestratified_120Ma.txt", append=T, col.names=F,row.names=F, sep="\t", quote = FALSE)
write.table(disp_probs_020MA,"Outputs/connectivity_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/connectivity_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
write.table(disp_probs_040MA,"Outputs/connectivity_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/connectivity_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
write.table(disp_probs_060MA,"Outputs/connectivity_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/connectivity_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
write.table(disp_probs_080MA,"Outputs/connectivity_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/connectivity_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
write.table(disp_probs_100MA,"Outputs/connectivity_table_timestratified_120Ma.txt", col.names=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"), append=T,row.names=F, sep="\t", quote = FALSE)
write.table(" ","Outputs/connectivity_table_timestratified_120Ma.txt", append=T,col.names=F,row.names=F, sep="\t", quote = FALSE)
write.table("","Outputs/connectivity_table_timestratified_120Ma.txt", col.names=F,append=T,sep="\t", row.names=F, quote = FALSE)
write.table("END","Outputs/connectivity_table_timestratified_120Ma.txt", col.names=F,append=T,sep="\t", row.names=F, quote = FALSE)


####---- Run a time-stratified DEC on modified Koppen-Geiger biomes with Modified States Space  ----####

time_stratified <- TRUE
state_space_modified <- TRUE
plus_w <- TRUE
plus_n <- TRUE
plus_x <- TRUE
plus_j <- TRUE

# get geography and phylogeny paths
geogfn   <- file.path("Outputs/koppen_regional_geog.txt")
trfn <- file.path( "Dataset_S10.tre")
tr <- read.tree(trfn)

# Max range
max_range_size <- 3
n_cores <- 10 # number of cores to use - change to suit machine

# set up a BGB object
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$force_sparse=FALSE    # sparse=FALSE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$speedup=TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE    # get ancestral states from optim run

if(state_space_modified == TRUE){
  BioGeoBEARS_run_object$lists_of_states_lists_0based[[1]] = ranges_list[[1]]
  BioGeoBEARS_run_object$lists_of_states_lists_0based[[2]] = ranges_list[[2]]
  BioGeoBEARS_run_object$lists_of_states_lists_0based[[3]] = ranges_list[[3]]
  BioGeoBEARS_run_object$lists_of_states_lists_0based[[4]] = ranges_list[[4]]
  BioGeoBEARS_run_object$lists_of_states_lists_0based[[5]] = ranges_list[[5]]
  BioGeoBEARS_run_object$lists_of_states_lists_0based[[6]] = ranges_list[[6]]
}

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing if desired
BioGeoBEARS_run_object$num_cores_to_use=n_cores

# Sparse Matrix? Speeds up big analyses but prone to bugs/errors
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

if(time_stratified == TRUE){
  # Set up a time-stratified analysis
  BioGeoBEARS_run_object$timesfn = "Outputs/timeperiods_table_120Ma.txt"
  
}

if(time_stratified == TRUE & plus_w == TRUE){
  # add dispersal probability multipliers
  BioGeoBEARS_run_object$dispersal_multipliers_fn = "Outputs/connectivity_table_timestratified_120Ma.txt"
}

if(time_stratified == TRUE & plus_x == TRUE){
  # add geographic distances
  BioGeoBEARS_run_object$distsfn = "Outputs/geographic_distances_table_timestratified_120Ma.txt"
}

if(time_stratified == TRUE & plus_n == TRUE){
  # add environmental distances
  BioGeoBEARS_run_object$envdistsfn = "Outputs/environmental_distances_table_timestratified_120Ma.txt"
}

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

if(time_stratified == TRUE & plus_w == TRUE){
  # freely estimate w for dispersal probabilities
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","type"] <- "free"
}

if(time_stratified == TRUE & plus_x == TRUE){
  # freely estimate x for geographic distances
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] <- "free"
}

if(time_stratified == TRUE & plus_n == TRUE){
  # freely estimate n for environmental distances
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["n","type"] <- "free"
}

if(plus_j == TRUE){
  # Add j as a free parameter
  jstart <- 0.0001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
}

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# set output path
resfn = paste0("Outputs/", model,"/", model, ".Rdata")

#  Run the analysis
res = try(bears_optim_run(BioGeoBEARS_run_object))
save(res, file=resfn)


####---- Stochastic Mapping ----####

#  Run the analysis
BSM_inputs_fn = paste0("Outputs/", model,"/", model, "_BSM_inputs_file.Rdata") 
stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)

# Run mapping
BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, 
                    savedir = paste0("Outputs/", model),
                    maxnum_maps_to_try=100, 
                    nummaps_goal=50, 
                    maxtries_per_branch=40000, 
                    save_after_every_try=TRUE,  
                    seedval=12345, 
                    wait_before_save=0.01)




