####---- Libraries ----####
library(ncdf4)
library(terra)

####---- Functions ----####

# Function to swap northern hemisphere summer with winter and visa-versa so we can have a global "summer" raster and a global "winter" raster
swapHemispheres <- function(north_summer, south_summer){
  template_north_summer <- north_summer
  template_mat_north_summer <- as.data.frame(template_north_summer, xy=T, na.rm=F)
  mat_all_summer <- template_mat_north_summer
  
  template_south_summer <- south_summer
  template_mat_south_summer <- as.data.frame(template_south_summer, xy=T, na.rm=F)
  mat_all_winter <- template_mat_south_summer
  
  mat_all_summer[which(mat_all_summer$y <0), 3:ncol(mat_all_summer)] <- template_mat_south_summer[which(template_mat_south_summer$y <0), 3:ncol(template_mat_south_summer)] # all summer
  mat_all_winter[which(mat_all_winter$y <0), 3:ncol(mat_all_winter)] <- template_mat_north_summer[which(template_mat_north_summer$y <0), 3:ncol(template_mat_north_summer)] # replacing southern summer with southern values during northern summer
  
  ras_s <- rast(mat_all_summer, type="xyz")
  ras_w <- rast(mat_all_winter, type="xyz")
  
  return(list(ras_s, ras_w))
}

# function to make longitude between -180 and 180, rather than 0 and 360
readjustLon <- function(ras, t=F){
  template_raster <- ras
  if(t==T){
    template_mat <- as.data.frame(template_raster, xy=T, na.rm=F)
    colnames(template_mat) <- c("y", "x", colnames(template_mat)[3])
    template_mat <- template_mat[, c(2,1,3)]
  } else{
    template_mat <- as.data.frame(template_raster, xy=T, na.rm=F)
    
  }
  template_mat$x[which(template_mat$x >180)] <- template_mat$x[which(template_mat$x >180)] - 360
  ras_fix <- rast( template_mat, type="xyz")
  return(ras_fix)
}

# function to turn monthly temp and prec values from Li et al. into modified Koppen-Geiger biomes
koppenizerLi <- function(year,li_year, li_etal_rast, scotese){
  
  #~~ 1. Elevation Values
  scotese_ras <- rast(read.csv(file.path("../Paleoclimate/PaleoDEMS_long_lat_elev_csv_v2", "PaleoDEMS_long_lat_elev_csv_v2.csv",scotese[year] )), type="xyz")
  crs(scotese_ras) <- "+init=EPSG:4326"
  
  #~~ 1. Monthly TEMP and PREC values
  # Monthly
  li_etal_rast_i <- li_etal_rast[[which(grepl(paste0("simulation=", li_year),names(li_etal_rast)))]]
  crs(li_etal_rast_i) <- "+init=EPSG:4326"
  li_etal_rast_i <- project(li_etal_rast_i, scotese_ras) 
  monthlies <- names(li_etal_rast_i)
  temp_months <- vector("list", 12)
  prec_months <- vector("list", 12)
  for(month in 0:11){
    temp_months[[month+1]] <- li_etal_rast_i[[grepl(paste0("T_month=", month), monthlies)]]
    prec_months[[month+1]] <- li_etal_rast_i[[grepl(paste0("P_month=", month), monthlies)]]
    ext(temp_months[[month+1]]) <- ext(li_etal_rast_i)
    ext(prec_months[[month+1]]) <- ext(li_etal_rast_i)
  }
  
  temp_months <- rast( temp_months)
  prec_months <- rast( prec_months)
  
  #~~ 2. Seasonality TEMP and PREC values
  TSE <- app(temp_months, fun= function(x) sd(x, na.rm=T)) 
  PSE <- app(prec_months, fun= function(x) sd(x, na.rm=T)/mean(x, na.rm=T)) 
  
  #~~ 3. Annual TEMP and PREC values 
  
  MAT <- mean(temp_months)
  MAP <- sum(prec_months)
  
  #~~ 4. Koppen Zones
  
  KOPPEN <- MAT
  values(KOPPEN) <- NA
  
  # seasonality
  # seasonal mins
  october_march_min   <- min(prec_months[[which(grepl(c("P_month=9|P_month=10|P_month=11|P_month=0|P_month=1|P_month=2"), names(prec_months)))]])
  april_september_min <- min(prec_months[[which(grepl(c("P_month=3|P_month=4|P_month=5|P_month=6|P_month=7|P_month=8"), names(prec_months)))]])
  
  monthly_mins <- swapHemispheres(north_summer =april_september_min, south_summer=october_march_min )
  summer_min <- monthly_mins[[1]]
  winter_min <- monthly_mins[[2]]
  
  # seasonal maxs
  october_march_max <- max(prec_months[[which(grepl(c("P_month=9|P_month=10|P_month=11|P_month=0|P_month=1|P_month=2"), names(prec_months)))]])
  april_september_max <- max(prec_months[[which(grepl(c("P_month=3|P_month=4|P_month=5|P_month=6|P_month=7|P_month=8"), names(prec_months)))]])
  
  monthly_maxs <- swapHemispheres(north_summer =april_september_max, south_summer=october_march_max )
  summer_max <- monthly_maxs[[1]]
  winter_max <- monthly_maxs[[2]]
  
  # seasonal total
  october_march_sum <- sum(prec_months[[which(grepl(c("P_month=9|P_month=10|P_month=11|P_month=0|P_month=1|P_month=2"), names(prec_months)))]])
  april_september_sum <- sum(prec_months[[which(grepl(c("P_month=3|P_month=4|P_month=5|P_month=6|P_month=7|P_month=8"), names(prec_months)))]])
  
  monthly_sum <- swapHemispheres(north_summer =april_september_sum, south_summer=october_march_sum )
  summer_sum <- monthly_sum[[1]]
  winter_sum <- monthly_sum[[2]]
  
  #~ 4 .1 TROPICAL
  # tropical all months have temps > 18 degrees
  values(KOPPEN)[which(apply(values(temp_months), 1, FUN=function(x){min(x, na.rm=T)>=18 }))] <- 1
  
  #~ 4 .2 TEMPERATE and Cold
  # Temperates have coldest month between -3 and 18 and at least one month above 10
  # Cool temperate (Cfb/Cfc/Cwb/Cwc)
  values(KOPPEN)[which(apply(values(temp_months), 1, FUN=function(x){min(x, na.rm=T) < 18 & any(x > 10)}))] <- 2
  
  # Warm Temperate - Cfb have at least 4 months above 10 
  values(KOPPEN)[which(apply(values(temp_months), 1, FUN=function(x){min(x, na.rm=T)> 0 & min(x, na.rm=T) < 18 & length(which(x > 10)) >= 4}))] <- 2
  
  # Cfa - subtropical humid Temperate with one month above 22 
  values(KOPPEN)[which(apply(values(temp_months), 1, FUN=function(x){min(x, na.rm=T)> 0 & min(x, na.rm=T) < 18 & length(which(x > 22)) >= 1}))] <- 2
  
  # Csa and Csb and Csc 
  Mediterranean <- summer_min < 40 & summer_min < (winter_max/3)
  values(KOPPEN)[which(apply(values(temp_months), 1, FUN=function(x){min(x, na.rm=T)> 0 & min(x, na.rm=T) < 18}) &
                         apply(values(Mediterranean), 1, FUN=function(x){x==1}) )] <- 3
  
  # Cwa - Subtropical monsoon temperate, one month above 22
  values(KOPPEN)[which(values(KOPPEN) == 2 & values(winter_min < (summer_max/10))==1)] <- 1
  
  #~ 4 .3 ARID
  # determine precipitation threshold (which depends on precipitation seasonality)
  threshold <- MAT*20
  prop_summer <- summer_sum / (summer_sum + winter_sum)
  prop_winter <- winter_sum / (summer_sum + winter_sum)
  
  # set threshold
  values(threshold )[which(values(prop_summer)>= 0.7)] <- values(threshold )[which(values(prop_summer)>= 0.7)] + 280
  values(threshold )[which(values(prop_winter)>= 0.7)] <- values(threshold )[which(values(prop_winter)>= 0.7)] + 140
  semi_arid <- MAP < (threshold)
  arid <- MAP < (threshold/2)
  
  # arid is semi arid + arid and warm
  values(KOPPEN)[which(values(MAP) < (values(threshold)))] <- 4
  values(KOPPEN)[which(values(MAP) < (values(threshold)/2))] <- 5
  
  #~ 4 .4 POLAR
  values(KOPPEN)[which(apply(values(temp_months), 1, FUN=function(x){max(x, na.rm=T) <= 10}))] <- 6 #
  
  # fix longitude
  KOPPEN <- readjustLon(KOPPEN)
  
  # crop elevation
  values(KOPPEN)[which(values(scotese_ras <=0))] <- NA
  values(MAT)[which(values(scotese_ras <=0))] <- NA
  values(MAP)[which(values(scotese_ras <=0))] <- NA
  values(TSE)[which(values(scotese_ras <=0))] <- NA
  values(PSE)[which(values(scotese_ras <=0))] <- NA
  return(list(MAT=MAT, MAP=MAP, TSE=TSE, PSE=PSE, KOP=KOPPEN, DEM=scotese_ras))
}

####---- Data ----####

# get paleoclimate variables Li et al. 2022 Sci Data 9, 371 (2022): https://doi.org/10.1038/s41597-022-01490-4 
li_etal_rast <- rast(file.path("../Paleoclimate/Li_etal_2022/High_Resolution_Climate_Simulation_Dataset_540_Myr.nc"))

# get paleoelevations Scotese and Wright (2018): https://www.earthbyte.org/paleodem-resource-scotese-and-wright-2018/
scotese_and_wright_rast <- list.files("../Paleoclimate/PaleoDEMS_long_lat_elev_csv_v2/PaleoDEMS_long_lat_elev_csv_v2.csv/")[1:30]

# Present-day Koppen from Beck et al
beck_etal_rast <-rast("data/koppen_Beck/Beck_KG_V1_present_0p5.tif")

 
# koppenise all of the Li et al data
# select years we are interested in
li_years <- c(54:41)

# create a list
raster_list  <- list()

#loop over list and fill with Koppen-ized paleoclimate 
for(year in 1:length(li_years)){
  print((54-li_years[year] )*10)
  raster_list [[year]] <- koppenizerLi(year, li_years[year], li_etal_rast, scotese_and_wright_rast[seq(from=1, to=29, by=2)])
}

# take just the koppen data for now
koppen_raster <- raster_list[[1]]$KOP

# combine the Li koppen data with the Beck et al data for the present day
beck_etal_df <- as.data.frame(beck_etal_rast, xy=TRUE)

# simplify biomes according to BGB
beck_etal_df$biome <- NA
beck_etal_df$biome[which(beck_etal_df$Beck_KG_V1_present_0p5 %in% c(0))] <- NA
beck_etal_df$biome[which(beck_etal_df$Beck_KG_V1_present_0p5 %in% c(1,2,3, 11))] <- "TropicalRainforest"
beck_etal_df$biome[which(beck_etal_df$Beck_KG_V1_present_0p5 %in% c(7, 6))] <- "SemiArid"
beck_etal_df$biome[which(beck_etal_df$Beck_KG_V1_present_0p5 %in% c(5, 4))] <- "Arid"
beck_etal_df$biome[which(beck_etal_df$Beck_KG_V1_present_0p5 %in% c(15, 16, 14, 12, 19, 26, 25))] <- "Subtropical"
beck_etal_df$biome[which(beck_etal_df$Beck_KG_V1_present_0p5 %in% c(8, 9, 10))] <- "Mediterranean"
beck_etal_df$biome <- factor(beck_etal_df$biome, levels=c("TropicalRainforest",  "Subtropical", "Mediterranean", "SemiArid", "Arid", "NA"))
beck_etal_df$koppen <- as.numeric(beck_etal_df$biome)

# turn into a raster
beck_etal_rast  <- rast(beck_etal_df[,c("x", "y", "koppen")], type="xyz")
beck_etal_rast   <- resample(beck_etal_rast, koppen_raster, method="near" )

# seperate ou different paleoclimate variables (for climate space analysis)
MAT <- rast(lapply(raster_list, FUN=function(x)x[["MAT"]]))
MAP <- rast(lapply(raster_list, FUN=function(x)x[["MAP"]]))
TSE <- rast(lapply(raster_list, FUN=function(x)x[["TSE"]]))
PSE <- rast(lapply(raster_list, FUN=function(x)x[["PSE"]]))
DEM <- rast(lapply(raster_list, FUN=function(x)x[["DEM"]]))
KOP <- rast(lapply(raster_list, FUN=function(x)x[["KOP"]]))
time(KOP) <- c(130,120,110, 100,90, 80,70, 60,50, 40, 30,20, 10, 0)

# make the present-day Beck et al.
KOP[[1]] <- beck_etal_rast

# save as netCDF
writeCDF(KOP, "Dataset_S6.nc",zname=, overwrite=TRUE, varname="kop", 
         longname="modified Koppen-Geiger climate zones", unit="m")
