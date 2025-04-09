####---- Libraries ----####
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
koppenizerLi <- function(year,li_year, clim_rast, scotese){
  
  #~~ 1. Elevation Values
  scotese_ras <- rast(read.csv(file.path("External_datasets/PaleoDEMS_long_lat_elev_csv_v2", "PaleoDEMS_long_lat_elev_csv_v2.csv",scotese[year] )), type="xyz")
  crs(scotese_ras) <- "+init=EPSG:4326"
  
  #~~ 1. Monthly TEMP and PREC values
  # Monthly
  clim_rast_i <- clim_rast[[which(grepl(paste0("simulation=", li_year),names(clim_rast)))]]
  crs(clim_rast_i) <- "+init=EPSG:4326"
  clim_rast_i <- project(clim_rast_i, scotese_ras) 
  monthlies <- names(clim_rast_i)
  temp_months <- vector("list", 12)
  prec_months <- vector("list", 12)
  for(month in 0:11){
    temp_months[[month+1]] <- clim_rast_i[[grepl(paste0("T_month=", month), monthlies)]]
    prec_months[[month+1]] <- clim_rast_i[[grepl(paste0("P_month=", month), monthlies)]]
    ext(temp_months[[month+1]]) <- ext(clim_rast_i)
    ext(prec_months[[month+1]]) <- ext(clim_rast_i)
  }
  
  temp_months <- rast( temp_months)
  prec_months <- rast( prec_months)
  
  #~~ 2. Seasonality TEMP and PREC values
  TSE <- app(rast(temp_months), fun= function(x) sd(x, na.rm=T)) 
  PSE <- app(rast(prec_months), fun= function(x) sd(x, na.rm=T)/mean(x, na.rm=T)) 
  
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
  
  plot(KOPPEN)
  return(KOPPEN)
}

####---- DATA ----####

# get paleoclimate variables Li et al. 2022 Sci Data 9, 371 (2022): https://doi.org/10.1038/s41597-022-01490-4 
clim_rast <- rast(file.path("High_Resolution_Climate_Simulation_Dataset_540_Myr.nc"))

# get paleoelevations Scotese and Wright (2018): https://www.earthbyte.org/paleodem-resource-scotese-and-wright-2018/
scotese <- list.files("PaleoDEMS_long_lat_elev_csv_v2/PaleoDEMS_long_lat_elev_csv_v2.csv/")[1:30]

# index the years we would like to extract from Li et al
li_years_index <- c(54:41)

####---- Paleobiome reconstruction ----####

# run  or just load?
run_slow <- FALSE

# set up a list to put in the rasters
koppen_li  <- list()

# loop over each year and create raster of modified Koppen-Geiger biomes
print(paste("extracting modified Koppen-Geiger biome classification from paleotemperature and paleoclimate values"))

if(run_slow==TRUE){
  
  for(year in 1:length(li_years_index)){
    print((54-li_years_index[year] )*10)
    koppen_li[[year]] <- koppenizerLi(year, li_years_index[year], clim_rast, scotese[seq(from=1, to=29, by=2)])
  }
  
  writeRaster(rast(koppen_li), "Li_Koppenised.tif", overwrite=T)
  
} else {
  
  koppen_li <- rast("Li_Koppenised.tif")
}


# get years
li_years <- seq(from=0, to=130, by=10)

print(paste("seperating biomes from different continents through time"))

# set up data structure
RegionalPolygons <- vector("list", 7)
RegionalPolygons <- lapply(RegionalPolygons, FUN=function(x){vector("list", 12)})
RegionalPolygons <- lapply(RegionalPolygons, FUN=function(x){names(x) <- c(
  "AustraloPapuan_Tropical",
  "AustraloPapuan_Subtropical",
  "AustraloPapuan_Mediterranean",
  "AustraloPapuan_SemiArid",
  "AustraloPapuan_Arid",
  "SthAmerica",
  "Madagascar",
  "Cape",
  "TropicalAsia",
  "NewCaledonia",
  "NewZealand",
  "Antarctica"); return(x)})

### ~~~~ Lower Cretaceous 140-120 ~~~###

# load roughly drawn continent outlines from QGIS
TropicalAsia  <- vect("External_datasets/QGIS_regional_polygons/120Ma_TropicalAsia.shp")
Afrotropics   <- vect("External_datasets/QGIS_regional_polygons/120Ma_Afrotropics.shp")
Neotropics    <- vect("External_datasets/QGIS_regional_polygons/120Ma_SouthAmerica.shp")
AustraloPapua <- vect("External_datasets/QGIS_regional_polygons/120Ma_AustraloPapuan.shp")
Antarctica    <- vect("External_datasets/QGIS_regional_polygons/120Ma_Ant.shp")
Madagascar    <- vect("External_datasets/QGIS_regional_polygons/120Ma_Madagascar.shp")
NewZealand    <- vect("External_datasets/QGIS_regional_polygons/120Ma_NZ.shp")
NewCaledonia  <- vect("External_datasets/QGIS_regional_polygons/120Ma_NC.shp")

# subset suitable habitat in each
li_years
ras_timestep <- koppen_li[[which(li_years==120)]]
names(ras_timestep) <- "koppen"

# Tropical Asia
TropicalAsia_ras <- mask(ras_timestep, TropicalAsia)
values(TropicalAsia_ras)[which(values(TropicalAsia_ras$koppen %in% c(3,4,5,6)))] <- NA
values(TropicalAsia_ras)[which(!is.na(values(TropicalAsia_ras$koppen)))] <- 1
TropicalAsia_poly <- as.polygons(TropicalAsia_ras, dissolve=T)

# Madagascar
Madagascar_ras <- mask(ras_timestep,Madagascar)
values(Madagascar_ras)[which(values(Madagascar_ras$koppen %in% c(6)))] <- NA
values(Madagascar_ras)[which(!is.na(values(Madagascar_ras$koppen)))] <- 1
Madagascar_poly <- as.polygons(Madagascar_ras, dissolve=T)

# Tropical South America
SthAmerica_ras <- mask(ras_timestep, Neotropics)
values(SthAmerica_ras)[which(values(SthAmerica_ras$koppen %in% c(6)))] <- NA
values(SthAmerica_ras)[which(!is.na(values(SthAmerica_ras$koppen)))] <- 1
SthAmerica_poly <- as.polygons(SthAmerica_ras, dissolve=T)

# AfroTropics
Afrotropics_ras <- mask(ras_timestep,Afrotropics)
values(Afrotropics_ras)[which(values(Afrotropics_ras$koppen %in% c(6)))] <- NA
values(Afrotropics_ras)[which(!is.na(values(Afrotropics_ras$koppen)))] <- 1
Afrotropics_poly <- as.polygons(Afrotropics_ras, dissolve=T)

# New Zealand
NewZealand_ras <- mask(ras_timestep,NewZealand)
values(NewZealand_ras)[which(values(NewZealand_ras$koppen %in% c(6)))] <- NA
values(NewZealand_ras)[which(!is.na(values(NewZealand_ras$koppen)))] <- 1
NewZealand_poly <- as.polygons(NewZealand_ras, dissolve=T)

# New Caledonia
NewCaledonia_ras <- mask(ras_timestep,NewCaledonia)
values(NewCaledonia_ras)[which(values(NewCaledonia_ras$koppen %in% c(6)))] <- NA
values(NewCaledonia_ras)[which(!is.na(values(NewCaledonia_ras$koppen)))] <- 1
NewCaledonia_poly <- as.polygons(NewCaledonia_ras, dissolve=T)

# Ausralo-Papua
AustraloPapuan_ras <- mask(ras_timestep, AustraloPapua)
values(AustraloPapuan_ras)[which(values(AustraloPapuan_ras$koppen) %in% c(6))] <- NA

AustraloPapuan_ras_subtrop <- AustraloPapuan_ras
AustraloPapuan_ras_temp    <- AustraloPapuan_ras

values(AustraloPapuan_ras_subtrop )[which(!values(AustraloPapuan_ras_subtrop )==2)] <- NA

AustraloPapuan_poly <- as.polygons(AustraloPapuan_ras,    dissolve=T)
AustraloPapuan_poly <- erase(AustraloPapuan_poly   , NewCaledonia_poly)

AustraloPapuan_poly_subtrop <- as.polygons(AustraloPapuan_ras_subtrop, dissolve=T)
AustraloPapuan_poly_subtrop <- erase(AustraloPapuan_poly_subtrop , NewCaledonia_poly)

# Antarctica
Antarctica_ras <- mask(ras_timestep, Antarctica)
values(Antarctica_ras)[which(!is.na(values(Antarctica_ras$koppen)))] <- 1
Antarctica_poly <- as.polygons(Antarctica_ras, dissolve=T)
Antarctica_poly <- erase(Antarctica_poly, AustraloPapuan_poly)
Antarctica_poly <- erase(Antarctica_poly, Afrotropics_poly)
Antarctica_poly <- erase(Antarctica_poly, SthAmerica_poly)
Antarctica_poly <- erase(Antarctica_poly, Madagascar_poly)
Antarctica_poly <- erase(Antarctica_poly, NewZealand_poly )
Antarctica_poly <- erase(Antarctica_poly, NewCaledonia_poly)

plot_it <- T
if(plot_it == TRUE){
  par(mfrow=c(1,1))
  plot(ras_timestep)
  plot(TropicalAsia_poly, add=T)
  plot(Afrotropics_poly, add=T)
  plot(AustraloPapuan_poly, add=T)
  plot(Antarctica_poly, add=T)
  plot(SthAmerica_poly, add=T)
  plot(Madagascar_poly, add=T)
  plot(NewZealand_poly, add=T)
  plot(NewCaledonia_poly, add=T)
}

RegionalPolygons[[1]][["TropicalAsia"]]                <- TropicalAsia_poly
RegionalPolygons[[1]][["Cape"]]                        <- Afrotropics_poly
RegionalPolygons[[1]][["Madagascar"]]                  <- Madagascar_poly
RegionalPolygons[[1]][["SthAmerica"]]                  <- SthAmerica_poly
RegionalPolygons[[1]][["Antarctica"]]                  <- Antarctica_poly
RegionalPolygons[[1]][["NewZealand"]]                  <- NewZealand_poly
RegionalPolygons[[1]][["NewCaledonia"]]                <- NewCaledonia_poly
RegionalPolygons[[1]][["AustraloPapuan_Tropical"]]     <- NA
RegionalPolygons[[1]][["AustraloPapuan_Subtropical"]]  <- AustraloPapuan_poly_subtrop
RegionalPolygons[[1]][["AustraloPapuan_SemiArid"]]     <- NA
RegionalPolygons[[1]][["AustraloPapuan_Mediterranean"]]<- NA
RegionalPolygons[[1]][["AustraloPapuan_Arid"]]         <- NA

### ~~~~ Lower Cretaceous 120-100 ~~~###

# load roughly drawn continent outlines from QGIS
TropicalAsia <-   vect("External_datasets/QGIS_regional_polygons/100Ma_TropicalAsia.shp")
Afrotropics <-    vect("External_datasets/QGIS_regional_polygons/100Ma_Afrotropics.shp")
Neotropics  <-    vect("External_datasets/QGIS_regional_polygons/100Ma_SouthAmerica.shp")
AustraloPapua <-  vect("External_datasets/QGIS_regional_polygons/100Ma_AustraloPapua.shp")
Antarctica <-     vect("External_datasets/QGIS_regional_polygons/100Ma_Ant.shp")
Madagascar<-      vect("External_datasets/QGIS_regional_polygons/100Ma_Madagascar.shp")
NewZealand <-     vect(ext(c(146.644518272425, 158.604651162791, -72.4585188861313, -63.8212955479859)))
NewCaledonia <-   vect(ext(c(144.651162790698, 149.302325581395, -51.1976614383887, -48.5400542574209)))


# subset suitable habitat in each
li_years
ras_timestep <- koppen_li[[which(li_years==100)]]
names(ras_timestep) <- "koppen"

# Tropical Asia
TropicalAsia_ras <- mask(ras_timestep, TropicalAsia)
values(TropicalAsia_ras)[which(values(TropicalAsia_ras$koppen %in% c(3,4,5,6)))] <- NA
values(TropicalAsia_ras)[which(!is.na(values(TropicalAsia_ras$koppen)))] <- 1
TropicalAsia_poly <- as.polygons(TropicalAsia_ras, dissolve=T)

# Madagascar
Madagascar_ras <- mask(ras_timestep,Madagascar)
values(Madagascar_ras)[which(values(Madagascar_ras$koppen %in% c(6)))] <- NA
values(Madagascar_ras)[which(!is.na(values(Madagascar_ras$koppen)))] <- 1
Madagascar_poly <- as.polygons(Madagascar_ras, dissolve=T)


# Tropical South America
SthAmerica_ras <- mask(ras_timestep, Neotropics)
values(SthAmerica_ras)[which(values(SthAmerica_ras$koppen %in% c(6)))] <- NA
values(SthAmerica_ras)[which(!is.na(values(SthAmerica_ras$koppen)))] <- 1
SthAmerica_poly <- as.polygons(SthAmerica_ras, dissolve=T)

# AfroTropics
Afrotropics_ras <- mask(ras_timestep,Afrotropics)
values(Afrotropics_ras)[which(values(Afrotropics_ras$koppen %in% c(6)))] <- NA
values(Afrotropics_ras)[which(!is.na(values(Afrotropics_ras$koppen)))] <- 1
Afrotropics_poly <- as.polygons(Afrotropics_ras, dissolve=T)

# New Zealand
NewZealand_ras <- mask(ras_timestep,NewZealand)
values(NewZealand_ras)[which(values(NewZealand_ras$koppen %in% c(6)))] <- NA
values(NewZealand_ras)[which(!is.na(values(NewZealand_ras$koppen)))] <- 1
NewZealand_poly <- as.polygons(NewZealand_ras, dissolve=T)

# New Caledonia
NewCaledonia_ras <- crop(ras_timestep,NewCaledonia)
values(NewCaledonia_ras) <- 1
NewCaledonia_poly <- as.polygons(NewCaledonia_ras, dissolve=T)

# Ausralo-Papua
AustraloPapuan_ras <- mask(ras_timestep, AustraloPapua)
values(AustraloPapuan_ras)[which(values(AustraloPapuan_ras$koppen) %in% c(6))] <- NA

AustraloPapuan_ras_trop    <- AustraloPapuan_ras
AustraloPapuan_ras_subtrop <- AustraloPapuan_ras
AustraloPapuan_ras_med     <- AustraloPapuan_ras

values(AustraloPapuan_ras_trop    )[which(!values(AustraloPapuan_ras_trop    )==1)] <- NA
values(AustraloPapuan_ras_subtrop )[which(!values(AustraloPapuan_ras_subtrop )==2)] <- NA
values(AustraloPapuan_ras_med     )[which(!values(AustraloPapuan_ras_med     )==3)] <- NA

AustraloPapuan_poly <- as.polygons(AustraloPapuan_ras,    dissolve=T)
AustraloPapuan_poly <- erase(AustraloPapuan_poly   , NewCaledonia_poly)
AustraloPapuan_poly <- erase(AustraloPapuan_poly   , NewZealand_poly)

AustraloPapuan_poly_trop    <- as.polygons(AustraloPapuan_ras_trop,    dissolve=T)
AustraloPapuan_poly_subtrop <- as.polygons(AustraloPapuan_ras_subtrop, dissolve=T)
AustraloPapuan_poly_med     <- as.polygons(AustraloPapuan_ras_med ,    dissolve=T)

AustraloPapuan_poly_trop    <- erase(AustraloPapuan_poly_trop    , NewCaledonia_poly)
AustraloPapuan_poly_subtrop <- erase(AustraloPapuan_poly_subtrop , NewCaledonia_poly)
AustraloPapuan_poly_med    <- erase(AustraloPapuan_poly_med   , NewCaledonia_poly)

AustraloPapuan_poly_trop    <- erase(AustraloPapuan_poly_trop    , NewZealand_poly)
AustraloPapuan_poly_subtrop <- erase(AustraloPapuan_poly_subtrop , NewZealand_poly)
AustraloPapuan_poly_med    <- erase(AustraloPapuan_poly_med    , NewZealand_poly)

# Antarctica
Antarctica_ras <- mask(ras_timestep, Antarctica)
values(Antarctica_ras)[which(!is.na(values(Antarctica_ras$koppen)))] <- 1
Antarctica_poly <- as.polygons(Antarctica_ras, dissolve=T)
Antarctica_poly <- erase(Antarctica_poly, AustraloPapuan_poly)
Antarctica_poly <- erase(Antarctica_poly, Afrotropics_poly)
Antarctica_poly <- erase(Antarctica_poly, SthAmerica_poly)
Antarctica_poly <- erase(Antarctica_poly, Madagascar_poly)
Antarctica_poly <- erase(Antarctica_poly, NewZealand_poly )
Antarctica_poly <- erase(Antarctica_poly, NewCaledonia_poly)

plot_it <- T
if(plot_it == TRUE){
  par(mfrow=c(1,1))
  plot(ras_timestep)
  plot(TropicalAsia_poly, add=T)
  plot(Afrotropics_poly, add=T)
  plot(AustraloPapuan_poly, add=T)
  plot(Antarctica_poly, add=T)
  plot(SthAmerica_poly, add=T)
  plot(Madagascar_poly, add=T)
  plot(NewZealand_poly, add=T)
  plot(NewCaledonia_poly, add=T)
}


RegionalPolygons[[2]][["TropicalAsia"]]                <- TropicalAsia_poly
RegionalPolygons[[2]][["Cape"]]                        <- Afrotropics_poly
RegionalPolygons[[2]][["Madagascar"]]                  <- Madagascar_poly
RegionalPolygons[[2]][["SthAmerica"]]                  <- SthAmerica_poly
RegionalPolygons[[2]][["Antarctica"]]                  <- Antarctica_poly
RegionalPolygons[[2]][["NewZealand"]]                  <- NewZealand_poly
RegionalPolygons[[2]][["NewCaledonia"]]                <- NewCaledonia_poly
RegionalPolygons[[2]][["AustraloPapuan_Tropical"]]     <- AustraloPapuan_poly_trop
RegionalPolygons[[2]][["AustraloPapuan_Subtropical"]]  <- AustraloPapuan_poly_subtrop
RegionalPolygons[[2]][["AustraloPapuan_SemiArid"]]     <- NA
RegionalPolygons[[2]][["AustraloPapuan_Mediterranean"]]<- AustraloPapuan_poly_med
RegionalPolygons[[2]][["AustraloPapuan_Arid"]]         <- NA

### ~~~~ Lower Cretaceous 100-80 ~~~###

# load roughly drawn continent outlines from QGIS
TropicalAsia  <- vect("External_datasets/QGIS_regional_polygons/80Ma_TropicalAsia.shp")
Afrotropics   <- vect("External_datasets/QGIS_regional_polygons/80Ma_Afrotropics.shp")
Neotropics    <- vect("External_datasets/QGIS_regional_polygons/80Ma_SouthAmerica.shp")
AustraloPapua <- vect("External_datasets/QGIS_regional_polygons/80Ma_AustraloPapua.shp")
Antarctica    <- vect("External_datasets/QGIS_regional_polygons/80Ma_Ant.shp")
Madagascar    <- vect("External_datasets/QGIS_regional_polygons/80Ma_Madagascar.shp")
NewZealand    <- vect("External_datasets/QGIS_regional_polygons/80Ma_NZ.shp")
NewCaledonia  <- vect("External_datasets/QGIS_regional_polygons/80Ma_NC.shp")

# subset suitable habitat in each
li_years
ras_timestep <- koppen_li[[which(li_years==80)]]
names(ras_timestep) <- "koppen"

# Tropical Asia
TropicalAsia_ras <- mask(ras_timestep, TropicalAsia)
values(TropicalAsia_ras)[which(values(TropicalAsia_ras$koppen %in% c(3,4,5,6)))] <- NA
values(TropicalAsia_ras)[which(!is.na(values(TropicalAsia_ras$koppen)))] <- 1
TropicalAsia_poly <- as.polygons(TropicalAsia_ras, dissolve=T)

# Madagascar
Madagascar_ras <- mask(ras_timestep,Madagascar)
values(Madagascar_ras)[which(values(Madagascar_ras$koppen %in% c(6)))] <- NA
values(Madagascar_ras)[which(!is.na(values(Madagascar_ras$koppen)))] <- 1
Madagascar_poly <- as.polygons(Madagascar_ras, dissolve=T)


# Tropical South America
SthAmerica_ras <- mask(ras_timestep, Neotropics)
values(SthAmerica_ras)[which(values(SthAmerica_ras$koppen %in% c(6)))] <- NA
values(SthAmerica_ras)[which(!is.na(values(SthAmerica_ras$koppen)))] <- 1
SthAmerica_poly <- as.polygons(SthAmerica_ras, dissolve=T)

# AfroTropics
Afrotropics_ras <- mask(ras_timestep,Afrotropics)
values(Afrotropics_ras)[which(values(Afrotropics_ras$koppen %in% c(6)))] <- NA
values(Afrotropics_ras)[which(!is.na(values(Afrotropics_ras$koppen)))] <- 1
Afrotropics_poly <- as.polygons(Afrotropics_ras, dissolve=T)

# New Zealand
NewZealand_ras <- mask(ras_timestep,NewZealand)
values(NewZealand_ras)[which(values(NewZealand_ras$koppen %in% c(6)))] <- NA
values(NewZealand_ras)[which(!is.na(values(NewZealand_ras$koppen)))] <- 1
NewZealand_poly <- as.polygons(NewZealand_ras, dissolve=T)

# New Caledonia
NewCaledonia_ras <- crop(ras_timestep,NewCaledonia)
values(NewCaledonia_ras)[which(!is.na(values(NewCaledonia_ras$koppen)))] <- 1
NewCaledonia_poly <- as.polygons(NewCaledonia_ras, dissolve=T)

# Ausralo-Papua
AustraloPapuan_ras <- mask(ras_timestep, AustraloPapua)
values(AustraloPapuan_ras)[which(values(AustraloPapuan_ras$koppen) %in% c(6))] <- NA

AustraloPapuan_ras_trop    <- AustraloPapuan_ras
AustraloPapuan_ras_subtrop <- AustraloPapuan_ras
AustraloPapuan_ras_med    <- AustraloPapuan_ras

values(AustraloPapuan_ras_trop    )[which(!values(AustraloPapuan_ras_trop    )==1)] <- NA
values(AustraloPapuan_ras_subtrop )[which(!values(AustraloPapuan_ras_subtrop )==2)] <- NA
values(AustraloPapuan_ras_med     )[which(!values(AustraloPapuan_ras_med     )==3)] <- NA

AustraloPapuan_poly <- as.polygons(AustraloPapuan_ras,    dissolve=T)
AustraloPapuan_poly <- erase(AustraloPapuan_poly   , NewCaledonia_poly)
AustraloPapuan_poly <- erase(AustraloPapuan_poly   , NewZealand_poly)

AustraloPapuan_poly_trop    <- as.polygons(AustraloPapuan_ras_trop,    dissolve=T)
AustraloPapuan_poly_subtrop <- as.polygons(AustraloPapuan_ras_subtrop, dissolve=T)
AustraloPapuan_poly_med     <- as.polygons(AustraloPapuan_ras_med ,    dissolve=T)

AustraloPapuan_poly_trop    <- erase(AustraloPapuan_poly_trop    , NewCaledonia_poly)
AustraloPapuan_poly_subtrop <- erase(AustraloPapuan_poly_subtrop , NewCaledonia_poly)
AustraloPapuan_poly_med    <- erase(AustraloPapuan_poly_med   , NewCaledonia_poly)

AustraloPapuan_poly_trop    <- erase(AustraloPapuan_poly_trop    , NewZealand_poly)
AustraloPapuan_poly_subtrop <- erase(AustraloPapuan_poly_subtrop , NewZealand_poly)
AustraloPapuan_poly_med    <- erase(AustraloPapuan_poly_med    , NewZealand_poly)

# Antarctica
Antarctica_ras <- mask(ras_timestep, Antarctica)
values(Antarctica_ras)[which(!is.na(values(Antarctica_ras$koppen)))] <- 1
Antarctica_poly <- as.polygons(Antarctica_ras, dissolve=T)
Antarctica_poly <- erase(Antarctica_poly, AustraloPapuan_poly)
Antarctica_poly <- erase(Antarctica_poly, Afrotropics_poly)
Antarctica_poly <- erase(Antarctica_poly, SthAmerica_poly)
Antarctica_poly <- erase(Antarctica_poly, Madagascar_poly)
Antarctica_poly <- erase(Antarctica_poly, NewZealand_poly )
Antarctica_poly <- erase(Antarctica_poly, NewCaledonia_poly)

plot_it <- T
if(plot_it == TRUE){
  par(mfrow=c(1,1))
  plot(ras_timestep)
  plot(TropicalAsia_poly, add=T)
  plot(Afrotropics_poly, add=T)
  plot(AustraloPapuan_poly, add=T)
  plot(Antarctica_poly, add=T)
  plot(SthAmerica_poly, add=T)
  plot(Madagascar_poly, add=T)
  plot(NewZealand_poly, add=T)
  plot(NewCaledonia_poly, add=T)
}


RegionalPolygons[[3]][["TropicalAsia"]]                <- TropicalAsia_poly
RegionalPolygons[[3]][["Cape"]]                        <- Afrotropics_poly
RegionalPolygons[[3]][["Madagascar"]]                  <- Madagascar_poly
RegionalPolygons[[3]][["SthAmerica"]]                  <- SthAmerica_poly
RegionalPolygons[[3]][["Antarctica"]]                  <- Antarctica_poly
RegionalPolygons[[3]][["NewZealand"]]                  <- NewZealand_poly
RegionalPolygons[[3]][["NewCaledonia"]]                <- NewCaledonia_poly
RegionalPolygons[[3]][["AustraloPapuan_Tropical"]]     <- AustraloPapuan_poly_trop
RegionalPolygons[[3]][["AustraloPapuan_Subtropical"]]  <- AustraloPapuan_poly_subtrop
RegionalPolygons[[3]][["AustraloPapuan_SemiArid"]]     <- NA
RegionalPolygons[[3]][["AustraloPapuan_Mediterranean"]]<- AustraloPapuan_poly_med
RegionalPolygons[[3]][["AustraloPapuan_Arid"]]         <- NA

### ~~~~ Lower Cretaceous - Palecoene 80-60 ~~~###
# load continents from QGIS
TropicalAsia  <- vect("External_datasets/QGIS_regional_polygons/60Ma_TropicalAsia.shp")
Afrotropics   <- vect("External_datasets/QGIS_regional_polygons/60Ma_Afrotropics.shp")
Neotropics    <- vect("External_datasets/QGIS_regional_polygons/60Ma_SouthAmerica.shp")
AustraloPapua <- vect("External_datasets/QGIS_regional_polygons/60Ma_AustraloPapua.shp")
Antarctica    <- vect("External_datasets/QGIS_regional_polygons/60Ma_Ant.shp")
Madagascar    <- vect("External_datasets/QGIS_regional_polygons/60Ma_Madagascar.shp")
NewZealand    <- vect("External_datasets/QGIS_regional_polygons/60Ma_NZ.shp")
NewCaledonia  <- vect("External_datasets/QGIS_regional_polygons/60Ma_NC.shp")

# subset suitable habitat in each
li_years
ras_timestep <- koppen_li[[which(li_years==60)]]
names(ras_timestep) <- "koppen"

# Tropical Asia
TropicalAsia_ras <- mask(ras_timestep, TropicalAsia)
values(TropicalAsia_ras)[which(values(TropicalAsia_ras$koppen %in% c(3,4,5,6)))] <- NA
values(TropicalAsia_ras)[which(!is.na(values(TropicalAsia_ras$koppen)))] <- 1
TropicalAsia_poly <- as.polygons(TropicalAsia_ras, dissolve=T)

# Madagascar
Madagascar_ras <- mask(ras_timestep,Madagascar)
values(Madagascar_ras)[which(values(Madagascar_ras$koppen %in% c(6)))] <- NA
values(Madagascar_ras)[which(!is.na(values(Madagascar_ras$koppen)))] <- 1
Madagascar_poly <- as.polygons(Madagascar_ras, dissolve=T)


# Tropical South America
SthAmerica_ras <- mask(ras_timestep, Neotropics)
values(SthAmerica_ras)[which(values(SthAmerica_ras$koppen %in% c(6)))] <- NA
values(SthAmerica_ras)[which(!is.na(values(SthAmerica_ras$koppen)))] <- 1
SthAmerica_poly <- as.polygons(SthAmerica_ras, dissolve=T)

# AfroTropics
Afrotropics_ras <- mask(ras_timestep,Afrotropics)
values(Afrotropics_ras)[which(values(Afrotropics_ras$koppen %in% c(6)))] <- NA
values(Afrotropics_ras)[which(!is.na(values(Afrotropics_ras$koppen)))] <- 1
Afrotropics_poly <- as.polygons(Afrotropics_ras, dissolve=T)

# New Zealand
NewZealand_ras <- mask(ras_timestep,NewZealand)
values(NewZealand_ras)[which(values(NewZealand_ras$koppen %in% c(6)))] <- NA
values(NewZealand_ras)[which(!is.na(values(NewZealand_ras$koppen)))] <- 1
NewZealand_poly <- as.polygons(NewZealand_ras, dissolve=T)

# New Caledonia
NewCaledonia_ras <- crop(ras_timestep,NewCaledonia)
values(NewCaledonia_ras)[which(!is.na(values(NewCaledonia_ras$koppen)))] <- 1
NewCaledonia_poly <- as.polygons(NewCaledonia_ras, dissolve=T)


# Ausralo-Papua
AustraloPapuan_ras <- mask(ras_timestep, AustraloPapua)
values(AustraloPapuan_ras)[which(values(AustraloPapuan_ras$koppen) %in% c(6))] <- NA

AustraloPapuan_ras_trop    <- AustraloPapuan_ras
AustraloPapuan_ras_subtrop <- AustraloPapuan_ras
AustraloPapuan_ras_med    <- AustraloPapuan_ras

values(AustraloPapuan_ras_trop    )[which(!values(AustraloPapuan_ras_trop    )==1)] <- NA
values(AustraloPapuan_ras_subtrop )[which(!values(AustraloPapuan_ras_subtrop )==2)] <- NA
values(AustraloPapuan_ras_med     )[which(!values(AustraloPapuan_ras_med     )==3)] <- NA

AustraloPapuan_poly <- as.polygons(AustraloPapuan_ras,    dissolve=T)
AustraloPapuan_poly <- erase(AustraloPapuan_poly   , NewCaledonia_poly)
AustraloPapuan_poly <- erase(AustraloPapuan_poly   , NewZealand_poly)

AustraloPapuan_poly_trop    <- as.polygons(AustraloPapuan_ras_trop,    dissolve=T)
AustraloPapuan_poly_subtrop <- as.polygons(AustraloPapuan_ras_subtrop, dissolve=T)
AustraloPapuan_poly_med     <- as.polygons(AustraloPapuan_ras_med ,    dissolve=T)

AustraloPapuan_poly_trop    <- erase(AustraloPapuan_poly_trop    , NewCaledonia_poly)
AustraloPapuan_poly_subtrop <- erase(AustraloPapuan_poly_subtrop , NewCaledonia_poly)
AustraloPapuan_poly_med    <- erase(AustraloPapuan_poly_med   , NewCaledonia_poly)

AustraloPapuan_poly_trop    <- erase(AustraloPapuan_poly_trop    , NewZealand_poly)
AustraloPapuan_poly_subtrop <- erase(AustraloPapuan_poly_subtrop , NewZealand_poly)
AustraloPapuan_poly_med    <- erase(AustraloPapuan_poly_med    , NewZealand_poly)

# Antarctica
Antarctica_ras <- mask(ras_timestep, Antarctica)
values(Antarctica_ras)[which(!is.na(values(Antarctica_ras$koppen)))] <- 1
Antarctica_poly <- as.polygons(Antarctica_ras, dissolve=T)
Antarctica_poly <- erase(Antarctica_poly, AustraloPapuan_poly)
Antarctica_poly <- erase(Antarctica_poly, Afrotropics_poly)
Antarctica_poly <- erase(Antarctica_poly, SthAmerica_poly)
Antarctica_poly <- erase(Antarctica_poly, Madagascar_poly)
Antarctica_poly <- erase(Antarctica_poly, NewZealand_poly )
Antarctica_poly <- erase(Antarctica_poly, NewCaledonia_poly)

plot_it <- T
if(plot_it == TRUE){
  par(mfrow=c(1,1))
  plot(ras_timestep)
  plot(TropicalAsia_poly, add=T)
  plot(Afrotropics_poly, add=T)
  plot(AustraloPapuan_poly, add=T)
  plot(Antarctica_poly, add=T)
  plot(SthAmerica_poly, add=T)
  plot(Madagascar_poly, add=T)
  plot(NewZealand_poly, add=T)
  plot(NewCaledonia_poly, add=T)
}


RegionalPolygons[[4]][["TropicalAsia"]]                <- TropicalAsia_poly
RegionalPolygons[[4]][["Cape"]]                        <- Afrotropics_poly
RegionalPolygons[[4]][["Madagascar"]]                  <- Madagascar_poly
RegionalPolygons[[4]][["SthAmerica"]]                  <- SthAmerica_poly
RegionalPolygons[[4]][["Antarctica"]]                  <- Antarctica_poly
RegionalPolygons[[4]][["NewZealand"]]                  <- NewZealand_poly
RegionalPolygons[[4]][["NewCaledonia"]]                <- NewCaledonia_poly
RegionalPolygons[[4]][["AustraloPapuan_Tropical"]]     <- AustraloPapuan_poly_trop
RegionalPolygons[[4]][["AustraloPapuan_Subtropical"]]  <- AustraloPapuan_poly_subtrop
RegionalPolygons[[4]][["AustraloPapuan_SemiArid"]]     <- NA
RegionalPolygons[[4]][["AustraloPapuan_Mediterranean"]]<- AustraloPapuan_poly_med
RegionalPolygons[[4]][["AustraloPapuan_Arid"]]         <- NA

### ~~~~ Palecoene - Eocene 60-40 ~~~###
# load roughly drawn continent outlines from QGIS

TropicalAsia  <- vect("External_datasets/QGIS_regional_polygons/40Ma_TropicalAsia.shp")
Afrotropics   <- vect("External_datasets/QGIS_regional_polygons/40Ma_Afrotropics.shp")
Neotropics    <- vect("External_datasets/QGIS_regional_polygons/40Ma_SouthAmerica.shp")
AustraloPapua <- vect("External_datasets/QGIS_regional_polygons/40Ma_AustraloPapua.shp")
Antarctica    <- vect("External_datasets/QGIS_regional_polygons/40Ma_Ant.shp")
Madagascar    <- vect("External_datasets/QGIS_regional_polygons/40Ma_Madagascar.shp")
NewZealand    <- vect("External_datasets/QGIS_regional_polygons/40Ma_NZ.shp")
NewCaledonia  <- vect("External_datasets/QGIS_regional_polygons/40Ma_NC.shp")

# subset suitable habitat in each
li_years
ras_timestep <- koppen_li[[which(li_years==40)]]
names(ras_timestep) <- "koppen"

# Tropical Asia
TropicalAsia_ras <- mask(ras_timestep, TropicalAsia)
values(TropicalAsia_ras)[which(values(TropicalAsia_ras$koppen %in% c(3,4,5,6)))] <- NA
values(TropicalAsia_ras)[which(!is.na(values(TropicalAsia_ras$koppen)))] <- 1
TropicalAsia_poly <- as.polygons(TropicalAsia_ras, dissolve=T)

# Madagascar
Madagascar_ras <- mask(ras_timestep,Madagascar)
values(Madagascar_ras)[which(values(Madagascar_ras$koppen %in% c(6)))] <- NA
values(Madagascar_ras)[which(!is.na(values(Madagascar_ras$koppen)))] <- 1
Madagascar_poly <- as.polygons(Madagascar_ras, dissolve=T)


# Tropical South America
SthAmerica_ras <- mask(ras_timestep, Neotropics)
values(SthAmerica_ras)[which(values(SthAmerica_ras$koppen %in% c(6)))] <- NA
values(SthAmerica_ras)[which(!is.na(values(SthAmerica_ras$koppen)))] <- 1
SthAmerica_poly <- as.polygons(SthAmerica_ras, dissolve=T)

# AfroTropics
Afrotropics_ras <- mask(ras_timestep,Afrotropics)
values(Afrotropics_ras)[which(values(Afrotropics_ras$koppen %in% c(6)))] <- NA
values(Afrotropics_ras)[which(!is.na(values(Afrotropics_ras$koppen)))] <- 1
Afrotropics_poly <- as.polygons(Afrotropics_ras, dissolve=T)

# New Zealand
NewZealand_ras <- mask(ras_timestep,NewZealand)
values(NewZealand_ras)[which(values(NewZealand_ras$koppen %in% c(6)))] <- NA
values(NewZealand_ras)[which(!is.na(values(NewZealand_ras$koppen)))] <- 1
NewZealand_poly <- as.polygons(NewZealand_ras, dissolve=T)

# New Caledonia
NewCaledonia_ras <- crop(ras_timestep,NewCaledonia)
values(NewCaledonia_ras)[which(!is.na(values(NewCaledonia_ras$koppen)))] <- 1
NewCaledonia_poly <- as.polygons(NewCaledonia_ras, dissolve=T)


# Ausralo-Papua
AustraloPapuan_ras <- mask(ras_timestep, AustraloPapua)
values(AustraloPapuan_ras)[which(values(AustraloPapuan_ras$koppen) %in% c(6))] <- NA

AustraloPapuan_ras_trop    <- AustraloPapuan_ras
AustraloPapuan_ras_subtrop <- AustraloPapuan_ras
AustraloPapuan_ras_semi    <- AustraloPapuan_ras
AustraloPapuan_ras_med    <- AustraloPapuan_ras
AustraloPapuan_ras_arid    <- AustraloPapuan_ras

values(AustraloPapuan_ras_trop    )[which(!values(AustraloPapuan_ras_trop    )==1)] <- NA
values(AustraloPapuan_ras_subtrop )[which(!values(AustraloPapuan_ras_subtrop )==2)] <- NA
values(AustraloPapuan_ras_semi    )[which(!values(AustraloPapuan_ras_semi    )==4)] <- NA
values(AustraloPapuan_ras_med     )[which(!values(AustraloPapuan_ras_med     )==3)] <- NA

AustraloPapuan_poly <- as.polygons(AustraloPapuan_ras,    dissolve=T)

AustraloPapuan_poly_trop    <- as.polygons(AustraloPapuan_ras_trop,    dissolve=T)
AustraloPapuan_poly_subtrop <- as.polygons(AustraloPapuan_ras_subtrop, dissolve=T)
AustraloPapuan_poly_semi    <- as.polygons(AustraloPapuan_ras_semi,    dissolve=T)
AustraloPapuan_poly_med     <- as.polygons(AustraloPapuan_ras_med ,    dissolve=T)

# Antarctica
Antarctica_ras <- mask(ras_timestep, Antarctica)
values(Antarctica_ras)[which(!is.na(values(Antarctica_ras$koppen)))] <- 1
Antarctica_poly <- as.polygons(Antarctica_ras, dissolve=T)
Antarctica_poly <- erase(Antarctica_poly, AustraloPapuan_poly)
Antarctica_poly <- erase(Antarctica_poly, Afrotropics_poly)
Antarctica_poly <- erase(Antarctica_poly, SthAmerica_poly)
Antarctica_poly <- erase(Antarctica_poly, Madagascar_poly)
Antarctica_poly <- erase(Antarctica_poly, NewZealand_poly )
Antarctica_poly <- erase(Antarctica_poly, NewCaledonia_poly)

plot_it <- T
if(plot_it == TRUE){
  par(mfrow=c(1,1))
  plot(ras_timestep)
  plot(TropicalAsia_poly, add=T)
  plot(Afrotropics_poly, add=T)
  plot(AustraloPapuan_poly, add=T)
  plot(Antarctica_poly, add=T)
  plot(SthAmerica_poly, add=T)
  plot(Madagascar_poly, add=T)
  plot(NewZealand_poly, add=T)
  plot(NewCaledonia_poly, add=T)
}


RegionalPolygons[[5]][["TropicalAsia"]]                <- TropicalAsia_poly
RegionalPolygons[[5]][["Cape"]]                        <- Afrotropics_poly
RegionalPolygons[[5]][["Madagascar"]]                  <- Madagascar_poly
RegionalPolygons[[5]][["SthAmerica"]]                  <- SthAmerica_poly
RegionalPolygons[[5]][["Antarctica"]]                  <- Antarctica_poly
RegionalPolygons[[5]][["NewZealand"]]                  <- NewZealand_poly
RegionalPolygons[[5]][["NewCaledonia"]]                <- NewCaledonia_poly
RegionalPolygons[[5]][["AustraloPapuan_Tropical"]]     <- AustraloPapuan_poly_trop
RegionalPolygons[[5]][["AustraloPapuan_Subtropical"]]  <- AustraloPapuan_poly_subtrop
RegionalPolygons[[5]][["AustraloPapuan_SemiArid"]]     <- AustraloPapuan_poly_semi
RegionalPolygons[[5]][["AustraloPapuan_Mediterranean"]]<- AustraloPapuan_poly_med
RegionalPolygons[[5]][["AustraloPapuan_Arid"]]         <- NA

### ~~~~ Eocene - Miocene 40-20 ~~~###

# load roughly drawn continent outlines from QGIS
TropicalAsia  <- vect("External_datasets/QGIS_regional_polygons/20Ma_TropicalAsia.shp")
Afrotropics   <- vect("External_datasets/QGIS_regional_polygons/20Ma_Afrotropics.shp")
Neotropics    <- vect("External_datasets/QGIS_regional_polygons/20Ma_SouthAmerica.shp")
AustraloPapua <- vect("External_datasets/QGIS_regional_polygons/20Ma_AustraloPapua.shp")
Antarctica    <- vect("External_datasets/QGIS_regional_polygons/20Ma_Ant.shp")
Madagascar    <- vect("External_datasets/QGIS_regional_polygons/20Ma_Madagascar.shp")
NewZealand    <- vect("External_datasets/QGIS_regional_polygons/20Ma_NZ.shp")
NewCaledonia  <- vect("External_datasets/QGIS_regional_polygons/20Ma_NC.shp")

# subset suitable habitat in each
li_years
ras_timestep <- koppen_li[[which(li_years==20)]]
names(ras_timestep) <- "koppen"

# Tropical Asia
TropicalAsia_ras <- mask(ras_timestep, TropicalAsia)
values(TropicalAsia_ras)[which(values(TropicalAsia_ras$koppen %in% c(3,4,5,6)))] <- NA
values(TropicalAsia_ras)[which(!is.na(values(TropicalAsia_ras$koppen)))] <- 1
TropicalAsia_poly <- as.polygons(TropicalAsia_ras, dissolve=T)

# Madagascar
Madagascar_ras <- mask(ras_timestep,Madagascar)
values(Madagascar_ras)[which(values(Madagascar_ras$koppen %in% c(6)))] <- NA
values(Madagascar_ras)[which(!is.na(values(Madagascar_ras$koppen)))] <- 1
Madagascar_poly <- as.polygons(Madagascar_ras, dissolve=T)


# Tropical South America
SthAmerica_ras <- mask(ras_timestep, Neotropics)
values(SthAmerica_ras)[which(values(SthAmerica_ras$koppen %in% c(6)))] <- NA
values(SthAmerica_ras)[which(!is.na(values(SthAmerica_ras$koppen)))] <- 1
SthAmerica_poly <- as.polygons(SthAmerica_ras, dissolve=T)

# AfroTropics
Afrotropics_ras <- mask(ras_timestep,Afrotropics)
values(Afrotropics_ras)[which(values(Afrotropics_ras$koppen %in% c(6)))] <- NA
values(Afrotropics_ras)[which(!is.na(values(Afrotropics_ras$koppen)))] <- 1
Afrotropics_poly <- as.polygons(Afrotropics_ras, dissolve=T)

# New Zealand
NewZealand_ras <- mask(ras_timestep,NewZealand)
values(NewZealand_ras)[which(values(NewZealand_ras$koppen %in% c(6)))] <- NA
values(NewZealand_ras)[which(!is.na(values(NewZealand_ras$koppen)))] <- 1
NewZealand_poly <- as.polygons(NewZealand_ras, dissolve=T)

# New Caledonia
NewCaledonia_ras <- crop(ras_timestep,NewCaledonia)
values(NewCaledonia_ras)[which(!is.na(values(NewCaledonia_ras$koppen)))] <- 1
NewCaledonia_poly <- as.polygons(NewCaledonia_ras, dissolve=T)


# Ausralo-Papua
AustraloPapuan_ras <- mask(ras_timestep, AustraloPapua)
values(AustraloPapuan_ras)[which(values(AustraloPapuan_ras$koppen) %in% c(6))] <- NA

AustraloPapuan_ras_trop    <- AustraloPapuan_ras
AustraloPapuan_ras_subtrop <- AustraloPapuan_ras
AustraloPapuan_ras_med    <- AustraloPapuan_ras

values(AustraloPapuan_ras_trop    )[which(!values(AustraloPapuan_ras_trop    )==1)] <- NA
values(AustraloPapuan_ras_subtrop )[which(!values(AustraloPapuan_ras_subtrop )==2)] <- NA
values(AustraloPapuan_ras_med     )[which(!values(AustraloPapuan_ras_med     )==3)] <- NA

AustraloPapuan_poly <- as.polygons(AustraloPapuan_ras,    dissolve=T)


AustraloPapuan_poly_trop    <- as.polygons(AustraloPapuan_ras_trop,    dissolve=T)
AustraloPapuan_poly_subtrop <- as.polygons(AustraloPapuan_ras_subtrop, dissolve=T)
AustraloPapuan_poly_med     <- as.polygons(AustraloPapuan_ras_med ,    dissolve=T)

# Antarctica
Antarctica_ras <- mask(ras_timestep, Antarctica)
values(Antarctica_ras)[which(!is.na(values(Antarctica_ras$koppen)))] <- 1
Antarctica_poly <- as.polygons(Antarctica_ras, dissolve=T)
Antarctica_poly <- erase(Antarctica_poly, AustraloPapuan_poly)
Antarctica_poly <- erase(Antarctica_poly, Afrotropics_poly)
Antarctica_poly <- erase(Antarctica_poly, SthAmerica_poly)
Antarctica_poly <- erase(Antarctica_poly, Madagascar_poly)
Antarctica_poly <- erase(Antarctica_poly, NewZealand_poly )
Antarctica_poly <- erase(Antarctica_poly, NewCaledonia_poly)

plot_it <- T
if(plot_it == TRUE){
  par(mfrow=c(1,1))
  plot(ras_timestep)
  plot(TropicalAsia_poly, add=T)
  plot(Afrotropics_poly, add=T)
  plot(AustraloPapuan_poly, add=T)
  plot(Antarctica_poly, add=T)
  plot(SthAmerica_poly, add=T)
  plot(Madagascar_poly, add=T)
  plot(NewZealand_poly, add=T)
  plot(NewCaledonia_poly, add=T)
}


RegionalPolygons[[6]][["TropicalAsia"]]                <- TropicalAsia_poly
RegionalPolygons[[6]][["Cape"]]                        <- Afrotropics_poly
RegionalPolygons[[6]][["Madagascar"]]                  <- Madagascar_poly
RegionalPolygons[[6]][["SthAmerica"]]                  <- SthAmerica_poly
RegionalPolygons[[6]][["Antarctica"]]                  <- Antarctica_poly
RegionalPolygons[[6]][["NewZealand"]]                  <- NewZealand_poly
RegionalPolygons[[6]][["NewCaledonia"]]                <- NewCaledonia_poly
RegionalPolygons[[6]][["AustraloPapuan_Tropical"]]     <- AustraloPapuan_poly_trop
RegionalPolygons[[6]][["AustraloPapuan_Subtropical"]]  <- AustraloPapuan_poly_subtrop
RegionalPolygons[[6]][["AustraloPapuan_SemiArid"]]     <- NA
RegionalPolygons[[6]][["AustraloPapuan_Mediterranean"]]<- AustraloPapuan_poly_med
RegionalPolygons[[6]][["AustraloPapuan_Arid"]]         <- NA


### ~~~~ Miocene - Present-day 20-0 Ma ~~~###

# Do a final time step based on the Present-day high resolution from Beck

# koppen data
koppen_083 <- rast("Beck_KG_V1_present_0p083.tif")

# as data frame
koppen_df <- as.data.frame(koppen_083, xy=TRUE)

# simplify biomes according to BGB
koppen_df$biome <- NA
koppen_df$biome[which(koppen_df$Beck_KG_V1_present_0p083 %in% c(0))] <- NA
koppen_df$biome[which(koppen_df$Beck_KG_V1_present_0p083 %in% c(29,30))] <- "Polar"
koppen_df$biome[which(koppen_df$Beck_KG_V1_present_0p083 %in% c(1,2,3,11))] <- "Tropical"
koppen_df$biome[which(koppen_df$Beck_KG_V1_present_0p083 %in% c(14, 15, 12, 16, 17, 18,19, 25, 26, 27))] <- "Subtropical"
koppen_df$biome[which(koppen_df$Beck_KG_V1_present_0p083 %in% c(8, 9, 10))] <- "Mediterranean"
koppen_df$biome[which(koppen_df$Beck_KG_V1_present_0p083 %in% c(6, 7))] <- "SemiArid"
koppen_df$biome[which(koppen_df$Beck_KG_V1_present_0p083 %in% c(4,5))] <- "Arid"
koppen_df$biome <- factor(koppen_df$biome, levels=c("Tropical", "Subtropical","Mediterranean","SemiArid", "Arid",   "Polar"))
koppen_df$koppen <- as.numeric(koppen_df$biome)

ras_timestep <- rast(koppen_df[,c("x", "y", "koppen")], type="xyz")

plot(ras_timestep)

# draw bounding boxes
TropicalAsia  <- vect(ext(c(87.7138018628281, 128.865368331922, -12.3268512974667, 23.8276836329068)))
Afrotropics   <- vect(ext(c(3.86138613861386, 42.4752475247525, -39.5256498134822, 5.7168176334353)))
Neotropics    <- vect(ext(c(-83.7582970931563, -33.4950789654383, -60.0106284372932, 12.4129010948191)))
NorthAmerica  <- vect(ext(c(-128.665598535134, -68.9265278095674, 21.8773396132201, 48.6246658608752)))
AustraloPapua <- vect(ext(c(109.0547035935, 160.141908903639, -49.3116979382312, 3)))
Antarctica    <- vect(ext(c(-176.943521594684, 174.551495016611, -86.3787375415282, -69.1029900332226)))
NewZealand    <- vect(ext(c(161.377889677272, 178.681620508125, -55.0726605146492, -31.6173128820901)))
NewCaledonia  <- vect(ext(c(162.591362126246, 171.229235880399, -23.9202657807309, -17.9401993355482)))
Madagascar    <- vect(ext(c(44.7029702970297, 52.8712871287129, -28.4004529003058, -11.3418176334352)))

# subset suitable habitat in each
# Tropical Asia
TropicalAsia_ras <- mask(ras_timestep, TropicalAsia)
values(TropicalAsia_ras)[which(values(TropicalAsia_ras$koppen %in% c(3,4,5,6)))] <- NA
values(TropicalAsia_ras)[which(!is.na(values(TropicalAsia_ras$koppen)))] <- 1
TropicalAsia_poly <- as.polygons(TropicalAsia_ras, dissolve=T)

# Madagascar
Madagascar_ras <- mask(ras_timestep,Madagascar)
values(Madagascar_ras)[which(values(Madagascar_ras$koppen %in% c(6)))] <- NA
values(Madagascar_ras)[which(!is.na(values(Madagascar_ras$koppen)))] <- 1
Madagascar_poly <- as.polygons(Madagascar_ras, dissolve=T)


# Tropical South America
SthAmerica_ras <- mask(ras_timestep, Neotropics)
values(SthAmerica_ras)[which(values(SthAmerica_ras$koppen %in% c(6)))] <- NA
values(SthAmerica_ras)[which(!is.na(values(SthAmerica_ras$koppen)))] <- 1
SthAmerica_poly <- as.polygons(SthAmerica_ras, dissolve=T)

# AfroTropics
Afrotropics_ras <- mask(ras_timestep,Afrotropics)
values(Afrotropics_ras)[which(values(Afrotropics_ras$koppen %in% c(6)))] <- NA
values(Afrotropics_ras)[which(!is.na(values(Afrotropics_ras$koppen)))] <- 1
Afrotropics_poly <- as.polygons(Afrotropics_ras, dissolve=T)

# New Zealand
NewZealand_ras <- mask(ras_timestep,NewZealand)
values(NewZealand_ras)[which(values(NewZealand_ras$koppen %in% c(6)))] <- NA
values(NewZealand_ras)[which(!is.na(values(NewZealand_ras$koppen)))] <- 1
NewZealand_poly <- as.polygons(NewZealand_ras, dissolve=T)

# New Caledonia
NewCaledonia_ras <- crop(ras_timestep,NewCaledonia)
values(NewCaledonia_ras)[which(!is.na(values(NewCaledonia_ras$koppen)))] <- 1
NewCaledonia_poly <- as.polygons(NewCaledonia_ras, dissolve=T)

# Australo-Papua
AustraloPapuan_ras <- mask(ras_timestep, AustraloPapua)
values(AustraloPapuan_ras)[which(values(AustraloPapuan_ras$koppen) %in% c(6))] <- NA

AustraloPapuan_ras_trop    <- AustraloPapuan_ras
AustraloPapuan_ras_subtrop <- AustraloPapuan_ras
AustraloPapuan_ras_semi    <- AustraloPapuan_ras
AustraloPapuan_ras_med    <- AustraloPapuan_ras
AustraloPapuan_ras_arid    <- AustraloPapuan_ras

values(AustraloPapuan_ras_trop    )[which(!values(AustraloPapuan_ras_trop    )==1)] <- NA
values(AustraloPapuan_ras_subtrop )[which(!values(AustraloPapuan_ras_subtrop )==2)] <- NA
values(AustraloPapuan_ras_semi    )[which(!values(AustraloPapuan_ras_semi    )==4)] <- NA
values(AustraloPapuan_ras_med     )[which(!values(AustraloPapuan_ras_med     )==3)] <- NA
values(AustraloPapuan_ras_arid     )[which(!values(AustraloPapuan_ras_arid  )==5)] <- NA

AustraloPapuan_poly <- as.polygons(AustraloPapuan_ras,    dissolve=T)
AustraloPapuan_poly_trop    <- as.polygons(AustraloPapuan_ras_trop,    dissolve=T)
AustraloPapuan_poly_subtrop <- as.polygons(AustraloPapuan_ras_subtrop, dissolve=T)
AustraloPapuan_poly_semi    <- as.polygons(AustraloPapuan_ras_semi,    dissolve=T)
AustraloPapuan_poly_med     <- as.polygons(AustraloPapuan_ras_med ,    dissolve=T)
AustraloPapuan_poly_arid     <- as.polygons(AustraloPapuan_ras_arid ,    dissolve=T)

# Antarctica
Antarctica_ras <- mask(ras_timestep, Antarctica)
values(Antarctica_ras)[which(!is.na(values(Antarctica_ras$koppen)))] <- 1
Antarctica_poly <- as.polygons(Antarctica_ras, dissolve=T)
Antarctica_poly <- erase(Antarctica_poly, AustraloPapuan_poly)
Antarctica_poly <- erase(Antarctica_poly, Afrotropics_poly)
Antarctica_poly <- erase(Antarctica_poly, SthAmerica_poly)
Antarctica_poly <- erase(Antarctica_poly, Madagascar_poly)
Antarctica_poly <- erase(Antarctica_poly, NewZealand_poly )
Antarctica_poly <- erase(Antarctica_poly, NewCaledonia_poly)

plot_it <- T
if(plot_it == TRUE){
  par(mfrow=c(1,1))
  plot(ras_timestep)
  plot(TropicalAsia_poly, add=T)
  plot(Afrotropics_poly, add=T)
  plot(AustraloPapuan_poly, add=T)
  plot(Antarctica_poly, add=T)
  plot(SthAmerica_poly, add=T)
  plot(Madagascar_poly, add=T)
  plot(NewZealand_poly, add=T)
  plot(NewCaledonia_poly, add=T)
}


RegionalPolygons[[7]][["TropicalAsia"]]                <- TropicalAsia_poly
RegionalPolygons[[7]][["Cape"]]                        <- Afrotropics_poly
RegionalPolygons[[7]][["Madagascar"]]                  <- Madagascar_poly
RegionalPolygons[[7]][["SthAmerica"]]                  <- SthAmerica_poly
RegionalPolygons[[7]][["Antarctica"]]                  <- Antarctica_poly
RegionalPolygons[[7]][["NewZealand"]]                  <- NewZealand_poly
RegionalPolygons[[7]][["NewCaledonia"]]                <- NewCaledonia_poly
RegionalPolygons[[7]][["AustraloPapuan_Tropical"]]     <- AustraloPapuan_poly_trop
RegionalPolygons[[7]][["AustraloPapuan_Subtropical"]]  <- AustraloPapuan_poly_subtrop
RegionalPolygons[[7]][["AustraloPapuan_SemiArid"]]     <- AustraloPapuan_poly_semi
RegionalPolygons[[7]][["AustraloPapuan_Mediterranean"]]<- AustraloPapuan_poly_med
RegionalPolygons[[7]][["AustraloPapuan_Arid"]]         <- AustraloPapuan_poly_arid


# paleo-area


areas_df <- data.frame(do.call(rbind, lapply(RegionalPolygons, getArea)))
areas_df$time <- seq(from=120, to=0, by=-10)
area_df_melt <- melt(areas_df, id="time")

write.csv(area_df_melt, "biome_paleoarea.csv")
