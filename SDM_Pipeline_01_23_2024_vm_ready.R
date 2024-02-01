# Habitat Suitability Models (HSMs)
# Nikol Dimitrov, Richard Scuster and Juan Zuloaga
# Pipeline - May 2023 /
# set Java options here before running any code or packages
options(java.parameters = "-Xmx64g" ) #for 256GB ram machine 
gc()
######### START PIPELINE ######################

# 0. Packages required -----------

list.of.packages <- c("spatialEco", "devtools", "rJava", "raster", "rasterVis", "devtools", "ggplot2", "countrycode", "CoordinateCleaner", "rgbif", "raster", 'dismo', 'ENMeval', 'dplyr', "ecospat",
                      'adehabitatHR', "sf", "spThin", "WorldClimTiles", "virtualspecies", "rnaturalearthdata")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

devtools::install_github("kapitzas/WorldClimTiles") # to load and merge BIOCLIM tiles
# (https://github.com/kapitzas/WorldClimTiles/blob/master/README.md)


if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only =TRUE) 

system.file("java", package="dismo")

# 1. Settings ----------

# 1.1. Creating main folder

# Checking of folder exists
if(!dir.exists("C:/HSM_NCC")){
  dir.create("C:/HSM_NCC")
  print("dir created")
}else{
  print("dir exists")
}

# Setting main folder as a default
setwd("C:/HSM_NCC")

# 1.2. Creating a directory to save temporary files that will be deleted after process

if(!dir.exists("./temp_HMs_to_REMOVE")){
  dir.create("./temp_HMs_to_REMOVE")
  print("dir created")
}else{
  print("dir exists")
}

# Creating directories to save results 

if(!dir.exists("./Results/Variable Importance")){
  dir.create("./Results/Variable Importance")
  print("dir created")
}else{
  print("dir exists")
}


if(!dir.exists("./Results")){
  dir.create("./Results")
  print("dir created")
}else{
  print("dir exists")
}

if(!dir.exists("./Results/Continuous_Results")){
  dir.create("./Results/Continuous_Results")
  print("dir created")
}else{
  print("dir exists")
}

if(!dir.exists("./Results/Binary_Results")){
  dir.create("./Results/Binary_Results")
  print("dir created")
}else{
  print("dir exists")
}

if(!dir.exists("./Results/Results_for_calculations")){
  dir.create("./Results/Results_for_calculations")
  print("dir created")
}else{
  print("dir exists")
}

if(!dir.exists("./Results/Model_stats")){
  dir.create("./Results/Model_stats")
  print("dir created")
}else{
  print("dir exists")
} 

if(!dir.exists("./Results/Model_issues")){
  dir.create("./Results/Model_stats")
  print("dir created")
}else{
  print("dir exists")
} 


# 1.3. Setting temp directory
rasterOptions(tmpdir=file.path("./temp_HMs_to_REMOVE"))

species_to_run <- c("Acronicta albarufa", "Adiantum aleuticum var. subpumilum", "Abronia latifolia")

#for (i in species_to_run){
  
# set start time 
  
start_time <- Sys.time()
  
# version.string R version 
R.Version()

# 1.4. Setting Projection to preserve area for all files 

# Lon-Lat 
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"

# Albers Equal Area Conic
#[See:](https://spatialreference.org/ref/esri/102001/proj4/) 
aeac="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Sourcing species cleanup functions 

source("Species_data_cleanup_functions.R")

# 2.Occurrences ------------

# 2.1. Loading data from GBIF (https://www.gbif.org/)

myspecies <- "Anaxyrus cognatus"


sp <- export_gbif_data(myspecies)

if(is.null(sp)){
  dat <- data.frame()
  species_name <- myspecies
  date_run <- as.character(Sys.Date())
  data_point_count <- "NULL"
  issue <- "No species data on GBIF"
  stage <- "GBIF export"
  n <- data.frame(cbind(species_name, date_run, data_point_count, issue, stage))
  dat <- rbind(dat, n)
  write.table(dat, file= "./Results/Model_issues/Reported_model_issues.csv", append = T, sep=',', row.names = F, col.names = F)
  #write_csv(dat, paste0("./Results/Model_issues/", species_name, ".csv"))
  next
}

if (nrow(sp) < 5){
  dat <- data.frame()
  species_name <- myspecies
  date_run <- as.character(Sys.Date())
  data_point_count <- nrow(sp)
  issue <- "Less than 5 observations - unable to run"
  stage <- "GBIF export"
  n <- data.frame(cbind(species_name, date_run, data_point_count, issue, stage))
  dat <- rbind(dat, n)
  write.table(dat, file= "./Results/Model_issues/Reported_model_issues.csv", append = T, sep=',', row.names = F, col.names = F)
  #write_csv(dat, paste0("./Results/Model_issues/", species_name, ".csv"))
  next
}

# 2.2. plotting uncleaned data 

plot_gbif_data(sp)

# 2.3. clean data 

cl_species <- clean_gbif_data(sp)

if (nrow(cl_species) < 5){
  dat <- data.frame()
  species_name <- myspecies
  date_run <- as.character(Sys.Date())
  data_point_count <- nrow(cl_species)
  issue <- "Less than 5 observations - unable to run"
  stage <- "Cleanup of species data"
  n <- data.frame(cbind(species_name, date_run, data_point_count, issue, stage))
  dat <- rbind(dat, n)
  write.table(dat, file= "./Results/Model_issues/Reported_model_issues.csv", append = T, sep=',', row.names = F, col.names = F)
  #write_csv(dat, paste0("./Results/Model_issues/", species_name, ".csv"))
  next
}

plot_gbif_data(cl_species)

# 2.4. Spatial thining

obs_thin <- thin_data(cl_species)

if (nrow(obs_thin) < 5){
  dat <- data.frame()
  species_name <- myspecies
  date_run <- as.character(Sys.Date())
  data_point_count <- nrow(obs_thin)
  issue <- "Less than 5 observations - unable to run"
  stage <- "Thinning of species data"
  n <- data.frame(cbind(species_name, date_run, data_point_count, issue, stage))
  dat <- rbind(dat, n)
  write.table(dat, file= "./Results/Model_issues/Reported_model_issues.csv", append = T, sep=',', row.names = F, col.names = F)
  #write_csv(dat, paste0("./Results/Model_issues/", species_name, ".csv"))
  next
}

# 2.5. spatial points 

Obs_gbif_data_sp <- convert_spatial_points(obs_thin)

Obs_gbif_data_sp

# 3. Area of study -------------

# Creating box extent to download predictors

# Several options here:
# Geographic range: not for many species
# Minimum convex polygon (MVP) with buffer

# We will implement two options; (i) if range map exist read .shp file;
#                               (ii) if not create MCP with buffer     

# 3.1. Loading species list table

species_rangemap <- read.csv("./Data/species_list/updated_species_list_taxonomy.csv")
species_rangemap
myspecies_rangemap <- species_rangemap[species_rangemap$Sci.Name == myspecies, ]

if(myspecies_rangemap$Range == "Yes"){
  
  # 3.2. Range map
  
  # Box extent
  box_extent_analysis <- st_read(paste0("./Data/ECCC_Polygons/",myspecies, ".shp"))
  
  box_extent_analysis_bf <- box_extent_analysis %>% 
    st_as_sf() %>%
    st_transform(crs = wgs84)
  
  # Transform for plotting
  box_extent_analysis_bf_aeac <- st_as_sf(box_extent_analysis_bf)%>%
    st_transform(aeac)
  
}else{
  
  # 3.3 MinImum convex polygon (mcp) with buffer  
  
  # Box extent
  box_extent_analysis <- mcp(Obs_gbif_data_sp, percent = 100)%>%
    st_as_sf()%>%
    st_transform(crs=aeac)
  
  # Buffering box extent
  box_extent_analysis_bf <-  st_buffer(box_extent_analysis, dist =  100000)%>%
    st_transform(crs=wgs84)
  
  # Transform for plotting
  box_extent_analysis_bf_aeac <- st_as_sf(box_extent_analysis_bf)%>%
    st_transform(aeac)
}

# 4. Predictors ------------

# 4.1. Topographic heterogeneity

# Vector_Ruggedness_Measure
vrm <- raster("./Data/Topographic_Index/vrm_1KMmn_GMTEDmd.tif") %>% 
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(vrm) <- "Vector_Ruggedness_Measure"

# Roughness
roughness <- raster("./Data/Topographic_Index/roughness_1KMmn_GMTEDmd.tif") %>%
  crop(box_extent_analysis_bf)   %>%
  resample(vrm, method = "ngb") %>% 
  mask(box_extent_analysis_bf)
names(roughness) <- "Roughness"

# Slope
Slope <- raster("./Data/Topographic_Index/slope_1KMmn_GMTEDmd.tif") %>%
  crop(box_extent_analysis_bf)  %>%
  resample(vrm, method = "ngb") %>%
  mask(box_extent_analysis_bf) 
names(Slope) <- "Slope"

# Eastness
Eastness <- raster("./Data/Topographic_Index/eastness_1KMmn_GMTEDmd.tif") %>%
  crop(box_extent_analysis_bf)  %>%
  resample(vrm, method = "ngb") %>% 
  mask(box_extent_analysis_bf)
names(Eastness) <- "Eastness"

# Northness
Northness <- raster("./Data/Topographic_Index/northness_1KMmn_GMTEDmd.tif") %>%
  crop(box_extent_analysis_bf)  %>%
  resample(vrm, method = "ngb") %>% 
  mask(box_extent_analysis_bf)
names(Northness) <- "Northness"


# 4.2. Bioclim

# Idenfitying tiles based on Area of study

box_extent_bioclim <- tile_name(as(box_extent_analysis_bf, 'Spatial'), "worldclim") # determine which WorldClim tiles your study area intersects with.

# Creating folder
if(dir.exists("./bioclim_t")){
}else{
  out_bioclim <-   dir.create("./bioclim_t")}

# Loading tiles
clim_tiles <- tile_get(box_extent_bioclim, name =  "worldclim", var="bio", path = "./bioclim_t") # for 0.5 arcmin worldclim tiles of 

# Merging tiles
clim_tiles_merge <- tile_merge(clim_tiles) %>%
   resample(vrm, method = "ngb") %>% 
   mask(box_extent_analysis_bf)


# 4.3. Productivity

# Loading Dynamic Habitat Index  (all .tiff files)
fpar_all <- list.files("./Data/fpar_can", pattern=".tif$", full.names = T)

#  Band 1
fpar_b1_mean <- lapply(fpar_all, raster, band=1) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf) %>%  # Crop to extent
  calc(mean)%>%                                                # Apply function to calculate mean
  resample(vrm, method = "ngb") %>% 
  mask(box_extent_analysis_bf)                      # resample using previous raster
names(fpar_b1_mean) <- "Cummulative_annual_productivity_b1"    # rename variable

#  Band 2
fpar_b2_mean <- lapply(fpar_all, raster, band=2) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf)%>%  # Crop to extent
  calc(mean)%>%                                                # Apply function to calculate mean
  resample(vrm, method = "ngb") %>%
  mask(box_extent_analysis_bf)                                 # resample using previous raster
names(fpar_b2_mean) <- "Minimum_annual_productivity_b2"        # rename variable

#  Band 3
fpar_b3_mean <- lapply(fpar_all, raster, band=3) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf)%>%                                # Crop to extent
  calc(mean)%>%                                                # Apply function to calculate mean
  resample(vrm, method = "ngb") %>%
  mask(box_extent_analysis_bf)                                 # resample using previous raster
names(fpar_b3_mean) <- "Variation_annual_productivity_b3"      # rename variable

# 4.4. Water proximity

# Percentage of lakes

lakes_100m <- raster("./Data/Lakes/Lakes_100_wgs84.tif")%>% # Lakes North America 100m resolution
  crop(box_extent_analysis_bf) 
  

lakes_100m_agg <- aggregate(lakes_100m, fact=10, fun=sum)%>% # aggregating 100m cells into 1000m cells, using 'sum' function
  resample(vrm, method = "ngb") %>%
  mask(box_extent_analysis_bf)
names(lakes_100m_agg) <- "Lakes_percentage"

# 4.5 Landcover data 

# Bare areas
bare_areas <- raster("./Data/Landcover Change/Bare areas reclass/bare_areas_1km_change_1992_2020.tif") %>% 
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(bare_areas) <- "Bare_areas"

# Cropland 
Cropland <- raster("./Data/Landcover Change/Cropland reclass/cropland_1km_change_1992_2020.tif") %>%
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(Cropland) <- "Cropland"

# Flooded vegetation 
Flooded_vegetation <- raster("./Data/Landcover Change/Flooded vegetation reclass/flooded_vegetation_1km_change_1992_2020.tif") %>%
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(Flooded_vegetation) <- "Flooded_vegetation"

# Grassland 
Grassland <- raster("./Data/Landcover Change/Grassland reclass/grassland_1km_change_1992_2020.tif") %>%                                         # Stack them
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(Grassland) <- "Grassland"

# Other vegetation 
Other_vegetation <- raster("./Data/Landcover Change/Other vegetation reclass/other_vegetation_1km_change_1992_2020.tif") %>%
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(Other_vegetation) <- "Other_vegetation"

# Snow and ice 
Snow_ice <- raster("./Data/Landcover Change/Snow and ice reclass/snow_ice_1km_change_1992_2020.tif") %>%                                        # Stack them
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(Snow_ice) <- "Snow_and_ice"

# Tree cover 
Tree_cover <- raster("./Data/Landcover Change/Forests reclass/forests_1km_change_1992_2020.tif") %>%
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(Tree_cover) <- "Tree_cover"

# Urban areas 
Urban_areas<- raster("./Data/Landcover change/Urban reclass/urban_1km_change_1992_2020.tif")  %>%
  crop(box_extent_analysis_bf) %>% 
  mask(box_extent_analysis_bf)
names(Urban_areas) <- "Urban_areas"

# 4.6.  Stacking Predictors ----------

predictors <- stack(vrm,
                    roughness,
                    Slope,
                    Eastness,
                    Northness,
                    lakes_100m_agg,
                    if(maxValue(fpar_b1_mean) != 0){fpar_b1_mean}else{},
                    if(maxValue(fpar_b2_mean) != 0){fpar_b2_mean}else{},
                    if(maxValue(fpar_b3_mean) != 0){fpar_b3_mean}else{},
                    clim_tiles_merge, 
                    if(maxValue(bare_areas) != 0 ){bare_areas}else{}, 
                    if(maxValue(Cropland) != 0 ){Cropland}else{}, 
                    if(maxValue(Grassland) != 0 ){Grassland}else{},
                    if(maxValue(Flooded_vegetation) != 0 ){Flooded_vegetation}else{}, 
                    if(maxValue(Other_vegetation) != 0 ){Other_vegetation}else{},
                    if(maxValue(Snow_ice) != 0){Snow_ice}else{},
                    if(maxValue(Tree_cover) != 0 ){Tree_cover}else{},
                    if(maxValue(Urban_areas)!=0){Urban_areas}else{}, na.rm=T)

rm(vrm, roughness, Slope, Eastness, Northness, lakes_100m_agg, lakes_100m, fpar_b1_mean, fpar_b2_mean,
   fpar_b3_mean, clim_tiles_merge, clim_tiles, bare_areas, Cropland, Grassland, Flooded_vegetation, Other_vegetation,
   Snow_ice, Tree_cover, Urban_areas)

# 4.7. Removing collinear variables

# Calculating collinearity
collinearity_test <- removeCollinearity(predictors,
                                        multicollinearity.cutoff = 0.70, 
                                        plot = F, select.variables = T, sample.points = FALSE)

# Sub-setting variables
noncollinear_predictors <- stack(subset(predictors, collinearity_test))

noncollinear_predictors[[2]]

rm(predictors, box_extent_analysis, box_extent_analysis_bf_aeac)

# 5. Creating background points ----------
# Identify number of background points (bg) for HSMs within the geographic range (EBARs), testing three bg.

# 5.1. Calculating number of grid cells in predictor variable
pixel_1000Km2 <- (1000*1000)/1000000
study_area_size <- cellStats((((noncollinear_predictors[[1]]*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
study_area_size

# 5.2. Selecting bg (two options: in case small geographic ranges)

small_study_area <- 0.40 # take 40% if area of study contains less than 100,000 grid cells

huge_study_area <- 0.20 # take 20% if area of study contains  100,000 - 1,000,000 grid cells

very_huge_study_area <- 0.05 # take 5% if area of study contains more than 1,000,000 grid cells

extremely_huge_study_area <- 0.001 #take 1% if area of study contains more than 2,000,000 grid cells 

study_area_size

if(study_area_size < 100000){
  
  background_points <- c(
    round(study_area_size*small_study_area, 0)
  )
}else if(study_area_size >= 100000 & study_area_size < 1000000){
  
  background_points <- c(
    round(study_area_size*huge_study_area, 0)
  )
}else if(study_area_size >= 1000000 & study_area_size < 2000000){
  
  background_points <- c(
    round(study_area_size*very_huge_study_area, 0)
  )
}else if (study_area_size >= 2000000){
  
  background_points <- c(
    round(study_area_size*extremely_huge_study_area, 0)
  )} 

# 6. Creating sampling bias layer-----------

# Creating a sampling bias layer to select background point from these areas
# UsIng MASS::kde2d

#sub step - only for obs less than 25 


if (nrow(obs_thin) <= 25){
 x <- raster::extract(noncollinear_predictors, obs_thin)
 x
 z <- data.frame(obs_thin, x)
 z <- na.omit(z)
 if (nrow(z) < 5){
   dat <- data.frame()
   species_name <- myspecies
   date_run <- as.character(Sys.Date())
   data_point_count <- nrow(z)
   issue <- "Less than 5 observations - unable to run"
   stage <- "Removal of NA observations"
   n <- data.frame(cbind(species_name, date_run, data_point_count, issue, stage))
   dat <- rbind(dat, n)
   write.table(dat, file= "./Results/Model_issues/Reported_model_issues.csv", append = T, sep=',', row.names = F, col.names = F)
   #write_csv(dat, paste0("./Results/Model_issues/", species_name, ".csv"))
   next
 }

 obs_thin <- z %>% 
  dplyr::select("Longitude", "Latitude")
}

points_thin_sf <- st_as_sf(data.frame(obs_thin), coords=c("Longitude","Latitude"), crs=wgs84)


obs_density <- MASS::kde2d(st_coordinates(st_as_sf(points_thin_sf))[,1],
                           st_coordinates(st_as_sf(points_thin_sf))[,2],
                           lims = c(range(noncollinear_predictors[[1]]@extent@xmin,
                                          noncollinear_predictors[[1]]@extent@xmax),
                                    range(noncollinear_predictors[[1]]@extent@ymin,
                                          noncollinear_predictors[[1]]@extent@ymax))) %>%
  
  raster()

obs_density_a_wgs84 <- projectRaster(obs_density, crs = wgs84, res = 0.008333333)%>%
  resample(noncollinear_predictors[[1]])



# 7. Model fitting ------

# 7.1. Renaming observations columns
occs <- st_as_sf(x = obs_thin,                         
                 coords =c("Longitude","Latitude"), crs=wgs84)

occs_df <-  data.frame(st_coordinates(occs))%>%
  dplyr::rename(Longitude = X,  Latitude =  Y)

# 7.2. creating pseudo-absences

bg_points <- dismo::randomPoints(obs_density_a_wgs84, n = background_points) %>% as.data.frame()
colnames(bg_points) <- colnames(occs)

# 7.3. background using bias
bg_bias <- xyFromCell(!is.na(obs_density_a_wgs84),
                      sample(ncell(!is.na(obs_density_a_wgs84)),
                             nrow(bg_points),
                             prob =  values(!is.na(obs_density_a_wgs84))))
colnames(bg_bias) <- colnames(occs_df)

# 8. Modeling (using Maxent in ENMeval) basic parameters

# 8.1. Model settings

# 8.1.1. Selecting features based on number of occurrences ---------------
if(nrow(occs_df) <=10){
  meth = 'jackknife'
  features <- "L"
}else if(nrow(occs_df) >10 && nrow(occs_df) <=15){
  meth = 'jackknife'
  features <- c("L", "Q", "LQ")
}else if(nrow(occs_df) >15 && nrow(occs_df) <=25){
  meth = 'jackknife'
  features <- c("L", "Q", "H", "LQ", "QH")
}else if(nrow(occs_df) >25 && nrow(occs_df) <=80){
  meth = 'randomkfold'
  features <- c("L", "Q", "H", "LQ", "QH")
}else if(nrow(occs_df) > 80){
  meth = 'randomkfold'
  features <- c("Q",  "H", "LQ", "LQP", "QPT")
}
features
meth

# 8.1.2. Partition: select method based on number of occurrences

if(!is.null(occs_df)){
  
  if(nrow(occs_df) <=25){
    user_partition  <- get.jackknife(occs_df, bg_bias)
    
  }else if(nrow(occs_df) > 25){
    user_partition  <- get.randomkfold(occs_df, bg_bias, 10)
  }
}

# 8.1.3. Run model, using ENMeval package (Maxent algorithm)
gc()
parallel::detectCores()

model_species <- ENMeval::ENMevaluate(occs = occs_df,
                                      envs = noncollinear_predictors,
                                      bg = bg_bias, 
                                      algorithm = "maxent.jar",
                                      method = meth,
                                      partitions = 'user',
                                      user.grp = (list(occs.grp = user_partition$occs.grp,
                                                       bg.grp = user_partition$bg.grp)),
                                      tune.args = list(fc = features, rm = c(0.05,0.5, 1)), # let' test 3 values
                                      parallel =  FALSE,
                                      doClamp = T,
                                      updateProgress = TRUE,
                                      parallelType = "doParallel")


# 8.2. Model results 

# 8.2.1. All models
res_model_species <- model_species@results

# 8.2.2. let's add background point and observations
res_model_species$bg_points <- length(model_species@bg[,1])
res_model_species$occs_points <- length(model_species@occs.grp)

# 8.2.3. Let's now choose the optimal model settings based on the sequential criteria and examine it.
# lowest value of omission rate (or.10p.avg)
# to resolve tied or.10p.avg, we used highest AUC value (auc.val.avg)

optimal_model <- eval.results(model_species)%>%
  filter(!is.na(or.10p.avg))%>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))

optimal_model

# 8.2.4.Variable importance
mod.seq <- eval.models(model_species)[[optimal_model$tune.args]]

plot(mod.seq, type = "cloglog") 
variable_importance <- as.data.frame(model_species@variable.importance[[as.vector(optimal_model$tune.args[1])]])[,-2]
v <- variable_importance %>% 
  arrange(desc(permutation.importance))

write.csv(v, file = paste0("./Results/Variable Importance/", myspecies, "_variable_importance.csv"))
# 9. Predictions -------

all_model_prediction <- eval.predictions(model_species)
best_model_prediction <- eval.predictions(model_species)[[optimal_model$tune.args]]

plot(best_model_prediction) 
plot(occs, add=T)

# 10.  Model uncertainty ---------

# 10.1.Modelling
model_10 <- list()

for(i in 1:3){
  cat(paste0("Testing background points_model_", i), '\n')
  
  bg_points <- dismo::randomPoints(obs_density_a_wgs84, n = background_points) %>% as.data.frame()
  colnames(bg_points) <- colnames(obs_thin)
  
  # 10.2.background using bias
  bg_bias <- xyFromCell(!is.na(obs_density_a_wgs84),
                        sample(ncell(!is.na(obs_density_a_wgs84)),
                               nrow(bg_points),
                               prob =  values(!is.na(obs_density_a_wgs84))))
  colnames(bg_bias) <- colnames(occs_df)
  
  model_10[[i]] <- ENMevaluate(occs = occs_df,
                               envs = noncollinear_predictors,
                               bg = bg_bias, 
                               algorithm = 'maxent.jar',
                               partitions = 'block',
                               tune.args = list(fc = "L", rm = 1),
                               parallel =  FALSE,
                               updateProgress = TRUE,
                               parallelType = "doParallel"
  )
  
}

# 10.3. Calculating coefficient of variance
model_10_predictions <- stack(lapply(model_10, eval.predictions))
uncertainty <- raster::cv(model_10_predictions, na.rm=F)

# 11. Projecting outputs -----------
# Projecting to aeac to preserve areas and allow calculations 
model_species_prediction_p <- projectRaster(best_model_prediction, crs = aeac, res = 1000, method = "bilinear")
uncertainty_p <- projectRaster(uncertainty, crs = aeac, res = 1000, method = "bilinear")

par(mfrow=c(1,2))
plot(model_species_prediction_p)
plot(uncertainty_p)

par(mfrow=c(2,1))
plot(model_species_prediction_p)
plot(best_model_prediction)
plot(uncertainty_p)

plot(best_model_prediction)

plot(model_species_prediction_p)
dev.off()
# 12. Thresholding values:

# 90% of records 

pred_vals_presence <- raster::extract(best_model_prediction, occs)

pred_vals_presence <- pred_vals_presence[!is.na(pred_vals_presence)]

n10 <- ceiling(length(pred_vals_presence) * 0.1)

threshold <- sort(pred_vals_presence)[n10]

binarize <- function(x){
  ifelse(x <= threshold, 0,
  ifelse(x > threshold, 1, NA))
}

or.10.threshold.imposed <- calc(best_model_prediction, fun=binarize)
model_species_prediction_binary_p <- projectRaster(or.10.threshold.imposed, crs = aeac, res = 1000, method = "bilinear")

plot(or.10.threshold.imposed)
# 12. Write outputs -----------

# 12.1. Load and project into NCC national grid 

nat_grid <- raster("./Data/national_grid/constant_grid.tif")
species_prediction_nat_grid <- projectRaster(best_model_prediction, crs = proj4string(nat_grid), res = 1000, method = "bilinear")
species_uncertainty_nat_grid <- projectRaster(uncertainty, crs = proj4string(nat_grid), res = 1000, method = "bilinear")
binary_species_prediction_nat_grid <- projectRaster(or.10.threshold.imposed, crs = proj4string(nat_grid), res = 1000, method = "bilinear")
# 12.1.1 Write NCC projected rasters 
writeRaster(species_prediction_nat_grid, filename = paste0("./Results/Continuous_Results/", myspecies,"_continuous.tif"), options = c('TFW = YES'), overwrite = TRUE)
writeRaster(species_uncertainty_nat_grid, filename = paste0("./Results/Continuous_Results/",myspecies, "_uncertainty.tif"), options = c('TFW = YES'), overwrite = TRUE)
writeRaster(binary_species_prediction_nat_grid, filename = paste0("./Binary_Results/",myspecies, "_binary.tif"), options = c('TFW = YES'), overwrite = TRUE)

# 12.2. Write aeac projected rasters for calculations 
writeRaster(model_species_prediction_p, filename = paste0("./Results/Results_for_calculations/", myspecies,"_continuous.tif"), options = c('TFW = YES'), overwrite = TRUE)
writeRaster(uncertainty_p, filename = paste0("./Results/Results_for_calculations/",myspecies, "_uncertainty.tif"), options = c('TFW = YES'), overwrite = TRUE)
writeRaster(model_species_prediction_binary_p, filename= paste0("./Results/Results_for_calculations/",myspecies, "_binary.tif"), options = c('TFW = YES'), overwrite = TRUE)

# Report runtime, AUC, OR: 
end_time <- Sys.time()


dat <- data.frame()
species_name <- myspecies
date_run <- as.character(Sys.Date())
data_point_count <- nrow(obs_thin)
run_time <- difftime(end_time, start_time , units= 'hours')
auc_value <- optimal_model$auc.val.avg 
or_value <- optimal_model$or.10p.avg

if (optimal_model$auc.val.avg >= 0.75 & optimal_model$or.10p.avg <= 0.25){
  model_result <- "model informative, no issues"
}
if (optimal_model$auc.val.avg >= 0.75 & optimal_model$or.10p.avg > 0.25 & optimal_model$or.10p.avg <= 0.5){
  model_result <- "AUC normal, OR value between 0.25-0.5 - some concern in model predictions"
}
if (optimal_model$auc.val.avg < 0.75 & optimal_model$or.10p.avg > 0.25 & optimal_model$or.10p.avg <= 0.5){
  model_result <- "AUC below 0.75 - not informative, OR value between 0.25-0.5 - some concern in model predictions"
}
if (optimal_model$auc.val.avg >= 0.75 & optimal_model$or.10p.avg > 0.5){
  model_result <- "AUC normal, OR value greater than 0.5 - high concern in model predictions"
}
if (optimal_model$auc.val.avg < 0.75 & optimal_model$or.10p.avg > 0.5){
  model_result <- "AUC below 0.75 - not informative, OR value greater than 0.5 - high concern in model predictions"
}
if (optimal_model$auc.val.avg < 0.75 & optimal_model$or.10p.avg <= 0.25){
  model_result <- "AUC below 0.75 - not informative, OR normal" 
}

n <- data.frame(cbind(species_name, date_run, data_point_count, run_time, auc_value, or_value, model_result))
dat <- rbind(dat, n)
write.table(dat, file= "./Results/Model_stats/Reported_model_results.csv", append = T, sep=',', row.names = F, col.names = F)
#}

# remove temporary files 
temp_fold <- './temp_HMs_to_REMOVE'

# get all files in the directories, recursively
f <- list.files(temp_fold, include.dirs = F, full.names = T, recursive = T)
# remove the files
file.remove(f)
#############END PIPELINE ####################################
