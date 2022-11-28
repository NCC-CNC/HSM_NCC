
# Habitat Suitability Models (HSMs)
# Nikol Dimitrov, Richard Scuster and Juan Zuloaga
# Basic script

# set Java options here before running any code or packages

options(java.parameters = '-Xmx8g') 

# version.string R version 
R.Version()

######### START PIPELINE ######################

# 0. Packages required -----------

list.of.packages <- c("spatialEco", "rgdal", "devtools", "rJava", "rasterVis", "devtools", "ggplot2", "countrycoude",
                      "CoordinateCleaner", "rgbif", "raster", 'dismo', 'ENMeval', 'dplyr', "ecospat",
                      'adehabitatHR', "rgeos", "sf", "spThin", "WorldClimTiles", "virtualspecies", "rnaturalearthdata")
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

# 1.3. Setting temp directory
  rasterOptions(tmpdir=file.path("./temp_HMs_to_REMOVE"))

  
# 1.4. Setting Projection to preserve area for all files 
  
  # Lon-Lat 
   wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
   
  # Albers Equal Area Conic
  #[See:](https://spatialreference.org/ref/esri/102001/proj4/) 
    aeac="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"


# 2. Occurrences -----------

# 2.1. Loading data from GBIF (https://www.gbif.org/)

myspecies <- "Ammodramus bairdii"
sp <- export_gbif_data(myspecies)

# 2.2. plotting uncleaned data 

plot_gbif_data(sp)

# 2.3. clean data 

cl_species <- clean_gbif_data(sp)
             
plot_gbif_data(cl_species)

# 2.3.1. (sub-step for species where range exists) read range in and remove points that are outside of it 

range_eccc <- readOGR(paste0("./Data/range_maps/ECCC_Ranges_Polygons/", chartr(" ", "_", myspecies), ".shp"))
range_eccc$species <- paste0(myspecies)

range_flags <- cc_iucn(x = cl_species, range = range_eccc, lon = "decimalLongitude", lat = "decimalLatitude", value = "flagged")

cl_species_f <- cl_species[range_flags,]

# 2.4. Spatial thining

obs_thin <- thin_data(cl_species_f)

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

# 3.1 Loading species list table
  
  species_rangemap <- read.csv("./Data/species_list/species_list.csv")
  myspecies_rangemap <- species_rangemap[species_rangemap$Species.Name == myspecies, ]
  
  
if(myspecies_rangemap$Range == "Yes"){
  
 # 3.2 Range map
  
  # Box extent
  box_extent_analysis <- st_read(paste0("./Data/range_maps/", chartr(" ", "_", myspecies), ".shp"))

    
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

### Let's start with WORLDCLIM (we are going to use CHELSA)


# 4.1. Topographic heterogeneity
  
  # Vector_Ruggedness_Measure
  vrm <- raster("./Data/Topographic_Index/vrm_1KMmn_GMTEDmd.tif") %>%
    crop(box_extent_analysis_bf) 
  names(vrm) <- "Vector_Ruggedness_Measure"
  
  # Roughness
  roughness <- raster("./Data/Topographic_Index/roughness_1KMmn_GMTEDmd.tif") %>%
    crop(box_extent_analysis_bf)   %>%
    resample(vrm, method = "ngb")
  names(roughness) <- "Roughness"
  
  # Slope
  Slope <- raster("./Data/Topographic_Index/slope_1KMmn_GMTEDmd.tif") %>%
    crop(box_extent_analysis_bf)  %>%
    resample(vrm, method = "ngb")
  names(Slope) <- "Slope"
  
  
  # Eastness
  Eastness <- raster("./Data/Topographic_Index/eastness_1KMmn_GMTEDmd.tif") %>%
    crop(box_extent_analysis_bf)  %>%
    resample(vrm, method = "ngb")
  names(Eastness) <- "Eastness"
  
  # Northness
  Northness <- raster("./Data/Topographic_Index/northness_1KMmn_GMTEDmd.tif") %>%
    crop(box_extent_analysis_bf)  %>%
    resample(vrm, method = "ngb")
  names(Northness) <- "Northness"
  
  
# 4.2. Bioclim
  
  # Idenfitying tiles based on Area of study
  
  box_extent_bioclim <- tile_name(box_extent_analysis_bf, "worldclim") # determine which WorldClim tiles your study area intersects with.
  
  # Creating folder
  if(dir.exists("./bioclim_t")){
  }else{
    out_bioclim <-   dir.create("./bioclim_t")}
  
  # Loading tiles
  clim_tiles <- tile_get(box_extent_bioclim, name =  "worldclim", var="bio", path = "./bioclim_t") # for 0.5 arcmin worldclim tiles of 
  
  # Merging tiles
  clim_tiles_merge <- tile_merge(clim_tiles) %>%
      resample(vrm, method = "ngb")
 
# 4.3. Productivity
  
  # Loading Dynamic Habitat Index  (all .tiff files)
    fpar_all <- list.files("./Data/fpar_can", pattern=".tif$", full.names = T)
  
  #  Band 1
      fpar_b1_mean <- lapply(fpar_all, raster, band=1) %>%           # Create raster
        stack()  %>%                                                 # Stack them
        crop(box_extent_analysis_bf)%>%                              # Crop to extent
        calc(mean)%>%                                                # Apply function to calculate mean
        resample(vrm, method = "ngb")                                # resample using previous raster
      names(fpar_b1_mean) <- "Cummulative_annual_productivity_b1"    # rename variable
  
  
  #  Band 2
      fpar_b2_mean <- lapply(fpar_all, raster, band=2) %>%           # Create raster
        stack()  %>%                                                 # Stack them
        crop(box_extent_analysis_bf)%>%                              # Crop to extent
        calc(mean)%>%                                                # Apply function to calculate mean
        resample(vrm, method = "ngb")                                # resample using previous raster
      names(fpar_b2_mean) <- "Minimum_annual_productivity_b2"        # rename variable
          
      
  #  Band 3
      fpar_b3_mean <- lapply(fpar_all, raster, band=3) %>%           # Create raster
        stack()  %>%                                                 # Stack them
        crop(box_extent_analysis_bf)%>%                              # Crop to extent
        calc(mean)%>%                                                # Apply function to calculate mean
        resample(vrm, method = "ngb")                                # resample using previous raster
      names(fpar_b3_mean) <- "Variation_annual_productivity_b3"      # rename variable

# 4.4. Water proximity
      
  # Percentage of lakes
      
    lakes_100m <- raster("./Data/Lakes/Lakes_100_wgs84.tif")%>% # Lakes North America 100m resolution
      crop(box_extent_analysis_bf)
    lakes_100m_agg <- aggregate(lakes_100m, fact=10, fun=sum)%>% # aggregating 100m cells into 1000m cells, using 'sum' function
      resample(vrm, method = "ngb")
    names(lakes_100m_agg) <- "Lakes_percentage"

# 4.5 Landcover data 

# Bare areas
bare_areas_all <- list.files("./Data/Bare_areas_mean", pattern=".tif$", full.names = T)
bare_areas <- lapply(bare_areas_all, raster) %>% 
  stack() %>% 
  crop(box_extent_analysis_bf)%>%                              
  calc(mean)
names(bare_areas) <- "Bare_areas"

# Cropland 
Cropland_all <- list.files("./Data/Cropland_mean", pattern=".tif$", full.names = T)
Cropland <- lapply(Cropland_all, raster) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf) %>% 
  calc(mean)
names(Cropland) <- "Cropland"

# Flooded vegetation 
Flooded_vegetation_all <- list.files("./Data/Flooded_vegetation_mean", pattern=".tif$", full.names = T)
Flooded_vegetation <- lapply(Flooded_vegetation_all, raster) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf) %>% 
  calc(mean)
names(Flooded_vegetation) <- "Flooded_vegetation"

# Grassland 
Grassland_all <- list.files("./Data/Grassland_mean", pattern=".tif$", full.names = T)
Grassland <- lapply(Grassland_all, raster) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf) %>% 
  calc(mean)
names(Grassland) <- "Grassland"

# Other vegetation 
Other_vegetation_all <- list.files("./Data/Other_vegetation_mean", pattern=".tif$", full.names = T)
Other_vegetation <- lapply(Other_vegetation_all, raster, band=1) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf) %>% 
  calc(mean)
names(Other_vegetation) <- "Other_vegetation"

# Snow and ice 
Snow_ice_all <- list.files("./Data/Snow_ice_mean", pattern=".tif$", full.names = T)
Snow_ice <- lapply(Snow_ice_all, raster, band=1) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf) %>% 
  calc(mean)
names(Snow_ice) <- "Snow_and_ice"

# Tree cover 
Tree_cover_all <- list.files("./Data/Tree_cover_mean", pattern=".tif$", full.names = T)
Tree_cover <- lapply(Tree_cover_all, raster, band=1) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf) %>% 
  calc(mean)
names(Tree_cover) <- "Tree_cover"

# Urban areas 
Urban_areas_all <- list.files("./Data/Urban_mean", pattern=".tif$", full.names = T)
Urban_areas <- lapply(Urban_areas_all, raster, band=1) %>%           # Create raster
  stack()  %>%                                                 # Stack them
  crop(box_extent_analysis_bf) %>% 
  calc(mean)
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
        
# 4.7. Removing collinear variables

# Calculating collinearity
    collinearity_test <- removeCollinearity(predictors,
                                             multicollinearity.cutoff = 0.70, 
                                             plot = F, select.variables = T, sample.points = FALSE)

  # Sub-setting variables
    noncollinear_predictors <- stack(subset(predictors, collinearity_test))
    noncollinear_predictors[[2]]
     
# 5. Creating background points ----------
    # Identify number of background points (bg) for HSMs within the geographic range (EBARs), testing three bg.
    
   
    # Calculating number of grid cells in predictor variable
    pixel_1000Km2 <- (1000*1000)/1000000
    study_area_size <- cellStats((((noncollinear_predictors[[1]]*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
    
    # Selecting bg (two options: in case small geographic ranges)
    
    small_study_area <- 0.40 # take 40% if area of study contains less than 100,000 grid cells
    
    huge_study_area <- 0.20 # take 20% if area of study contains  100,000 - 1,000,000 grid cells
    
    very_huge_study_area <- 0.05 # take 20% if area of study contains less than 1,000,000 grid cells

    extremely_huge_study_area <- 0.01 #take 1% if area of study contains more than 2,000,000 grid cells 
    
    
    if(study_area_size < 100000){
     
      background_points <- c(
        round(study_area_size*small_study_area, 0)
      )
    }else if(study_area_size >= 100000 & study_area_size < 1000000){
     
      background_points <- c(
       round(study_area_size*huge_study_area, 0)
      )
    }else if(study_area_size >= 1000000){
      
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

# 7.1. Checking for intersection of predictors stack and observations
    
    points_predictors_overlap <- raster::extract(noncollinear_predictors, points_thin_sf, df=T)%>%
    na.omit()
    obs_thin_noNA <- obs_thin[as.vector(points_predictors_overlap$ID),]
    
    # Creating final occs (occurrences) object
    
    occs <- st_as_sf(x = obs_thin_noNA,                         
                        coords =c("Longitude","Latitude"), crs=wgs84) #Renaming observations columns
      
    occs_df <-  data.frame(st_coordinates(occs))%>%
     dplyr::rename(Longitude = X,  Latitude =  Y)
 
   
# 7.2 creating background points (pseudo-absences) using bias layer
    
    bg_points <- dismo::randomPoints(obs_density_a_wgs84, n = background_points) %>% as.data.frame()
    colnames(bg_points) <- colnames(obs_thin_noNA)

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
    
    if(nrow(occs_df) <=25){
      user_partition  <- get.jackknife(occs_df, bg_bias)
      
    }else if(nrow(occs_df) > 25){
      user_partition  <- get.randomkfold(occs_df, bg_bias, 10)
      
    }
    
    
    
# 8.1.3. Run model, using ENMeval package (Maxent algorithm)
    
  model_species <- ENMeval::ENMevaluate(occs = occs_df,
                                        envs = noncollinear_predictors,
                                        bg = bg_bias, 
                                        algorithm = 'maxent.jar',
                                        method = meth,
                                        partitions = 'user',
                                        user.grp = list(occs.grp = user_partition$occs.grp,
                                                        bg.grp = user_partition$bg.grp),
                                        tune.args = list(fc = features, rm = c(0.05,0.5, 1)), # let' test 3 values
                                        parallel =  TRUE,
                                        doClamp = T,
                                        updateProgress = TRUE,
                                        parallelType = "doParallel"
                                        
  )

 # 8.1.4. Model results 
  
  # All models
  res_model_species <- model_species@results
  
  # let's add background point and observations
    res_model_species$bg_points <- length(model_species@bg[,1])
    res_model_species$occs_points <- length(model_species@occs.grp)
  
  #best_model_species <- model_species@models[model_species@results$or.10p.avg == min(res_model_species$or.10p.avg)]
  
  # Let's now choose the optimal model settings based on the sequential criteria and examine it.
    # lowest value of omission rate (or.10p.avg)
    # to resolve tied or.10p.avg, we used highest AUC value (auc.val.avg)
  
    optimal_model <- eval.results(model_species)%>%
      filter(!is.na(or.10p.avg))%>% 
      filter(or.10p.avg == min(or.10p.avg)) %>% 
      filter(auc.val.avg == max(auc.val.avg))
  
  # Variable importance
    mod.seq <- eval.models(model_species)[[optimal_model$tune.args]]
    variable_importance <- as.data.frame(model_species@variable.importance[[as.vector(optimal_model$tune.args[1])]])[,-2]
    variable_importance
    
# 9. Predictions -------

  all_model_prediction <- eval.predictions(model_species)
  
  best_model_prediction <- eval.predictions(model_species)[[optimal_model$tune.args]]
  plot(best_model_prediction)
  plot(occs, add=T)
  
# 10.  Model uncertainty ---------

# let's run 10 models and calculate the coefficient of variance (the only think that will change is background points)

  # Modelling
    model_10 <- list()
    
    for(i in 1:3){
      cat(paste0("Testing background points_model_", i), '\n')
      
      
      bg_points <- dismo::randomPoints(obs_density_a_wgs84, n = background_points) %>% as.data.frame()
      colnames(bg_points) <- colnames(obs_thin_noNA)
      
      # background using bias
      bg_bias_b <- xyFromCell(!is.na(obs_density_a_wgs84),
                            sample(ncell(!is.na(obs_density_a_wgs84)),
                                   nrow(bg_points),
                                   prob =  values(!is.na(obs_density_a_wgs84))))
      colnames(bg_bias_b) <- colnames(occs_df)
      
      
      # 8.3.2 Partition: select method based on number of occurrences
      
      if(nrow(occs_df) <=25){
        user_partition_b  <- get.jackknife(occs_df, bg_bias_b)
        
      }else if(nrow(occs_df) > 25){
        user_partition_b  <- get.randomkfold(occs_df, bg_bias_b, 10)
        
      }
      
      
      user_partition_b$bg.grp
      
      model_10[[i]] <- ENMevaluate(occs = occs_df,
                                   envs = noncollinear_predictors, bg = bg_bias_b,
                                   partitions='user',
                                   user.grp = list(occs.grp = user_partition_b$occs.grp,
                                                   bg.grp = user_partition_b$bg.grp),
                                   method = meth,
                                   algorithm = 'maxent.jar',
                                   tune.args = list(fc = optimal_model$fc, rm = as.numeric(as.vector(optimal_model$rm))),
                                   parallel =  TRUE,
                                   updateProgress = TRUE,
                                   parallelType = "doParallel"
      )
      
    }

  # Calculating coefficient of variance
    model_10_predictions <- stack(lapply(model_10, eval.predictions))
    uncertainty <- raster::cv(model_10_predictions, na.rm=F)

# 11. Projecting outputs -----------
  # Projecting to aeac to preserve areas and allow calculations 
    model_species_prediction_p <- projectRaster(best_model_prediction, crs = aeac, res = 1000, method = "bilinear")
    uncertainty_p <- projectRaster(uncertainty, crs = aeac, res = 1000, method = "bilinear")
  
    par(mfrow=c(1,2))
    plot(model_species_prediction_p)
    plot(box_extent_analysis_bf_aeac, add=T)
    plot(uncertainty_p)
    plot(box_extent_analysis_bf_aeac, add=T)

plot(best_model_prediction, ylim = c(40, 50))

levelplot(best_model_prediction, margin = FALSE)
levelplot(uncertainty_p, margin = FALSE)

# 12. Write outputs -----------

# 12.1. Load and project into NCC national grid 

nat_grid <- raster("./Data/national_grid/constant_grid.tif")
species_prediction_nat_grid <- projectRaster(best_model_prediction, crs = proj4string(nat_grid), res = 1000, method = "bilinear")
species_uncertainty_nat_grid <- projectRaster(uncertainty, crs = proj4string(nat_grid), res = 1000, method = "bilinear")

# 12.1.1 Write NCC projected rasters 
writeRaster(species_prediction_nat_grid, filename = paste0("./Results/SAR_from_ECCC/", myspecies,".tif"), options = c('TFW = YES'))
writeRaster(species_uncertainty_nat_grid, filename = paste0("./Results/SAR_from_ECCC/",myspecies, "_uncertainty.tif"), options = c('TFW = YES'))

# 12.2. Write aeac projected rasters for calculations 
writeRaster(model_species_prediction_p, filename = paste0("./Results/SAR_from_ECCC_for_calculations/", myspecies,".tif"), options = c('TFW = YES'))
writeRaster(uncertainty_p, filename = paste0("./Results/SAR_from_ECCC_for_calculations/",myspecies, "_uncertainty.tif"), options = c('TFW = YES'))

#############END PIPELINE ####################################
    
#############END PIPELINE ####################################

