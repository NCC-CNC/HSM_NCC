# ***********************************************************************************************************
# utf8 encoding
# File Name:     Species_data_cleanup_functions.R
# Author:        Nikol Dimitrov & Juan Zuloaga
# Notes:         Species occurrences data cleanup 
# ***********************************************************************************************************

# ***********************************************************************************************************
# utf8 encoding
# File Name:     Species_data_cleanup_functions.R
# Author:        Nikol Dimitrov
# Notes:         Functions to clean up species data to be fed into SDM Pipeline
# ***********************************************************************************************************

# load needed packages for species data cleanup 
list.of.packages <- c("readr", "rgbif", "raster", 'dplyr', "ggplot2", "countrycode", "CoordinateCleaner", "sp", "spThin", "rnaturalearthdata")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only =TRUE) 


# Setting main folder as a default
setwd("C:/HSM_NCC")


# 1.5. Setting Projection to preserve area for all files 

# Lon-Lat 
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

### function format - to speed up processes and scale up to many species 

# FUNCTION 1 - Export data from gbif and store data into object 

## Parameters: 
#### speciesName = scientific name of species of interest 
#;US
export_gbif_data <- function(speciesName){
  Obs_data <- read_csv(paste0("C:/HSM_NCC/Data/2026_01_21_CA_GBIF_List/", speciesName, ".csv"))
  if (is.null(Obs_data)) {
    stop("no species data available")
  } else {
    return(Obs_data)
  }
}

# FUNCTION 2 - Clean up species dataset 

## Parameters: 
### data = takes in loaded dataset from GBIF
### flagged.data = set to FALSE by default (meaning only not flagged or "clean" data is returned, 
### if "TRUE" flagged data is returned)

clean_gbif_data <- function(data, flagged.data = FALSE){ 
  # remove occurrences with incomplete coordinate information
  dat <- data %>% 
    filter(!is.na(decimalLatitude)) %>% 
    filter(!is.na(decimalLongitude))
  
  # convert country code from iso2c to iso3c
  dat$countryCode <- countrycode(dat$countryCode, origin = "iso2c", destination = "iso3c")
  
  # flag suspicious data 
  flag_data <- clean_coordinates(x = dat, 
                                 lon = "decimalLongitude", 
                                 lat = "decimalLatitude", 
                                 countries = "countryCode", 
                                 species = "scientificName", 
                                 tests = c("capitals", "centroids", "equal","gbif", 
                                           "institutions", "zeros", "countries", "outliers"),
                                 verbose = TRUE)
  # if parameter for flagged data is true return only flagged data (set to false by default)
  if (flagged.data == TRUE){                 
    # flagged records dataframe 
    data_fl <- dat[!flag_data$.summary,]
    return(data_fl)
  } else {
    # cleaned dataframe 
    data_cl <- dat[flag_data$.summary,]
    
    #remove occurrences with coordinate undertainty greater than 1km
    if ("coordinateUncertaintyInMeters" %in% colnames(data_cl)){
      data_cl <- data_cl %>% 
        filter(coordinateUncertaintyInMeters / 1000 <= 1 | is.na(coordinateUncertaintyInMeters))
    }
    #filter for records that are human observations or occurrences
    if ("basisOfRecord" %in% colnames(data_cl)){
      data_cl <- data_cl %>% 
        filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "OCCURRENCE")) 
    }
    #remove suspicious individual counts (0 counts or really large counts > 99)
    if ("individualCount" %in% colnames(data_cl)){
      data_cl <- data_cl %>% 
        filter(individualCount > 0 | is.na(individualCount))%>%
        filter(individualCount < 99 | is.na(individualCount)) 
    } 
    return(data_cl)
  } 
}

# FUNCTION 3 -  plot the data 

## Parameters: 
#### data = takes in loaded dataset from GBIF 

plot_gbif_data <- function(data) { 
  if (dim(data)[1] == 0) {
    stop("no species data available")
  } else {
    w_b <- borders("world", colour = "gray50", fill = "gray50")
    ggplot() + coord_fixed() + w_b + 
      geom_point(data = data, aes(decimalLongitude, decimalLatitude),
                 colour = "red", size = 0.5) + 
      theme_bw() +
      xlim(min(data$decimalLongitude) - 20, max(data$decimalLongitude) + 20) + 
      ylim(min(data$decimalLatitude) - 20, max(data$decimalLatitude) + 20)
  }
}

# FUNCTION 4 - Spatial thinning alogrithm 

## Parameters: 
#### data = takes in loaded dataset from GBIF 
#### thinning_par = thinning paramter set by default to 1km 

thin_data <- function(data, thinning_par = 1) {
  
  # run algorithm 
  Obs_data_thinned <- thin(data, lat.col = 'decimalLatitude', long.col = 'decimalLongitude', spec.col = 'scientificName', thin.par = thinning_par, 
                           reps = 1,
                           locs.thinned.list.return = TRUE,
                           write.files = FALSE,
                           write.log.file = FALSE)%>%
    data.frame() 
  
  return(Obs_data_thinned)
}

# FUNCTION - Remove NA predictor variable observations 

# FUNCTION 5 - Convert to spatial points 

## Parameters: 
#### data = takes in loaded dataset from GBIF 

convert_spatial_points <- function(data){
  
  Obs_data_sp <- data.frame(cbind(data$Longitude, data$Latitude))%>% 
    sp::SpatialPoints(proj4string = CRS(wgs84))
  
  return(Obs_data_sp)
  
}
