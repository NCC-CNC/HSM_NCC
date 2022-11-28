setwd("C:/HSM_NCC")
library(terra)

all <- list.files("./Data/Water_bodies_km", pattern=".tif$", full.names = T)


pclass <- function(x, y=c(0,1)) {
  return( length(which(x %in% y)) / length(x) )
}

pclass(0)


for (i in all){
  
  m <- terra::rast(i)
  
  b <- focal(m, w=matrix(1/9, nc=3, nr=3), fun=pclass)
  
  writeRaster(b, paste("./Data/Water_bodies_prop/", names(m), ".tif"), overwrite = TRUE)} 

