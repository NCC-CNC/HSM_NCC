library(readr)
library(tidyverse)
library(rgbif)
s <- read_csv("C:/Results_RDS/2026_01_21_sp_list.csv")

gbif_taxon_keys <- s %>% 
  name_backbone_checklist()  %>% # match to backbone 
  filter(!matchType == "NONE") %>% # get matched names
  pull(speciesKey) 

Sys.setenv(GBIF_USER = "nikoldim_97",
           GBIF_PWD = "Formosa@57",
           GBIF_EMAIL = "nikol.dimitrov@natureconservancy.ca")

Obs_gbif_data <- occ_download(pred_in("speciesKey", gbif_taxon_keys),
                              pred("hasCoordinate", TRUE),
                              pred("hasGeospatialIssue", FALSE),
                              pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "OCCURRENCE")),
                              pred_lte("coordinateUncertaintyInMeters", 1000),
                              pred_in("country", c("CA", "US")),
                              format = "SIMPLE_CSV")

Obs_gbif_data

occ_download_wait('0002713-260120142942310')

d <- occ_download_get('0002713-260120142942310') %>%
  occ_download_import()


l <- d %>%
  split(f = as.factor(.$species))

lapply(names(l), function(x) {write.csv(l[[x]], file = paste0("C:/HSM_NCC/Data/GBIF_List/", x, ".csv"))})