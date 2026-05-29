# Habitat Suitability Models (HSMs)
# Nikol Dimitrov, Richard Scuster and Juan Zuloaga
# Pipeline - May 2023
# set Java options here before running any code or packages
options(java.parameters = "-Xmx64g") # for 256GB ram machine
gc()

######### START PIPELINE ######################

# 0. Packages required -----------

list.of.packages <- c("usethis", "spatialEco", "devtools", "rJava", "raster", "rasterVis",
                      "ggplot2", "countrycode", "CoordinateCleaner", "rgbif",
                      "dismo", "ENMeval", "dplyr", "ecospat",
                      "adehabitatHR", "sf", "spThin", "WorldClimTiles", "virtualspecies",
                      "rnaturalearthdata", "geodata")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

devtools::install_github("kapitzas/WorldClimTiles") # https://github.com/kapitzas/WorldClimTiles

if (length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only = TRUE)

system.file("java", package = "dismo")

# 1. Settings ----------

# Helper: create directory (and parents) if absent
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# 1.1. Main folder
ensure_dir("C:/HSM_NCC")
setwd("C:/HSM_NCC")

# 1.2. Output directories
ensure_dir("./temp_HMs_to_REMOVE")
ensure_dir("./Results/Variable Importance")
ensure_dir("./Results/Continuous_Results")
ensure_dir("./Results/Binary_Results")
ensure_dir("./Results/Results_for_calculations")
ensure_dir("./Results/Model_stats")
ensure_dir("./Results/Model_issues")  # fixed: was creating Model_stats here instead

# 1.3. Temp directory for raster scratch files
rasterOptions(tmpdir = file.path("./temp_HMs_to_REMOVE"))

# Helper: append a model issue record to the issues log
log_issue <- function(species_name, data_point_count, issue, stage) {
  n <- data.frame(
    species_name     = species_name,
    date_run         = as.character(Sys.Date()),
    data_point_count = as.character(data_point_count),
    issue            = issue,
    stage            = stage,
    stringsAsFactors = FALSE
  )
  write.table(n, file = "C:/Results_RDS/Model_issues/Reported_model_issues.csv",
              append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
}

# 1.4. Projections (constant across all species — defined once outside the loop)

# Lon-Lat
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"

# Albers Equal Area Conic (https://spatialreference.org/ref/esri/102001/proj4/)
aeac <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

source("Species_data_cleanup_functions_new.R")
R.Version()

species_to_run <- c(
  "Lagopus lagopus",    "Lasiurus borealis",  "Lepus americanus",
  "Lontra canadensis",  "Lynx canadensis",    "Mergus merganser",
  "Ondatra zibethicus", "Phoca vitulina",     "Sterna paradisaea",
  "Ursus americanus",   "Salix vestita",      "Salix reticulata",
  "Salix brachycarpa",  "Salix arctica"
)

for (myspecies in species_to_run) {

  start_time <- Sys.time()

  # 2. Occurrences ------------

  species_rangemap   <- read.csv("C:/Results_RDS/updated_species_list.csv")
  myspecies_rangemap <- species_rangemap[species_rangemap$Sci.Name == myspecies, ]

  sp <- export_gbif_data(myspecies)

  if (is.null(sp)) {
    log_issue(myspecies, "NULL", "No species data on GBIF", "GBIF export")
    next
  }

  if (nrow(sp) < 5) {
    log_issue(myspecies, nrow(sp), "Less than 5 observations - unable to run", "GBIF export")
    next
  }

  plot_gbif_data(sp)
  cl_species <- clean_gbif_data(sp)

  if (myspecies_rangemap$TAX_GROUP == "Birds") {
    cl_species <- cl_species %>% filter(month %in% 4:8)
  }

  if (nrow(cl_species) < 5) {
    log_issue(myspecies, nrow(cl_species), "Less than 5 observations - unable to run", "Cleanup of species data")
    next
  }

  plot_gbif_data(cl_species)

  mode_year <- getmode(cl_species$year)
  land_year <- max(1992, min(2022, mode_year))  # clamp to valid landcover year range

  obs_thin <- thin_data(cl_species)
  gc()

  if (nrow(obs_thin) < 5) {
    log_issue(myspecies, nrow(obs_thin), "Less than 5 observations - unable to run", "Thinning of species data")
    next
  }

  Obs_gbif_data_sp <- convert_spatial_points(obs_thin)

  # 3. Area of study -------------

  if (myspecies_rangemap$Range == "Yes") {

    box_extent_analysis <- st_read(paste0("./Data/ECCC_Polygons/", myspecies, ".shp"))
    box_extent_analysis_bf <- box_extent_analysis %>%
      st_as_sf() %>%
      st_transform(crs = wgs84)
    box_extent_analysis_bf_aeac <- st_as_sf(box_extent_analysis_bf) %>% st_transform(aeac)

  } else {

    box_extent_analysis <- mcp(Obs_gbif_data_sp, percent = 100) %>%
      st_as_sf() %>%
      st_transform(crs = aeac)
    box_extent_analysis_bf <- st_buffer(box_extent_analysis, dist = 100000) %>%
      st_transform(crs = wgs84)
    box_extent_analysis_bf_aeac <- st_as_sf(box_extent_analysis_bf) %>% st_transform(aeac)
  }

  # 4. Predictors ------------

  # 4.1. Topographic heterogeneity
  vrm <- raster("./Data/Topographic_Index/vrm_1KMmn_GMTEDmd.tif") %>%
    crop(box_extent_analysis_bf) %>%
    mask(box_extent_analysis_bf)
  names(vrm) <- "Vector_Ruggedness_Measure"

  load_topo <- function(filename) {
    raster(paste0("./Data/Topographic_Index/", filename)) %>%
      crop(box_extent_analysis_bf) %>%
      resample(vrm, method = "ngb") %>%
      mask(box_extent_analysis_bf)
  }

  roughness <- load_topo("roughness_1KMmn_GMTEDmd.tif"); names(roughness) <- "Roughness"
  Slope     <- load_topo("slope_1KMmn_GMTEDmd.tif");     names(Slope)     <- "Slope"
  Eastness  <- load_topo("eastness_1KMmn_GMTEDmd.tif");  names(Eastness)  <- "Eastness"
  Northness <- load_topo("northness_1KMmn_GMTEDmd.tif"); names(Northness) <- "Northness"
  gc()

  # 4.2. Bioclim (WorldClim 2.1 30s, bands 1-19)
  # FIX: original code had names(bioclim_1) <- ... for all 19 layers,
  #      so only bioclim_1 got renamed (to the last assigned value).
  bioclim_layer_names <- c(
    "Annual Mean Temperature",
    "Mean Diurnal Range",
    "Isothermality",
    "Temperature Seasonality",
    "Max Temperature of Warmest Month",
    "Min Temperature of Coldest Month",
    "Temperature Annual Range",           # was mislabeled as "Min Temperature of Coldest Month"
    "Mean Temperature of Wettest Quarter",
    "Mean Temperature of Driest Quarter",
    "Mean Temperature of Warmest Quarter",
    "Mean Temperature of Coldest Quarter",
    "Annual Precipitation",
    "Precipitation of Wettest Month",
    "Precipitation of Driest Month",
    "Precipitation Seasonality",
    "Precipitation of Wettest Quarter",
    "Precipitation of Driest Quarter",
    "Precipitation of Warmest Quarter",
    "Precipitation of Coldest Quarter"
  )

  bioclim_layers <- lapply(seq_along(bioclim_layer_names), function(b) {
    r <- raster(sprintf("C:/HSM_NCC/bioclim_t/climate/wc2.1_30s/wc2.1_30s_bio_%d.tif", b)) %>%
      crop(box_extent_analysis_bf) %>%
      mask(box_extent_analysis_bf) %>%
      resample(vrm, method = "ngb")
    names(r) <- bioclim_layer_names[b]
    r
  })
  gc()

  # 4.3. Productivity (fpar bands 1-3)
  fpar_all <- list.files("./Data/fpar_can", pattern = "\\.tif$", full.names = TRUE)

  load_fpar_band <- function(band_num, band_name) {
    r <- lapply(fpar_all, raster, band = band_num) %>%
      stack() %>%
      crop(box_extent_analysis_bf) %>%
      calc(mean) %>%
      resample(vrm, method = "ngb") %>%
      mask(box_extent_analysis_bf)
    names(r) <- band_name
    r
  }

  fpar_b1_mean <- load_fpar_band(1, "Cummulative_annual_productivity_b1")
  fpar_b2_mean <- load_fpar_band(2, "Minimum_annual_productivity_b2")
  fpar_b3_mean <- load_fpar_band(3, "Variation_annual_productivity_b3")
  gc()

  # 4.4. Water proximity
  lakes_100m <- raster("./Data/Lakes/Lakes_100_wgs84.tif") %>%
    crop(box_extent_analysis_bf)
  lakes_100m_agg <- aggregate(lakes_100m, fact = 10, fun = sum) %>%
    resample(vrm, method = "ngb") %>%
    mask(box_extent_analysis_bf)
  names(lakes_100m_agg) <- "Lakes_percentage"
  gc()

  # 4.5. Landcover
  lc_config <- list(
    list(subdir = "bare_areas_reclass",         suffix = "bare_areas",         layer_name = "Bare_areas"),
    list(subdir = "cropland_reclass",           suffix = "cropland",           layer_name = "Cropland"),
    list(subdir = "flooded_vegetation_reclass", suffix = "flooded_vegetation", layer_name = "Flooded_vegetation"),
    list(subdir = "grassland_reclass",          suffix = "grassland",          layer_name = "Grassland"),
    list(subdir = "other_vegetation_reclass",   suffix = "other_vegetation",   layer_name = "Other_vegetation"),
    list(subdir = "snow_ice_reclass",           suffix = "snow_ice",           layer_name = "Snow_and_ice"),
    list(subdir = "forests_reclass",            suffix = "forests",            layer_name = "Tree_cover"),
    list(subdir = "urban_reclass",              suffix = "urban",              layer_name = "Urban_areas")
  )

  lc_layers <- lapply(lc_config, function(cfg) {
    path <- sprintf("./Data/Landcover/%s/landcover_%d_%s_resample.tif",
                    cfg$subdir, land_year, cfg$suffix)
    r <- raster(path) %>%
      crop(box_extent_analysis_bf) %>%
      mask(box_extent_analysis_bf)
    names(r) <- cfg$layer_name
    r
  })
  gc()

  # 4.6. Stack predictors — drop any layer whose max == 0
  fpar_nonempty <- Filter(function(r) maxValue(r) != 0, list(fpar_b1_mean, fpar_b2_mean, fpar_b3_mean))
  lc_nonempty   <- Filter(function(r) maxValue(r) != 0, lc_layers)

  all_layers <- c(
    list(vrm, roughness, Slope, Eastness, Northness, lakes_100m_agg),
    fpar_nonempty,
    bioclim_layers,
    lc_nonempty
  )
  predictors <- stack(all_layers)

  rm(vrm, roughness, Slope, Eastness, Northness, lakes_100m_agg, lakes_100m,
     fpar_b1_mean, fpar_b2_mean, fpar_b3_mean, bioclim_layers, lc_layers,
     fpar_nonempty, lc_nonempty)
  gc()

  # 4.7. Remove collinear variables
  collinearity_test <- removeCollinearity(predictors,
                                          multicollinearity.cutoff = 0.70,
                                          plot = FALSE, select.variables = TRUE,
                                          sample.points = FALSE)
  gc()

  noncollinear_predictors <- stack(subset(predictors, collinearity_test))
  rm(predictors, box_extent_analysis, box_extent_analysis_bf_aeac)

  # 5. Background points ----------

  pixel_1000Km2   <- (1000 * 1000) / 1000000
  study_area_size <- cellStats(((noncollinear_predictors[[1]] * 0) + 1) * pixel_1000Km2, "sum", na.rm = TRUE)
  gc()

  # FIX: original had an unreachable >= 3000000 branch (caught by >= 2000000 first)
  bg_fraction <- if      (study_area_size >= 3000000) 0.00025
                 else if (study_area_size >= 2000000) 0.001
                 else if (study_area_size >= 1000000) 0.05
                 else if (study_area_size >= 100000)  0.20
                 else                                 0.40

  background_points <- round(study_area_size * bg_fraction, 0)

  # 6. Sampling bias layer -----------

  gc()

  if (nrow(obs_thin) <= 50) {
    x <- raster::extract(noncollinear_predictors, obs_thin)
    z <- na.omit(data.frame(obs_thin, x))
    if (nrow(z) < 5) {
      log_issue(myspecies, nrow(z), "Less than 5 observations - unable to run", "Removal of NA observations")
      next
    }
    obs_thin <- z %>% dplyr::select("Longitude", "Latitude")
  }

  points_thin_sf <- st_as_sf(data.frame(obs_thin), coords = c("Longitude", "Latitude"), crs = wgs84)
  pts_coords     <- st_coordinates(points_thin_sf)  # compute once, reuse below

  ext <- noncollinear_predictors[[1]]@extent
  obs_density <- MASS::kde2d(
    pts_coords[, 1], pts_coords[, 2],
    lims = c(range(ext@xmin, ext@xmax), range(ext@ymin, ext@ymax))
  ) %>% raster()

  obs_density_a_wgs84 <- projectRaster(obs_density, crs = wgs84, res = 0.008333333) %>%
    resample(noncollinear_predictors[[1]])

  gc()

  # 7. Model fitting ------

  occs <- st_as_sf(obs_thin, coords = c("Longitude", "Latitude"), crs = wgs84)
  occs_df <- data.frame(st_coordinates(occs)) %>%
    dplyr::rename(Longitude = X, Latitude = Y)

  bg_points <- dismo::randomPoints(obs_density_a_wgs84, n = background_points) %>% as.data.frame()
  colnames(bg_points) <- colnames(occs_df)

  # Cache the non-NA mask — original computed !is.na(obs_density_a_wgs84) three separate times
  density_mask <- !is.na(obs_density_a_wgs84)
  bg_bias <- xyFromCell(density_mask,
                        sample(ncell(density_mask), nrow(bg_points), prob = values(density_mask)))
  colnames(bg_bias) <- colnames(occs_df)
  gc()

  # 8. Modeling (Maxent via ENMeval) ------

  # 8.1. Feature classes and partition method, scaled to occurrence count
  if (nrow(occs_df) <= 10) {
    meth <- "jackknife"; features <- "L"
  } else if (nrow(occs_df) <= 15) {
    meth <- "jackknife"; features <- c("L", "Q", "LQ")
  } else if (nrow(occs_df) <= 25) {
    meth <- "jackknife"; features <- c("L", "Q", "H", "LQ", "QH")
  } else if (nrow(occs_df) <= 80) {
    meth <- "randomkfold"; features <- c("L", "Q", "H", "LQ", "QH")
  } else {
    meth <- "randomkfold"; features <- c("Q", "H", "LQ", "LQP", "QPT")
  }

  user_partition <- if (nrow(occs_df) <= 25) {
    get.jackknife(occs_df, bg_bias)
  } else {
    get.randomkfold(occs_df, bg_bias, 10)
  }

  gc()
  parallel::detectCores()

  model_species <- ENMeval::ENMevaluate(
    occs       = occs_df,
    envs       = noncollinear_predictors,
    bg         = bg_bias,
    algorithm  = "maxent.jar",
    method     = meth,
    partitions = "user",
    user.grp   = list(occs.grp = user_partition$occs.grp, bg.grp = user_partition$bg.grp),
    tune.args  = list(fc = features, rm = c(0.05, 0.5, 1)),
    parallel   = FALSE,
    doClamp    = TRUE,
    updateProgress = TRUE,
    parallelType = "doParallel"
  )

  # 8.2. Model results
  res_model_species             <- model_species@results
  res_model_species$bg_points   <- length(model_species@bg[, 1])
  res_model_species$occs_points <- length(model_species@occs.grp)

  optimal_model <- eval.results(model_species) %>%
    filter(!is.na(or.10p.avg)) %>%
    filter(or.10p.avg == min(or.10p.avg)) %>%
    filter(auc.val.avg == max(auc.val.avg))

  # 8.3. Variable importance
  mod.seq <- eval.models(model_species)[[optimal_model$tune.args]]
  plot(mod.seq, type = "cloglog")

  variable_importance <- as.data.frame(
    model_species@variable.importance[[as.vector(optimal_model$tune.args[1])]]
  )[, -2]
  v <- variable_importance %>% arrange(desc(permutation.importance))
  write.csv(v, file = paste0("C:/Results_RDS/Variable Importance/", myspecies, "_variable_importance.csv"))

  # 9. Predictions -------

  best_model_prediction <- eval.predictions(model_species)[[optimal_model$tune.args]]
  plot(best_model_prediction)
  plot(occs, add = TRUE)

  # 10. Model uncertainty (3 replicate background samples) ---------

  model_10 <- vector("list", 3)
  for (j in seq_along(model_10)) {  # FIX: was `for(i in 1:3)` which shadowed the outer loop variable
    cat(paste0("Testing background points_model_", j, "\n"))

    bg_pts_j <- dismo::randomPoints(obs_density_a_wgs84, n = background_points) %>% as.data.frame()
    colnames(bg_pts_j) <- colnames(obs_thin)
    bg_bias_j <- xyFromCell(density_mask,
                            sample(ncell(density_mask), nrow(bg_pts_j), prob = values(density_mask)))
    colnames(bg_bias_j) <- colnames(occs_df)

    model_10[[j]] <- ENMevaluate(
      occs       = occs_df,
      envs       = noncollinear_predictors,
      bg         = bg_bias_j,
      algorithm  = "maxent.jar",
      partitions = "block",
      tune.args  = list(fc = "L", rm = 1),
      parallel   = FALSE,
      updateProgress = TRUE,
      parallelType = "doParallel"
    )
  }

  model_10_predictions <- stack(lapply(model_10, eval.predictions))
  uncertainty          <- raster::cv(model_10_predictions, na.rm = FALSE)

  # 11. Project outputs to AEAC -----------

  model_species_prediction_p <- projectRaster(best_model_prediction, crs = aeac, res = 1000, method = "bilinear")
  uncertainty_p              <- projectRaster(uncertainty,            crs = aeac, res = 1000, method = "bilinear")

  # 12. Threshold (10th-percentile omission) and write outputs -----------

  write_rds(occs,                  paste0("C:/Results_RDS/RDS/Occs/", myspecies))
  write_rds(best_model_prediction, paste0("C:/Results_RDS/RDS/Model Predictions/", myspecies))

  pred_vals_presence <- raster::extract(best_model_prediction, occs)
  pred_vals_presence <- pred_vals_presence[!is.na(pred_vals_presence)]
  n10       <- ceiling(length(pred_vals_presence) * 0.1)
  threshold <- sort(pred_vals_presence)[n10]

  binarize <- function(x) ifelse(x > threshold, 1, 0)

  or.10.threshold.imposed           <- calc(best_model_prediction, fun = binarize)
  model_species_prediction_binary_p <- projectRaster(or.10.threshold.imposed, crs = aeac, res = 1000, method = "ngb")

  plot(or.10.threshold.imposed)

  # 12.1. Project to NCC national grid and write
  nat_grid <- raster("./Data/national_grid/constant_grid.tif")
  nat_crs  <- proj4string(nat_grid)

  species_prediction_nat_grid        <- projectRaster(best_model_prediction,   crs = nat_crs, res = 1000, method = "bilinear")
  species_uncertainty_nat_grid       <- projectRaster(uncertainty,              crs = nat_crs, res = 1000, method = "bilinear")
  binary_species_prediction_nat_grid <- projectRaster(or.10.threshold.imposed, crs = nat_crs, res = 1000, method = "ngb")

  writeRaster(species_prediction_nat_grid,
              paste0("C:/Results_RDS/Continuous_Results/", myspecies, "_continuous.tif"),
              options = c("TFW=YES"), overwrite = TRUE)
  writeRaster(species_uncertainty_nat_grid,
              paste0("C:/Results_RDS/Continuous_Results/", myspecies, "_uncertainty.tif"),
              options = c("TFW=YES"), overwrite = TRUE)
  writeRaster(binary_species_prediction_nat_grid,
              paste0("C:/Results_RDS/Binary_Results/", myspecies, ".tif"),
              options = c("TFW=YES"), datatype = "INT4U", overwrite = TRUE)

  # 12.2. AEAC projected rasters for area calculations
  writeRaster(model_species_prediction_p,
              paste0("C:/Results_RDS/Results_for_calculations/", myspecies, "_continuous.tif"),
              options = c("TFW=YES"), overwrite = TRUE)
  writeRaster(uncertainty_p,
              paste0("C:/Results_RDS/Results_for_calculations/", myspecies, "_uncertainty.tif"),
              options = c("TFW=YES"), overwrite = TRUE)
  writeRaster(model_species_prediction_binary_p,
              paste0("C:/Results_RDS/Results_for_calculations/", myspecies, "_binary.tif"),
              options = c("TFW=YES"), datatype = "INT4U", overwrite = TRUE)

  # Report runtime, AUC, omission rate
  end_time <- Sys.time()

  auc_val <- optimal_model$auc.val.avg
  or_val  <- optimal_model$or.10p.avg

  model_result <- dplyr::case_when(
    auc_val >= 0.75 & or_val <= 0.25                 ~ "model informative, no issues",
    auc_val >= 0.75 & or_val > 0.25 & or_val <= 0.5 ~ "AUC normal, OR value between 0.25-0.5 - some concern in model predictions",
    auc_val <  0.75 & or_val > 0.25 & or_val <= 0.5 ~ "AUC below 0.75 - not informative, OR value between 0.25-0.5 - some concern in model predictions",
    auc_val >= 0.75 & or_val > 0.5                  ~ "AUC normal, OR value greater than 0.5 - high concern in model predictions",
    auc_val <  0.75 & or_val > 0.5                  ~ "AUC below 0.75 - not informative, OR value greater than 0.5 - high concern in model predictions",
    auc_val <  0.75 & or_val <= 0.25                 ~ "AUC below 0.75 - not informative, OR normal"
  )

  stats_row <- data.frame(
    species_name          = myspecies,
    common_landcover_year = mode_year,
    landcover_year_used   = land_year,
    date_run              = as.character(Sys.Date()),
    data_point_count      = nrow(obs_thin),
    run_time              = as.numeric(difftime(end_time, start_time, units = "hours")),
    auc_value             = auc_val,
    or_value              = or_val,
    model_result          = model_result,
    stringsAsFactors      = FALSE
  )
  write.table(stats_row, file = "C:/Results_RDS/Model_stats/Reported_model_results.csv",
              append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)

  # Clean up temp raster scratch files
  f <- list.files("./temp_HMs_to_REMOVE", include.dirs = FALSE, full.names = TRUE, recursive = TRUE)
  file.remove(f)
}
