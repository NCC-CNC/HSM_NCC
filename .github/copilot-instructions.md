# HSM_NCC Copilot Instructions

## Project Overview
This is a Species Distribution Modeling (SDM) pipeline for **Habitat Suitability Models (HSMs)** targeting Canadian species. The project uses **MaxEnt** algorithm via the `ENMeval` R package to predict species distributions based on environmental predictors and species occurrence data from GBIF.

**Key partners**: The Nature Conservancy of Canada, McGill University

## Architecture & Data Flow

### Core Pipeline Stages
1. **Species Data Acquisition** ([gbif_backbone_species_pull.R](gbif_backbone_species_pull.R))
   - Downloads species occurrence data from GBIF via `rgbif` package
   - Stores individual species CSVs in `Data/GBIF_List/` directory
   - Filters: Canada/US only, coordinate uncertainty ≤1000m, human observations

2. **Data Cleaning & Preparation** ([HSM_species_cleanup_functions.R](HSM_species_cleanup_functions.R))
   - Removes incomplete coordinates, flags suspicious records (outliers, capitals, centroids)
   - Thins observations to 1km resolution using `spThin` to reduce spatial bias
   - Minimum threshold: 5+ observations per species to proceed
   - Outputs: Cleaned spatial point layers (sp/sf objects)

3. **Study Area Definition** ([SDM_Pipeline_01_23_2024_vm_ready.R](SDM_Pipeline_01_23_2024_vm_ready.R), lines ~240)
   - Uses species range map shapefile (if available) OR creates Minimum Convex Polygon (MCP)
   - Buffers MCP by 10km; reprojects to WGS84 for predictor download
   - Stored in ECCC_Polygons directory or dynamically generated

4. **Predictor Assembly** (lines ~260-450)
   - **Topographic**: VRM, Roughness, Eastness, Northness, Slope (EarthEnv)
   - **Climate**: 19 BIOCLIM variables (WorldClim 1970-2000, 30-second resolution)
   - **Productivity**: Dynamic Habitat Index (DHI) - 3 bands from fPAR satellite imagery
   - **Landcover**: ESA-CCI aggregated to 8 types, converted to 9x9 moving window proportions
   - **Lakes**: Calculated from HydroLakes v1 within 1000-pixel window
   - All predictors: cropped to study area, masked, resampled to 1km resolution

5. **MaxEnt Modeling** (lines ~600+)
   - Feature selection by observation count: L, LQ, LQH, LQHP
   - Regularization multiplier (RM): 0.05, 0.50, 1.00 based on observation count
   - Partitioning: k-fold cross-validation (n>25) or jackknife (n<25)
   - Model selection: lowest omission rate, ties broken by highest AUC

6. **Uncertainty Quantification**
   - Runs best model 10x with variable background point selection
   - Outputs: coefficient of variation map across 10 runs

7. **Output Conversion**
   - Continuous suitability map (logistic output)
   - Binary map using 10th percentile training presence threshold
   - Projects to NCC national grid (Albers Equal Area Conic, 1km²)

## Critical Patterns & Conventions

### Environment Configuration
- **Java memory**: Set at script start: `options(java.parameters = "-Xmx64g")`
- **Working directory**: Hardcoded to `C:/HSM_NCC` throughout (shared structure across VMs)
- **File structure**: Assumes strict directory organization: `Data/GBIF_List/`, `Data/Topographic_Index/`, `Data/fpar_can/`, `Results/Continuous_Results/`, etc.

### Species-Level Loop Pattern
```r
# Main processing loop iterates over species (typically from Species_List.csv)
# For each species:
# 1. Load cleaned GBIF data
# 2. Check observation count; skip if <5 after thinning
# 3. Load predictors (expensive operation - crops, masks, resamples each)
# 4. Run ENMeval with model selection
# 5. Generate uncertainty map
# 6. Save continuous + binary outputs
# 7. Log errors to "Reported_model_issues.csv"
```

### Predictor Stacking & Multicollinearity
- Uses `virtualspecies::removeCollinearity()` to filter predictors at r=0.7 correlation threshold
- This ensures uncorrelated feature set for each species

### Key Decision Points (Adaptive Logic)
- **Landcover year selection**: Matches most common year of species observations
- **Background points**: 40%, 20%, 5%, or 1% based on pixel count in study area
- **Feature complexity**: Tied to observation count (lines ~150-160 in README)

## Key Dependencies & External Data

### R Packages
- **Core SDM**: `ENMeval`, `dismo`, `ecospat`
- **Spatial**: `raster`, `sf`, `sp`, `adehabitatHR`, `spThin`
- **GBIF**: `rgbif`, `CoordinateCleaner`
- **Predictors**: `WorldClimTiles` (GitHub custom), `geodata`

### External Datasets (Pre-downloaded)
- **EarthEnv topography** (worldenv.org)
- **WorldClim BIOCLIM 30-second** tiles (~2GB)
- **Silvis Lab DHI** (fpar_can, 3-band GeoTIFF)
- **HydroLakes v1** (lakes dataset)
- **ESA-CCI Landcover** (1992 & 2020)
- **ECCC Species Range Maps** (shapefiles, optional)

### Important Note on Scaling
- Pipeline designed for **Google Cloud Compute** (16 VMs parallel processing)
- Currently supports 1000+ species in scaled production environment

## Common Workflows

### Running Single Species
1. Set `myspecies` variable to scientific name (e.g., "Cervus elaphus")
2. Ensure species CSV exists in `Data/2026_01_21_CA_GBIF_List/`
3. Run pipeline sequentially; monitor `Reported_model_issues.csv` for failures

### Debugging Model Failures
- Check `Reported_model_issues.csv` for stage-specific errors
- Common issues: <5 observations after thinning, missing predictors in study area
- Verify predictor raster alignment: all must resample to VRM template raster

### Performance Considerations
- **Bottleneck**: predictor raster I/O and resampling (~minutes per species)
- **Memory**: Set Java to 64GB for machine with sufficient RAM
- **Garbage collection**: `gc()` called strategically to free raster memory

## Coordinate Systems
- **Working projection**: WGS84 (EPSG:4326) for predictor download
- **Output projection**: Albers Equal Area Conic (aeac) for final NCC grid
- Always check `st_transform()` calls when modifying study area definition

## File Naming Conventions
- Species data: Scientific name (e.g., `Cervus_elaphus.csv`)
- Results: `{species}_continuous.tif`, `{species}_binary.tif`, `{species}_uncertainty.tif`
- Temporary files: stored in `temp_HMs_to_REMOVE/` (cleaned up post-run)

## When to Contact Developers
This pipeline is maintained by: Nikol Dimitrov, Richard Schuster (TNC), Juan Zuloaga (McGill). Key contacts:
- **Species data issues**: Check GBIF download and CoordinateCleaner flags
- **Predictor alignment**: Verify all rasters conform to VRM template
- **Model selection logic**: Refer to Kass et al. 2020 and ENMeval documentation
