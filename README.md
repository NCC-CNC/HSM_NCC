# Habitat Suitability Models (HSMs)
## Collaborators
- Juan Zuloaga and Andrew Gonzalez (McGill University)
- Nikol Dimitrov and Richard Schuster (The Nature Conservancy of Canada)

## Objectives
Develop HSMs using citizen science data as a tool to help The Nature Conservancy of Canada in conservation planning.

## Scope
This is an exploratory analysis using coarse resolution datasets (~1km2).  

### Data
- Observations (points): Citizen science data
- Predictors (rasters):
  - Topographic heterogeneity developed by <a href="https://www.nature.com/articles/sdata201840" target="_blank">Amatulli et al 2018</a>. Access date (September 17, 20)  <a href="http://www.earthenv.org/topography" target="_blank">on this website</a>. `Vector Ruggedness Measure (VRM)`, `Roughness`, `Eastness`, `Northness`, `Slope`.
  - The Dynamic Habitat Index (DHI) is an integrated metric of vegetation production from satellite imagery that measures the fraction of photosynthetically active radiation (or fPAR) intercepted by vegetation (<a href="https://www.sciencedirect.com/science/article/pii/S1470160X08000071?casa_token=r7JKpy2f-ocAAAAA:MxkcwYeyPJx-n8_i4efA3gqAWuXOcebBwILc_faNT1oP2otQFxFiF_Zvzcq9As0n0wTBnW2ATA#bib53" target="_blank">Coops et al 2008</a>). `Cummulative annual productivity`, `the minimum level of vegetation cover`, and `the degree of seasonality`.
  - BIOCLIM. Nineteen (19) bioclimatic variables from <a href="https://www.worldclim.org/data/worldclim21.html" target="_blank">Worldclim</a> (<a href="https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/joc.1276" target="_blank">Hijmans et al 2005</a>) were used in this analysis for both 'Historical' and 'Projected' climate. 


### Methods
- We used <a href="https://www.sciencedirect.com/science/article/pii/S030438000500267X" target="_blank">Maxent</a>  algorithm implemented in the `ENMeval` v.2.0.0 package (<a href="https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13628?campaign=woletoc" target="_blank">Kass et al 2021</a>) in R.
- Spatial thining: To avoid issues of spatial sampling and reduce spatial autocorrelation in our HSMs estimates we thinned the occurrence data set of each species using a thinning algorithm (implemented in the <a href="https://cran.r-project.org/web/packages/spThin/index.html" target="_blank">spThin</a> package by <a href="https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.01132" target="_blank">Aiello-Lammens et al 2015</a>), using a distance of 5 Km.
- 
### Outputs (raster 1km2)
- Habitat suitability map 
- Uncertainty map (coeficient of variation 10 runs best model)

### Tasks
- Script basic structure (completed)
- Observations: for every single species include long and lat coordinates (csv file )or shapefile. GBIF dataset???? If so, cleaning module 
- Predictors (Juan will transfer data sets to Nikol)
- Add to script (Juan)
  - Base maps: Canada boundary map (with provinces), Protected Areas,  Cities, Roads
  - Spatial thining
  - Add predictors: Bioclim(19), productivity (3), Topographic heterogeneity (5), distance to lakes(1)
  - Stack predictors
  - Remove collinear predictors
  - Check intersection observation vs predictos (to evaluate observations with NAs)
  - Define background points
  - Create sampling bias layer
  - Model settings
  - Select best model
  - Identify potencial high suitable areas
  - Variable importance, performance metrics
  - Unceratinty map
- Markdown document (Juan generates it and Nikol develop content)
  
 
