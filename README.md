# [Joint species distribution models with imperfect detection for high-dimensional spatial data](https://arxiv.org/abs/2204.02707)

### [Jeffrey W. Doser](https://www.jeffdoser.com/), [Andrew O. Finley](https://www.finley-lab.com/), [Sudipto Banerjee](http://sudipto.bol.ucla.edu/)

### Ecology 

### Code/Data DOI: [![DOI](https://zenodo.org/badge/472435954.svg)](https://zenodo.org/badge/latestdoi/472435954)

### spOccupancy Package [Website](https://www.jeffdoser.com/files/spoccupancy-web/) and [Repository](https://github.com/doserjef/spOccupancy/)

### Please contact the first author for questions about the code or data used in the manuscript: Jeffrey W. Doser (doserjef@msu.edu)

---------------------------------

## Abstract

Determining spatial distributions of species and communities are key objectives of ecology and conservation. Joint species distribution models use multi-species detection-nondetection data to estimate species and community distributions. The analysis of such data is complicated by residual correlations between species, imperfect detection, and spatial autocorrelation. While methods exist to accommodate each of these complexities, there are few examples in the literature that address and explore all three complexities simultaneously. Here we developed a spatial factor multi-species occupancy model to explicitly account for species correlations, imperfect detection, and spatial autocorrelation. The proposed model uses a spatial factor dimension reduction approach and Nearest Neighbor Gaussian Processes to ensure computational efficiency for data sets with both a large number of species (e.g., > 100) and spatial locations (e.g., 100,000). We compare the proposed model performance to five candidate models, each addressing a subset of the three complexities. We implemented the proposed and competing models in the spOccupancy software, designed to facilitate application via an accessible, well-documented, and open-source R package. Using simulations, we found ignoring the three complexities when present leads to inferior model predictive performance, and the impacts of failing to account for one or more complexities will depend on the objectives of a given study. Using a case study on 98 bird species across the continental US, the spatial factor multi-species occupancy model had the highest predictive performance among the candidate models. Our proposed framework, together with its implementation in spOccupancy, serves as a user-friendly tool to understand spatial variation in species distributions and biodiversity while addressing common complexities in multi-species detection-nondetection data.     

## Repository Directory

All code and resulting model objects were created and saved using spOccupancy v0.5.0.

### [code/case-study](./code/case-study)

Contains all code used in the BBS case study

+ `bbs-data-prep.R`: preps raw BBS data for analysis with the spatial factor multi-species occupancy model.
+ `bbs-pred-data-prep.R`: code to extract the prediction grid for predicting species richness of the two bird communities across the continental US. 
+ `cross-val`: contains scripts to perform cross-validation using all five candidate models in the case study.
+ `lulc-data-prep.R`: extracts land-use and land-cover variables from USGS EROS for use as occupancy predictors in the case study when fitting the model.
+ `lulc-pred-data-prep.R`: extracts land-use and land-cover variables from USGS EROS for use in prediction of species richness across the continental US.
+ `main-lfJSDM.R`: script to run a non-spatial joint species distribution model with BBS data using `lfJSDM()`. 
+ `main-lfMsPGOcc.R`: script to run a non-spatial latent factor multi-species occupancy model with BBS data using `lfMsPGOcc()`.
+ `main-msPGOcc.R`: script to run a non-spatial multi-species occupancy model with BBS data using `msPGOcc()`. 
+ `main-sfJSDM.R`: script to run a spatial joint species distribution model with BBS data using `sfJSDM()`.
+ `main-sfMsPGOcc.R`: script to run a spatial factor multi-species occupancy model with BBS data using `sfMsPGOcc()`. 
+ `ppt-data-prep.R`: script to prepare the precipitation from PRISM for calculation of the bioclim variables that were used as predictors in the occupancy portion of the model
+ `predict-extract.R`: summarizes the full posterior predictive distributions from the prediction results into means and standard deviations to produce smaller model objects for plotting.
+ `predict-lfMsPGOcc.R`: script to predict species richness for the two bird communities across the continental US using `lfMsPGOcc()`. 
+ `predict-sfJSDM.R`: script to predict species richness for the two bird communities across the continental US using `sfJSDM()`.
+ `predict-sfMsPGOcc.R`: script to predict species richness for the two bird communities across the continental US using `sfMsPGOcc()`. 
+ `predict-spOccupancy-data-prep.R`: combines all covariate data into the format necessary for doing prediction in `spOccupancy`. Also calculates bioclim variables from the precipitation and temperature data.
+ `spOccupancy-data-prep.R`: combine all covariate data and detection-nondetection data into the format for fitting models in `spOccupancy`. Also calculated bioclim variables from the precipitation and temperature data.
+ `summary.R`: script to summarize results from all model fits and generate all figures in the manuscript.
+ `tmax-data-prep.R`: script to prepare the maximum temperature data from PRISM for calculation of the bioclim variables that were used as predictors in the occupancy portion of the model.
+ `tmin-data-prep.R`: script to prepare the minimum temperature data from PRISM for calculation of the bioclim variables that were used as predictors in the occupancy portion of the model.

### [code/sims](./code/sims)

Contains all code used in the simulation study. 

+ `main-sim-lfJSDM.R`: runs a simulation study comparing the six candidate models when data are generated from Simulation Scenario 1.
+ `main-sim-sfJSDM.R`: runs a simulation study comparing the six candidate models when data are generated from Simulation Scenario 2.
+ `main-sim-msPGOcc.R`: runs a simulation study comparing the six candidate models when data are generated from Simulation Scenario 3.
+ `main-sim-spMsPGOcc.R`: runs a simulation study comparing the six candidate models when data are generated from Simulation Scenario 4.
+ `main-sim-lfMsPGOcc.R`: runs a simulation study comparing the six candidate models when data are generated from Simulation Scenario 5.
+ `main-sim-sfMsPGocc.R`: runs a simulation study comparing the six candidate models when data are generated from Simulation Scenario 6.
+ `summary.R`: script that summarizes the simulation results and generates all tables reported in the manuscript. 

### [data](./data/)

Contains data used for BBS case study.

+ `BBS`: directory containing the raw BBS data from the USGS. 
+ `bbs-data-formatted.R`: a temporary object that contains the the formatted BBS data and detection covariates for use in fitting the models.
+ `bird-species-table-bateman.csv`: CSV file from Bateman et al. (2020) containing information on species classifications that was used to select the two communities of birds.
+ `climate-data/`: contains R data file objects of the raw precipitation and temperature data extracted from PRISM at both the locations of the BBS route used to fit the model and across the prediction grid. 
+ `data-bundle.rda`: the data object in the format necessary for fitting models with `spOccupancy`. 
+ `lulc-covs.rda`: the land-use and land-cover variables from USGS EROS extracted at the locations of the BBS routes used to fit the model. 
+ `lulc-pred-covs.rda`: the land-use and land-cover variables from USGS EROS extracted at a prediction grid across the US.
+ `pred-coords.rda`: the coordinates of the prediction locations across the continental US.
+ `pred-data-bundle.rda`: the full prediction data object used for predicting occupancy and species richness across the continental US.


### [results](.results/)

Directory containing results from `spOccupancy` model fits. Many resulting objects are too large to include on Github, so please contact the first author (doserjef@msu.edu) if these files are desired without having to run the models yourself using the scripts in this repository.

+ `bbs-pred-lfMsPGOcc-summary.R`: R data object consisting of posterior summaries of prediction results from `lfMsPGOcc()`. 
+ `bbs-pred-sfJSDM-summary.R`: R data object consisting of posterior summaries of prediction results from `sfJSDM()`.
+ `bbs-pred-sfMsPGOcc-summary.R`: R data object consisting of posterior summaries of prediction results from `sfMsPGOcc()`.
+ `out-of-sample-deviance.R`: table of out-of-sample deviance metrics from the cross-validation used to assess model predictive performance.
+ `sim-lfJSDM-2022-10-20.R`: simulation results when generating data from Simulation Scenario 1.
+ `sim-sfJSDM-2022-10-20.R`: simulation results when generating data from Simulation Scenario 2.
+ `sim-msPGOcc-2022-10-21.R`: simulation results when generating data from Simulation Scenario 3.
+ `sim-spMsPGOcc-2022-10-21.R`: simulation results when generating data from Simulation Scenario 4.
+ `sim-lfMsPGOcc-2022-10-21.R`: simulation results when generating data from Simulation Scenario 5.
+ `sim-sfMsPGOcc-2022-10-21.R`: simulation results when generating data from Simulation Scenario 6.
