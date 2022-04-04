# Joint species distribution models with imperfect detection for high-dimensional spatial data

### In Prep

### Jeffrey W. Doser, Andrew O. Finley, Sudipto Banerjee

### spOccupancy Package [Website](https://www.jeffdoser.com/files/spoccupancy-web/) and [Repository](https://github.com/doserjef/spOccupancy/)

### Please contact the first author for questions about the code or data used in the empirical case studies: Jeffrey W. Doser (doserjef@msu.edu)

---------------------------------

## Abstract

Determining spatial distributions of species and communities are key objectives of ecology and conservation. Joint species distribution models use multi-species detection-nondetection data to estimate species and community distributions. The analysis of such data is complicated by residual correlations between species, imperfect detection, and spatial autocorrelation. While methods exist to accommodate each of these complexities, there are few examples in the literature that address and explore all three complexities simultaneously. Here we developed a spatial factor multi-species occupancy model to explicitly account for species correlations, imperfect detection, and spatial autocorrelation. The proposed model uses a spatial factor dimension reduction approach and Nearest Neighbor Gaussian Processes to ensure computational efficiency for data sets with both a large number of species (e.g., > 100) and spatial locations (e.g., 100,000). We compare the proposed model performance to five candidate models, each addressing a subset of the three complexities. We implemented the proposed and competing models in the `spOccupancy` software, designed to facilitate application via an accessible, well-documented, and open-source R package. Using simulations, we found ignoring the three complexities when present leads to inferior model predictive performance, and the impacts of failing to account for one or more complexities will depend on the objectives of a given study. Using a case study on 98 bird species across the continental US, the spatial factor multi-species occupancy model had the highest predictive performance among the candidate models. Further, our model successfully distinguished between two biogeographical species groups within the 98 species, indicating the potential of our framework as a model-based ordination technique. Our proposed framework, together with its implementation in `spOccupancy`, serves as a user-friendly tool to understand spatial variation in species distributions and biodiversity metrics while addressing common complexities in multi-species detection-nondetection data.    

## Repository Directory

All code and resulting model objects were created and saved using spOccupancy v0.3.0.

### [code/case-study](./code/case-study)

Contains all code used in the BBS case study

+ `bbs-data-prep.R`: preps raw BBS data for analysis with the spatial factor multi-species occupancy model.
+ `bbs-pred-data-prep.R`: prepares the covariate data for predictions across the continental US.
+ `cross-val.R`: contains scripts to perform cross-validation using all five candidate models in the case study.
+ `main-lfJSDM.R`: script to run a non-spatial joint species distribution model with BBS data using `lfJSDM()`. 
+ `main-lfMsPGOcc.R`: script to run a non-spatial latent factor multi-species occupancy model with BBS data using `lfMsPGOcc()`.
+ `main-msPGOcc.R`: script to run a non-spatial multi-species occupancy model with BBS data using `msPGOcc()`. 
+ `main-sfJSDM.R`: script to run a spatial joint species distribution model with BBS data using `sfJSDM()`.
+ `main-sfMsPGOcc.R`: script to run a spatial factor multi-species occupancy model with BBS data using `sfMsPGOcc()`. 
+ `predict-extract.R`: summarizes the full posterior predictive distributions from the prediction results into means and standard deviations to produce smaller model objects for plotting.
+ `predict-lfMsPGOcc.R`: script to predict species richness for the two bird communities across the continental US using `lfMsPGOcc()`. 
+ `predict-sfJSDM.R`: script to predict species richness for the two bird communities across the continental US using `sfJSDM()`.
+ `predict-sfMsPGOcc.R`: script to predict species richness for the two bird communities across the continental US using `sfMsPGOcc()`. 
+ `summary.R`: script to summarize results from all model fits and generate all figures in the manuscript. 

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
+ `bird-species-table-bateman.csv`: CSV file from Bateman et al. (2020) containing information on species classifications that was used to select the two communities of birds.
+ `data-bundle.R`: the data object in the format necessary for fitting models with `spOccupancy`. 
+ `full-bbs-pred-dat.rda`: R data object containing the covariate prediction data to predict species richness across the continental US.
+ `inits-lfJSDM.rda`: R data object containing initial values used to fit `lfJSDM()`. 
+ `inits-lfMsPGOcc.rda`: R data object containing initial values used to fit `lfMsPGOcc()`.
+ `inits-msPGOcc.rda`: R data object containing initial values used to fit `msPGOcc()`. 
+ `inits-sfJSDM.rda`: R data object containing initial values used to fit `sfJSDM()`.
+ `inits-sfMsPGOcc.rda`: R data object containing initial values used to fit `sfMsPGOcc()`.

### [results](.results/)

Directory containing results from `spOccupancy` model fits. Many resulting objects are too large to include on Github, so please contact the first author (doserjef@msu.edu) if these files are desired without having to run the models yourself using the scripts in this repository.

+ `bbs-pred-lfMsPGOcc-summary.R`: R data object consisting of posterior summaries of prediction results from `lfMsPGOcc()`. 
+ `bbs-pred-sfJSDM-summary.R`: R data object consisting of posterior summaries of prediction results from `sfJSDM()`.
+ `bbs-pred-sfMsPGOcc-summary.R`: R data object consisting of posterior summaries of prediction results from `sfMsPGOcc()`.
+ `out-of-sample-deviance.R`: table of out-of-sample deviance metrics from the cross-validation used to assess model predictive performance.
+ `sim-lfJSDM-2022-03-30.R`: simulation results when generating data from Simulation Scenario 1.
+ `sim-sfJSDM-2022-03-30.R`: simulation results when generating data from Simulation Scenario 2.
+ `sim-msPGOcc-2022-03-30.R`: simulation results when generating data from Simulation Scenario 3.
+ `sim-spMsPGOcc-2022-03-30.R`: simulation results when generating data from Simulation Scenario 4.
+ `sim-lfMsPGOcc-2022-03-30.R`: simulation results when generating data from Simulation Scenario 5.
+ `sim-sfMsPGOcc-2022-03-30.R`: simulation results when generating data from Simulation Scenario 6.




