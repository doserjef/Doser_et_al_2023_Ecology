# predict-spOccupancy-data-prep.R: this script prepares the climate and 
#                                  LULC covariates for use in spOccupancy
#                                  prediction functions. It also extracts
#                                  the five bioclim variables from the 
#                                  downloaded climate data.
# Author: Jeffrey W. Doser
rm(list = ls())
# For calculating bioclim variables. 
library(dismo)

# Read in prediction coordinates ------------------------------------------
load("data/pred-coords.rda")

# Read in climate data ----------------------------------------------------
# Load 2017 and 2018 climate data
# Need to extract June 2017 - May 2018 data for calculation of bioclim 
# variables following Rushing et al. (2019) and Clement et al. (2016)
load("data/climate-data/ppt-2017-pred.rda")
ppt.2017 <- ppt.pred.vals
load("data/climate-data/ppt-2018-pred.rda")
ppt.2018 <- ppt.pred.vals
# Final ppt values for time period of interest
ppt.vals <- cbind(ppt.2017[, 6:12], ppt.2018[, 1:5])
load("data/climate-data/tmax-2017-pred.rda")
tmax.2017 <- tmax.pred.vals
load("data/climate-data/tmax-2018-pred.rda")
tmax.2018 <- tmax.pred.vals
# Final tmax values for time period of interest
tmax.vals <- cbind(tmax.2017[, 6:12], tmax.2018[, 1:5])
load("data/climate-data/tmin-2017-pred.rda")
tmin.2017 <- tmin.pred.vals
load("data/climate-data/tmin-2018-pred.rda")
tmin.2018 <- tmin.pred.vals
# Final tmin values for time period of interest
tmin.vals <- cbind(tmin.2017[, 6:12], tmin.2018[, 1:5])
# Calculate bioclim variables
bioclim.vars <- biovars(ppt.vals, tmin.vals, tmax.vals)
# Extract 5 bioclim variables that have low correlation and are good for 
# modeling species ranges. 
# bio1: mean annual temperature
# bio2: mean diurnal temperature range
# bio8: mean tempearature of the wettest quarter
# bio12: annual precipitation
# bio18: precipitation of warmest quarter. 
climate.vars <- bioclim.vars[, c('bio1', 'bio2', 'bio8', 'bio12', 'bio18')]

# Load LULC data ----------------------------------------------------------
load("data/lulc-pred-covs.rda")
# Don't use crop as this has the highest correlation with another variable
# (forest). 

# Get all variables in a data frame ---------------------------------------
pred.covs <- data.frame(climate.vars, 
		       lulc.pred.covs[, -which(colnames(lulc.pred.covs) == 'crop')])
# Format data for design matrix for use in predict() ----------------------
X.0 <- matrix(NA, nrow(pred.covs), ncol(pred.covs) + 1)
colnames(X.0) <- c('(Intercept)', colnames(pred.covs))
# Load data from model fitting
load("data/data-bundle.rda")
# Standardize prediction values by those used to fit the model
X.0[, "(Intercept)"] <- 1
X.0[, "bio1"] <- (pred.covs$bio1 - mean(data.list$occ.covs$bio1)) / sd(data.list$occ.covs$bio1)
X.0[, "bio2"] <- (pred.covs$bio2 - mean(data.list$occ.covs$bio2)) / sd(data.list$occ.covs$bio2)
X.0[, "bio8"] <- (pred.covs$bio8 - mean(data.list$occ.covs$bio8)) / sd(data.list$occ.covs$bio8)
X.0[, "bio12"] <- (pred.covs$bio12 - mean(data.list$occ.covs$bio12)) / sd(data.list$occ.covs$bio12)
X.0[, "bio18"] <- (pred.covs$bio18 - mean(data.list$occ.covs$bio18)) / sd(data.list$occ.covs$bio18)
X.0[, "water"] <- (pred.covs$water - mean(data.list$occ.covs$water)) / sd(data.list$occ.covs$water)
X.0[, "barren"] <- (pred.covs$barren - mean(data.list$occ.covs$barren)) / sd(data.list$occ.covs$barren)
X.0[, "forest"] <- (pred.covs$forest - mean(data.list$occ.covs$forest)) / sd(data.list$occ.covs$forest)
X.0[, "grass"] <- (pred.covs$grass - mean(data.list$occ.covs$grass)) / sd(data.list$occ.covs$grass)
X.0[, "shrub"] <- (pred.covs$shrub - mean(data.list$occ.covs$shrub)) / sd(data.list$occ.covs$shrub)
X.0[, "hay"] <- (pred.covs$hay - mean(data.list$occ.covs$hay)) / sd(data.list$occ.covs$hay)
X.0[, "wet"] <- (pred.covs$wet - mean(data.list$occ.covs$wet)) / sd(data.list$occ.covs$wet)
X.0[, "devel"] <- (pred.covs$devel - mean(data.list$occ.covs$devel)) / sd(data.list$occ.covs$devel)

save(pred.coords, pred.covs, X.0, file = 'data/pred-data-bundle.rda')
