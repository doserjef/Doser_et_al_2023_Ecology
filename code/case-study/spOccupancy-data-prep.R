# spOccupancy-data-prep.R: this script prepares the climate and 
#                          LULC covariates for use in spOccupancy
#                          model-fitting functions. It also extracts
#                          the five bioclim variables from the 
#                          downloaded climate data. This script is for
#                          extracting the data at the BBS locations used 
#                          when fitting the candidate models.
# Author: Jeffrey W. Doser
rm(list = ls())
# For calculating bioclim variables. 
library(dismo)

# Read in the BBS formatted data ------------------------------------------
# Brings in coords (spatial coordinates in lat/long), det.covs (detection 
# covariates), sp.codes (species names), and y (detection-nondetection data array)
load("data/bbs-data-formatted.R")

# Read in climate data ----------------------------------------------------
# Load 2017 and 2018 climate data
# Need to extract June 2017 - May 2018 data for calculation of bioclim 
# variables following Rushing et al. (2019) and Clement et al. (2016)
load("data/climate-data/ppt-2017.rda")
ppt.2017 <- ppt.vals
load("data/climate-data/ppt-2018.rda")
ppt.2018 <- ppt.vals
# Final ppt values for time period of interest
ppt.vals <- cbind(ppt.2017[, 6:12], ppt.2018[, 1:5])
load("data/climate-data/tmax-2017.rda")
tmax.2017 <- tmax.vals
load("data/climate-data/tmax-2018.rda")
tmax.2018 <- tmax.vals
# Final tmax values for time period of interest
tmax.vals <- cbind(tmax.2017[, 6:12], tmax.2018[, 1:5])
load("data/climate-data/tmin-2017.rda")
tmin.2017 <- tmin.vals
load("data/climate-data/tmin-2018.rda")
tmin.2018 <- tmin.vals
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
load("data/lulc-covs.rda")
# Don't use crop as this has the highest correlation with another variable
# (forest). 

# Format data for spOccupancy ---------------------------------------------
occ.covs <- data.frame(climate.vars, 
		       lulc.covs[, -which(colnames(lulc.covs) == 'crop')])
data.list <- list(y = y, 
		  det.covs = det.covs, 
		  occ.covs = occ.covs, 
                  coords = as.matrix(coords))
save(data.list, sp.codes, file = 'data/data-bundle.rda')
