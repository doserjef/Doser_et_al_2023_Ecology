# cv-main-sfJSDM.R: this script runs a spatially-explicit joint species distribution
#                   model using the sfJSDM function in the spOccupancy R package. 
#                   This model accounts for species correlations and residual 
#                   spatial autocorrelation, but fails to account for imperfect
#                   detection. This script fits the model with 75% of the locations
#                   for use in cross-validation. 
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(spOccupancy)
library(coda)
library(sf)

# Get chain number from command line run ----------------------------------
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(chain) == 0) base::stop('Need to tell spOccupancy the chain number')

# Read in the data --------------------------------------------------------
load("data/data-bundle.R")

# Subset the data for cross-validation ------------------------------------
n.hold.out <- round(nrow(data.list$occ.covs) * 0.25)
# For consistency across different model fits
set.seed(100)
pred.indx <- sample(1:nrow(data.list$occ.covs), n.hold.out, replace = FALSE)
# Subset the data 
data.list$y <- data.list$y[, -pred.indx, ]
data.list$det.covs <- lapply(data.list$det.covs, function(a) a[-pred.indx])
data.list$occ.covs <- data.list$occ.covs[-pred.indx, ]
data.list$coords <- data.list$coords[-pred.indx, ]

# Get coordinates in a projection instead of lat-long ---------------------
coords.sf <- st_as_sf(data.frame(data.list$coords),
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Get coordinates in Albers Equal Area
coords.albers <- st_coordinates(coords.sf.albers)
# Convert coordinates to km in Albers equal area. 
data.list$coords <- coords.albers / 1000

# Prep the model ----------------------------------------------------------
# Reformat data for p-ignorant model --
jsdm.data.list <- list()
jsdm.data.list$covs <- data.frame(elev = data.list$occ.covs$elev, 
			          forest = data.list$occ.covs$forest, 
			          day = data.list$det.covs[[1]], 
			          tod = data.list$det.covs[[2]], 
			          obs = data.list$det.covs[[3]]) 
jsdm.data.list$y <- apply(data.list$y, c(1, 2), max, na.rm = TRUE)
jsdm.data.list$coords <- data.list$coords
# Priors ------------------------------
# Priors for spatial range parameters
dist.bbs <- dist(data.list$coords)
mean.dist <- mean(dist.bbs)
min.dist <- min(dist.bbs)
max.dist <- max(dist.bbs)
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist))
# Use the default initial values, which will grab random values. 
# Run the model -----------------------------------------------------------
n.batch <- 2000
batch.length <- 25
n.neighbors <- 5
n.factors <- 5
cov.model <- "exponential"
n.samples <- n.batch * batch.length
n.burn <- 20000
n.thin <- 30
n.chains <- 1
out <- sfJSDM(formula = ~ scale(elev) + I(scale(elev)^2) + scale(forest) +
		          scale(day) + I(scale(day)^2) + scale(tod) + (1 | obs),
	      data = jsdm.data.list, priors = prior.list, 
	      n.neighbors = n.neighbors, cov.model = cov.model, NNGP = TRUE,
	      n.factors = n.factors, n.batch = n.batch, batch.length = batch.length, 
	      n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	      n.chains = n.chains, n.report = 40, n.omp.threads = 1)

# Save results ------------------------------------------------------------
save(out, file = paste("results/bbs-cv-sfJSDM-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
