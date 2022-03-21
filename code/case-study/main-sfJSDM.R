# main-sfJSDM.R: this script runs a spatially-explicit joint species distribution
#                model using the sfJSDM function in the spOccupancy R package. 
#                This model accounts for species correlations and residual 
#                spatial autocorrelation, but fails to account for imperfect
#                detection. 
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
# Reorder species to help with mixing
# Putting these five species first after exploratory analysis
# REVI, GRSP, PIWO, EAME, BTNW
start.sp <- c('REVI', 'GRSP', 'PIWO', 'EAME', 'BTNW')
# Other species code
indices <- rep(NA, 5)
for (i in 1:5) {
  indices[i] <- which(sp.codes == start.sp[i])
}
indices.other <- 1:nrow(data.list$y)
indices.other <- indices.other[-indices]
# Ordered y
y.ordered <- data.list$y[c(indices, indices.other), , ]
# Update the new data. 
data.list$y <- y.ordered
sp.codes <- sp.codes[c(indices, indices.other)]

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
# Use default initial values for everything except phi and lambda. The spatial range
# parameters and factor loadings are quite sensitive to initial values, and so I fix them to 
# the following values based on a preliminary fit of the model. I then 
# assess convergence of all other parameters in the model using the Gelman-Rubin 
# diagnostic, and assess convergence/adequate mixing of the spatial range
# parameters and spatial factors using traceplots, Gewke Diagnostic, and ESS. 
# Load lambda initial values
n.factors <- 5
N <- nrow(data.list$y)
load("data/inits-sfJSDM.rda")
inits.sfJSDM$z <- NULL
inits.sfJSDM$lambda <- matrix(inits.sfJSDM$lambda, N, n.factors)
save(inits.sfJSDM, file = "data/inits-sfJSDM.rda")
# Run the model -----------------------------------------------------------
n.batch <- 6000
batch.length <- 25
n.neighbors <- 5
n.factors <- 5
cov.model <- "exponential"
n.samples <- n.batch * batch.length
n.burn <- 100000
n.thin <- 50
n.chains <- 1
# TODO: testing
n.batch <- 80
n.burn <- 1000
n.thin <- 1
out <- sfJSDM(formula = ~ scale(elev) + I(scale(elev)^2) + scale(forest) +
		          scale(day) + I(scale(day)^2) + scale(tod) + (1 | obs),
	      data = jsdm.data.list, priors = prior.list, inits = inits.sfJSDM,
	      n.neighbors = n.neighbors, cov.model = cov.model, NNGP = TRUE,
	      n.factors = n.factors, n.batch = n.batch, batch.length = batch.length, 
	      n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	      n.chains = n.chains, n.report = 40, n.omp.threads = 3)
save(out, file = paste("results/bbs-sfJSDM-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
