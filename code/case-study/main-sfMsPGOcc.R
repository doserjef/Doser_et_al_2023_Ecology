# main-sfMsPGOcc.R: this script runs a spatially-explicit joint species distribution model
#                   with imperfect detection using sfMsPGOcc from the spOccupancy 
#                   package for the BBS data set. Alternatively, this model can be
#                   viewed as a multispecies occupancy model with species correlations.
#                   This model accounts for (1) imperfect detection;
#                   (2) spatial autocorrelation; and (3) residual species correlations.
# Author: Jeffrey W. Doser

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
# Priors for spatial range parameters
dist.bbs <- dist(data.list$coords)
mean.dist <- mean(dist.bbs)
min.dist <- min(dist.bbs)
max.dist <- max(dist.bbs)
# Priors ------------------------------
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
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
load("data/inits-sfMsPGOcc.rda")
inits.sfMsPGOcc$lambda <- matrix(inits.sfMsPGOcc$lambda, N, n.factors)
save(inits.sfMsPGOcc, file = "data/inits-sfMsPGOcc.rda")
# Run the model -----------------------------------------------------------
n.batch <- 6000
batch.length <- 25
n.neighbors <- 5
cov.model <- "exponential"
n.samples <- n.batch * batch.length
n.burn <- 100000
n.thin <- 50
n.chains <- 1
n.batch <- 400
n.thin <- 5
n.burn <- 5000
out <- sfMsPGOcc(occ.formula = ~ scale(elev) + I(scale(elev)^2) + scale(forest), 
		 det.formula = ~ scale(day) + I(scale(day)^2) + scale(tod) + (1 | obs), 
		 data = data.list, priors = prior.list, inits = inits.sfMsPGOcc,
		 n.neighbors = n.neighbors, cov.model = cov.model, NNGP = TRUE,
		 n.factors = n.factors, n.batch = n.batch, batch.length = batch.length, 
		 n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, n.omp.threads = 3,
		 n.chains = n.chains)
save(out, file = paste("results/bbs-sfMsPGOcc-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
