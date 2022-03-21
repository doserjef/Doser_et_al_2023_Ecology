# main-msPGOcc.R: this script runs a non-spatial multi-species occupancy model
#                 msPGOcc from the spOccupancy package for the BBS data set. 
#                 This model accounts for imperfect detection, but does not 
#                 account for residual species correlations or residual 
#                 spatial autocorrelation. 
# Author: Jeffrey W. Doser

rm(list = ls())
library(spOccupancy)
library(coda)

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

# Prep the model ----------------------------------------------------------
# Priors ------------------------------
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Initial values ----------------------
load("data/inits-msPGOcc.rda")
# Run the model -----------------------------------------------------------
# n.samples <- 150000
# n.burn <- 100000
# n.thin <- 50
# n.chains <- 1
n.samples <- 10000
n.burn <- 5000
n.thin <- 5
n.chains <- 1
out <- msPGOcc(occ.formula = ~ scale(elev) + I(scale(elev)^2) + scale(forest), 
	       det.formula = ~ scale(day) + I(scale(day)^2) + scale(tod) + (1 | obs), 
	       data = data.list, priors = prior.list, inits = inits.msPGOcc,
	       n.samples = n.samples, n.burn = n.burn, 
	       n.thin = n.thin, n.chains = n.chains, n.report = 100)
# Save results ------------------------------------------------------------
save(out, file = paste("results/bbs-msPGOcc-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
