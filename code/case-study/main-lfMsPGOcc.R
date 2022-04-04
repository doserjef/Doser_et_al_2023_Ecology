# main-lfMsPGOcc.R: this script runs a non-spatial joint species distribution model
#                   with imperfect detection using lfMsPGOcc from the spOccupancy 
#                   package for the BBS data set. Alternatively, this model can be
#                   viewed as a multispecies occupancy model with species correlations.
#                   This function accounts for (1) imperfect detection and
#                   (2) residual species correlations. 
# Author: Jeffrey W. Doser

rm(list = ls())
library(spOccupancy)
library(coda)

# Get chain number from command line run ----------------------------------
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# Alternatively, if not running the script from the command line:
# chain <- 1
# Or, can use the n.chains function in spOccupancy (for sequential runs of
# chains).
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
# Updated species codes
sp.codes <- sp.codes[c(indices, indices.other)]

# Prep the model ----------------------------------------------------------
# Priors ------------------------------
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Load initial values to help with convergence and mixing.
load("data/inits-lfMsPGOcc.rda")
# Run the model -----------------------------------------------------------
n.factors <- 5
n.samples <- 150000
n.burn <- 100000
n.thin <- 50
n.chains <- 1
out <- lfMsPGOcc(occ.formula = ~ scale(elev) + I(scale(elev)^2) + scale(forest), 
		 det.formula = ~ scale(day) + I(scale(day)^2) + scale(tod) + (1 | obs), 
		 data = data.list, priors = prior.list, inits = inits.lfMsPGOcc,
		 n.factors = n.factors, n.samples = n.samples, n.burn = n.burn, 
		 n.thin = n.thin, n.chains = n.chains, n.report = 100)
# Save results ------------------------------------------------------------
save(out, file = paste("results/bbs-lfMsPGOcc-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
