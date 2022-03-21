# cv-main-msPGOcc.R: this script runs a non-spatial multi-species occupancy model
#                    msPGOcc from the spOccupancy package for the BBS data set. 
#                    This model accounts for imperfect detection, but does not 
#                    account for residual species correlations or residual 
#                    spatial autocorrelation. This script uses only 75% of the
#                    data for use in cross-validation. 
# Author: Jeffrey W. Doser

rm(list = ls())
library(spOccupancy)
library(coda)

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

# Prep the model ----------------------------------------------------------
# Priors ------------------------------
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   alpha.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
		   tau.sq.alpha.ig = list(a = 0.1, b = 0.1))
# Run the model -----------------------------------------------------------
n.samples <- 50000
n.burn <- 20000
n.thin <- 30
n.chains <- 1
out <- msPGOcc(occ.formula = ~ scale(elev) + I(scale(elev)^2) + scale(forest), 
	       det.formula = ~ scale(day) + I(scale(day)^2) + scale(tod) + (1 | obs), 
	       data = data.list, priors = prior.list, 
	       n.samples = n.samples, n.burn = n.burn, 
	       n.thin = n.thin, n.chains = n.chains, n.report = 100)
# Save results ------------------------------------------------------------
save(out, file = paste("results/bbs-cv-msPGOcc-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
