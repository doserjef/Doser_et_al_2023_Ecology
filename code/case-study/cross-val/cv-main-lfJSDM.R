# cv-main-lfJSDM.R: this script runs a non-spatial joint species distribution model
#                   for the BBS data set using the lfJSDM function from the 
#                   spOccupancy R package. This model accounts for
#                   species correlations, but fails to account for residual
#                   spatial autocorrelation and imperfect detection. This script
#                   only uses 75% of the available data locations for use in 
#                   cross-validation
# Author: Jeffrey W. Doser
# Citation: 

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
# Reformat data for p-ignorant model --
jsdm.data.list <- list()
jsdm.data.list$covs <- data.frame(elev = data.list$occ.covs$elev, 
			          forest = data.list$occ.covs$forest, 
			          day = data.list$det.covs[[1]], 
			          tod = data.list$det.covs[[2]], 
			          obs = data.list$det.covs[[3]]) 
jsdm.data.list$y <- apply(data.list$y, c(1, 2), max, na.rm = TRUE)
jsdm.data.list$coords <- data.list$coords
# The occ.covs and det.covs portion of the data will be ignored in lfJSDM()
# Priors ------------------------------
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1))
# Run the model -----------------------------------------------------------
n.factors <- 5
n.samples <- 50000
n.burn <- 20000
n.thin <- 30
n.chains <- 1
out <- lfJSDM(formula = ~ scale(elev) + I(scale(elev)^2) + scale(forest) +
		          scale(day) + I(scale(day)^2) + scale(tod) + (1 | obs), 
	      data = jsdm.data.list, priors = prior.list, 
	      n.factors = n.factors, n.samples = n.samples, n.burn = n.burn, 
	      n.thin = n.thin, n.chains = n.chains, n.report = 100)
# Save results ------------------------------------------------------------
save(out, file = paste("results/bbs-cv-lfJSDM-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
