# cv-main-lfJSDM.R: this script runs a non-spatial joint species distribution model
#                   for the BBS data set using the lfJSDM function from the 
#                   spOccupancy R package. This model accounts for
#                   species correlations, but fails to account for residual
#                   spatial autocorrelation and imperfect detection. This script
#                   only uses 75% of the available data locations for use in 
#                   cross-validation.
# Author: Jeffrey W. Doser
# Citation: 

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
load("data/data-bundle.rda")
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
# Update species codes
sp.codes <- sp.codes[c(indices, indices.other)]

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
jsdm.data.list$covs <- data.frame(data.list$occ.covs,
			          day = data.list$det.covs[[1]], 
			          tod = data.list$det.covs[[2]], 
			          obs = data.list$det.covs[[3]]) 
jsdm.data.list$y <- apply(data.list$y, c(1, 2), max, na.rm = TRUE)
jsdm.data.list$coords <- data.list$coords
# The occ.covs and det.covs portion of the data will be ignored in lfJSDM()
# Priors ------------------------------
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
		   tau.sq.beta.ig = list(a = 0.1, b = 0.1))
# Load initial values
# load("data/inits-lfJSDM.rda")
n.factors <- 5
lambda.inits <- matrix(0, nrow = nrow(data.list$y), ncol = n.factors)
diag(lambda.inits) <- 1
inits.list <- list(beta = 0, tau.sq.beta = 1, beta.comm = 0,
		   sigma.sq.psi = 4, lambda = lambda.inits)
# Run the model -----------------------------------------------------------
n.samples <- 50000
n.burn <- 20000
n.thin <- 30
n.chains <- 1
out <- lfJSDM(formula = ~ scale(bio1) + scale(bio2) + scale(bio8) + scale(bio12) +
	                       scale(bio18) + scale(water) + scale(barren) + scale(forest) +
			       scale(grass) + scale(shrub) + scale(hay) + scale(wet) +
			       scale(devel) + scale(day) + I(scale(day)^2) + 
			       scale(tod) + (1 | obs),
	      data = jsdm.data.list, priors = prior.list, inits = inits.list, 
	      n.factors = n.factors, n.samples = n.samples, n.burn = n.burn, 
	      n.thin = n.thin, n.chains = n.chains, n.report = 100)
# Save results ------------------------------------------------------------
save(out, file = paste("results/bbs-cv-lfJSDM-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
