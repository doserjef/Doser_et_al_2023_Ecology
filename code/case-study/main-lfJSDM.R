# main-lfJSDM.R: this script runs a non-spatial joint species distribution model
#                for the BBS data set using the lfJSDM function from the 
#                spOccupancy R package. This model accounts for
#                species correlations, but fails to account for residual
#                spatial autocorrelation and imperfect detection.
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
# Use default initial values for everything except lambda. The factor loadings are quite sensitive 
# to initial values, and so I fix them to the following values based on a preliminary fit of the model. I then
# assess convergence of all other parameters in the model using the Gelman-Rubin
# diagnostic, and assess convergence/adequate mixing of the latent factors 
# using traceplots, Gewke Diagnostic, and ESS.
# Load lambda initial values
n.factors <- 5
N <- nrow(data.list$y)
# Load initial values
load("data/inits-lfJSDM.rda")
# inits.lfJSDM$z <- NULL
# inits.lfJSDM$lambda <- matrix(inits.lfJSDM$lambda, N, n.factors)
# save(inits.lfJSDM, file = "data/inits-lfJSDM.rda")
# Run the model -----------------------------------------------------------
# n.samples <- 150000
# n.burn <- 100000
# n.thin <- 50
# n.chains <- 1
n.samples <- 2000
n.burn <- 1000
n.thin <- 1
n.chains <- 1
out <- lfJSDM(formula = ~ scale(elev) + I(scale(elev)^2) + scale(forest) +
		          scale(day) + I(scale(day)^2) + scale(tod) + (1 | obs), 
	      data = jsdm.data.list, priors = prior.list, inits = inits.lfJSDM, 
	      n.factors = n.factors, n.samples = n.samples, n.burn = n.burn, 
	      n.thin = n.thin, n.chains = n.chains, n.report = 100)
# Save results ------------------------------------------------------------
save(out, file = paste("results/bbs-lfJSDM-", chain, "-chain-", 
		       Sys.Date(), ".R", sep = ''))
