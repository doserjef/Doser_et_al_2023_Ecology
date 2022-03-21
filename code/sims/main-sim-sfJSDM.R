# main-sim-sfJSDM.R: this script runs a simulation study comparing the 
#                    six candidate models when data are generated from a
#                    multi-species occupancy model with species interactions
#                    and spatial autocorrelation, under a condition 
#                    of constant and high detection when a JSDM with no
#                    imperfect detection may be considered adequate. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(coda)

# Parameters for simulations ----------------------------------------------
# Number of data sets simulated
n.sims <- 100
# Number of models compared
n.models <- 6
# Random seeds for each data set
set.seed(438)
my.seeds <- sample(1:100000, n.sims, replace = FALSE)
# Spatial locations
J.x <- 15
J.y <- 15
# Total number of spatial locations
J <- J.x * J.y
# Number of replicates
n.rep <- rep(3, J)
# Number of species
N <- 10
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, 1.0)
p.occ <- length(beta.mean)
tau.sq.beta <- c(1.5, .5)
# Detection assumed to be constant, high, and consistent across species.
alpha.mean <- c(plogis(1.386294))
tau.sq.alpha <- c(0.0001)
p.det <- length(alpha.mean)
# Random effects (none for simulations)
psi.RE <- list()
p.RE <- list()
# No factor model
factor.model <- TRUE
n.factors <- 3
# Spatial autocorrelation
phi <- runif(n.factors, 3/1, 3/.1)
sp <- TRUE
# Indices of hold out locations
pred.indx <- seq(1, J, by = 5)

# Simulation setup --------------------------------------------------------
psi.true <- array(NA, dim = c(n.sims, N, J))
psi.mean.samples <- array(NA, dim = c(n.sims, n.models, N, J))
psi.low.samples <- array(NA, dim = c(n.sims, n.models, N, J))
psi.high.samples <- array(NA, dim = c(n.sims, n.models, N, J))
beta.true <- array(NA, dim = c(n.sims, N, p.occ))
beta.mean.samples <- array(NA, dim = c(n.sims, n.models, N * p.occ))
beta.low.samples <- array(NA, dim = c(n.sims, n.models, N * p.occ))
beta.high.samples <- array(NA, dim = c(n.sims, n.models, N * p.occ))
run.time.samples <- matrix(NA, n.sims, n.models)
# MCMC Info ---------------------------
n.samples <- 15000
batch.length <- 25
n.batch <- n.samples / batch.length
n.burn <- 10000
n.thin <- 5
n.chains <- 3
accept.rate <- 0.43

# Simulate data -----------------------------------------------------------
for (j in 1:n.sims) {
  print(paste("Currently on iteration ", j, " out of ", n.sims, sep = ''))
  set.seed(my.seeds[j])
  # Draw species-level effects from community means.
  beta <- matrix(NA, nrow = N, ncol = p.occ)
  alpha <- matrix(NA, nrow = N, ncol = p.det)
  for (i in 1:p.occ) {
    beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
  }
  for (i in 1:p.det) {
    alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
  }
  
  dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
  	        psi.RE = psi.RE, p.RE = p.RE, sp = TRUE, factor.model = factor.model, 
                n.factors = 3, phi = phi, cov.model = 'exponential')
  
  X <- dat$X[-pred.indx, ]
  X.full <- dat$X
  X.p <- dat$X.p[-pred.indx, , ]
  y <- dat$y[, -pred.indx, ]
  y.pred <- dat$y[, pred.indx, ]
  coords <- dat$coords[-pred.indx, ]
  coords.full <- dat$coords
  colnames(X) <- c('int', 'occ.cov.1')
  # Save the data true parameter values
  psi.true[j, , ]  <- dat$psi
  beta.true[j, , ] <- beta
  data.list <- list(y = y, 
  		    occ.covs = X, 
                    coords = coords) 

  # msPGOcc --------------------------- 
  print("Running msPGOcc")
  out <- msPGOcc(occ.formula = ~ occ.cov.1, 
  	         det.formula = ~ 1, 
  	         data = data.list, 
  	         n.samples = n.samples,
  	         n.burn = n.burn,
  	         n.thin = n.thin,
  	         n.chains = n.chains, 
                 verbose = FALSE) 
  beta.mean.samples[j, 1, ] <- apply(out$beta.samples, 2, mean)
  beta.low.samples[j, 1, ] <- apply(out$beta.samples, 2, quantile, 0.025)
  beta.high.samples[j, 1, ] <- apply(out$beta.samples, 2, quantile, 0.975)
  out.pred <- predict(out, X.full)
  psi.mean.samples[j, 1, , ] <- apply(out.pred$psi.0.samples, c(2, 3), mean)
  psi.low.samples[j, 1, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.025)
  psi.high.samples[j, 1, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.975)
  run.time.samples[j, 1] <- out$run.time[3]
  # spMsPGOcc -------------------------
  print("Running spMsPGOcc")
  out <- spMsPGOcc(occ.formula = ~ occ.cov.1, 
        	   det.formula = ~ 1, 
        	   data = data.list,
        	   batch.length = batch.length, 
        	   n.batch = n.batch, 
        	   n.burn = n.burn, 
        	   n.thin = n.thin, 
        	   n.chains = n.chains, 
        	   NNGP = TRUE, 
        	   n.neighbors = 5, 
        	   verbose = FALSE)
  beta.mean.samples[j, 2, ] <- apply(out$beta.samples, 2, mean)
  beta.low.samples[j, 2, ] <- apply(out$beta.samples, 2, quantile, 0.025)
  beta.high.samples[j, 2, ] <- apply(out$beta.samples, 2, quantile, 0.975)
  out.pred <- predict(out, X.full, coords.full, verbose = FALSE)
  psi.mean.samples[j, 2, , ] <- apply(out.pred$psi.0.samples, c(2, 3), mean)
  psi.low.samples[j, 2, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.025)
  psi.high.samples[j, 2, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.975)
  run.time.samples[j, 2] <- out$run.time[3]
  # lfMsPGOcc --------------------------- 
  print("Running lfMsPGOcc")
  out <- lfMsPGOcc(occ.formula = ~ occ.cov.1, 
  	         det.formula = ~ 1, 
  	         data = data.list, 
  	         n.samples = n.samples,
        	 n.factors = 3,
  	         n.burn = n.burn,
  	         n.thin = n.thin,
  	         n.chains = n.chains, 
                 verbose = FALSE) 
  beta.mean.samples[j, 3, ] <- apply(out$beta.samples, 2, mean)
  beta.low.samples[j, 3, ] <- apply(out$beta.samples, 2, quantile, 0.025)
  beta.high.samples[j, 3, ] <- apply(out$beta.samples, 2, quantile, 0.975)
  out.pred <- predict(out, X.full, coords.full)
  psi.mean.samples[j, 3, , ] <- apply(out.pred$psi.0.samples, c(2, 3), mean)
  psi.low.samples[j, 3, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.025)
  psi.high.samples[j, 3, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.975)
  run.time.samples[j, 3] <- out$run.time[3]
  # sfMsPGOcc -------------------------
  print("Running sfMsPGOcc")
  out <- sfMsPGOcc(occ.formula = ~ occ.cov.1, 
        	   det.formula = ~ 1, 
        	   data = data.list,
        	   batch.length = batch.length, 
        	   n.batch = n.batch, 
        	   n.burn = n.burn, 
        	   n.thin = n.thin, 
        	   n.chains = n.chains, 
        	   NNGP = TRUE, 
        	   n.neighbors = 5, 
        	   n.factors = 3,
        	   verbose = FALSE)
  beta.mean.samples[j, 4, ] <- apply(out$beta.samples, 2, mean)
  beta.low.samples[j, 4, ] <- apply(out$beta.samples, 2, quantile, 0.025)
  beta.high.samples[j, 4, ] <- apply(out$beta.samples, 2, quantile, 0.975)
  out.pred <- predict(out, X.full, coords.full, verbose = FALSE)
  psi.mean.samples[j, 4, , ] <- apply(out.pred$psi.0.samples, c(2, 3), mean)
  psi.low.samples[j, 4, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.025)
  psi.high.samples[j, 4, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.975)
  run.time.samples[j, 4] <- out$run.time[3]
  # lfJSDM ----------------------------
  print("Running lfJSDM")
  y.jsdm <- apply(y, c(1, 2), max, na.rm = TRUE)
  data.list <- list(y = y.jsdm, 
        	    covs = X,
        	    coords = coords)
  out <- lfJSDM(formula = ~ occ.cov.1, 
        	data = data.list, 
        	n.samples = n.samples, 
        	n.factors = 3, 
        	n.burn = n.burn, 
        	n.thin = n.thin, 
        	n.chains = n.chains, 
        	verbose = FALSE)
  beta.mean.samples[j, 5, ] <- apply(out$beta.samples, 2, mean)
  beta.low.samples[j, 5, ] <- apply(out$beta.samples, 2, quantile, 0.025)
  beta.high.samples[j, 5, ] <- apply(out$beta.samples, 2, quantile, 0.975)
  out.pred <- predict(out, X.full, coords.full)
  psi.mean.samples[j, 5, , ] <- apply(out.pred$psi.0.samples, c(2, 3), mean)
  psi.low.samples[j, 5, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.025)
  psi.high.samples[j, 5, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.975)
  run.time.samples[j, 5] <- out$run.time[3]
  # sfJSDM ----------------------------
  out <- sfJSDM(formula = ~ occ.cov.1, 
        	data = data.list,
        	batch.length = batch.length, 
        	n.batch = n.batch, 
        	n.burn = n.burn, 
        	n.thin = n.thin, 
        	n.chains = n.chains, 
        	NNGP = TRUE, 
        	n.neighbors = 5, 
        	n.factors = 3,
        	verbose = FALSE)
  beta.mean.samples[j, 6, ] <- apply(out$beta.samples, 2, mean)
  beta.low.samples[j, 6, ] <- apply(out$beta.samples, 2, quantile, 0.025)
  beta.high.samples[j, 6, ] <- apply(out$beta.samples, 2, quantile, 0.975)
  out.pred <- predict(out, X.full, coords.full, verbose = FALSE)
  psi.mean.samples[j, 6, , ] <- apply(out.pred$psi.0.samples, c(2, 3), mean)
  psi.low.samples[j, 6, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.025)
  psi.high.samples[j, 6, , ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.975)
  run.time.samples[j, 6] <- out$run.time[3]
} # j

# Save results ------------------------------------------------------------
save(beta.mean.samples, beta.low.samples, beta.high.samples, 
     psi.mean.samples, psi.low.samples, psi.high.samples, 
     run.time.samples, coords.full, beta.true, psi.true, 
     file = paste("results/sim-sfJSDM-", Sys.Date(), ".R", sep = ''))
