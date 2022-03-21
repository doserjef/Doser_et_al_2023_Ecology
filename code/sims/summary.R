# summary.R: this script summarizes the simulation results comparing the 
#            six different models that account for subsets of the three 
#            major complexities when modeling species distributions and 
#            communities: (1) imperfect detection, (2) spatial
#            autocorrelation, and (3) residual species correlations.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(coda)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(sf)

# Read in the results -----------------------------------------------------
# Result files are sorted by the model that generated the simulated the data. 
# msPGOcc ------------------------------
load("results/sim-msPGOcc-2022-02-12.R")
beta.high.msPGOcc <- beta.high.samples
beta.low.msPGOcc <- beta.low.samples
beta.mean.msPGOcc <- beta.mean.samples
beta.true.msPGOcc <- beta.true
psi.mean.msPGOcc <- psi.mean.samples
psi.low.msPGOcc <- psi.low.samples
psi.high.msPGOcc <- psi.high.samples
psi.true.msPGOcc <- psi.true
run.time.msPGOcc <- run.time.samples
# spMsPGOcc ------------------------------
load("results/sim-spMsPGOcc-2022-02-12.R")
beta.high.spMsPGOcc <- beta.high.samples
beta.low.spMsPGOcc <- beta.low.samples
beta.mean.spMsPGOcc <- beta.mean.samples
beta.true.spMsPGOcc <- beta.true
psi.mean.spMsPGOcc <- psi.mean.samples
psi.low.spMsPGOcc <- psi.low.samples
psi.high.spMsPGOcc <- psi.high.samples
psi.true.spMsPGOcc <- psi.true
run.time.spMsPGOcc <- run.time.samples
# lfMsPGOcc ------------------------------
load("results/sim-lfMsPGOcc-2022-02-12.R")
beta.high.lfMsPGOcc <- beta.high.samples
beta.low.lfMsPGOcc <- beta.low.samples
beta.mean.lfMsPGOcc <- beta.mean.samples
beta.true.lfMsPGOcc <- beta.true
psi.mean.lfMsPGOcc <- psi.mean.samples
psi.low.lfMsPGOcc <- psi.low.samples
psi.high.lfMsPGOcc <- psi.high.samples
psi.true.lfMsPGOcc <- psi.true
run.time.lfMsPGOcc <- run.time.samples
# sfMsPGOcc ------------------------------
load("results/sim-sfMsPGOcc-2022-02-12.R")
beta.high.sfMsPGOcc <- beta.high.samples
beta.low.sfMsPGOcc <- beta.low.samples
beta.mean.sfMsPGOcc <- beta.mean.samples
beta.true.sfMsPGOcc <- beta.true
psi.mean.sfMsPGOcc <- psi.mean.samples
psi.low.sfMsPGOcc <- psi.low.samples
psi.high.sfMsPGOcc <- psi.high.samples
psi.true.sfMsPGOcc <- psi.true
run.time.sfMsPGOcc <- run.time.samples
# lfJSDM ------------------------------
load("results/sim-lfJSDM-2022-02-12.R")
beta.high.lfJSDM <- beta.high.samples
beta.low.lfJSDM <- beta.low.samples
beta.mean.lfJSDM <- beta.mean.samples
beta.true.lfJSDM <- beta.true
psi.mean.lfJSDM <- psi.mean.samples
psi.low.lfJSDM <- psi.low.samples
psi.high.lfJSDM <- psi.high.samples
psi.true.lfJSDM <- psi.true
run.time.lfJSDM <- run.time.samples
# sfJSDM ------------------------------
load("results/sim-sfJSDM-2022-02-12.R")
beta.high.sfJSDM <- beta.high.samples
beta.low.sfJSDM <- beta.low.samples
beta.mean.sfJSDM <- beta.mean.samples
beta.true.sfJSDM <- beta.true
psi.mean.sfJSDM <- psi.mean.samples
psi.low.sfJSDM <- psi.low.samples
psi.high.sfJSDM <- psi.high.samples
psi.true.sfJSDM <- psi.true
run.time.sfJSDM <- run.time.samples

# Summarize bias in Predicted Occupancy Probability -----------------------
# Calculate bias in psi for each scenario
n.models <- 6
# msPGOcc
psi.bias.msPGOcc <- array(NA, dim = dim(psi.mean.msPGOcc))
for (i in 1:n.models) {
  psi.bias.msPGOcc[, i, , ] <- (psi.mean.msPGOcc[, i, , ] - psi.true.msPGOcc)^2
}
# spMsPGOcc
psi.bias.spMsPGOcc <- array(NA, dim = dim(psi.mean.spMsPGOcc))
for (i in 1:n.models) {
  psi.bias.spMsPGOcc[, i, , ] <- (psi.mean.spMsPGOcc[, i, , ] - psi.true.spMsPGOcc)^2
}
# lfMsPGOcc
psi.bias.lfMsPGOcc <- array(NA, dim = dim(psi.mean.lfMsPGOcc))
for (i in 1:n.models) {
  psi.bias.lfMsPGOcc[, i, , ] <- (psi.mean.lfMsPGOcc[, i, , ] - psi.true.lfMsPGOcc)^2
}
# sfMsPGOcc
psi.bias.sfMsPGOcc <- array(NA, dim = dim(psi.mean.sfMsPGOcc))
for (i in 1:n.models) {
  psi.bias.sfMsPGOcc[, i, , ] <- (psi.mean.sfMsPGOcc[, i, , ] - psi.true.sfMsPGOcc)^2
}
# lfJSDM
psi.bias.lfJSDM <- array(NA, dim = dim(psi.mean.lfJSDM))
for (i in 1:n.models) {
  psi.bias.lfJSDM[, i, , ] <- (psi.mean.lfJSDM[, i, , ] - psi.true.lfJSDM)^2
}
# sfJSDM
psi.bias.sfJSDM <- array(NA, dim = dim(psi.mean.sfJSDM))
for (i in 1:n.models) {
  psi.bias.sfJSDM[, i, , ] <- (psi.mean.sfJSDM[, i, , ] - psi.true.sfJSDM)^2
}

# Average all RMSEs for prediction locations only
psi.pred.bias.msPGOcc <- apply(psi.bias.msPGOcc[, , , pred.indx], c(2, 4), mean)
psi.pred.bias.spMsPGOcc <- apply(psi.bias.spMsPGOcc[, , , pred.indx], c(2, 4), mean)
psi.pred.bias.lfMsPGOcc <- apply(psi.bias.lfMsPGOcc[, , , pred.indx], c(2, 4), mean)
psi.pred.bias.sfMsPGOcc <- apply(psi.bias.sfMsPGOcc[, , , pred.indx], c(2, 4), mean)
psi.pred.bias.lfJSDM <- apply(psi.bias.lfJSDM[, , , pred.indx], c(2, 4), mean)
psi.pred.bias.sfJSDM <- apply(psi.bias.sfJSDM[, , , pred.indx], c(2, 4), mean)

# Average RMSEs -----------------------
apply(psi.pred.bias.msPGOcc, 1, mean)
apply(psi.pred.bias.spMsPGOcc, 1, mean)
apply(psi.pred.bias.lfMsPGOcc, 1, mean)
apply(psi.pred.bias.sfMsPGOcc, 1, mean)
apply(psi.pred.bias.lfJSDM, 1, mean)
apply(psi.pred.bias.sfJSDM, 1, mean)

# Average all RMSEs across data simulations and species
psi.bias.msPGOcc <- apply(psi.bias.msPGOcc, c(2, 4), mean)
psi.bias.spMsPGOcc <- apply(psi.bias.spMsPGOcc, c(2, 4), mean)
psi.bias.lfMsPGOcc <- apply(psi.bias.lfMsPGOcc, c(2, 4), mean)
psi.bias.sfMsPGOcc <- apply(psi.bias.sfMsPGOcc, c(2, 4), mean)
psi.bias.lfJSDM <- apply(psi.bias.lfJSDM, c(2, 4), mean)
psi.bias.sfJSDM <- apply(psi.bias.sfJSDM, c(2, 4), mean)

# Average RMSEs -----------------------
apply(psi.bias.msPGOcc, 1, mean)
apply(psi.bias.spMsPGOcc, 1, mean)
apply(psi.bias.lfMsPGOcc, 1, mean)
apply(psi.bias.sfMsPGOcc, 1, mean)
apply(psi.bias.lfJSDM, 1, mean)
apply(psi.bias.sfJSDM, 1, mean)


# Combine into a data frame for plotting
# psi.bias.df <- data.frame(bias = c(c(psi.bias.msPGOcc), c(psi.bias.spMsPGOcc), 
# 				   c(psi.bias.lfMsPGOcc), c(psi.bias.sfMsPGOcc), 
# 				   c(psi.bias.lfJSDM), c(psi.bias.sfJSDM)), 
# 			  model.sim = rep(c('msPGOcc', 'spMsPGOcc', 
# 					    'lfMsPGOcc', 'sfMsPGOcc', 
# 					    'lfJSDM', 'sfJSDM'), each = J * n.models), 
# 			  model.fit = rep(c('msPGOcc', 'spMsPGOcc', 
# 					    'lfMsPGOcc', 'sfMsPGOcc', 
# 					    'lfJSDM', 'sfJSDM'), times = J * n.models))
# psi.bias.df$model.sim <- factor(psi.bias.df$model.sim, 
# 				levels = c('msPGOcc', 'spMsPGOcc', 
# 					   'lfMsPGOcc', 'sfMsPGOcc', 
# 					   'lfJSDM', 'sfJSDM'))
# psi.bias.df$model.fit <- factor(psi.bias.df$model.fit, 
# 				levels = c('msPGOcc', 'spMsPGOcc', 
# 					   'lfMsPGOcc', 'sfMsPGOcc', 
# 					   'lfJSDM', 'sfJSDM'))
# # Overall bias (but this likely masks a lot of info)
# psi.bias.df %>%
#   group_by(model.sim, model.fit) %>%
#   summarize(mean.bias = mean(abs(bias))) %>%
#   print(n = nrow(.))
# # Add coordinates to psi.bias.df
# psi.bias.df$x <- rep(rep(coords.full[, 1], each = n.models), times = n.models)
# psi.bias.df$y <- rep(rep(coords.full[, 2], each = n.models), times = n.models)
# 
# psi.bias.df %>%
#   #filter(!(model.sim %in% c('lfJSDM', 'sfJSDM')), !(model.fit %in% c('lfJSDM', 'sfJSDM'))) %>%
#   ggplot(aes(x = x, y = y, fill = bias)) + 
#     geom_raster() + 
#     scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
# 		         na.value = 'gray') + 
#     facet_grid(vars(model.sim), vars(model.fit)) +
#     theme_bw(base_size = 18) +
#     theme(axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(), 
#           axis.ticks.y = element_blank(), 
#           axis.text.y = element_blank()) +
#     labs(x = "Data Simulating Model", y = "Model", fill = "Bias")

# Calculate Coverage Rates ------------------------------------------------
J <- ncol(psi.bias.msPGOcc)
n.models <- 6
n.sims <- dim(beta.mean.msPGOcc)[1]
N <- dim(psi.low.msPGOcc)[3]
# msPGOcc
psi.covered.msPGOcc <- array(NA, dim = dim(psi.low.msPGOcc))
for (i in 1:n.sims) {
  for (k in 1:N) {
    for (j in 1:n.models) {
      psi.covered.msPGOcc[i, j, k, ] <- ifelse((psi.true.msPGOcc[i, k, ] > 
          				   psi.low.msPGOcc[i, j, k, ]) & 
          			          (psi.true.msPGOcc[i, k, ] < 
          			           psi.high.msPGOcc[i, j, k, ]), 
          		                   1, 0)
    } # j
  } # k
} # i
# spMsPGOcc
psi.covered.spMsPGOcc <- array(NA, dim = dim(psi.low.spMsPGOcc))
for (i in 1:n.sims) {
  for (k in 1:N) {
    for (j in 1:n.models) {
      psi.covered.spMsPGOcc[i, j, k, ] <- ifelse((psi.true.spMsPGOcc[i, k, ] > 
          				   psi.low.spMsPGOcc[i, j, k, ]) & 
          			          (psi.true.spMsPGOcc[i, k, ] < 
          			           psi.high.spMsPGOcc[i, j, k, ]), 
          		                   1, 0)
    } # j
  } # k
} # i
# lfMsPGOcc
psi.covered.lfMsPGOcc <- array(NA, dim = dim(psi.low.lfMsPGOcc))
for (i in 1:n.sims) {
  for (k in 1:N) {
    for (j in 1:n.models) {
      psi.covered.lfMsPGOcc[i, j, k, ] <- ifelse((psi.true.lfMsPGOcc[i, k, ] > 
          				   psi.low.lfMsPGOcc[i, j, k, ]) & 
          			          (psi.true.lfMsPGOcc[i, k, ] < 
          			           psi.high.lfMsPGOcc[i, j, k, ]), 
          		                   1, 0)
    } # j
  } # k
} # i
# sfMsPGOcc
psi.covered.sfMsPGOcc <- array(NA, dim = dim(psi.low.sfMsPGOcc))
for (i in 1:n.sims) {
  for (k in 1:N) {
    for (j in 1:n.models) {
      psi.covered.sfMsPGOcc[i, j, k, ] <- ifelse((psi.true.sfMsPGOcc[i, k, ] > 
          				   psi.low.sfMsPGOcc[i, j, k, ]) & 
          			          (psi.true.sfMsPGOcc[i, k, ] < 
          			           psi.high.sfMsPGOcc[i, j, k, ]), 
          		                   1, 0)
    } # j
  } # k
} # i
# lfJSDM
psi.covered.lfJSDM <- array(NA, dim = dim(psi.low.lfJSDM))
for (i in 1:n.sims) {
  for (k in 1:N) {
    for (j in 1:n.models) {
      psi.covered.lfJSDM[i, j, k, ] <- ifelse((psi.true.lfJSDM[i, k, ] > 
          				   psi.low.lfJSDM[i, j, k, ]) & 
          			          (psi.true.lfJSDM[i, k, ] < 
          			           psi.high.lfJSDM[i, j, k, ]), 
          		                   1, 0)
    } # j
  } # k
} # i
# sfJSDM
psi.covered.sfJSDM <- array(NA, dim = dim(psi.low.sfJSDM))
for (i in 1:n.sims) {
  for (k in 1:N) {
    for (j in 1:n.models) {
      psi.covered.sfJSDM[i, j, k, ] <- ifelse((psi.true.sfJSDM[i, k, ] > 
          				   psi.low.sfJSDM[i, j, k, ]) & 
          			          (psi.true.sfJSDM[i, k, ] < 
          			           psi.high.sfJSDM[i, j, k, ]), 
          		                   1, 0)
    } # j
  } # k
} # i
# Coverage for each species within each data sim
psi.cov.by.sp.msPGOcc <- apply(psi.covered.msPGOcc, c(1, 2, 3), function(a) sum(a) / J)
psi.cov.by.sp.spMsPGOcc <- apply(psi.covered.spMsPGOcc, c(1, 2, 3), function(a) sum(a) / J)
psi.cov.by.sp.lfMsPGOcc <- apply(psi.covered.lfMsPGOcc, c(1, 2, 3), function(a) sum(a) / J)
psi.cov.by.sp.sfMsPGOcc <- apply(psi.covered.sfMsPGOcc, c(1, 2, 3), function(a) sum(a) / J)
psi.cov.by.sp.lfJSDM <- apply(psi.covered.lfJSDM, c(1, 2, 3), function(a) sum(a) / J)
psi.cov.by.sp.sfJSDM <- apply(psi.covered.sfJSDM, c(1, 2, 3), function(a) sum(a) / J)
# Average Coverage Rates --------------
apply(psi.cov.by.sp.msPGOcc, 2, mean) * 100
apply(psi.cov.by.sp.spMsPGOcc, 2, mean) * 100
apply(psi.cov.by.sp.lfMsPGOcc, 2, mean) * 100
apply(psi.cov.by.sp.sfMsPGOcc, 2, mean) * 100
apply(psi.cov.by.sp.lfJSDM, 2, mean) * 100
apply(psi.cov.by.sp.sfJSDM, 2, mean) * 100

# Look at Occurrence Coefficients -----------------------------------------
# Bias --------------------------------
p.occ <- dim(beta.mean.msPGOcc)[3]
# msPGOcc
beta.bias.msPGOcc <- array(NA, dim = dim(beta.mean.msPGOcc))
beta.true.small.msPGOcc <- matrix(beta.true.msPGOcc, n.sims, p.occ)
for (i in 1:n.models) {
  beta.bias.msPGOcc[, i, ] <- (beta.mean.msPGOcc[, i, ] - beta.true.small.msPGOcc)^2
}
# spMsPGOcc
beta.bias.spMsPGOcc <- array(NA, dim = dim(beta.mean.spMsPGOcc))
beta.true.small.spMsPGOcc <- matrix(beta.true.spMsPGOcc, n.sims, p.occ)
for (i in 1:n.models) {
  beta.bias.spMsPGOcc[, i, ] <- (beta.mean.spMsPGOcc[, i, ] - beta.true.small.spMsPGOcc)^2
}
# lfMsPGOcc
beta.bias.lfMsPGOcc <- array(NA, dim = dim(beta.mean.lfMsPGOcc))
beta.true.small.lfMsPGOcc <- matrix(beta.true.lfMsPGOcc, n.sims, p.occ)
for (i in 1:n.models) {
  beta.bias.lfMsPGOcc[, i, ] <- (beta.mean.lfMsPGOcc[, i, ] - beta.true.small.lfMsPGOcc)^2
}
# sfMsPGOcc
beta.bias.sfMsPGOcc <- array(NA, dim = dim(beta.mean.sfMsPGOcc))
beta.true.small.sfMsPGOcc <- matrix(beta.true.sfMsPGOcc, n.sims, p.occ)
for (i in 1:n.models) {
  beta.bias.sfMsPGOcc[, i, ] <- (beta.mean.sfMsPGOcc[, i, ] - beta.true.small.sfMsPGOcc)^2
}
# lfJSDM
beta.bias.lfJSDM <- array(NA, dim = dim(beta.mean.lfJSDM))
beta.true.small.lfJSDM <- matrix(beta.true.lfJSDM, n.sims, p.occ)
for (i in 1:n.models) {
  beta.bias.lfJSDM[, i, ] <- (beta.mean.lfJSDM[, i, ] - beta.true.small.lfJSDM)^2
}
# sfJSDM
beta.bias.sfJSDM <- array(NA, dim = dim(beta.mean.sfJSDM))
beta.true.small.sfJSDM <- matrix(beta.true.sfJSDM, n.sims, p.occ)
for (i in 1:n.models) {
  beta.bias.sfJSDM[, i, ] <- (beta.mean.sfJSDM[, i, ] - beta.true.small.sfJSDM)^2
}

# Average all biases across data simulations and species
# Note for future this is hardcoded. 
# Intercept
# Means
beta.0.bias.msPGOcc <- apply(beta.bias.msPGOcc[, , 1:10], 2, mean)
beta.0.bias.spMsPGOcc <- apply(beta.bias.spMsPGOcc[, , 1:10], 2, mean)
beta.0.bias.lfMsPGOcc <- apply(beta.bias.lfMsPGOcc[, , 1:10], 2, mean)
beta.0.bias.sfMsPGOcc <- apply(beta.bias.sfMsPGOcc[, , 1:10], 2, mean)
beta.0.bias.lfJSDM <- apply(beta.bias.lfJSDM[, , 1:10], 2, mean)
beta.0.bias.sfJSDM <- apply(beta.bias.sfJSDM[, , 1:10], 2, mean)
# # Lower
# beta.0.bias.low.msPGOcc <- apply(beta.bias.msPGOcc[, , 1:10], 2, quantile, 0.025)
# beta.0.bias.low.spMsPGOcc <- apply(beta.bias.spMsPGOcc[, , 1:10], 2, quantile, 0.025)
# beta.0.bias.low.lfMsPGOcc <- apply(beta.bias.lfMsPGOcc[, , 1:10], 2, quantile, 0.025)
# beta.0.bias.low.sfMsPGOcc <- apply(beta.bias.sfMsPGOcc[, , 1:10], 2, quantile, 0.025)
# beta.0.bias.low.lfJSDM <- apply(beta.bias.lfJSDM[, , 1:10], 2, quantile, 0.025)
# beta.0.bias.low.sfJSDM <- apply(beta.bias.sfJSDM[, , 1:10], 2, quantile, 0.025)
# # Upper
# beta.0.bias.high.msPGOcc <- apply(beta.bias.msPGOcc[, , 1:10], 2, quantile, 0.975)
# beta.0.bias.high.spMsPGOcc <- apply(beta.bias.spMsPGOcc[, , 1:10], 2, quantile, 0.975)
# beta.0.bias.high.lfMsPGOcc <- apply(beta.bias.lfMsPGOcc[, , 1:10], 2, quantile, 0.975)
# beta.0.bias.high.sfMsPGOcc <- apply(beta.bias.sfMsPGOcc[, , 1:10], 2, quantile, 0.975)
# beta.0.bias.high.lfJSDM <- apply(beta.bias.lfJSDM[, , 1:10], 2, quantile, 0.975)
# beta.0.bias.high.sfJSDM <- apply(beta.bias.sfJSDM[, , 1:10], 2, quantile, 0.975)
# Covariate
# Mean bias
beta.1.bias.msPGOcc <- apply(beta.bias.msPGOcc[, , 11:20], 2, mean)
beta.1.bias.spMsPGOcc <- apply(beta.bias.spMsPGOcc[, , 11:20], 2, mean)
beta.1.bias.lfMsPGOcc <- apply(beta.bias.lfMsPGOcc[, , 11:20], 2, mean)
beta.1.bias.sfMsPGOcc <- apply(beta.bias.sfMsPGOcc[, , 11:20], 2, mean)
beta.1.bias.lfJSDM <- apply(beta.bias.lfJSDM[, , 11:20], 2, mean)
beta.1.bias.sfJSDM <- apply(beta.bias.sfJSDM[, , 11:20], 2, mean)
# Lower
# beta.1.bias.low.msPGOcc <- apply(beta.bias.msPGOcc[, , 11:20], 2, quantile, 0.025)
# beta.1.bias.low.spMsPGOcc <- apply(beta.bias.spMsPGOcc[, , 11:20], 2, quantile, 0.025)
# beta.1.bias.low.lfMsPGOcc <- apply(beta.bias.lfMsPGOcc[, , 11:20], 2, quantile, 0.025)
# beta.1.bias.low.sfMsPGOcc <- apply(beta.bias.sfMsPGOcc[, , 11:20], 2, quantile, 0.025)
# beta.1.bias.low.lfJSDM <- apply(beta.bias.lfJSDM[, , 11:20], 2, quantile, 0.025)
# beta.1.bias.low.sfJSDM <- apply(beta.bias.sfJSDM[, , 11:20], 2, quantile, 0.025)
# # Upper
# beta.1.bias.high.msPGOcc <- apply(beta.bias.msPGOcc[, , 11:20], 2, quantile, 0.975)
# beta.1.bias.high.spMsPGOcc <- apply(beta.bias.spMsPGOcc[, , 11:20], 2, quantile, 0.975)
# beta.1.bias.high.lfMsPGOcc <- apply(beta.bias.lfMsPGOcc[, , 11:20], 2, quantile, 0.975)
# beta.1.bias.high.sfMsPGOcc <- apply(beta.bias.sfMsPGOcc[, , 11:20], 2, quantile, 0.975)
# beta.1.bias.high.lfJSDM <- apply(beta.bias.lfJSDM[, , 11:20], 2, quantile, 0.975)
# beta.1.bias.high.sfJSDM <- apply(beta.bias.sfJSDM[, , 11:20], 2, quantile, 0.975)

# Put everything in a data frame
# beta.0.plot.df <- data.frame(mean = c(beta.0.bias.msPGOcc, beta.0.bias.spMsPGOcc, 
# 				      beta.0.bias.lfMsPGOcc, beta.0.bias.sfMsPGOcc, 
# 				      beta.0.bias.lfJSDM, beta.0.bias.sfJSDM), 
# 			     abs.mean = c(beta.0.abs.bias.msPGOcc, beta.0.abs.bias.spMsPGOcc, 
# 				          beta.0.abs.bias.lfMsPGOcc, beta.0.abs.bias.sfMsPGOcc, 
# 				          beta.0.abs.bias.lfJSDM, beta.0.abs.bias.sfJSDM), 
# 			     low = c(beta.0.bias.low.msPGOcc, beta.0.bias.low.spMsPGOcc, 
# 				     beta.0.bias.low.lfMsPGOcc, beta.0.bias.low.sfMsPGOcc, 
# 				     beta.0.bias.low.lfJSDM, beta.0.bias.low.sfJSDM),
# 			     high = c(beta.0.bias.high.msPGOcc, beta.0.bias.high.spMsPGOcc, 
# 				     beta.0.bias.high.lfMsPGOcc, beta.0.bias.high.sfMsPGOcc, 
# 				     beta.0.bias.high.lfJSDM, beta.0.bias.high.sfJSDM), 
# 			     model.sim = factor(rep(c("msPGOcc", "spMsPGOcc", "lfMsPGOcc", 
# 					     "sfMsPGOcc", "lfJSDM", "sfJSDM"), each = n.models), 
# 						levels= c('msPGOcc', 'spMsPGOcc', 
# 							  'lfMsPGOcc', 'sfMsPGOcc', 
# 							  'lfJSDM', 'sfJSDM')), 
# 			     model.fit = factor(rep(c("msPGOcc", "spMsPGOcc", "lfMsPGOcc", 
# 					     "sfMsPGOcc", "lfJSDM", "sfJSDM"), times = n.models), 
# 						levels= c('msPGOcc', 'spMsPGOcc', 
# 							  'lfMsPGOcc', 'sfMsPGOcc', 
# 							  'lfJSDM', 'sfJSDM'))) 
# # Bias plot
# ggplot(data = beta.0.plot.df) +
#   geom_point(aes(x = model.fit, y = mean), size = 4.2) +
#   geom_segment(aes(x = model.fit, xend = model.fit, y = low,
#  		   yend = high), size = 1.2, lineend = 'round') +
#   theme_bw(base_size = 18) + 
#   facet_wrap(vars(model.sim)) +
#   geom_hline(yintercept = 0) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
#   labs(x = 'Model', 
#        y = 'Bias in Intercept') 
# 
# beta.1.plot.df <- data.frame(mean = c(beta.1.bias.msPGOcc, beta.1.bias.spMsPGOcc, 
# 				      beta.1.bias.lfMsPGOcc, beta.1.bias.sfMsPGOcc, 
# 				      beta.1.bias.lfJSDM, beta.1.bias.sfJSDM), 
# 			     abs.mean = c(beta.1.abs.bias.msPGOcc, beta.1.abs.bias.spMsPGOcc, 
# 				          beta.1.abs.bias.lfMsPGOcc, beta.1.abs.bias.sfMsPGOcc, 
# 				          beta.1.abs.bias.lfJSDM, beta.1.abs.bias.sfJSDM), 
# 			     low = c(beta.1.bias.low.msPGOcc, beta.1.bias.low.spMsPGOcc, 
# 				     beta.1.bias.low.lfMsPGOcc, beta.1.bias.low.sfMsPGOcc, 
# 				     beta.1.bias.low.lfJSDM, beta.1.bias.low.sfJSDM),
# 			     high = c(beta.1.bias.high.msPGOcc, beta.1.bias.high.spMsPGOcc, 
# 				     beta.1.bias.high.lfMsPGOcc, beta.1.bias.high.sfMsPGOcc, 
# 				     beta.1.bias.high.lfJSDM, beta.1.bias.high.sfJSDM), 
# 			     model.sim = factor(rep(c("msPGOcc", "spMsPGOcc", "lfMsPGOcc", 
# 					     "sfMsPGOcc", "lfJSDM", "sfJSDM"), each = n.models), 
# 						levels= c('msPGOcc', 'spMsPGOcc', 
# 							  'lfMsPGOcc', 'sfMsPGOcc', 
# 							  'lfJSDM', 'sfJSDM')), 
# 			     model.fit = factor(rep(c("msPGOcc", "spMsPGOcc", "lfMsPGOcc", 
# 					     "sfMsPGOcc", "lfJSDM", "sfJSDM"), times = n.models), 
# 						levels= c('msPGOcc', 'spMsPGOcc', 
# 							  'lfMsPGOcc', 'sfMsPGOcc', 
# 							  'lfJSDM', 'sfJSDM'))) 
# # Bias
# ggplot(data = beta.1.plot.df) +
#   geom_point(aes(x = model.fit, y = mean), size = 4.2) +
#   geom_segment(aes(x = model.fit, xend = model.fit, y = low,
#  		   yend = high), size = 1.2, lineend = 'round') +
#   theme_bw(base_size = 18) + 
#   facet_wrap(vars(model.sim)) +
#   geom_hline(yintercept = 0) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
#   labs(x = 'Model', 
#        y = 'Bias in Covariate Effect') 

# Average RMSEs across all simulations for intercept 
beta.0.bias.msPGOcc
beta.0.bias.spMsPGOcc
beta.0.bias.lfMsPGOcc
beta.0.bias.sfMsPGOcc
beta.0.bias.lfJSDM
beta.0.bias.sfJSDM
# Average RMSEs across all simulations for covariate effect
beta.1.bias.msPGOcc
beta.1.bias.spMsPGOcc
beta.1.bias.lfMsPGOcc
beta.1.bias.sfMsPGOcc
beta.1.bias.lfJSDM
beta.1.bias.sfJSDM

# Calculate Coverage Rates ------------------------------------------------
# Intercept ---------------------------
# msPGOcc
beta.0.covered.msPGOcc <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.0.covered.msPGOcc[i, j, ] <- ifelse((beta.true.small.msPGOcc[i,1:10] > 
					   beta.low.msPGOcc[i, j, 1:10]) & 
				          (beta.true.small.msPGOcc[i, 1:10] < 
				           beta.high.msPGOcc[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# spMsPGOcc
beta.0.covered.spMsPGOcc <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.0.covered.spMsPGOcc[i, j, ] <- ifelse((beta.true.small.spMsPGOcc[i,1:10] > 
					   beta.low.spMsPGOcc[i, j, 1:10]) & 
				          (beta.true.small.spMsPGOcc[i, 1:10] < 
				           beta.high.spMsPGOcc[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# lfMsPGOcc
beta.0.covered.lfMsPGOcc <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.0.covered.lfMsPGOcc[i, j, ] <- ifelse((beta.true.small.lfMsPGOcc[i,1:10] > 
					   beta.low.lfMsPGOcc[i, j, 1:10]) & 
				          (beta.true.small.lfMsPGOcc[i, 1:10] < 
				           beta.high.lfMsPGOcc[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# sfMsPGOcc
beta.0.covered.sfMsPGOcc <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.0.covered.sfMsPGOcc[i, j, ] <- ifelse((beta.true.small.sfMsPGOcc[i,1:10] > 
					   beta.low.sfMsPGOcc[i, j, 1:10]) & 
				          (beta.true.small.sfMsPGOcc[i, 1:10] < 
				           beta.high.sfMsPGOcc[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# lfJSDM
beta.0.covered.lfJSDM <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.0.covered.lfJSDM[i, j, ] <- ifelse((beta.true.small.lfJSDM[i,1:10] > 
					   beta.low.lfJSDM[i, j, 1:10]) & 
				          (beta.true.small.lfJSDM[i, 1:10] < 
				           beta.high.lfJSDM[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# sfJSDM
beta.0.covered.sfJSDM <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.0.covered.sfJSDM[i, j, ] <- ifelse((beta.true.small.sfJSDM[i,1:10] > 
					   beta.low.sfJSDM[i, j, 1:10]) & 
				          (beta.true.small.sfJSDM[i, 1:10] < 
				           beta.high.sfJSDM[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# Coverage rates for intercept
apply(apply(beta.0.covered.msPGOcc, c(1, 2), function(a) sum(a) / (p.occ / 2)), 2, mean)
apply(apply(beta.0.covered.spMsPGOcc, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean)
apply(apply(beta.0.covered.lfMsPGOcc, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean)
apply(apply(beta.0.covered.sfMsPGOcc, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean)
apply(apply(beta.0.covered.lfJSDM, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean)
apply(apply(beta.0.covered.sfJSDM, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean)
# Covariate Effect --------------------
# msPGOcc
beta.1.covered.msPGOcc <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.1.covered.msPGOcc[i, j, ] <- ifelse((beta.true.small.msPGOcc[i,1:10] > 
					   beta.low.msPGOcc[i, j, 1:10]) & 
				          (beta.true.small.msPGOcc[i, 1:10] < 
				           beta.high.msPGOcc[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# spMsPGOcc
beta.1.covered.spMsPGOcc <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.1.covered.spMsPGOcc[i, j, ] <- ifelse((beta.true.small.spMsPGOcc[i,1:10] > 
					   beta.low.spMsPGOcc[i, j, 1:10]) & 
				          (beta.true.small.spMsPGOcc[i, 1:10] < 
				           beta.high.spMsPGOcc[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# lfMsPGOcc
beta.1.covered.lfMsPGOcc <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.1.covered.lfMsPGOcc[i, j, ] <- ifelse((beta.true.small.lfMsPGOcc[i,1:10] > 
					   beta.low.lfMsPGOcc[i, j, 1:10]) & 
				          (beta.true.small.lfMsPGOcc[i, 1:10] < 
				           beta.high.lfMsPGOcc[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# sfMsPGOcc
beta.1.covered.sfMsPGOcc <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.1.covered.sfMsPGOcc[i, j, ] <- ifelse((beta.true.small.sfMsPGOcc[i,1:10] > 
					   beta.low.sfMsPGOcc[i, j, 1:10]) & 
				          (beta.true.small.sfMsPGOcc[i, 1:10] < 
				           beta.high.sfMsPGOcc[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# lfJSDM
beta.1.covered.lfJSDM <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.1.covered.lfJSDM[i, j, ] <- ifelse((beta.true.small.lfJSDM[i,1:10] > 
					   beta.low.lfJSDM[i, j, 1:10]) & 
				          (beta.true.small.lfJSDM[i, 1:10] < 
				           beta.high.lfJSDM[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i
# sfJSDM
beta.1.covered.sfJSDM <- array(NA, dim = c(n.sims, n.models, p.occ / 2))
for (i in 1:n.sims) {
  for (j in 1:n.models) {
    beta.1.covered.sfJSDM[i, j, ] <- ifelse((beta.true.small.sfJSDM[i,1:10] > 
					   beta.low.sfJSDM[i, j, 1:10]) & 
				          (beta.true.small.sfJSDM[i, 1:10] < 
				           beta.high.sfJSDM[i, j, 1:10]), 
			                   1, 0)
  } # j
} # i

# Coverage rates for occurrence covariate effect
apply(apply(beta.1.covered.msPGOcc, c(1, 2), function(a) sum(a) / (p.occ / 2)), 2, mean) * 100
apply(apply(beta.1.covered.spMsPGOcc, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean) * 100
apply(apply(beta.1.covered.lfMsPGOcc, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean) * 100
apply(apply(beta.1.covered.sfMsPGOcc, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean) * 100
apply(apply(beta.1.covered.lfJSDM, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean) * 100
apply(apply(beta.1.covered.sfJSDM, c(1, 2), function(a) sum(a) / (p.occ/2)), 2, mean) * 100

# Run Time Comparisons ----------------
run.time.all <- rbind(run.time.msPGOcc, run.time.spMsPGOcc,
		      run.time.lfMsPGOcc, run.time.sfMsPGOcc, 
		      run.time.lfJSDM, run.time.sfJSDM)
# Average run time in minutes 
apply(run.time.all, 2, mean) / 60
