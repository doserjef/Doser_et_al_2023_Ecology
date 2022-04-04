# predict-out-of-sample.R: this script uses the model fits from the five 
#                          candidate models using 75% of the data, predicts
#                          occurrence at the remaining 25% of the data, and 
#                          then summarizes the results into a model deviance
#                          metric that is used to compare the performance of 
#                          the models for prediction purposes.
# Author: Jeffrey W. Doser

rm(list = ls())
library(coda)
library(spOccupancy)
library(sf)

# Read in the data --------------------------------------------------------
load("data/data-bundle.R")
# Reorder species according to how the model was fit
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
# Load in the indices of the hold-out sites
load("results/pred-indx.R")

# Prepare the prediction data ---------------------------------------------
y.pred <- apply(data.list$y[, pred.indx, ], c(1, 2), max, na.rm = TRUE)
occ.covs.pred <- data.list$occ.covs[pred.indx, ]
coords.pred <- data.list$coords[pred.indx, ]
det.covs.pred <- data.frame(day = data.list$det.covs$day[pred.indx], 
			    tod = data.list$det.covs$tod[pred.indx], 
			    obs = data.list$det.covs$obs[pred.indx])
# Get predicted z-values from other models
# msPGOcc predictions
load('results/bbs-msPGOcc-2-chain-2022-03-29.R')
z.pred.samples.msPGOcc <- out$z.samples[, , pred.indx]
# lfMsPGOcc predictions
load('results/bbs-lfMsPGOcc-2-chain-2022-03-28.R')
z.pred.samples.lfMsPGOcc <- out$z.samples[, , pred.indx]
# sfMsPGOcc predictions
load('results/bbs-sfMsPGOcc-2-chain-2022-03-29.R')
z.pred.samples.sfMsPGOcc <- out$z.samples[, , pred.indx]

# Convert coordinates to albers equal area
coords.sf <- st_as_sf(data.frame(coords.pred),
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Get coordinates in Albers Equal Area
coords.albers <- st_coordinates(coords.sf.albers)
# Convert coordinates to km in Albers equal area. 
coords.pred <- coords.albers / 1000
# Get covariates scaled to the right values based on the values used to fit the model.
occ.covs.fit <- data.list$occ.covs[-pred.indx, ]
det.covs.fit <- lapply(data.list$det.covs, function(a) a[-pred.indx])
elev.pred <- (occ.covs.pred$elev - mean(occ.covs.fit$elev)) / sd(occ.covs.fit$elev)
forest.pred <- (occ.covs.pred$forest - mean(occ.covs.fit$forest)) / sd(occ.covs.fit$forest)
# Occurrence prediction design matrix
X.0 <- cbind(1, elev.pred, elev.pred^2, forest.pred)
day.pred <- (det.covs.pred$day - mean(det.covs.fit$day)) / sd(det.covs.fit$day)
tod.pred <- (det.covs.pred$tod - mean(det.covs.fit$tod)) / sd(det.covs.fit$tod)
# Detection prediction design matrix.
X.p.0 <- cbind(1, day.pred, day.pred^2, tod.pred, det.covs.pred$obs)
colnames(X.p.0) <- c('int', 'day', 'day2', 'tod', 'obs')
# sfMsPGOcc ---------------------------------------------------------------
load("results/bbs-cv-sfMsPGOcc-1-chain-2022-03-30.R")
# Predict occurrence ------------------
out.pred <- predict(out, X.0, coords.pred, n.omp.threads = 10, 
		    verbose = TRUE, n.report = 10, 
                    type = 'occupancy')
# Predict detection -------------------
out.p.pred <- predict(out, X.p.0, type = 'detection')
p.0.samples <- out.p.pred$p.0.samples

# Compute hold-out value deviance -----
# Recover detection samples 
N <- nrow(y.pred)
like.samples <- array(NA, c(N, nrow(X.p.0)))
like.msPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.lfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.sfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
for (i in 1:N) {
  print(paste("Currently on species ", i, " out of ", N, sep = ''))
  for (j in 1:nrow(X.p.0)) {
    like.samples[i, j] <- mean(dbinom(y.pred[i, j], 1,
      			           p.0.samples[, i, j] * out.pred$z.0.samples[, i, j]))
    like.msPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.msPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.lfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.lfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.sfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.sfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
  }
}
# Set 0 values to NA
like.samples[like.samples == 0] <- NA
like.msPGOcc.samples[like.msPGOcc.samples == 0] <- NA
like.lfMsPGOcc.samples[like.lfMsPGOcc.samples == 0] <- NA
like.sfMsPGOcc.samples[like.sfMsPGOcc.samples == 0] <- NA
# Compute deviance metrics
deviance.sfMsPGOcc <- apply(like.samples, 1, function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.ms.sfMsPGOcc <- apply(like.msPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.lf.sfMsPGOcc <- apply(like.lfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.sf.sfMsPGOcc <- apply(like.sfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))

# lfMsPGOcc ---------------------------------------------------------------
load("results/bbs-cv-lfMsPGOcc-1-chain-2022-03-30.R")
# Predict occurrence ------------------
out.pred <- predict(out, X.0, coords.pred, type = 'occupancy')
# Predict detection -------------------
out.p.pred <- predict(out, X.p.0, type = 'detection')
p.0.samples <- out.p.pred$p.0.samples
# Compute hold-out value deviance -----
# Recover detection samples 
N <- nrow(y.pred)
like.samples <- array(NA, c(N, nrow(X.p.0)))
like.msPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.lfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.sfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
for (i in 1:N) {
  print(paste("Currently on species ", i, " out of ", N, sep = ''))
  for (j in 1:nrow(X.p.0)) {
    like.samples[i, j] <- mean(dbinom(y.pred[i, j], 1,
      			           p.0.samples[, i, j] * out.pred$z.0.samples[, i, j]))
    like.msPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.msPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.lfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.lfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.sfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.sfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
  }
}
like.samples[like.samples == 0] <- NA
like.msPGOcc.samples[like.msPGOcc.samples == 0] <- NA
like.lfMsPGOcc.samples[like.lfMsPGOcc.samples == 0] <- NA
like.sfMsPGOcc.samples[like.sfMsPGOcc.samples == 0] <- NA
# Compute deviance metrics
deviance.lfMsPGOcc <- apply(like.samples, 1, function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.ms.lfMsPGOcc <- apply(like.msPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.lf.lfMsPGOcc <- apply(like.lfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.sf.lfMsPGOcc <- apply(like.sfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))

# msPGOcc ---------------------------------------------------------------
load("results/bbs-cv-msPGOcc-1-chain-2022-03-30.R")
# Predict occurrence ------------------
out.pred <- predict(out, X.0, coords.pred, type = 'occupancy')
# Predict detection -------------------
out.p.pred <- predict(out, X.p.0, type = 'detection')
p.0.samples <- out.p.pred$p.0.samples
# Compute hold-out value deviance -----
# Recover detection samples 
N <- nrow(y.pred)
like.samples <- array(NA, c(N, nrow(X.p.0)))
like.msPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.lfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.sfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
for (i in 1:N) {
  print(paste("Currently on species ", i, " out of ", N, sep = ''))
  for (j in 1:nrow(X.p.0)) {
    like.samples[i, j] <- mean(dbinom(y.pred[i, j], 1,
      			           p.0.samples[, i, j] * out.pred$z.0.samples[, i, j]))
    like.msPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.msPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.lfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.lfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.sfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.sfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
  }
}
like.samples[like.samples == 0] <- NA
like.msPGOcc.samples[like.msPGOcc.samples == 0] <- NA
like.lfMsPGOcc.samples[like.lfMsPGOcc.samples == 0] <- NA
like.sfMsPGOcc.samples[like.sfMsPGOcc.samples == 0] <- NA
# Compute deviance metrics
deviance.msPGOcc <- apply(like.samples, 1, function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.ms.msPGOcc <- apply(like.msPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.lf.msPGOcc <- apply(like.lfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.sf.msPGOcc <- apply(like.sfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))

# lfJSDM ---------------------------------------------------------------
load("results/bbs-cv-lfJSDM-1-chain-2022-03-29.R")
X.jsdm.0 <- cbind(X.0, X.p.0[, -1])
colnames(X.jsdm.0) <- c('int', 'elev', 'elev.2', 'forest', 'day', 'day.2', 'tod', 'obs')
# Predict at the hold-out locations ---------------------------------------
out.pred <- predict(out, X.jsdm.0, coords.pred, ignore.RE = FALSE)
# Compute hold-out value deviance -----------------------------------------
like.samples <- array(NA, c(N, nrow(X.p.0)))
like.msPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.lfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.sfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
for (i in 1:N) {
  for (j in 1:nrow(X.p.0)) {
    like.samples[i, j] <- mean(dbinom(y.pred[i, j], 1,
      			           out.pred$psi.0.samples[, i, j]))
    like.msPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.msPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.lfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.lfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.sfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.sfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
  }
}
like.samples[like.samples == 0] <- NA
like.msPGOcc.samples[like.msPGOcc.samples == 0] <- NA
like.lfMsPGOcc.samples[like.lfMsPGOcc.samples == 0] <- NA
like.sfMsPGOcc.samples[like.sfMsPGOcc.samples == 0] <- NA
# Compute deviance metrics
deviance.lfJSDM <- apply(like.samples, 1, function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.ms.lfJSDM <- apply(like.msPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.lf.lfJSDM <- apply(like.lfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.sf.lfJSDM <- apply(like.sfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))

# sfJSDM ---------------------------------------------------------------
load("results/bbs-cv-sfJSDM-1-chain-2022-03-30.R")
X.jsdm.0 <- cbind(X.0, X.p.0[, -1])
colnames(X.jsdm.0) <- c('int', 'elev', 'elev.2', 'forest', 'day', 'day.2', 'tod', 'obs')
# Predict at the hold-out locations ---------------------------------------
out.pred <- predict(out, X.jsdm.0, coords.pred)
# Compute hold-out value deviance -----------------------------------------
like.samples <- array(NA, c(N, nrow(X.p.0)))
like.msPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.lfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
like.sfMsPGOcc.samples <- array(NA, c(N, nrow(X.p.0)))
for (i in 1:N) {
  for (j in 1:nrow(X.p.0)) {
    like.samples[i, j] <- mean(dbinom(y.pred[i, j], 1,
      			           out.pred$psi.0.samples[, i, j]))
    like.msPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.msPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.lfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.lfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
    like.sfMsPGOcc.samples[i, j] <- mean(dbinom(z.pred.samples.sfMsPGOcc[, i, j], 1,
					      out.pred$psi.0.samples[, i, j]))
  }
}
like.samples[like.samples == 0] <- NA
like.msPGOcc.samples[like.msPGOcc.samples == 0] <- NA
like.lfMsPGOcc.samples[like.lfMsPGOcc.samples == 0] <- NA
like.sfMsPGOcc.samples[like.sfMsPGOcc.samples == 0] <- NA
# Compute deviance metrics
deviance.sfJSDM <- apply(like.samples, 1, function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.ms.sfJSDM <- apply(like.msPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.lf.sfJSDM <- apply(like.lfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))
deviance.z.sf.sfJSDM <- apply(like.sfMsPGOcc.samples, 1, 
				 function(a) -2 * sum(log(a), na.rm = TRUE))

# Store all deviances together --------------------------------------------
n.deviances <- 5 * 4
deviance.df <- data.frame(deviance = c(deviance.msPGOcc, deviance.lfMsPGOcc, 
				       deviance.sfMsPGOcc, deviance.lfJSDM, 
				       deviance.sfJSDM, deviance.z.ms.msPGOcc, 
				       deviance.z.ms.lfMsPGOcc, 
				       deviance.z.ms.sfMsPGOcc, deviance.z.ms.lfJSDM, 
				       deviance.z.ms.sfJSDM, deviance.z.lf.msPGOcc, 
				       deviance.z.lf.lfMsPGOcc, 
				       deviance.z.lf.sfMsPGOcc, deviance.z.lf.lfJSDM, 
				       deviance.z.lf.sfJSDM, deviance.z.sf.msPGOcc, 
				       deviance.z.sf.lfMsPGOcc, 
				       deviance.z.sf.sfMsPGOcc, deviance.z.sf.lfJSDM, 
				       deviance.z.sf.sfJSDM), 
			  sp = rep(sp.codes, times = n.deviances), 
			  model.fit = rep(rep(c('msPGOcc', 'lfMsPGOcc', 
					    'sfMsPGOcc', 'lfJSDM', 'sfJSDM'), 
					  each = N), times = 4), 
                          type.deviance = rep(c('y', 'z.msPGOcc', 
					            'z.lfMsPGOcc', 'z.sfMsPGOcc'), 
				             each = N * 5))

save(deviance.df, file = "results/out-of-sample-deviance.R")
