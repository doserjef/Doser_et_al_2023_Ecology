# predict-sfMsPGOcc.R: this script predicts species richness for the two 
#                      communities (Eastern Forest Birds and Grasslands)
#                      across the continental US. Predictions are performed
#                      using a spatial factor multispecies occupancy model.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(coda)
library(tidyverse)

# Read in model output object ---------------------------------------------
load("results/bbs-sfMsPGOcc-2-chain-2022-03-28.R")

# Data Prep ---------------------------------------------------------------
# Read in data used for model fitting
load("data/data-bundle.R")
# Get species in correct order used for model fitting.  
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

# Read in prediction values
# These come from one state at a time, so this code reads in each states 
# prediction values, then binds them all together
# pred.names <- list.files(path = "data/", pattern = '^bbs', full.names = TRUE)
# pred.list <- list()
# for (i in 1:length(pred.names)) {
#    load(pred.names[[i]])
#    pred.list[[i]] <- pred.df
# }
# pred.dat <- bind_rows(pred.list)
# # Remove the few missing values
# pred.dat <- pred.dat[!is.na(pred.dat$elev), ]
# Save the full data set
# save(pred.dat, file = 'data/full-bbs-pred-dat.rda')

# Load the full data set for all prediction locations
load("data/full-bbs-pred-dat.rda")
# Albers equal area in KM
coords.0 <- pred.dat[, c('X', 'Y')] / 1000
# Standardize by values used to fit the model
elev.pred <- (pred.dat$elev - mean(data.list$occ.covs$elev)) / sd(data.list$occ.covs$elev)
forest.pred <- (pred.dat$forest - mean(data.list$occ.covs$forest)) / 
	        sd(data.list$occ.covs$forest)
X.0 <- cbind(1, elev.pred, elev.pred^2, forest.pred)
names(X.0) <- c('int', 'elev', 'elev.2', 'pf')

# Get info on bird communities --------------------------------------------
# Read in community classification data from Bateman et al. (2020)
comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
my.sp.code <- comm.group.dat %>%
  filter(Group %in% c("Eastern.Forests", "Grasslands")) %>%
  filter(Code %in% sp.codes, Season == 'Breeding')
# Classify birds into different communities
ef.birds <- my.sp.code %>%
  filter(Group == 'Eastern.Forests')
ef.indices <- which(sp.codes %in% ef.birds$Code)
n.forest <- length(ef.indices)
grass.birds <- my.sp.code %>%
  filter(Group == 'Grasslands')
grass.indices <- which(sp.codes %in% grass.birds$Code)
n.grass <- length(grass.indices)

# Predict piece by piece across the continental US. -----------------------
J.str <- nrow(X.0)
# Split up the data set.
vals <- split(1:J.str, ceiling(seq_along(1:J.str) / 600))
rich.grass.samples <- array(NA, dim = c(out$n.post, J.str))
rich.forest.samples <- array(NA, dim = c(out$n.post, J.str))
w.samples <- array(NA, dim = c(out$n.post, out$q, J.str))
for (j in 1:length(vals)) {
  print(paste("Currently on set ", j, " out of ", length(vals), sep = ''))
  curr.indx <- vals[[j]]
  out.pred <- predict(out, X.0[curr.indx, ], coords.0[curr.indx, ], n.omp.threads = 5, 
		      verbose = TRUE)
  rich.grass.samples[, curr.indx] <- apply(out.pred$z.0.samples[, grass.indices, ], 
					   c(1, 3), sum)
  rich.forest.samples[, curr.indx] <- apply(out.pred$z.0.samples[, ef.indices, ], 
					   c(1, 3), sum)
  w.samples[, , curr.indx] <- out.pred$w.0.samples
}

# Save results ------------------------------------------------------------
save(rich.forest.samples, 
     rich.grass.samples, w.samples, coords.0, file = 'results/bbs-pred-sfMsPGOcc-rich-results.R')
