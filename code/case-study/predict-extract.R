# pred-extract.R: this script summarizes the full posteriors from the 
#                 prediction results into the means and standard deviations 
#                 for plotting. This drastically reduces the sizes of the output
#                 objects that predict across the continental US.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(coda)
# sfMsPGOcc ---------------------------
load("results/bbs-pred-sfMsPGOcc-rich-results.R")

# Eastern forest bird community richness
rich.forest.mean <- apply(rich.forest.samples, 2, mean)
rich.forest.sd <- apply(rich.forest.samples, 2, sd)
# Grassland bird community richness
rich.grass.mean <- apply(rich.grass.samples, 2, mean)
rich.grass.sd <- apply(rich.grass.samples, 2, sd)
# Spatial factors
w.mean <- apply(w.samples, c(2, 3), mean)
w.sd <- apply(w.samples, c(2, 3), sd)

save(coords.0, rich.forest.mean, rich.forest.sd, rich.grass.mean, 
     rich.grass.sd, w.mean, w.sd, file = 'results/bbs-pred-sfMsPGOcc-summary.R')

# lfMsPGOcc ---------------------------
load("results/bbs-pred-lfMsPGOcc-rich-results.R")
# Eastern forest bird community richness
rich.forest.mean <- apply(rich.forest.samples, 2, mean)
rich.forest.sd <- apply(rich.forest.samples, 2, sd)
# Grassland bird community richness
rich.grass.mean <- apply(rich.grass.samples, 2, mean)
rich.grass.sd <- apply(rich.grass.samples, 2, sd)
save(coords.0, rich.forest.mean, rich.forest.sd, rich.grass.mean, 
     rich.grass.sd, file = 'results/bbs-pred-lfMsPGOcc-summary.R')

# sfJSDM ---------------------------
load("results/bbs-pred-sfJSDM-rich-results.R")

# Eastern forest bird community richness
rich.forest.mean <- apply(rich.forest.samples, 2, mean)
rich.forest.sd <- apply(rich.forest.samples, 2, sd)
# Grassland bird community richness
rich.grass.mean <- apply(rich.grass.samples, 2, mean)
rich.grass.sd <- apply(rich.grass.samples, 2, sd)
# Spatial factors
w.mean <- apply(w.samples, c(2, 3), mean)
w.sd <- apply(w.samples, c(2, 3), sd)
save(coords.0, rich.forest.mean, rich.forest.sd, rich.grass.mean, 
     rich.grass.sd, w.mean, w.sd, file = 'results/bbs-pred-sfJSDM-summary.R')
