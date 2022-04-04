# summary.R: this script summarizes results from the full spatial factor
#            multi-species occupancy model for the BBS case study of 98 
#            bird species across the continental US. 
# Author: Jeffrey W. Doser

# Clear out the workspace
rm(list = ls())
# Load the packages
library(tidyverse)
library(spOccupancy)
library(coda)
library(sf)
library(viridis)
library(corrplot)
library(ggpubr)
library(stars)

# Read in files -----------------------------------------------------------
load("results/bbs-sfMsPGOcc-2-chain-2022-03-29.R")

# Read in raw data --------------------------------------------------------
load("data/data-bundle.R")
# Get data in the order used to fit the model. 
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

# Get coordinates in a projection instead of lat-long ---------------------
coords.sf <- st_as_sf(data.frame(data.list$coords),
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Get coordinates in Albers Equal Area
coords.albers <- st_coordinates(coords.sf.albers)
data.list$coords <- coords.albers

# Exploratory plots of model results --------------------------------------
# Chan change species code to different things as below, either by uncommenting 
# the loop and going across all species, or just doing it manually
#for (i in 1:length(sp.codes)) {
  #name.sp <- sp.codes[i]
  name.sp <- "BLJA"
  curr.sp <- which(sp.codes == name.sp)
  curr.psi.mean <- apply(out$psi.samples[, curr.sp, ], 2, mean)
  curr.df <- coords.sf.albers
  curr.df$psi <- curr.psi.mean
  # Occurrence
  plot(ggplot(curr.df) + 
    geom_sf(aes(col = psi)) + 
    scale_color_viridis() + 
    theme_bw(base_size = 18) +
    labs(title = name.sp) )
  #Sys.sleep(2)
#}

# Predicted Species Richness ----------------------------------------------
# Load in summary statistics of prediction results
load("results/bbs-pred-sfMsPGOcc-summary.R")
# Load in prediction data. 
load("data/full-bbs-pred-dat.rda")
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")

# Create data frame for plotting
plot.stars.df <- data.frame(x = coords.0[, 1], y = coords.0[, 2],
			    ef.rich.mean = rich.forest.mean,
			    ef.rich.sd = rich.forest.sd, 
			    grass.rich.mean = rich.grass.mean, 
                            grass.rich.sd = rich.grass.sd,
                            elev = pred.dat$elev, 
                            forest = pred.dat$forest, 
                            w.1.mean = w.mean[1, ], 
                            w.2.mean = w.mean[2, ], 
                            w.3.mean = w.mean[3, ], 
                            w.4.mean = w.mean[4, ], 
                            w.5.mean = w.mean[5, ], 
                            w.1.sd = w.sd[1, ], 
                            w.2.sd = w.sd[2, ], 
                            w.3.sd = w.sd[3, ], 
                            w.4.sd = w.sd[4, ], 
                            w.5.sd = w.sd[5, ]) 

# Load and add in the non-spatial prediction summary statistics (lfMsPGOcc)
load('results/bbs-pred-lfMsPGOcc-summary.R')
plot.stars.df$ef.rich.mean.non.sp <- rich.forest.mean
plot.stars.df$ef.rich.sd.non.sp <- rich.forest.sd
plot.stars.df$grass.rich.mean.non.sp <- rich.grass.mean
plot.stars.df$grass.rich.sd.non.sp <- rich.grass.sd
plot.stars.df$ef.rich.diff <- plot.stars.df$ef.rich.mean - plot.stars.df$ef.rich.mean.non.sp
plot.stars.df$grass.rich.diff <- plot.stars.df$grass.rich.mean - plot.stars.df$grass.rich.mean.non.sp

# Load and add in the spatial JSDM estimates (sfJSDM)
load('results/bbs-pred-sfJSDM-summary.R')
plot.stars.df$ef.rich.mean.sfJSDM <- rich.forest.mean
plot.stars.df$ef.rich.sd.sfJSDM <- rich.forest.sd
plot.stars.df$grass.rich.mean.sfJSDM <- rich.grass.mean
plot.stars.df$grass.rich.sd.sfJSDM <- rich.grass.sd
plot.stars.df$ef.rich.diff.sfJSDM <- plot.stars.df$ef.rich.mean - plot.stars.df$ef.rich.mean.sfJSDM
plot.stars.df$grass.rich.diff.sfJSDM <- plot.stars.df$grass.rich.mean - plot.stars.df$grass.rich.mean.sfJSDM

# Convert to a stars object for plotting. 
pred.stars <- st_as_stars(plot.stars.df, dims = c('x', 'y'))

# Spatial richness plot ---------------------------------------------------
# Eastern forest bird mean richness
ef.rich.pred.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))
# Eastern forest bird sd richness
ef.rich.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.sd),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

# Grassland bird mean richness
grass.rich.pred.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

# Create Figure 1
ggarrange(ef.rich.pred.plot, ef.rich.sd.plot, 
	  grass.rich.pred.plot, grass.rich.sd.plot, 
	  nrow = 2, ncol = 2, labels = c("(A) Eastern Forest Mean Richness", 
					 "(B) Eastern Forest SD Richness", 
					 "(C) Grassland Mean Richness", 
					 "(D) Grassland SD Richness"),
	  font.label = list(size = 25), hjust = -0.15)

 # Grassland bird sd richness
grass.rich.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.sd),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

# Save to hard drive if you like.
ggsave(device = 'pdf', filename = 'figures/Fig1.pdf', height = 14, width = 20)

# Mean difference between sfMsPGOcc and lfMsPGocc estimates
# Eastern forests
ef.rich.diff.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.diff),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))
# Grassland
grass.rich.diff.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.diff),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

# Mean differences between sfMsPGOcc and sfJSDM
# Eastern forests
ef.rich.diff.sfJSDM.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.diff.sfJSDM),
	     interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))
# Grassland
grass.rich.diff.sfJSDM.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.diff.sfJSDM),
	     interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))
# Figure 2
ggarrange(ef.rich.diff.plot, grass.rich.diff.plot,
	  ef.rich.diff.sfJSDM.plot, grass.rich.diff.sfJSDM.plot, 
	  nrow = 2, ncol = 2, labels = c("(A) Eastern Forest sfMsPGOcc - lfMsPGOcc", 
					 "(B) Grassland sfMsPGOcc - lfMsPGOcc", 
					 "(C) Eastern Forest sfMsPGOcc - sfJSDM", 
					 "(D) Grassland sfMsPGOcc - sfJSDM"), 
	  font.label = list(size = 25), hjust = -0.15)
ggsave(device = 'pdf', filename = 'figures/Fig2.pdf', height = 14, width = 20)

# Spatial Factors ---------------------------------------------------------
# First spatial process ---------------
w.1.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.1.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.1.plot
# Part of Figure 3
ggsave(device = 'pdf', filename = 'figures/w-1-fig.pdf', height = 7, width = 10)
# Second spatial process --------------
w.2.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.2.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.2.plot
# Part of Figure 3 
ggsave(device = 'pdf', filename = 'figures/w-2-fig.pdf', height = 7, width = 10)
# Third spatial process ---------------
w.3.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.3.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.3.plot
# Part of Figure S1
ggsave(device = 'pdf', filename = 'figures/w-3-fig.pdf', height = 7, width = 10)
# Fourth spatial process --------------
w.4.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.4.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.4.plot
# Part of Figure S2
ggsave(device = 'pdf', filename = 'figures/w-4-fig.pdf', height = 7, width = 10)
# Fifth spatial process ---------------
w.5.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.5.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.5.plot
# Part of Figure S3
ggsave(device = 'pdf', filename = 'figures/w-5-fig.pdf', height = 7, width = 10)

# Factor loading densities ------------------------------------------------
# Species guild information from Bateman et al. (2020)
comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
my.sp.code <- comm.group.dat %>%
  filter(Group %in% c("Eastern.Forests", "Grasslands")) %>%
  filter(Code %in% sp.codes, Season == 'Breeding')
ef.birds <- my.sp.code %>%
  filter(Group == 'Eastern.Forests')
ef.indices <- which(sp.codes %in% ef.birds$Code)
grass.birds <- my.sp.code %>%
  filter(Group == 'Grasslands')
grass.indices <- which(sp.codes %in% grass.birds$Code)

# Means of the species-specific factor loadings
lambda.means <- apply(out$lambda.samples, 2, mean)
N <- dim(data.list$y)[1]
q <- dim(data.list$y)[3]
lambda.means <- matrix(lambda.means, N, q)

lambda.ef.means <- lambda.means[ef.indices, ]
lambda.grass.means <- lambda.means[grass.indices, ]
lambda.plot.df <- data.frame(rbind(lambda.ef.means, lambda.grass.means), 
			     comm = c(rep('Eastern Forest', nrow(lambda.ef.means)), 
				      rep('Grassland', nrow(lambda.grass.means))))
names(lambda.plot.df) <- c('Factor1', 'Factor2', 'Factor3', 'Factor4', 
			   'Factor5', 'Community')
# Display a summary of the factor loadings in the two communities
lambda.plot.df %>%
  group_by(Community) %>%
  summarize(sum.1.positive = sum(Factor1 > 0) / n(), 
	    sum.2.positive = sum(Factor2 > 0) / n(), 
            sum.3.positive = sum(Factor3 > 0) / n(), 
            sum.4.positive = sum(Factor4 > 0) / n(), 
            sum.5.positive = sum(Factor5 > 0) / n())
# Factor 1
lambda.1 <- ggplot(lambda.plot.df, aes(x = Factor1, fill = Community)) + 
  geom_density(alpha = 0.5, col = 'black') + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 1", y = "Density") +
  theme(legend.position = 'top')
# Part of Figure 3
lambda.1
ggsave(device = 'pdf', filename = 'figures/lambda-1-fig.pdf', height = 5, width = 6)
# Factor 2
lambda.2 <- ggplot(lambda.plot.df, aes(x = Factor2, fill = Community)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 2", y = "Density")+ 
  theme(legend.position = 'top')
# Part of Figure 3
lambda.2
ggsave(device = 'pdf', filename = 'figures/lambda-2-fig.pdf', height = 5, width = 6)
# Factor 3
lambda.3 <- ggplot(lambda.plot.df, aes(x = Factor3, fill = Community)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 3", y = "Density") + 
  theme(legend.position = 'top')
# Part of Figure S1
lambda.3
ggsave(device = 'pdf', filename = 'figures/lambda-3-fig.pdf', height = 5, width = 6)
# Factor 4
lambda.4 <- ggplot(lambda.plot.df, aes(x = Factor4, fill = Community)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 4", y = "Density") + 
  theme(legend.position = 'top')
# Part of Figure S2
lambda.4
ggsave(device = 'pdf', filename = 'figures/lambda-4-fig.pdf', height = 5, width = 6)
# Factor 5
lambda.5 <- ggplot(lambda.plot.df, aes(x = Factor5, fill = Community)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 5", y = "Density") + 
  theme(legend.position = 'top')
lambda.5
# Part of Figure S3
ggsave(device = 'pdf', filename = 'figures/lambda-5-fig.pdf', height = 5, width = 6)

# Compare WAIC across models ----------------------------------------------
# lfJSDM ---------------------------
load("results/bbs-lfJSDM-2-chain-2022-03-27.R")
out.2.lfJSDM <- out
waicOcc(out.2.lfJSDM)
# sfJSDM ---------------------------
load("results/bbs-sfJSDM-2-chain-2022-03-28.R")
out.2.sfJSDM <- out
waicOcc(out.2.sfJSDM)
# msPGOcc ---------------------------
load("results/bbs-msPGOcc-2-chain-2022-03-29.R")
out.2.msPGOcc <- out
waicOcc(out.2.msPGOcc)
# lfMsPGOcc ---------------------------
load("results/bbs-lfMsPGOcc-2-chain-2022-03-28.R")
out.2.lfMsPGOcc <- out
waicOcc(out.2.lfMsPGOcc)
# sfMsPGOcc ---------------------------
load("results/bbs-sfMsPGOcc-2-chain-2022-03-29.R")
out.2.sfMsPGOcc <- out
waicOcc(out.2.sfMsPGOcc)

# Out-of-sample model-validation results ----------------------------------
# This loads an object deviance.df that consists of all the out-of-sample
# deviance metrics.
load("results/out-of-sample-deviance.R")

# Average deviance for the raw data across species
deviance.df %>%
  filter(type.deviance == 'y') %>%
  group_by(model.fit) %>%
  summarize(dev.comm = mean(deviance))
# Average deviance across species and the model-generated z values from 
# three different occupancy models
deviance.df %>%
  filter(type.deviance != 'y') %>%
  group_by(model.fit) %>%
  summarize(dev.comm = mean(deviance))
