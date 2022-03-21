# summary.R: this script summarizes results from the full spatially explicit
#            joint species distribution model that explicitly accounts for 
#            imperfect detection. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(spOccupancy)
library(coda)
library(sf)
library(viridis)
library(corrplot)
library(ggpubr)
library(stars)

# Read in files -----------------------------------------------------------
# Chains were run separately to allow for parallel runs instead of using the
# n.chains argument in sfMsPGOcc. 
load("results/bbs-sfMsPGOcc-1-chain-2022-02-09.R")
out.1 <- out
# load("results/bbs-sfMsPGOcc-2-chain-2022-02-04.R")
# out.2 <- out
# load("results/bbs-sfMsPGOcc-3-chain-2022-02-04.R")
# out.3 <- out
# load("results/bbs-lfMsPGOcc-1-chain-2022-02-04.R")
# out.no.sp <- out
# out.no.sp$psiRE <- FALSE
# # Necessary to work with current version of spOccupancy.
# out.1$psiRE <- FALSE
# out.2$psiRE <- FALSE
# out.3$psiRE <- FALSE

# For now
out <- out.1
# Read in raw data --------------------------------------------------------
load("data/data-bundle.R")

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

# Create some basic plots to compare spatial vs. non-spatial maps. 
# for (i in 1:length(sp.codes)) {
  # Common grassland: DICK, EAKI, EAME, GRSP (moderate), BOBO (moderate)
  # Common forest: BLJA, REVI, PIWO (less so), BTNW
  # Ordering for Factors: REVI, GRSP, PIWO, EAME, BTNW
  name.sp <- 'BTNW'
  curr.sp <- which(sp.codes == name.sp)
  # name.sp <- sp.codes[i]
  # curr.sp <- i	
  curr.psi.mean <- apply(out$psi.samples[, curr.sp, ], 2, mean)
  # curr.lf.psi.mean <- apply(out.no.sp$psi.samples[, curr.sp, ], 2, mean)
  curr.df <- coords.sf.albers
  curr.df$psi <- curr.psi.mean
  # curr.df$psi.no.sp <- curr.lf.psi.mean
  # Occurrence
  psi.sp.plot <- ggplot(curr.df) + 
    geom_sf(aes(col = psi)) + 
    scale_color_viridis() + 
    theme_bw(base_size = 18) +
    labs(title = name.sp) 
  # psi.no.sp.plot <- ggplot(curr.df) + 
  #   geom_sf(aes(col = psi.no.sp)) + 
  #   scale_color_viridis() + 
  #   theme_bw(base_size = 18) 
  # ggarrange(psi.sp.plot, psi.no.sp.plot, labels = c('Spatial', 'Non-spatial'), ncol = 2)
  psi.sp.plot
# }

# Factor loadings ---------------------------------------------------------
lambda.means <- apply(out$lambda.samples, 2, mean)
N <- dim(data.list$y)[1]
q <- dim(data.list$y)[3]
lambda.means <- matrix(lambda.means, N, q)
lambda.means[curr.sp, ]

# Predicted Species Richness ----------------------------------------------
load("results/bbs-pred-sfMsPGOcc-summary.R")
load("data/full-bbs-pred-dat.rda")
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")

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

# Add in the non-spatial estimates
load('results/bbs-pred-lfMsPGOcc-summary.R')
plot.stars.df$ef.rich.mean.non.sp <- rich.forest.mean
plot.stars.df$ef.rich.sd.non.sp <- rich.forest.sd
plot.stars.df$grass.rich.mean.non.sp <- rich.grass.mean
plot.stars.df$grass.rich.sd.non.sp <- rich.grass.sd
plot.stars.df$ef.rich.diff <- plot.stars.df$ef.rich.mean - plot.stars.df$ef.rich.mean.non.sp
plot.stars.df$grass.rich.diff <- plot.stars.df$grass.rich.mean - plot.stars.df$grass.rich.mean.non.sp

# Add in the spatial JSDM estimates
load('results/bbs-pred-sfJSDM-summary.R')
plot.stars.df$ef.rich.mean.sfJSDM <- rich.forest.mean
plot.stars.df$ef.rich.sd.sfJSDM <- rich.forest.sd
plot.stars.df$grass.rich.mean.sfJSDM <- rich.grass.mean
plot.stars.df$grass.rich.sd.sfJSDM <- rich.grass.sd
plot.stars.df$ef.rich.diff.sfJSDM <- plot.stars.df$ef.rich.mean - plot.stars.df$ef.rich.mean.sfJSDM
plot.stars.df$grass.rich.diff.sfJSDM <- plot.stars.df$grass.rich.mean - plot.stars.df$grass.rich.mean.sfJSDM

pred.stars <- st_as_stars(plot.stars.df, dims = c('x', 'y'))

# Spatial richness plot ---------------------------------------------------
ef.rich.pred.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

ef.rich.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.sd),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

grass.rich.pred.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

grass.rich.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.sd),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

# Difference between sfMsPGOcc and lfMsPGocc estimates
ef.rich.diff.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.diff),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  # scale_fill_gradient2(midpoint = 0, low = 'red', mid = 'white', high = 'blue', 
  # 	               na.value = NA) + 
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

grass.rich.diff.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.diff),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  # scale_fill_gradient2(midpoint = 0, low = 'red', mid = 'white', high = 'blue', 
  # 	               na.value = NA) + 
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))

# Full plot of all four maps
ggarrange(ef.rich.pred.plot, ef.rich.sd.plot, 
	  grass.rich.pred.plot, grass.rich.sd.plot, 
	  ef.rich.diff.plot, grass.rich.diff.plot, 
	  nrow = 3, ncol = 2, labels = c("(A) Eastern Forest Mean Richness", "(B) Eastern Forest SD Richness", "(C) Grassland Mean Richness", "(D) Grassland SD Richness", "(E) Eastern Forest Spatial - Nonspatial", "(F) Grassland Spatial - Nonspatial"), 
	  font.label = list(size = 25), hjust = -0.15)
ggsave(device = 'pdf', filename = 'figures/sfMsPGOcc-richness-plot.pdf', height = 20.75, width = 20)

# Non-spatial richness plot -----------------------------------------------
ef.rich.non.sp.pred.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.mean.non.sp),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))

ef.rich.non.sp.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.sd.non.sp),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "SD") +
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))

grass.rich.non.sp.pred.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.mean.non.sp),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))

grass.rich.non.sp.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.sd.non.sp),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "SD") +
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))

# Full plot of all four maps
ggarrange(ef.rich.non.sp.pred.plot, ef.rich.non.sp.sd.plot, 
	  grass.rich.non.sp.pred.plot, grass.rich.non.sp.sd.plot, 
	  nrow = 2, ncol = 2, labels = c("(A)", "(B)", "(C)", "(D)"), 
	  font.label = list(size = 25))
ggsave(device = 'pdf', filename = 'figures/lfMsPGOcc-richness-plot.pdf', height = 15, width = 20)

# Plot the differences between sfMsPGOcc and sfJSDM -----------------------
ef.rich.diff.sfJSDM.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = ef.rich.diff.sfJSDM),
	     interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  # scale_fill_gradient2(midpoint = 0, low = 'red', mid = 'white', high = 'blue', 
  # 	               na.value = NA) + 
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))
grass.rich.diff.sfJSDM.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = grass.rich.diff.sfJSDM),
	     interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  # scale_fill_gradient2(midpoint = 0, low = 'red', mid = 'white', high = 'blue', 
  # 	               na.value = NA) + 
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.92, 0.3), 
        legend.background = element_rect(fill = NA))
ggarrange(ef.rich.diff.sfJSDM.plot, grass.rich.diff.sfJSDM.plot, 
	  nrow = 1, ncol = 2, labels = c("(A) Eastern Forest", "(B) Grassland"), 
	  font.label = list(size = 25))
ggsave(device = 'pdf', filename = 'figures/richness-sfMsPGOcc-sfJSDM-diff-plot.pdf', height = 7, width = 20)

# Spatial Factors ---------------------------------------------------------
# First spatial process
# A pretty clear longitudinal gradient, geared towards Eastern US
w.1.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.1.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.1.plot
ggsave(device = 'pdf', filename = 'figures/w-1-fig.pdf', height = 7, width = 10)
# Second spatial process
# Kind of a weird one: NE + west. This makes total sense. These are pretty much 
# high latitude birds. Those that have ranges that span into forests in Canada, 
# then fall down lower into the PNW.  
w.2.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.2.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.2.plot
ggsave(device = 'pdf', filename = 'figures/w-2-fig.pdf', height = 7, width = 10)
# Third spatial process
# Latitudinal gradient
w.3.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.3.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.3.plot
ggsave(device = 'pdf', filename = 'figures/w-3-fig.pdf', height = 7, width = 10)
# Fourth spatial process
# Prairie pothole region. 
w.4.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.4.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.4.plot
ggsave(device = 'pdf', filename = 'figures/w-4-fig.pdf', height = 7, width = 10)
# Fifth spatial process
w.5.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.5.mean),interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0) +
  #scale_fill_viridis(na.value = NA) +
  scale_fill_gradientn(colors = magma(10), na.value = NA) +
  labs(x = "Longitude", y = "Latitude", fill = "Mean") +
  theme_bw(base_size = 25) + 
  theme(legend.position = c(0.92, 0.205),
	legend.background = element_rect(fill = 'white', color = 'gray', size = 0.3))
w.5.plot
ggsave(device = 'pdf', filename = 'figures/w-5-fig.pdf', height = 7, width = 10)

# Densities of factor loadings by species community -----------------------
# Species guild information
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

# Assess raw data prevalence
raw.prop <- apply(data.list$y, 1, mean)
names(raw.prop) <- sp.codes
# Grassland
raw.prop[grass.indices]
# Common grassland: DICK, EAKI, EAME, GRSP (moderate), BOBO (moderate)
# Forest
raw.prop[ef.indices]
# Common forest: BLJA, REVI, PIWO (less so), BTNW

# Ordering: REVI, DICK


lambda.ef.means <- lambda.means[ef.indices, ]
lambda.grass.means <- lambda.means[grass.indices, ]
lambda.plot.df <- data.frame(rbind(lambda.ef.means, lambda.grass.means), 
			     comm = c(rep('Eastern Forest', nrow(lambda.ef.means)), 
				      rep('Grassland', nrow(lambda.grass.means))))
names(lambda.plot.df) <- c('Factor1', 'Factor2', 'Factor3', 'Factor4', 
			   'Factor5', 'Community')
lambda.plot.df %>%
  group_by(Community) %>%
  summarize(sum.1.positive = sum(Factor1 > 0) / n(), 
	    sum.4.positive = sum(Factor4 > 0) / n())
# Factor 1
lambda.1 <- ggplot(lambda.plot.df, aes(x = Factor1, fill = Community)) + 
  geom_density(alpha = 0.5, col = 'black') + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 1", y = "Density") +
  theme(legend.position = 'top')
lambda.1
ggsave(device = 'pdf', filename = 'figures/lambda-1-fig.pdf', height = 5, width = 6)
# Factor 2
lambda.2 <- ggplot(lambda.plot.df, aes(x = Factor2, fill = Community)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 2", y = "Density")+ 
  theme(legend.position = 'top')
lambda.2
ggsave(device = 'pdf', filename = 'figures/lambda-2-fig.pdf', height = 5, width = 6)
# Factor 3
lambda.3 <- ggplot(lambda.plot.df, aes(x = Factor3, fill = Community)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 3", y = "Density") + 
  theme(legend.position = 'top')
lambda.3
ggsave(device = 'pdf', filename = 'figures/lambda-3-fig.pdf', height = 5, width = 6)
# Factor 4
lambda.4 <- ggplot(lambda.plot.df, aes(x = Factor4, fill = Community)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_viridis_d() + 
  theme_bw(base_size = 18) + 
  labs(x = "Factor 4", y = "Density") + 
  theme(legend.position = 'top')
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
ggsave(device = 'pdf', filename = 'figures/lambda-5-fig.pdf', height = 5, width = 6)

# Out-of-sample Performance -----------------------------------------------
load("results/out-of-sample-deviance.R")

# Average deviance for each type across species
deviance.df %>%
  group_by(model.fit, type.deviance) %>%
  summarize(dev.comm = mean(deviance)) %>%
  arrange(type.deviance)
