# lulc-data-prep.R: this script extracts the LULC data from the USGS EROS
#                   Center in 2018 at the locations of the BBS routes that 
#                   were surveyed in 2018.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(stars)
library(sf)

# Grab filename to process ------------------------------------------------
# The file CONUS_A1B_y2018.tif comes from the USGS EROS data center and 
# can be downloaded at https://www.usgs.gov/centers/eros/data. 
file.name <- 'CONUS_A1B_y2018.tif'

# Read in BBS data --------------------------------------------------------
load("data/bbs-data-formatted.R")
# Predict in 1km cells across area of interest ----------------------------
# Get states of interest as an sf object
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Extract the covariates one state at a time ------------------------------
coords.sf <- st_as_sf(coords, 
		      coords = c('Longitude', 'Latitude'), 
		      crs = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
coords.albers <- coords.sf %>%
  st_transform(crs = st_crs(usa))
# 5 km buffer around the starting point of each route. 
buffer.radius <- 5000

# Land Cover Classes ------------------
# 1: Water
# 2: Developed
# 3: Mechanically Disturbed National Forests
# 4: Mechanically Disturbed Other Public Lands
# 5: Mechanically Disturbed Private
# 6: Mining
# 7: Barren
# 8: Deciduous Forest
# 9: Evergreen Forest
# 10: Mixed Forest
# 11: Grassland
# 12: Shrubland
# 13: Cropland
# 14: Hay/Pasture Land
# 15: Herbaceous Wetland
# 16: Woody Wetland
# 17: Perennial Ice/Snow

# Load rasters ------------------------------------------------------------
# NOTE: this line of code reads in the USGS EROS file to extract. You will need
#       to change this depending on where you read in your code.
lulc.curr <- read_stars(paste("~/Dropbox/DSFBGWZ22/data/eros-lulc/", file.name, sep = ''))
# Convert coordinates to the coordinates of the raster. 
coords.ACEA <- coords.albers %>%
  st_transform(crs = st_crs(lulc.curr))
# Buffer the coordinates
coords.ACEA.buffered <- st_buffer(coords.ACEA, dist = buffer.radius)
# Functions to calculate land use proportions -----------------------------
prop.water <- function(a) {
  mean(a %in% 1, na.rm = TRUE)
}
prop.barren <- function(a) {
  mean(a %in% 7, na.rm = TRUE)
}
prop.forest <- function(a) {
  mean(a %in% 8:10, na.rm = TRUE)
}
prop.grass <- function(a) {
  mean(a %in% 11, na.rm = TRUE)
}
prop.shrub <- function(a) {
  mean(a %in% 12, na.rm = TRUE)
}
prop.crop <- function(a) {
  mean(a %in% 13, na.rm = TRUE)
}
prop.hay <- function(a) {
  mean(a %in% 14, na.rm = TRUE)
}
prop.wet <- function(a) {
  mean(a %in% c(15, 16), na.rm = TRUE)
}
prop.devel <- function(a) {
  mean(a == 2, na.rm = TRUE)
}
# Extract the covariate values --------------------------------------------
J <- nrow(coords.albers)
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/100))
lulc.covs <- matrix(NA, nrow(coords.ACEA), 9)
colnames(lulc.covs) <- c('water', 'barren', 'forest', 'grass', 'shrub', 
			 'crop', 'hay', 'wet', 'devel')
for (j in 1:length(vals)) {
  print(j)
  coords.curr <- coords.ACEA.buffered[vals[[j]], ]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.water)
  lulc.covs[vals[[j]], 'water'] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.barren)
  lulc.covs[vals[[j]], 'barren'] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.forest)
  lulc.covs[vals[[j]], 'forest'] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.grass)
  lulc.covs[vals[[j]], 'grass'] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.shrub)
  lulc.covs[vals[[j]], 'shrub'] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.crop)
  lulc.covs[vals[[j]], 'crop'] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.hay)
  lulc.covs[vals[[j]], 'hay'] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.wet)
  lulc.covs[vals[[j]], 'wet'] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.devel)
  lulc.covs[vals[[j]], 'devel'] <- tmp[[1]]
}

# Save extracted covariate values -----------------------------------------
save(lulc.covs, file = 'data/lulc-covs.rda')

