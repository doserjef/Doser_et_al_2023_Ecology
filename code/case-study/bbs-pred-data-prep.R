# bbs-pred-data-prep.R: code to extract the prediction grid across the US. 
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(tidyverse)
library(sp)
library(raster)
library(sf)
library(stars)

# Read in BBS data --------------------------------------------------------
load("data/data-bundle.rda")
# Predict in 1km cells across area of interest ----------------------------
# Get states of interest as an sf object
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Get state number from command line --------------------------------------
for (i in 1:nrow(usa)) {
  print(i)
  curr.state <- i
  # Do for one state at a time to speed things up. 
  state.curr <- usa[curr.state, ]
  
  # Grid the area for predictions. 
  # This is 12 x 12 km. 
  grid.pred <- st_as_stars(st_bbox(usa), dx = 12000, dy = 12000)
  # Convert to data frame
  coords.pred <- as.data.frame(grid.pred, center = TRUE)
  # Convert coordinates to an sf object
  coords.pred.sf <- st_as_sf(coords.pred, 
  			   coords = c('x', 'y'), 
  			   crs = st_crs(usa))
  
  # Intersect with current state of interest
  coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(state.curr))
  if (nrow(coords.pred.sf) == 0) {
    print("Current state has zero rows")
    next
  }
  
  # Convert to Lat-Long
  coords.pred.lat.long <- coords.pred.sf %>%
    st_transform(crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
  # Convert to sp object for use in FedData
  coords.pred.lat.long <- as_Spatial(coords.pred.lat.long)
  J <- nrow(coords.pred.sf)
  
  # Save in a data frame for running predict functions in spOccupancy
  coords.df <- st_coordinates(coords.pred.sf)
  pred.df <- data.frame(coords.df) 
  save(pred.df, file = paste('data/bbs-pred-data-', curr.state, '.rda', sep = ''))
}

# Combine all states together ---------------------------------------------
# Read in prediction values
# These come from one state at a time, so this code reads in each states 
# prediction values, then binds them all together
pred.names <- list.files(path = "data/", pattern = '^bbs', full.names = TRUE)
pred.list <- list()
for (i in 1:length(pred.names)) {
   load(pred.names[[i]])
   pred.list[[i]] <- pred.df
}
pred.coords <- bind_rows(pred.list)
# Save the full data set
save(pred.coords, file = 'data/pred-coords.rda')

