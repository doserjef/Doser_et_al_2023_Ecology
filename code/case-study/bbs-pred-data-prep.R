# bbs-pred-data-prep.R: prepare covariate data for predictions across
#                       the continental US.
# Author: Jeffrey W. Doser
# Citation: 

rm(list = ls())
library(tidyverse)
library(sp)
library(raster)
library(FedData)
library(sf)
library(stars)

# Read in BBS BTNW data ---------------------------------------------------
load("data/data-bundle.R")
# Predict in 1km cells across area of interest ----------------------------
# Get states of interest as an sf object
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Get state number from command line --------------------------------------
curr.state <- as.numeric(commandArgs(trailingOnly = TRUE))
if(length(curr.state) == 0) base::stop('Need to tell me which state to grab data from')
# Do for one state at a time to speed things up. 
state.curr <- usa[curr.state, ]

# Grid the area for predictions. 
# This is currently 11 x 11 km. 
grid.pred <- st_as_stars(st_bbox(usa), dx = 12000, dy = 12000)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred, 
			   coords = c('x', 'y'), 
			   crs = st_crs(usa))

# Intersect with current state of interest
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(state.curr))
if (nrow(coords.pred.sf) == 0) stop("Current state has zero rows")

# Convert to Lat-Long
coords.pred.lat.long <- coords.pred.sf %>%
  st_transform(crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
# Convert to sp object for use in FedData
coords.pred.lat.long <- as_Spatial(coords.pred.lat.long)
J <- nrow(coords.pred.sf)
vals <- split(1:J, ceiling(seq_along(1:J) / 24))
# Function to calculate percent forested area
props <- function(a, na.rm = TRUE) {
  my.sum <- sum(!is.na(a))	
  prop.for <- sum(a %in% c(41, 42, 43), na.rm = na.rm) / my.sum
  return(prop.for)
}
# Merge rasters all together
elev.pred <- rep(NA, J)
for.pred <- rep(NA, J)
ned.dat <- list()
nlcd.dat <- list()
for (i in 1:length(vals)) {
  print(paste("Currently on iteration ", i, " out of ", length(vals), sep = ''))
  tryCatch({
    ned.dat[[i]] <- get_ned(template = coords.pred.lat.long[vals[[i]], ], 
          		  label = paste('btnw', i))
    tmp <- extract(ned.dat[[i]], coords.pred.lat.long[vals[[i]], ])
    nlcd.dat[[i]] <- get_nlcd(template = coords.pred.lat.long[vals[[i]], ], 
          		    label = paste('btnw', i), year = 2016)
    tmp.nlcd <- extract(nlcd.dat[[i]], coords.pred.lat.long[vals[[i]], ], buffer = 1000, fun = props)
    elev.pred[vals[[i]]] <- tmp
    for.pred[vals[[i]]] <- tmp.nlcd
  }, error = function(e){})
}

coords.df <- st_coordinates(coords.pred.sf)
pred.df <- data.frame(coords.df, 
		      elev = elev.pred, 
		      forest = for.pred)
save(pred.df, file = paste('data/bbs-pred-data-', curr.state, '.rda', sep = ''))

