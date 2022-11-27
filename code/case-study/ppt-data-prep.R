# ppt-data-prep.R: this script extracts precipitation data from PRISM
#                  at the locations of the 2619 BBS routes in 2018/2017 and
#                  across the continental US for subsequent prediction.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(stars)
library(prism)

# Load formatted BBS data -------------------------------------------------
# Loads data object formatted for spOccupancy. Here only using the coords 
# portion of the object to extract the spatial coordinates. 
load("data/bbs-data-formatted.R")
# Convert coords to Albers Equal Area. 
coords.sf <- st_as_sf(data.frame(coords),
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
# Get coordinates in Albers Equal Area
coords.albers <- st_coordinates(coords.sf.albers)

# Set PRISM directory -----------------------------------------------------
prism_set_dl_dir("data/prism-climate")

# Extract data at BBS routes ----------------------------------------------
# 12 months in a year. 
ppt.vals <- matrix(NA, nrow = nrow(coords.albers), ncol = 12)
# NOTE: can change this to extract values for different years. For the case 
#       study in the associated manuscript, we extracted the values for 2017 
#       and 2018, and then used data spanning from June 2017 - May 2018.
curr.year <- 2018
get_prism_monthlys(type = "ppt", year = curr.year, mon = 1:12, keepZip = FALSE)
# Get file name to data of interest
# Loop through all the months
for (i in 1:ncol(ppt.vals)) {
  ppt.curr.path <- prism_archive_subset("ppt", "monthly", years = curr.year, mon = i)
  # Get absolute file path
  ppt.curr.abs <- pd_to_file(ppt.curr.path)
  # Download the raster
  ppt.curr <- read_stars(ppt.curr.abs)
  # Get coordinates with buffer for extracting the climate data
  coords.proj <- coords.sf.albers %>%
    st_transform(crs = st_crs(ppt.curr))
  # Extract precipitation value at the starting location of the route. 
  tmp <- st_extract(ppt.curr, at = coords.proj)
  ppt.vals[, i] <- tmp[[1]]
}

# Extract data at prediction locations ------------------------------------
load("data/pred-coords.rda")
# Albers equal area in KM
coords.0 <- pred.coords / 1000
coords.0.sf <- st_as_sf(coords.0, 
			coords = c('X', 'Y'), 
			crs = st_crs(coords.sf.albers))
ppt.pred.vals <- matrix(NA, nrow = nrow(coords.0.sf), ncol = 12)
for (i in 1:ncol(ppt.vals)) {
  ppt.curr.path <- prism_archive_subset("ppt", "monthly", years = curr.year, mon = i)
  # Get absolute file path
  ppt.curr.abs <- pd_to_file(ppt.curr.path)
  # Download the raster
  ppt.curr <- read_stars(ppt.curr.abs)
  # Get coordinates with buffer for extracting the climate data
  coords.proj <- coords.0.sf %>%
    st_transform(crs = st_crs(ppt.curr))
  # Extract precipitation value at the starting location of the route. 
  tmp <- st_extract(ppt.curr, at = coords.proj)
  ppt.pred.vals[, i] <- tmp[[1]]
}

# Save values to hard drive -----------------------------------------------
save(ppt.vals, file = paste("data/climate-data/ppt-", curr.year, ".rda", sep = ''))
save(ppt.pred.vals, file = paste("data/climate-data/ppt-", curr.year, "-pred.rda", sep = ''))

