# bbs-data-prep.R: this script takes the RAW BBS data and preps it for analysis
#                  with the spatial factor multi-species occupancy model. The 
#                  script also extracts the covariates for use in the detection
#                  portion of the occupancy model.
# Author: Jeffrey W. Doser

rm(list = ls())
library(tidyverse)
library(lubridate)
library(sp)
library(raster)
library(FedData)
library(sf)
library(stars)
# For linking bird codes with species info.
library(wildlifeR)

# Read in BBS Data --------------------------------------------------------
# These are the 10-stop summary data
bbs.dat <- list.files(path = "data/BBS/States/", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()
# Get associated route data
route.dat <- read.csv("data/BBS/routes.csv")
# Get associated weather data
weather.dat <- read.csv("data/BBS/weather.csv")
# Note that route id is nested within state. 
# Join BBS data with route data
bbs.dat <- left_join(bbs.dat, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# Only grab data from 2018
bbs.2018 <- bbs.dat %>%
  filter(Year == 2018)
# Select columns of interest
bbs.2018 <- bbs.2018 %>%
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) %>%
  dplyr::select(-CountryNum)
# Fill in implicit zeros. For this situation, gives a value for each species in 
# each existing combination of RouteDataID, Latitude, and Longitude
# RouteDataID is a unique identifier for a route within a given year of sampling. 
bbs.2018 <- bbs.2018 %>%
  complete(AOU, nesting(RouteDataID, Latitude, Longitude))
# Replace NAs with 0s for all columns at once.
bbs.2018 <- bbs.2018 %>%
  replace(is.na(.), 0)
# Filter for community of interest. 
aou.info <- AOU_species_codes %>%
  mutate(AOU = spp.num)
bbs.2018 <- left_join(bbs.2018, aou.info, by = c("AOU"))
# For now, try it out with eastern forest and grassland forest species. 
comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
my.sp.code <- comm.group.dat %>%
  filter(Group %in% c("Eastern.Forests", "Grasslands"))
curr.dat <- bbs.2018 %>%
  mutate(alpha.code = as.character(alpha.code)) %>%	
  filter(alpha.code %in% my.sp.code$Code)

# All detection level variables are at the route level, so don't
# need to switch over to wide format. 
# Join data with Weather data
# Get date in proper format
weather.dat <- weather.dat %>% 
  unite('date', sep = '-', Year, Month, Day, remove = FALSE)
weather.dat$date <- as.Date(weather.dat$date, tz = "America/New_York")
# Get julian date of each survey
weather.dat$julian <- as.numeric(format(weather.dat$date, '%j'))
weather.covs <- weather.dat %>%
  filter(RouteDataID %in% unique(curr.dat$RouteDataID))

# Get data in format for spOccupancy --------------------------------------
# Ensure ordering is the same. 
weather.covs <- weather.covs %>%
  arrange(RouteDataID)
curr.dat <- curr.dat %>%
  arrange(alpha.code, RouteDataID)
# Detection Covariates ----------------
det.covs <- list(day = weather.covs$julian, 
		 tod = weather.covs$StartTime, 
		 obs = as.numeric(factor(weather.covs$ObsN)))

# Detection-nondetection data ---------
N <- n_distinct(curr.dat$alpha.code)
sp.codes <- unique(curr.dat$alpha.code)
J <- n_distinct(curr.dat$RouteDataID)
# Number of stop replicates (here reduced to 5, not using the full 50 stop data)
K <- 5
y <- array(NA, dim = c(N, J, K))
for (i in 1:N) {
  tmp <- curr.dat %>%
    filter(alpha.code == sp.codes[i])
  y[i, , ] <- as.matrix(tmp[, c('Count10', 'Count20', 'Count30', 'Count40',
				'Count50')])
  y[i, , ] <- ifelse(y[i, , ] > 0, 1, 0)
}

coords <- curr.dat[1:J, c('Latitude', 'Longitude')]
save(y, det.covs, coords, sp.codes, file = "data/bbs-data-formatted.R")



