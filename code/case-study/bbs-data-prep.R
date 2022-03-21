# bbs-data-prep.R: this script takes the RAW BBS data and preps it for analysis
#                  with the binary spatial factor NNGP model. 
# Author: Jeffrey W. Doser
# Citation: 

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
# TODO: need to go back here and make sure this is all correct. 
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

# Occurrence covariates
# Since curr.dat is ordered by species, then site, can just grab the first J values
coords <- curr.dat[1:J, c('Latitude', 'Longitude')]
# Need to break this up into pieces to make it work, otherwise it times out
coords.sp <- data.frame(coords)
coords.sp <- coords.sp %>% arrange(Longitude, Latitude)
coordinates(coords.sp) <- ~Longitude + Latitude
proj4string(coords.sp) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/20))
elev.cov <- rep(0, J)
for.cov <- rep(0, J)
# Get proportion of forest
props <- function(a, na.rm = TRUE) {
  my.sum <- sum(!is.na(a))	
  prop.for <- sum(a %in% c(41, 42, 43), na.rm = na.rm) / my.sum
  return(prop.for)
}
ned.dat <- list()
nlcd.dat <- list()
for (i in 1:length(vals)) {
  print(paste("Currently on iteration ", i, " out of ", length(vals), sep = ''))
  ned.dat[[i]] <- get_ned(template = coords.sp[vals[[i]], ], label = paste('btbw', i))
  nlcd.dat[[i]] <- get_nlcd(template = coords.sp[vals[[i]], ], label = paste('btbw', i), year = 2016)
  elev.cov[vals[[i]]] <- extract(ned.dat[[i]], coords.sp[vals[[i]], ])
  for.cov[vals[[i]]] <- extract(nlcd.dat[[i]], coords.sp[vals[[i]], ], buffer = 1000, fun = props)
}

occ.covs <- data.frame(elev = elev.cov, 
		       forest = for.cov)

data.list <- list(y = y, det.covs = det.covs, 
		  occ.covs = occ.covs, coords = coords)

save(data.list, sp.codes, file = "data/data-bundle.R")



