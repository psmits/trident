# because of course
library(tidyverse)

# geolocation stuff
library(sp)
library(dismo)
library(raster)
library(rgdal)
library(XML)
library(maptools)
library(geosphere)
library(fields)


# read in data file and clean column names
nano <- read_tsv(file = '../data/nannotax_2018-02-27_18-00-36.csv')
names(nano) <- str_to_lower(names(nano))
names(nano) <- str_replace_all(names(nano), ' ', '.')
names(nano)[19] <- 'age'


# assign million year bins
nano <- nano %>%
  mutate(mybin = ntile(age, 65))


# assign everything a geographic cell
eq <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=m")
globe.map <- readOGR('../data/ne_10m_coastline.shp')  # from natural earth
proj4string(globe.map) <- eq
spatialref <- SpatialPoints(coords = nano[, c('longitude', 'latitude')],
                            proj4string = eq)  # wgs1984.proj
r <- raster(globe.map, nrows = 50, ncols = 50)
ras <- rasterize(spatialref, r)
# get cell # for each observation
nano$cell <- cellFromXY(ras, 
                        xy = as.data.frame(nano[, c('longitude', 'latitude')]))



# get all the important geographic range information
sprange <- nano %>%
  group_by(genus, species, mybin) %>%
  dplyr::summarise(nocc = n(),  # number of occurrences in bin
                   ncell = n_distinct(cell),  # number of unique cells in bin
                   # latitidinal extent in bin
                   latext = abs(max(latitude) - min(latitude)),  
                   # great circle distance (on ellipsoid) in km in bin
                   maxgcd = max(distm(cbind(longitude, latitude), 
                                      fun = distGeo)) / 1000) 

