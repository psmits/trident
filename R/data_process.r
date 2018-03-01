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
neptune <- list.files(path = '../data', pattern = 'nannotax', full.names = TRUE)
nano <- purrr::map(neptune, read_tsv)
tt <- purrr::map(nano, ~ str_to_lower(names(.)))
tt <- purrr::map(tt, ~ str_replace_all(., ' ', '.'))
nano <- purrr::map2(nano, tt, function(x, y) {names(x) <- y; x})
nano <- purrr::reduce(nano, rbind)
names(nano)[19] <- 'age'
nano[nano$fossil.group == 'FALSE', 'fossil.group'] <- 'F'


nano[, c('latitude', 'longitude')] <- apply(nano[, c('latitude', 'longitude')], 
                                            2, as.numeric)

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

# make a genus_species combo
nano <- nano %>% dplyr::mutate(fullname = str_c(genus, '_', species))


# get all the important geographic range information
sprange <- nano %>%
  group_by(fullname, mybin) %>%
  dplyr::summarise(nocc = n(),  # number of occurrences in bin
                   ncell = n_distinct(cell),  # number of unique cells in bin
                   # latitidinal extent in bin
                   latext = abs(max(latitude) - min(latitude)),  
                   # great circle distance (on ellipsoid) in km in bin
                   maxgcd = max(distm(cbind(longitude, latitude), 
                                      fun = distGeo)) / 1000,
                   fg = unique(fossil.group))


