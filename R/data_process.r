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


# read in data file and clean column names
neptune <- list.files(path = '../data', pattern = 'nannotax', full.names = TRUE)
nano <- purrr::map(neptune, read_tsv)
tt <- purrr::map(nano, ~ str_to_lower(names(.)))
tt <- purrr::map(tt, ~ str_replace_all(., ' ', '.'))
nano <- purrr::map2(nano, tt, function(x, y) {names(x) <- y; x})
nano <- purrr::reduce(nano, rbind)

# important but impossible to call variable name
names(nano)[19] <- 'age'
names(nano)[23] <- 'plat'
names(nano)[24] <- 'plng'

# one of the fossil group tibble's miscoded fossil group as logical
nano[nano$fossil.group == 'FALSE', 'fossil.group'] <- 'F'

# rbind miscoded as strings
nano$longitude <- parse_double(nano$longitude)
nano$latitude <- parse_double(nano$latitude)

# some don't have paleolat or long
nano <- nano %>%
  filter(!is.na(plat), !is.na(plng))

# assign million year bins (decreasing; youngest lowest)
nano <- nano %>%
  dplyr::arrange(desc(age)) %>%
  dplyr::mutate(mybin = ntile(age, 65))

# assign everything a geographic cell
# uses paleocoordinates
eq <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=m")
globe.map <- readOGR('../data/ne_10m_coastline.shp')  # from natural earth
proj4string(globe.map) <- eq
spatialref <- SpatialPoints(coords = nano[, c('plng', 'plat')],
                            proj4string = eq)  # wgs1984.proj
r <- raster(globe.map, nrows = 50, ncols = 50)
ras <- rasterize(spatialref, r)
# get cell # for each observation
nano$cell <- cellFromXY(ras, 
                        xy = as.data.frame(nano[, c('plng', 'plat')]))

# make a genus_species combo
nano <- nano %>% dplyr::mutate(fullname = str_c(genus, '_', species))


# get all the important geographic range information
sprange <- nano %>%
  group_by(fullname, mybin) %>%
  dplyr::summarise(nocc = n(),  # number of occurrences in bin
                   ncell = n_distinct(cell),  # number of unique cells in bin
                   # latitidinal extent in bin
                   latext = abs(max(plat) - min(plat)),  
                   # great circle distance (on ellipsoid) in km in bin
                   maxgcd = max(distm(cbind(plng, plat), 
                                      fun = distGeo)) / 1000,
                   nprov = n_distinct(longhurst.code),
                   fg = unique(fossil.group))
# still grouped by fullname!

# relative age in bins
longi <- sprange %>%
  dplyr::mutate(relage = abs(mybin - max(mybin))) %>%
  ungroup() %>%
  dplyr::mutate(id = as.numeric(as.factor(fullname)))

# write to file
write_rds(longi, path = '../data/longitude.rds')

# survival dataset
survi <- longi %>%
  group_by(fullname) %>%
  dplyr::summarise(duration = max(relage),
                   cohort = max(mybin),
                   fg = unique(fg),
                   dead = ifelse(min(mybin) == 1, 0, 1)) %>%
  ungroup() %>%
  dplyr::mutate(id = as.numeric(as.factor(fullname)))

# write to file
write_rds(survi, path = '../data/survival.rds')
