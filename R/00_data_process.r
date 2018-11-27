# because of course
library(tidyverse)
library(janitor)

# geolocation stuff
library(sp)
library(dismo)
library(raster)
library(rgdal)
library(XML)
library(maptools)
library(geosphere)

# plotting stuff
library(scales)

# misc
source('../R/helper01_process_foo.r')
source('../R/helper03_misc_foo.r')

# process the raw data into an anlyzable format

bin_width <- 1
age_max <- 63
form <- 'counti'

# occurrence information
neptune <- list.files(path = '../data', pattern = 'nannotax', full.names = TRUE) %>%
  purrr::map(., read_tsv) %>%
  purrr::map(., clean_names) %>%
  purrr::reduce(., rbind)

# map
sp <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=m")
globe.map <- readOGR('../data/ne_10m_coastline.shp')  # from natural earth
proj4string(globe.map) <- sp

# climate data
mgca <- read_tsv('../data/cramer/cramer_temp.txt') %>% 
  clean_names

# keep something constant across datasets
partial_raw_to_clean <- partial(raw_to_clean,
                                nano = neptune,
                                form = 'count',
                                bin_width = bin_width, 
                                age_max = age_max,
                                sp = sp,
                                mgca = mgca)
counti <- partial_raw_to_clean(restrict = FALSE)
counti_restrict <- partial_raw_to_clean(restrict = TRUE)

# write to file
write_rds(counti, path = '../data/counting.rds')
write_rds(counti_restrict, path = '../data/counting_restrict.rds')

