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
source('../R/process_foo.r')


# constants
bin_width <- 1
age_max <- 63
bin_number <- age_max / bin_width

# read in data file and clean column names
neptune <- list.files(path = '../data', pattern = 'nannotax', full.names = TRUE)
nano <- purrr::map(neptune, read_tsv)
nano <- purrr::map(nano, clean_names)
nano <- purrr::reduce(nano, rbind)

# important but impossible to call variable name
names(nano)[19] <- 'age'
names(nano)[23] <- 'plat'
names(nano)[24] <- 'plng'

# one of the fossil group tibble's miscoded fossil group as logical
nano[nano$fossil_group == 'FALSE', 'fossil_group'] <- 'F'

# miscoded as strings
nano <- nano %>%
  mutate_at(c('longitude', 'latitude'), .funs = parse_double)

nano <- nano %>%
  filter(#fossil_group == 'F',
         fossil_group != 'DN')

nano <- nano %>% 
  filter(!is.na(plat), 
         !is.na(plng)) %>%
  group_by(taxon_id) %>%
  filter(max(age) < age_max) %>%
  ungroup() %>%
  # some don't have paleolat or long
  # also need to limit to line up with mgca
  # then
  # assign million year bins (decreasing; youngest lowest)
  dplyr::arrange(desc(age)) %>%
  dplyr::mutate(mybin = break_my(age, by = 1)) %>%
  # then
  # make a species genus combo
  dplyr::mutate(fullname = str_c(genus, '_', species)) %>% 
  arrange(fullname)

# assign everything a geographic cell using paleocoordinates
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

# ecological information
idinfo <- read_csv('../data/ezard2011/ezard_id2.csv')
idinfo <- idinfo %>% clean_names
idinfo <- idinfo %>% 
  dplyr::mutate(fullname = str_replace_all(species_in_lineage, ' ', '_'),
                code = str_to_lower(lineage_code))
idinfo <- idinfo %>% separate_rows(fullname, sep = '-')
# combine with nano with left join

# trait information
trait <- read_csv('../data/2010-09-06_aLext.csv')
trait <- trait %>% clean_names
trait$ec <- fct_collapse(factor(trait$ec),
                         'mixed' = c(1, 2),
                         'thermocline' = '3',
                         'subthermocline' = '4')
trait <- trait %>%
  dplyr::mutate(label = str_to_lower(label),
                code = str_to_lower(nm),
                pn = str_to_lower(pn)) %>%
separate_rows(code, sep = '-')

# get the most core version of the ecological information
eco_foram <- inner_join(idinfo, trait, by = 'code')
# haven't combined with microfossil data yet
# left join means i keep all rows but add info to where exists
nano <- 
  left_join(nano, eco_foram, by = 'fullname')


# get all the important geographic range information
# have to summarize the big matrix nano
sprange <- nano %>%
  group_by(fullname, mybin) %>%
  dplyr::summarise(nocc = n(),  # number of occurrences in bin
                   ncell = n_distinct(cell),  # number of unique cells in bin
                   # latitidinal extent in bin
                   latext = abs(max(plat) - min(plat)),  
                   # great circle distance (on ellipsoid) in km in bin
                   area = areaPolygon(cbind(plng, plat)),
                   maxgcd = dist_gcd(plng, plat),
                   nprov = n_distinct(longhurst_code),
                   fossil_group = plurality(fossil_group)) %>%
  filter(maxgcd > 0,
         latext > 0)


# want to add in the lear mgca data
mgca <- read_tsv('../data/cramer/cramer_temp.txt')
mgca <- mgca %>% clean_names
mgca <- mgca %>%
  mutate(mybin = break_my(age)) %>%
  group_by(mybin) %>%
  summarize(temp = mean(temperature, na.rm = TRUE),
            temp_sd = sd(temperature, na.rm = TRUE)) %>%
  ungroup()

# put the climate data into the dataframe
# easy because i've binned them the say way using break_my
sprange <- left_join(sprange, mgca, by = 'mybin')

# exclude the last cohort because artifact re. never possible to die
# another example of me not knowing the tidy solution
# because i'm not excluding rows within groups
# i'm removing groups based on row value
r <- split(sprange, sprange$fullname)
r <- r[map_lgl(r, ~ max(.x$mybin) > 1)]
sprange <- reduce(r, rbind)

# longitudinal dataset
# relative age in bins
longi <- sprange %>%
  dplyr::mutate(relage = abs(mybin - max(mybin))) %>%
  ungroup() %>%
  dplyr::mutate(id = as.factor(fullname))
# longi has a weird data sorting issue i need to figure out
# so this is currently "inelegant" (not tidy)
ff <- split(longi, longi$fullname)
ff <- purrr::map(ff, function(x) {
                   x <- x[order(x$relage), ]
                   x})
longi <- purrr::reduce(ff, bind_rows)


# write to file
write_rds(longi, path = '../data/longitude.rds')


# survival dataset
# continuous data form
survi <- longi %>%
  group_by(fullname) %>%
  dplyr::summarise(duration = max(relage),
                   cohort = max(mybin),
                   dead = ifelse(min(mybin) == 1, 0, 1)) %>%
  #eco = unique(eco),
  #morph = unique(morph),
  #keel = unique(keel),
  #symb = unique(symb),
  #spin = unique(spin)) %>%
  ungroup() %>%
  dplyr::mutate(id = as.factor(fullname))

survi <- survi %>%
  dplyr::mutate(cohort = as.character(cohort),
                cc = fct_drop(cohort))
survi$cc.rescale <- plyr::mapvalues(survi$cc, 
                                    from = sort(unique(survi$cc)), 
                                    to = seq(length(unique(survi$cc))))


# counting process form
counti <- longi %>%
  group_by(fullname) %>%
  dplyr::mutate(time1 = relage,
                time2 = relage + 1,
                event = ifelse(relage == max(relage) &
                               min(mybin) != 1, 1, 0),
                cohort = as.character(max(mybin)),
                cc = fct_drop(cohort)) %>%
  ungroup()


# gaps, names, etc.
counti$cc.rescale <- plyr::mapvalues(counti$cc, 
                                     from = sort(unique(counti$cc)), 
                                     to = seq(length(unique(counti$cc))))
# necessary data transform to handle factor
counti$cc.rescale <- with(counti, {
                            factor(cc.rescale, 
                                   levels = sort(unique(as.numeric(cc.rescale))))})

# write to file
write_rds(survi, path = '../data/survival.rds')
write_rds(counti, path = '../data/counting.rds')
