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

# survival
library(survival)
library(flexsurv)
library(survminer)

# plotting stuff
library(ggplot2)
library(scales)


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
  filter(!is.na(plat), !is.na(plng), fossil.group == 'F')

# assign million year bins (decreasing; youngest lowest)
nano <- nano %>%
  dplyr::arrange(desc(age)) %>%
  dplyr::mutate(mybin = ntile(age, 65))

# make a genus_species combo
nano <- nano %>% 
  dplyr::mutate(fullname = str_c(genus, '_', species)) %>% 
  arrange(fullname)



# join with metadata
# species ids
idinfo <- read_csv('../data/ezard2011/ezard_id2.csv')
names(idinfo) <- str_to_lower(names(idinfo))
names(idinfo) <- str_replace_all(names(idinfo), ' ', '.')
idinfo <- idinfo %>% 
  dplyr::mutate(fullname = str_replace_all(species.in.lineage, ' ', '_'),
                code = str_to_lower(lineage.code))

idinfo <- idinfo %>% separate_rows(fullname, sep = '-')

# combine with nano with left join
nano <- left_join(nano, idinfo, by = 'fullname')

# trait information
trait <- read_csv('../data/2010-09-06_aLext.csv')
names(trait) <- str_to_lower(names(trait))
names(trait) <- str_replace_all(names(trait), ' ', '.')

trait$ec <- fct_collapse(factor(trait$ec),
                         'mixed' = c(1, 2),
                         'thermocline' = '3',
                         'subthermocline' = '4')
trait <- trait %>%
  dplyr::mutate(label = str_to_lower(label),
                code = str_to_lower(nm),
                pn = str_to_lower(pn)) %>%
  separate_rows(code, sep = '-')


# combine all the data with an inner joing
nano <- inner_join(nano, trait, by = 'code')


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
                   eco = names(which.max(table(ec))),
                   morph = names(which.max(table(mp))),
                   keel = names(which.max(table(kl))),
                   symb = names(which.max(table(sy))),
                   spin = names(which.max(table(sp))))




# longitudinal dataset
# relative age in bins
longi <- sprange %>%
  dplyr::mutate(relage = abs(mybin - max(mybin))) %>%
  ungroup() %>%
  dplyr::mutate(id = as.numeric(as.factor(fullname)))

# make some data plots
gg <- unique(longi$fullname)[21:40]
lls <- longi[longi$fullname %in% gg, ]
glt <- ggplot(lls, aes(x = relage, y = maxgcd))
glt <- glt + geom_point() + geom_line()
glt <- glt + facet_wrap(~ fullname, ncol = 4)


# write to file
write_rds(longi, path = '../data/longitude.rds')




# survival dataset
# continuous data form
survi <- longi %>%
  group_by(fullname) %>%
  dplyr::summarise(duration = max(relage),
                   cohort = max(mybin),
                   dead = ifelse(min(mybin) == 1, 0, 1),
                   eco = unique(eco),
                   morph = unique(morph),
                   keel = unique(keel),
                   symb = unique(symb),
                   spin = unique(spin)) %>%
  ungroup() %>%
  dplyr::mutate(id = as.numeric(as.factor(fullname)))

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

counti$cc.rescale <- plyr::mapvalues(counti$cc, 
                                     from = sort(unique(counti$cc)), 
                                     to = seq(length(unique(counti$cc))))

sc <- with(counti, {survfit(Surv(time = time1, time2 = time2, 
                                 event = event, type = 'counting') ~ cc.rescale,
                            data = counti)})
gsc <- ggsurvplot(fit = sc, data = counti)


# make some data plots
sf <- with(survi, {survfit(Surv(duration, dead) ~ 1, survi)})
# survival plot (K-M est)
gsf <- ggsurvplot(fit = sf, data = survi, fun = 'pct')  
# event plot
gsfe <- ggsurvevents(fit = sf, data = survi)  

# write to file
write_rds(survi, path = '../data/survival.rds')
write_rds(counti, path = '../data/counting.rds')
