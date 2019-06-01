#' Sample N times within groups for dplyr
#'
#' Solution is from https://github.com/tidyverse/dplyr/issues/361.
#' Solution written by Kendon Bell
#'
#' @param tbl a tbl that is already grouped
#' @param size an integer scalar indicating the number of samples per group
#' @param replace sample with replacement? logical
#' @param weight weighting of samples? default NULL/all equal
#' @return tbl with each group subsampled to size 
sample_n_groups = function(tbl, size, replace = FALSE, weight = NULL) {

  # regroup when done
  grps = tbl %>% groups %>% lapply(as.character) %>% unlist
  # check length of groups non-zero
  keep = tbl %>% summarise() %>% ungroup() %>% sample_n(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping variable
  tbl %>% right_join(keep, by=grps) %>% group_by_(.dots = grps)
}

#' Prepare tibble for analysis
#'
#' This is mostly an internal function. Takes a tibble and returns a tibble with some variables transformed.
#'
#' @param x tibble
#' @param fg vector of characters defining which fossil_groups to include
#' @return a tibble ready for analysis
prepare_analysis <- function(x, fg = NULL) {
  # do i want to restrict the taxonomic group?
  tb <- x

  if(!is.null(fg)) {
    tb <- tb %>%
      filter(fossil_group %in% fg)
  }

  # prep the data
  tb <- tb %>%
    arrange(fullname, desc(mybin)) %>%
    mutate_at(.vars = vars(ncell, latext, maxgcd, area, temp),
              .funs = ~ scale(log1p(.x))) %>%
    group_by(fullname) %>%
    mutate(diff_maxgcd = c(0, diff(maxgcd)),
           lag1_maxgcd = lag(maxgcd, default = 0),
           lag2_maxgcd = lag(maxgcd, n = 2, default = 0),
           lead1_maxgcd = lead(maxgcd, default = 0),
           lead2_maxgcd = lead(maxgcd, n = 2, default = 0),
           diff_latext = c(0, diff(latext)),
           lag1_latext = lag(latext, default = 0),
           lag2_latext = lag(latext, n = 2, default = 0),
           lead1_latext = lead(latext, default = 0),
           lead2_latext = lead(latext, n = 2, default = 0),
           diff_area = c(0, diff(area)),
           lag1_area = lag(area, default = 0),
           lag2_area = lag(area, n = 2, default = 0),
           lead1_area = lead(area, default = 0),
           lead2_area = lead(area, n = 2, default = 0),
           diff_temp = c(0, diff(temp)),
           lag1_temp = lag(temp, default = 0),
           lag2_temp = lag(temp, n = 2, default = 0),
           lead1_temp = lead(temp, default = 0),
           lead2_temp = lead(temp, n = 2, default = 0)) %>%
    ungroup() %>%
    mutate(fact_mybin = as.factor(mybin),
           fact_relage = as.factor(relage))

    tb$fact_mybin <- as.factor(tb$mybin)
    tb$fact_relage <- as.factor(tb$relage)

    tb
}


#' Break time data up into bins
#' 
#' Have fun with this. basic rules. greater than equal to base, less than top.
#' 
#' @param x vector of ages
#' @param by bin width
#' @param age logical bin age returned, not number (default FALSE, return bin number)
#' @return vector of bin memberships
#' @author Peter D Smits <peterdavidsmits@gmail.com>
break_my <- function(x, by = NULL, number = NULL, age = FALSE) {

  if(is.null(by) & is.null(number)) {
    return('no scheme given. specify either bin width or number of bins.')
  }

  if(!is.null(by) & !is.null(number)) {
    return('too much information. specify either bin width OR number of bins, not both.')
  }

  # range to bin
  top <- ceiling(max(x))
  bot <- floor(min(x))

  # create bins
  if(!is.null(by)) {
    unt <- seq(from = bot, to = top, by = by)
  } else if(!is.null(number)) {
    unt <- seq(from = bot, to = top, length.out = number + 1)
  }

  # bin top and bottom
  unt1 <- unt[-length(unt)]
  unt2 <- unt[-1]

  # assign memberships
  uu <- map2(unt1, unt2, ~ which(between(x, left = .x, right = .y)))

  # what if we want the "age" of the bin, not just number?
  if(age == TRUE) {
    unt_age <- map2_dbl(unt1, unt2, ~ median(c(.x, .y)))
  }

  # create output vector
  y <- x
  for(ii in seq(length(uu))) {
    if(age == FALSE) {
      y[uu[[ii]]] <- ii
    } else if(age == TRUE) {
      y[uu[[ii]]] <- unt_age[ii]
    }
  }
  y
}



#' Translate raw neptune reads into analyzable forms
#'
#' Using neptune database, a map, a climate record, and some information, give the counting form dataframe for further analysis.
#'
#' @param nano tibble of neptune data
#' @param form output format (count is counting format). can also longi(tude) (default count)
#' @param bin_width numeric scalar of how many million years each bin is (default 1 my)
#' @param age_max numeric scalar oldest fossil occurrences allowed
#' @param restrict logical scalar if restricting the duration of individual species (default false)
#' @param .width length 2 vector with upper and lower percentage for truncating temporal ranges (default 0.01 and 0.99)
#' @param sp sp object of globe
#' @param mgca tibble of temperature estiamtes
#' @param prov character vector of Longhurst codes to KEEP in the dataset (default NULL)
#' @return a tibble
raw_to_clean <- function(nano, 
                         form = 'count',
                         bin_width = 1, 
                         age_max, 
                         restrict = FALSE, 
                         .width = c(0.01, 0.99),
                         sp,
                         mgca,
                         prov = NULL) {

  #bin_number <- age_max / bin_width

  # important but impossible to call variable name
  names(nano)[19] <- 'age'
  names(nano)[23] <- 'plat'
  names(nano)[24] <- 'plng'

  # one of the fossil group tibble's miscoded fossil group as logical
  nano[nano$fossil_group == 'FALSE', 'fossil_group'] <- 'F'

  nano <- nano %>%
    filter(#fossil_group == 'F',
           fossil_group != 'DN')

  # filter and give ages
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
    dplyr::mutate(mybin = break_my(age, by = bin_width)) %>%
    # then
    # make a species genus combo
    dplyr::mutate(fullname = str_c(genus, '_', species)) %>% 
    arrange(fullname)


  # restrict occurrences to core sequence
  if(restrict == TRUE) {
    nano <- nano %>%
      group_by(fullname) %>%
      mutate(lower = quantile(age, probs = min(.width)),
             upper = quantile(age, probs = max(.width))) %>%
      filter(age > lower,
             age < upper) %>%
      ungroup()
  }

  # remove points not in selected region
  if(!is.null(prov)) {
    nano <- nano %>%
      filter(longhurst_code %in% prov)
  }


  # assign everything a geographic cell using paleocoordinates
  spatialref <- SpatialPoints(coords = nano[, c('plng', 'plat')],
                              proj4string = sp)  # wgs1984.proj
  r <- raster(globe.map, nrows = 50, ncols = 50)
  ras <- rasterize(spatialref, r)
  # get cell # for each observation
  nano$cell <- cellFromXY(ras, 
                          xy = as.data.frame(nano[, c('plng', 'plat')]))
  
  
  # get all the important geographic range information
  # have to summarize the big matrix nano
  sprange <- 
    nano %>%
    group_by(fullname, mybin) %>%
    dplyr::summarise(nocc = n(),  # number of occurrences in bin
                     ncell = n_distinct(cell),  # number of unique cells in bin
                     # latitidinal extent in bin
                     latext = abs(max(plat) - min(plat)),  
                     # great circle distance (on ellipsoid) in km in bin
                     area = areaPolygon(cbind(plng, plat)),
                     maxgcd = dist_gcd(plng, plat),
                     nprov = n_distinct(longhurst_code),
                     # NEW minimum spanning tree, only need measure in km
                     mst = MSTDist(longs = plng, lats = plat)$MST_km,
                     fossil_group = plurality(fossil_group)) %>%
    filter(maxgcd > 0,
           latext > 0)


  # want to add in the lear mgca data
  mgca <- mgca %>%
    mutate(mybin = break_my(age, by = bin_width)) %>%
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
  sprange <- split(sprange, sprange$fullname) %>%
    .[map_lgl(., ~ max(.x$mybin) > 1)] %>%
    reduce(., rbind)


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

  if(form == 'long') {
    return(longi) 
  } else if(form == 'count') {
    return(counti)
  }
}
