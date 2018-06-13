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

#' Measure maximum great circle from a geocoordinates
#'
#' Uses functions from geosphere to measure maximum great circle distance amoungst a cloud of points. 
#' Coordinates are in longitude, and latitude. Output is in meters. 
#'
#' @param x vector of longitudes
#' @param y vector of latitudes
#' @return maximum greater distance measured in meters
dist_gcd <- function(x, y) max(distm(cbind(x, y), fun = distGeo))

#' Get most common entry of a vector of characters or factors
#' 
#' Given a vector of characters or factors, count number of occurrences of each unique entry.
#' Return most common (using which.max rules).
#' 
#' @param x vector of characters or factors
#' @return most common entry in vector x
plurality <- function(x) names(which.max(table(x)))




#' Prepare tibble for analysis
#'
#' This is mostly an internal function. Takes a tibble and returns a tibble with some variables transformed.
#'
#' @param x tibble
#' @param fossil.group vector of characters defining which fossil groups to include
#' @return a tibble ready for analysis
prepare_analysis <- function(x, fossil.group = NULL) {
  # do i want to restrict the taxonomic group?
  tb <- x

  if(!is.null(fossil.group)) {
    tb <- tb %>%
      filter(fossil.group %in% fossil.group)
  }

  # prep the data
  tb <- tb %>%
    arrange(fullname, desc(mybin)) %>%
    mutate_at(.vars = vars(ncell, latext, maxgcd, area, temp),
              .funs = ~ arm::rescale(log1p(.x))) %>%
    group_by(fullname) %>%
    mutate(diff_maxgcd = c(0, diff(maxgcd)),
           lag1_maxgcd = lag(maxgcd, default = 0),
           lag2_maxgcd = lag(maxgcd, n = 2, default = 0),
           lead1_maxgcd = lead(maxgcd, default = 0),
           lead2_maxgcd = lead(maxgcd, n = 2, default = 0),
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
#' Have fun.
#' 
#' @param x vector of ages
#' @param by bin width
#' @return vector of bin memberships
break_my <- function(x, by = 1, number = NULL) {
  top <- ceiling(max(x))
  bot <- floor(min(x))
  if(!is.null(by)) {
    unt <- seq(from = bot, to = top, by = by)
  } else {
    unt <- seq(from = bot, to = top, length.out = number)
  }
  unt1 <- unt[-length(unt)]
  unt2 <- unt[-1]
  uu <- map(map2(unt1, unt2, ~ x >= .x & x < .y), which)

  y <- x
  for(ii in seq(length(uu))) {
    y[uu[[ii]]] <- ii
  }
  y
}
