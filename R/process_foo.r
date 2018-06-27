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
              .funs = ~ arm::rescale(log1p(.x))) %>%
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

#' Did function return an error?
#' 
#' Given result, did it return a try-error? Assumes results wrapped in try(...). 
#' This is used with purrr::keep to process vectors.
#' TRUE means there is no error because we're checking for clean. 
#' FALSE means an error was returned.
#' This is a clean-ish way of using try(...) when you didn't or can't use purrr::safely.
#'
#' @param x element.
#' @return logical TRUE for no error, FALSE for error.
check_class <- function(x) class(x) != 'try-error'
