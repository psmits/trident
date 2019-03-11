#' calculate in sample roc for taxon X time
#'
#' what it says on the tin
#'
#' @param fit_list list of model fits
#' @param .data tibble of neptune data
#' @param key vector of model names
#' @return a tibble
ins_roc_taxon_time <- function(fit_list, .data, key) {
  # posterior predictive distribution for ins_fit
  ppd <- future_map(fit_list, 
                    ~ posterior_linpred(.x, transform = TRUE, draws = 100))

  # setup breaking by taxonXtime
  bysplit <- .data %>%
    dplyr::select(mybin, fossil_group, event) %>%
    mutate(type = paste0(fossil_group, ':', mybin))

  bb <- split(bysplit, bysplit$type)

  # break by taxonXtime

  ppd_split <- map(ppd, ~ as_tibble(t(.x)) %>%
                   split(., bysplit$type))

  # calculate roc
  safe_roc <- safely(roc)
  safe_auc <- safely(auc)
  foo <- function(event, prob) {
    safe_auc(safe_roc(event, prob)$result)
  }
  split_auc <- map(ppd_split, function(x) 
                   map2(x, bb, ~ apply(.x, 2, function(a)  
                                       foo(.y$event, a)))) %>%
  map(., ~ map(.x, ~ reduce(map(.x, 'result'), c)))

auc_taxon_time <- map(split_auc, function(x) map(x, ~ as_tibble(x = .x)))%>%
  map(., ~ bind_rows(.x, .id = 'phyla_time')) %>%
  bind_rows(., .id = 'model') %>%
  separate(., col = phyla_time, into = c('phyla', 'time'), sep = ':') %>%
  mutate(time = as.numeric(time),
         fossil_group = case_when(phyla == 'D' ~ 'Dinoflagellates',
                                  phyla == 'R' ~ 'Radiolaria',
                                  phyla == 'F' ~ 'Foraminifera',
                                  phyla == 'N' ~ 'Calc. nanno.'),
         model = case_when(model == 1 ~ key[1],
                           model == 2 ~ key[2],
                           model == 3 ~ key[3],
                           model == 4 ~ key[4]),
         model = factor(model, levels = rev(key)))

  auc_taxon_time
}


#' calculate in sample roc for taxon X time
#'
#' what it says on the tin
#'
#' @param fit list of model fits
#' @param .data tibble of neptune data
#' @param key vector of model names
#' @return a tibble
oos_roc_taxon_time <- function(fit, .data, key) {

  # predict the test data
  pred <- get_pred(fit, .data)

  # break up the posterior predictive estimates by taxon and time
  # break up the data to match nicely
  
  # prepare the data in format
  bysplit <- .data %>%
    map(., 
        ~ .x %>%
          dplyr::select(mybin, fossil_group, event) %>%
          mutate(type = paste0(fossil_group, ':', mybin)))

  bb <- bysplit %>%
    map(., ~ split(.x, .x$type))         # by taxon and time
  
  # break up probs
  tt <- map(pred, ~ as_tibble(t(.x))) %>%
    map2(., bysplit, ~ split(.x, .y$type))  
  
  
  # calculate auc for each taxonXtime combo
  safe_roc <- safely(roc)
  safe_auc <- safely(auc)
  foo <- function(event, prob) {
    safe_auc(safe_roc(event, prob)$result)
  }
  split_auc <- map2(tt, bb, function(xx, yy)
                    map2(xx, yy, ~ apply(.x, 2, function(aa)
                                         foo(.y$event, aa)))) %>%
    map(., ~ map(.x, ~ reduce(map(.x, 'result'), c)))
  
  # recombine
  fold_auc_taxon_time <- map(split_auc, function(x) map(x, ~ as_tibble(x = .x))) %>%
    map(., ~ bind_rows(.x, .id = 'phyla_time')) %>%
    bind_rows(., .id = 'model') %>%
    separate(., col = phyla_time, into = c('phyla', 'time'), sep = ':') %>%
    separate(., col = model, into = c('model', 'fold'), sep = '_') %>%
    mutate(time = as.numeric(time),
           fossil_group = case_when(phyla == 'D' ~ 'Dinoflagellates',
                                    phyla == 'R' ~ 'Radiolaria',
                                    phyla == 'F' ~ 'Foraminifera',
                                    phyla == 'N' ~ 'Calc. nanno.'),
           model = case_when(model == 'mod1' ~ key[1],
                             model == 'mod2' ~ key[2],
                             model == 'mod3' ~ key[3],
                             model == 'mod4' ~ key[4]),
           model = factor(model, levels = rev(key)))
 
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


