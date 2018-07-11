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
break_my <- function(x, by = NULL, number = NULL) {
  top <- ceiling(max(x))
  bot <- floor(min(x))
  if(!is.null(by)) {
    unt <- seq(from = bot, to = top, by = by)
  } else if(!is.null(number)) {
    unt <- seq(from = bot, to = top, length.out = number + 1)
  }
  unt1 <- unt[-length(unt)]
  unt2 <- unt[-1]
  uu <- map2(unt1, unt2, ~ which(between(x, left = .x, right = .y)))

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



#' Get optimal cut from ROC results
#' 
#' Uses ROCR. I got this from somewhere on the internet.
#'
#' @param perf performance object
#' @param pred prediction object
#' @return vector with sensitivity, specificity, and cutpoint
opt_cut <- function(perf, pred) {
  cut_ind <- mapply(FUN = function(x, y, p) {
                      d <- (x - 0)^2 + (y - 1)^2
                      ind <- which(d == min(d))
                      c(sensitivity = y[[ind]], 
                        specificty = 1 - x[[ind]],
                        cutoff = p[[ind]])
               }, perf@x.values, perf@y.values, pred@cutoffs)
  cut_ind
}

#' Get optimal cutpoints for a lot of ROC results
#'
#' Looped wrapper around opt_cut to get optimal cutpoint for a given posterior predictive distribution
#'
#' @param pp_prob list of posterior predictive draws, where each list element is the PPD for one model.
#' @param target the actual values that are being estimated
get_cutpoints <- function(pp_prob, target) {
  pred <- plyr::alply(pp_prob, 1, function(x) prediction(x, target))
  perf <- map(pred, ~ performance(.x, measure = 'tpr', x.measure = 'fpr'))
  cutpoints <- map2(perf, pred, opt_cut)
}


#' Reassign class membership based on new cutpoints
#'
#' Give probabilities and new cutpoint level, reassigns classes.
#'
#' @param pp_prob list of posterior predictive draws, where each list element is the PPD for one model.
#' @param list_cutpoint list of new cutpoints as estimated from get_cutpoints
#' @return list of reassigned class predictions
cut_newpoint <- function(pp_prob, list_cutpoint) {
  pp_est_new <- list()
  for(jj in seq(length(pp_prob))) {
    mm <- matrix(ncol = ncol(pp_prob[[jj]]), nrow = nrow(pp_prob[[jj]]))
    for(ii in seq(nrow(mm))) {
      mm[ii, ] <- pp_prob[[jj]][ii, ] > list_cutpoint[[jj]]
    }
    pp_est_new[[jj]] <- mm * 1
  }
  pp_est_new
}

#' Posterior roc
#' 
#' This is mostly a convenient function to quickly get ROC estimates for entire posterior distribution (awkward matrix form).
#' 
#' @param x matrix of posterior 
#' @param y matrix with true data. named element event is vector of 0/1
post_roc <- function(x, y) apply(x, 1, function(a) roc(y$event, a))

#' Extract posterior AUC, safely
#'
#' Wrapped with safely to prevent response has only 1 level stuff
#'
#' @param pred matrix of predicted 0/1 values
#' @param counti_fold testing dataset
#' @return list w/ two elements result and error
get_auc_time <- function(pred, counti_fold) {
  fold_time <- split(counti_fold, counti_fold$mybin)

  pred_time <- split(t(pred), counti_fold$mybin) %>%
    map2(., fold_time, ~ matrix(.x, nrow = nrow(.y))) %>%
    map(., t)

  safe_pr <- safely(post_roc)
  auc_time <- map2(pred_time, fold_time, safe_pr) %>%
    map(., function(x) map_dbl(x$result, ~ auc(.x)))
  auc_time
}



