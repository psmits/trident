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


