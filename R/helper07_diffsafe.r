#' Safe differencing with uncertain vector length
#'
#' What it says on the tin.
#'
#' @param x numerical vector
#' @param lag difference lag
#' @return vector of differences of appropriate length for tibble
diff_safe <- function(x, lag) {
  n <- length(x)
  if(n < lag) {
    out <- rep(0, n)
  } else {
    out <- c(rep(0, lag), diff(x, lag = lag))
  }
  out
}


