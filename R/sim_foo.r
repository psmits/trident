#' Simulate a species geographic range trajectory
#'
#' This is a bounded random walk. Starts at 0, ends when back to 0. Uses discrete time steps.
#'
#' Simulation can end before the maximum number of discete time steps.
#' 
#' Change to value is determined by draws from a standard normal distribution.
#'
#' @param max_length scalar int for maximum number of discrete time steps.
#' @return vector of of real values
sim_range <- function(max_length = 10) {

  # this is easily solved by a fold across a series 
  # that series begins at 0, has max_length changes
  # accumulate changes
  run <- accumulate(c(0, rnorm(max_length)), `+`)

  # first time we hit negative is extinction
  # this could be right away (index == 2)
  death <- run %>%
    detect_index(~ . < 0)
  if(death == 0) {
    out <- run 
  } else if(death > 0) {
    out <- run[1:(death - 1)]
  }

  out
}

