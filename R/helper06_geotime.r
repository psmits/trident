#a goal is grey rectangles

#' Highlight timeseries using geological time scale
get_geotime_box <- function(time_range, type = 'age', alternate = 'even') {
  gts <- read_csv(here('data', 'geologic_time_scale.csv'))

  lts <- 
    gts %>%
    filter(type == !!type,
           end >= min(time_range),
           start <= max(time_range)) %>%
    rowwise() %>%
    mutate(xmin = max(start, min(time_range)),
           xmax = min(end, max(time_range))) %>%
    ungroup() %>%
    mutate(ymin = -Inf,
           ymax = Inf)
  
  # how are the bands alternating?
  rr <- seq(nrow(lts))
  if(alternate == 'even') {
    lts <- lts[(rr %% 2) == 0, ]
  } else if(alternate == 'odd') {
    lts <- lts[(rr %% 2) == 1, ]
  }

  lts
}
