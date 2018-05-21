#' Species' geographic range over time, compared
#' 
#' Make a plot comparing species geographic ranges over time, where all species are scaled to duration 1.
#' This inherently means species with duration of 1 time unit cannot be plotted.
#' 
#' @param data tibble of longitudinal data
#' @return ggplot object
plot_georange_compare <- function(data) {
  scale_01 <- function(x) (x - min(x)) / (max(x) - min(x))
  
  data %>%
    group_by(fullname) %>%
    mutate(scaleage = scale_01(mybin)) %>%
    ggplot(., aes(x = scaleage, y = log(maxgcd), group = fullname)) + 
    geom_line(alpha = 0.1)

}


