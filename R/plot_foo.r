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

#' Taxonomic group by Covariate over time plot
#'
#' Given an rstanarm object conforming to this analysis, produce a ggplot object where each coviariate is plotted over time for each of the observed taxonomic groups.
#'
#' @param disc_best rstanarm object
#' @return a ggplot object
plot_taxon_covariate_time <- function(disc_best) {
  comp_const <- disc_best %>%
    spread_samples(maxgcd, diff_maxgcd)

  comp_var <- disc_best %>%
    spread_samples(b[i, f]) %>%
    spread(i, b) %>%
    filter(str_detect(f, 'fact_mybin')) %>%
    mutate(type = str_remove_all(f, '[0-9]'),
           age = as.numeric(str_extract(f, '[0-9]+'))) %>%
    ungroup

  # need to confirm correct line-up
  cv <- split(comp_var, comp_var$type) %>%
    map(., ~ .x %>% arrange(.chain, .iteration, age))

  # first element is just temporal effect
  core <- full_join(comp_const, cv[[1]], by = c('.chain', '.iteration')) %>%
    mutate(eff_maxgcd = maxgcd.x + maxgcd.y,
           eff_diff_maxgcd = diff_maxgcd.x + diff_maxgcd.y) %>%
  dplyr::select(.chain, .iteration, f, age, 
                eff_maxgcd, eff_diff_maxgcd) %>%
  arrange(.chain, .iteration, age)

  # need to confirm correct line-up
  by_taxon <- map(cv[-1], ~ full_join(core, .x, by = c('.chain', '.iteration', 'age'))) %>%
    map(., ~ .x %>%
        mutate(taxon_eff_maxgcd = eff_maxgcd + maxgcd,
               taxon_eff_diff_maxgcd = eff_diff_maxgcd + diff_maxgcd))
  by_taxon <- reduce(by_taxon, bind_rows)
  
  by_taxon <- by_taxon %>%
    gather(key, value, -.chain, -.iteration, -f.x, 
           -age, -f.y, -`(Intercept)`, -type,
           -eff_maxgcd, -eff_diff_maxgcd,
           -diff_maxgcd, -maxgcd) %>%
  filter(!is.na(type)) %>%
  ggplot(aes(x = age, y = value)) +
  stat_lineribbon() +
  scale_fill_brewer() +
  facet_grid(key ~ type, scales = 'free_y')
  
  by_taxon
}


#' Plot risk of extinction for selection of species over time
#'
#' Generate a plot of species risk of extinction over duration at each observation.
#' Risk of extinction is calculated given all covariate information for a species over time.
#' 
#' @param data tibble to make plots from
#' @param model rstanarm object to calculated from
#' @param nsp number of species to grab
#' @return ggplot object
plot_risk_time <- function(data, model, nsp = 4) {
  temp <- data %>%
    group_by(fullname) %>%
    sample_n_groups(size = nsp)

  # estimate log-odds extinction for new species at each time
  temp_est <- posterior_linpred(object = model, 
                                newdata = temp, 
                                transform = TRUE)
  temp_est <- reshape2::melt(temp_est)
  names(temp_est) <- c('iterations', 'row', 'value')
  temp_est <- as.tibble(temp_est) %>% 
    mutate(fullname = row,
           relage = row)
  temp_est$fullname <- plyr::mapvalues(temp_est$fullname, 
                                       from = unique(temp_est$fullname), 
                                       to = temp$fullname)
  temp_est$relage <- plyr::mapvalues(temp_est$relage, 
                                     from = unique(temp_est$relage), 
                                     to = temp$relage)

  # now we can combine the empirical data with the linpred values
  full_est <- temp %>%
    dplyr::select(fullname, maxgcd, relage) %>% 
    left_join(., temp_est, by = c('fullname', 'relage'))

  ext_plot <- full_est %>%
    ggplot(aes(x = relage, y = value)) + 
    stat_pointinterval(.prob = c(0.5, 0.8)) +
    facet_grid(~ fullname) + 
    labs(x = 'Age (My)', y = 'P(T = t | T >= t, x)')

  range_plot <- full_est %>%
    group_by(fullname, relage) %>%
    summarise(maxgcd = mean(maxgcd)) %>%
    ggplot(aes(x = relage, y = maxgcd)) + 
    geom_line() +
    facet_grid(~ fullname) + 
    labs(x = 'Age (My)', y = 'P(T = t | T >= t, x)')

  full_plot <- full_est %>%
    gather(key, value, -fullname, -iterations, -relage, -row) %>%
    ggplot(aes(x = relage, y = value)) +
    geom_point(alpha = 0.01) + 
    facet_grid(key ~ fullname, scales = 'free_y') +
    labs(x = 'Age (My)', y = 'P(T = t | T >= t, x)')

  out <- list(ext_plot, range_plot, full_plot)
  out
}



#' Faceted plot of baseline hazard by fossil group
#'
#' Given model fit, produce a faceted plot of discrete time hazard. 
#' Hazard is presented as a probability of going extinct. 
#' 
#' @param model rstanarm object
#' @return ggplot2 object
plot_taxon_hazard <- function(model) {
  # the average
  haz_avg <- model %>%
    spread_samples(b[i, f], `(Intercept)`) %>%
    filter(str_detect(f, pattern = 'fact_relage'),
           !str_detect(f, pattern = 'fossil_group')) %>%
    mutate(cr = `(Intercept)` + b) %>%    # keep log-odds scale
    mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
      arrange(age)

  # components from fossil_groups
  db <- model %>%
    spread_samples(b[i, f]) %>%
    filter(str_detect(f, pattern = 'fact_relage'),
           str_detect(f, pattern = 'fossil_group')) %>%
    mutate(type = str_remove_all(f, '[0-9]'),
           age = as.numeric(str_extract(f, '[0-9]+'))) %>%
    arrange(age) %>%
    split(., .$type) %>%
    map(., ~ left_join(.x, haz_avg, by = c('.chain', '.iteration', 'age'))) %>%
    map(., ~ .x %>% 
        ungroup() %>%
        mutate(effect = b.x + cr) %>%
        dplyr::select(-.chain, -.iteration,
                      -i.x, -f.x, -b.x, 
                      -i.y, -f.y, -b.y, 
                      -`(Intercept)`, -cr) %>%
        mutate(type = str_extract(type, pattern = '[A-Z]'))) %>%
    bind_rows %>%
    mutate(effect_prob = invlogit(effect)) %>%  # put on prob scale
    ggplot(aes(x = age, y = effect_prob)) + 
      stat_lineribbon() +
      scale_fill_brewer() +
      facet_grid(type ~ .) +
      labs(x = 'Age (My)', y = 'P(T = t | T >= t, x)')

  db
}


#' Faceted plot of AUC time series
#' 
#' This spews a bunch of error messages by they are handeled through exception handling, so the error messages are actually warnings :P.
#' 
#' @param data tibble of data used to fit models
#' @param model_pp list of model posterior predictive draws
#' @param model_key vector characters naming each model
#' @return ggplot object
plot_roc_series <- function(data, model_pp, model_key) {
#data <- counti_trans
#model_pp <- pp_est_new
#model_key <- model_key
  # important exception functions
  safe_roc <- safely(roc)
  safe_auc <- safely(auc)
  
  out <- list()
  for(kk in seq(length(model_pp))) {
    est <- list()
    for(ii in seq(max(data$mybin))) {
      tw <- model_pp[[kk]][, data$fact_mybin == ii]
      ew <- data %>%
        filter(fact_mybin == ii) %>%
        dplyr::select(event)
      ew <- as.data.frame(ew)[, 1]
      oo <- list()
      for(jj in seq(nrow(tw))) {
        # see here for the big source of exceptions...
        oo[[jj]] <- safe_roc(ew, tw[jj, ])
      }
      est[[ii]] <- oo
    }
    est <- set_names(est, seq(length(est)))
    out[[kk]] <- est
  }
  # of those that aren't errors, try to get AUC
  es <- map(out, ~ map(.x, ~ map(.x, ~ safe_auc(.x$result)))) %>%
    map(., ~ map(.x, ~ reduce(map(.x, 'result'), c)))

  # zero out the error-d entries
  tt <- map(es, ~ map_lgl(.x, ~ !is.null(.x)))
  slate <- map2(es, tt, ~ keep(.x, .y))

  # massage into graph
  roc_ts <- 
    map(slate, bind_rows) %>%
    map(., ~ .x %>% gather(key, value)) %>%
    bind_rows(., .id = 'model') %>%
    mutate(key = parse_integer(key),
           model = plyr::mapvalues(model, 
                                   from = seq(length(model_key)), 
                                   to = model_key)) %>%
    ggplot(aes(x = key, y = value)) +
    stat_lineribbon() +
    scale_fill_brewer() +
    facet_grid(model ~ .) +
    labs(x = 'Time (My)', y = 'AUC') +
    NULL
  roc_ts
}
