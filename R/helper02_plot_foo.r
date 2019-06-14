#' Plot of stacked Bayes R2 values
#' 
#' Given a list of posterior fits, computes bayes r2. returns a plot
#' @param data list of tibbles
#' @param model_key model types
#' @return ggplot object
plot_bayesr2 <- function(data, model_key) {
  br2 <- map(data, bayes_R2)
  br2_gg <- reshape2::melt(br2) %>%
    rename(model = L1) %>%
    mutate(model_name = plyr::mapvalues(model, unique(model), model_key),
           model_name = factor(model_name, levels = model_key)) %>%
    ggplot(aes(x = value, y = model_name)) +
    geom_halfeyeh(width = c(0.5, 0.8)) +
    labs(x = 'Bayesian R^2', y = 'Model')
  br2_gg
}


#' Plot of stacked Bayes R2 values
#' 
#' Given a list of posterior fits, computes bayes r2. returns a plot
#' @param data an appropriate tibble
#' @param model_key model types
#' @return ggplot object
plot_roc_curve <- function(data, model_key) {
  cur <- data %>%
    ggplot(aes(x = fpr, 
               y = tpr, 
               group = sim,
               colour = model)) +
    geom_line(alpha = 0.01) +
    facet_wrap(model ~ ., nrow = 2) +
    coord_equal(ratio = 1) +
    labs(x = 'False Positive Rate',
         y = 'True Positive Rate') +
    theme(legend.position = 'none')
  cur
}


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
    spread_draws(mst, diff1_mst, diff2_mst, diff3_mst, temp, lag1_temp)

  comp_var <- disc_best %>%
    spread_draws(b[i, f]) %>%
    spread(i, b) %>%
    filter(str_detect(f, 'fact_mybin')) %>%
    mutate(type = str_remove_all(f, '[0-9]'),
           age = as.numeric(str_extract(f, '[0-9]+'))) %>%
    ungroup

  # need to confirm correct line-up
  cv <- split(comp_var, comp_var$type) %>%
    map(., ~ .x %>% arrange(.chain, .iteration, age))

  # first element is just temporal effect
  core <- full_join(comp_const, cv[[1]], 
                    by = c('.chain', '.iteration', '.draw')) %>%
    mutate(eff_mst = mst.x + mst.y,
           eff_diff1_mst = diff1_mst.x + diff1_mst.y,
           eff_diff2_mst = diff2_mst.x + diff2_mst.y,
           eff_diff3_mst = diff3_mst.x + diff3_mst.y,
           eff_temp = temp.x + temp.y,
           eff_lag1_temp = lag1_temp.x + lag1_temp.y) %>%
    dplyr::select(.chain, .iteration, .draw, f, age, 
                  eff_mst, eff_diff1_mst, eff_diff2_mst, eff_diff3_mst,
                  eff_temp, eff_lag1_temp) %>%
    arrange(.chain, .iteration, .draw, age)

  # need to confirm correct line-up
  by_taxon <- map(cv[-1], ~ full_join(core, .x, by = c('.chain', '.iteration', 'age'))) %>%
    map(., ~ .x %>%
        mutate(taxon_eff_mst = eff_mst + mst,
               taxon_eff_diff1_mst = eff_diff1_mst + diff1_mst,
               taxon_eff_diff2_mst = eff_diff2_mst + diff2_mst,
               taxon_eff_diff3_mst = eff_diff3_mst + diff3_mst,
               taxon_eff_temp = eff_temp + temp,
               taxon_eff_lag1_temp = eff_lag1_temp + lag1_temp)) %>%
    reduce(bind_rows) %>%
    gather(key, value, 
           -.chain, -.iteration, 
           -.draw.y, -.draw.x,
           -f.x, -f.y, 
           -age, 
           -`(Intercept)`, -type,
           -eff_mst, -eff_diff1_mst, -eff_diff2_mst, -eff_diff3_mst, 
           -eff_temp, -eff_lag1_temp,
           -diff1_mst, -diff2_mst, -diff3_mst, -mst,
           -lag1_temp, -temp) %>%
    filter(!is.na(type)) %>%
    mutate(key = fct_recode(key, 
                            geo_range = 'taxon_eff_mst',
                            geo_change1 = 'taxon_eff_diff1_mst',
                            geo_change2 = 'taxon_eff_diff2_mst',
                            geo_change3 = 'taxon_eff_diff3_mst',
                            temp_now = 'taxon_eff_temp',
                            temp_lag = 'taxon_eff_lag1_temp'),
           key = fct_relevel(key, 
                             'geo_range', 
                             'geo_change1', 'geo_change2', 'geo_change3', 
                             'temp_now', 'temp_lag'),
           type = fct_recode(type,
                             diatoms = 'fossil_group:fact_mybin:D:',
                             forams = 'fossil_group:fact_mybin:F:',
                             nannofossil_misc = 'fossil_group:fact_mybin:N:',
                             radiolarians = 'fossil_group:fact_mybin:R:'))

  out <- by_taxon %>%
    ggplot(aes(x = age, y = value)) +
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
    stat_lineribbon(size = 0.5) +
    scale_fill_brewer() +
    facet_grid(type ~ key, scales = 'free_y') +
    theme(legend.position = 'bottom') +
    NULL
  
  out
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
  temp_est <- as_tibble(temp_est) %>% 
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
    dplyr::select(fullname, mst, relage) %>% 
    left_join(., temp_est, by = c('fullname', 'relage'))

  ext_plot <- full_est %>%
    ggplot(aes(x = relage, y = value)) + 
    stat_pointinterval(.width = c(0.5, 0.8)) +
    facet_grid(~ fullname) + 
    labs(x = 'Age (My)', y = 'P(T = t | T >= t, x)')

  range_plot <- full_est %>%
    group_by(fullname, relage) %>%
    summarise(mst = mean(mst)) %>%
    ggplot(aes(x = relage, y = mst)) + 
    geom_line() +
    facet_grid(~ fullname) + 
    labs(x = 'Age (My)', y = 'Geographic range (standardized)')

  out <- list(ext_plot, range_plot)
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
    spread_draws(b[i, f], `(Intercept)`) %>%
    filter(str_detect(f, pattern = 'fact_relage'),
           !str_detect(f, pattern = 'fossil_group')) %>%
    mutate(cr = `(Intercept)` + b) %>%    # keep log-odds scale
    mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
      arrange(age)

  # components from fossil_groups
  db <- model %>%
    spread_draws(b[i, f]) %>%
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
    stat_lineribbon(size = 0.5) +
    scale_fill_brewer() +
    facet_grid(type ~ .) +
    labs(x = 'Age (My)', y = 'P(T = t | T >= t, x)') +
    theme(legend.position = 'bottom') +
    NULL

  db
}


#' Faceted plot of AUC time series
#' 
#' This spews a bunch of error messages by they are handeled through exception handling, so the error messages are actually warnings :P.
#' 
#' @param data tibble of data used to fit models
#' @param model_pp list of model posterior estimates
#' @param model_key vector characters naming each model
#' @return ggplot object
plot_roc_series <- function(data, model_pp, model_key) {
  # important exception functions
  safe_roc <- safely(roc)
  safe_auc <- safely(auc)
  
  out <- list()
  for(kk in seq(length(model_pp))) {
    est <- list()

    # gather the predicted vs observed event
    for(ii in seq(max(data$mybin))) {
      tw <- model_pp[[kk]][, data$fact_mybin == ii]
      ew <- data %>%
        filter(fact_mybin == ii) %>%
        dplyr::select(event)
      ew <- as.data.frame(ew)[, 1]
      oo <- list()
   
      # calculate roc for each posterior draw
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
                                   to = model_key),
           model = factor(model, levels = rev(model_key)))

  rects <- get_geotime_box(range(roc_ts$key))


  # calculate breaks
  brks <- seq(min(roc_ts$key), max(roc_ts$key), by = 5) %>%
    round(., -1) %>%
    unique(.)
  
  roc_ts <- roc_ts %>%
    ggplot() +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = 'gray80', alpha = 0.8) +
    stat_lineribbon(aes(x = key, y = value), size = 0.5) +
    scale_fill_brewer(name = 'CI', aes(x = key, y = value)) + 
    scale_x_reverse(breaks = brks) +
    theme(legend.position = 'bottom') +
    NULL
    
  roc_ts
}


#' Visualize aspects of the neptune database
#'
#' This function is only called for its side effects!
#'
#' @param .data tibble of neptune data
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
view_neptune <- function(.data, name = 'full', path = '../results/figure/') {
  # occurrences through time, labeled if LAD
  octg <-  
    .data %>%
    mutate(state = case_when(event == 0 ~ 'Standard',
                             event == 1 ~ 'Last'),
           fossil_group = case_when(fossil_group == 'D' ~ 'Diatoms',
                                    fossil_group == 'R' ~ 'Radiolaria',
                                    fossil_group == 'F' ~ 'Foraminifera',
                                    fossil_group == 'N' ~ 'Calc. nanno.')) %>%
    ggplot(aes(x = mybin, fill = state)) +
    stat_bin() +
    facet_grid(fossil_group ~ ., switch = 'y') +
    scale_fill_manual(name = 'Occurrence type',
                      values = c('goldenrod', 'skyblue')) +
    theme(legend.position = 'bottom') +
    labs(title = 'Occurrences', x = 'Time (My before present)', y = 'Count')

  filename <- paste0(path, '/occ_time_label_', name, '.png')
  ggsave(filename = filename,
         plot = octg, 
         width = 4, height = 6)
  
  
  # relative "abundance" of microfossil_groups over time
  ocag <- .data %>%
    ggplot(aes(x = mybin, fill = fossil_group)) +
    geom_histogram(position = 'fill')
  
  filename <- paste0(path, '/abn_time_stack_', name, '.png')
  ggsave(filename = filename,
         plot = ocag, 
         width = 6, height = 6)
  
  # occurrences by relage
  ocrg <- .data %>% 
    group_by(fullname) %>%
    dplyr::summarize(maxage = max(relage),
                     fossil_group = plurality(fossil_group),
                     died = any(event == 1)) %>%
    mutate(state = case_when(died == 0 ~ 'Extant',
                             died == 1 ~ 'Extinct'),
           fossil_group = case_when(fossil_group == 'D' ~ 'Diatoms',
                                    fossil_group == 'R' ~ 'Radiolaria',
                                    fossil_group == 'F' ~ 'Foraminifera',
                                    fossil_group == 'N' ~ 'Calc. nanno.')) %>%
    ggplot(aes(x = maxage, fill = state)) +
    stat_bin() +
    facet_grid(fossil_group ~ ., switch = 'y') +
    scale_fill_manual(name = 'State', 
                      values = c('goldenrod', 'skyblue')) +
    theme(legend.position = 'bottom') +
    labs(title = 'Age distribution', x = 'Age (My)', y = 'Count')
  
  filename <- paste0(path, '/age_label_', name, '.png')
  ggsave(filename = filename,
         plot = ocrg, 
         width = 4, height = 6)
  
  
  
  # make a plot of a random selection of species
  set.seed(100)
  srg <- .data %>% 
    group_by(fullname) %>%
    sample_n_groups(size = 8) %>%
    ungroup %>%
    ggplot(aes(x = relage, y = maxgcd, group = fullname, colour = fullname)) +
    geom_line() +
    geom_point() +
    theme(legend.position = 'bottom')

  filename <- paste0(path, '/range_time_', name, '.png')
  ggsave(filename = filename,
         plot = srg, 
         width = 8, height = 6)
  
  # lots of little code here
  # for FAD/LAD accumulation curves
  ft <- .data %>%
    group_by(fullname) %>%
    summarize(fad = max(mybin),
              fossil_group = plurality(fossil_group)) %>%
    group_by(fossil_group, fad) %>%
    summarize(n = n()) %>%
    arrange(desc(fad)) %>%
    mutate(nsum = cumsum(n),
           time = fad)
  
  # LADs over time
  lt <- .data %>%
    group_by(fullname) %>%
    filter(!all(event == 0)) %>%
    summarize(lad = max(mybin),
              fossil_group = plurality(fossil_group)) %>%
    group_by(fossil_group, lad) %>%
    summarize(n = n()) %>%
    arrange(desc(lad)) %>%
    mutate(nsum = cumsum(n),
           time = lad)
  
  # put them together
  ccg <- bind_rows(ft, lt, .id = 'type') %>%
    mutate(type = plyr::mapvalues(type, 1:2, c('FAD', 'LAD'))) %>%
    ggplot(aes(x = time, y = nsum, colour = type, group = type)) +
    geom_line() +
    facet_grid(~ fossil_group) +
    labs(x = 'Time (My)', y = 'Cummulative count')

  filename <- paste0(path, '/fad_lad_count_wide_', name, '.png')
  ggsave(filename = filename,
         plot = ccg, 
         width = 6, height = 3)
  
  ccg2 <- ccg + facet_grid(fossil_group ~ ., switch = 'y', scales = 'free_y')
  filename <- paste0(path, '/fad_lad_count_tall_', name, '.png')
  ggsave(filename = filename,
         plot = ccg2, 
         width = 4, height = 6)
}


#' Visualize aspects of model fit
#'
#' This function is only called for its side effects!
#'
#' @param fit_list list of model fits
#' @param .data tibble of neptune data
#' @param key vector of model names
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
p_model <- function(fit_list, .data, key, name, path) {
  # posterior probability of observation surviving
  pp_prob <- future_map(fit_list, 
                        ~ posterior_linpred(.x, transform = TRUE, draws = 100))
 
  eroc <- map(pp_prob, ~ apply(.x, 1, function(y) roc(.data$event, y)))
  
  # ROC curve plot
  roc_df <- map(eroc, function(x) 
             imap(x, ~ tibble(sim = .y, 
                              fpr = 1 - .x$specificities,
                              tpr = .x$sensitivities))) %>%
    map(., ~ reduce(.x, rbind)) %>%
    imap(., ~ add_column(.x, mod = .y)) %>%
    reduce(., rbind) %>%
    as_tibble(.) %>%
    # long hand way that doesn't need plyr
    mutate(model = case_when(mod == 1 ~ key[1],
                             mod == 2 ~ key[2],
                             mod == 3 ~ key[3],
                             mod == 4 ~ key[4])) 
  cur <- plot_roc_curve(roc_df)
  fn <- paste0(path, '/roc_curve_', name, '.png')
  ggsave(filename = fn,
         plot = cur,
         width = 5, height = 8)
  
  # get AUC values for the above
  # plot as histogram
  auc_hist <- map(eroc, ~ map(.x, function(y) auc(y)[[1]])) %>%
    reshape2::melt(.) %>%
    as_tibble %>%
    rename(model = L1,
           draw = L2) %>%
    mutate(model_name = plyr::mapvalues(model, unique(model), key),
           model_name = factor(model_name, levels = key)) %>%
    ggplot(aes(x = value, y = model_name)) +
    geom_halfeyeh(.width = c(0.5, 0.8)) +
    labs(x = 'AUC', y = NULL) +
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15)) +
    NULL
  fn <- paste0(path, '/auc_hist_', name, '.png')
  ggsave(filename = fn,
         plot = auc_hist,
         width = 6, height = 6)
  fn <- paste0(path, '/auc_hist_zoom_', name, '.png')
  ggsave(filename = fn,
         plot = auc_hist + 
           xlim(0.5, 1) +
           geom_vline(xintercept = 0.5, colour = 'red'),
         width = 6, height = 6)
}


#' Visualize aspects of model fit
#'
#' This function is only called for its side effects!
#'
#' @param fit_list list of model fits
#' @param fit_list list of model fits
#' @param .data tibble of neptune data
#' @param key vector of model names
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
p_model_time <- function(fit_list, .data, key, name, path) {
  # posterior probability of observation surviving
  pp_prob <- future_map(fit_list, 
                        ~ posterior_linpred(.x, transform = TRUE, draws = 100))
  eroc <- map(pp_prob, ~ apply(.x, 1, function(y) roc(.data$event, y)))
  
  # roc as timeseries to see best and worst times
  roc_ts <- plot_roc_series(.data, pp_prob, key) +
    NULL
    #coord_cartesian(ylim = c(0.4, 1), xlim = c(0, 62)) +
  fn <- paste0(path, '/auc_ts_tiny_', name, '.png')

  ## kludge to generate pastable geotime scale
  #ggsave(filename = fn,
  #       plot = gggeo_scale(roc_ts, 
  #                          dat = 'epochs', 
  #                          size = 3, 
  #                          rot = 90,
  #                          height = 0.2),
  #       width = 8, height = 6)
  
  # facet
  roc_ts <- roc_ts +
    facet_grid(model ~ .) +
    labs(y = 'AUC', x = 'Time (Mya)') +
    geom_hline(yintercept = 0.5, colour = 'red', linetype = 'dashed') +
    NULL

  #srts <- gggeo_scale(roc_ts, dat = 'epochs', size = 3, rot = 90, height = 0.2)

  fn <- paste0(path, '/auc_ts_', name, '.png')
  ggsave(filename = fn,
         plot = roc_ts, #srts
         width = 11, height = 8.5)
}

#' Visualize aspects of model fit
#'
#' This function is only called for its side effects!
#'
#' @param fit_list list of model fits
#' @param fit_list list of model fits
#' @param .data tibble of neptune data
#' @param key vector of model names
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
p_model_taxon <- function(fit_list, .data, key, name, path) {
  # posterior probability of observation surviving
  pp_prob <- future_map(fit_list, 
                        ~ posterior_linpred(.x, transform = TRUE, draws = 100))
  # view through taxonomic window
  pp_taxon <- map(pp_prob, ~ split(data.frame(t(.x)), 
                                   .data$fossil_group)) %>%
    map(., function(x) map(x, ~ t(.x)))
  counti_taxon <- split(.data, .data$fossil_group)
  
  
  auc_taxon <- map(pp_taxon, function(x) 
                   map2(x, counti_taxon, ~ 
                        apply(.x, 1, function(a) roc(.y$event, a)))) %>%
    map(., function(a) 
        map(a, function(d)
            map_dbl(d, ~ auc(.x)[[1]]))) %>%
    set_names(key) %>%
    reshape2::melt(.) %>%                # clean up nested list
    rename(taxon = L2,
           model = L1) %>%
    mutate(model = factor(model, levels = rev(key)),
           taxon = case_when(taxon == 'D' ~ 'Diatoms',
                             taxon == 'R' ~ 'Radiolaria',
                             taxon == 'F' ~ 'Foraminifera',
                             taxon == 'N' ~ 'Calc. nanno.')) %>%
    ggplot(aes(x = value, y = model)) +
    geom_halfeyeh(.width = c(0.5, 0.8)) +
    facet_wrap(~ taxon) +
    labs(x = 'AUC', y = NULL) 
  fn <- paste0(path, '/auc_taxon_', name, '.png')
  ggsave(filename = fn,
         plot = auc_taxon,
         width = 8, height = 8)
}

#' Visualize aspects of model fit
#'
#' This function is only called for its side effects!
#'
#' @param fit_list list of model fits
#' @param .data tibble of neptune data
#' @param key vector of model names
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
p_model_taxon_time <- function(fit_list, .data, key, name, path) {
  # posterior probability of observation surviving

  auc_taxon_time <- ins_roc_taxon_time(fit_list, .data, key)
  
  rects <- get_geotime_box(range(auc_taxon_time$time))
  
  brks <- seq(min(auc_taxon_time$time), max(auc_taxon_time$time), by = 5) %>%
    round(., -1) %>%
    unique(.)
  
  auc_taxon_time <- auc_taxon_time %>%
    ggplot() +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = 'gray80', alpha = 0.8) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    stat_lineribbon(aes(x = time, y = value), size = 0.5) + 
    scale_fill_brewer() +
    scale_x_reverse(name = 'CI', breaks = brks) +
    facet_grid(fossil_group ~ model) +
    theme(legend.position = 'bottom') +
    labs(x = 'Time (Mya)', y = 'AUC') +
    NULL
  fn <- paste0(path, '/auc_taxon_time_', name, '.png')
  ggsave(filename = fn,
         plot = auc_taxon_time,
         width = 11, height = 8.5)
}





#' convenience function to get posterior predictive distribution
#' 
#' convenience function to get posterior predictive distribution
#'
#' @param fit list of model fits
#' @param .data tibble of neptune data
#' @return list
get_pred <- function(fit, .data) {
  future_map2(fit, .data,
              ~ posterior_linpred(object = .x, 
                                  newdata = .y,
                                  draws = 100))
}

#' Visualize aspects of cross-validation fit
#'
#' This function is only called for its side effects!
#'
#' @param fit list of model fits
#' @param .data list of tibbles
#' @param key vector of model names
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
cv_model <- function(fit, .data, key, name, path) {
  # predict the test data
  pred <- get_pred(fit, .data)
  
  # calculate ROC/AUC for the predictions of test data
  # to determine how good our out of sample predictions are
  pred_auc <- map2(pred, .data, post_roc) %>%
    map(., function(x) map_dbl(x, ~ auc(.x))) %>%
    set_names(., names(fit))
  
  # make a plot about out-of-sample predictive accuracy
  oos_auc <- bind_cols(pred_auc) %>% 
    gather() %>%
    separate(key, into = c('mod', 'fold'), sep = '\\_') %>%
    mutate(mod = plyr::mapvalues(mod, unique(mod), key),
           mod = factor(mod, levels = key)) %>%
    ggplot(aes(x = value, y = mod)) +
    geom_halfeyeh(.width = c(0.5, 0.8)) +
    labs(x = 'AUC', y = NULL) +
    scale_colour_brewer() +
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15)) +
    NULL
  fn <- paste0(path, '/fold_auc_', name, '.png')
  ggsave(filename = fn,
         plot = oos_auc,
         width = 6, height = 6)
  fn <- paste0(path, '/fold_auc_zoom_', name, '.png')
  ggsave(filename = fn,
         plot = oos_auc + 
           xlim(0.5, 1) +
           geom_vline(xintercept = 0.5, colour = 'red'),
         width = 6, height = 6)
}





#' Visualize aspects of cross-validation fit -- by time
#'
#' This function is only called for its side effects!
#'
#' @param fit list of model fits
#' @param .data list of tibbles
#' @param key vector of model names
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
cv_model_time <- function(fit, .data, key, name, path) {
  # predict the test data
  pred <- get_pred(fit, .data)
  
  time_auc <- map2(pred, .data, get_auc_time) %>%
    set_names(., names(fit))
  
  # could this be done with enframe?
  ta <- reshape2::melt(time_auc) %>% 
    as_tibble %>%
    mutate(time = parse_double(L2)) %>%
    separate(L1, into = c('mod', 'fold'), sep = '\\_') %>%
    mutate(mod = plyr::mapvalues(mod, unique(mod), key),
           mod = factor(mod, levels = rev(key)))
  
  rects <- get_geotime_box(range(ta$time))
  
  brks <- seq(min(ta$time), max(ta$time), by = 5) %>%
    round(., -1) %>%
    unique(.)
  
  ta <- ta %>%
    ggplot() +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = 'gray80', alpha = 0.8) +
    stat_lineribbon(aes(x = time, y = value), size = 0.5) +
    scale_fill_brewer() +
    scale_x_reverse(breaks = brks) +
    theme(legend.position = 'bottom') +
    labs(y = 'AUC', x = 'Time (Mya)') +
    NULL
    #coord_cartesian(ylim = c(0.4, 1), xlim = c(0, 50)) +
  fn <- paste0(path, '/fold_auc_time_tiny_', name, '.png')

  ## kludge to generate pastable geotime scale
  #ggsave(filename = fn,
  #       plot = gggeo_scale(ta, 
  #                          dat = 'epochs', 
  #                          size = 3, 
  #                          rot = 90,
  #                          height = 0.2),
  #       width = 11, height = 8.5)

  ta <- ta +
    facet_grid(mod ~ .) +
    geom_hline(yintercept = 0.5, colour = 'red', linetype = 'dashed') +
    NULL
  
  #sta <- gggeo_scale(roc_ts, dat = 'epochs', size = 3, rot = 90, height = 0.2)

  fn <- paste0(path, '/fold_auc_time_', name, '.png')
  ggsave(filename = fn,
         plot = ta, #sta,
         width = 11, height = 8.5)
}



#' Visualize aspects of cross-validation fit -- by taxon
#'
#' This function is only called for its side effects!
#'
#' @param fit list of model fits
#' @param .data list of tibbles
#' @param key vector of model names
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
cv_model_taxon <- function(fit, .data, key, name, path) {
  # predict the test data
  pred <- get_pred(fit, .data)
 
  pp_taxon <- map2(pred, .data, 
                   ~ split(data.frame(t(.x)), 
                           .y$fossil_group)) %>%
    map(., function(x) map(x, ~ t(.x)))
  
  counti_fold_taxon <- map(.data, ~ split(.x, .x$fossil_group))
  
  fold_auc_taxon <- map2(pp_taxon, counti_fold_taxon, 
       ~ map2(.x, .y, 
              ~ apply(.x, 1, function(a) fast_auc(a, .y$event)))) %>%
    reshape2::melt(.) %>%
    as_tibble(.) %>%
    rename(taxon = L2,
           model = L1) %>%
    separate(., col = model, into = c('model', 'fold'), sep = '_') %>%
    mutate(taxon = case_when(taxon == 'D' ~ 'Diatoms',
                             taxon == 'R' ~ 'Radiolaria',
                             taxon == 'F' ~ 'Foraminifera',
                             taxon == 'N' ~ 'Calc. nanno.'),
           model = case_when(model == 'mod1' ~ key[1],
                             model == 'mod2' ~ key[2],
                             model == 'mod3' ~ key[3],
                             model == 'mod4' ~ key[4]),
           model = factor(model, levels = rev(key))) %>%
    ggplot(aes(x = value, y = model)) +
    geom_halfeyeh(.width = c(0.5, 0.8)) +
    facet_wrap(~ taxon) +
    labs(x = 'AUC', y = NULL)
  fn <- paste0(path, '/fold_auc_taxon_', name, '.png')
  ggsave(filename = fn,
         plot = fold_auc_taxon,
         width = 8, height = 8)
}


#' Visualize aspects of cross-validation fit -- by taxon/time
#'
#' This function is only called for its side effects!
#'
#' @param fit list of model fits
#' @param .data list of tibbles
#' @param key vector of model names
#' @param name character length 1 vector describing data
#' @param path character length 1 vector describing directory to drop figures into -- needs trailing slash!
#' @return NULL
cv_model_taxon_time <- function(fit, .data, key, name, path) {

  fold_auc_taxon_time <- oos_roc_taxon_time(fit, .data, key)
  
  rects <- get_geotime_box(range(fold_auc_taxon_time$time))
 
  brks <- seq(min(fold_auc_taxon_time$time), 
              max(fold_auc_taxon_time$time), by = 5) %>%
    round(., -1) %>%
    unique(.)

  fold_auc_taxon_time <- fold_auc_taxon_time %>%
    ggplot() +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = 'gray80', alpha = 0.8) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    stat_lineribbon(aes(x = time, y = value), size = 0.5) +
    scale_fill_brewer() +
    scale_x_reverse(breaks = brks) +
    facet_grid(fossil_group ~ model) +
    theme(legend.position = 'bottom') +
    labs(x = 'Time (Mya)', y = 'AUC') +
    NULL
  fn <- paste0(path, '/fold_auc_taxon_time_', name, '.png')
  ggsave(filename = fn,
         plot = fold_auc_taxon_time,
         width = 11, height = 8.5)
}
