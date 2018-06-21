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
    spread_samples(maxgcd, diff_maxgcd, `maxgcd:diff_maxgcd`)

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
           eff_diff_maxgcd = diff_maxgcd.x + diff_maxgcd.y,
           eff_interaction = `maxgcd:diff_maxgcd.x` + `maxgcd:diff_maxgcd.y`) %>%
  dplyr::select(.chain, .iteration, f, age, 
                eff_maxgcd, eff_diff_maxgcd, eff_interaction) %>%
  arrange(.chain, .iteration, age)

  # need to confirm correct line-up
  by_taxon <- map(cv[-1], ~ full_join(core, .x, by = c('.chain', '.iteration', 'age'))) %>%
    map(., ~ .x %>%
        mutate(taxon_eff_maxgcd = eff_maxgcd + maxgcd,
               taxon_eff_diff_maxgcd = eff_diff_maxgcd + diff_maxgcd,
               taxon_eff_interaction = eff_interaction + `maxgcd:diff_maxgcd`))
  by_taxon <- reduce(by_taxon, bind_rows)
  
  by_taxon <- by_taxon %>%
    gather(key, value, -.chain, -.iteration, -f.x, 
           -age, -f.y, -`(Intercept)`, -type,
           -eff_maxgcd, -eff_diff_maxgcd, -eff_interaction,
           -diff_maxgcd, -maxgcd, -`maxgcd:diff_maxgcd`) %>%
  filter(!is.na(type)) %>%
  ggplot(aes(x = age, y = value)) +
  stat_lineribbon() +
  scale_fill_brewer() +
  facet_grid(key ~ type)
  
  by_taxon
}
