library(pacman)

p_load(tidyverse, janitor, tidybayes, furrr, arm, rstanarm, 
       bayesplot, pROC, here)

source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper02_plot_foo.r'))
source(here('R', 'helper03_misc_foo.r'))
source(here('R', 'helper04_stan_utility.r'))
source(here('R', 'helper05_roc_utility.r'))
source(here('R', 'helper07_diffsafe.r'))

# important constants
plan(multiprocess)

# useful graphical constants
theme_set(theme_bw(base_size = 20))

# get data in
counti <- read_rds(here('data', 'counting.rds'))

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)

# read in model fits
disc_fit <- read_rds(here('results', 'disc_fit.rds'))
disc_best <- disc_fit[[1]]

## check how much information learned
#priors <- prior_summary(disc_best)
#prmn <- priors$prior$location
#prsc <- priors$prior$scale
#
#posc <- disc_best %>% 
#  spread_draws(maxgcd,
#               diff_maxgcd,
#               temp,
#               lag1_temp) %>% 
#clean_names %>% 
#dplyr::summarize_at(vars(-chain, -iteration), funs(sd(.)))
#
## 1 = lot of shrinkage, 0 = no shrinkage
#posterior_shrinkage <- 1 - (posc / prsc) ^ 2

# posterior intervals
effect_group <- 
  disc_best %>%
  gather_draws(`(Intercept)`, 
               mst,
               diff1_mst,
               diff2_mst,
               diff3_mst,
               temp,
               lag1_temp) %>%
  ungroup() %>%
  mutate(.ids = case_when(.variable == '(Intercept)' ~ 
                            'Average log-odds \nof extinction',
                          .variable == 'mst' ~
                            'Geographic range',
                          .variable == 'diff1_mst' ~
                            'Change in geographic \nrange (1 lag)',
                          .variable == 'diff2_mst' ~
                            'Change in geographic \nrange (2 lag)',
                          .variable == 'diff3_mst' ~
                            'Change in geographic \nrange (3 lag)',
                          .variable == 'temp' ~
                            'Global temperature',
                          .variable == 'lag1_temp' ~
                            'Previous global \ntemperature'),
         .ids = fct_inorder(.ids),
         .ids = fct_rev(.ids)) %>%
  ggplot(aes(y = .ids, x = .value)) +
  geom_halfeyeh(.width = c(0.8, 0.5)) +
  labs(x = 'Posterior estimate for parameter',
       y = 'Group-level average covariate effects')
ggsave(filename = here('results', 'figure', 'effect_est.png'),
       plot = effect_group, width = 10, height = 8)

# variance components
# between taxonmic groups
vary_between <- 
  disc_best %>%
  spread_draws(Sigma[var1, var2]) %>%
  filter(str_count(var1, ':') == 1) %>%
  separate(., var1, into = c('group_time', 'var1'), sep = ':') %>%
  filter(var1 == var2) %>%
  ungroup() %>%
  mutate(.variable = case_when(var1 == '(Intercept)' ~ 
                                 'Average log-odds of \nextinction',
                               var1 == 'mst' ~
                                 'Geographic range',
                               var1 == 'diff1_mst' ~
                                 'Change in geographic \nrange (1 lag)',
                               var1 == 'diff2_mst' ~
                                 'Change in geographic \nrange (2 lag)',
                               var1 == 'diff3_mst' ~
                                 'Change in geographic \nrange (3 lag)',
                               var1 == 'temp' ~
                                 'Global temperature',
                               var1 == 'lag1_temp' ~
                                 'Previous global \ntemperature'),
         .variable = fct_relevel(.variable,
                                 c('Geographic range'),
                                 after = 1),
         .variable = fct_rev(.variable))
  
vb_gg <- 
  vary_between %>%
  ggplot(aes(x = Sigma, y = .variable)) +
  geom_halfeyeh(.width = c(0.8, 0.5)) +
  labs(x = 'Variance of effect between taxonomic groups',
       y = 'Covariate effect') 
ggsave(filename = here::here('results', 'figure', 'variance_between.png'),
       plot = vb_gg)

# within taxonomic groups
vary_within <- 
  disc_best %>%
  spread_draws(Sigma[var1, var2]) %>%
  filter(str_count(var1, ':') == 2) %>%
  separate(., var1, into = c('group_taxon', 'group_time', 'var1'), sep = ':') %>%
  filter(var1 == var2) %>%
  ungroup() %>%
  mutate(.variable = case_when(var1 == '(Intercept)' ~ 
                                 'Average log-odds of \nextinction',
                               var1 == 'mst' ~
                                 'Geographic range',
                               var1 == 'diff1_mst' ~
                                 'Change in geographic \nrange (1 lag)',
                               var1 == 'diff2_mst' ~
                                 'Change in geographic \nrange (2 lag)',
                               var1 == 'diff3_mst' ~
                                 'Change in geographic \nrange (3 lag)',
                               var1 == 'temp' ~
                                 'Global temperature',
                               var1 == 'lag1_temp' ~
                                 'Previous global \ntemperature'),
         .variable = fct_relevel(.variable,
                                 c('Geographic range'),
                                 after = 1),
         .variable = fct_rev(.variable))

vw_gg <- 
  vary_within %>%
  ggplot(aes(x = Sigma, y = .variable)) +
  geom_halfeyeh(.width = c(0.8, 0.5)) +
  labs(x = 'Variance of effect over time',
       y = 'Covariate effect')
ggsave(filename = here::here('results', 'figure', 'variance_within.png'),
       plot = vw_gg)



# base-line hazard plot
hazard_plot <- 
  disc_best %>%
  spread_draws(b[i, f], `(Intercept)`) %>%
  filter(str_detect(f, pattern = 'fact_relage'),
         !str_detect(f, pattern = 'fossil_group')) %>%
  mutate(cr = invlogit(`(Intercept)` + b)) %>%
  mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
  ggplot(aes(y = cr, x = age)) + 
  stat_lineribbon() +
  scale_fill_brewer() +
  labs(x = 'Age (My)', y = 'P(T = t | T >= t, x)')
ggsave(filename = here('results', 'figure', 'hazard_baseline.png'),
       plot = hazard_plot, width = 6, height = 4)

# can also then compare between taxonomic groups because i have that calculated
# faceted baseline hazard plot by fossil_group
db <- plot_taxon_hazard(disc_best)
ggsave(filename = here('results', 'figure', 'hazard_bygroup.png'), 
       plot = db,
       width = 6, height = 8)


# make plots of effect change over time for
by_taxon <- plot_taxon_covariate_time(disc_best) + 
  scale_x_reverse()
ggsave(filename = here('results', 'figure', 'eff_time_group.png'), 
       plot = by_taxon,
       height = 6, width = 8)


## look at hazard over duration for a few examples
## grab a random selection of species
#set.seed(100)
#risk_plots <- plot_risk_time(counti_trans, disc_best, nsp = 4)
#ggsave(filename = here('results', 'figure', 'relrisk_ext.png'), 
#       plot = risk_plots[[1]],
#       height = 6, width = 8)
#ggsave(filename = here('results', 'figure', 'relrisk_range.png'), 
#       plot = risk_plots[[2]],
#       height = 6, width = 8)
