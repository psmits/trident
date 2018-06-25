# data manipulation
library(tidyverse)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/stan_utility.R')

# misc
library(pROC)
source('../R/process_foo.r')
source('../R/plot_foo.r')

# important constants
options(mc.cores = parallel::detectCores())

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)

# read in model fits
disc_fit <- read_rds('../data/disc_fit.rds')
disc_best <- disc_fit[[1]]

# posterior intervals
interval_est <- disc_best %>%
  spread_samples(`(Intercept)`, 
                 maxgcd, 
                 diff_maxgcd,
                 `maxgcd:diff_maxgcd`,
                 temp,
                 lag1_temp,
                 `Sigma[fossil.group:fact_mybin:(Intercept),(Intercept)]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd,(Intercept)]`,
                 `Sigma[fossil.group:fact_mybin:diff_maxgcd,(Intercept)]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd:diff_maxgcd,(Intercept)]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd,maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:diff_maxgcd,maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd:diff_maxgcd,maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:diff_maxgcd,diff_maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd:diff_maxgcd,diff_maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd:diff_maxgcd,maxgcd:diff_maxgcd]`,
                 `Sigma[fossil.group:fact_relage:(Intercept),(Intercept)]`,
                 `Sigma[fact_mybin:(Intercept),(Intercept)]`,
                 `Sigma[fact_mybin:maxgcd,(Intercept)]`,
                 `Sigma[fact_mybin:diff_maxgcd,(Intercept)]`,
                 `Sigma[fact_mybin:maxgcd:diff_maxgcd,(Intercept)]`,
                 `Sigma[fact_mybin:maxgcd,maxgcd]`,
                 `Sigma[fact_mybin:diff_maxgcd,maxgcd]`,
                 `Sigma[fact_mybin:maxgcd:diff_maxgcd,maxgcd]`,
                 `Sigma[fact_mybin:diff_maxgcd,diff_maxgcd]`,
                 `Sigma[fact_mybin:maxgcd:diff_maxgcd,diff_maxgcd]`,
                 `Sigma[fact_mybin:maxgcd:diff_maxgcd,maxgcd:diff_maxgcd]`,
                 `Sigma[fact_relage:(Intercept),(Intercept)]`) %>%
  median_qi(.prob = c(0.9, 0.5))

effect_eye <- disc_best %>%
  gather_samples(`(Intercept)`, 
                 maxgcd, 
                 diff_maxgcd,
                 `maxgcd:diff_maxgcd`,
                 temp,
                 lag1_temp) %>%
  ggplot(aes(y = term, x = estimate)) +
  geom_halfeyeh(.prob = c(0.9, 0.5))
ggsave(filename = '../doc/figure/effect_est.png',
       plot = effect_eye, width = 4, height = 6)

vary_eye <- disc_best %>%
  gather_samples(`Sigma[fossil.group:fact_mybin:(Intercept),(Intercept)]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd,(Intercept)]`,
                 `Sigma[fossil.group:fact_mybin:diff_maxgcd,(Intercept)]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd:diff_maxgcd,(Intercept)]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd,maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:diff_maxgcd,maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd:diff_maxgcd,maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:diff_maxgcd,diff_maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd:diff_maxgcd,diff_maxgcd]`,
                 `Sigma[fossil.group:fact_mybin:maxgcd:diff_maxgcd,maxgcd:diff_maxgcd]`,
                 `Sigma[fossil.group:fact_relage:(Intercept),(Intercept)]`,
                 `Sigma[fact_mybin:(Intercept),(Intercept)]`,
                 `Sigma[fact_mybin:maxgcd,(Intercept)]`,
                 `Sigma[fact_mybin:diff_maxgcd,(Intercept)]`,
                 `Sigma[fact_mybin:maxgcd:diff_maxgcd,(Intercept)]`,
                 `Sigma[fact_mybin:maxgcd,maxgcd]`,
                 `Sigma[fact_mybin:diff_maxgcd,maxgcd]`,
                 `Sigma[fact_mybin:maxgcd:diff_maxgcd,maxgcd]`,
                 `Sigma[fact_mybin:diff_maxgcd,diff_maxgcd]`,
                 `Sigma[fact_mybin:maxgcd:diff_maxgcd,diff_maxgcd]`,
                 `Sigma[fact_mybin:maxgcd:diff_maxgcd,maxgcd:diff_maxgcd]`,
                 `Sigma[fact_relage:(Intercept),(Intercept)]`) %>%
  ggplot(aes(y = term, x = estimate)) +
  geom_halfeyeh(.prob = c(0.9, 0.5))
ggsave(filename = '../doc/figure/variance_components.png',
       plot = vary_eye, width = 2, height = 4)

# base-line hazard plot
hazard_plot <- 
  disc_best %>%
  spread_samples(b[i, f], `(Intercept)`) %>%
  filter(str_detect(f, pattern = 'fact_mybin'),
         !str_detect(f, pattern = 'fossil.group')) %>%
  mutate(cr = invlogit(`(Intercept)` + b)) %>%
  mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
  ggplot(aes(y = cr, x = age)) + 
  stat_lineribbon() +
  scale_fill_brewer() +
  labs(x = 'age (My)', y = 'P(T = t | T >= t, x)')
ggsave(filename = '../doc/figure/hazard_baseline.png',
       plot = hazard_plot, width = 6, height = 4)
# can also then compare between taxonomic groups because i have that calculated



# make plots of effect change over time for
#   3 covariates
#   4 taxonomic groups
by_taxon <- plot_taxon_covariate_time(disc_best)
ggsave(filename = '../doc/figure/eff_time_group.png', plot = by_taxon,
       height = 8, width = 8)



# look at hazard over duration for a few examples
# grab a random selection of species
set.seed(100)
temp <- counti_trans %>%
  group_by(fullname) %>%
  sample_n_groups(size = 4)

# estimate log-odds extinction for new species at each time
temp_est <- posterior_linpred(object = disc_best, newdata = temp)
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
  stat_pointinterval() +
  facet_grid(~ fullname) + 
  labs(x = 'Age (My)', y = 'log-odds extinction')
ggplot(filename = '../doc/figure/relrisk_ext.png', plot = ext_plot,
       height = 8, width = 8)

range_plot <- full_est %>%
  group_by(fullname, relage) %>%
  summarise(maxgcd = mean(maxgcd)) %>%
  ggplot(aes(x = relage, y = maxgcd)) + 
  geom_line() +
  facet_grid(~ fullname) + 
  labs(x = 'Age (My)', y = 'log-odds extinction')
ggplot(filename = '../doc/figure/relrisk_range.png', plot = range_plot,
       height = 8, width = 8)

full_plot <- full_est %>%
  gather(key, value, -fullname, -iterations, -relage, -row) %>%
  ggplot(aes(x = relage, y = value)) +
  geom_point() + 
  facet_grid(key ~ fullname, scales = 'free_y') +
  labs(x = 'Age (My)', y = 'log-odds extinction')
ggsave(filename = '../doc/figure/relrisk_full.png', plot = full_plot,
       height = 8, width = 8)
