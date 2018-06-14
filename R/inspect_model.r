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
                 `Sigma[fact_mybin:(Intercept),(Intercept)]`,
                 `Sigma[fact_relage:(Intercept),(Intercept)]`,
                 `Sigma[fossil.group:(Intercept),(Intercept)]`) %>%
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
  gather_samples(`Sigma[fact_mybin:(Intercept),(Intercept)]`,
                 `Sigma[fact_relage:(Intercept),(Intercept)]`,
                 `Sigma[fossil.group:(Intercept),(Intercept)]`) %>%
  ggplot(aes(y = term, x = estimate)) +
  geom_halfeyeh(.prob = c(0.9, 0.5))
ggsave(filename = '../doc/figure/variance_components.png',
       plot = vary_eye, width = 2, height = 4)


# base-line hazard plot
hazard_plot <- 
  disc_best %>%
  spread_samples(b[i, f], `(Intercept)`) %>%
  filter(str_detect(f, pattern = 'fact_mybin')) %>%
  mutate(cr = invlogit(`(Intercept)` + b)) %>%
  mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
  ggplot(aes(y = cr, x = age)) + 
  stat_lineribbon() +
  scale_fill_brewer() +
  labs(x = 'age (My)', y = 'P(T = t | T >= t, x)')
ggsave(filename = '../doc/figure/hazard_baseline.png',
       plot = hazard_plot, width = 6, height = 4)
