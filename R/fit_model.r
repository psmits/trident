# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/stan_utility.R')

# misc
source('../R/process_foo.r')

# important constants
options(mc.cores = parallel::detectCores())

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# transform the data for analysis
counti_trans <- prepare_analysis(counti)

#pp fit the model

# past and vary
form <- formula(event ~ temp + lag1_temp + maxgcd + diff_maxgcd 
                + (1 | fact_relage/fossil_group)
                + (1 + temp + lag1_temp + maxgcd + diff_maxgcd | 
                   fact_mybin/fossil_group))

# past but no vary
form2 <- update(form, ~ . - (1 + temp + lag1_temp + maxgcd + diff_maxgcd | 
                             fact_mybin/fossil_group) 
                      + (1 | fact_mybin/fossil_group))

# no past but vary
form3 <- update(form, ~ . - diff_maxgcd - lag1_temp 
                      - (1 + temp + lag1_temp + maxgcd + diff_maxgcd | 
                         fact_mybin/fossil_group) 
                      + (1 + temp + maxgcd | fact_mybin/fossil_group))

# no past or vary
form4 <- update(form3, ~ . - (1 + temp + maxgcd | fact_mybin/fossil_group) 
                       + (1 | fact_mybin/fossil_group))

forms <- list(form, form2, form3, form4)

# partial evaluation to fill in priors and ancillary details
part_glmer <- partial(stan_glmer, 
                      family = 'binomial',
                      data = counti_trans,
                      prior_intercept = normal(0, 10, autoscale = FALSE),
                      prior_aux = cauchy(0, 5, autoscale = FALSE),
                      thin = 4, 
                      adapt_delta = 0.99999999)

# some models have differen't priors
list_part_glmer <- 
  list(partial(part_glmer, prior = normal(c(0, 0, -1, 0), 
                                          rep(3, 4), 
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(c(0, 0, -1, 0), 
                                          rep(3, 4), 
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(c(0, -1), 
                                          rep(3, 2), 
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(c(0, -1), 
                                          rep(3, 2), 
                                          autoscale = FALSE)))

# map the data to the models
disc_fit <- map2(forms, list_part_glmer, ~ .y(.x))

write_rds(disc_fit, path = '../data/disc_fit.rds')
