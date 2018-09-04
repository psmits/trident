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
library(pROC)
library(ROCR)
source('../R/process_foo.r')

# important constants
options(mc.cores = parallel::detectCores())

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# transform the data for analysis
counti_trans <- prepare_analysis(counti)

# split the data into temporal folds
# there are only 63 temporal units
# assuming 5 folds, that's 12-13 time units per fold
# train [1], test [2]
# train [1,2], test [3]
# train [1,2,3], test [4]
# train [1,2,3,4], test [5]

counti_trans <- counti_trans %>%
  mutate(fold = break_my(mybin, number = 5))

# 1 is young, 5 is old; reverse!
counti_fold <- rev(split(counti_trans, counti_trans$fold))  
# combine some folds to make training sets
counti_accum <- accumulate(counti_fold, bind_rows) 
counti_accum <- counti_accum[-1]

form <- formula(event ~ temp + lag1_temp + maxgcd + diff_maxgcd 
                + (1 | fact_relage/fossil_group)
                + (1 + temp + lag1_temp + maxgcd + diff_maxgcd | 
                   fact_mybin/fossil_group))

glmer_part <- 
  partial(stan_glmer, 
          prior = normal(c(0, 0, -1, 0), rep(3, 4), autoscale = FALSE),
          prior_intercept = normal(0, 10, autoscale = FALSE),
          prior_aux = cauchy(0, 5, autoscale = FALSE))
fit <- map(counti_accum, ~ glmer_part(form, 
                                      family = 'binomial', 
                                      data = .x,
                                      adapt_delta = 0.999999, 
                                      thin = 4))

write_rds(fit, path = '../data/training_fit.rds')
