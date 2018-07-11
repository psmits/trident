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
# assuming 5 folds, that's 12 time units per fold
# train [1], test [2]
# train [1,2], test [3]
# train [1,2,3], test [4]
# train [1,2,3,4], test [5]

counti_trans <- counti_trans %>%
  mutate(fold = break_my(mybin, number = 5))

# 1 is young, 5 is old
counti_fold <- rev(split(counti_trans, counti_trans$fold))  # each fold
counti_accum <- accumulate(counti_fold, bind_rows) # combine some folds
counti_accum <- counti_accum[-1]

form1 <- formula(event ~ temp + lag1_temp + maxgcd + diff_maxgcd + 
                 (1 + maxgcd + diff_maxgcd | fact_mybin/fossil_group) + 
                 (1 | fact_relage/fossil_group))

fit <- map(counti_accum, ~ stan_glmer(form1, 
                                      family = 'binomial', 
                                      data = .x,
                                      adapt_delta = 0.999, 
                                      thin = 4))

write_rds(fit, path = '../data/training_fit.rds')
