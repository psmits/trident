# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(brms)
library(bayesplot)

# misc
library(pROC)
library(ROCR)
source('../R/helper01_process_foo.r')
source('../R/helper03_misc_foo.r')
source('../R/helper04_stan_utility.r')

# important constants
options(mc.cores = parallel::detectCores())

# get data in
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

# all of the forms
forms <- list(form, form2, form3, form4)

# partial evaluation to fill in big details once
part_glmer <- partial(stan_glmer, 
                      family = 'binomial',
                      prior_intercept = normal(-2, 5, autoscale = FALSE),
                      prior_aux = cauchy(0, 1, autoscale = FALSE),
                      thin = 4, 
                      adapt_delta = 0.999999)

# some models have differen't priors
priors <- list(prior = normal(c(0, 0, -1, 0), 
                              rep(1, 4), 
                              autoscale = FALSE),
               prior = normal(c(0, 0, -1, 0), 
                              rep(1, 4), 
                              autoscale = FALSE),
               prior = normal(c(0, -1), 
                              rep(1, 2), 
                              autoscale = FALSE),
               prior = normal(c(0, -1), 
                              rep(1, 2), 
                              autoscale = FALSE))

# combine formula, prior, and model
by_formula <- map2(forms, priors, ~ partial(part_glmer, 
                                            formula = .x, 
                                            prior = .y))

# by_formula list of (partial) functions
# counti_accum list of folds
to_fit <- cross2(by_formula, counti_accum)

# all that is missing is the data
fit <- map(to_fit, ~ .x[[1]](data = .x[[2]]))
# everything is multiples of 4
#   1:4 first form applied to each fold
#   5:8 second form applied to each fold
#   9:12 third form applied to each fold
#   13:16 third form applied to each fold

# help us remember this fact
#   folds numbered from oldest (1) to youngest (4)
nn <- c('mod1_', 'mod2_', 'mod3_', 'mod4_') %>%
  map(., function(a) map(1:4, ~ paste0(a, 'fold', .x)))
fit <- set_names(fit, flatten(nn))


write_rds(fit, path = '../data/training_fit.rds')
