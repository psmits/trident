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
library(here)
source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper03_misc_foo.r'))
source(here('R', 'helper04_stan_utility.r'))

# important constants
options(mc.cores = parallel::detectCores())

# the function to do this
fit_models <- function(.data) {
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
  # i have to admit that this is really cool -- beginning me would never have figured this out!
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

  fit
}


# given above function, fit cross-validation models to the three data sets

# get data in
# split the data into temporal folds
# there are only 63 temporal units
# assuming 5 folds, that's 12-13 time units per fold
# train [1], test [2]
# train [1,2], test [3]
# train [1,2,3], test [4]
# train [1,2,3,4], test [5]
# 1 is young, 5 is old; reverse!
counti_fold <- read_rds(here('data', 'counting.rds')) %>%
  prepare_analysis(.) %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)

counti_fold_rt <- read_rds(here('data', 'counting_restrict_time.rds')) %>%
  prepare_analysis(.) %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)

counti_fold_rl <- read_rds(here('data', 'counting_restrict_local.rds')) %>%
  prepare_analysis(.) %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)

# combine some folds to make training sets
counti_accum <- accumulate(counti_fold, bind_rows) [-1]
counti_accum_rt <- accumulate(counti_fold_rt, bind_rows) [-1]
counti_accum_rl <- accumulate(counti_fold_rl, bind_rows) [-1]

# fit models to the folds
fit_full <- fit_models(.data = counti_accum)
fit_rt <- fit_models(.data = counti_accum_rt)
fit_rl <- fit_models(.data = counti_accum_rl)

# write out cross-validation fits
write_rds(fit_full, path = here('results', 'training_fit.rds'))
write_rds(fit_rt, path = here('results', 'training_fit_rt.rds'))
write_rds(fit_rl, path = here('results', 'training_fit_rl.rds'))
