# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)

# misc
library(here)
source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper03_misc_foo.r'))
source(here('R', 'helper04_stan_utility.r'))

# important constants
options(mc.cores = parallel::detectCores())

#priors <- c(set_prior('normal(-2, 1)', class= 'Intercept'),
#            set_prior('normal(0, 1)', class = 'b'),
#            set_prior('normal(-1, 1)', class = 'b', coef = 'maxgcd'),
#            set_prior('cauchy(0, 1)', class = 'sd'),
#            set_prior('lkj(1)', class = 'cor'))
#
#form <- bf(event ~ temp + lag1_temp + maxgcd + diff_maxgcd 
#           + (1 | fact_relage/fossil_group)
#           + (1 + temp + lag1_temp + maxgcd + diff_maxgcd | 
#              fact_mybin/fossil_group))
#
#brmfit <- brm(formula = form,
#              data = counti_trans, 
#              family = bernoulli(), 
#              chains = 4, 
#              cores = 4, 
#              prior = priors)

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
                      prior_intercept = normal(-2, 5, autoscale = FALSE),
                      prior_aux = cauchy(0, 1, autoscale = FALSE),
                      prior_covariance = decov(regularization = 2),
                      thin = 4, 
                      adapt_delta = 0.999999)

# some models have differen't priors
list_part_glmer <- 
  list(partial(part_glmer, prior = normal(c(0, 0, -1, 0), 
                                          rep(1, 4), 
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(c(0, 0, -1, 0), 
                                          rep(1, 4), 
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(c(0, -1), 
                                          rep(1, 2), 
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(c(0, -1), 
                                          rep(1, 2), 
                                          autoscale = FALSE)))

# get data in
counti <- read_rds(here('data', 'counting.rds'))
counti_restrict_time <- read_rds(here('data', 'counting_restrict_time.rds'))
counti_restrict_local <- read_rds(here('data', 'counting_restrict_local.rds'))

# transform the data for analysis
counti_trans <- prepare_analysis(counti)
counti_rt_trans <- prepare_analysis(counti_restrict_time)
counti_rl_trans <- prepare_analysis(counti_restrict_local)

# map the data to the models
disc_fit <- map2(forms, list_part_glmer, ~ .y(.x, data = counti_trans))
disc_fit_rt <- map2(forms, list_part_glmer, ~ .y(.x, data = counti_rt_trans))
disc_fit_rl <- map2(forms, list_part_glmer, ~ .y(.x, data = counti_rl_trans))

# write out the model fits
write_rds(disc_fit, path = here('results', 'disc_fit.rds'))
write_rds(disc_fit_rt, path = here('results', 'disc_fit_rt.rds'))
write_rds(disc_fit_rl, path = here('results', 'disc_fit_rl.rds'))
