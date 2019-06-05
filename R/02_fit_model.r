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
source(here('R', 'helper07_diffsafe.r'))

# important constants
options(mc.cores = parallel::detectCores())


# get data in
counti <- read_rds(here('data', 'counting.rds'))
counti_restrict_time <- read_rds(here('data', 'counting_restrict_time.rds'))
counti_restrict_local <- read_rds(here('data', 'counting_restrict_local.rds'))

# transform the data for analysis
counti_trans <- prepare_analysis(counti)
counti_rt_trans <- prepare_analysis(counti_restrict_time)
counti_rl_trans <- prepare_analysis(counti_restrict_local)

## test eff
#test_maxgcd <- 
#  stan_glmer(formula = event ~ maxgcd + (1 | fact_mybin / fossil_group),
#             data = counti_trans,
#             family = 'binomial')
#
#test_mst <- 
#  stan_glmer(formula = event ~ mst + (1 | fact_mybin / fossil_group),
#             data = counti_trans,
#             family = 'binomial')
#
#loo(test_maxgcd)
#loo(test_mst)

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
form <- 
  formula(event ~ 
          temp + lag1_temp + mst + diff1_mst + diff2_mst + diff3_mst +
          (temp + lag1_temp + mst + diff1_mst + diff2_mst + diff3_mst | fact_mybin/fossil_group) +
          (1 | fact_relage/fossil_group))

# past but no vary
form2 <- 
  formula(event ~ 
          temp + lag1_temp + mst + diff1_mst + diff2_mst + diff3_mst +
          (1 | fact_mybin/fossil_group) +
          (1 | fact_relage/fossil_group))

# no past but vary
form3 <- 
  formula(event ~ 
          temp + mst + 
          (temp + mst | fact_mybin/fossil_group) +
          (1 | fact_relage/fossil_group))

# no past or vary
form4 <- 
  formula(event ~ 
          temp + mst + 
          (1 | fact_mybin/fossil_group) +
          (1 | fact_relage/fossil_group))


forms <- list(form, form2, form3, form4)

# partial evaluation to fill in priors and ancillary details?
part_glmer <- partial(stan_glmer, 
                      family = 'binomial',
                      prior_intercept = normal(-2, 5, autoscale = FALSE),
                      prior_aux = cauchy(0, 1, autoscale = FALSE),
                      prior_covariance = decov(regularization = 2),
                      thin = 4, 
                      adapt_delta = 0.999999)

# some models have differen't priors
list_part_glmer <- 
  list(partial(part_glmer, prior = normal(location = c(0, 0, -1, 0, 0, 0),
                                          scale = rep(1, 6),
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(location = c(0, 0, -1, 0, 0, 0),
                                          scale = rep(1, 6),
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(location = c(0, -1),
                                          scale = rep(1, 2),
                                          autoscale = FALSE)),
       partial(part_glmer, prior = normal(location = c(0, -1),
                                          scale = rep(1, 2),
                                          autoscale = FALSE)))


# map the data to the models
disc_fit <- map2(forms, list_part_glmer, ~ .y(.x, data = counti_trans))
disc_fit_rt <- map2(forms, list_part_glmer, ~ .y(.x, data = counti_rt_trans))
disc_fit_rl <- map2(forms, list_part_glmer, ~ .y(.x, data = counti_rl_trans))

# write out the model fits
write_rds(disc_fit, path = here('results', 'disc_fit.rds'))
write_rds(disc_fit_rt, path = here('results', 'disc_fit_rt.rds'))
write_rds(disc_fit_rl, path = here('results', 'disc_fit_rl.rds'))
