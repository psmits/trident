# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# analysis packages
library(survival)
library(splines)
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/stan_utility.R')


# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# data transforms
longi <- longi %>%
  dplyr::mutate_at(.vars = vars(ncell, latext, maxgcd), 
                   .funs = ~ arm::rescale(log(.x)))
counti$cc.rescale <- with(counti, {
                          factor(cc.rescale, 
                          levels = sort(unique(as.numeric(cc.rescale))))})

# this probably isn't necessary but i find it helpful to be explicit
longi <- longi %>% 
  mutate_at(.vars = vars(keel, symb, spin),
            .funs = ~ factor(.x, levels = c(0, 1)))
counti <- counti %>%
  mutate_at(.vars = vars(keel, symb, spin),
            .funs = ~ factor(.x, levels = c(0, 1)))
# should i model.matrix before sticking in function?


# set up model fit
#   could include extrinsic climate as time-dependent covariate of survival


# no covariates
fit1 <- stan_jm(formulaLong = maxgcd ~ ns(relage, 2) + (ns(relage, 2) | id),
             dataLong = longi,
             formulaEvent = survival::Surv(time = time1,
                                           time2 = time2,
                                           event = event,
                                           type = 'counting') ~ cc.rescale, 
           dataEvent = counti,
           time_var = 'relage',
           id_var = 'id', 
           assoc = c('etavalue', 'etaslope'),
           basehaz = 'weibull',
           chains = 4,
           cores = detectCores(),
           iter = 8000)
write_rds(fit1, path = '../data/jm_fit1.rds')

# covariates just for survival
fit2 <- stan_jm(formulaLong = maxgcd ~ ns(relage, 2) + (ns(relage, 2) | id),
             dataLong = longi,
             formulaEvent = survival::Surv(time = time1,
                                           time2 = time2,
                                           event = event,
                                           type = 'counting') ~ 
             cc.rescale + keel + symb + spin + eco,
           dataEvent = counti,
           time_var = 'relage',
           id_var = 'id', 
           assoc = c('etavalue', 'etaslope'),
           basehaz = 'weibull',
           chains = 4,
           cores = detectCores(),
           iter = 8000)
write_rds(fit2, path = '../data/jm_fit2.rds')

# all covariates
fit3 <- stan_jm(formulaLong = maxgcd ~ 
               ns(relage, 2) + keel + symb + spin + eco + (ns(relage, 2) | id),
             dataLong = longi,
             formulaEvent = survival::Surv(time = time1,
                                           time2 = time2,
                                           event = event,
                                           type = 'counting') ~ 
             cc.rescale + keel + symb + spin + eco,
           dataEvent = counti,
           time_var = 'relage',
           id_var = 'id', 
           assoc = c('etavalue', 'etaslope'),
           basehaz = 'weibull',
           chains = 4,
           cores = detectCores(),
           iter = 8000)
write_rds(fit3, path = '../data/jm_fit3.rds')
