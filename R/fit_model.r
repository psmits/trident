# data manipulation
library(tidyverse)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)

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

# data transforms
longi <- longi %>%
  dplyr::mutate_at(.vars = vars(ncell, latext, maxgcd), 
                   .funs = ~ arm::rescale(log1p(.x)))
longi <- longi %>%
  dplyr::filter(fg == 'F')

survi <- survi %>%
  dplyr::mutate(cc = parse_factor(cohort, levels = sort(unique(cohort))))
survi <- survi %>%
  dplyr::filter(fg == 'F')

# set up model fit
fit <- stan_jm(formulaLong = maxgcd ~ relage + (relage | id),
               dataLong = longi,
               formulaEvent = survival::Surv(duration, dead) ~ 1,
               dataEvent = survi,
               assoc = c('etavalue', 'etaslope'),
               time_var = 'relage',
               id_var = 'id', 
               chains = 4,
               cores = detectCores(),
               iter = 8000)

# check the model
pp <- posterior_predict(fit)



# what are we really interested in?
#   association parameters
#   ns(relage) because nonlinear?

