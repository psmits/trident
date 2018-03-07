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
counti <- read_rds('../data/counting.rds')

# data transforms
longi <- longi %>%
  dplyr::mutate_at(.vars = vars(ncell, latext, maxgcd), 
                   .funs = ~ arm::rescale(log1p(.x)))


# set up model fit
fit <- stan_jm(formulaLong = maxgcd ~ relage + (relage | id),
               dataLong = longi,
               formulaEvent = survival::Surv(time = time1,
                                             time2 = time2,
                                             event = event,
                                             type = 'counting') ~ cc.rescale,
               dataEvent = counti,
               assoc = c('etavalue', 'etaslope'),
               time_var = 'relage',
               id_var = 'id', 
               chains = 4,
               cores = detectCores(),
               iter = 8000)

# check the model
pp <- posterior_predict(fit)           # posterior predictive distribution
# worry is what exactly am a simulating?

# use bayes plot to compare


# what are we really interested in?
#   association parameters
#   ns(relage) because nonlinear?

