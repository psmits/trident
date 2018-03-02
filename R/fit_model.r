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

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')

# data transforms
longi <- longi %>%
  dplyr::mutate_at(.vars = vars(ncell, latext, maxgcd), 
                   .funs = ~ arm::rescale(log1p(.x)))

survi <- survi %>%
  dplyr::mutate(cc = parse_factor(cohort, levels = sort(unique(cohort))))

# set up model fit
fit <- stan_jm(formulaLong = maxgcd ~ ns(relage, 2) + (ns(relage, 2) | id),
               dataLong = longi,
               formulaEvent = survival::Surv(duration, dead) ~ fg,
               dataEvent = survi,
               assoc = c('etavalue'),
               time_var = 'relage',
               id_var = 'id', 
               chains = 4,
               cores = detectCores(),
               iter = 500)

# what are we really interested in?
#   association parameters
#   ns(relage) because nonlinear?
