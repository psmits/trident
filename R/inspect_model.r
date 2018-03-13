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
fit <- read_rds('../data/jm_fit.rds')

# data transforms
longi <- longi %>%
  dplyr::mutate_at(.vars = vars(ncell, latext, maxgcd), 
                   .funs = ~ arm::rescale(log1p(.x)))
counti$cc.rescale <- with(counti, {
                            factor(cc.rescale, 
                            levels = sort(unique(as.numeric(cc.rescale))))})

# check the model...
# apparently the wrapper doesn't work for counting process formulation
# only one row per observation...key is i need to provide cc.rescale for each
# worry is what exactly am a simulating?
nsu <- counti %>% 
  group_by(fullname) %>%
  dplyr::summarize(cc.rescale = unique(cc.rescale),
                   time1 = max(time1),
                   time2 = max(time2),
                   event = max(event),
                   id = unique(id),
                   eco = unique(eco),
                   keep = unique(keel),
                   symb = unique(symb),
                   spin = unique(spin))

pp <- posterior_survfit(fit, newdataLong = longi, newdataEvent = nsu)
# use bayes plot to compare


# what are we really interested in?
#   association parameters
#   ns(relage) because nonlinear?
