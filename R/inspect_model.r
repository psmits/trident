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
                   .funs = ~ log1p(.x)) %>%
dplyr::mutate_at(.vars = vars(ncell, latext, maxgcd), 
                 .funs = ~ .x - mean(.x))
longi <- longi %>% 
  mutate_at(.vars = vars(keel, symb, spin),
            .funs = ~ factor(.x, levels = c(0, 1)))
counti <- counti %>%
  mutate_at(.vars = vars(keel, symb, spin),
            .funs = ~ factor(.x, levels = c(0, 1)))
            #.funs = ~ as.numeric(.x))

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
                   time = max(time1),
                   time1 = max(time1),
                   time2 = max(time2),
                   event = max(event),
                   id = unique(id),
                   eco = unique(eco),
                   keel = unique(keel),
                   symb = unique(symb),
                   spin = unique(spin))

pp.s <- posterior_survfit(fit, newdataLong = longi, newdataEvent = nsu) 
pp.l <- posterior_predict(fit)

# simuates survival probability for each species for up to max observed age
# use bayes plot to compare
pmean <- ppc_stat(longi$maxgcd, pp.l, 'mean')
psd <- ppc_stat(longi$maxgcd, pp.l, 'sd')
pmed <- ppc_stat(longi$maxgcd, pp.l, 'median')
pmeang <- ppc_stat_grouped(longi$maxgcd, pp.l, 
                           group = longi$id, 'mean')
psdg <- ppc_stat_grouped(longi$maxgcd, pp.l, 
                         group = longi$id, 'sd')
pmedg <- ppc_stat_grouped(longi$maxgcd, pp.l, 
                          group = longi$id, 'median')
pe <- ppc_ecdf_overlay(longi$maxgcd, pp.l[1:50, ])
pd <- ppc_dens_overlay(longi$maxgcd, pp.l[1:50, ])


# check posterior predict for species through time
# get residuals from that


# what are we really interested in?
#   association parameters
#   ns(relage) because nonlinear?
