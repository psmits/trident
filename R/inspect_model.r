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
fit1 <- read_rds('../data/jm_fit1.rds')
fit2 <- read_rds('../data/jm_fit2.rds')
fit3 <- read_rds('../data/jm_fit3.rds')

# data transforms
longi <- longi %>%
  dplyr::mutate_at(.vars = vars(ncell, latext, maxgcd), 
                   .funs = ~ arm::rescale(log(.x)))
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

# survival process
#   keeps breaking for some reason -> issues with factors but i can't find it
pp_surv1 <- posterior_survfit(fit1, newdataLong = longi, newdataEvent = nsu)
pp_surv2 <- posterior_survfit(fit2, newdataLong = longi, newdataEvent = nsu)
pp_surv3 <- posterior_survfit(fit3, newdataLong = longi, newdataEvent = nsu)

# longitudinal process
pp_long1 <- posterior_predict(fit1)
pp_long2 <- posterior_predict(fit2)
pp_long3 <- posterior_predict(fit3)


# simuates survival probability for each species for up to max observed age
# use bayes plot to compare

pp_checks <- function(value, ppc, group) {
  pmean <- ppc_stat(value, ppc, 'mean')
  psd <- ppc_stat(value, ppc, 'sd')
  pmed <- ppc_stat(value, ppc, 'median')
  pmeang <- ppc_stat_grouped(value, ppc, 
                             group = group, 'mean')
  psdg <- ppc_stat_grouped(value, ppc,
                           group = group, 'sd')
  pmedg <- ppc_stat_grouped(value, ppc,
                            group = group, 'median')
  pe <- ppc_ecdf_overlay(value, ppc[1:50, ])
  pd <- ppc_dens_overlay(value, ppc[1:50, ])
  out <- list(pmean = pmean, psd = psd, pmed = pmed,
              pmeang = pmeang, psdg = psdg, pmedg = pmedg,
              pe = pe, pd = pd)
  out
}

surv_checks1a <- pp_checks(nsu$time1, pp_surv1, survi$cc.rescale)
surv_checks1b <- pp_checks(nsu$time2, pp_surv1, survi$cc.rescale)

long_checks1 <- pp_checks(longi$maxgcd, pp_long1, longi$id)
long_checks2 <- pp_checks(longi$maxgcd, pp_long2, longi$id)
long_checks3 <- pp_checks(longi$maxgcd, pp_long3, longi$id)



# i need a plot average hazard
#   really interested in knowing its shape
#   want to know if hazard (non)monotonic and/or age-dependent
#     these things don't necessarily go hand-in-hand
#     hazard can be nonmonotonic but not age-dependent if something else induces it (e.g. geographic range change)
