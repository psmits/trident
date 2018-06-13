# data manipulation
library(tidyverse)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)

# parallel processing
library(parallel)
library(future)

# survival
library(survival)
library(LTRCtrees)
library(rpart.plot)
library(partykit)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/stan_utility.R')

# misc
library(splines)
library(pROC)
source('../R/process_foo.r')

# important constants
options(mc.cores = parallel::detectCores())
future::plan(strategy = multicore)

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)

# read in model fits
disc_fit <- read_rds('../data/disc_fit.rds')
# estimate of out-of-sample performance
# these should be parallelized
cl <- map(disc_fit, ~ future::future(loo(.x)))
cl <- map(cl, ~ future::value(.x))
compare_loo_tab <- loo::compare(x = cl)
wl <- map(disc_fit, ~ future::future(waic(.x)))
wl <- map(wl, ~ future::value(.x))
compare_waic_tab <- loo::compare(x = wl)


# posterior probability of observation surviving
pe <- map(disc_fit, ~ future::future(posterior_predict(.x)))
pp_est <- map(pe, ~ future::value(.x))
pp <- map(disc_fit, ~ future::future(posterior_linpred(.x, transform = TRUE)))
pp_prob <- map(pe, ~ future::value(.x))

# adequacy of fit ROC
pp_roc <- map(pp_est, ~ apply(.x, 1, function(y) roc(counti_trans$event, y)))
pp_auc <- map(pp_roc, function(y) map_dbl(y, ~ auc(.x)))


# break up by point in time
#tt <- map(pp_est, ~ split(.x, counti_trans$fact_mybin))
#tt <- map(tt, ~ matrix(.x, nrow = 4000))
#ee <- split(counti_trans$event, counti_trans$fact_mybin)
#pp_auc_time <- mcMap(function(x, y) apply(x, 1, function(a) auc(roc(a, y))),
#                     x = tt, y = ee)
#pat <- map2(tt, ee, ~ future::future(apply(.x, 1, function(a) auc(roc(a, .y)))))
#pat <- map(pat, ~ future::value(.x))
