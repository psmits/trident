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
model_key <- c('Past and vary', 
               'No past but vary', 
               'No past or vary')
disc_fit <- read_rds('../data/disc_fit.rds')

# diagnostics
np <- map(disc_fit, nuts_params)
ll <- map(disc_fit, log_posterior)
# mcmc
rhat <- map(disc_fit, ~ mcmc_rhat(rhat(.x)))
neff <- map(disc_fit, ~ mcmc_neff(neff_ratio(.x)))
# hmc/nuts
accept <- map2(np, ll, ~ mcmc_nuts_acceptance(.x, .y))
diverg <- map2(np, ll, ~ mcmc_nuts_divergence(.x, .y))
steps <- map2(np, ll, ~ mcmc_nuts_stepsize(.x, .y))
treed <- map2(np, ll, ~ mcmc_nuts_treedepth(.x, .y))
energy <- map2(np, ll, ~ mcmc_nuts_energy(.x, .y))


# estimate of out-of-sample performance
# these should be parallelized; done using future
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

#pp_proc <- map(pp_prob, ~ apply(.x, 1, function(y) 
#                                roc(counti_trans$event, y, 
#                                    partial.auc = c(0, max(y))))

pp_auc <- map(pp_roc, function(y) map_dbl(y, ~ auc(.x)))

roc_hist <- bind_rows(imap(pp_auc, ~ data.frame(model = .y, roc = .x))) %>%
  mutate(model = recode(model, !!!model_key),
         model = fct_relevel(model, !!model_key[2], after = 1)) %>%
  ggplot(aes(x = roc, y = model)) +
  geom_eyeh()
ggsave(filename = '../doc/figure/roc_hist.png', plot = roc_hist,
       width = 6, height = 6)

# break up by point in time and plot
tt <- map(pp_est, ~ split(.x, counti_trans$fact_mybin))
tt <- map(tt, ~ map(.x, function(y) matrix(y, nrow = 1000)))
ee <- split(counti_trans$event, counti_trans$fact_mybin)
row_auc <- function(x, y) apply(x, 1, function(a) auc(roc(a, y)))
pat <- map(tt, ~ map2(.x, ee, function(a, b) try(row_auc(a, b))))


rr <- map(pat, ~ keep(.x, check_class)) %>%
  map(., ~ bind_rows(.x)) %>%
  bind_rows(.id = 'model') %>%
  gather(when, value, -model) %>%
  mutate(when = as.numeric(when)) %>%
  ggplot(aes(x = when, y = value)) +
  geom_eye() +
  facet_grid(~ model)
ggsave(filename = '../doc/figure/roc_time.png', plot = rr,
       width = 8, height = 8)
