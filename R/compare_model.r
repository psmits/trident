# data manipulation
library(tidyverse)
library(janitor)
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
library(ROCR)
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


## estimate of out-of-sample performance
## these should be parallelized; done using future
#cl <- map(disc_fit, ~ future::future(loo(.x)))
#cl <- map(cl, ~ future::value(.x))
#compare_loo_tab <- loo::compare(x = cl)
#wl <- map(disc_fit, ~ future::future(waic(.x)))
#wl <- map(wl, ~ future::value(.x))
#compare_waic_tab <- loo::compare(x = wl)

# posterior probability of observation surviving
pe <- map(disc_fit, ~ future::future(posterior_predict(.x)))
pp_est <- map(pe, ~ future::value(.x))
pp <- map(disc_fit, ~ future::future(posterior_linpred(.x, transform = TRUE)))
pp_prob <- map(pp, ~ future::value(.x))

# adequacy of fit ROC
opt_cut <- function(perf, pred) {
  cut_ind <- mapply(FUN = function(x, y, p) {
                      d <- (x - 0)^2 + (y - 1)^2
                      ind <- which(d == min(d))
                      c(sensitivity = y[[ind]], 
                        specificty = 1 - x[[ind]],
                        cutoff = p[[ind]])
               }, perf@x.values, perf@y.values, pred@cutoffs)
  cut_ind
}
get_cutpoints <- function(pp_prob, target) {
  pred <- plyr::alply(pp_prob, 1, function(x) prediction(x, target))
  perf <- map(pred, ~ performance(.x, measure = 'tpr', x.measure = 'fpr'))
  cutpoints <- map2(perf, pred, opt_cut)
}
pgc <- partial(get_cutpoints, target = counti_trans$event)
list_cutpoint <- map(map(pp_prob, pgc), ~ map(.x, ~ .x[3]))

pp_est_new <- list()
for(jj in seq(length(pp_prob))) {
  mm <- matrix(ncol = ncol(pp_prob[[jj]]), nrow = nrow(pp_prob[[jj]]))
  for(ii in seq(nrow(mm))) {
    mm[ii, ] <- pp_prob[[jj]][ii, ] > list_cutpoint[[jj]][[ii]]
  }
  pp_est_new[[jj]] <- mm * 1
}
# this resets the cutpoint for determining 0 vs 1 by setting it very low.
# this is a product of really high class imbalance. 
# cutpoint based on maximizing both sensitivity and specificity.
# this is calculated for every posterior predictive simulation for each model.
# these new, rescaled 0-1 results are then fed through the ROC/AUC machine.

pp_roc <- map(pp_est_new, ~ apply(.x, 1, function(y) roc(counti_trans$event, y)))
pp_auc <- map(pp_roc, function(y) map_dbl(y, ~ auc(.x)))

roc_hist <- bind_rows(imap(pp_auc, ~ data.frame(model = .y, roc = .x))) %>%
  mutate(model = recode(model, !!!model_key),
         model = fct_relevel(model, !!model_key[2], after = 1)) %>%
  ggplot(aes(x = roc, y = model)) +
  geom_eyeh()
ggsave(filename = '../doc/figure/roc_hist.png', plot = roc_hist,
       width = 6, height = 6)

# roc as timeseries to see best and worst times
roc_ts <- plot_roc_series(counti_trans, pp_est_new, model_key)
ggsave(filename = '../doc/figure/roc_ts.png', plot = roc_ts,
       width = 8, height = 8)
