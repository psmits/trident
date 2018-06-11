# data manipulation
library(tidyverse)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)

# parallel processing
library(parallel)

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

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)
inspect <- counti_trans %>%            # nice to see the super relevant columns
  dplyr::select(event, relage, mybin, maxgcd, diff_maxgcd)

# read in model fits
disc_fit <- read_rds('../data/disc_fit.rds')

# estimate of out-of-sample performance
compare_loo <- map(disc_fit, loo)
compare_loo_tab <- loo::compare(x = compare_loo)
compare_waic <- map(disc_fit, waic)
compare_waic_tab <- loo::compare(x = compare_waic)

# limit
disc_best <- disc_fit[[4]]
# posterior probability of observation surviving
pp_est <- posterior_predict(disc_best)
pp_prob <- posterior_linpred(disc_best, transform = TRUE)
# adequacy of fit ROC
pp_auc <- apply(pp_est, 1, function(x) roc(counti_trans$event, x))


# break up by point in time
tt <- split(pp_est, counti_trans$fact_mybin)
tt <- lapply(tt, function(x) matrix(x, nrow = 4000))
ee <- split(counti_trans$event, counti_trans$fact_mybin)
pp_auc_time <- map2(tt, ee, ~ apply(.x, 1, function(a) auc(roc(a, .y))))


# posterior intervals
interval_est <- disc_best %>%
  spread_samples(`(Intercept)`, 
                 maxgcd, 
                 diff_maxgcd,
                 `maxgcd:diff_maxgcd`,
                 `Sigma[fact_mybin:(Intercept),(Intercept)]`,
                 `Sigma[fact_relage:(Intercept),(Intercept)]`,
                 `Sigma[fossil.group:(Intercept),(Intercept)]`) %>%
  median_qi(.prob = c(0.9, 0.5))

effect_eye <- disc_best %>%
  gather_samples(`(Intercept)`, 
                 maxgcd, 
                 diff_maxgcd,
                 `maxgcd:diff_maxgcd`) %>%
  ggplot(aes(y = term, x = estimate)) +
  geom_halfeyeh(.prob = c(0.9, 0.5))

vary_eye <- disc_best %>%
  gather_samples(`Sigma[fact_mybin:(Intercept),(Intercept)]`,
                 `Sigma[fact_relage:(Intercept),(Intercept)]`,
                 `Sigma[fossil.group:(Intercept),(Intercept)]`) %>%
  ggplot(aes(y = term, x = estimate)) +
  geom_halfeyeh(.prob = c(0.9, 0.5))


# base-line hazard plot
hazard_plot <- 
  disc_best %>%
  spread_samples(b[i, f], `(Intercept)`) %>%
  filter(str_detect(f, pattern = 'fact_mybin')) %>%
  mutate(cr = invlogit(`(Intercept)` + b)) %>%
  mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
  ggplot(aes(y = cr, x = age)) + 
  stat_lineribbon() +
  scale_fill_brewer() +
  labs(x = 'age (My)', y = 'P(T = t | T >= t, x)')
