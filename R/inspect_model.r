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
tree_fit <- read_rds('../data/tree_fit.rds')

# estimate of out-of-sample performance
compare_loo <- map(disc_fit, loo)
compare_loo_tab <- loo::compare(x = compare_loo)
compare_waic <- map(disc_fit, waic)
compare_waic_tab <- loo::compare(x = compare_waic)

# posterior probability of observation surviving
pp_est <- map(disc_fit, posterior_predict)
#fitted_samples(disc_fit[[1]])
pp_prob <- map(disc_fit, ~ posterior_linpred(.x, transform = TRUE))
#predicted_samples(disc_fit[[1]])

# adequacy of fit ROC
pp_roc <- map(pp_est, ~ apply(.x, 1, function(y) auc(roc(counti_trans$event, y))))

# for each time bin
#pp_est
#counti_trans$fact_mybin


# pic a model
#disc_best <- disc_fit[[2]]
disc_best <- disc_fit[[4]]
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

interval_eye <- disc_best %>%
  gather_samples(`(Intercept)`, 
                 maxgcd, 
                 diff_maxgcd) %>%
  ggplot(aes(y = term, x = estimate)) +
  geom_halfeyeh(.prob = c(0.9, 0.5))

interval_range <- disc_best %>%
  gather_samples(`(Intercept)`, 
                 maxgcd, 
                 diff_maxgcd) %>%
  median_qi(.prob = c(0.9, 0.5)) %>%
  ggplot(aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_pointintervalh()


# base-line hazard plot
hazard_plot <- 
  disc_best %>%
  spread_samples(b[i, f], `(Intercept)`) %>%
  filter(str_detect(f, pattern = 'fact_mybin')) %>%
  mutate(cr = invlogit(`(Intercept)` + b)) %>%
  median_qi(.prob = 0.9) %>%
  mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
  ggplot(aes(y = cr, x = age, 
             ymin = cr.low, ymax = cr.high)) +
  geom_pointrange(fatten = 2) +
  geom_line()




# tree-based model
#rpart.plot(tree_fit, type = 2)
#tree_fit_party <- as.party(tree_fit)
#tree_fit_party$fitted[["(response)"]]<- Surv(counti$time1, counti$time2, counti$event)
#plot(tree_fit_party)
