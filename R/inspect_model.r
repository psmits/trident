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
  dplyr::select(event, relage, mybin, maxgcd, diff_maxgcd, lag1_maxgcd)

# read in model fits
disc_fit <- read_rds('../data/disc_fit.rds')
tree_fit <- read_rds('../data/tree_fit.rds')



# glmer model
# posterior intervals
interval_est <- disc_fit %>%
  spread_samples(`(Intercept)`, 
                 maxgcd, 
                 diff_maxgcd, 
                 lag1_maxgcd, 
                 `Sigma[fact_relage:(Intercept),(Intercept)]`) %>%
  median_qi(.prob = c(0.9, 0.5))

interval_eye <- disc_fit %>%
  gather_samples(maxgcd, diff_maxgcd, lag1_maxgcd) %>%
  ggplot(aes(y = term, x = estimate)) +
  geom_halfeyeh(.prob = c(0.9, 0.5))

interval_range <- disc_fit %>%
  gather_samples(`(Intercept)`, 
                 maxgcd, 
                 diff_maxgcd, 
                 lag1_maxgcd) %>%
  median_qi(.prob = c(0.9, 0.5)) %>%
  ggplot(aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_pointintervalh()


# base-line hazard plot, measured as continuation ratio
hazard_plot <- 
  disc_fit %>%
  spread_samples(b[i, f], `(Intercept)`) %>%
  mutate(cr = exp(`(Intercept)` + b)) %>%
  median_qi(.prob = 0.9) %>%
  mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
  ggplot(aes(y = cr, x = age, 
             ymin = cr.low, ymax = cr.high)) +
  geom_pointrange(fatten = 2) +
  geom_line()




# posterior probability of observation surviving
pp_est <- posterior_predict(disc_fit)
pp_prob <- posterior_linpred(disc_fit, transform = TRUE)





# tree-based model
rpart.plot(tree_fit, type = 2)
tree_fit_party <- as.party(tree_fit)
tree_fit_party$fitted[["(response)"]]<- Surv(counti$time1, counti$time2, counti$event)
plot(tree_fit_party)
