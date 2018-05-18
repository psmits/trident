# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# analysis packages
library(survival)
library(LTRCtrees)
library(rpart.plot)
library(partykit)
library(splines)
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/stan_utility.R')

# important constants
options(mc.cores = parallel::detectCores())

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# necessary data transform to handle factor
counti$cc.rescale <- with(counti, {
                          factor(cc.rescale, 
                          levels = sort(unique(as.numeric(cc.rescale))))})

# prep the data
counti <- counti %>%
  arrange(fullname, desc(mybin)) %>%
  mutate_at(.vars = vars(ncell, latext, maxgcd, area),
            .funs = ~ arm::rescale(log1p(.x))) %>%
  group_by(fullname) %>%
  mutate(diff_maxgcd = c(0, diff(maxgcd)),
         lag1_maxgcd = lag(maxgcd, default = 0),
         lag2_maxgcd = lag(maxgcd, n = 2, default = 0),
         lead1_maxgcd = lead(maxgcd, default = 0),
         lead2_maxgcd = lead(maxgcd, n = 2, default = 0)) %>%
  ungroup()

counti$fact_mybin <- as.factor(counti$mybin)
counti$fact_relage <- as.factor(counti$relage)


# fit the model
form <- formula(event ~ relage + maxgcd + diff_maxgcd + lag1_maxgcd + 
                (1 | fact_mybin) + (1 | id))
disc_fit <- stan_glmer(form, family = 'binomial', data = counti)


# tree-based methods
# survival tree with time-varying covariates
counti <- as.data.frame(counti)  # this is fucking annoying but necessary
form <- formula(Surv(time1, time2, event) ~ 
                maxgcd + diff_maxgcd + lag1_maxgcd + relage)
tree_fit <- LTRCART(formula = form, data = counti)
#rpart.plot(tree_fit, type = 2)
