# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/stan_utility.R')

# misc
source('../R/process_foo.r')

# important constants
options(mc.cores = parallel::detectCores())

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# transform the data for analysis
counti <- prepare_analysis(counti)

# fit the model
form <- formula(event ~ maxgcd * diff_maxgcd + 
                temp + lag1_temp +
                (1 | fact_relage) + 
                (1 | fossil.group) + 
                (1 | fact_mybin))
form2 <- update(form, 
                ~ . - diff_maxgcd - maxgcd:diff_maxgcd - lag1_temp)
forms <- list(form, form2)

disc_fit <- map(forms, ~ stan_glmer(.x, family = 'binomial', data = counti,
                                    adapt_delta = 0.999, thin = 4))
write_rds(disc_fit, path = '../data/disc_fit.rds')
