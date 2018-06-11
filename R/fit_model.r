# data manipulation
library(tidyverse)

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

# transform the data for analysis
counti <- prepare_analysis(counti)


# fit the model
form <- formula(event ~ maxgcd * diff_maxgcd +
                (1 | fact_relage))
forms <- list(form,
             update(form, ~ . + (1 | fact_mybin)),
             update(form, ~ . + (1 | fossil.group)),
             update(form, ~ . + (1 | fact_mybin) + (1 | fossil.group)))
forms2 <- map(forms, ~ update(.x, ~ . - diff_maxgcd - maxgcd:diff_maxgcd))
forms <- append(forms, forms2)

disc_fit <- map(forms, ~ stan_glmer(.x, family = 'binomial', data = counti,
                                    adapt_delta = 0.999))
write_rds(disc_fit, path = '../data/disc_fit.rds')



# tree-based methods
# survival tree with time-varying covariates
#counti <- as.data.frame(counti)  # this is fucking annoying but necessary
#form <- formula(Surv(time1, time2, event) ~ 
#                maxgcd + diff_maxgcd + relage + mybin + fossil.group)
#tree_fit <- LTRCART(formula = form, data = counti)
#write_rds(tree_fit, path = '../data/tree_fit.rds')
