# data manipulation
library(tidyverse)
library(janitor)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)
library(deeptime)

# parallel processing
library(parallel)
library(future)
library(furrr)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)

# misc
library(pROC)
library(ROCR)
library(ggridges)

# my function set
source('../R/helper01_process_foo.r')
source('../R/helper02_plot_foo.r')
source('../R/helper03_misc_foo.r')
source('../R/helper04_stan_utility.r')
source('../R/helper05_roc_utility.r')
source('../R/helper06_geotime.r')

# important constants
future::plan(strategy = multicore)

# get data in
counti_fold <- read_rds('../data/counting.rds') %>%
  prepare_analysis(.) %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)

counti_fold_rt <- read_rds('../data/counting_restrict_time.rds') %>%
  prepare_analysis(.) %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)

counti_fold_rl <- read_rds('../data/counting_restrict_local.rds') %>%
  prepare_analysis(.) %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)

counti_fold_match <- rep(counti_fold[-1], 4)
counti_fold_match_rt <- rep(counti_fold_rt[-1], 4)
counti_fold_match_rl <- rep(counti_fold_rl[-1], 4)

# get in fits and posterior work
model_key <- c('Past and vary', 
               'Past but no vary',
               'No past but vary', 
               'No past or vary')
fit <- read_rds('../data/training_fit.rds')
fit_rt <- read_rds('../data/training_fit_rt.rds')
fit_rl <- read_rds('../data/training_fit_rl.rds')

# general summary
cv_model(fit, counti_fold_match, 
         model_key, 'full', '../results/figure/')
cv_model(fit_rt, counti_fold_match_rt, 
         model_key, 'restrict_time', '../results/figure/')
cv_model(fit_rl, counti_fold_match_rl, 
         model_key, 'restrict_local', '../results/figure/')

# over time
cv_model_time(fit, counti_fold_match, 
              model_key, 'full', '../results/figure/')
cv_model_time(fit_rt, counti_fold_match_rt, 
              model_key, 'restrict_time', '../results/figure/')
cv_model_time(fit_rl, counti_fold_match_rl, 
              model_key, 'restrict_local', '../results/figure/')

# by taxon
cv_model_taxon(fit, counti_fold_match, 
               model_key, 'full', '../results/figure/')
cv_model_taxon(fit_rt, counti_fold_match_rt, 
               model_key, 'restrict_time', '../results/figure/')
cv_model_taxon(fit_rl, counti_fold_match_rl, 
               model_key, 'restrict_local', '../results/figure/')

# by taxon and time
cv_model_taxon_time(fit, counti_fold_match, 
                    model_key, 'full', '../results/figure/')
cv_model_taxon_time(fit_rt, counti_fold_match_rt, 
                    model_key, 'restrict_time', '../results/figure/')
cv_model_taxon_time(fit_rl, counti_fold_match_rl, 
                    model_key, 'restrict_local', '../results/figure/')
