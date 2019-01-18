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
library(here)

# my function set
source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper02_plot_foo.r'))
source(here('R', 'helper03_misc_foo.r'))
source(here('R', 'helper04_stan_utility.r'))
source(here('R', 'helper05_roc_utility.r'))
source(here('R', 'helper06_geotime.r'))

# important constants
future::plan(strategy = multicore)

# get data in
counti_fold <- read_rds(here('data', 'counting.rds')) %>%
  prepare_analysis(.) %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)

counti_fold_rt <- read_rds(here('data', 'counting_restrict_time.rds')) %>%
  prepare_analysis(.) %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)

counti_fold_rl <- read_rds(here('data', 'counting_restrict_local.rds')) %>%
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
fit <- read_rds(here('results', 'training_fit.rds'))
fit_rt <- read_rds(here('results', 'training_fit_rt.rds'))
fit_rl <- read_rds(here('results', 'training_fit_rl.rds'))

# the following is called entirely for its side effects
# general summary
cv_model(fit, counti_fold_match, 
         model_key, 'full', 
         here('results', 'figure'))
cv_model(fit_rt, counti_fold_match_rt, 
         model_key, 'restrict_time',
         here('results', 'figure'))
cv_model(fit_rl, counti_fold_match_rl, 
         model_key, 'restrict_local',
         here('results', 'figure'))

# over time
cv_model_time(fit, counti_fold_match, 
              model_key, 'full',
              here('results', 'figure'))
cv_model_time(fit_rt, counti_fold_match_rt, 
              model_key, 'restrict_time',
              here('results', 'figure'))
cv_model_time(fit_rl, counti_fold_match_rl, 
              model_key, 'restrict_local',
              here('results', 'figure'))

# by taxon
cv_model_taxon(fit, counti_fold_match, 
               model_key, 'full',
               here('results', 'figure'))
cv_model_taxon(fit_rt, counti_fold_match_rt, 
               model_key, 'restrict_time',
               here('results', 'figure'))
cv_model_taxon(fit_rl, counti_fold_match_rl, 
               model_key, 'restrict_local',
               here('results', 'figure'))

# by taxon and time
cv_model_taxon_time(fit, counti_fold_match, 
                    model_key, 'full',
                    here('results', 'figure'))
cv_model_taxon_time(fit_rt, counti_fold_match_rt, 
                    model_key, 'restrict_time',
                    here('results', 'figure'))
cv_model_taxon_time(fit_rl, counti_fold_match_rl, 
                    model_key, 'restrict_local',
                    here('results', 'figure'))
