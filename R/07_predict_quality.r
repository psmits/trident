library(pacman)

p_load(tidyverse, magrittr, janitor, ggridges, furrr, tidybayes, 
       arm, rstanarm, bayesplot, splines, pROC, ROCR, here)

p_load_gh('willgearty/deeptime')

source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper02_plot_foo.r'))
source(here('R', 'helper03_misc_foo.r'))
source(here('R', 'helper04_stan_utility.r'))
source(here('R', 'helper05_roc_utility.r'))
source(here('R', 'helper06_geotime.r'))

# important constants
plan(multiprocess)

# get data in
counti <- 
  read_rds(here('data', 'counting.rds')) %>%
  prepare_analysis()
# cv data set
counti_fold <- 
  counti %>%
  mutate(fold = break_my(mybin, number = 5)) %>%
  split(., .$fold) %>%
  rev(.)
counti_fold_match <- rep(counti_fold[-1], 4)

# read in model fits
model_key <- c('VP', 'CP', 'V', 'C')
ins_fit <- read_rds(here('results', 'disc_fit.rds'))
oos_fit <- read_rds(here('results', 'training_fit.rds'))

# in-sample posterior predictive distribution
ins_roc <- ins_roc_taxon_time(ins_fit, counti, model_key)

# out-of-sample posterior predictive distribution
oos_roc <- oos_roc_taxon_time(oos_fit, counti_fold_match, model_key)


# compare auc estimate at each time point





