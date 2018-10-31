# data manipulation
library(tidyverse)
library(janitor)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/helper03_stan_utility.r')

# misc
library(pROC)
library(ROCR)
library(ggridges)
source('../R/helper01_process_foo.r')
source('../R/helper02_plot_foo.r')

# important constants
options(mc.cores = parallel::detectCores())

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)

# fold data prep
counti_trans <- counti_trans %>%
  mutate(fold = break_my(mybin, number = 5))

# 1 is young, 5 is old
counti_fold <- rev(split(counti_trans, counti_trans$fold))  # each fold


# get in fits and posterior work
model_key <- c('Past and vary', 
               'Past but no vary',
               'No past but vary', 
               'No past or vary')
fit <- read_rds('../data/training_fit.rds')

# given trained models for the rolled-out folds, estimate the next fold

# 16 models, 1:4 as formula for each fold etc
counti_fold_match <- rep(counti_fold[-1], 4)

# given new cutpoints, try predicting the test data
pred <- map2(fit, counti_fold_match,
             ~ posterior_linpred(object = .x, 
                                 newdata = .y,
                                 draws = 1000))


# calculate ROC/AUC for the predictions of test data
# to determine how good our out of sample predictions are
pred_auc <- map2(pred, counti_fold_match, post_roc) %>%
  map(., function(x) map_dbl(x, ~ auc(.x))) %>%
  set_names(., names(fit))


# now to make a plot about out-of-sample predictive accuracy
# because 1000s of numbers are hard to visualize

oos_auc <- bind_cols(pred_auc) %>% 
  gather() %>%
  separate(key, into = c('mod', 'fold'), sep = '\\_') %>%
  mutate(mod = plyr::mapvalues(mod, unique(mod), model_key),
         mod = factor(mod, levels = model_key)) %>%
  ggplot(aes(x = value, y = mod)) +
  geom_density_ridges(rel_min_height = 0.01) +
  labs(x = 'AUC', y = 'NULL') +
  scale_colour_brewer() +
  NULL
ggsave(filename = '../results/figure/fold_auc.png', plot = oos_auc,
       width = 6, height = 6)


# then do it for time
time_auc <- map2(pred, counti_fold_match, get_auc_time) %>%
  set_names(., names(fit))

# could this be done with enframe?
ta <- reshape2::melt(time_auc) %>% 
  as.tibble %>%
  mutate(time = parse_double(L2)) %>%
  separate(L1, into = c('mod', 'fold'), sep = '\\_') %>%
  mutate(mod = plyr::mapvalues(mod, unique(mod), model_key),
         mod = factor(mod, levels = rev(model_key))) %>%
  ggplot(aes(x = time, y = value)) +
  stat_lineribbon() +
  facet_grid(mod ~ .) +
  scale_fill_brewer() +
  scale_x_reverse()
ggsave(filename = '../results/figure/fold_auc_time.png', plot = ta,
       width = 8, height = 6)
