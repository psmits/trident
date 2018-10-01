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
source('../R/helper03_stan_utility.R')

# misc
library(pROC)
library(ROCR)
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
counti_accum <- accumulate(counti_fold, bind_rows) # combine some folds
counti_accum <- counti_accum[-1]


# get in fits and posterior work
fit <- read_rds('../data/training_fit.rds')

# given trained models for the rolled-out folds, estimate the next fold

# first: get the cutpoint from the training data set b/c class imbalance
plin <- map2(fit, counti_accum, ~ posterior_linpred(object = .x, newdata = .y))
cut_points <- map2(plin, counti_accum, ~ get_cutpoints(.x, .y$event)) %>%
  map(., ~ .x[[1]][3])

# second: given new cutpoints, predict test data class
pred <- map2(fit, counti_fold[-1], 
             ~ posterior_linpred(object = .x, newdata = .y)) %>%
  cut_newpoint(., cut_points)

# third: calculate ROC/AUC for the predictions of test data
# to determine how good our out of sample predictions are
pred_auc <- map2(pred, counti_fold[-1], post_roc) %>%
  map(., function(x) map_dbl(x, ~ auc(.x)))

# now to make a plot about out-of-sample predictive accuracy
# because 1000s of numbers are hard to visualize
oos_auc <- bind_cols(pred_auc) %>% 
  gather() %>%
  mutate(key = str_extract(key, '[0-9]')) %>%
  ggplot(aes(x = value, fill = key)) + 
  stat_density(alpha = 0.5) +
  labs(x = 'AUC', y = 'density')
ggsave(filename = '../doc/figure/fold_auc.png', plot = oos_auc,
       width = 6, height = 6)

# then do it for time
time_auc <- map2(pred, counti_fold[-1], get_auc_time)

# could this be done with enframe?
ta <- reshape2::melt(time_auc) %>% 
  as.tibble %>%
  mutate(fold = factor(L1),
         time = parse_double(L2)) %>%
  ggplot(aes(x = time, y = value, colour = fold)) +
  stat_interval(alpha = 0.5, .width = c(0.5, 0.8))
ggsave(filename = '../doc/figure/fold_auc_time.png', plot = ta,
       width = 8, height = 6)
