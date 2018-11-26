# data manipulation
library(tidyverse)
library(janitor)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)
library(deeptime)

# parallel processing
library(parallel)

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
source('../R/helper03_stan_utility.r')
source('../R/helper04_roc_utility.r')

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
  geom_halfeyeh(.width = c(0.5, 0.8)) +
  labs(x = 'AUC ROC', y = NULL) +
  scale_colour_brewer() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15)) +
  NULL
ggsave(filename = '../results/figure/fold_auc.png', plot = oos_auc,
       width = 6, height = 6)
ggsave(filename = '../results/figure/fold_auc_zoom.png', 
       plot = oos_auc + 
         xlim(0.5, 1) +
         geom_vline(xintercept = 0.5, colour = 'red'),
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
  scale_fill_brewer() +
  scale_x_reverse() +
  coord_cartesian(ylim = c(0.4, 1), xlim = c(0, 50)) +
  labs(y = 'AUC ROC', x = 'Time (Mya)') +
  NULL
ggsave(filename = '../results/figure/fold_auc_time_tiny.png', 
       plot = gggeo_scale(ta, 
                          dat = 'epochs', 
                          size = 3, 
                          rot = 90,
                          height = 0.2),
       width = 8, height = 6)
ta <- ta +
  facet_grid(mod ~ .) +
  geom_hline(yintercept = 0.5, colour = 'red', linetype = 'dashed') +
  NULL
ggsave(filename = '../results/figure/fold_auc_time.png', 
       plot = ta,
       width = 8, height = 6)


# by taxonomic group
pp_taxon <- map2(pred, counti_fold_match, 
                 ~ split(data.frame(t(.x)), 
                         .y$fossil_group)) %>%
  map(., function(x) map(x, ~ t(.x)))

counti_fold_taxon <- map(counti_fold_match, ~ split(.x, .x$fossil_group))

fold_auc_taxon <- map2(pp_taxon, counti_fold_taxon, 
     ~ map2(.x, .y, 
            ~ apply(.x, 1, function(a) fast_auc(a, .y$event)))) %>%
  reshape2::melt(.) %>%
  as.tibble(.) %>%
  rename(taxon = L2,
         model = L1) %>%
  separate(., col = model, into = c('model', 'fold'), sep = '_') %>%
  mutate(taxon = case_when(taxon == 'D' ~ 'Dinoflagellates',
                           taxon == 'R' ~ 'Radiolaria',
                           taxon == 'F' ~ 'Foraminifera',
                           taxon == 'N' ~ 'Calc. nanno.'),
         model = case_when(model == 'mod1' ~ model_key[1],
                           model == 'mod2' ~ model_key[2],
                           model == 'mod3' ~ model_key[3],
                           model == 'mod4' ~ model_key[4]),
         model = factor(model, levels = rev(model_key))) %>%
  ggplot(aes(x = value, y = model)) +
  geom_halfeyeh(.width = c(0.5, 0.8)) +
  facet_wrap(~ taxon) +
  labs(x = 'AUC ROC', y = NULL)
ggsave(filename = '../results/figure/fold_auc_taxon.png',
       plot = fold_auc_taxon,
       width = 8, height = 8)


# by taxon and time

