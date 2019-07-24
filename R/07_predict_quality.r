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
source(here('R', 'helper07_diffsafe.r'))

# important constants
plan(multiprocess)

# useful graphical constants
theme_set(theme_bw(base_size = 20))

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
# group_by taxon, time
# nest list columns
# somehow join by taxon/time
# map2 list columns
#  oos - ins
# plot that difference

ins_grp <- 
  ins_roc %>% 
  group_by_at(vars(model, phyla, time, fossil_group)) %>%
  group_nest()

oos_grp <- 
  oos_roc %>%
  group_by_at(vars(model, phyla, time, fossil_group)) %>%
  group_nest()

# join oos and ins dropping non-forecasted times
jnt <- inner_join(ins_grp, oos_grp, 
                  by = c('model', 'phyla', 'time', 'fossil_group'),
                  suffix = c('_ins', '_oos'))

# calculate difference in AUC estimate between oos and ins
auc_diff <- jnt %>% 
  mutate(diff = map2(data_ins, data_oos, ~ .y$value - .x$value)) %>%
  dplyr::select(-data_ins, -data_oos) %>%
  unnest()

# chart junk
rects <- get_geotime_box(range(auc_diff$time))

brks <- seq(min(auc_diff$time), max(auc_diff$time), by = 5) %>%
  round(., -1) %>%
  unique(.)


# plot those differences
diff_gg <- 
  auc_diff %>%
  ggplot() +
  geom_rect(data = rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = 'gray80', alpha = 0.8) +
  stat_lineribbon(mapping = aes(x = time, y = diff), 
                  size = 0, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
  scale_fill_brewer() +
  facet_grid(fossil_group ~ model) +
  scale_x_reverse(name = 'CI', breaks = brks) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(legend.position = 'bottom') +
  labs(x = 'Time (Mya)', y = 'Difference between out-of-sample and in-sample AUC') +
  NULL

ggsave(filename = here::here('results', 'figure', 'auc_diff.png'),
       plot = diff_gg,
       width = 11, height = 11)


