library(pacman)

p_load(tidyverse, furrr, here,
       arm, rstanarm, bayesplot,
       scales, ggridges,
       pROC)

source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper02_plot_foo.r'))
source(here('R', 'helper03_misc_foo.r'))
source(here('R', 'helper04_stan_utility.r'))

# important constants
theme_set(theme_bw())
plan(multiprocess)


# get data in
counti <- read_rds(here('data', 'counting.rds'))
counti_rt <- read_rds(here('data', 'counting_restrict_time.rds'))
counti_rl <- read_rds(here('data', 'counting_restrict_local.rds'))

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)
counti_rt_trans <- prepare_analysis(counti_rt)
counti_rl_trans <- prepare_analysis(counti_rl)

view_neptune(.data = counti_trans, 
             name = 'full',
             path = here('results', 'figure'))
view_neptune(.data = counti_rt_trans, 
             name = 'restrict_time',
             path = here('results', 'figure'))
view_neptune(.data = counti_rl_trans, 
             name = 'restrict_local',
             path = here('results', 'figure'))


# looking at the temperature data
mgca <- read_tsv(here('data', 'cramer', 'cramer_temp.txt'))
names(mgca) <- str_to_lower(names(mgca))
names(mgca) <- str_remove_all(names(mgca), '[^[[:alnum:]]]')
mgca <- mgca %>%
  mutate(bin = break_my(age, by = 1)) %>%
  group_by(bin) %>%
  mutate(temp_bin = mean(temperature, na.rm = TRUE)) %>%
  ungroup

# this is stupid but i'm bad at tidyr gather/spread
a <- mgca %>% 
  dplyr::select(age, temperature)
b <- mgca %>% 
  dplyr::select(bin, temp_bin) %>% 
  mutate(bin = bin - 0.5) %>%
  rename(age = bin,
         temperature = temp_bin) 
mg <- bind_rows(a, b, .id = 'scheme')

# temperature over time based on Mg/Ca
tpg <- mg %>%
  mutate(scheme = plyr::mapvalues(scheme, 1:2, c('raw', 'bin'))) %>%
  ggplot(aes(x = age, y = temperature, colour = scheme)) + 
  geom_line() +
  scale_colour_manual(values = c('goldenrod', 'skyblue')) +
  labs(x = 'Time (My)', y = 'Temperature diff. from modern (C)')
ggsave(filename = here('results', 'figure', 'cramer_temp.png'),
       plot = tpg, width = 6, height = 3)






# number of observations per time bin as we change binning scheme
# calculate
#  number of observations
#  number of observations per species
#  geographic range of species
vary_width <- read_rds(here('data', 'counting_vary_binwidth.rds'))

inplace <- function(x) {
  xx <- log(x)
  out <- xx / length(xx)
  out
}

width_effect <- 
  vary_width %>%
  dplyr::select(-latext) %>%
  group_by(width) %>%
  mutate_at(vars(-width), funs(inplace(.))) %>%
  gather(key = 'key', value = 'value', -width) %>%
  ggplot(aes(x = width, y = value)) +
  geom_jitter(height = 0, width = 0.05, alpha = 0.01) +
  facet_grid(key ~ ., scales = 'free', switch = 'y') + 
  labs(x = 'Timebin width', y = 'log(value) / n observations') +
  NULL
ggsave(filename = here('results', 'figure', 'binwidth_effect.png'),
       plot = width_effect, width = 8, height = 6)
