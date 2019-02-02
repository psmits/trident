# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)

# misc
library(scales)
library(ggridges)
library(pROC)
library(here)

source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper02_plot_foo.r'))
source(here('R', 'helper03_misc_foo.r'))
source(here('R', 'helper04_stan_utility.r'))

theme_set(theme_bw())

# important constants
options(mc.cores = parallel::detectCores())

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
