library(tidyverse)
library(zoo)

source('../R/process_foo.r')
source('../R/sim_foo.r')

nsim <- 1000
nstep <- 100

simul <- rerun(.n = nsim, sim_range(nstep)) %>%
  imap(~ tibble(range = .x, 
                time = seq(length(.x)), 
                sim = .y)) %>%
  reduce(bind_rows)

# get rid of things that 
#   only extisted for their first moment sofilter by group size
#   never died
# standardize time
# time bins
simul <- simul %>%
  group_by(sim) %>%
  nest() %>%
  mutate(n = map(data, n_distinct)) %>%
  unnest(n = n) %>%
  filter(n > 1,
         n < nstep) %>%
  select(-n) %>%
  unnest() %>%
  group_by(sim) %>%
  mutate(time_standard = (time - min(time)) / (max(time) - min(time))) %>%
  ungroup() %>%
  mutate(time_bin = break_my(time_standard, number = 25)) %>%
  group_by(time_bin) %>%
  mutate(range_avg = mean(range)) %>%
  ungroup()

simul_avg <- simul %>%
  group_by(time_bin) %>%
  dplyr::summarize(time_standard = mean(time_standard),
                   range_avg = unique(range_avg)) %>%
  mutate(range_roll = rollmean(range_avg, 3, na.pad = TRUE))


# basic plot of runs
graph_basic <- simul %>%
  ggplot(aes(x = time, y = range, group = sim)) + geom_line()

# normalize by duration
graph_norm <- simul %>%
  ggplot(aes(x = time_standard, y = range, group = sim)) +
  geom_line() +
  geom_line(mapping = aes(x = time_standard, y = range_avg, group = NULL), 
            colour = 'blue', size = 1.2) +
  geom_line(mapping = aes(x = time_standard, y = range_roll, group = NULL),
            data = simul_avg, colour = 'red')

