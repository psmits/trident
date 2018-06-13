# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/stan_utility.R')

# misc
library(scales)
library(ggridges)
library(pROC)
source('../R/process_foo.r')

# important constants
options(mc.cores = parallel::detectCores())

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)


# occurrences through time, labeled if LAD
octg <- counti %>%
  ggplot(aes(x = mybin, fill = factor(event))) +
  stat_bin() +
  facet_grid(fossil.group ~ ., switch = 'y')

# relative "abundance" of microfossil groups over time
ocag <- counti %>%
  ggplot(aes(x = mybin, fill = fossil.group)) +
  geom_histogram(position = 'fill')

# occurrences by relage
ocrg <- counti %>% 
  group_by(fullname) %>%
  dplyr::summarize(maxage = max(relage),
                   fossil.group = plurality(fossil.group),
                   died = any(event == 1)) %>%
  ggplot(aes(x = maxage, fill = died)) +
  stat_bin() +
  facet_grid(fossil.group ~ ., switch = 'y')

# make a plot of a random selection of species
srg <- counti %>% 
  mutate(logmaxgcd = log1p(maxgcd)) %>%
  group_by(fullname) %>%
  sample_n_groups(size = 8) %>%
  ungroup %>%
  ggplot(aes(x = relage, y = logmaxgcd, group = fullname, colour = fullname)) +
  geom_line() +
  geom_point() +
  theme(legend.position = 'bottom')

# lots of little code here
# for FAD/LAD accumulation curves
ft <- counti %>%
  group_by(fullname) %>%
  summarize(fad = max(mybin),
            fossil.group = plurality(fossil.group)) %>%
  group_by(fossil.group, fad) %>%
  summarize(n = n()) %>%
  arrange(desc(fad)) %>%
  mutate(nsum = cumsum(n),
         time = fad)

# LADs over time
lt <- counti %>%
  group_by(fullname) %>%
  filter(!all(event == 0)) %>%
  summarize(lad = max(mybin),
            fossil.group = plurality(fossil.group)) %>%
  group_by(fossil.group, lad) %>%
  summarize(n = n()) %>%
  arrange(desc(lad)) %>%
  mutate(nsum = cumsum(n),
         time = lad)

# put them together
ccg <- bind_rows(ft, lt, .id = 'type') %>%
  mutate(type = plyr::mapvalues(type, 1:2, c('FAD', 'LAD'))) %>%
  ggplot(aes(x = time, y = nsum, colour = type, group = type)) +
  geom_line() +
  facet_grid(~ fossil.group) +
  labs(x = 'Time (My)', y = 'Cummulative count')

