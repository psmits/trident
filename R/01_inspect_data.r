# data manipulation
library(tidyverse)

# parallel processing
library(parallel)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/helper03_stan_utility.r')

# misc
library(scales)
library(ggridges)
library(pROC)
source('../R/helper01_process_foo.r')

theme_set(theme_bw())

# important constants
options(mc.cores = parallel::detectCores())

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)
#counti_trans %>% dplyr::select(fullname, mybin, relage, maxgcd, diff_maxgcd)


# occurrences through time, labeled if LAD
octg <-  
  counti %>%
  mutate(state = case_when(event == 0 ~ 'Standard',
                           event == 1 ~ 'Last'),
         fossil_group = case_when(fossil_group == 'D' ~ 'Dinoflagellates',
                                  fossil_group == 'R' ~ 'Radiolaria',
                                  fossil_group == 'F' ~ 'Foraminifera',
                                  fossil_group == 'N' ~ 'Calc. nanno.')) %>%
  ggplot(aes(x = mybin, fill = state)) +
  stat_bin() +
  facet_grid(fossil_group ~ ., switch = 'y') +
  scale_fill_manual(name = 'Occurrence type',
                    values = c('skyblue', 'goldenrod')) +
  theme(legend.position = 'bottom') +
  labs(title = 'Occurrences', x = 'Time (My before present)', y = 'Count')
ggsave(filename = '../results/figure/occ_time_label.png',
       plot = octg, width = 4, height = 6)


# relative "abundance" of microfossil_groups over time
ocag <- counti %>%
  ggplot(aes(x = mybin, fill = fossil_group)) +
  geom_histogram(position = 'fill')
ggsave(filename = '../results/figure/abn_time_stack.png',
       plot = ocag, width = 6, height = 6)

# occurrences by relage
ocrg <- counti %>% 
  group_by(fullname) %>%
  dplyr::summarize(maxage = max(relage),
                   fossil_group = plurality(fossil_group),
                   died = any(event == 1)) %>%
  mutate(state = case_when(died == 0 ~ 'Extant',
                           died == 1 ~ 'Extinct'),
         fossil_group = case_when(fossil_group == 'D' ~ 'Dinoflagellates',
                                  fossil_group == 'R' ~ 'Radiolaria',
                                  fossil_group == 'F' ~ 'Foraminifera',
                                  fossil_group == 'N' ~ 'Calc. nanno.')) %>%
  ggplot(aes(x = maxage, fill = state)) +
  stat_bin() +
  facet_grid(fossil_group ~ ., switch = 'y') +
  scale_fill_manual(name = 'State', 
                    values = c('skyblue', 'goldenrod')) +
  theme(legend.position = 'bottom') +
  labs(title = 'Age distribution', x = 'Age (My)', y = 'Count')
ggsave(filename = '../results/figure/age_label.png',
       plot = ocrg, width = 4, height = 6)



# make a plot of a random selection of species
set.seed(100)
srg <- counti %>% 
  mutate(logmaxgcd = log1p(maxgcd)) %>%
  group_by(fullname) %>%
  sample_n_groups(size = 8) %>%
  ungroup %>%
  ggplot(aes(x = relage, y = logmaxgcd, group = fullname, colour = fullname)) +
  geom_line() +
  geom_point() +
  theme(legend.position = 'bottom')
ggsave(filename = '../results/figure/range_time.png',
       plot = srg, width = 8, height = 6)

# lots of little code here
# for FAD/LAD accumulation curves
ft <- counti %>%
  group_by(fullname) %>%
  summarize(fad = max(mybin),
            fossil_group = plurality(fossil_group)) %>%
  group_by(fossil_group, fad) %>%
  summarize(n = n()) %>%
  arrange(desc(fad)) %>%
  mutate(nsum = cumsum(n),
         time = fad)

# LADs over time
lt <- counti %>%
  group_by(fullname) %>%
  filter(!all(event == 0)) %>%
  summarize(lad = max(mybin),
            fossil_group = plurality(fossil_group)) %>%
  group_by(fossil_group, lad) %>%
  summarize(n = n()) %>%
  arrange(desc(lad)) %>%
  mutate(nsum = cumsum(n),
         time = lad)

# put them together
ccg <- bind_rows(ft, lt, .id = 'type') %>%
  mutate(type = plyr::mapvalues(type, 1:2, c('FAD', 'LAD'))) %>%
  ggplot(aes(x = time, y = nsum, colour = type, group = type)) +
  geom_line() +
  facet_grid(~ fossil_group) +
  labs(x = 'Time (My)', y = 'Cummulative count')
ggsave(filename = '../results/figure/fad_lad_count_wide.png',
       plot = ccg, width = 6, height = 3)
ccg2 <- ccg + facet_grid(fossil_group ~ ., switch = 'y', scales = 'free_y')
ggsave(filename = '../results/figure/fad_lad_count_tall.png',
       plot = ccg2, width = 4, height = 6)


# looking at the temperature data
mgca <- read_tsv('../data/cramer/cramer_temp.txt')
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
ggsave(filename = '../results/figure/cramer_temp.png',
       plot = tpg, width = 6, height = 3)
