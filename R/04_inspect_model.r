# data manipulation
library(tidyverse)
library(janitor)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)

# parallel processing
library(parallel)
library(furrr)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)

# misc
library(pROC)
library(here)

source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper02_plot_foo.r'))
source(here('R', 'helper03_misc_foo.r'))
source(here('R', 'helper04_stan_utility.r'))
source(here('R', 'helper05_roc_utility.r'))

# important constants
future::plan(strategy = multicore)

# get data in
counti <- read_rds(here('data', 'counting.rds'))

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)

# read in model fits
disc_fit <- read_rds(here('results', 'disc_fit.rds')

# check how much information learned
priors <- prior_summary(disc_best)
prmn <- priors$prior$location
prsc <- priors$prior$scale

posc <- disc_best %>% 
  spread_draws(maxgcd,
               diff_maxgcd,
               temp,
               lag1_temp) %>% 
clean_names %>% 
dplyr::summarize_at(vars(-chain, -iteration), funs(sd(.)))

# 1 = lot of shrinkage, 0 = no shrinkage
posterior_shrinkage <- 1 - (posc / prsc) ^ 2


# posterior intervals
interval_est <- disc_best %>%
  spread_draws(`(Intercept)`, 
               maxgcd, 
               diff_maxgcd,
               #`maxgcd:diff_maxgcd`,
               temp,
               lag1_temp,
               `Sigma[fossil_group:fact_mybin:(Intercept),(Intercept)]`,
               `Sigma[fossil_group:fact_mybin:maxgcd,(Intercept)]`,
               `Sigma[fossil_group:fact_mybin:diff_maxgcd,(Intercept)]`,
               #`Sigma[fossil_group:fact_mybin:maxgcd:diff_maxgcd,(Intercept)]`,
               `Sigma[fossil_group:fact_mybin:maxgcd,maxgcd]`,
               `Sigma[fossil_group:fact_mybin:diff_maxgcd,maxgcd]`,
               #`Sigma[fossil_group:fact_mybin:maxgcd:diff_maxgcd,maxgcd]`,
               `Sigma[fossil_group:fact_mybin:diff_maxgcd,diff_maxgcd]`,
               #`Sigma[fossil_group:fact_mybin:maxgcd:diff_maxgcd,diff_maxgcd]`,
               #`Sigma[fossil_group:fact_mybin:maxgcd:diff_maxgcd,maxgcd:diff_maxgcd]`,
               `Sigma[fossil_group:fact_relage:(Intercept),(Intercept)]`,
               `Sigma[fact_mybin:(Intercept),(Intercept)]`,
               `Sigma[fact_mybin:maxgcd,(Intercept)]`,
               `Sigma[fact_mybin:diff_maxgcd,(Intercept)]`,
               #`Sigma[fact_mybin:maxgcd:diff_maxgcd,(Intercept)]`,
               `Sigma[fact_mybin:maxgcd,maxgcd]`,
               `Sigma[fact_mybin:diff_maxgcd,maxgcd]`,
               #`Sigma[fact_mybin:maxgcd:diff_maxgcd,maxgcd]`,
               `Sigma[fact_mybin:diff_maxgcd,diff_maxgcd]`,
               #`Sigma[fact_mybin:maxgcd:diff_maxgcd,diff_maxgcd]`,
               #`Sigma[fact_mybin:maxgcd:diff_maxgcd,maxgcd:diff_maxgcd]`,
               `Sigma[fact_relage:(Intercept),(Intercept)]`) %>%
median_qi(.width = c(0.8, 0.5))

effect_eye <- disc_best %>%
  gather_draws(`(Intercept)`, 
               maxgcd, 
               diff_maxgcd,
               #`maxgcd:diff_maxgcd`,
               temp,
               lag1_temp) %>%
ggplot(aes(y = .variable, x = .value)) +
geom_halfeyeh(.width = c(0.8, 0.5))
ggsave(filename = here('results', 'figure', 'effect_est.png'),
       plot = effect_eye, width = 4, height = 6)


# variance component
interp <- c('time: variance intercept, between',
            'time: variance maxgcd, between',
            'time: variance diff_maxgcd, between',
            'time: variance intercept, within',
            'time: variance maxgcd, within',
            'time: variance diff_maxgcd, within',
            'age: variance intercept, overall',
            'age: variance intercept, within')

vary_eye <- disc_best %>%
  gather_draws(`Sigma[fact_mybin:(Intercept),(Intercept)]`,
               `Sigma[fact_mybin:maxgcd,maxgcd]`,
               `Sigma[fact_mybin:diff_maxgcd,diff_maxgcd]`,
               `Sigma[fossil_group:fact_mybin:(Intercept),(Intercept)]`,
               `Sigma[fossil_group:fact_mybin:maxgcd,maxgcd]`,
               `Sigma[fossil_group:fact_mybin:diff_maxgcd,diff_maxgcd]`,
               `Sigma[fact_relage:(Intercept),(Intercept)]`,
               `Sigma[fossil_group:fact_relage:(Intercept),(Intercept)]`) %>%
  ungroup() %>%
  mutate(term = plyr::mapvalues(.variable, 
                                from = unique(.variable), 
                                to = interp)) %>%
  ggplot(aes(y = term, x = .value)) +
  geom_halfeyeh(.width = c(0.8, 0.5)) +
  labs(x = 'Estimate', y = 'Variance component')
ggsave(filename = here('results', 'figure', 'variance_components.png'),
       plot = vary_eye, width = 6, height = 6)


# probability of overall average ests being greater than 0
get_percent <- function(x) sum(x > 0) / length(x)
effect_prob <- disc_best %>%
  spread_draws(maxgcd, 
               diff_maxgcd, 
               #`maxgcd:diff_maxgcd`, 
               temp, 
               lag1_temp) %>%
summarize_at(.vars = c('maxgcd', 
                       'diff_maxgcd', 
                       # 'maxgcd:diff_maxgcd', 
                       'temp', 
                       'lag1_temp'), 
             .funs = get_percent)

# base-line hazard plot
hazard_plot <- disc_best %>%
  spread_draws(b[i, f], `(Intercept)`) %>%
  filter(str_detect(f, pattern = 'fact_relage'),
         !str_detect(f, pattern = 'fossil_group')) %>%
  mutate(cr = invlogit(`(Intercept)` + b)) %>%
  mutate(age = as.numeric(str_extract(f, '[0-9]+'))) %>%
  ggplot(aes(y = cr, x = age)) + 
  stat_lineribbon() +
  scale_fill_brewer() +
  labs(x = 'Age (My)', y = 'P(T = t | T >= t, x)')
ggsave(filename = here('results', 'figure', 'hazard_baseline.png'),
       plot = hazard_plot, width = 6, height = 4)

# can also then compare between taxonomic groups because i have that calculated
# faceted baseline hazard plot by fossil_group
db <- plot_taxon_hazard(disc_best)
ggsave(filename = here('results', 'figure', 'hazard_bygroup.png'), 
       plot = db,
       width = 6, height = 8)


# make plots of effect change over time for
#   3 covariates
#   4 taxonomic groups
by_taxon <- plot_taxon_covariate_time(disc_best)
ggsave(filename = here('results', 'figure', 'eff_time_group.png'), 
       plot = by_taxon,
       height = 6, width = 8)


# look at hazard over duration for a few examples
# grab a random selection of species
set.seed(100)
risk_plots <- plot_risk_time(counti_trans, disc_best, nsp = 4)
ggsave(filename = here('results', 'figure', 'relrisk_ext.png'), 
       plot = risk_plots[[1]],
       height = 6, width = 8)
ggsave(filename = here('results', 'figure', 'relrisk_range.png'), 
       plot = risk_plots[[2]],
       height = 6, width = 8)
