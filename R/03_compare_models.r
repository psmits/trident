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

# get data in
counti <- read_rds(here('data', 'counting.rds'))
#counti_restrict_time <- read_rds(here('data', 'counting_restrict_time.rds'))
#counti_restrict_local <- read_rds(here('data', 'counting_restrict_local.rds'))

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)
#counti_trans_rt <- prepare_analysis(counti_restrict_time)
#counti_trans_rl <- prepare_analysis(counti_restrict_local)

# read in model fits
#model_key <- c('Past and vary', 
#               'Past but no vary',
#               'No past but vary', 
#               'No past or vary')
model_key <- c('VP', 'CP', 'V', 'C')
disc_fit <- read_rds(here('results', 'disc_fit.rds'))
#disc_fit_rt <- read_rds(here('results', 'disc_fit_rt.rds'))
#disc_fit_rl <- read_rds(here('results', 'disc_fit_rl.rds'))


# bayes R2
br2 <- plot_bayesr2(disc_fit, model_key)
ggsave(filename = here('results', 'figure', 'bayes_r2_full.png'),
       plot = br2,
       width = 6, height = 6)
#br2_rt <- plot_bayesr2(disc_fit_rt, model_key)
#ggsave(filename = here('results', 'figure', 'bayes_r2_restrict_time.png'),
#       plot = br2_rt,
#       width = 6, height = 6)
#br2_rl <- plot_bayesr2(disc_fit_rl, model_key)
#ggsave(filename = here('results', 'figure', 'bayes_r2_restrict_local.png'),
#       plot = br2_rl,
#       width = 6, height = 6)


# look at posterior predictive distribution for...
# model
p_model(fit_list = disc_fit,
        .data = counti_trans,
        key = model_key,
        name = 'full',
        path = here('results', 'figure'))
#p_model(fit_list = disc_fit_rt,
#        .data = counti_trans_rt,
#        key = model_key,
#        name = 'restrict_time',
#        path = here('results', 'figure'))
#p_model(fit_list = disc_fit_rl,
#        .data = counti_trans_rl,
#        key = model_key,
#        name = 'restrict_local',
#        path = here('results', 'figure'))


# model X time
p_model_time(fit_list = disc_fit,
             .data = counti_trans,
             key = model_key,
             name = 'full',
             path = here('results', 'figure'))
#p_model_time(fit_list = disc_fit_rt,
#             .data = counti_trans_rt,
#             key = model_key,
#             name = 'restrict_time',
#             path = here('results', 'figure'))
#p_model_time(fit_list = disc_fit_rl,
#             .data = counti_trans_rl,
#             key = model_key,
#             name = 'restrict_local',
#             path = here('results', 'figure'))

# model X taxon
p_model_taxon(fit_list = disc_fit,
              .data = counti_trans,
              key = model_key,
              name = 'full',
              path = here('results', 'figure'))
#p_model_taxon(fit_list = disc_fit_rt,
#              .data = counti_trans_rt,
#              key = model_key,
#              name = 'restrict_time',
#              path = here('results', 'figure'))
#p_model_taxon(fit_list = disc_fit_rl,
#              .data = counti_trans_rl,
#              key = model_key,
#              name = 'restrict_local',
#              path = here('results', 'figure'))

# model X taxon/time
p_model_taxon_time(fit_list = disc_fit,
                   .data = counti_trans,
                   key = model_key,
                   name = 'full',
                   path = here('results', 'figure'))
#p_model_taxon_time(fit_list = disc_fit_rt,
#                   .data = counti_trans_rt,
#                   key = model_key,
#                   name = 'restrict_time',
#                   path = here('results', 'figure'))
#p_model_taxon_time(fit_list = disc_fit_rl,
#                   .data = counti_trans_rl,
#                   key = model_key,
#                   name = 'restrict_local',
#                   path = here('results', 'figure'))
